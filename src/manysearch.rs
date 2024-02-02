/// multisearch: massively search of many queries against many large subjects.
/// the OG MAGsearch, branchwater, etc.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.
use anyhow::Result;
use rayon::prelude::*;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use crate::utils::{csvwriter_thread, load_collection, load_sketches, ReportType, SearchResult};
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;

pub fn manysearch(
    query_filepath: String,
    against_filepath: String,
    selection: &Selection,
    threshold: f64,
    output: Option<String>,
    allow_failed_sigpaths: bool,
) -> Result<()> {
    // Load query collection
    let query_collection = load_collection(
        &query_filepath,
        selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;
    // load all query sketches into memory, downsampling on the way
    let query_sketchlist = load_sketches(query_collection, selection, ReportType::Query).unwrap();

    // Against: Load all _paths_, not signatures, into memory.
    let against_collection = load_collection(
        &against_filepath,
        selection,
        ReportType::Against,
        allow_failed_sigpaths,
    )?;

    // set up a multi-producer, single-consumer channel.
    let (send, recv) = std::sync::mpsc::sync_channel::<SearchResult>(rayon::current_num_threads());

    // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = csvwriter_thread(recv, output);

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_sigs = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let send = against_collection
        .par_iter()
        .filter_map(|(_idx, record)| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 {
                eprintln!("Processed {} search sigs", i);
            }

            let mut results = vec![];

            // against downsampling happens here
            match against_collection.sig_from_record(record) {
                Ok(against_sig) => {
                    if let Some(against_mh) = against_sig.minhash() {
                        for query in query_sketchlist.iter() {
                            let overlap =
                                query.minhash.count_common(against_mh, false).unwrap() as f64;
                            let query_size = query.minhash.size() as f64;
                            let target_size = against_mh.size() as f64;
                            let containment_query_in_target = overlap / query_size;
                            let containment_in_target = overlap / target_size;
                            let max_containment =
                                containment_query_in_target.max(containment_in_target);
                            let jaccard = overlap / (target_size + query_size - overlap);

                            if containment_query_in_target > threshold {
                                results.push(SearchResult {
                                    query_name: query.name.clone(),
                                    query_md5: query.md5sum.clone(),
                                    match_name: against_sig.name(),
                                    containment: containment_query_in_target,
                                    intersect_hashes: overlap as usize,
                                    match_md5: Some(against_sig.md5sum()),
                                    jaccard: Some(jaccard),
                                    max_containment: Some(max_containment),
                                });
                            }
                        }
                    } else {
                        eprintln!(
                            "WARNING: no compatible sketches in path '{}'",
                            record.internal_location()
                        );
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                }
                Err(err) => {
                    eprintln!("Sketch loading error: {}", err);
                    eprintln!(
                        "WARNING: no compatible sketches in path '{}'",
                        record.internal_location()
                    );
                    let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
            }

            Some(results)
        })
        .flatten()
        .try_for_each_with(send, |s, m| s.send(m));

    // do some cleanup and error handling -
    if let Err(e) = send {
        eprintln!("Unable to send internal data: {:?}", e);
    }

    if let Err(e) = thrd.join() {
        eprintln!("Unable to join internal thread: {:?}", e);
    }

    // done!
    let i: usize = processed_sigs.fetch_max(0, atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} search sigs", i);

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    if skipped_paths > 0 {
        eprintln!(
            "WARNING: skipped {} search paths - no compatible signatures.",
            skipped_paths
        );
    }
    if failed_paths > 0 {
        eprintln!(
            "WARNING: {} search paths failed to load. See error messages above.",
            failed_paths
        );
    }

    Ok(())
}
