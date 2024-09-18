/// multisearch: massively search of many queries against many large subjects.
/// the OG MAGsearch, branchwater, etc.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.
use anyhow::Result;
use rayon::prelude::*;
use stats::{median, stddev};
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use crate::utils::{csvwriter_thread, load_collection, load_sketches, ReportType, SearchResult};
use sourmash::ani_utils::ani_from_containment;
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
    let allow_empty_collection = false;
    // Load query collection
    let query_collection = load_collection(
        &query_filepath,
        selection,
        ReportType::Query,
        allow_failed_sigpaths,
        allow_empty_collection,
    )?;
    // load all query sketches into memory, downsampling on the way
    let query_sketchlist = load_sketches(query_collection, selection, ReportType::Query).unwrap();

    // Against: Load all _paths_, not signatures, into memory.
    let against_collection = load_collection(
        &against_filepath,
        selection,
        ReportType::Against,
        allow_failed_sigpaths,
        allow_empty_collection,
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
            if i % 1000 == 0 && i > 0 {
                eprintln!("Processed {} search sigs", i);
            }

            let mut results = vec![];

            // against downsampling happens here
            match against_collection.sig_from_record(record) {
                Ok(against_sig) => {
                    if let Some(against_mh) = against_sig.minhash() {
                        for query in query_sketchlist.iter() {
                            // to do - let user choose?
                            let calc_abund_stats = against_mh.track_abundance();

                            let against_mh_ds = against_mh.downsample_scaled(query.minhash.scaled()).unwrap();
                            let overlap =
                                query.minhash.count_common(&against_mh_ds, false).unwrap() as f64;

                            // only calculate results if we have shared hashes
                            if overlap > 0.0 {
                                let query_size = query.minhash.size() as f64;
                                let containment_query_in_target = overlap / query_size;
                                if containment_query_in_target > threshold {
                                    let target_size = against_mh.size() as f64;
                                    let containment_target_in_query = overlap / target_size;

                                    let max_containment =
                                        containment_query_in_target.max(containment_target_in_query);
                                    let jaccard = overlap / (target_size + query_size - overlap);

                                    let qani = ani_from_containment(
                                        containment_query_in_target,
                                        against_mh.ksize() as f64,
                                    );
                                    let mani = ani_from_containment(
                                        containment_target_in_query,
                                        against_mh.ksize() as f64,
                                    );
                                    let query_containment_ani = Some(qani);
                                    let match_containment_ani = Some(mani);
                                    let average_containment_ani = Some((qani + mani) / 2.);
                                    let max_containment_ani = Some(f64::max(qani, mani));

                                    let (total_weighted_hashes, n_weighted_found, average_abund, median_abund, std_abund) = if calc_abund_stats {
                                        match query.minhash.inflated_abundances(&against_mh_ds) {
                                            Ok((abunds, sum_weighted_overlap)) => {
                                                let sum_all_abunds = against_mh_ds.sum_abunds() as usize;
                                                let average_abund = sum_weighted_overlap as f64 / abunds.len() as f64;
                                                let median_abund = median(abunds.iter().cloned()).unwrap();
                                                let std_abund = stddev(abunds.iter().cloned());
                                                (Some(sum_all_abunds), Some(sum_weighted_overlap as usize), Some(average_abund), Some(median_abund), Some(std_abund))
                                            }
                                            Err(e) => {
                                                eprintln!("Error calculating abundances for query: {}, against: {}; Error: {}", query.name, against_sig.name(), e);
                                                continue;
                                            }
                                        }
                                    } else {
                                        (None, None, None, None, None)
                                    };

                                    results.push(SearchResult {
                                        query_name: query.name.clone(),
                                        query_md5: query.md5sum.clone(),
                                        match_name: against_sig.name(),
                                        containment: containment_query_in_target,
                                        intersect_hashes: overlap as usize,
                                        match_md5: Some(against_sig.md5sum()),
                                        jaccard: Some(jaccard),
                                        max_containment: Some(max_containment),
                                        average_abund,
                                        median_abund,
                                        std_abund,
                                        query_containment_ani,
                                        match_containment_ani,
                                        average_containment_ani,
                                        max_containment_ani,
                                        n_weighted_found,
                                        total_weighted_hashes,
                                    });
                                }
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
