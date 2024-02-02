/// mastiff_manygather: mastiff-indexed version of fastmultigather.
use anyhow::Result;
use camino::Utf8PathBuf as PathBuf;
use rayon::prelude::*;
use sourmash::index::revindex::{RevIndex, RevIndexOps};
use sourmash::prelude::*;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use crate::utils::{
    csvwriter_thread, is_revindex_database, load_collection, BranchwaterGatherResult, ReportType,
};

pub fn mastiff_manygather(
    queries_file: String,
    index: PathBuf,
    selection: &Selection,
    threshold_bp: usize,
    output: Option<String>,
    allow_failed_sigpaths: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    if !is_revindex_database(&index) {
        bail!("'{}' is not a valid RevIndex database", index);
    }
    // Open database once
    let db = RevIndex::open(index, true)?;
    println!("Loaded DB");

    let query_collection = load_collection(
        &queries_file,
        selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;

    // set up a multi-producer, single-consumer channel.
    let (send, recv) =
        std::sync::mpsc::sync_channel::<BranchwaterGatherResult>(rayon::current_num_threads());

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

    let send = query_collection
        .par_iter()
        .filter_map(|(_idx, record)| {
            let threshold = threshold_bp / selection.scaled()? as usize;

            // query downsampling happens here
            match query_collection.sig_from_record(record) {
                Ok(query_sig) => {
                    let mut results = vec![];
                    if let Some(query_mh) = query_sig.minhash() {
                        // Gather!
                        let (counter, query_colors, hash_to_color) =
                            db.prepare_gather_counters(query_mh);

                        let matches = db.gather(
                            counter,
                            query_colors,
                            hash_to_color,
                            threshold,
                            query_mh,
                            Some(selection.clone()),
                        );
                        // extract results TODO: ADD REST OF GATHER COLUMNS
                        if let Ok(matches) = matches {
                            for match_ in &matches {
                                results.push(BranchwaterGatherResult {
                                    query_name: query_sig.name().clone(),
                                    query_md5: query_sig.md5sum().clone(),
                                    match_name: match_.name().clone(),
                                    match_md5: match_.md5().clone(),
                                    f_match_query: match_.f_match(),
                                    intersect_bp: match_.intersect_bp(),
                                });
                            }
                        } else {
                            eprintln!("Error gathering matches: {:?}", matches.err());
                        }
                    } else {
                        eprintln!(
                            "WARNING: no compatible sketches in path '{}'",
                            query_sig.filename()
                        );
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }

                    if results.is_empty() {
                        None
                    } else {
                        Some(results)
                    }
                }
                Err(err) => {
                    eprintln!("Error loading sketch: {}", err);
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    None
                }
            }
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
            "WARNING: skipped {} query paths - no compatible signatures.",
            skipped_paths
        );
    }
    if failed_paths > 0 {
        eprintln!(
            "WARNING: {} query paths failed to load. See error messages above.",
            failed_paths
        );
    }

    Ok(())
}
