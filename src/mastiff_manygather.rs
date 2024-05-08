/// mastiff_manygather: mastiff-indexed version of fastmultigather.
use anyhow::Result;
use camino::Utf8PathBuf as PathBuf;
use rayon::prelude::*;
use sourmash::index::revindex::{RevIndex, RevIndexOps};
use sourmash::prelude::*;
use sourmash::signature::SigsTrait;
use std::sync::atomic;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::sync::Mutex;

use crate::utils::{
    csvwriter_thread, is_revindex_database, load_collection, BranchwaterGatherResult, ReportType,
    ThreadManager,
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
    let db = RevIndex::open(index, true, None)?;
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

    // set up manager to allow for ctrl-c handling
    let manager = ThreadManager::new(send, thrd);

    // Wrap ThreadManager in Arc<Mutex> for safe sharing across threads
    let manager_shared = Arc::new(Mutex::new(manager));

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_sigs = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    query_collection
        .par_iter()
        .filter_map(|(_idx, record)| {
            let threshold = threshold_bp / selection.scaled()? as usize;
            let ksize = selection.ksize()?;
            if manager_shared
                .lock()
                .unwrap()
                .interrupted
                .load(Ordering::SeqCst)
            {
                println!("Ctrl-C received, signaling shutdown...");
                return None; // Early exit if interrupted
            }

            // query downsampling happens here
            match query_collection.sig_from_record(record) {
                Ok(query_sig) => {
                    let mut results = vec![];
                    if let Some(query_mh) = query_sig.minhash() {
                        let _ = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
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
                        if let Ok(matches) = matches {
                            for match_ in &matches {
                                results.push(BranchwaterGatherResult {
                                    intersect_bp: match_.intersect_bp(),
                                    f_orig_query: match_.f_orig_query(),
                                    f_match: match_.f_match(),
                                    f_unique_to_query: match_.f_unique_to_query(),
                                    f_unique_weighted: match_.f_unique_weighted(),
                                    average_abund: match_.average_abund(),
                                    median_abund: match_.median_abund(),
                                    std_abund: match_.std_abund(),
                                    match_filename: match_.filename().clone(),
                                    match_name: match_.name().clone(),
                                    match_md5: match_.md5().clone(),
                                    f_match_orig: match_.f_match_orig(),
                                    unique_intersect_bp: match_.unique_intersect_bp(),
                                    gather_result_rank: match_.gather_result_rank(),
                                    remaining_bp: match_.remaining_bp(),
                                    query_filename: query_sig.filename(),
                                    query_name: query_sig.name().clone(),
                                    query_md5: query_sig.md5sum().clone(),
                                    query_bp: query_mh.n_unique_kmers() as usize,
                                    ksize: ksize as usize,
                                    moltype: query_mh.hash_function().to_string(),
                                    scaled: query_mh.scaled() as usize,
                                    query_n_hashes: query_mh.size(),
                                    query_abundance: query_mh.track_abundance(),
                                    query_containment_ani: match_.query_containment_ani(),
                                    match_containment_ani: match_.match_containment_ani(),
                                    average_containment_ani: match_.average_containment_ani(),
                                    max_containment_ani: match_.max_containment_ani(),
                                    n_unique_weighted_found: match_.n_unique_weighted_found(),
                                    sum_weighted_found: match_.sum_weighted_found(),
                                    total_weighted_hashes: match_.total_weighted_hashes(),

                                    query_containment_ani_ci_low: match_
                                        .query_containment_ani_ci_low(),
                                    query_containment_ani_ci_high: match_
                                        .query_containment_ani_ci_high(),
                                    match_containment_ani_ci_low: match_
                                        .match_containment_ani_ci_low(),
                                    match_containment_ani_ci_high: match_
                                        .match_containment_ani_ci_high(),
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
        .try_for_each_with(manager_shared.clone(), |manager, result| {
            manager.lock().unwrap().send(result)
        })?;

    // do some cleanup and error handling -
    manager_shared.lock().unwrap().perform_cleanup();

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
