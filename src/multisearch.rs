/// multisearch: massively parallel in-memory sketch search.
use anyhow::Result;
use rayon::prelude::*;
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use crate::utils::{
    csvwriter_thread, load_collection, load_sketches, MultiSearchResult, ReportType,
};
use sourmash::ani_utils::ani_from_containment;

/// Search many queries against a list of signatures.
///
/// Note: this function loads all _queries_ into memory, and iterates over
/// database once.

pub fn multisearch(
    query_filepath: String,
    against_filepath: String,
    threshold: f64,
    selection: &Selection,
    allow_failed_sigpaths: bool,
    estimate_ani: bool,
    output: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    let allow_empty_collection = false;
    // Load all queries into memory at once.

    let query_collection = load_collection(
        &query_filepath,
        selection,
        ReportType::Query,
        allow_failed_sigpaths,
        allow_empty_collection,
    )?;
    let queries = load_sketches(query_collection, selection, ReportType::Query).unwrap();

    // Load all against sketches into memory at once.
    let against_collection = load_collection(
        &against_filepath,
        selection,
        ReportType::Against,
        allow_failed_sigpaths,
        allow_empty_collection,
    )?;
    let against = load_sketches(against_collection, selection, ReportType::Against).unwrap();

    // set up a multi-producer, single-consumer channel.
    let (send, recv) =
        std::sync::mpsc::sync_channel::<MultiSearchResult>(rayon::current_num_threads());

    // // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = csvwriter_thread(recv, output);

    //
    // Main loop: iterate (in parallel) over all search signature paths,
    // loading them individually and searching them. Stuff results into
    // the writer thread above.
    //

    let processed_cmp = AtomicUsize::new(0);
    let ksize = selection.ksize().unwrap() as f64;

    let send = against
        .par_iter()
        .filter_map(|against| {
            let mut results = vec![];
            // search for matches & save containment.
            for query in queries.iter() {
                let i = processed_cmp.fetch_add(1, atomic::Ordering::SeqCst);
                if i % 100000 == 0 && i > 0 {
                    eprintln!("Processed {} comparisons", i);
                }

                let overlap = query.minhash.count_common(&against.minhash, false).unwrap() as f64;
                // use downsampled sizes
                let query_size = query.minhash.size() as f64;
                let target_size = against.minhash.size() as f64;

                let containment_query_in_target = overlap / query_size;

                if containment_query_in_target > threshold {
                    let containment_target_in_query = overlap / target_size;
                    let max_containment =
                        containment_query_in_target.max(containment_target_in_query);
                    let jaccard = overlap / (target_size + query_size - overlap);
                    let mut query_containment_ani = None;
                    let mut match_containment_ani = None;
                    let mut average_containment_ani = None;
                    let mut max_containment_ani = None;

                    // estimate ANI values
                    if estimate_ani {
                        let qani = ani_from_containment(containment_query_in_target, ksize);
                        let mani = ani_from_containment(containment_target_in_query, ksize);
                        query_containment_ani = Some(qani);
                        match_containment_ani = Some(mani);
                        average_containment_ani = Some((qani + mani) / 2.);
                        max_containment_ani = Some(f64::max(qani, mani));
                    }

                    results.push(MultiSearchResult {
                        query_name: query.name.clone(),
                        query_md5: query.md5sum.clone(),
                        match_name: against.name.clone(),
                        match_md5: against.md5sum.clone(),
                        containment: containment_query_in_target,
                        max_containment,
                        jaccard,
                        intersect_hashes: overlap,
                        query_containment_ani,
                        match_containment_ani,
                        average_containment_ani,
                        max_containment_ani,
                    })
                }
            }
            if results.is_empty() {
                None
            } else {
                Some(results)
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
    let i: usize = processed_cmp.fetch_max(0, atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} comparisons", i);

    Ok(())
}
