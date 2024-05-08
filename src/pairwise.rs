/// pairwise: massively parallel in-memory pairwise comparisons.
use anyhow::Result;
use rayon::prelude::*;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;
use std::sync::Arc;
use std::sync::Mutex;

use crate::utils::{
    csvwriter_thread, load_collection, load_sketches, MultiSearchResult, ReportType, ThreadManager,
};
use sourmash::ani_utils::ani_from_containment;
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;

/// Perform pairwise comparisons of all signatures in a list.
///
/// Note: this function loads all _signatures_ into memory.

pub fn pairwise(
    siglist: String,
    threshold: f64,
    selection: &Selection,
    allow_failed_sigpaths: bool,
    estimate_ani: bool,
    write_all: bool,
    output: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Load all sigs into memory at once.
    let collection = load_collection(
        &siglist,
        selection,
        ReportType::General,
        allow_failed_sigpaths,
    )?;

    if collection.len() <= 1 {
        bail!(
            "Pairwise requires two or more sketches. Check input: '{:?}'",
            &siglist
        )
    }
    let sketches = load_sketches(collection, selection, ReportType::General).unwrap();

    // set up a multi-producer, single-consumer channel.
    let (send, recv) =
        std::sync::mpsc::sync_channel::<MultiSearchResult>(rayon::current_num_threads());

    // // & spawn a thread that is dedicated to printing to a buffered output
    let thrd = csvwriter_thread(recv, output);

    // set up manager to allow for ctrl-c handling
    let manager = ThreadManager::new(send, thrd);

    // Wrap ThreadManager in Arc<Mutex> for safe sharing across threads
    let manager_shared = Arc::new(Mutex::new(manager));

    //
    // Main loop: iterate (in parallel) over all signature,
    // Results written to the writer thread above.

    let processed_cmp = AtomicUsize::new(0);
    let ksize = selection.ksize().unwrap() as f64;

    sketches.par_iter().enumerate().for_each(|(idx, query)| {
        // Clone the Arc to get a new reference for this thread
        let manager_clone = manager_shared.clone();
        let mut has_written_comparison = false;
        for against in sketches.iter().skip(idx + 1) {
            if manager_shared
                .lock()
                .unwrap()
                .interrupted
                .load(atomic::Ordering::SeqCst)
            {
                println!("Ctrl-C received, signaling shutdown...");
                return; // Early return to stop processing further
            }
            let overlap = query.minhash.count_common(&against.minhash, false).unwrap() as f64;
            let query1_size = query.minhash.size() as f64;
            let query2_size = against.minhash.size() as f64;

            let containment_q1_in_q2 = overlap / query1_size;
            let containment_q2_in_q1 = overlap / query2_size;

            if containment_q1_in_q2 > threshold || containment_q2_in_q1 > threshold {
                has_written_comparison = true;
                let max_containment = containment_q1_in_q2.max(containment_q2_in_q1);
                let jaccard = overlap / (query1_size + query2_size - overlap);
                let mut query_containment_ani = None;
                let mut match_containment_ani = None;
                let mut average_containment_ani = None;
                let mut max_containment_ani = None;

                // estimate ANI values
                if estimate_ani {
                    let qani = ani_from_containment(containment_q1_in_q2, ksize);
                    let mani = ani_from_containment(containment_q2_in_q1, ksize);
                    query_containment_ani = Some(qani);
                    match_containment_ani = Some(mani);
                    average_containment_ani = Some((qani + mani) / 2.);
                    max_containment_ani = Some(f64::max(qani, mani));
                }
                manager_clone
                    .lock()
                    .unwrap()
                    .send(MultiSearchResult {
                        query_name: query.name.clone(),
                        query_md5: query.md5sum.clone(),
                        match_name: against.name.clone(),
                        match_md5: against.md5sum.clone(),
                        containment: containment_q1_in_q2,
                        max_containment,
                        jaccard,
                        intersect_hashes: overlap,
                        query_containment_ani,
                        match_containment_ani,
                        average_containment_ani,
                        max_containment_ani,
                    })
                    .unwrap()
            }

            let i = processed_cmp.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 100000 == 0 && i > 0 {
                eprintln!("Processed {} comparisons", i);
            }
        }
        if write_all & !has_written_comparison {
            let mut query_containment_ani = None;
            let mut match_containment_ani = None;
            let mut average_containment_ani = None;
            let mut max_containment_ani = None;

            if estimate_ani {
                query_containment_ani = Some(1.0);
                match_containment_ani = Some(1.0);
                average_containment_ani = Some(1.0);
                max_containment_ani = Some(1.0);
            }

            manager_clone
                .lock()
                .unwrap()
                .send(MultiSearchResult {
                    query_name: query.name.clone(),
                    query_md5: query.md5sum.clone(),
                    match_name: query.name.clone(),
                    match_md5: query.md5sum.clone(),
                    containment: 1.0,
                    max_containment: 1.0,
                    jaccard: 1.0,
                    intersect_hashes: query.minhash.size() as f64,
                    query_containment_ani,
                    match_containment_ani,
                    average_containment_ani,
                    max_containment_ani,
                })
                .unwrap()
        }
    });

    // do some cleanup and error handling -
    manager_shared.lock().unwrap().perform_cleanup();

    // done!
    let i: usize = processed_cmp.load(atomic::Ordering::SeqCst);
    eprintln!("DONE. Processed {} comparisons", i);

    Ok(())
}
