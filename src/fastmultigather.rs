/// fastmultigather: Run gather for multiple queries against a list of files.
use anyhow::Result;
use rayon::prelude::*;

use sourmash::prelude::ToWriter;
use sourmash::{selection::Selection, signature::SigsTrait};

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use camino::Utf8Path as PathBuf;

use std::collections::HashSet;
use std::fs::File;

use sourmash::signature::Signature;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::sketch::Sketch;

use crate::utils::{
    consume_query_by_gather, load_collection, load_sketches, write_prefetch, PrefetchResult,
    ReportType,
};

pub fn fastmultigather(
    query_filepath: String,
    against_filepath: String,
    threshold_bp: usize,
    scaled: usize,
    selection: &Selection,
    allow_failed_sigpaths: bool,
    save_matches: bool,
    create_empty_results: bool,
) -> Result<()> {
    let allow_empty_collection = false;
    // load query collection
    let query_collection = load_collection(
        &query_filepath,
        selection,
        ReportType::Query,
        allow_failed_sigpaths,
        allow_empty_collection,
    )?;

    let threshold_hashes: u64 = {
        let x = threshold_bp / scaled;
        if x > 0 {
            x
        } else {
            1
        }
    }
    .try_into()
    .unwrap();

    println!("threshold overlap: {} {}", threshold_hashes, threshold_bp);

    // load against collection
    let against_collection = load_collection(
        &against_filepath,
        selection,
        ReportType::Against,
        allow_failed_sigpaths,
        allow_empty_collection,
    )?;
    // load against sketches into memory, downsampling on the way
    let against = load_sketches(against_collection, selection, ReportType::Against).unwrap();

    // Iterate over all queries => do prefetch and gather!
    let processed_queries = AtomicUsize::new(0);
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    query_collection.par_iter().for_each(|(_idx, record)| {
        // increment counter of # of queries. q: could we instead use the _idx from par_iter(), or will it vary based on thread?
        let _i = processed_queries.fetch_add(1, atomic::Ordering::SeqCst);
        // Load query sig (downsampling happens here)
        match query_collection.sig_from_record(record) {
            Ok(query_sig) => {
                let name = query_sig.name();
                let prefix = name.split(' ').next().unwrap_or_default().to_string();
                let location = PathBuf::new(&prefix).file_name().unwrap();
                if let Some(query_mh) = query_sig.minhash() {
                    let mut matching_hashes = if save_matches { Some(Vec::new()) } else { None };
                    let matchlist: BinaryHeap<PrefetchResult> = against
                        .iter()
                        .filter_map(|against| {
                            let mut mm: Option<PrefetchResult> = None;
                            if let Ok(overlap) = against.minhash.count_common(query_mh, false) {
                                if overlap >= threshold_hashes {
                                    if save_matches {
                                        if let Ok(intersection) =
                                            against.minhash.intersection(query_mh)
                                        {
                                            matching_hashes
                                                .as_mut()
                                                .unwrap()
                                                .extend(intersection.0);
                                        }
                                    }
                                    let result = PrefetchResult {
                                        name: against.name.clone(),
                                        md5sum: against.md5sum.clone(),
                                        minhash: against.minhash.clone(),
                                        location: against.location.clone(),
                                        overlap,
                                    };
                                    mm = Some(result);
                                }
                            }
                            mm
                        })
                        .collect();
                    if !matchlist.is_empty() {
                        let prefetch_output = format!("{}.prefetch.csv", location);
                        let gather_output = format!("{}.gather.csv", location);

                        // Save initial list of matches to prefetch output
                        write_prefetch(&query_sig, Some(prefetch_output), &matchlist).ok();

                        // Now, do the gather!
                        consume_query_by_gather(
                            query_sig.clone(),
                            scaled as u64,
                            matchlist,
                            threshold_hashes,
                            Some(gather_output),
                        )
                        .ok();

                        // Save matching hashes to .sig file if save_matches is true
                        if save_matches {
                            if let Some(hashes) = matching_hashes {
                                let sig_filename = format!("{}.matches.sig", name);
                                if let Ok(mut file) = File::create(&sig_filename) {
                                    let unique_hashes: HashSet<u64> = hashes.into_iter().collect();
                                    let mut new_mh = KmerMinHash::new(
                                        query_mh.scaled().try_into().unwrap(),
                                        query_mh.ksize().try_into().unwrap(),
                                        query_mh.hash_function().clone(),
                                        query_mh.seed(),
                                        false,
                                        query_mh.num(),
                                    );
                                    new_mh
                                        .add_many(&unique_hashes.into_iter().collect::<Vec<_>>())
                                        .ok();
                                    let mut signature = Signature::default();
                                    signature.push(Sketch::MinHash(new_mh));
                                    signature.set_filename(&name);
                                    if let Err(e) = signature.to_writer(&mut file) {
                                        eprintln!("Error writing signature file: {}", e);
                                    }
                                } else {
                                    eprintln!("Error creating signature file: {}", sig_filename);
                                }
                            }
                        }
                    } else {
                        println!("No matches to '{}'", location);
                        if create_empty_results {
                            let prefetch_output = format!("{}.prefetch.csv", location);
                            let gather_output = format!("{}.gather.csv", location);
                            // touch output files
                            match std::fs::File::create(&prefetch_output) {
                                Ok(_) => {}
                                Err(e) => {
                                    eprintln!("Failed to create empty prefetch output: {}", e)
                                }
                            }
                            match std::fs::File::create(&gather_output) {
                                Ok(_) => {}
                                Err(e) => eprintln!("Failed to create empty gather output: {}", e),
                            }
                        }
                    }
                } else {
                    // different warning here? Could not load sig from record??
                    eprintln!(
                        "WARNING: no compatible sketches in path '{}'",
                        record.internal_location()
                    );
                    let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                }
            }
            Err(_) => {
                // different warning here? Could not load sig from record??
                eprintln!(
                    "WARNING: no compatible sketches in path '{}'",
                    record.internal_location()
                );
                let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
            }
        }
    });

    println!(
        "DONE. Processed {} queries total.",
        processed_queries.into_inner()
    );

    let skipped_paths = skipped_paths.into_inner();
    let failed_paths = failed_paths.into_inner();

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
