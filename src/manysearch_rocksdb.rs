/// manysearch_rocksdb: rocksdb-indexed version of manysearch.
use anyhow::Result;
use camino::Utf8PathBuf as PathBuf;
use rayon::prelude::*;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use sourmash::ani_utils::ani_from_containment;
use sourmash::index::revindex::{RevIndex, RevIndexOps};
use sourmash::selection::Selection;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::SigStore;

use crate::utils::{
    csvwriter_thread, is_revindex_database, load_collection, ManySearchResult, MultiCollection,
    ReportType,
};

pub fn manysearch_rocksdb(
    queries_path: String,
    index: PathBuf,
    selection: Selection,
    minimum_containment: f64,
    output: Option<String>,
    allow_failed_sigpaths: bool,
    output_all_comparisons: bool,
) -> Result<()> {
    if !is_revindex_database(&index) {
        bail!("'{}' is not a valid RevIndex database", index);
    }

    // Open database once
    let db = match RevIndex::open(index, true, None) {
        Ok(db) => db,
        Err(e) => {
            return Err(anyhow::anyhow!(
                "cannot open RocksDB database. Error is: {}",
                e
            ))
        }
    };

    println!("Loaded DB");

    // grab scaled from the database.
    let (_, max_db_scaled) = db
        .collection()
        .min_max_scaled()
        .expect("no records in db?!");

    let selection_scaled: u32 = match selection.scaled() {
        Some(scaled) => {
            if *max_db_scaled > scaled {
                return Err(anyhow::anyhow!(
                    "Error: database scaled is higher than requested scaled"
                ));
            }
            scaled
        }
        None => {
            eprintln!("Setting scaled={} from the database", *max_db_scaled);
            *max_db_scaled
        }
    };

    let mut set_selection = selection;
    set_selection.set_scaled(selection_scaled);

    // Load query paths
    let query_collection = load_collection(
        &queries_path,
        &set_selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;

    let (n_processed, skipped_paths, failed_paths) = manysearch_rocksdb_obj(
        &query_collection,
        &db,
        minimum_containment,
        output,
        output_all_comparisons,
    )?;

    // done!
    eprintln!("DONE. Processed {} search sigs", n_processed);

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

pub(crate) fn manysearch_rocksdb_obj(
    query_collection: &MultiCollection,
    db: &RevIndex,
    minimum_containment: f64,
    output: Option<String>,
    output_all_comparisons: bool,
) -> Result<(usize, usize, usize)> {
    // set up a multi-producer, single-consumer channel.
    let (send, recv) =
        std::sync::mpsc::sync_channel::<ManySearchResult>(rayon::current_num_threads());

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

    let send_result = query_collection
        .par_iter()
        .filter_map(|(coll, _idx, record)| {
            let i = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
            if i % 1000 == 0 && i > 0 {
                eprintln!("Processed {} search sigs", i);
            }

            let mut results = vec![];
            match coll.sig_from_record(record) {
                Ok(query_sig) => {
                    let query_name = query_sig.name().clone();
                    let query_md5 = query_sig.md5sum().clone();
                    let query_file = query_sig.filename().clone();

                    if let Ok(query_mh) = <SigStore as TryInto<KmerMinHash>>::try_into(query_sig) {
                        let query_size = query_mh.size();
                        let counter = db.counter_for_query(&query_mh);
                        let matches =
                            db.matches_from_counter(counter, minimum_containment as usize);

                        // filter the matches for containment
                        for (path, overlap) in matches {
                            let containment = overlap as f64 / query_size as f64;
                            if containment >= minimum_containment || output_all_comparisons {
                                let query_containment_ani = Some(ani_from_containment(
                                    containment,
                                    query_mh.ksize() as f64,
                                ));

                                results.push(ManySearchResult {
                                    query_name: query_name.clone(),
                                    query_md5: query_md5.clone(),
                                    match_name: path.clone(),
                                    containment,
                                    intersect_hashes: overlap as u64,
                                    ksize: query_mh.ksize() as u16,
                                    scaled: query_mh.scaled(),
                                    moltype: query_mh.hash_function().to_string(),
                                    match_md5: None,
                                    jaccard: None,
                                    max_containment: None,
                                    // can't calculate from here -- need to get these from w/in sourmash
                                    average_abund: None,
                                    median_abund: None,
                                    std_abund: None,
                                    query_containment_ani,
                                    match_containment_ani: None,
                                    average_containment_ani: None,
                                    max_containment_ani: None,
                                    n_weighted_found: None,
                                    total_weighted_hashes: None,
                                });
                            }
                        }
                    } else {
                        eprintln!("WARNING: no compatible sketches in path '{}'", query_file);
                        let _ = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                    if results.is_empty() {
                        None
                    } else {
                        Some(results)
                    }
                }
                Err(err) => {
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    eprintln!("Sketch loading error: {}", err);
                    eprintln!(
                        "WARNING: could not load sketches from path '{}'",
                        record.internal_location()
                    );
                    None
                }
            }
        })
        .flatten()
        .try_for_each_with(send, |s, results| {
            if let Err(e) = s.send(results) {
                Err(format!("Unable to send internal data: {:?}", e))
            } else {
                Ok(())
            }
        });

    send_result.expect("Error during parallel processing");
    thrd.join().expect("Unable to join internal thread.");

    let i = processed_sigs.load(atomic::Ordering::SeqCst);

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    Ok((i, skipped_paths, failed_paths))
}
