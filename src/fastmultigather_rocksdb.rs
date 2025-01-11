/// fastmultigather_rocksdb: rocksdb-indexed version of fastmultigather.
use anyhow::Result;
use camino::Utf8PathBuf as PathBuf;
use rayon::prelude::*;
use sourmash::index::revindex::{RevIndex, RevIndexOps};
use sourmash::prelude::*;
use sourmash::signature::SigsTrait;
use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::SigStore;

use crate::utils::{
    csvwriter_thread, is_revindex_database, load_collection, BranchwaterGatherResult,
    MultiCollection, ReportType,
};

pub fn fastmultigather_rocksdb(
    queries_file: String,
    index: PathBuf,
    selection: Selection,
    threshold_bp: u32,
    output: Option<String>,
    allow_failed_sigpaths: bool,
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

    let query_collection = load_collection(
        &queries_file,
        &set_selection,
        ReportType::Query,
        allow_failed_sigpaths,
    )?;

    let (n_processed, skipped_paths, failed_paths) =
        fastmultigather_rocksdb_obj(&query_collection, &db, &set_selection, threshold_bp, output)?;

    println!("DONE. Processed {} queries total.", n_processed);

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

pub(crate) fn fastmultigather_rocksdb_obj(
    query_collection: &MultiCollection,
    db: &RevIndex,
    selection: &Selection,
    threshold_bp: u32,
    output: Option<String>,
) -> Result<(usize, usize, usize)> {
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
    let failed_gathers = AtomicUsize::new(0);

    let send = query_collection
        .par_iter()
        .filter_map(|(coll, _idx, record)| {
            let threshold = threshold_bp / selection.scaled().expect("scaled is not set!?");
            let ksize = selection.ksize().expect("ksize not set!?");

            // query downsampling happens here
            match coll.sig_from_record(record) {
                Ok(query_sig) => {
                    let query_filename = query_sig.filename();
                    let query_name = query_sig.name();
                    let query_md5 = query_sig.md5sum();

                    let mut results = vec![];
                    if let Ok(query_mh) = <SigStore as TryInto<KmerMinHash>>::try_into(query_sig) {
                        let _ = processed_sigs.fetch_add(1, atomic::Ordering::SeqCst);
                        // Gather!
                        let (counter, query_colors, hash_to_color) =
                            db.prepare_gather_counters(&query_mh);

                        let matches = db.gather(
                            counter,
                            query_colors,
                            hash_to_color,
                            threshold as usize,
                            &query_mh,
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
                                    query_filename: query_filename.clone(),
                                    query_name: query_name.clone(),
                                    query_md5: query_md5.clone(),
                                    query_bp: query_mh.n_unique_kmers(),
                                    ksize: ksize as u16,
                                    moltype: query_mh.hash_function().to_string(),
                                    scaled: query_mh.scaled(),
                                    query_n_hashes: query_mh.size() as u64,
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
                            let _ = failed_gathers.fetch_add(1, atomic::Ordering::SeqCst);
                        }
                    } else {
                        eprintln!(
                            "WARNING: no compatible sketches in path '{}'",
                            query_filename
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
    send.expect("Unable to send internal data");
    thrd.join().expect("Unable to join CSV writing thread.");

    // done!
    let n_processed: usize = processed_sigs.fetch_max(0, atomic::Ordering::SeqCst);
    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);
    let failed_gathers = failed_gathers.load(atomic::Ordering::SeqCst);

    if n_processed == 0 {
        return Err(anyhow::anyhow!("no search sigs found!?"));
    }

    if failed_gathers > 0 {
        return Err(anyhow::anyhow!(
            "{} failed gathers. See error messages above.",
            failed_gathers
        ));
    }

    Ok((n_processed, skipped_paths, failed_paths))
}
