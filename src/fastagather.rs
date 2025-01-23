use anyhow::{bail, Result};
use needletail::parse_fastx_file;
use sourmash::selection::Selection;
// use sourmash::encodings::HashFunctions;
use camino::Utf8PathBuf;
use sourmash::index::revindex::RevIndex;
use sourmash::index::revindex::RevIndexOps;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::{KmerMinHash, KmerMinHashBTree};
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::utils::{
    build_selection, csvwriter_thread, is_revindex_database, BranchwaterGatherResult,
};

// to start, implement straightforward record --> sketch --> gather
// other ideas/to do:
// - add full-file (lower resolution) prefetch first, to reduce search space
// - parallelize and/or batch records?
// - write function to filter fasta entries for those with matches (or those without)
// - could use that with this structure for charcoal decontam or other functions
// - add rocksdb search -- only way this will make sense.
// to do -- use input moltype to check that we can build desired moltype
// let _input_moltype = input_moltype.to_ascii_lowercase();

#[allow(clippy::too_many_arguments)]
pub fn fastagather(
    query_filename: String,
    index: String,
    _input_moltype: String,
    threshold_hashes: u64,
    ksize: u32,
    scaled: u32,
    moltype: String,
    output: String,
) -> Result<()> {
    let index_path = Utf8PathBuf::from(index.clone());

    if !is_revindex_database(&index_path) {
        bail!("'{}' is not a valid RevIndex database", index_path);
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

    let selection = build_selection(ksize as u8, Some(scaled), &moltype);

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

    let (n_processed, n_failed_records, n_matched_gathers) = fastmultigather_fasta_rocksdb_obj(
        query_filename.as_str(),
        &db,
        &set_selection,
        threshold_hashes,
        output,
    )?;

    println!("DONE. Processed {} records total.", n_processed);

    if n_failed_records > 0 {
        eprintln!(
            "WARNING: {} records could not be queried. See error messages above.",
            n_failed_records
        );
    }

    eprintln!(
        "Found gather results for {}/{} records.",
        n_matched_gathers, n_processed,
    );

    Ok(())
}

pub(crate) fn fastmultigather_fasta_rocksdb_obj(
    query_filename: &str,
    db: &RevIndex,
    selection: &Selection,
    threshold_hashes: u64,
    output: String,
) -> Result<(usize, usize, usize)> {
    // Set up a multi-producer, single-consumer channel.
    let (send, recv) =
        std::sync::mpsc::sync_channel::<BranchwaterGatherResult>(rayon::current_num_threads());

    // Spawn a thread for writing CSV output.
    let thrd = csvwriter_thread(recv, Some(output));

    // Atomic counters for tracking progress and failures.
    let processed_records = AtomicUsize::new(0);
    let failed_records = AtomicUsize::new(0);
    let matched_gathers = AtomicUsize::new(0);

    // Extract and validate the required fields from `Selection` using `expect`.
    let ksize = selection.ksize().expect("ksize is not set in selection");
    let scaled = selection.scaled().expect("scaled is not set in selection");
    let hash_function = selection
        .moltype()
        .expect("moltype is not set in selection");

    // Build the minhash template.
    let mh_template = KmerMinHashBTree::new(
        scaled,
        ksize,
        hash_function,
        42,                                // Example seed value
        selection.abund().unwrap_or(true), // Default to tracking abundance if not specified
        0,                                 // Default value for num, can be changed as needed
    );

    // Open the FASTA file reader.
    let mut reader = parse_fastx_file(query_filename)
        .map_err(|err| anyhow::anyhow!("Error opening file {}: {:?}", query_filename, err))?;

    // Main loop to process each record in the FASTA file.
    while let Some(record_result) = reader.next() {
        // Clone the template minhash for this record.
        let mut build_mh = mh_template.clone();
        match record_result {
            Ok(record) => {
                let query_name = std::str::from_utf8(record.id())
                    .expect("record.id() contains invalid UTF-8")
                    .to_string();
                if let Err(err) = build_mh.add_sequence(&record.seq(), true) {
                    eprintln!(
                        "Error building minhash for record '{}' in {}: {:?}",
                        query_name, query_filename, err
                    );
                    failed_records.fetch_add(1, Ordering::SeqCst);
                    continue;
                }
                let query_mh = KmerMinHash::from(build_mh);
                let query_md5 = query_mh.md5sum();
                let n_unique_kmers = query_mh.n_unique_kmers();
                processed_records.fetch_add(1, Ordering::SeqCst);

                // Prepare gather counters.
                let (counter, query_colors, hash_to_color) = db.prepare_gather_counters(&query_mh);

                // Perform the gather operation.
                match db.gather(
                    counter,
                    query_colors,
                    hash_to_color,
                    threshold_hashes as usize,
                    &query_mh,
                    Some(selection.clone()),
                ) {
                    Ok(matches) => {
                        let mut results = vec![];
                        for match_ in matches {
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
                                query_filename: query_filename.to_string(),
                                query_name: query_name.clone(),
                                query_md5: query_md5.clone(),
                                query_bp: n_unique_kmers,
                                ksize: selection.ksize().expect("ksize not set!?") as u16,
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
                                query_containment_ani_ci_low: match_.query_containment_ani_ci_low(),
                                query_containment_ani_ci_high: match_
                                    .query_containment_ani_ci_high(),
                                match_containment_ani_ci_low: match_.match_containment_ani_ci_low(),
                                match_containment_ani_ci_high: match_
                                    .match_containment_ani_ci_high(),
                            });
                        }

                        if !results.is_empty() {
                            eprintln!("Matches found for record: {}", query_name);
                            matched_gathers.fetch_add(1, Ordering::SeqCst);
                            for result in results {
                                if let Err(err) = send.send(result) {
                                    eprintln!("Error sending result: {:?}", err);
                                }
                            }
                        }
                    }
                    Err(err) => {
                        eprintln!(
                            "Error gathering matches for query '{}': {:?}",
                            query_name, err
                        );
                    }
                }
            }
            Err(err) => {
                eprintln!("Error processing record in {}: {:?}", query_filename, err);
                failed_records.fetch_add(1, Ordering::SeqCst);
            }
        }
    }

    // Cleanup and join the writer thread.
    drop(send); // Close the sender to signal the thread to finish.
    thrd.join().expect("Unable to join CSV writer thread.");

    // Gather final stats.
    let n_processed = processed_records.load(Ordering::SeqCst);
    let n_failed_records = failed_records.load(Ordering::SeqCst);
    let n_matched_gathers = matched_gathers.load(Ordering::SeqCst);

    if n_processed == 0 {
        Err(anyhow::anyhow!(
            "No records were processed from the FASTA file."
        ))
    } else {
        Ok((n_processed, n_failed_records, n_matched_gathers))
    }
}
