/// fastgather: Run gather with a query against a list of files.
use anyhow::Result;
use sourmash::prelude::Select;
use sourmash::selection::Selection;

use crate::utils::{
    consume_query_by_gather, load_collection, load_sketches_above_threshold, write_prefetch,
    ReportType,
};

#[allow(clippy::too_many_arguments)]
pub fn fastgather(
    query_filepath: String,
    against_filepath: String,
    threshold_bp: usize,
    scaled: usize,
    selection: &Selection,
    gather_output: Option<String>,
    prefetch_output: Option<String>,
    allow_failed_sigpaths: bool,
) -> Result<()> {
    let allow_empty_collection = false;
    let query_collection = load_collection(
        &query_filepath,
        selection,
        ReportType::Query,
        allow_failed_sigpaths,
        allow_empty_collection,
    )?;

    if query_collection.len() != 1 {
        bail!(
            "Fastgather requires a single query sketch. Check input: '{:?}'",
            &query_filepath
        )
    }
    // get single query sig and minhash
    let query_sig = query_collection.sig_for_dataset(0)?; // need this for original md5sum
    let query_sig_ds = query_sig.clone().select(selection)?; // downsample
    let query_mh = match query_sig_ds.minhash() {
        Some(query_mh) => query_mh,
        None => {
            bail!("No query sketch matching selection parameters.");
        }
    };
    // load collection to match against.
    let against_collection = load_collection(
        &against_filepath,
        selection,
        ReportType::Against,
        allow_failed_sigpaths,
        allow_empty_collection,
    )?;

    // calculate the minimum number of hashes based on desired threshold
    let threshold_hashes: u64 = {
        let x = threshold_bp / scaled;
        if x > 0 {
            x
        } else {
            1
        }
    }
    .try_into()?;

    eprintln!(
        "using threshold overlap: {} {}",
        threshold_hashes, threshold_bp
    );

    // load a set of sketches, filtering for those with overlaps > threshold
    let result = load_sketches_above_threshold(against_collection, query_mh, threshold_hashes)?;
    let matchlist = result.0;
    let skipped_paths = result.1;
    let failed_paths = result.2;
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

    if matchlist.is_empty() {
        eprintln!("No search signatures loaded, exiting.");
        return Ok(());
    }

    if prefetch_output.is_some() {
        write_prefetch(&query_sig, prefetch_output, &matchlist).ok();
    }

    // run the gather!
    consume_query_by_gather(
        query_sig,
        scaled as u64,
        matchlist,
        threshold_hashes,
        gather_output,
    )
    .ok();
    Ok(())
}
