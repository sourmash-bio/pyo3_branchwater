//! Rust-to-Python interface code for sourmash_plugin_branchwater, using pyo3.
//!
//! If you're using Rust, you're probably most interested in
//! [utils](utils/index.html)

use pyo3::prelude::*;

#[macro_use]
extern crate simple_error;

mod utils;
use crate::utils::build_selection;
use crate::utils::is_revindex_database;
mod check;
mod cluster;
mod fastagather;
mod fastgather;
mod fastmultigather;
mod fastmultigather_rocksdb;
mod index;
mod manysearch;
mod manysearch_rocksdb;
mod manysketch;
mod multisearch;
mod pairwise;
mod search_significance;
mod singlesketch;

use camino::Utf8PathBuf as PathBuf;

#[pyfunction]
#[pyo3(signature = (querylist_path, siglist_path, threshold, ksize, scaled, moltype, output_path=None, ignore_abundance=false, output_all_comparisons=false))]
#[allow(clippy::too_many_arguments)]
fn do_manysearch(
    querylist_path: String,
    siglist_path: String,
    threshold: f64,
    ksize: u8,
    scaled: Option<u32>,
    moltype: String,
    output_path: Option<String>,
    ignore_abundance: Option<bool>,
    output_all_comparisons: Option<bool>,
) -> anyhow::Result<u8> {
    let againstfile_path: PathBuf = siglist_path.clone().into();
    let selection = build_selection(ksize, scaled, &moltype);
    eprintln!("selection scaled: {:?}", selection.scaled());
    let allow_failed_sigpaths = true;

    let ignore_abundance = ignore_abundance.unwrap_or(false);
    let output_all_comparisons = output_all_comparisons.unwrap_or(false);

    // if siglist_path is revindex, run rocksdb manysearch; otherwise run manysearch
    if is_revindex_database(&againstfile_path) {
        // note: manysearch_rocksdb ignores abundance automatically.
        match manysearch_rocksdb::manysearch_rocksdb(
            querylist_path,
            againstfile_path,
            selection,
            threshold,
            output_path,
            allow_failed_sigpaths,
            output_all_comparisons,
        ) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    } else {
        match manysearch::manysearch(
            querylist_path,
            siglist_path,
            selection,
            threshold,
            output_path,
            allow_failed_sigpaths,
            ignore_abundance,
            output_all_comparisons,
        ) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    }
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (query_filename, siglist_path, threshold_bp, ksize, scaled, moltype, output_path_prefetch=None, output_path_gather=None))]
fn do_fastgather(
    query_filename: String,
    siglist_path: String,
    threshold_bp: u64,
    ksize: u8,
    scaled: Option<u32>,
    moltype: String,
    output_path_prefetch: Option<String>,
    output_path_gather: Option<String>,
) -> anyhow::Result<u8> {
    let selection = build_selection(ksize, scaled, &moltype);
    let allow_failed_sigpaths = true;

    match fastgather::fastgather(
        query_filename,
        siglist_path,
        threshold_bp,
        selection,
        output_path_prefetch,
        output_path_gather,
        allow_failed_sigpaths,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (query_filenames, siglist_path, threshold_bp, ksize, scaled, moltype, output_path=None, save_matches=false, create_empty_results=false))]
fn do_fastmultigather(
    query_filenames: String,
    siglist_path: String,
    threshold_bp: u64,
    ksize: u8,
    scaled: Option<u32>,
    moltype: String,
    output_path: Option<String>,
    save_matches: bool,
    create_empty_results: bool,
) -> anyhow::Result<u8> {
    let againstfile_path: camino::Utf8PathBuf = siglist_path.clone().into();
    let selection = build_selection(ksize, scaled, &moltype);
    let allow_failed_sigpaths = true;

    // if a siglist path is a revindex, run rocksdb fastmultigather. If not, run multigather
    if is_revindex_database(&againstfile_path) {
        match fastmultigather_rocksdb::fastmultigather_rocksdb(
            query_filenames,
            againstfile_path,
            selection.clone(),
            threshold_bp as u32,
            output_path,
            allow_failed_sigpaths,
        ) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    } else {
        match fastmultigather::fastmultigather(
            query_filenames,
            siglist_path,
            threshold_bp as u32,
            scaled,
            selection,
            allow_failed_sigpaths,
            save_matches,
            output_path,
            create_empty_results,
        ) {
            Ok(_) => Ok(0),
            Err(e) => {
                eprintln!("Error: {e}");
                Ok(1)
            }
        }
    }
}

#[pyfunction]
fn set_global_thread_pool(num_threads: usize) -> PyResult<usize> {
    if std::panic::catch_unwind(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
    })
    .is_ok()
    {
        Ok(rayon::current_num_threads())
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            "Could not set the number of threads. Global thread pool might already be initialized.",
        ))
    }
}

#[pyfunction]
#[pyo3(signature = (siglist, ksize, scaled, moltype, output, colors, use_internal_storage))]
fn do_index(
    siglist: String,
    ksize: u8,
    scaled: Option<u32>,
    moltype: String,
    output: String,
    colors: bool,
    use_internal_storage: bool,
) -> anyhow::Result<u8> {
    let selection = build_selection(ksize, scaled, &moltype);
    let allow_failed_sigpaths = false;
    match index::index(
        siglist,
        selection,
        output,
        colors,
        allow_failed_sigpaths,
        use_internal_storage,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
fn do_check(index: String, quick: bool, rw: bool) -> anyhow::Result<u8> {
    let idx: PathBuf = index.into();
    match check::check(idx, quick, rw) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
#[pyo3(signature = (querylist_path, siglist_path, threshold, ksize, scaled, moltype, estimate_ani, estimate_prob_overlap, output_all_comparisons, output_path=None))]
#[allow(clippy::too_many_arguments)]
fn do_multisearch(
    querylist_path: String,
    siglist_path: String,
    threshold: f64,
    ksize: u8,
    scaled: Option<u32>,
    moltype: String,
    estimate_ani: bool,
    estimate_prob_overlap: bool,
    output_all_comparisons: bool,
    output_path: Option<String>,
) -> anyhow::Result<u8> {
    let _ = env_logger::try_init();

    let selection = build_selection(ksize, scaled, &moltype);
    let allow_failed_sigpaths = true;

    match multisearch::multisearch(
        querylist_path,
        siglist_path,
        threshold,
        selection,
        allow_failed_sigpaths,
        estimate_ani,
        estimate_prob_overlap,
        output_all_comparisons,
        output_path,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (siglist_path, threshold, ksize, scaled, moltype, estimate_ani, write_all, output_all_comparisons, output_path=None))]
fn do_pairwise(
    siglist_path: String,
    threshold: f64,
    ksize: u8,
    scaled: Option<u32>,
    moltype: String,
    estimate_ani: bool,
    write_all: bool,
    output_all_comparisons: bool,
    output_path: Option<String>,
) -> anyhow::Result<u8> {
    let selection = build_selection(ksize, scaled, &moltype);
    let allow_failed_sigpaths = true;
    match pairwise::pairwise(
        siglist_path,
        threshold,
        selection,
        allow_failed_sigpaths,
        estimate_ani,
        write_all,
        output_all_comparisons,
        output_path,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
fn do_manysketch(
    filelist: String,
    param_str: String,
    output: String,
    singleton: bool,
    force: bool,
) -> anyhow::Result<u8> {
    match manysketch::manysketch(filelist, param_str, output, singleton, force) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
#[pyo3(signature = (input_filenames, input_moltype, param_str, output, name))]
fn do_singlesketch(
    input_filenames: Vec<String>,
    input_moltype: String,
    param_str: String,
    output: String,
    name: String,
) -> anyhow::Result<u8> {
    match singlesketch::singlesketch(input_filenames, input_moltype, param_str, output, name) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
#[pyo3(signature = (pairwise_csv, output_clusters, similarity_column, similarity_threshold, cluster_sizes=None))]
fn do_cluster(
    pairwise_csv: String,
    output_clusters: String,
    similarity_column: String,
    similarity_threshold: f64,
    cluster_sizes: Option<String>,
) -> anyhow::Result<u8> {
    match cluster::cluster(
        pairwise_csv,
        output_clusters,
        similarity_column,
        similarity_threshold,
        cluster_sizes,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (query_filename, index, input_moltype, threshold_hashes, ksize, scaled, moltype, output))]
fn do_fastagather(
    query_filename: String,
    index: String,
    input_moltype: String,
    threshold_hashes: u64,
    ksize: u32,
    scaled: u32,
    moltype: String,
    output: String,
) -> anyhow::Result<u8> {
    match fastagather::fastagather(
        query_filename,
        index,
        input_moltype,
        threshold_hashes,
        ksize,
        scaled,
        moltype,
        output,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

/// Module interface for the `sourmash_plugin_branchwater` extension module.

#[pymodule]
fn sourmash_plugin_branchwater(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_manysearch, m)?)?;
    m.add_function(wrap_pyfunction!(do_fastgather, m)?)?;
    m.add_function(wrap_pyfunction!(do_fastmultigather, m)?)?;
    m.add_function(wrap_pyfunction!(do_index, m)?)?;
    m.add_function(wrap_pyfunction!(do_check, m)?)?;
    m.add_function(wrap_pyfunction!(do_manysketch, m)?)?;
    m.add_function(wrap_pyfunction!(set_global_thread_pool, m)?)?;
    m.add_function(wrap_pyfunction!(do_multisearch, m)?)?;
    m.add_function(wrap_pyfunction!(do_pairwise, m)?)?;
    m.add_function(wrap_pyfunction!(do_cluster, m)?)?;
    m.add_function(wrap_pyfunction!(do_singlesketch, m)?)?;
    m.add_function(wrap_pyfunction!(do_fastagather, m)?)?;

    Ok(())
}
