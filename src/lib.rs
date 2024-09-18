/// Python interface Rust code for sourmash_plugin_branchwater.
use pyo3::prelude::*;

#[macro_use]
extern crate simple_error;

mod utils;
use crate::utils::build_selection;
use crate::utils::is_revindex_database;
mod check;
mod cluster;
mod fastgather;
mod fastmultigather;
mod index;
mod manysearch;
mod manysketch;
mod mastiff_manygather;
mod mastiff_manysearch;
mod multisearch;
mod pairwise;
mod sigcat;

use camino::Utf8PathBuf as PathBuf;

#[pyfunction]
#[pyo3(signature = (querylist_path, siglist_path, threshold, ksize, scaled, moltype, output_path=None))]
fn do_manysearch(
    querylist_path: String,
    siglist_path: String,
    threshold: f64,
    ksize: u8,
    scaled: usize,
    moltype: String,
    output_path: Option<String>,
) -> anyhow::Result<u8> {
    let againstfile_path: PathBuf = siglist_path.clone().into();
    let selection = build_selection(Some(ksize), Some(scaled), Some(&moltype));
    eprintln!("selection scaled: {:?}", selection.scaled());
    let allow_failed_sigpaths = true;

    // if siglist_path is revindex, run mastiff_manysearch; otherwise run manysearch
    if is_revindex_database(&againstfile_path) {
        match mastiff_manysearch::mastiff_manysearch(
            querylist_path,
            againstfile_path,
            &selection,
            threshold,
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
        match manysearch::manysearch(
            querylist_path,
            siglist_path,
            &selection,
            threshold,
            output_path,
            allow_failed_sigpaths,
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
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
    moltype: String,
    output_path_prefetch: Option<String>,
    output_path_gather: Option<String>,
) -> anyhow::Result<u8> {
    let selection = build_selection(Some(ksize), Some(scaled), Some(&moltype));
    let allow_failed_sigpaths = true;

    match fastgather::fastgather(
        query_filename,
        siglist_path,
        threshold_bp,
        scaled,
        &selection,
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
#[pyo3(signature = (query_filenames, siglist_path, threshold_bp, ksize, scaled, moltype, output_path=None, save_matches=false, create_empty_results=false))]
fn do_fastmultigather(
    query_filenames: String,
    siglist_path: String,
    threshold_bp: usize,
    ksize: u8,
    scaled: usize,
    moltype: String,
    output_path: Option<String>,
    save_matches: bool,
    create_empty_results: bool,
) -> anyhow::Result<u8> {
    let againstfile_path: camino::Utf8PathBuf = siglist_path.clone().into();
    let selection = build_selection(Some(ksize), Some(scaled), Some(&moltype));
    let allow_failed_sigpaths = true;

    // if a siglist path is a revindex, run mastiff_manygather. If not, run multigather
    if is_revindex_database(&againstfile_path) {
        match mastiff_manygather::mastiff_manygather(
            query_filenames,
            againstfile_path,
            &selection,
            threshold_bp,
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
        if output_path.is_some() {
            bail!("output path specified, but not running fastmultigather against a rocksdb. See issue #239");
        }
        match fastmultigather::fastmultigather(
            query_filenames,
            siglist_path,
            threshold_bp,
            scaled,
            &selection,
            allow_failed_sigpaths,
            save_matches,
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
fn do_index(
    siglist: String,
    ksize: u8,
    scaled: usize,
    moltype: String,
    output: String,
    colors: bool,
    use_internal_storage: bool,
) -> anyhow::Result<u8> {
    let selection = build_selection(Some(ksize), Some(scaled), Some(&moltype));
    let allow_failed_sigpaths = false;
    match index::index(
        siglist,
        &selection,
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
fn do_check(index: String, quick: bool) -> anyhow::Result<u8> {
    let idx: PathBuf = index.into();
    match check::check(idx, quick) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pyfunction]
#[pyo3(signature = (querylist_path, siglist_path, threshold, ksize, scaled, moltype, estimate_ani, output_path=None))]
#[allow(clippy::too_many_arguments)]
fn do_multisearch(
    querylist_path: String,
    siglist_path: String,
    threshold: f64,
    ksize: u8,
    scaled: usize,
    moltype: String,
    estimate_ani: bool,
    output_path: Option<String>,
) -> anyhow::Result<u8> {
    let selection = build_selection(Some(ksize), Some(scaled), Some(&moltype));
    let allow_failed_sigpaths = true;

    match multisearch::multisearch(
        querylist_path,
        siglist_path,
        threshold,
        &selection,
        allow_failed_sigpaths,
        estimate_ani,
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
#[pyo3(signature = (siglist_path, threshold, ksize, scaled, moltype, estimate_ani, write_all, output_path=None))]
fn do_pairwise(
    siglist_path: String,
    threshold: f64,
    ksize: u8,
    scaled: usize,
    moltype: String,
    estimate_ani: bool,
    write_all: bool,
    output_path: Option<String>,
) -> anyhow::Result<u8> {
    let selection = build_selection(Some(ksize), Some(scaled), Some(&moltype));
    let allow_failed_sigpaths = true;
    match pairwise::pairwise(
        siglist_path,
        threshold,
        &selection,
        allow_failed_sigpaths,
        estimate_ani,
        write_all,
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
fn do_sigcat(
    sigfiles: String,
    output_path: String,
    force: bool,
    ksize: Option<u8>,
    scaled: Option<usize>,
    moltype: Option<String>,
) -> anyhow::Result<u8> {
    let selection = build_selection(ksize, scaled, moltype.as_deref());
    let allow_failed_sigpaths = true;
    match sigcat::sig_cat(
        sigfiles,
        &selection,
        allow_failed_sigpaths,
        output_path,
        force,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

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
    m.add_function(wrap_pyfunction!(do_sigcat, m)?)?;
    Ok(())
}
