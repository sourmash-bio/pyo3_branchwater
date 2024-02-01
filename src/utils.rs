/// Utility functions for sourmash_plugin_branchwater.
use rayon::prelude::*;
use sourmash::encodings::HashFunctions;
use sourmash::manifest::Manifest;
use sourmash::selection::Select;

use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::panic;
use std::path::{Path, PathBuf};

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use anyhow::{anyhow, Result};
use std::cmp::{Ordering, PartialOrd};

use sourmash::collection::Collection;
use sourmash::manifest::Record;
use sourmash::selection::Selection;
use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::sketch::Sketch;
use sourmash::storage::{FSStorage, InnerStorage, SigStore};

/// Structure to hold overlap information from comparisons.

pub struct PrefetchResult {
    pub name: String,
    pub md5sum: String,
    pub minhash: KmerMinHash,
    pub overlap: u64,
}

impl Ord for PrefetchResult {
    fn cmp(&self, other: &PrefetchResult) -> Ordering {
        self.overlap.cmp(&other.overlap)
    }
}

impl PartialOrd for PrefetchResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for PrefetchResult {
    fn eq(&self, other: &Self) -> bool {
        self.overlap == other.overlap
    }
}

impl Eq for PrefetchResult {}

/// Find sketches in 'sketchlist' that overlap with 'query' above
/// specified threshold.

pub fn prefetch(
    query_mh: &KmerMinHash,
    sketchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
) -> BinaryHeap<PrefetchResult> {
    sketchlist
        .into_par_iter()
        .filter_map(|result| {
            let mut mm = None;
            let searchsig = &result.minhash;
            // TODO: fix Select so we can go back to downsample: false here
            let overlap = searchsig.count_common(query_mh, true);
            if let Ok(overlap) = overlap {
                if overlap >= threshold_hashes {
                    let result = PrefetchResult { overlap, ..result };
                    mm = Some(result);
                }
            }
            mm
        })
        .collect()
}

/// Write list of prefetch matches.
pub fn write_prefetch(
    query: &SigStore,
    prefetch_output: Option<String>,
    matchlist: &BinaryHeap<PrefetchResult>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Define the writer to stdout by default
    let mut writer: Box<dyn Write> = Box::new(std::io::stdout());

    if let Some(output_path) = &prefetch_output {
        // Account for potential missing dir in output path
        let directory_path = Path::new(output_path).parent();

        // If a directory path exists in the filename, create it if it doesn't already exist
        if let Some(dir) = directory_path {
            create_dir_all(dir)?;
        }

        let file = File::create(output_path)?;
        writer = Box::new(BufWriter::new(file));
    }

    writeln!(
        &mut writer,
        "query_filename,query_name,query_md5,match_name,match_md5,intersect_bp"
    )
    .ok();

    for m in matchlist.iter() {
        writeln!(
            &mut writer,
            "{},\"{}\",{},\"{}\",{},{}",
            query.filename(),
            query.name(),
            query.md5sum(),
            m.name,
            m.md5sum,
            m.overlap
        )
        .ok();
    }

    Ok(())
}

/// Load a list of filenames from a file. Exits on bad lines.
pub fn load_sketchlist_filenames<P: AsRef<Path>>(sketchlist_filename: &P) -> Result<Vec<PathBuf>> {
    let sketchlist_file = BufReader::new(File::open(sketchlist_filename)?);

    let mut sketchlist_filenames: Vec<PathBuf> = Vec::new();
    for line in sketchlist_file.lines() {
        let line = match line {
            Ok(v) => v,
            Err(_) => {
                return {
                    let filename = sketchlist_filename.as_ref().display();
                    let msg = format!("invalid line in fromfile '{}'", filename);
                    Err(anyhow!(msg))
                }
            }
        };

        if !line.is_empty() {
            let mut path = PathBuf::new();
            path.push(line);
            sketchlist_filenames.push(path);
        }
    }
    Ok(sketchlist_filenames)
}

pub fn load_fasta_fromfile<P: AsRef<Path>>(
    sketchlist_filename: &P,
) -> Result<Vec<(String, PathBuf, String)>> {
    let mut rdr = csv::Reader::from_path(sketchlist_filename)?;

    // Check for right header
    let headers = rdr.headers()?;
    if headers.len() != 3
        || headers.get(0).unwrap() != "name"
        || headers.get(1).unwrap() != "genome_filename"
        || headers.get(2).unwrap() != "protein_filename"
    {
        return Err(anyhow!(
            "Invalid header. Expected 'name,genome_filename,protein_filename', but got '{}'",
            headers.iter().collect::<Vec<_>>().join(",")
        ));
    }

    let mut results = Vec::new();

    let mut row_count = 0;
    let mut genome_count = 0;
    let mut protein_count = 0;
    // Create a HashSet to keep track of processed rows.
    let mut processed_rows = std::collections::HashSet::new();
    let mut duplicate_count = 0;

    for result in rdr.records() {
        let record = result?;

        // Skip duplicated rows
        let row_string = record.iter().collect::<Vec<_>>().join(",");
        if processed_rows.contains(&row_string) {
            duplicate_count += 1;
            continue;
        }
        processed_rows.insert(row_string.clone());
        row_count += 1;
        let name = record
            .get(0)
            .ok_or_else(|| anyhow!("Missing 'name' field"))?
            .to_string();

        let genome_filename = record
            .get(1)
            .ok_or_else(|| anyhow!("Missing 'genome_filename' field"))?;
        if !genome_filename.is_empty() {
            results.push((
                name.clone(),
                PathBuf::from(genome_filename),
                "dna".to_string(),
            ));
            genome_count += 1;
        }

        let protein_filename = record
            .get(2)
            .ok_or_else(|| anyhow!("Missing 'protein_filename' field"))?;
        if !protein_filename.is_empty() {
            results.push((name, PathBuf::from(protein_filename), "protein".to_string()));
            protein_count += 1;
        }
    }
    // Print warning if there were duplicated rows.
    if duplicate_count > 0 {
        println!("Warning: {} duplicated rows were skipped.", duplicate_count);
    }
    println!(
        "Loaded {} rows in total ({} genome and {} protein files)",
        row_count, genome_count, protein_count
    );
    Ok(results)
}

pub fn load_mh_with_name_and_md5(
    collection: Collection,
    selection: &Selection,
    report_type: ReportType,
) -> Result<Vec<(KmerMinHash, String, String)>> {
    let mut sketchinfo: Vec<(KmerMinHash, String, String)> = Vec::new();
    for (_idx, record) in collection.iter() {
        if let Ok(sig) = collection.sig_from_record(record) {
            if let Some(ds_mh) = sig.clone().select(&selection)?.minhash().cloned() {
                sketchinfo.push((ds_mh, record.name().to_string(), record.md5().to_string()));
            }
        } else {
            bail!(
                "Error: Failed to load {} record: {}",
                report_type,
                record.name()
            );
        }
    }
    Ok(sketchinfo)
}

/// Load a collection of sketches from a file, filtering to keep only
/// those with a minimum overlap.

pub fn load_sketches_above_threshold(
    against_collection: Collection,
    selection: &Selection,
    query: &KmerMinHash,
    threshold_hashes: u64,
) -> Result<(BinaryHeap<PrefetchResult>, usize, usize)> {
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let matchlist: BinaryHeap<PrefetchResult> = against_collection
        .par_iter()
        .filter_map(|(_idx, against_record)| {
            let mut results = Vec::new();
            // Load against into memory
            if let Ok(against_sig) = against_collection.sig_from_record(against_record) {
                for sketch in against_sig.sketches() {
                    if let Sketch::MinHash(against_mh) = sketch {
                        // currently downsampling here to avoid changing md5sum
                        if let Ok(overlap) = against_mh.count_common(query, true) {
                            if overlap >= threshold_hashes {
                                let result = PrefetchResult {
                                    name: against_record.name().to_string(),
                                    md5sum: against_mh.md5sum(),
                                    minhash: against_mh.clone(),
                                    overlap,
                                };
                                results.push(result);
                            }
                        }
                    } else {
                        eprintln!(
                            "WARNING: no compatible sketches in path '{}'",
                            against_sig.filename()
                        );
                        let _i = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                }
            } else {
                // this shouldn't happen here anymore -- likely would happen at load_collection
                eprintln!(
                    "WARNING: could not load sketches for record '{}'",
                    against_record.internal_location()
                );
                let _i = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
            }
            if results.is_empty() {
                None
            } else {
                Some(results)
            }
        })
        .flatten()
        .collect();

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    Ok((matchlist, skipped_paths, failed_paths))
}

pub enum ReportType {
    Query,
    Against,
    Pairwise,
}

impl std::fmt::Display for ReportType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let description = match self {
            ReportType::Query => "query",
            ReportType::Against => "search",
            ReportType::Pairwise => "signature",
        };
        write!(f, "{}", description)
    }
}

pub fn load_collection(
    sigpath: &camino::Utf8PathBuf,
    selection: &Selection,
    report_type: ReportType,
) -> Result<Collection> {
    if !sigpath.exists() {
        bail!("No such file or directory: '{}'", sigpath);
    }

    let mut n_failed = 0;
    let collection = if sigpath.extension().map_or(false, |ext| ext == "zip") {
        match Collection::from_zipfile(&sigpath) {
            Ok(collection) => collection,
            Err(_) => bail!("failed to load {} zipfile: '{}'", report_type, sigpath),
        }
    } else {
        // if pathlist is just a signature path, load it into a collection
        match Signature::from_path(sigpath) {
            Ok(signatures) => {
                // Load the collection from the signature
                match Collection::from_sigs(signatures) {
                    Ok(collection) => collection,
                    Err(_) => bail!(
                        "loaded {} signatures but failed to load as collection: '{}'",
                        report_type,
                        sigpath
                    ),
                }
            }
            // if not, try to load file as list of sig paths
            Err(_) => {
                //             // using core fn doesn't allow us to ignore failed paths; I reimplement loading here to allow
                let sketchlist_file = BufReader::new(File::open(sigpath)?);
                let records: Vec<Record> = sketchlist_file
                    .lines()
                    .filter_map(|line| {
                        let path = line.ok()?;
                        match Signature::from_path(&path) {
                            Ok(signatures) => {
                                let recs: Vec<Record> = signatures
                                    .into_iter()
                                    .flat_map(|v| Record::from_sig(&v, &path))
                                    .collect();
                                Some(recs)
                            }
                            Err(err) => {
                                eprintln!("Sketch loading error: {}", err);
                                eprintln!("WARNING: could not load sketches from path '{}'", path);
                                n_failed += 1;
                                None
                            }
                        }
                    })
                    .flatten()
                    .collect();

                let manifest: Manifest = records.into();
                Collection::new(
                    manifest,
                    InnerStorage::new(
                        FSStorage::builder()
                            .fullpath("".into())
                            .subdir("".into())
                            .build(),
                    ),
                )
            }
        }
    };

    let n_total = collection.len();
    let selected = collection.select(selection)?;
    let n_skipped = n_total - selected.len();
    report_on_collection_loading(&selected, n_skipped, n_failed, report_type)?;
    Ok(selected)
}

/// Uses the output of collection loading function to report the
/// total number of sketches loaded, as well as the number of files,
/// if any, that failed to load or contained no compatible sketches.
/// If no sketches were loaded, bail.
///
/// # Arguments
///
/// * `sketchlist` - A slice of loaded `SmallSignature` sketches.
/// * `skipped_paths` - # paths that contained no compatible sketches.
/// * `failed_paths` - # paths that failed to load.
/// * `report_type` - ReportType Enum (Query or Against). Used to specify
///                   which sketch input this information pertains to.
///
/// # Returns
///
/// Returns `Ok(())` if at least one signature was successfully loaded.
/// Returns an error if no signatures were loaded.
///
/// # Errors
///
/// Returns an error if:
/// * No signatures were successfully loaded.
pub fn report_on_collection_loading(
    collection: &Collection,
    skipped_paths: usize,
    failed_paths: usize,
    report_type: ReportType,
) -> Result<()> {
    if failed_paths > 0 {
        eprintln!(
            "WARNING: {} {} paths failed to load. See error messages above.",
            failed_paths, report_type
        );
    }
    if skipped_paths > 0 {
        eprintln!(
            "WARNING: skipped {} {} paths - no compatible signatures.",
            skipped_paths, report_type
        );
    }

    // Validate sketches
    if collection.is_empty() {
        bail!("No {} signatures loaded, exiting.", report_type);
    }
    eprintln!("Loaded {} {} signature(s)", collection.len(), report_type);
    Ok(())
}

/// Execute the gather algorithm, greedy min-set-cov, by iteratively
/// removing matches in 'matchlist' from 'query'.

pub fn consume_query_by_gather(
    query: SigStore,
    matchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
    gather_output: Option<String>,
) -> Result<()> {
    // Define the writer to stdout by default
    let mut writer: Box<dyn Write> = Box::new(std::io::stdout());

    if let Some(output_path) = &gather_output {
        // Account for potential missing dir in output path
        let directory_path = Path::new(output_path).parent();

        // If a directory path exists in the filename, create it if it doesn't already exist
        if let Some(dir) = directory_path {
            create_dir_all(dir)?;
        }

        let file = File::create(output_path)?;
        writer = Box::new(BufWriter::new(file));
    }
    writeln!(
        &mut writer,
        "query_filename,rank,query_name,query_md5,match_name,match_md5,intersect_bp"
    )
    .ok();

    let mut matching_sketches = matchlist;
    let mut rank = 0;

    let mut last_matches = matching_sketches.len();

    // let location = query.location;
    let location = query.filename(); // this is different (original fasta filename) than query.location was (sig name)!!

    let sketches = query.sketches();
    let orig_query_mh = match sketches.get(0) {
        Some(Sketch::MinHash(mh)) => Ok(mh),
        _ => Err(anyhow::anyhow!("No MinHash found")),
    }?;
    let mut query_mh = orig_query_mh.clone();
    let mut last_hashes = orig_query_mh.size();

    eprintln!(
        "{} iter {}: start: query hashes={} matches={}",
        location,
        rank,
        orig_query_mh.size(),
        matching_sketches.len()
    );

    while !matching_sketches.is_empty() {
        let best_element = matching_sketches.peek().unwrap();

        // remove!
        query_mh.remove_from(&best_element.minhash)?;

        writeln!(
            &mut writer,
            "{},{},\"{}\",{},\"{}\",{},{}",
            location,
            rank,
            query.name(),
            query.md5sum(),
            best_element.name,
            best_element.md5sum,
            best_element.overlap
        )
        .ok();

        // recalculate remaining overlaps between query and all sketches.
        // note: this is parallelized.
        matching_sketches = prefetch(&query_mh, matching_sketches, threshold_hashes);
        rank += 1;

        let sub_hashes = last_hashes - query_mh.size();
        let sub_matches = last_matches - matching_sketches.len();

        eprintln!(
            "{} iter {}: remaining: query hashes={}(-{}) matches={}(-{})",
            location,
            rank,
            query_mh.size(),
            sub_hashes,
            matching_sketches.len(),
            sub_matches
        );

        last_hashes = query_mh.size();
        last_matches = matching_sketches.len();
    }
    Ok(())
}

pub fn build_selection(ksize: u8, scaled: usize, moltype: &str) -> Selection {
    let hash_function = match moltype {
        "dna" => HashFunctions::Murmur64Dna,
        "protein" => HashFunctions::Murmur64Protein,
        "dayhoff" => HashFunctions::Murmur64Dayhoff,
        "hp" => HashFunctions::Murmur64Hp,
        _ => panic!("Unknown molecule type: {}", moltype),
    };
    // let hash_function = HashFunctions::try_from(moltype)
    //     .map_err(|_| panic!("Unknown molecule type: {}", moltype))
    //     .unwrap();

    Selection::builder()
        .ksize(ksize.into())
        .scaled(scaled as u32)
        .moltype(hash_function)
        .build()
}

pub fn is_revindex_database(path: &camino::Utf8PathBuf) -> bool {
    // quick file check for Revindex database:
    // is path a directory that contains a file named 'CURRENT'?
    if path.is_dir() {
        let current_file = path.join("CURRENT");
        current_file.exists() && current_file.is_file()
    } else {
        false
    }
}

pub struct SearchResult {
    pub query_name: String,
    pub query_md5: String,
    pub match_name: String,
    pub containment: f64,
    pub intersect_hashes: usize,
    pub match_md5: Option<String>,
    pub jaccard: Option<f64>,
    pub max_containment: Option<f64>,
}

impl ResultType for SearchResult {
    fn header_fields() -> Vec<&'static str> {
        vec![
            "query_name",
            "query_md5",
            "match_name",
            "containment",
            "intersect_hashes",
            "match_md5",
            "jaccard",
            "max_containment",
        ]
    }

    fn format_fields(&self) -> Vec<String> {
        vec![
            format!("\"{}\"", self.query_name), // Wrap query_name with quotes
            self.query_md5.clone(),
            format!("\"{}\"", self.match_name), // Wrap match_name with quotes
            self.containment.to_string(),
            self.intersect_hashes.to_string(),
            match &self.match_md5 {
                Some(md5) => md5.clone(),
                None => "".to_string(),
            },
            match &self.jaccard {
                Some(jaccard) => jaccard.to_string(),
                None => "".to_string(),
            },
            match &self.max_containment {
                Some(max_containment) => max_containment.to_string(),
                None => "".to_string(),
            },
        ]
    }
}

pub struct ManifestRow {
    pub md5: String,
    pub md5short: String,
    pub ksize: u32,
    pub moltype: String,
    pub num: u32,
    pub scaled: u64,
    pub n_hashes: usize,
    pub with_abundance: bool,
    pub name: String,
    pub filename: String,
    pub internal_location: String,
}

pub fn bool_to_python_string(b: bool) -> String {
    match b {
        true => "True".to_string(),
        false => "False".to_string(),
    }
}

impl ResultType for ManifestRow {
    fn header_fields() -> Vec<&'static str> {
        vec![
            "internal_location",
            "md5",
            "md5short",
            "ksize",
            "moltype",
            "num",
            "scaled",
            "n_hashes",
            "with_abundance",
            "name",
            "filename",
        ]
    }

    fn format_fields(&self) -> Vec<String> {
        vec![
            self.internal_location.clone(),
            self.md5.clone(),
            self.md5short.clone(),
            self.ksize.to_string(),
            self.moltype.clone(),
            self.num.to_string(),
            self.scaled.to_string(),
            self.n_hashes.to_string(),
            bool_to_python_string(self.with_abundance),
            format!("\"{}\"", self.name), // Wrap name with quotes
            self.filename.clone(),
        ]
    }
}

pub fn make_manifest_row(
    sig: &Signature,
    filename: &Path,
    internal_location: &str,
    scaled: u64,
    num: u32,
    abund: bool,
    is_dna: bool,
    is_protein: bool,
) -> ManifestRow {
    if is_dna && is_protein {
        panic!("Both is_dna and is_protein cannot be true at the same time.");
    } else if !is_dna && !is_protein {
        panic!("Either is_dna or is_protein must be true.");
    }
    let moltype = if is_dna {
        "DNA".to_string()
    } else {
        "protein".to_string()
    };
    let sketch = &sig.sketches()[0];
    ManifestRow {
        internal_location: internal_location.to_string(),
        md5: sig.md5sum(),
        md5short: sig.md5sum()[0..8].to_string(),
        ksize: sketch.ksize() as u32,
        moltype,
        num,
        scaled,
        n_hashes: sketch.size(),
        with_abundance: abund,
        name: sig.name().to_string(),
        // filename: filename.display().to_string(),
        filename: filename.to_str().unwrap().to_string(),
    }
}

pub fn open_stdout_or_file<P: AsRef<Path>>(output: Option<P>) -> Box<dyn Write + Send + 'static> {
    // if output is a file, use open_output_file
    if let Some(path) = output {
        Box::new(open_output_file(&path))
    } else {
        Box::new(std::io::stdout())
    }
}

pub fn open_output_file<P: AsRef<Path>>(output: &P) -> BufWriter<File> {
    let file = File::create(output).unwrap_or_else(|e| {
        eprintln!("Error creating output file: {:?}", e);
        std::process::exit(1);
    });
    BufWriter::new(file)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Params {
    pub ksize: u32,
    pub track_abundance: bool,
    pub num: u32,
    pub scaled: u64,
    pub seed: u32,
    pub is_protein: bool,
    pub is_dna: bool,
}
use std::hash::Hash;
use std::hash::Hasher;

impl Hash for Params {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ksize.hash(state);
        self.track_abundance.hash(state);
        self.num.hash(state);
        self.scaled.hash(state);
        self.seed.hash(state);
        self.is_protein.hash(state);
        self.is_dna.hash(state);
    }
}

pub enum ZipMessage {
    SignatureData(Vec<Signature>, Vec<Params>, PathBuf),
    WriteManifest,
}

pub fn sigwriter<P: AsRef<Path> + Send + 'static>(
    recv: std::sync::mpsc::Receiver<ZipMessage>,
    output: String,
) -> std::thread::JoinHandle<Result<()>> {
    std::thread::spawn(move || -> Result<()> {
        let file_writer = open_output_file(&output);

        let options = zip::write::FileOptions::default()
            .compression_method(zip::CompressionMethod::Stored)
            .large_file(true);
        let mut zip = zip::ZipWriter::new(file_writer);
        let mut manifest_rows: Vec<ManifestRow> = Vec::new();
        // keep track of md5sum occurrences to prevent overwriting duplicates
        let mut md5sum_occurrences: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();

        while let Ok(message) = recv.recv() {
            match message {
                ZipMessage::SignatureData(sigs, params, filename) => {
                    if sigs.len() != params.len() {
                        bail!("Mismatched lengths of signatures and parameters");
                    }
                    for (sig, param) in sigs.iter().zip(params.iter()) {
                        let md5sum_str = sig.md5sum();
                        let count = md5sum_occurrences.entry(md5sum_str.clone()).or_insert(0);
                        *count += 1;
                        let sig_filename = if *count > 1 {
                            format!("signatures/{}_{}.sig.gz", md5sum_str, count)
                        } else {
                            format!("signatures/{}.sig.gz", md5sum_str)
                        };
                        write_signature(sig, &mut zip, options, &sig_filename);
                        manifest_rows.push(make_manifest_row(
                            sig,
                            &filename,
                            &sig_filename,
                            param.scaled,
                            param.num,
                            param.track_abundance,
                            param.is_dna,
                            param.is_protein,
                        ));
                    }
                }
                ZipMessage::WriteManifest => {
                    println!("Writing manifest");
                    // Start the CSV file inside the zip
                    zip.start_file("SOURMASH-MANIFEST.csv", options).unwrap();

                    // write manifest version line
                    writeln!(&mut zip, "# SOURMASH-MANIFEST-VERSION: 1.0").unwrap();
                    // Write the header
                    let header = ManifestRow::header_fields();
                    if let Err(e) = writeln!(&mut zip, "{}", header.join(",")) {
                        eprintln!("Error writing header: {:?}", e);
                    }

                    // Write each manifest row
                    for row in &manifest_rows {
                        let formatted_fields = row.format_fields(); // Assuming you have a format_fields method on ManifestRow
                        if let Err(e) = writeln!(&mut zip, "{}", formatted_fields.join(",")) {
                            eprintln!("Error writing item: {:?}", e);
                        }
                    }
                    // finalize the zip file writing.
                    zip.finish().unwrap();
                }
            }
        }
        Ok(())
    })
}

pub trait ResultType {
    fn header_fields() -> Vec<&'static str>;
    fn format_fields(&self) -> Vec<String>;
}

pub fn csvwriter_thread<T: ResultType + Send + 'static, P>(
    recv: std::sync::mpsc::Receiver<T>,
    output: Option<P>,
) -> std::thread::JoinHandle<()>
where
    T: ResultType,
    P: Clone + std::convert::AsRef<std::path::Path>,
{
    // create output file
    let out = open_stdout_or_file(output.as_ref());
    // spawn a thread that is dedicated to printing to a buffered output
    std::thread::spawn(move || {
        let mut writer = out;

        let header = T::header_fields();
        if let Err(e) = writeln!(&mut writer, "{}", header.join(",")) {
            eprintln!("Error writing header: {:?}", e);
        }
        writer.flush().unwrap();

        for item in recv.iter() {
            let formatted_fields = item.format_fields();
            if let Err(e) = writeln!(&mut writer, "{}", formatted_fields.join(",")) {
                eprintln!("Error writing item: {:?}", e);
            }
            writer.flush().unwrap();
        }
    })
}

pub fn write_signature(
    sig: &Signature,
    zip: &mut zip::ZipWriter<BufWriter<File>>,
    zip_options: zip::write::FileOptions,
    sig_filename: &str,
) {
    let wrapped_sig = vec![sig];
    let json_bytes = serde_json::to_vec(&wrapped_sig).unwrap();

    let gzipped_buffer = {
        let mut buffer = std::io::Cursor::new(Vec::new());
        {
            let mut gz_writer = niffler::get_writer(
                Box::new(&mut buffer),
                niffler::compression::Format::Gzip,
                niffler::compression::Level::Nine,
            )
            .unwrap();
            gz_writer.write_all(&json_bytes).unwrap();
        }
        buffer.into_inner()
    };

    zip.start_file(sig_filename, zip_options).unwrap();
    zip.write_all(&gzipped_buffer).unwrap();
}
