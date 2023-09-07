/// Utility functions for pyo3_branchwater.

use rayon::prelude::*;

use std::fs::File;
use std::io::Read;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use zip::read::ZipArchive;
use tempfile::tempdir;

use std::sync::atomic;
use std::sync::atomic::AtomicUsize;

use std::collections::BinaryHeap;

use anyhow::{Result, anyhow};

use std::cmp::{PartialOrd, Ordering};

use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;
use sourmash::prelude::MinHashOps;
use sourmash::prelude::FracMinHashOps;

/// Track a name/minhash.

pub struct SmallSignature {
    pub location: String,
    pub name: String,
    pub md5sum: String,
    pub minhash: KmerMinHash,
}

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

/// check to see if two KmerMinHash are compatible.
///
/// CTB note: despite the name, downsampling is not performed?
/// Although it checks if they are compatible in one direction...

pub fn check_compatible_downsample(
    me: &KmerMinHash,
    other: &KmerMinHash,
) -> Result<(), sourmash::Error> {
    /* // ignore num minhashes.
    if self.num != other.num {
        return Err(Error::MismatchNum {
            n1: self.num,
            n2: other.num,
        }
        .into());
    }
    */
    use sourmash::Error;

    if me.ksize() != other.ksize() {
        return Err(Error::MismatchKSizes);
    }
    if me.hash_function() != other.hash_function() {
        // TODO: fix this error
        return Err(Error::MismatchDNAProt);
    }
    if me.max_hash() < other.max_hash() {
        return Err(Error::MismatchScaled);
    }
    if me.seed() != other.seed() {
        return Err(Error::MismatchSeed);
    }
    Ok(())
}


/// Given a vec of search Signatures, each containing one or more sketches,
/// and a template Sketch, return a compatible (& now downsampled)
/// Sketch from the search Signatures..
///
/// CTB note: this will return the first acceptable match, I think, ignoring
/// all others.


pub fn prepare_query(search_sigs: &[Signature], template: &Sketch, location: &str) -> Option<SmallSignature> {

    for search_sig in search_sigs.iter() {
        // find exact match for template?
        if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(template) {
            return Some(SmallSignature {
                location: location.to_string().clone(),
                name: search_sig.name(),
                md5sum: mh.md5sum(),
                minhash: mh.clone()
            });
        } else {
            // no - try to find one that can be downsampled
            if let Sketch::MinHash(template_mh) = template {
                for sketch in search_sig.sketches() {
                    if let Sketch::MinHash(ref_mh) = sketch {
                        if check_compatible_downsample(&ref_mh, template_mh).is_ok() {
                            let max_hash = max_hash_for_scaled(template_mh.scaled());
                            let mh = ref_mh.downsample_max_hash(max_hash).unwrap();
                            return Some(SmallSignature {
                                location: location.to_string().clone(),
                                name: search_sig.name(),
                                md5sum: ref_mh.md5sum(), // original
                                minhash: mh,             // downsampled
                            });
                        }
                    }
                }
            }
        }
    }
    None
}

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
            let overlap = searchsig.count_common(query_mh, false);
            if let Ok(overlap) = overlap {
                if overlap >= threshold_hashes {
                    let result = PrefetchResult {
                        overlap,
                        ..result
                    };
                    mm = Some(result);
                }
            }
            mm
        })
        .collect()
}

/// Write list of prefetch matches.

pub fn write_prefetch<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display + Clone>(
    query: &SmallSignature,
    prefetch_output: Option<P>,
    matchlist: &BinaryHeap<PrefetchResult>
) -> Result<()> {
    // Set up a writer for prefetch output
    let prefetch_out: Box<dyn Write> = match prefetch_output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(prefetch_out);
    writeln!(&mut writer, "query_filename,query_name,query_md5,match_name,match_md5,intersect_bp").ok();

    for m in matchlist.iter() {
        writeln!(&mut writer, "{},\"{}\",{},\"{}\",{},{}", query.location,
                 query.name, query.md5sum,
                 m.name, m.md5sum, m.overlap).ok();
    }

    Ok(())
}

/// Load a list of filenames from a file. Exits on bad lines.

pub fn load_sketchlist_filenames<P: AsRef<Path>>(sketchlist_filename: &P) ->
    Result<Vec<PathBuf>>
{
    let sketchlist_file = BufReader::new(File::open(sketchlist_filename)?);

    let mut sketchlist_filenames : Vec<PathBuf> = Vec::new();
    for line in sketchlist_file.lines() {
        let line = match line {
            Ok(v) => v,
            Err(_) => return {
                let filename = sketchlist_filename.as_ref().display();
                let msg = format!("invalid line in fromfile '{}'", filename);
                Err(anyhow!(msg))
            },
        };

        if !line.is_empty() {
            let mut path = PathBuf::new();
            path.push(line);
            sketchlist_filenames.push(path);
        }
    }
    Ok(sketchlist_filenames)
}

pub fn load_sketch_fromfile<P: AsRef<Path>>(sketchlist_filename: &P) -> Result<Vec<(String, PathBuf, String)>> {
    let mut rdr = csv::Reader::from_path(sketchlist_filename)?;

    // Check for right header
    let headers = rdr.headers()?;
    if headers.len() != 3 ||
    headers.get(0).unwrap() != "name" ||
    headers.get(1).unwrap() != "genome_filename" ||
    headers.get(2).unwrap() != "protein_filename" {
        return Err(anyhow!("Invalid header. Expected 'name,genome_filename,protein_filename', but got '{}'", headers.iter().collect::<Vec<_>>().join(",")));
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
        let name = record.get(0).ok_or_else(|| anyhow!("Missing 'name' field"))?.to_string();

        let genome_filename = record.get(1).ok_or_else(|| anyhow!("Missing 'genome_filename' field"))?;
        if !genome_filename.is_empty() {
            results.push((name.clone(), PathBuf::from(genome_filename), "dna".to_string()));
            genome_count += 1;
        }

        let protein_filename = record.get(2).ok_or_else(|| anyhow!("Missing 'protein_filename' field"))?;
        if !protein_filename.is_empty() {
            results.push((name, PathBuf::from(protein_filename), "protein".to_string()));
            protein_count += 1;
        }
    }
    // Print warning if there were duplicated rows.
    if duplicate_count > 0 {
        println!("Warning: {} duplicated rows were skipped.", duplicate_count);
    }
    println!("Loaded {} rows in total ({} genome and {} protein files)", row_count, genome_count, protein_count);
    Ok(results)
}


/// Load a collection of sketches from a file in parallel.
pub fn load_sketches(sketchlist_paths: Vec<PathBuf>, template: &Sketch) ->
    Result<(Vec<SmallSignature>, usize, usize)>
{
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let sketchlist : Vec<SmallSignature> = sketchlist_paths
        .par_iter()
        .filter_map(|m| {
            let filename = m.display().to_string();

            match Signature::from_path(m) {
                Ok(sigs) => {
                    let sm = prepare_query(&sigs, template, &filename);
                    if sm.is_none() {
                        // track number of paths that have no matching sigs
                        let _i = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                    sm
                },
                Err(err) => {
                    // failed to load from this path - print error & track.
                    eprintln!("Sketch loading error: {}", err);
                    eprintln!("WARNING: could not load sketches from path '{}'", filename);
                    let _i = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    None
                }
            }
        })
        .collect();

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);
    Ok((sketchlist, skipped_paths, failed_paths))
}

/// Load a collection of sketches from a file, filtering to keep only
/// those with a minimum overlap.

pub fn load_sketches_above_threshold(
    sketchlist_paths: Vec<PathBuf>,
    template: &Sketch,
    query: &KmerMinHash,
    threshold_hashes: u64
) ->
    Result<(BinaryHeap<PrefetchResult>, usize, usize)>
{
    let skipped_paths = AtomicUsize::new(0);
    let failed_paths = AtomicUsize::new(0);

    let matchlist: BinaryHeap<PrefetchResult> = sketchlist_paths
        .par_iter()
        .filter_map(|m| {
            let sigs = Signature::from_path(m);
            let location = m.display().to_string();

            match sigs {
                Ok(sigs) => {
                    let mut mm = None;

                    if let Some(sm) = prepare_query(&sigs, template,
                                                           &location) {
                        let mh = sm.minhash;
                        if let Ok(overlap) = mh.count_common(query, false) {
                            if overlap >= threshold_hashes {
                                let result = PrefetchResult {
                                    name: sm.name,
                                    md5sum: sm.md5sum,
                                    minhash: mh,
                                    overlap,
                                };
                                mm = Some(result);
                            }
                        }
                    } else {
                        eprintln!("WARNING: no compatible sketches in path '{}'",
                                  m.display());
                        let _i = skipped_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    }
                    mm
                }
                Err(err) => {
                    eprintln!("Sketch loading error: {}", err);
                    let _ = failed_paths.fetch_add(1, atomic::Ordering::SeqCst);
                    eprintln!("WARNING: could not load sketches from path '{}'",
                          m.display());
                    None
                }
            }
        })
        .collect();

    let skipped_paths = skipped_paths.load(atomic::Ordering::SeqCst);
    let failed_paths = failed_paths.load(atomic::Ordering::SeqCst);

    Ok((matchlist, skipped_paths, failed_paths))
}

/// Execute the gather algorithm, greedy min-set-cov, by iteratively
/// removing matches in 'matchlist' from 'query'.

pub fn consume_query_by_gather<P: AsRef<Path> + std::fmt::Debug + std::fmt::Display + Clone>(
    query: SmallSignature,
    matchlist: BinaryHeap<PrefetchResult>,
    threshold_hashes: u64,
    gather_output: Option<P>,
) -> Result<()> {
    // Set up a writer for gather output
    let gather_out: Box<dyn Write> = match gather_output {
        Some(path) => Box::new(BufWriter::new(File::create(path).unwrap())),
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(gather_out);
    writeln!(&mut writer, "query_filename,rank,query_name,query_md5,match_name,match_md5,intersect_bp").ok();

    let mut matching_sketches = matchlist;
    let mut rank = 0;

    let mut last_hashes = query.minhash.size();
    let mut last_matches = matching_sketches.len();

    let location = query.location;
    let mut query_mh = query.minhash;

    eprintln!("{} iter {}: start: query hashes={} matches={}", location, rank,
              query_mh.size(), matching_sketches.len());

    while !matching_sketches.is_empty() {
        let best_element = matching_sketches.peek().unwrap();

        // remove!
        query_mh.remove_from(&best_element.minhash)?;

        writeln!(&mut writer, "{},{},\"{}\",{},\"{}\",{},{}", location, rank,
                 query.name, query.md5sum,
                 best_element.name, best_element.md5sum,
                 best_element.overlap).ok();

        // recalculate remaining overlaps between query and all sketches.
        // note: this is parallelized.
        matching_sketches = prefetch(&query_mh, matching_sketches, threshold_hashes);
        rank += 1;

        let sub_hashes = last_hashes - query_mh.size();
        let sub_matches = last_matches - matching_sketches.len();

        eprintln!("{} iter {}: remaining: query hashes={}(-{}) matches={}(-{})", location, rank,
            query_mh.size(), sub_hashes, matching_sketches.len(), sub_matches);

        last_hashes = query_mh.size();
        last_matches = matching_sketches.len();

    }
    Ok(())
}
  

  // mastiff rocksdb functions

pub fn build_template(ksize: u8, scaled: usize) -> Sketch {
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    Sketch::MinHash(template_mh)
}

pub fn read_signatures_from_zip<P: AsRef<Path>>(
    zip_path: P,
) -> Result<(Vec<PathBuf>, tempfile::TempDir), Box<dyn std::error::Error>> {
    let mut signature_paths = Vec::new();
    let temp_dir = tempdir()?;
    let zip_file = File::open(&zip_path)?;
    let mut zip_archive = ZipArchive::new(zip_file)?;

    for i in 0..zip_archive.len() {
        let mut file = zip_archive.by_index(i)?;
        let mut sig = Vec::new();
        file.read_to_end(&mut sig)?;

        let file_name = Path::new(file.name()).file_name().unwrap().to_str().unwrap();
        if file_name.ends_with(".sig") || file_name.ends_with(".sig.gz") {
            println!("Found signature file: {}", file_name);
            let mut new_file = File::create(temp_dir.path().join(file_name))?;
            new_file.write_all(&sig)?;

            // Push the created path directly to the vector
            signature_paths.push(temp_dir.path().join(file_name));
        }
    }
    println!("wrote {} signatures to temp dir", signature_paths.len());
    Ok((signature_paths, temp_dir))
}

pub fn is_revindex_database(path: &Path) -> bool {
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
        vec!["query_name", "query_md5", "match_name", "containment", "intersect_hashes", "match_md5", "jaccard", "max_containment"]
    }

    fn format_fields(&self) -> Vec<String> {
        vec![
            format!("\"{}\"", self.query_name),  // Wrap query_name with quotes
            self.query_md5.clone(),
            format!("\"{}\"", self.match_name),  // Wrap match_name with quotes
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
            }
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
        vec!["internal_location", "md5", "md5short", "ksize", "moltype", "num", "scaled", "n_hashes", "with_abundance", "name", "filename"]
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
            format!("\"{}\"", self.name),  // Wrap name with quotes
            self.filename.clone(),
        ]
    }
}

pub fn make_manifest_row(sig: &Signature, filename: &Path, internal_location: &str, scaled: u64, num: u32, abund: bool, is_dna: bool, is_protein: bool) -> ManifestRow {
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

pub fn open_stdout_or_file<P: AsRef<Path>>(
    output: Option<P>
) -> Box<dyn Write + Send + 'static> {
    // if output is a file, use open_output_file
    if let Some(path) = output {
        Box::new(open_output_file(&path))
    } else {
        Box::new(std::io::stdout())
    }
}

pub fn open_output_file<P: AsRef<Path>>(
    output: &P
) -> BufWriter<File> {
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

        let options = zip::write::FileOptions::default().compression_method(zip::CompressionMethod::Stored);
        let mut zip = zip::ZipWriter::new(file_writer);
        let mut manifest_rows: Vec<ManifestRow> = Vec::new();
        // keep track of md5sum occurrences to prevent overwriting duplicates
        let mut md5sum_occurrences: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

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
                        manifest_rows.push(make_manifest_row(sig, &filename, &sig_filename, param.scaled, param.num, param.track_abundance, param.is_dna, param.is_protein));
                    }
                },
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
                        let formatted_fields = row.format_fields();  // Assuming you have a format_fields method on ManifestRow
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
            ).unwrap();
            gz_writer.write_all(&json_bytes).unwrap();
        }
        buffer.into_inner()
    };

    zip.start_file(sig_filename, zip_options).unwrap();
    zip.write_all(&gzipped_buffer).unwrap();
}