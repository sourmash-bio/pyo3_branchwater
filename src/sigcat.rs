/// sigcat: concatenate signatures into a single sourmash zip file
use anyhow::Result;
use sourmash::{collection::Collection, selection::Selection};

use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use camino::Utf8PathBuf as PathBuf;
use sourmash::manifest::{Manifest, Record};
use sourmash::prelude::Select;
use zip::write::{FileOptions, ZipWriter};
use zip::CompressionMethod;

use crate::utils::{load_collection, open_output_file, write_signature, ReportType};

pub fn sig_cat(
    input_sigs: String,
    selection: &Selection,
    allow_failed_sigpaths: bool,
    output_path: String,
    allow_empty_collection: bool,
) -> Result<()> {
    // Check if output_path ends with ".zip"
    if !output_path.ends_with(".zip") {
        return Err(anyhow::anyhow!("Output file must end with '.zip'"));
    }

    // split input_sigs string into list
    let sig_inputs: Vec<&str> = input_sigs.split(' ').collect();

    let mut collection_list: Vec<Collection> = Vec::new();
    // read each input collection
    for sigs in sig_inputs {
        let collection = load_collection(
            &sigs.to_string(),
            selection,
            ReportType::General,
            allow_failed_sigpaths,
            allow_empty_collection,
        )?;
        if !collection.is_empty() {
            collection_list.push(collection);
        }
    }
    if collection_list.is_empty() {
        bail!("No signatures to concatenate, exiting")
    }
    let collected_sigs = AtomicUsize::new(0);

    // open new zipfile and write all sigs
    let outpath: PathBuf = output_path.into();

    let file_writer = open_output_file(&outpath);

    let options = FileOptions::default()
        .compression_method(CompressionMethod::Stored)
        .unix_permissions(0o644)
        .large_file(true);

    let mut zip = ZipWriter::new(file_writer);
    let mut manifest_rows: Vec<Record> = Vec::new();
    // keep track of MD5 sum occurrences to prevent overwriting duplicates
    let mut md5sum_occurrences: HashMap<String, usize> = HashMap::new();

    for collection in collection_list {
        collection.iter().for_each(|(_idx, record)| {
            // todo: count the number we're adding? or count failures?
            let _i = collected_sigs.fetch_add(1, Ordering::SeqCst);
            // if scaled doesn't match sig scaled, we need to downsample.
            let sig = collection
                .sig_from_record(record)
                .expect("failed to get sig from record")
                .select(selection)
                .unwrap();

            let md5sum_str = sig.md5sum(); // this is now the downsampled md5sum -- okay?

            let count = md5sum_occurrences.entry(md5sum_str.clone()).or_insert(0);
            *count += 1;
            let sig_filename = if *count > 1 {
                format!("signatures/{}_{}.sig.gz", md5sum_str, count)
            } else {
                format!("signatures/{}.sig.gz", md5sum_str)
            };
            write_signature(&sig, &mut zip, options.clone(), &sig_filename);
            // need to build new record to reflect new sig location
            let records: Vec<Record> = Record::from_sig(&sig, sig_filename.as_str());
            manifest_rows.extend(records);
        });
    }
    // now write manifest
    println!("Writing manifest");
    zip.start_file("SOURMASH-MANIFEST.csv", options)?;
    let manifest: Manifest = manifest_rows.clone().into();
    manifest.to_writer(&mut zip)?;
    zip.finish()?;

    // report
    let i: usize = collected_sigs.load(Ordering::SeqCst);
    println!("Concatenated {} signatures into '{}'.", i, outpath);

    Ok(())
}
