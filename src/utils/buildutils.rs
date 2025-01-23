//! sketching utilities

use anyhow::{anyhow, Context, Result};
use camino::Utf8PathBuf;
use getset::{Getters, Setters};
use needletail::parser::SequenceRecord;
use needletail::{parse_fastx_file, parse_fastx_reader, parse_fastx_stdin};
use serde::Serialize;
use sourmash::cmd::ComputeParameters;
use sourmash::encodings::{HashFunctions, Idx};
use sourmash::errors::SourmashError;
use sourmash::manifest::Record;
use sourmash::selection::Selection;
use sourmash::signature::Signature;
use sourmash::signature::SigsTrait;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Display;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{Cursor, Seek, Write};
use std::num::ParseIntError;
use std::ops::Index;
use std::str::FromStr;
use zip::write::{FileOptions, ZipWriter};
use zip::CompressionMethod;

#[derive(Default, Debug, Clone)]
pub struct MultiSelection {
    pub selections: Vec<Selection>,
}

impl MultiSelection {
    pub fn from_moltypes(moltypes: Vec<&str>) -> Result<Self, SourmashError> {
        let selections: Result<Vec<Selection>, SourmashError> = moltypes
            .into_iter()
            .map(|moltype_str| {
                let moltype = HashFunctions::try_from(moltype_str)?;
                let mut new_selection = Selection::default(); // Create a default Selection
                new_selection.set_moltype(moltype); // Set the moltype
                Ok(new_selection)
            })
            .collect();

        Ok(MultiSelection {
            selections: selections?,
        })
    }

    pub fn from_input_moltype(input_moltype: &str) -> Result<Self, SourmashError> {
        // currently we don't allow translation. Will need to change this when we do.
        // is there a better way to do this?
        let mut moltypes = vec!["DNA", "skipm1n3", "skipm2n3"]; // change so default is just dna?
        if input_moltype == "protein" {
            moltypes = vec!["protein", "dayhoff", "hp"];
        }
        let selections: Result<Vec<Selection>, SourmashError> = moltypes
            .into_iter()
            .map(|moltype_str| {
                let moltype = HashFunctions::try_from(moltype_str)?;
                let mut new_selection = Selection::default(); // Create a default Selection
                new_selection.set_moltype(moltype); // Set the moltype
                Ok(new_selection)
            })
            .collect();

        Ok(MultiSelection {
            selections: selections?,
        })
    }

    pub fn from_selection(selection: Selection) -> Self {
        MultiSelection {
            selections: vec![selection],
        }
    }
}

pub trait MultiSelect {
    fn select(&mut self, multi_selection: &MultiSelection) -> Result<(), SourmashError>;
}

#[derive(Debug, Clone, Getters, Setters, Serialize)]
pub struct BuildRecord {
    // fields are ordered the same as Record to allow serialization to manifest
    // required fields are currently immutable once set
    #[getset(get = "pub", set = "pub")]
    internal_location: Option<Utf8PathBuf>,

    #[getset(get = "pub", set = "pub")]
    md5: Option<String>,

    #[getset(get = "pub", set = "pub")]
    md5short: Option<String>,

    #[getset(get_copy = "pub", set = "pub")]
    ksize: u32,

    moltype: String,

    #[getset(get = "pub")]
    num: u32,

    #[getset(get = "pub")]
    scaled: u32,

    #[getset(get = "pub", set = "pub")]
    n_hashes: Option<usize>,

    #[getset(get_copy = "pub", set = "pub")]
    #[serde(serialize_with = "intbool")]
    with_abundance: bool,

    #[getset(get = "pub", set = "pub")]
    name: Option<String>,

    #[getset(get = "pub", set = "pub")]
    filename: Option<String>,

    #[getset(get_copy = "pub")]
    #[serde(skip)]
    pub seed: u32,

    #[serde(skip)]
    pub hashed_params: u64,

    #[serde(skip)]
    pub sequence_added: bool,
}

// from sourmash (intbool is currently private there)
fn intbool<S>(x: &bool, s: S) -> std::result::Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    if *x {
        s.serialize_i32(1)
    } else {
        s.serialize_i32(0)
    }
}

impl BuildRecord {
    // no general default, but we have defaults for each moltype
    pub fn default_dna() -> Self {
        Self {
            internal_location: None,
            md5: None,
            md5short: None,
            ksize: 31,
            moltype: "DNA".to_string(),
            num: 0,
            scaled: 1000,
            n_hashes: None,
            with_abundance: false,
            name: None,
            filename: None,
            seed: 42,
            hashed_params: 0,
            sequence_added: false,
        }
    }

    pub fn default_protein() -> Self {
        Self {
            moltype: "protein".to_string(),
            ksize: 10,
            scaled: 200,
            ..Self::default_dna()
        }
    }

    pub fn default_dayhoff() -> Self {
        Self {
            moltype: "dayhoff".to_string(),
            ksize: 10,
            scaled: 200,
            ..Self::default_dna()
        }
    }

    pub fn default_hp() -> Self {
        Self {
            moltype: "hp".to_string(),
            ksize: 10,
            scaled: 200,
            ..Self::default_dna()
        }
    }

    pub fn default_skipm1n3() -> Self {
        Self {
            moltype: "skipm1n3".to_string(),
            ksize: 21,
            scaled: 1000,
            ..Self::default_dna()
        }
    }

    pub fn default_skipm2n3() -> Self {
        Self {
            moltype: "skipm2n3".to_string(),
            ksize: 21,
            scaled: 1000,
            ..Self::default_dna()
        }
    }

    pub fn moltype(&self) -> HashFunctions {
        self.moltype.as_str().try_into().unwrap()
    }

    pub fn from_record(record: &Record) -> Self {
        Self {
            ksize: record.ksize(),
            moltype: record.moltype().to_string(),
            num: *record.num(),
            scaled: *record.scaled() as u32,
            with_abundance: record.with_abundance(),
            ..Self::default_dna() // ignore remaining fields
        }
    }

    pub fn matches_selection(&self, selection: &Selection) -> bool {
        let mut valid = true;

        if let Some(ksize) = selection.ksize() {
            valid = valid && self.ksize == ksize;
        }

        if let Some(moltype) = selection.moltype() {
            valid = valid && self.moltype() == moltype;
        }

        if let Some(abund) = selection.abund() {
            valid = valid && self.with_abundance == abund;
        }

        if let Some(scaled) = selection.scaled() {
            // num sigs have self.scaled = 0, don't include them
            valid = valid && self.scaled != 0 && self.scaled <= scaled as u32;
        }

        if let Some(num) = selection.num() {
            valid = valid && self.num == num;
        }

        valid
    }

    pub fn params(&self) -> (u32, String, bool, u32, u32) {
        (
            self.ksize,
            self.moltype.clone(),
            self.with_abundance,
            self.num,
            self.scaled,
        )
    }
}

impl PartialEq for BuildRecord {
    fn eq(&self, other: &Self) -> bool {
        self.ksize == other.ksize
            && self.moltype == other.moltype
            && self.with_abundance == other.with_abundance
            && self.num == other.num
            && self.scaled == other.scaled
    }
}

impl Eq for BuildRecord {}

impl Hash for BuildRecord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ksize.hash(state);
        self.moltype.hash(state);
        self.scaled.hash(state);
        self.num.hash(state);
        self.with_abundance.hash(state);
    }
}

#[derive(Debug, Default, Clone)]
pub struct BuildManifest {
    records: Vec<BuildRecord>,
}

impl BuildManifest {
    pub fn new() -> Self {
        BuildManifest {
            records: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn size(&self) -> usize {
        self.records.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = &BuildRecord> {
        self.records.iter()
    }

    // clear all records
    pub fn clear(&mut self) {
        self.records.clear();
    }

    pub fn summarize_params(&self) -> HashSet<(u32, String, bool, u32, u32)> {
        self.iter().map(|record| record.params()).collect()
    }

    pub fn filter_manifest(&self, other: &BuildManifest) -> Self {
        // Create a HashSet of references to the `BuildRecord`s in `other`
        let pairs: HashSet<_> = other.records.iter().collect();

        // Filter `self.records` to retain only those `BuildRecord`s that are NOT in `pairs`
        let records = self
            .records
            .iter()
            .filter(|&build_record| !pairs.contains(build_record))
            .cloned()
            .collect();

        Self { records }
    }

    pub fn add_record(&mut self, record: BuildRecord) {
        self.records.push(record);
    }

    pub fn extend_records(&mut self, other: impl IntoIterator<Item = BuildRecord>) {
        self.records.extend(other);
    }

    pub fn extend_from_manifest(&mut self, other: &BuildManifest) {
        self.records.extend(other.records.clone()); // Clone the records from the other manifest
    }

    pub fn to_writer<W: Write>(&self, mut wtr: W) -> Result<()> {
        // Write the manifest version as a comment
        wtr.write_all(b"# SOURMASH-MANIFEST-VERSION: 1.0\n")?;

        // Use CSV writer to serialize records
        let mut csv_writer = csv::Writer::from_writer(wtr);

        for record in &self.records {
            // don't write empty records (empty template sigs aren't written from BuildCollection)
            if record.sequence_added {
                csv_writer.serialize(record)?; // Serialize each BuildRecord
            }
        }

        csv_writer.flush()?; // Ensure all data is written

        Ok(())
    }

    pub fn write_manifest_to_zip<W: Write + Seek>(
        &self,
        zip: &mut ZipWriter<W>,
        options: &FileOptions<()>,
    ) -> Result<()> {
        zip.start_file("SOURMASH-MANIFEST.csv", *options)?;
        self.to_writer(zip)?;
        Ok(())
    }
}

impl MultiSelect for BuildManifest {
    fn select(&mut self, multi_selection: &MultiSelection) -> Result<(), SourmashError> {
        // Retain only the records that match any selection
        self.records.retain(|record| {
            multi_selection
                .selections
                .iter()
                .any(|selection| record.matches_selection(selection))
        });
        Ok(())
    }
}

impl From<Vec<BuildRecord>> for BuildManifest {
    fn from(records: Vec<BuildRecord>) -> Self {
        BuildManifest { records }
    }
}

impl Index<usize> for BuildManifest {
    type Output = BuildRecord;

    fn index(&self, index: usize) -> &Self::Output {
        &self.records[index]
    }
}

impl<'a> IntoIterator for &'a BuildManifest {
    type Item = &'a BuildRecord;
    type IntoIter = std::slice::Iter<'a, BuildRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.records.iter()
    }
}

impl<'a> IntoIterator for &'a mut BuildManifest {
    type Item = &'a mut BuildRecord;
    type IntoIter = std::slice::IterMut<'a, BuildRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.records.iter_mut()
    }
}

#[derive(Debug, Default, Clone)]
pub struct BuildCollection {
    pub manifest: BuildManifest,
    pub sigs: Vec<Signature>,
}

impl BuildCollection {
    pub fn new() -> Self {
        BuildCollection {
            manifest: BuildManifest::new(),
            sigs: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.manifest.is_empty()
    }

    pub fn size(&self) -> usize {
        self.manifest.size()
    }

    pub fn dna_size(&self) -> Result<usize, SourmashError> {
        let multiselection = MultiSelection::from_moltypes(vec!["dna"])?;
        let mut mf = self.manifest.clone(); // temporary mutable copy
        mf.select(&multiselection)?;
        Ok(mf.records.len())
    }

    pub fn skipm1n3_size(&self) -> Result<usize, SourmashError> {
        let multiselection = MultiSelection::from_moltypes(vec!["skipm1n3"])?;
        let mut mf = self.manifest.clone();
        mf.select(&multiselection)?;
        Ok(mf.records.len())
    }

    pub fn skipm2n3_size(&self) -> Result<usize, SourmashError> {
        let multiselection = MultiSelection::from_moltypes(vec!["skipm2n3"])?;
        let mut mf = self.manifest.clone();
        mf.select(&multiselection)?;
        Ok(mf.records.len())
    }

    pub fn anydna_size(&self) -> Result<usize, SourmashError> {
        let multiselection = MultiSelection::from_moltypes(vec!["DNA", "skipm1n3", "skipm2n3"])?;
        let mut mf = self.manifest.clone();
        mf.select(&multiselection)?;
        Ok(mf.records.len())
    }

    pub fn protein_size(&self) -> Result<usize, SourmashError> {
        let multiselection = MultiSelection::from_moltypes(vec!["protein"])?;
        let mut mf = self.manifest.clone();
        mf.select(&multiselection)?;
        Ok(mf.records.len())
    }

    pub fn anyprotein_size(&self) -> Result<usize, SourmashError> {
        let multiselection = MultiSelection::from_moltypes(vec!["protein", "dayhoff", "hp"])?;
        let mut mf = self.manifest.clone(); // temporary mutable copy
        mf.select(&multiselection)?;
        Ok(mf.records.len())
    }

    pub fn parse_ksize(value: &str) -> Result<u32, String> {
        value
            .parse::<u32>()
            .map_err(|_| format!("cannot parse k='{}' as a valid integer", value))
    }

    pub fn parse_int_once<T>(
        value: &str,
        field: &str,
        current: &mut Option<T>,
    ) -> Result<(), String>
    where
        T: FromStr<Err = ParseIntError> + Display + Copy,
    {
        let parsed_value = value
            .parse::<T>()
            .map_err(|_| format!("cannot parse {}='{}' as a valid integer", field, value))?;

        // Check for conflicts; we don't allow multiple values for the same field.
        if let Some(old_value) = *current {
            return Err(format!(
                "Conflicting values for '{}': {} and {}",
                field, old_value, parsed_value
            ));
        }

        *current = Some(parsed_value);
        Ok(())
    }

    pub fn parse_moltype(item: &str, current: &mut Option<String>) -> Result<String, String> {
        let new_moltype = match item {
            "protein" | "dna" | "dayhoff" | "hp" | "skipm1n3" | "skipm2n3" => item.to_string(),
            _ => return Err(format!("unknown moltype '{}'", item)),
        };

        // Check for conflicts and update the moltype.
        if let Some(existing) = current {
            if *existing != new_moltype {
                return Err(format!(
                    "Conflicting moltype settings in param string: '{}' and '{}'",
                    existing, new_moltype
                ));
            }
        }

        *current = Some(new_moltype.clone());
        Ok(new_moltype)
    }

    pub fn parse_abundance(item: &str, current: &mut Option<bool>) -> Result<(), String> {
        let new_abundance = item == "abund";

        if let Some(existing) = *current {
            if existing != new_abundance {
                return Err(format!(
                    "Conflicting abundance settings in param string: '{}'",
                    item
                ));
            }
        }

        *current = Some(new_abundance);
        Ok(())
    }

    pub fn summarize_params(&self) -> HashSet<(u32, String, bool, u32, u32)> {
        let params: HashSet<_> = self.manifest.iter().map(|record| record.params()).collect();

        // Print a description of the summary
        eprintln!("Building {} sketch types:", params.len());

        for (ksize, moltype, with_abundance, num, scaled) in &params {
            eprintln!(
                "    {},k={},scaled={},num={},abund={}",
                moltype, ksize, scaled, num, with_abundance
            );
        }
        params
    }

    pub fn parse_params(p_str: &str) -> Result<(BuildRecord, Vec<u32>), String> {
        let mut ksizes = Vec::new();
        let mut moltype: Option<String> = None;
        let mut track_abundance: Option<bool> = None;
        let mut num: Option<u32> = None;
        let mut scaled: Option<u32> = None;
        let mut seed: Option<u32> = None;

        for item in p_str.split(',') {
            match item {
                _ if item.starts_with("k=") => {
                    ksizes.push(Self::parse_ksize(&item[2..])?);
                }
                "abund" | "noabund" => {
                    Self::parse_abundance(item, &mut track_abundance)?;
                }
                "protein" | "dna" | "DNA" | "dayhoff" | "hp" | "skipm1n3" | "skipm2n3" => {
                    Self::parse_moltype(item, &mut moltype)?;
                }
                _ if item.starts_with("num=") => {
                    Self::parse_int_once(&item[4..], "num", &mut num)?;
                }
                _ if item.starts_with("scaled=") => {
                    Self::parse_int_once(&item[7..], "scaled", &mut scaled)?;
                }
                _ if item.starts_with("seed=") => {
                    Self::parse_int_once(&item[5..], "seed", &mut seed)?;
                }
                _ => {
                    return Err(format!(
                        "Error parsing params string '{}': Unknown component '{}'",
                        p_str, item
                    ));
                }
            }
        }

        // Ensure that moltype was set
        let moltype = moltype.ok_or_else(|| {
            format!(
                "Error parsing params string '{}': No moltype provided",
                p_str
            )
        })?;

        // Create a moltype-specific default BuildRecord or return an error if unsupported.
        let mut base_record = match moltype.as_str() {
            "dna" | "DNA" => BuildRecord::default_dna(),
            "protein" => BuildRecord::default_protein(),
            "dayhoff" => BuildRecord::default_dayhoff(),
            "hp" => BuildRecord::default_hp(),
            "skipm1n3" => BuildRecord::default_skipm1n3(),
            "skipm2n3" => BuildRecord::default_skipm2n3(),
            _ => {
                return Err(format!(
                    "Error parsing params string '{}': Unsupported moltype '{}'",
                    p_str, moltype
                ));
            }
        };

        // Apply parsed values
        if let Some(track_abund) = track_abundance {
            base_record.with_abundance = track_abund;
        }
        if let Some(n) = num {
            base_record.num = n;
        }
        if let Some(s) = scaled {
            base_record.scaled = s;
        }
        if let Some(s) = seed {
            base_record.seed = s;
        }

        // Use the default ksize if none were specified.
        if ksizes.is_empty() {
            ksizes.push(base_record.ksize);
        }

        // Ensure that num and scaled are mutually exclusive unless num is 0.
        if let (Some(n), Some(_)) = (num, scaled) {
            if n != 0 {
                return Err(format!(
                    "Error parsing params string '{}': Cannot specify both 'num' (non-zero) and 'scaled' in the same parameter string",
                    p_str
                ));
            }
        }

        Ok((base_record, ksizes))
    }

    pub fn from_param_str(params_str: &str) -> Result<Self, String> {
        if params_str.trim().is_empty() {
            return Err("Parameter string cannot be empty.".to_string());
        }

        let mut coll = BuildCollection::new();
        let mut seen_records = HashSet::new();

        for p_str in params_str.split('_') {
            // Use `parse_params` to get the base record and ksizes.
            let (base_record, ksizes) = Self::parse_params(p_str)?;

            // Iterate over each ksize and add a signature to the collection.
            for k in ksizes {
                let mut record = base_record.clone();
                record.ksize = k;

                // Check if the record is already in the set.
                if seen_records.insert(record.clone()) {
                    // Add the record and its associated signature to the collection.
                    coll.add_template_sig_from_record(&record);
                }
            }
        }
        Ok(coll)
    }

    pub fn from_manifest(manifest: &BuildManifest) -> Self {
        let mut collection = BuildCollection::new();

        // Iterate over each `BuildRecord` in the provided `BuildManifest`.
        for record in &manifest.records {
            // Add a signature to the collection using the `BuildRecord` and `input_moltype`.
            collection.add_template_sig_from_record(record);
        }

        collection
    }

    pub fn from_selection(selection: &Selection) -> Result<Self, String> {
        let mut collection = BuildCollection::new();

        // Set a default ksize if none is provided
        let ksizes = if let Some(ksize) = selection.ksize() {
            vec![ksize]
        } else {
            vec![21] // Default ksize
        };

        // Default moltype if not provided
        let moltype = selection
            .moltype()
            .clone()
            .ok_or("Moltype must be specified in selection")?;

        for ksize in ksizes {
            let mut record = match moltype {
                HashFunctions::Murmur64Dna => BuildRecord::default_dna(),
                HashFunctions::Murmur64Protein => BuildRecord::default_protein(),
                HashFunctions::Murmur64Dayhoff => BuildRecord::default_dayhoff(),
                HashFunctions::Murmur64Hp => BuildRecord::default_hp(),
                _ => {
                    return Err(format!("Unsupported moltype '{:?}' in selection", moltype));
                }
            };

            // Apply selection parameters to the BuildRecord
            record.ksize = ksize;
            if let Some(track_abundance) = selection.abund() {
                record.with_abundance = track_abundance;
            }
            if let Some(num) = selection.num() {
                record.num = num;
            }
            if let Some(scaled) = selection.scaled() {
                record.scaled = scaled;
            }

            // Add the template signature and record to the collection
            collection.add_template_sig_from_record(&record);
        }

        Ok(collection)
    }

    pub fn add_template_sig_from_record(&mut self, record: &BuildRecord) {
        // Adjust ksize for protein, dayhoff, or hp, which require tripling the k-mer size.
        let adjusted_ksize = match record.moltype.as_str() {
            "protein" | "dayhoff" | "hp" => record.ksize * 3,
            _ => record.ksize,
        };

        // Construct ComputeParameters.
        let cp = ComputeParameters::builder()
            .ksizes(vec![adjusted_ksize])
            .scaled(record.scaled as u32)
            .protein(record.moltype == "protein")
            .dna(record.moltype == "DNA")
            .dayhoff(record.moltype == "dayhoff")
            .hp(record.moltype == "hp")
            .skipm1n3(record.moltype == "skipm1n3")
            .skipm2n3(record.moltype == "skipm2n3")
            .num_hashes(record.num)
            .track_abundance(record.with_abundance)
            .build();

        // Create a Signature from the ComputeParameters.
        let sig = Signature::from_params(&cp);

        // Clone the `BuildRecord` and use it directly.
        let template_record = record.clone();

        // Add the record and signature to the collection.
        self.manifest.records.push(template_record);
        self.sigs.push(sig);
    }

    pub fn filter_manifest(&mut self, other: &BuildManifest) {
        self.manifest = self.manifest.filter_manifest(other)
    }

    pub fn filter_by_manifest(&mut self, other: &BuildManifest) {
        // Create a HashSet for efficient filtering based on the `BuildRecord`s in `other`.
        let other_records: HashSet<_> = other.records.iter().collect();

        // Retain only the records that are not in `other_records`, filtering in place.
        let mut sig_index = 0;
        self.manifest.records.retain(|record| {
            let keep = !other_records.contains(record);
            if !keep {
                // Remove the corresponding signature at the same index.
                self.sigs.remove(sig_index);
            } else {
                sig_index += 1; // Only increment if we keep the record and signature.
            }
            keep
        });
    }

    pub fn filter(&mut self, params_set: &HashSet<u64>) {
        let mut index = 0;
        while index < self.manifest.records.len() {
            let record = &self.manifest.records[index];

            // filter records with matching Params
            if params_set.contains(&record.hashed_params) {
                self.manifest.records.remove(index);
                self.sigs.remove(index);
            } else {
                index += 1;
            }
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = (Idx, &BuildRecord)> {
        self.manifest.iter().enumerate().map(|(i, r)| (i as Idx, r))
    }

    pub fn record_for_dataset(&self, dataset_id: Idx) -> Result<&BuildRecord> {
        Ok(&self.manifest[dataset_id as usize])
    }

    pub fn sigs_iter_mut(&mut self) -> impl Iterator<Item = &mut Signature> {
        self.sigs.iter_mut()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = (&mut BuildRecord, &mut Signature)> {
        // zip together mutable iterators over records and sigs
        self.manifest.records.iter_mut().zip(self.sigs.iter_mut())
    }

    fn build_sigs_from_record(
        &mut self,
        input_moltype: &str,
        record: &SequenceRecord,
    ) -> Result<()> {
        // Optionally use `par_iter_mut` for parallel execution
        self.iter_mut().try_for_each(|(rec, sig)| {
            if input_moltype == "protein"
                && (rec.moltype() == HashFunctions::Murmur64Protein
                    || rec.moltype() == HashFunctions::Murmur64Dayhoff
                    || rec.moltype() == HashFunctions::Murmur64Hp)
            {
                sig.add_protein(&record.seq())
                    .context("Failed to add protein")?;
                if !rec.sequence_added {
                    rec.sequence_added = true;
                }
            } else if (input_moltype == "DNA" || input_moltype == "dna")
                && (rec.moltype() == HashFunctions::Murmur64Dna
                    || rec.moltype() == HashFunctions::Murmur64Skipm1n3
                    || rec.moltype() == HashFunctions::Murmur64Skipm2n3)
            {
                sig.add_sequence(&record.seq(), true)
                    .context("Failed to add sequence")?;
                if !rec.sequence_added {
                    rec.sequence_added = true;
                }
            }
            Ok(())
        })
    }

    pub fn build_sigs_from_data(
        &mut self,
        data: Vec<u8>,
        input_moltype: &str,
        name: String,
        filename: String,
    ) -> Result<()> {
        let cursor = Cursor::new(data);
        let mut fastx_reader =
            parse_fastx_reader(cursor).context("Failed to parse FASTA/FASTQ data")?;

        // Iterate over FASTA records and add sequences/proteins to sigs
        while let Some(record) = fastx_reader.next() {
            let record = record.context("Failed to read record")?;
            self.build_sigs_from_record(input_moltype, &record)?;
        }

        // After processing sequences, update sig, record information
        self.update_info(name, filename);

        Ok(())
    }

    pub fn build_sigs_from_file_or_stdin(
        &mut self,
        input_moltype: &str, // "protein" or "DNA"
        name: String,
        filename: String,
    ) -> Result<u64> {
        // Create a FASTX reader from the file or stdin
        let mut fastx_reader = if filename == "-" {
            parse_fastx_stdin().context("Failed to parse FASTA/FASTQ data from stdin")?
        } else {
            parse_fastx_file(&filename).context("Failed to open file for FASTA/FASTQ data")?
        };

        // Counter for the number of records processed
        let mut record_count: u64 = 0;

        // Parse records and add sequences to signatures
        while let Some(record_result) = fastx_reader.next() {
            let record = record_result.context("Failed to read a record from input")?;

            self.build_sigs_from_record(input_moltype, &record)?;

            record_count += 1;
        }

        // Update signature and record metadata
        self.update_info(name, filename);

        // Return the count of records parsed
        Ok(record_count)
    }

    pub fn build_singleton_sigs(
        &mut self,
        record: SequenceRecord,
        input_moltype: &str, // (protein/dna); todo - use hashfns?
        filename: String,
    ) -> Result<()> {
        self.build_sigs_from_record(input_moltype, &record)?;
        // After processing sequences, update sig, record information
        let record_name = std::str::from_utf8(record.id())
            .expect("could not get record id")
            .to_string();
        self.update_info(record_name, filename);

        Ok(())
    }

    pub fn update_info(&mut self, name: String, filename: String) {
        // update the records to reflect information the signature;
        for (record, sig) in self.iter_mut() {
            if record.sequence_added {
                // update signature name, filename
                sig.set_name(name.as_str());
                sig.set_filename(filename.as_str());

                // update record: set name, filename, md5sum, n_hashes
                record.set_name(Some(name.clone()));
                record.set_filename(Some(filename.clone()));
                record.set_md5(Some(sig.md5sum()));
                record.set_md5short(Some(sig.md5sum()[0..8].into()));
                record.set_n_hashes(Some(
                    sig.get_sketch().expect("cannot retrieve sketch").size(),
                ));

                // note, this needs to be set when writing sigs (not here)
                // record.set_internal_location("")
            }
        }
    }

    pub fn write_sigs(&mut self, output: &str) -> Result<()> {
        let gzip = output.ends_with(".gz");
        if output == "-" {
            // Write to stdout
            let stdout = std::io::stdout();
            let mut handle = stdout.lock();
            self.write_sigs_as_json(&mut handle, gzip)
                .context("Failed to write signatures to stdout")?;
            handle.flush().context("Failed to flush stdout")?;
        } else if output.ends_with(".zip") {
            let options = FileOptions::default()
                .compression_method(CompressionMethod::Stored)
                .unix_permissions(0o644)
                .large_file(true);
            // Write to a zip file
            let file =
                File::create(output).context(format!("Failed to create file: {}", output))?;
            let mut zip = ZipWriter::new(file);
            let mut md5sum_occurrences: HashMap<String, usize> = HashMap::new();
            self.write_sigs_to_zip(&mut zip, &mut md5sum_occurrences, &options)
                .context(format!(
                    "Failed to write signatures to zip file: {}",
                    output
                ))?;
            println!("Writing manifest");
            self.manifest.write_manifest_to_zip(&mut zip, &options)?;
            zip.finish()?;
        } else {
            // Write JSON to output file
            let file =
                File::create(output).context(format!("Failed to create file: {}", output))?;
            let mut writer = std::io::BufWriter::new(file);
            self.write_sigs_as_json(&mut writer, gzip)
                .context(format!("Failed to write signatures to file: {}", output))?;
        }
        Ok(())
    }

    pub fn write_sigs_to_zip<W: Write + Seek>(
        &mut self, // need mutable to update records
        zip: &mut ZipWriter<W>,
        md5sum_occurrences: &mut HashMap<String, usize>,
        options: &FileOptions<()>,
    ) -> Result<()> {
        // iterate over both records and signatures
        for (record, sig) in self.iter_mut() {
            // skip any empty sig templates (no sequence added)
            // TODO --> test that this is working
            if !record.sequence_added {
                continue;
            }
            let md5sum_str = sig.md5sum();
            let count = md5sum_occurrences.entry(md5sum_str.clone()).or_insert(0);
            *count += 1;

            // Generate the signature filename
            let sig_filename = if *count > 1 {
                format!("signatures/{}_{}.sig.gz", md5sum_str, count)
            } else {
                format!("signatures/{}.sig.gz", md5sum_str)
            };

            // Update record's internal_location with the signature filename
            record.internal_location = Some(sig_filename.clone().into());

            // Serialize signature to JSON
            let wrapped_sig = vec![sig.clone()];
            let json_bytes = serde_json::to_vec(&wrapped_sig)
                .map_err(|e| anyhow!("Error serializing signature: {}", e))?;

            // Gzip compress the JSON bytes
            let gzipped_buffer = {
                let mut buffer = Cursor::new(Vec::new());
                {
                    let mut gz_writer = niffler::get_writer(
                        Box::new(&mut buffer),
                        niffler::compression::Format::Gzip,
                        niffler::compression::Level::Nine,
                    )?;
                    gz_writer.write_all(&json_bytes)?;
                }
                buffer.into_inner()
            };

            zip.start_file(sig_filename, *options)?;
            zip.write_all(&gzipped_buffer)
                .map_err(|e| anyhow!("Error writing zip entry for signature: {}", e))?;
        }

        Ok(())
    }

    pub fn write_sigs_as_json<W: Write>(
        &mut self, // mutable to update records if needed
        writer: &mut W,
        gzip: bool,
    ) -> Result<()> {
        // Create a vector to store all signatures
        let mut all_signatures = Vec::new();

        // Iterate over both records and signatures
        for (record, sig) in self.iter_mut() {
            // Skip any empty sig templates (no sequence added)
            if !record.sequence_added {
                continue;
            }

            // Add the signature to the collection for JSON serialization
            all_signatures.push(sig.clone());
        }

        // Serialize all signatures to JSON
        let json_bytes = serde_json::to_vec(&all_signatures)
            .map_err(|e| anyhow!("Error serializing signatures to JSON: {}", e))?;

        if gzip {
            // Gzip compress the JSON bytes
            let mut gz_writer = niffler::get_writer(
                Box::new(writer),
                niffler::compression::Format::Gzip,
                niffler::compression::Level::Nine,
            )?;
            gz_writer.write_all(&json_bytes)?;
        } else {
            // Write uncompressed JSON to the writer
            writer.write_all(&json_bytes)?;
        }

        Ok(())
    }
}

impl<'a> IntoIterator for &'a mut BuildCollection {
    type Item = (&'a mut BuildRecord, &'a mut Signature);
    type IntoIter =
        std::iter::Zip<std::slice::IterMut<'a, BuildRecord>, std::slice::IterMut<'a, Signature>>;

    fn into_iter(self) -> Self::IntoIter {
        self.manifest.records.iter_mut().zip(self.sigs.iter_mut())
    }
}

impl MultiSelect for BuildCollection {
    // in sourmash core, we don't need to select sigs themselves. Is this due to the way that Idx/Storage work?
    fn select(&mut self, multi_selection: &MultiSelection) -> Result<(), SourmashError> {
        // Retain records and sigs in place
        let mut i = 0;
        self.manifest.records.retain(|record| {
            let keep = multi_selection
                .selections
                .iter()
                .any(|selection| record.matches_selection(selection));

            if !keep {
                self.sigs.remove(i); // Remove corresponding signature
            } else {
                i += 1;
            }
            keep
        });

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_params_str() {
        let params_str = "k=31,abund,dna";
        let result = BuildCollection::parse_params(params_str);

        assert!(
            result.is_ok(),
            "Expected 'k=31,abund,dna' to be valid, but got an error: {:?}",
            result
        );

        let (record, ksizes) = result.unwrap();

        // Verify that the Record, ksizes have the correct settings.
        assert_eq!(record.moltype, "DNA");
        assert_eq!(record.with_abundance, true);
        assert_eq!(ksizes, vec![31]);
        assert_eq!(record.scaled, 1000, "Expected default scaled value of 1000");
        assert_eq!(record.num, 0, "Expected default num value of 0");
    }

    #[test]
    fn test_from_param_str() {
        let params_str = "k=31,abund,dna_dna,k=21,k=31,k=51,abund_k=10,protein";
        let coll_result = BuildCollection::from_param_str(params_str);

        assert!(
            coll_result.is_ok(),
            "Param str '{}' is valid, but got an error: {:?}",
            params_str,
            coll_result
        );

        let coll = coll_result.unwrap();

        // Ensure the BuildCollection contains the expected number of records.
        // Note that "k=31,abund,dna" appears in two different parameter strings, so it should only appear once.
        assert_eq!(
            coll.manifest.records.len(),
            4,
            "Expected 4 unique BuildRecords in the collection, but found {}",
            coll.manifest.records.len()
        );

        // Define the expected BuildRecords for comparison.
        let expected_records = vec![
            BuildRecord {
                ksize: 31,
                moltype: "DNA".to_string(),
                with_abundance: true,
                ..BuildRecord::default_dna()
            },
            BuildRecord {
                ksize: 21,
                moltype: "DNA".to_string(),
                with_abundance: true,
                ..BuildRecord::default_dna()
            },
            BuildRecord {
                ksize: 51,
                moltype: "DNA".to_string(),
                with_abundance: true,
                ..BuildRecord::default_dna()
            },
            BuildRecord::default_protein(),
        ];

        // Verify that each expected BuildRecord is present in the collection.
        for expected_record in expected_records {
            assert!(
                coll.manifest.records.contains(&expected_record),
                "Expected BuildRecord with ksize: {}, moltype: {}, with_abundance: {} not found in the collection",
                expected_record.ksize,
                expected_record.moltype,
                expected_record.with_abundance
            );
        }

        // Optionally, check that the corresponding signatures are present.
        assert_eq!(
            coll.sigs.len(),
            4,
            "Expected 4 Signatures in the collection, but found {}",
            coll.sigs.len()
        );
    }

    #[test]
    fn test_invalid_params_str_conflicting_moltypes() {
        let params_str = "k=31,abund,dna,protein";
        let result = BuildCollection::from_param_str(params_str);

        assert!(
            result.is_err(),
            "Expected 'k=31,abund,dna,protein' to be invalid due to conflicting moltypes, but got a successful result"
        );

        // Check if the error message contains the expected conflict text.
        if let Err(e) = result {
            assert!(
                e.contains("Conflicting moltype settings"),
                "Expected error to contain 'Conflicting moltype settings', but got: {}",
                e
            );
        }
    }

    #[test]
    fn test_unknown_component_error() {
        // Test for an unknown component that should trigger an error.
        let result = BuildCollection::from_param_str("dna,k=31,notaparam");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Error parsing params string 'dna,k=31,notaparam': Unknown component 'notaparam'"
        );
    }

    #[test]
    fn test_unknown_component_error2() {
        // Test a common param string error (k=31,51 compared with valid k=31,k=51)
        let result = BuildCollection::from_param_str("dna,k=31,51,abund");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Error parsing params string 'dna,k=31,51,abund': Unknown component '51'"
        );
    }

    #[test]
    fn test_conflicting_num_and_scaled() {
        // Test for specifying both num and scaled, which should result in an error.
        let result = BuildCollection::from_param_str("dna,k=31,num=10,scaled=1000");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Error parsing params string 'dna,k=31,num=10,scaled=1000': Cannot specify both 'num' (non-zero) and 'scaled' in the same parameter string"
        );
    }

    #[test]
    fn test_conflicting_abundance() {
        // Test for providing conflicting abundance settings, which should result in an error.
        let result = BuildCollection::from_param_str("dna,k=31,abund,noabund");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Conflicting abundance settings in param string: 'noabund'"
        );
    }

    #[test]
    fn test_invalid_ksize_format() {
        // Test for an invalid ksize format that should trigger an error.
        let result = BuildCollection::from_param_str("dna,k=abc");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "cannot parse k='abc' as a valid integer"
        );
    }

    #[test]
    fn test_invalid_num_format() {
        // Test for an invalid number format that should trigger an error.
        let result = BuildCollection::from_param_str("dna,k=31,num=abc");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "cannot parse num='abc' as a valid integer"
        );
    }

    #[test]
    fn test_invalid_scaled_format() {
        // Test for an invalid scaled format that should trigger an error.
        let result = BuildCollection::from_param_str("dna,k=31,scaled=abc");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "cannot parse scaled='abc' as a valid integer"
        );
    }

    #[test]
    fn test_invalid_seed_format() {
        // Test for an invalid seed format that should trigger an error.
        let result = BuildCollection::from_param_str("dna,k=31,seed=abc");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "cannot parse seed='abc' as a valid integer"
        );
    }

    #[test]
    fn test_repeated_values() {
        // repeated scaled
        let result = BuildCollection::from_param_str("dna,k=31,scaled=1,scaled=1000");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Conflicting values for 'scaled': 1 and 1000"
        );

        // repeated num
        let result = BuildCollection::from_param_str("dna,k=31,num=1,num=1000");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Conflicting values for 'num': 1 and 1000"
        );

        // repeated seed
        let result = BuildCollection::from_param_str("dna,k=31,seed=1,seed=42");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Conflicting values for 'seed': 1 and 42"
        );
    }

    #[test]
    fn test_missing_ksize() {
        // Test for a missing ksize, using default should not result in an error.
        let result = BuildCollection::from_param_str("dna,abund");
        assert!(result.is_ok(), "Expected Ok but got an error.");
    }

    #[test]
    fn test_repeated_ksize() {
        // Repeated ksize settings should not trigger an error since it is valid to have multiple ksizes.
        let result = BuildCollection::from_param_str("dna,k=31,k=21");
        assert!(result.is_ok(), "Expected Ok but got an error.");
    }

    #[test]
    fn test_empty_string() {
        // Test for an empty parameter string, which should now result in an error.
        let result = BuildCollection::from_param_str("");
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(result.unwrap_err(), "Parameter string cannot be empty.");
    }

    #[test]
    fn test_filter_by_manifest_with_matching_records() {
        // Create a BuildCollection with some records and signatures.

        let rec1 = BuildRecord::default_dna();
        let rec2 = BuildRecord {
            ksize: 21,
            moltype: "DNA".to_string(),
            scaled: 1000,
            ..BuildRecord::default_dna()
        };
        let rec3 = BuildRecord {
            ksize: 31,
            moltype: "DNA".to_string(),
            scaled: 1000,
            with_abundance: true,
            ..BuildRecord::default_dna()
        };

        let bmanifest = BuildManifest {
            records: vec![rec1.clone(), rec2.clone(), rec3.clone()],
        };
        // let mut dna_build_collection = BuildCollection::from_manifest(&bmanifest, "DNA");
        let mut dna_build_collection = BuildCollection::from_manifest(&bmanifest);

        // Create a BuildManifest with records to filter out.
        let filter_manifest = BuildManifest {
            records: vec![rec1],
        };

        // Apply the filter.
        dna_build_collection.filter_by_manifest(&filter_manifest);

        // check that the default DNA sig remains
        assert_eq!(dna_build_collection.manifest.size(), 2);

        let remaining_records = &dna_build_collection.manifest.records;

        assert!(remaining_records.contains(&rec2));
        assert!(remaining_records.contains(&rec3));
    }

    #[test]
    fn test_add_template_sig_from_record() {
        // Create a BuildCollection.
        let mut build_collection = BuildCollection::new();

        // Create a DNA BuildRecord.
        let dna_record = BuildRecord {
            ksize: 31,
            moltype: "DNA".to_string(),
            scaled: 1000,
            with_abundance: true,
            ..BuildRecord::default_dna()
        };

        // Add the DNA record to the collection
        build_collection.add_template_sig_from_record(&dna_record);

        // Verify that the record was added.
        assert_eq!(build_collection.manifest.records.len(), 1);
        assert_eq!(build_collection.sigs.len(), 1);

        let added_record = &build_collection.manifest.records[0];
        assert_eq!(added_record.moltype, "DNA");
        assert_eq!(added_record.ksize, 31);
        assert_eq!(added_record.with_abundance, true);

        // Create a protein BuildRecord.
        let protein_record = BuildRecord {
            ksize: 10,
            moltype: "protein".to_string(),
            scaled: 200,
            with_abundance: false,
            ..BuildRecord::default_dna()
        };

        // Add the protein record to the collection
        build_collection.add_template_sig_from_record(&protein_record);

        // Verify that the protein record was added.
        assert_eq!(build_collection.manifest.records.len(), 2);
        assert_eq!(build_collection.sigs.len(), 2);

        let added_protein_record = &build_collection.manifest.records[1];
        assert_eq!(added_protein_record.moltype, "protein");
        assert_eq!(added_protein_record.ksize, 10);
        assert_eq!(added_protein_record.with_abundance, false);

        // Create a BuildRecord with a non-matching moltype.
        let dayhoff_record = BuildRecord {
            ksize: 10,
            moltype: "dayhoff".to_string(),
            scaled: 200,
            with_abundance: true,
            ..BuildRecord::default_dna()
        };

        // Add dayhoff record.
        build_collection.add_template_sig_from_record(&dayhoff_record);

        // Verify that the record was added.
        assert_eq!(build_collection.manifest.records.len(), 3);
        assert_eq!(build_collection.sigs.len(), 3);

        let added_dayhoff_record = &build_collection.manifest.records[2];
        assert_eq!(added_dayhoff_record.moltype, "dayhoff");
        assert_eq!(added_dayhoff_record.ksize, 10);
        assert_eq!(added_dayhoff_record.with_abundance, true);
    }

    #[test]
    fn test_from_selection_dna_with_defaults() {
        // Create a selection with DNA moltype and default parameters
        let selection = Selection::builder()
            .moltype(HashFunctions::Murmur64Dna)
            .build();

        // Call from_selection
        let build_collection = BuildCollection::from_selection(&selection)
            .expect("Failed to create BuildCollection from selection");

        // Validate that the collection is not empty
        assert!(
            !build_collection.is_empty(),
            "BuildCollection should not be empty"
        );

        // Validate that the manifest contains the correct record
        assert_eq!(
            build_collection.manifest.size(),
            1,
            "Expected one record in the manifest"
        );

        let record = &build_collection.manifest.records[0];
        assert_eq!(record.moltype, "dna", "Expected moltype to be 'dna'");
        assert_eq!(record.ksize, 21, "Expected default ksize to be 21");
        assert!(
            !record.with_abundance,
            "Expected default abundance to be false"
        );
    }

    #[test]
    fn test_from_selection_with_custom_parameters() {
        // Create a selection with custom parameters
        let selection = Selection::builder()
            .moltype(HashFunctions::Murmur64Protein)
            .ksize(31)
            .abund(true)
            .scaled(1000)
            .build();

        // Call from_selection
        let build_collection = BuildCollection::from_selection(&selection)
            .expect("Failed to create BuildCollection from selection");

        // Validate that the collection is not empty
        assert!(
            !build_collection.is_empty(),
            "BuildCollection should not be empty"
        );

        // Validate that the manifest contains the correct record
        assert_eq!(
            build_collection.manifest.size(),
            1,
            "Expected one record in the manifest"
        );

        let record = &build_collection.manifest.records[0];
        assert_eq!(
            record.moltype, "protein",
            "Expected moltype to be 'protein'"
        );
        assert_eq!(record.ksize, 31, "Expected ksize to be 31");
        assert!(record.with_abundance, "Expected abundance to be true");
        assert_eq!(record.scaled, 1000, "Expected scaled to be 1000");
    }

    #[test]
    fn test_from_selection_multiple_ksizes() {
        // Create a selection with multiple ksizes
        let selection = Selection::builder()
            .moltype(HashFunctions::Murmur64Dayhoff)
            .ksize(21) // Simulate multiple ksizes by changing test logic
            .build();

        // Call from_selection
        let build_collection = BuildCollection::from_selection(&selection)
            .expect("Failed to create BuildCollection from selection");

        // Validate that the collection contains the correct number of records
        assert!(
            !build_collection.is_empty(),
            "BuildCollection should not be empty"
        );

        assert_eq!(
            build_collection.manifest.size(),
            1,
            "Expected one record in the manifest"
        );

        let record = &build_collection.manifest.records[0];
        assert_eq!(
            record.moltype, "dayhoff",
            "Expected moltype to be 'dayhoff'"
        );
        assert_eq!(record.ksize, 21, "Expected ksize to be 21");
    }

    #[test]
    fn test_from_selection_missing_moltype() {
        // Create a selection without a moltype
        let selection = Selection::builder().ksize(31).build();

        // Call from_selection and expect an error
        let result = BuildCollection::from_selection(&selection);
        assert!(result.is_err(), "Expected an error due to missing moltype");
        assert_eq!(
            result.unwrap_err(),
            "Moltype must be specified in selection",
            "Unexpected error message"
        );
    }
}
