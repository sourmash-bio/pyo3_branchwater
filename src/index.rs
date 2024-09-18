use sourmash::index::revindex::RevIndex;
use sourmash::index::revindex::RevIndexOps;
use sourmash::prelude::*;
use std::path::Path;

use crate::utils::{load_collection, ReportType};

pub fn index<P: AsRef<Path>>(
    siglist: String,
    selection: &Selection,
    output: P,
    colors: bool,
    allow_failed_sigpaths: bool,
    use_internal_storage: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Loading siglist");
    let allow_empty_collection = false;

    let collection = load_collection(
        &siglist,
        selection,
        ReportType::General,
        allow_failed_sigpaths,
        allow_empty_collection,
    )?;

    let mut index = RevIndex::create(
        output.as_ref(),
        collection.select(selection)?.try_into()?,
        colors,
    )?;

    if use_internal_storage {
        index.internalize_storage()?;
    }

    Ok(())
}
