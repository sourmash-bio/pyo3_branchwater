use crate::utils::is_revindex_database;
use anyhow::Result;

use sourmash::index::revindex::{RevIndex, RevIndexOps};

pub fn check(index: camino::Utf8PathBuf, quick: bool, rw: bool) -> Result<()> {
    if !is_revindex_database(&index) {
        bail!("'{}' is not a valid RevIndex database", index);
    }

    println!("Opening DB (rw mode? {})", rw);
    let db = match RevIndex::open(index, !rw, None) {
        Ok(db) => db,
        Err(e) => {
            return Err(anyhow::anyhow!(
                "cannot open RocksDB database. Error is: {}",
                e
            ))
        }
    };

    println!("Starting check");
    db.check(quick);

    println!("Finished check");
    Ok(())
}
