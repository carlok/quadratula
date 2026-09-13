//! quadratula: the refutation floor of Schröder's 990 quasigroup laws.
//!
//! Exhaustive enumeration of small quasigroups up to isomorphism, their law
//! signatures, and the parastrophe (S3) action on tables and on laws.

pub mod enumerate;
pub mod law;
pub mod parastrophe;
pub mod qg;
pub mod search;
pub mod sig;

pub use sig::Signature;

/// The pinned law file, resolved against the crate root.
pub const LAWS_PATH: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/data/lefloch/quasigroup-equations.txt"
);

pub fn load_laws() -> Vec<law::Law> {
    let text = std::fs::read_to_string(LAWS_PATH).expect("read pinned law file");
    law::parse_laws(&text).expect("parse pinned law file")
}
