//! Explicit non-destructive annotation decisions and documentary contributions.
//! Attribution is source evidence, never biological verification or event authority.
mod assessment;
mod contract;
mod model;
mod provenance;
mod sources;
use super::library;
use crate::encoding::bounded_json;
use collagene_model::{Error, Feature, Molecule};
pub use model::*;
pub use provenance::{Contribution, Ledger, SourcePin};
pub const ENGINE_VERSION: &str = "annotation-merge-1";
pub const MAX_REQUEST_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_SELECTIONS: usize = 1000;
pub const MAX_LEDGER_BYTES: usize = 64 * 1024;
pub const MAX_CONTRIBUTIONS: usize = 32;
pub const QUALIFIER: &str = "annotation_merge:contributions";
fn bad(message: &str) -> Error {
    Error::Format(message.into())
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
}
fn text(value: &str, maximum: usize) -> bool {
    !value.trim().is_empty() && value.len() <= maximum && !value.chars().any(char::is_control)
}
pub fn revision(target: &Molecule, request: &Request) -> Result<String, Error> {
    request.validate(target)?;
    bounded_json::hash(&(ENGINE_VERSION, request))
}
pub fn preview(target: &Molecule, request: &Request) -> Result<Preview, Error> {
    assessment::preview(target, request)
}
pub fn apply(target: &Molecule, request: &Request) -> Result<Applied, Error> {
    assessment::apply(target, request).map(|(applied, _)| applied)
}
pub(crate) fn apply_with_projection(
    target: &Molecule,
    request: &Request,
) -> Result<(Applied, crate::sequence::editing::AnnotationIdentity), Error> {
    assessment::apply(target, request)
}
