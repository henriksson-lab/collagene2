//! Published V-domain display convention for CDR1/2, conditional on independently
//! established loop boundaries. This is not germline storage-label conversion.
use super::*;
pub const LOOP_ENGINE: &str = "imgt-cdr12-central-gaps-1";
data! {pub enum Kind {Cdr1,Cdr2}}
pub fn positions(kind: Kind, length: usize) -> Result<Vec<Position>> {
    let (first, last, minimum, maximum) = match kind {
        Kind::Cdr1 => (27_u16, 38_u16, 5, 12),
        Kind::Cdr2 => (56, 65, 0, 10),
    };
    if length < minimum || length > maximum {
        return Err(Error::Invalid("unsupported IMGT CDR1/2 length"));
    }
    let left = length.div_ceil(2);
    let right = length / 2;
    Ok((0..length)
        .map(|i| {
            position(if i < left {
                first + i as u16
            } else {
                last + 1 - right as u16 + (i - left) as u16
            })
        })
        .collect())
}
pub fn catalog_revision() -> String {
    hash(b"imgt-cdr12-central-gaps-1\0",&(LOOP_ENGINE,
        "https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html",
        "CDR1 27..38 length5..12; CDR2 56..65 length0..10; ceil half left/floor half right; source order; no insertion labels; conditional loop boundaries only"),65536).expect("fixed loop rule catalog")
}
