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

data! {pub struct Anchors {
    pub source_revision:String,
    /// Zero-based source indices of provided framework boundaries 26/39 or55/66.
    pub left_index:usize,pub right_index:usize,pub provenance:String,
}}
data! {pub struct ConditionalLoop {
    pub engine:String,pub source_revision:String,pub rule_revision:String,pub context_revision:String,
    pub kind:Kind,pub anchors:Anchors,pub start:usize,pub end:usize,pub residues:Vec<NumberedResidue>,
}}
/// Conditional on independently established boundary identities. Unlike CDR3,
/// these framework anchors do not have invariant amino-acid identities to check.
pub fn number_anchored_loop(
    source: &Molecule,
    kind: &Kind,
    anchors: &Anchors,
) -> Result<ConditionalLoop> {
    let revision = source_revision(source)?;
    number_validated_loop(source, kind, anchors, &revision)
}
// Caller has already checked this immutable source and its exact revision.
pub(crate) fn number_validated_loop(
    source: &Molecule,
    kind: &Kind,
    anchors: &Anchors,
    revision: &str,
) -> Result<ConditionalLoop> {
    if revision != anchors.source_revision {
        return Err(Error::Binding);
    }
    if anchors.provenance.trim().is_empty()
        || anchors.provenance.len() > 1024
        || anchors.provenance.chars().any(char::is_control)
    {
        return Err(Error::Invalid("loop boundary provenance"));
    }
    let bytes = source.sequence.as_str().as_bytes();
    let start = anchors
        .left_index
        .checked_add(1)
        .ok_or(Error::Invalid("loop boundary index"))?;
    let end = anchors.right_index;
    if start > end || end >= bytes.len() {
        return Err(Error::Invalid("loop boundary order/bounds"));
    }
    if bytes[anchors.left_index..=end].contains(&b'*') {
        return Err(Error::Invalid("stop in conditional loop or boundary"));
    }
    let labels = positions(kind.clone(), end - start)?;
    let rule_revision = catalog_revision();
    let context_revision = hash(
        b"imgt-cdr12-context-1\0",
        &(LOOP_ENGINE, revision, &rule_revision, kind, anchors),
        65536,
    )?;
    let residues = labels
        .into_iter()
        .enumerate()
        .map(|(i, position)| NumberedResidue {
            source_index: start + i,
            residue: bytes[start + i] as char,
            position,
        })
        .collect();
    Ok(ConditionalLoop {
        engine: LOOP_ENGINE.into(),
        source_revision: revision.into(),
        rule_revision,
        context_revision,
        kind: kind.clone(),
        anchors: anchors.clone(),
        start,
        end,
        residues,
    })
}
impl ConditionalLoop {
    pub fn validate_for(&self, source: &Molecule, kind: &Kind, anchors: &Anchors) -> Result<()> {
        if self.engine != LOOP_ENGINE
            || self.residues.len() > 12
            || *self != number_anchored_loop(source, kind, anchors)?
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
