//! Identities for worker requests whose sources span several scientific values.
use collagene_model::{Error, Molecule};
pub(super) fn profile_binding(
    revisions: &[String],
    reverse: &[bool],
    overlaps: &[usize],
    circular: bool,
    options: &crate::overlap_profiles::ProfileOptions,
) -> Result<String, Error> {
    if revisions.is_empty()
        || revisions.len() > 1000
        || reverse.len() != revisions.len()
        || overlaps.len() > 1000
    {
        return Err(Error::Format(
            "profile binding needs explicit bounded source/orientation lists".into(),
        ));
    }
    let bytes = serde_json::to_vec(&(revisions, reverse, overlaps, circular, options))
        .map_err(|e| Error::Format(e.to_string()))?;
    Ok(blake3::hash(&bytes).to_hex().to_string())
}
pub fn overlap_profile_revision(
    molecules: &[Molecule],
    reverse: &[bool],
    overlaps: &[usize],
    circular: bool,
    options: &crate::overlap_profiles::ProfileOptions,
) -> Result<String, Error> {
    if molecules.is_empty() || molecules.len() > 1000 {
        return Err(Error::TooLarge);
    }
    profile_binding(
        &molecules.iter().map(Molecule::revision).collect::<Vec<_>>(),
        reverse,
        overlaps,
        circular,
        options,
    )
}

/// Exact immutable alignment content identity, including original source evidence.
pub fn alignment_revision(
    document: &crate::alignment_formats::AlignmentDocument,
) -> Result<String, Error> {
    Ok(
        blake3::hash(&crate::alignment_formats::save_alignment(document)?)
            .to_hex()
            .to_string(),
    )
}
