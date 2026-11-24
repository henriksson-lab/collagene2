//! Exact-source digest derivation using the existing native workspace/history model.
//! Installation reconstructs the product afresh. Imported recipe metadata is evidence only.
pub use crate::run_lineage::EventBinding as SourceBinding;
use crate::{
    Error,
    analysis::Enzyme,
    history::DocumentHistory,
    restriction_inspection::{self, InspectionOptions},
    run_lineage,
    workspace::{self, Workspace},
};
use serde::{Deserialize, Serialize};
macro_rules! dto{($($x:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
dto! {pub struct InstalledDigest{pub archive:Vec<u8>,pub source:SourceBinding,pub document_id:String,pub event_id:String,pub product_revision:String,pub batch_revision:String}}
pub struct PreparedInstallation {
    pub workspace: Workspace,
    pub source: SourceBinding,
    pub document_id: String,
    pub event_id: String,
    pub product_revision: String,
    pub batch_revision: String,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("digest installation: {s}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn identifier(s: &str) -> bool {
    !s.is_empty()
        && s.len() <= 128
        && s.bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b':'))
}
fn binding_valid(b: &SourceBinding) -> Result<(), Error> {
    if !identifier(&b.document_id)
        || !b.event_id.strip_prefix("event-").is_some_and(digest)
        || !digest(&b.revision)
    {
        return Err(bad("invalid exact source-event binding"));
    }
    Ok(())
}
pub fn source_role(b: &SourceBinding) -> Result<String, Error> {
    binding_valid(b)?;
    Ok(blake3::hash(
        &postcard::to_allocvec(&("restriction-source-event-1", b)).map_err(|_| Error::TooLarge)?,
    )
    .to_hex()
    .to_string())
}
fn source<'a>(w: &'a Workspace, b: &SourceBinding) -> Result<&'a crate::Molecule, Error> {
    binding_valid(b)?;
    let document = w
        .documents
        .iter()
        .find(|d| d.id == b.document_id)
        .ok_or_else(|| bad("source document is unavailable"))?;
    if document.head != b.revision || document.head_event.as_deref() != Some(b.event_id.as_str()) {
        return Err(bad(
            "source document or event changed; rerun virtual digest",
        ));
    }
    let event = w
        .history
        .iter()
        .find(|e| e.id == b.event_id)
        .ok_or_else(|| bad("source event is unavailable"))?;
    if event.document_id != b.document_id || event.revision != b.revision {
        return Err(bad("source event belongs to another document/revision"));
    }
    w.revisions
        .get(&b.revision)
        .ok_or_else(|| bad("source revision is unavailable"))
}
/// Resolve only the current active exact source event. No producer metadata establishes chemistry.
pub(crate) fn current_source<'a>(
    w: &'a Workspace,
    b: &SourceBinding,
) -> Result<&'a crate::Molecule, Error> {
    w.validate()?;
    if w.active_document != b.document_id {
        return Err(bad("active source changed; rerun virtual digest"));
    }
    source(w, b)
}
pub fn installation_revision(
    w: &Workspace,
    b: &SourceBinding,
    enzymes: &[Enzyme],
    options: &InspectionOptions,
    index: usize,
) -> Result<String, Error> {
    let m = current_source(w, b)?;
    let product =
        restriction_inspection::product_revision(&source_role(b)?, m, enzymes, options, index)?;
    Ok(blake3::hash(
        &postcard::to_allocvec(&(
            "restriction-installation-1",
            run_lineage::workspace_fingerprint(w)?,
            b,
            product,
        ))
        .map_err(|_| Error::TooLarge)?,
    )
    .to_hex()
    .to_string())
}
/// Pure atomic preparation: failure never changes the supplied workspace. The native
/// archive records exact computational ancestry, not measured chemistry or goal verification.
pub fn prepare_install(
    w: &Workspace,
    b: &SourceBinding,
    enzymes: &[Enzyme],
    options: &InspectionOptions,
    index: usize,
) -> Result<PreparedInstallation, Error> {
    let batch_revision = installation_revision(w, b, enzymes, options, index)?;
    let original = source(w, b)?;
    let role = source_role(b)?;
    let product = restriction_inspection::product(&role, original, enzymes, options, index)?;
    product.validate_bindings_for(&role, original, enzymes, options, index)?;
    let enzyme_names = enzymes
        .iter()
        .map(|e| e.name.as_str())
        .collect::<Vec<_>>()
        .join(", ");
    prepare_fresh_install(
        w,
        b,
        product.molecule,
        batch_revision,
        format!(
            "Virtual digest fragment {} with {} · source event {}; complete-cleavage assumptions, computational producer evidence only",
            index + 1,
            enzyme_names,
            b.event_id
        ),
    )
}
/// Install a product reconstructed by a trusted in-crate scientific caller. All producer
/// metadata must be present before deriving the event; no revision is mutated afterward.
pub(crate) fn prepare_fresh_install(
    w: &Workspace,
    b: &SourceBinding,
    product: crate::Molecule,
    batch_revision: String,
    description: String,
) -> Result<PreparedInstallation, Error> {
    current_source(w, b)?;
    product.validate()?;
    let histories = DocumentHistory::restore_all(w)?;
    let history = histories
        .iter()
        .find(|d| d.id == b.document_id)
        .ok_or_else(|| bad("source history is unavailable"))?;
    if history.head_event() != b.event_id || history.molecule.revision() != b.revision {
        return Err(bad("restored source head changed"));
    }
    let prefix = format!("digest-{}", &product.revision()[..24]);
    let mut id = prefix.clone();
    let mut suffix = 1usize;
    while w.documents.iter().any(|d| d.id == id) {
        id = format!("{prefix}-{suffix}");
        suffix = suffix.checked_add(1).ok_or(Error::TooLarge)?;
    }
    let derived = history.derive_event(&b.event_id, id.clone(), product, description)?;
    let event_id = derived.head_event().to_owned();
    let product_revision = derived.molecule.revision();
    let mut next = w.clone();
    derived.contribute(&mut next)?;
    next.active_document = id.clone();
    next.validate()?;
    Ok(PreparedInstallation {
        workspace: next,
        source: b.clone(),
        document_id: id,
        event_id,
        product_revision,
        batch_revision,
    })
}
pub fn install(
    w: &Workspace,
    b: &SourceBinding,
    enzymes: &[Enzyme],
    options: &InspectionOptions,
    index: usize,
) -> Result<InstalledDigest, Error> {
    let p = prepare_install(w, b, enzymes, options, index)?;
    let archive = workspace::save_workspace(&p.workspace)?;
    let result = InstalledDigest {
        archive,
        source: p.source,
        document_id: p.document_id,
        event_id: p.event_id,
        product_revision: p.product_revision,
        batch_revision: p.batch_revision,
    };
    result.validate()?;
    Ok(result)
}
impl InstalledDigest {
    /// Envelope/binding shape only; decode the archive and validate the new UI state
    /// before committing it. This method does not establish scientific authority.
    pub fn validate(&self) -> Result<(), Error> {
        binding_valid(&self.source)?;
        if !identifier(&self.document_id)
            || self.document_id == self.source.document_id
            || !self.event_id.strip_prefix("event-").is_some_and(digest)
            || !digest(&self.product_revision)
            || !digest(&self.batch_revision)
            || self.archive.is_empty()
            || self.archive.len() > 128 * 1024 * 1024
            || !workspace::is_workspace_archive(&self.archive)
        {
            return Err(bad("invalid prepared archive/binding"));
        }
        Ok(())
    }
    pub fn validate_archive(&self, w: &Workspace) -> Result<(), Error> {
        self.validate()?;
        w.validate()?;
        if w.active_document != self.document_id {
            return Err(bad("prepared active document differs"));
        }
        let doc = w
            .documents
            .iter()
            .find(|d| d.id == self.document_id)
            .ok_or_else(|| bad("prepared product document missing"))?;
        if doc.head != self.product_revision
            || doc.head_event.as_deref() != Some(self.event_id.as_str())
        {
            return Err(bad("prepared product revision/event differs"));
        }
        let event = w
            .history
            .iter()
            .find(|e| e.id == self.event_id)
            .ok_or_else(|| bad("prepared product event missing"))?;
        if event.document_id != self.document_id
            || event.revision != self.product_revision
            || event.parents != vec![self.source.event_id.clone()]
        {
            return Err(bad(
                "prepared product must retain only its exact source-event parent",
            ));
        }
        source(w, &self.source)?;
        Ok(())
    }
}
