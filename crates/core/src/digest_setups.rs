//! Named exact-event digest recipes with retained historical inspection evidence.
//! Capture and reopening never scan DNA or restore a fresh scientific assessment.
use crate::{
    Error, Molecule,
    catalog_digest::{CatalogRestrictionInspection, CatalogSelection, InspectionOptions, catalog},
    restriction_history,
    run_lineage::EventBinding,
    workspace::{self, Workspace, WorkspaceAttachment},
};
use serde::{Deserialize, Serialize};
use std::io::{self, Write};

pub const KIND: &str = "named-digest-setup-v1";
pub const SCHEMA_VERSION: u16 = 1;
pub const MAX_SAVED_BYTES: usize = workspace::MAX_ENTRY_BYTES;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum CleavageMode {
    /// The existing ideal-cleavage inspector does not assess methylation,
    /// supplier conditions, partial digestion, or experimental cleavage.
    AssumedCompleteMethylationUnassessed,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Recipe {
    pub source: EventBinding,
    /// Exact catalog edition, content hash, canonical IDs and ordered geometry.
    pub selection: CatalogSelection,
    pub options: InspectionOptions,
    pub cleavage_mode: CleavageMode,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub schema_version: u16,
    pub name: String,
    pub recipe: Recipe,
    /// Historical producer evidence only, including incomplete scans/products.
    pub report: CatalogRestrictionInspection,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Summary {
    pub attachment_hash: String,
    pub name: String,
    pub source: EventBinding,
    pub source_name: String,
    pub catalog_edition: String,
    pub batch_revision: String,
    pub scan_complete: bool,
    pub products_determined: bool,
}

fn bad(message: &str) -> Error {
    Error::Format(format!("digest setup: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn name_valid(name: &str) -> bool {
    !name.is_empty()
        && name == name.trim()
        && name.len() <= 256
        && !name.chars().any(char::is_control)
}
impl Recipe {
    /// Structural historical validation; does not establish installed catalog
    /// membership or grant a measured methylation/cleavage interpretation.
    pub fn validate(&self) -> Result<(), Error> {
        restriction_history::source_role(&self.source)?;
        self.selection.validate()?;
        self.options.validate()
    }
}

// Follow the established native attachment guards. There is no new digest-count
// or document-count cap below the workspace's own capacities.
fn workspace_guard(frame: &Workspace) -> Result<(), Error> {
    if frame.documents.len() > workspace::MAX_DOCUMENTS
        || frame.revisions.len() > workspace::MAX_REVISIONS
        || frame.history.len() > workspace::MAX_HISTORY_EVENTS
        || frame.candidates.len() > workspace::MAX_CANDIDATES
        || frame.attachments.len() > workspace::MAX_ATTACHMENTS
    {
        return Err(Error::TooLarge);
    }
    let mut bytes = 0usize;
    for size in frame
        .revisions
        .values()
        .map(|m| m.sequence.len())
        .chain(frame.attachments.values().map(|asset| asset.bytes.len()))
    {
        if size > workspace::MAX_ENTRY_BYTES {
            return Err(Error::TooLarge);
        }
        bytes = bytes
            .checked_add(size)
            .filter(|n| *n <= workspace::MAX_EXPANDED_BYTES)
            .ok_or(Error::TooLarge)?;
    }
    for molecule in frame.revisions.values() {
        // Unrelated documents retain the workspace's native 100,000-feature
        // capacity; this attachment must not impose its source-analysis budget.
        molecule.validate()?;
    }
    for candidate in &frame.candidates {
        if candidate.plan.inputs.len() > 1000 || candidate.plan.steps.len() > 1000 {
            return Err(Error::TooLarge);
        }
        for input in &candidate.plan.inputs {
            input.molecule.validate()?;
        }
    }
    frame.validate()
}
fn resolve<'a>(frame: &'a Workspace, source: &EventBinding) -> Result<&'a Molecule, Error> {
    restriction_history::source_role(source)?;
    if !frame.documents.iter().any(|d| d.id == source.document_id) {
        return Err(bad("original source document is not retained"));
    }
    let event = frame
        .history
        .iter()
        .find(|event| event.id == source.event_id)
        .ok_or_else(|| bad("exact source event is not retained"))?;
    if event.document_id != source.document_id || event.revision != source.revision {
        return Err(bad("source document/event/revision mismatch"));
    }
    let molecule = frame
        .revisions
        .get(&source.revision)
        .ok_or_else(|| bad("original source revision is unavailable"))?;
    crate::model::bounds::model_preflight(&[molecule])?;
    Ok(molecule)
}
fn historical_binding(frame: &Workspace, record: &Record) -> Result<(), Error> {
    if record.schema_version != SCHEMA_VERSION || !name_valid(&record.name) {
        return Err(bad("unsupported record version or invalid name"));
    }
    record.recipe.validate()?;
    let source = resolve(frame, &record.recipe.source)?;
    record.report.validate()?;
    let inspection = &record.report.inspection;
    if record.report.selection != record.recipe.selection
        || inspection.options != record.recipe.options
        || inspection.source_revision != record.recipe.source.revision
        || inspection.source_role != restriction_history::source_role(&record.recipe.source)?
        || inspection.source_length != source.sequence.len()
        || inspection.topology != source.topology
    {
        return Err(bad("retained report differs from exact source/recipe"));
    }
    // Intentionally do not require the current engine or installed catalog here.
    // Old documentary pins remain inspectable; a new worker run must resolve the
    // exact catalog and apply current source/chemistry rules before any reuse.
    Ok(())
}
struct Limited(Vec<u8>);
impl Write for Limited {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.0.len().saturating_add(bytes.len()) > MAX_SAVED_BYTES {
            return Err(io::Error::other("digest setup byte limit"));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Retain a privately computed catalog inspection under a name. The host checks
/// the live source/draft/epoch before installation; capture validates the current
/// scientific envelope without rerunning the scan. Only attachments are changed.
pub fn capture(
    catalog: &catalog::Catalog,
    frame: &Workspace,
    source: &EventBinding,
    name: &str,
    report: &CatalogRestrictionInspection,
) -> Result<(Workspace, String), Error> {
    if !name_valid(name) {
        return Err(bad(
            "give the digest a trimmed 1–256 byte name without control characters",
        ));
    }
    workspace_guard(frame)?;
    let molecule = resolve(frame, source)?;
    report.validate_bindings_for(
        catalog,
        &report.selection,
        &restriction_history::source_role(source)?,
        molecule,
        &report.inspection.options,
    )?;
    let recipe = Recipe {
        source: source.clone(),
        selection: report.selection.clone(),
        options: report.inspection.options.clone(),
        cleavage_mode: CleavageMode::AssumedCompleteMethylationUnassessed,
    };
    #[derive(Serialize)]
    struct Borrowed<'a> {
        schema_version: u16,
        name: &'a str,
        recipe: &'a Recipe,
        report: &'a CatalogRestrictionInspection,
    }
    let mut encoded = Limited(Vec::new());
    serde_json::to_writer(
        &mut encoded,
        &Borrowed {
            schema_version: SCHEMA_VERSION,
            name,
            recipe: &recipe,
            report,
        },
    )
    .map_err(|_| Error::TooLarge)?;
    let hash = blake3::hash(&encoded.0).to_hex().to_string();
    if let Some(existing) = frame.attachments.get(&hash) {
        if existing.kind != KIND
            || existing.media_type != "application/json"
            || existing.name != name
            || existing.bytes != encoded.0
        {
            return Err(bad("same attachment bytes have inconsistent type/name"));
        }
        return Ok((frame.clone(), hash));
    }
    if frame.attachments.len() == workspace::MAX_ATTACHMENTS {
        return Err(Error::TooLarge);
    }
    let mut staged = frame.clone();
    staged.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: KIND.into(),
            name: name.into(),
            media_type: "application/json".into(),
            bytes: encoded.0,
        },
    );
    staged.validate()?;
    Ok((staged, hash))
}
fn load_guarded(frame: &Workspace, hash: &str) -> Result<Record, Error> {
    if !digest(hash) {
        return Err(bad("invalid setup attachment identity"));
    }
    let asset = frame
        .attachments
        .get(hash)
        .ok_or_else(|| bad("saved digest setup is unavailable"))?;
    if asset.kind != KIND
        || asset.media_type != "application/json"
        || asset.bytes.len() > MAX_SAVED_BYTES
        || blake3::hash(&asset.bytes).to_hex().as_str() != hash
    {
        return Err(bad("invalid setup attachment kind, bytes or size"));
    }
    let record: Record = serde_json::from_slice(&asset.bytes).map_err(|e| bad(&e.to_string()))?;
    historical_binding(frame, &record)?;
    if asset.name != record.name {
        return Err(bad("attachment name differs from recorded setup name"));
    }
    Ok(record)
}
/// Historical evidence only. Neither opening nor listing establishes freshness,
/// catalog availability, reaction conditions, or current digest validity.
pub fn load(frame: &Workspace, hash: &str) -> Result<Record, Error> {
    workspace_guard(frame)?;
    load_guarded(frame, hash)
}
pub fn list(frame: &Workspace) -> Result<Vec<Summary>, Error> {
    workspace_guard(frame)?;
    let mut summaries = Vec::new();
    for (hash, asset) in &frame.attachments {
        if asset.kind != KIND {
            continue;
        }
        let record = load_guarded(frame, hash)?;
        summaries.push(Summary {
            attachment_hash: hash.clone(),
            name: record.name,
            source_name: resolve(frame, &record.recipe.source)?.name.clone(),
            source: record.recipe.source,
            catalog_edition: record.recipe.selection.catalog_edition,
            batch_revision: record.report.batch_revision,
            scan_complete: record.report.inspection.scan_complete,
            products_determined: record.report.inspection.products_determined,
        });
    }
    Ok(summaries)
}
/// Delete only the exact validated setup attachment. No source, document head,
/// operation event or unrelated attachment is removed.
pub fn delete(frame: &Workspace, hash: &str) -> Result<Workspace, Error> {
    workspace_guard(frame)?;
    load_guarded(frame, hash)?;
    let mut staged = frame.clone();
    staged.attachments.remove(hash);
    staged.validate()?;
    Ok(staged)
}
/// Full recipe for a future fresh predicted-gel bridge. This is deliberately not
/// gel::Preparation::Digest: that names-only path loses catalog geometry/options.
/// A consumer must resolve this exact event/catalog and recompute before making
/// current gel predictions; retained fragment evidence alone cannot do that.
pub fn gel_recipe(frame: &Workspace, hash: &str) -> Result<Recipe, Error> {
    Ok(load(frame, hash)?.recipe)
}
