//! Named editable declarations stored separately from the one active draft.
//!
//! Source snapshots are replaced ONLY at the documented Draft typed paths:
//! slots[*].{captured,anneal_bottom} and
//! alternatives.bins[*].choices[*].slot.{captured,anneal_bottom}.
//! Each non-null snapshot becomes {source: SourceInfo}; molecules are recovered
//! only from that exact retained document/event/revision. Closed run-origin
//! capsules remain intact because their embedded evidence cannot be replaced by
//! equal-DNA references. No verification or candidate authority is restored.
use crate::{
    assembly_draft_model::{self, Draft},
    source_identity::SourceInfo,
};
use collagene_core::workspace::{Workspace, WorkspaceAttachment};
use collagene_model::Molecule;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};

pub(crate) const KIND: &str = "assembly-design-v1";
pub(crate) const MAX_DESIGNS: usize = 64;
const MAX_PAYLOAD: usize = 8 << 20;
const MAX_DRAFT: usize = 32 << 20;
const MAX_NAME: usize = 160;
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Entry {
    pub(crate) id: String,
    pub(crate) name: String,
    pub(crate) attachment_hash: String,
    pub(crate) draft_version: u32,
    /// Distinct exact original source events, not repeated library occurrences.
    pub(crate) source_count: usize,
}
/// Opaque saved-design attachment that cannot be interpreted as this schema.
/// Recovery identifies bytes by their validated workspace hash, not a guessed
/// design ID, and does not require any source documents to remain available.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MalformedEntry {
    pub(crate) attachment_hash: String,
    pub(crate) name: String,
    pub(crate) error: String,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Saved {
    schema: u16,
    id: String,
    name: String,
    draft: Value,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SourceRef {
    source: SourceInfo,
}
fn fail(text: &str) -> String {
    format!("Saved assembly design: {text}")
}
fn name(text: &str) -> Result<String, String> {
    let text = text.trim();
    if text.is_empty() || text.len() > MAX_NAME || text.chars().any(char::is_control) {
        return Err(fail(
            "provide a nonempty name of at most160 bytes without control characters.",
        ));
    }
    Ok(text.into())
}
fn valid_id(id: &str) -> bool {
    id.strip_prefix("design-").is_some_and(|s| {
        s.len() == 64
            && s.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    })
}
struct Bounded {
    bytes: Vec<u8>,
    max: usize,
}
impl Write for Bounded {
    fn write(&mut self, value: &[u8]) -> std::io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(value.len())
            .is_none_or(|n| n > self.max)
        {
            return Err(std::io::Error::other("saved design size limit"));
        }
        self.bytes.extend_from_slice(value);
        Ok(value.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encode(value: &impl Serialize, max: usize) -> Result<Vec<u8>, String> {
    let mut writer = Bounded {
        bytes: Vec::new(),
        max,
    };
    serde_json::to_writer(&mut writer, value).map_err(|e| e.to_string())?;
    Ok(writer.bytes)
}
/// Visit only declared snapshot positions, never arbitrary keys inside a closed
/// scientific capsule or user text. Counts are bounded before iterating children.
fn snapshots(
    value: &mut Value,
    mut apply: impl FnMut(&mut Value) -> Result<(), String>,
) -> Result<(), String> {
    fn slot(
        value: &mut Value,
        apply: &mut impl FnMut(&mut Value) -> Result<(), String>,
    ) -> Result<(), String> {
        let object = value
            .as_object_mut()
            .ok_or_else(|| fail("invalid source slot."))?;
        let captured = object
            .get_mut("captured")
            .ok_or_else(|| fail("source slot has no captured field."))?;
        if !captured.is_null() {
            apply(captured)?;
        }
        if let Some(bottom) = object.get_mut("anneal_bottom") {
            if !bottom.is_null() {
                apply(bottom)?;
            }
        }
        Ok(())
    }
    let root = value
        .as_object_mut()
        .ok_or_else(|| fail("invalid editable declaration."))?;
    let slots = root
        .get_mut("slots")
        .and_then(Value::as_array_mut)
        .ok_or_else(|| fail("missing source slots."))?;
    if slots.len() > 12 {
        return Err(fail("too many source slots."));
    }
    for value in slots {
        slot(value, &mut apply)?;
    }
    if let Some(alternatives) = root.get_mut("alternatives").filter(|v| !v.is_null()) {
        let bins = alternatives
            .get_mut("bins")
            .and_then(Value::as_array_mut)
            .ok_or_else(|| fail("missing library bins."))?;
        if bins.len() > 12 {
            return Err(fail("too many source-library bins."));
        }
        let mut count = 0usize;
        for bin in bins {
            let choices = bin
                .get_mut("choices")
                .and_then(Value::as_array_mut)
                .ok_or_else(|| fail("missing source choices."))?;
            count = count
                .checked_add(choices.len())
                .filter(|n| *n <= 256)
                .ok_or_else(|| fail("too many source choices."))?;
            if choices.len() > 64 {
                return Err(fail("too many alternatives in one bin."));
            }
            for choice in choices {
                slot(
                    choice
                        .get_mut("slot")
                        .ok_or_else(|| fail("missing library source slot."))?,
                    &mut apply,
                )?;
            }
        }
    }
    Ok(())
}
/// One borrowed index per save/load transaction. Workspace validation has
/// already established graph membership and revision hashes. Referenced sources
/// use logarithmic exact event lookup even in long histories/source libraries.
struct Resolver<'a> {
    documents: BTreeSet<&'a str>,
    events: BTreeMap<(&'a str, &'a str), &'a str>,
    revisions: &'a BTreeMap<String, Molecule>,
}
impl<'a> Resolver<'a> {
    fn new(frame: &'a Workspace) -> Self {
        Self {
            documents: frame.documents.iter().map(|d| d.id.as_str()).collect(),
            events: frame
                .history
                .iter()
                .map(|event| {
                    (
                        (event.document_id.as_str(), event.id.as_str()),
                        event.revision.as_str(),
                    )
                })
                .collect(),
            revisions: &frame.revisions,
        }
    }
    fn resolve(&self, source: &SourceInfo) -> Result<&'a Molecule, String> {
        if !self.documents.contains(source.binding.document_id.as_str()) {
            return Err(fail("the exact original source document is unavailable."));
        }
        let revision=self.events.get(&(
            source.binding.document_id.as_str(),source.binding.event_id.as_str(),
        )).copied().ok_or_else(||fail("the exact retained source event is unavailable; equal DNA or a current head cannot substitute."))?;
        if revision != source.binding.revision {
            return Err(fail("the original source event revision changed."));
        }
        let molecule = self
            .revisions
            .get(revision)
            .ok_or_else(|| fail("the retained source revision is missing."))?;
        if molecule.name != source.name
            || molecule.sequence.len() != source.length
            || molecule.sequence.alphabet() != source.alphabet
            || molecule.topology != source.topology
        {
            return Err(fail(
                "source description differs from its exact retained event.",
            ));
        }
        Ok(molecule)
    }
}
fn source_reference(value: &Value) -> Result<SourceRef, String> {
    let parsed: SourceRef = serde_json::from_value(value.clone()).map_err(|e| e.to_string())?;
    // SourceInfo's general-purpose decoder permits extension fields. This
    // persisted schema deliberately accepts exactly the declared reference shape.
    if serde_json::to_value(&parsed).map_err(|e| e.to_string())? != *value {
        return Err(fail("source reference has unknown or noncanonical fields."));
    }
    Ok(parsed)
}
fn version(value: &Value) -> Result<u32, String> {
    value
        .get("version")
        .and_then(Value::as_u64)
        .filter(|v| *v == u64::from(assembly_draft_model::CURRENT_VERSION))
        .map(|v| v as u32)
        .ok_or_else(|| fail("unsupported editable-draft version."))
}
fn source_values(value: &Value) -> Result<Vec<&Value>, String> {
    fn slot<'a>(value: &'a Value, output: &mut Vec<&'a Value>) -> Result<(), String> {
        let object = value
            .as_object()
            .ok_or_else(|| fail("invalid source slot."))?;
        let captured = object
            .get("captured")
            .ok_or_else(|| fail("source slot has no captured field."))?;
        if !captured.is_null() {
            output.push(captured);
        }
        if let Some(bottom) = object.get("anneal_bottom").filter(|v| !v.is_null()) {
            output.push(bottom);
        }
        Ok(())
    }
    let root = value
        .as_object()
        .ok_or_else(|| fail("invalid editable declaration."))?;
    let slots = root
        .get("slots")
        .and_then(Value::as_array)
        .ok_or_else(|| fail("missing source slots."))?;
    if slots.len() > 12 {
        return Err(fail("too many source slots."));
    }
    let mut output = Vec::new();
    for value in slots {
        slot(value, &mut output)?;
    }
    if let Some(alternatives) = root.get("alternatives").filter(|v| !v.is_null()) {
        let bins = alternatives
            .get("bins")
            .and_then(Value::as_array)
            .ok_or_else(|| fail("missing library bins."))?;
        if bins.len() > 12 {
            return Err(fail("too many source-library bins."));
        }
        let mut count = 0usize;
        for bin in bins {
            let choices = bin
                .get("choices")
                .and_then(Value::as_array)
                .ok_or_else(|| fail("missing source choices."))?;
            count = count
                .checked_add(choices.len())
                .filter(|n| *n <= 256)
                .ok_or_else(|| fail("too many source choices."))?;
            if choices.len() > 64 {
                return Err(fail("too many alternatives in one bin."));
            }
            for choice in choices {
                slot(
                    choice
                        .get("slot")
                        .ok_or_else(|| fail("missing library source slot."))?,
                    &mut output,
                )?;
            }
        }
    }
    Ok(output)
}
fn count_sources(draft: &Value) -> Result<usize, String> {
    let mut refs = BTreeSet::new();
    for value in source_values(draft)? {
        let source = source_reference(value)?.source;
        refs.insert((
            source.binding.document_id,
            source.binding.event_id,
            source.binding.revision,
        ));
    }
    Ok(refs.len())
}
fn stored(attachment: &WorkspaceAttachment) -> Result<Saved, String> {
    if attachment.kind != KIND
        || attachment.media_type != "application/json"
        || attachment.bytes.len() > MAX_PAYLOAD
    {
        return Err(fail("invalid attachment kind, type or size."));
    }
    let saved: Saved = serde_json::from_slice(&attachment.bytes).map_err(|e| e.to_string())?;
    if saved.schema != 1
        || !valid_id(&saved.id)
        || name(&saved.name)? != saved.name
        || saved.name != attachment.name
    {
        return Err(fail("saved identity/name/schema is invalid."));
    }
    version(&saved.draft)?;
    count_sources(&saved.draft)?;
    Ok(saved)
}
fn entries(frame: &Workspace) -> Result<Vec<(Entry, Saved)>, String> {
    let assets = frame
        .attachments
        .iter()
        .filter(|(_, a)| a.kind == KIND)
        .collect::<Vec<_>>();
    if assets.len() > MAX_DESIGNS {
        return Err(fail("at most64 named designs are supported."));
    }
    let mut ids = BTreeSet::new();
    let mut names = BTreeSet::new();
    let mut result = Vec::with_capacity(assets.len());
    for (hash, asset) in assets {
        let saved = stored(asset)?;
        if !ids.insert(saved.id.clone()) || !names.insert(saved.name.to_lowercase()) {
            return Err(fail("duplicate saved identity or name."));
        }
        let entry = Entry {
            id: saved.id.clone(),
            name: saved.name.clone(),
            attachment_hash: hash.clone(),
            draft_version: version(&saved.draft)?,
            source_count: count_sources(&saved.draft)?,
        };
        result.push((entry, saved));
    }
    result.sort_by(|(a, _), (b, _)| {
        a.name
            .to_lowercase()
            .cmp(&b.name.to_lowercase())
            .then_with(|| a.id.cmp(&b.id))
    });
    Ok(result)
}
pub(crate) fn list(frame: &Workspace) -> Result<Vec<Entry>, String> {
    frame.validate().map_err(|e| e.to_string())?;
    Ok(entries(frame)?
        .into_iter()
        .map(|(entry, _)| entry)
        .collect())
}
/// Keep normal listing strict while allowing an explicitly selected unreadable
/// named-design attachment to be removed without blocking the rest of the list.
pub(crate) fn malformed(frame: &Workspace) -> Result<Vec<MalformedEntry>, String> {
    frame.validate().map_err(|e| e.to_string())?;
    Ok(frame
        .attachments
        .iter()
        .filter(|(_, asset)| asset.kind == KIND)
        .filter_map(|(hash, asset)| {
            stored(asset).err().map(|error| MalformedEntry {
                attachment_hash: hash.clone(),
                name: asset.name.clone(),
                error: error.chars().take(256).collect(),
            })
        })
        .collect())
}
/// Remove only the exact current invalid named-design bytes. Valid designs use
/// the normal ID-bound delete path; unrelated attachments can never be selected.
pub(crate) fn discard_malformed(frame: &Workspace, hash: &str) -> Result<Workspace, String> {
    frame.validate().map_err(|e| e.to_string())?;
    let asset = frame
        .attachments
        .get(hash)
        .ok_or_else(|| fail("the unreadable saved-design attachment is unavailable."))?;
    if asset.kind != KIND || stored(asset).is_ok() {
        return Err(fail(
            "recovery only removes the exact unreadable named-design attachment.",
        ));
    }
    let mut next = frame.clone();
    next.attachments.remove(hash);
    next.validate().map_err(|e| e.to_string())?;
    Ok(next)
}
/// Save a new named snapshot. This never replaces the active autosave or an
/// existing named design. Names refer to saved snapshots; Draft.name stays intact.
pub(crate) fn save(
    frame: &Workspace,
    display_name: &str,
    draft: &Draft,
) -> Result<(Workspace, Entry), String> {
    frame.validate().map_err(|e| e.to_string())?;
    draft.validate()?;
    let name = name(display_name)?;
    let existing = entries(frame)?;
    if existing.len() >= MAX_DESIGNS
        || existing
            .iter()
            .any(|(entry, _)| entry.name.to_lowercase() == name.to_lowercase())
    {
        return Err(fail(
            "the name already exists or the64-design limit was reached.",
        ));
    }
    let raw = assembly_draft_model::encode_raw(draft)?;
    let mut value: Value = serde_json::from_slice(&raw).map_err(|e| e.to_string())?;
    let resolver = Resolver::new(frame);
    snapshots(&mut value, |value| {
        let captured: assembly_draft_model::CapturedSource =
            serde_json::from_value(value.clone()).map_err(|e| e.to_string())?;
        let actual = resolver.resolve(&captured.source)?;
        if actual != &captured.molecule {
            return Err(fail(
                "pinned source snapshot differs from its retained original event.",
            ));
        }
        *value = serde_json::to_value(SourceRef {
            source: captured.source,
        })
        .map_err(|e| e.to_string())?;
        Ok(())
    })?;
    let identity = encode(
        &(&name, &value, frame.attachments.keys().collect::<Vec<_>>()),
        MAX_DRAFT,
    )?;
    let id = (0..=MAX_DESIGNS)
        .map(|nonce| {
            let mut hash = blake3::Hasher::new();
            hash.update(b"collagene-saved-design-1");
            hash.update(&identity);
            hash.update(&(nonce as u64).to_le_bytes());
            format!("design-{}", hash.finalize().to_hex())
        })
        .find(|id| !existing.iter().any(|(entry, _)| &entry.id == id))
        .ok_or_else(|| fail("no unused design identity is available."))?;
    let saved = Saved {
        schema: 1,
        id,
        name,
        draft: value,
    };
    let bytes = encode(&saved, MAX_PAYLOAD)?;
    let hash = blake3::hash(&bytes).to_hex().to_string();
    if frame.attachments.contains_key(&hash) {
        return Err(fail("an attachment identity collision prevents saving."));
    }
    let entry = Entry {
        id: saved.id.clone(),
        name: saved.name.clone(),
        attachment_hash: hash.clone(),
        draft_version: version(&saved.draft)?,
        source_count: count_sources(&saved.draft)?,
    };
    let mut next = frame.clone();
    next.attachments.insert(
        hash,
        WorkspaceAttachment {
            kind: KIND.into(),
            name: saved.name,
            media_type: "application/json".into(),
            bytes,
        },
    );
    next.validate().map_err(|e| e.to_string())?;
    Ok((next, entry))
}
struct Count {
    bytes: usize,
    max: usize,
}
impl Write for Count {
    fn write(&mut self, value: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(value.len())
            .filter(|n| *n <= self.max)
            .ok_or_else(|| std::io::Error::other("saved design hydrated-size limit"))?;
        Ok(value.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn json_size(value: &impl Serialize, max: usize) -> Result<usize, String> {
    let mut counter = Count { bytes: 0, max };
    serde_json::to_writer(&mut counter, value).map_err(|e| e.to_string())?;
    Ok(counter.bytes)
}
fn hydration_preflight(resolver: &Resolver<'_>, draft: &Value) -> Result<(), String> {
    #[derive(Serialize)]
    struct BorrowedCapture<'a> {
        source: &'a SourceInfo,
        molecule: &'a Molecule,
    }
    let mut predicted = json_size(draft, MAX_DRAFT)?;
    for value in source_values(draft)? {
        let reference = source_reference(value)?;
        let molecule = resolver.resolve(&reference.source)?;
        collagene_model::bounds::model_preflight(&[molecule]).map_err(|e| e.to_string())?;
        let original = json_size(value, MAX_PAYLOAD)?;
        let full = json_size(
            &BorrowedCapture {
                source: &reference.source,
                molecule,
            },
            MAX_DRAFT,
        )?;
        predicted = predicted
            .checked_sub(original)
            .and_then(|n| n.checked_add(full))
            .filter(|n| *n <= MAX_DRAFT)
            .ok_or_else(|| fail("hydrated declaration exceeds32MiB before source copies."))?;
    }
    Ok(())
}

/// Hydrate editable state only. Source events must remain available exactly;
/// selecting this snapshot never installs a candidate or restores a positive check.
pub(crate) fn load(frame: &Workspace, id: &str) -> Result<Draft, String> {
    frame.validate().map_err(|e| e.to_string())?;
    let (_, mut saved) = entries(frame)?
        .into_iter()
        .find(|(entry, _)| entry.id == id)
        .ok_or_else(|| fail("the named design is unavailable."))?;
    let resolver = Resolver::new(frame);
    hydration_preflight(&resolver, &saved.draft)?;
    snapshots(&mut saved.draft, |value| {
        let reference = source_reference(value)?;
        let molecule = resolver.resolve(&reference.source)?;
        *value = serde_json::to_value(assembly_draft_model::CapturedSource {
            source: reference.source,
            molecule: molecule.clone(),
        })
        .map_err(|e| e.to_string())?;
        Ok(())
    })?;
    let hydrated = encode(&saved.draft, MAX_DRAFT)?;
    let draft = assembly_draft_model::decode(&hydrated)?;
    // Native decoding may reject/remove a stale run-origin declaration. Refuse
    // any silent change rather than dropping a closed capsule from a saved design.
    let canonical = assembly_draft_model::encode_raw(&draft)?;
    let mut expected = saved.draft;
    expected["version"] = Value::from(draft.version);
    if serde_json::from_slice::<Value>(&canonical).map_err(|e| e.to_string())? != expected {
        return Err(fail(
            "native draft recovery would change saved declarations or closed source evidence.",
        ));
    }
    Ok(draft)
}
pub(crate) fn rename(frame: &Workspace, id: &str, new_name: &str) -> Result<Workspace, String> {
    frame.validate().map_err(|e| e.to_string())?;
    let name = name(new_name)?;
    let existing = entries(frame)?;
    if existing
        .iter()
        .any(|(entry, _)| entry.id != id && entry.name.to_lowercase() == name.to_lowercase())
    {
        return Err(fail("the saved-design name already exists."));
    }
    let (entry, mut saved) = existing
        .into_iter()
        .find(|(entry, _)| entry.id == id)
        .ok_or_else(|| fail("the named design is unavailable."))?;
    if saved.name == name {
        return Ok(frame.clone());
    }
    saved.name = name.clone();
    let bytes = encode(&saved, MAX_PAYLOAD)?;
    let hash = blake3::hash(&bytes).to_hex().to_string();
    if frame.attachments.contains_key(&hash) {
        return Err(fail("an attachment identity collision prevents rename."));
    }
    let mut next = frame.clone();
    next.attachments.remove(&entry.attachment_hash);
    next.attachments.insert(
        hash,
        WorkspaceAttachment {
            kind: KIND.into(),
            name,
            media_type: "application/json".into(),
            bytes,
        },
    );
    next.validate().map_err(|e| e.to_string())?;
    Ok(next)
}
pub(crate) fn delete(frame: &Workspace, id: &str) -> Result<Workspace, String> {
    frame.validate().map_err(|e| e.to_string())?;
    let (entry, _) = entries(frame)?
        .into_iter()
        .find(|(entry, _)| entry.id == id)
        .ok_or_else(|| fail("the named design is unavailable."))?;
    let mut next = frame.clone();
    next.attachments.remove(&entry.attachment_hash);
    next.validate().map_err(|e| e.to_string())?;
    Ok(next)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        assembly_draft_model::CapturedSource,
        assembly_edit_model, assembly_substitution_model,
        cloning_methods_model::{Mode, OverlapPreparation},
    };
    use collagene_core::{history::DocumentHistory, run_lineage::EventBinding};
    use collagene_model::Topology;
    use std::sync::Arc;
    fn captured(frame: &Workspace, index: usize) -> Arc<CapturedSource> {
        let doc = &frame.documents[index];
        let molecule = frame.revisions[&doc.head].clone();
        Arc::new(CapturedSource {
            source: SourceInfo {
                binding: EventBinding {
                    document_id: doc.id.clone(),
                    event_id: doc.head_event.clone().unwrap(),
                    revision: doc.head.clone(),
                },
                name: molecule.name.clone(),
                length: molecule.sequence.len(),
                alphabet: molecule.sequence.alphabet(),
                topology: molecule.topology,
            },
            molecule,
        })
    }
    fn fixture() -> (Workspace, Draft) {
        let mut frame = Workspace::from_documents(
            "Exact source events",
            vec![
                Molecule::dna("Original", "AACCGGTT", Topology::Linear).unwrap(),
                Molecule::dna("Bottom", "TTGGCCAA", Topology::Linear).unwrap(),
            ],
            0,
        )
        .unwrap();
        let mut draft = Draft::default();
        draft.name = "Unchanged intended construct name".into();
        for slot in &mut draft.slots {
            slot.captured = Some(captured(&frame, 0));
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    start: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = 0;
            }
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    end: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = 8;
            }
        }
        draft.sync_rows();
        draft.methods.rows[0].forward = "ac gt".into();
        draft.methods.overlaps = " unfinished overlap ".into();
        draft.validate().unwrap();
        let bytes = assembly_draft_model::encode_raw(&draft).unwrap();
        frame.attachments.insert(
            blake3::hash(&bytes).to_hex().to_string(),
            WorkspaceAttachment {
                kind: crate::assembly_draft_model::ATTACHMENT_KIND.into(),
                name: "Active draft".into(),
                media_type: "application/json".into(),
                bytes,
            },
        );
        frame.validate().unwrap();
        (frame, draft)
    }
    fn replace_payload(
        frame: &Workspace,
        hash: &str,
        mutate: impl FnOnce(&mut Saved),
    ) -> Workspace {
        let mut next = frame.clone();
        let mut saved = stored(next.attachments.get(hash).unwrap()).unwrap();
        mutate(&mut saved);
        let bytes = encode(&saved, MAX_PAYLOAD).unwrap();
        let mut asset = next.attachments.remove(hash).unwrap();
        asset.bytes = bytes;
        asset.name = saved.name;
        next.attachments
            .insert(blake3::hash(&asset.bytes).to_hex().to_string(), asset);
        next.validate().unwrap();
        next
    }
    #[test]
    fn current_and_literal_declarations_roundtrip_without_source_copies() {
        let (frame, base) = fixture();
        let before = serde_json::to_value(&frame).unwrap();
        {
            let draft = base.clone();
            draft.validate().unwrap();
            let (saved, entry) = save(&frame, "Current draft", &draft).unwrap();
            assert_eq!(entry.source_count, 1);
            let restored = load(&saved, &entry.id).unwrap();
            assert!(restored == draft);
            let text =
                String::from_utf8(saved.attachments[&entry.attachment_hash].bytes.clone()).unwrap();
            assert!(!text.contains("\"molecule\""));
            assert!(!text.contains("AACCGGTT"));
            assert_eq!(restored.methods.rows[0].forward, "ac gt");
            assert_eq!(restored.methods.overlaps, " unfinished overlap ");
            assert!(saved.candidates.is_empty());
            assert!(saved.active_candidate.is_none());
        }
        let substituted = assembly_substitution_model::declare(&base, 0, "3", "T", None).unwrap();
        let edited = assembly_edit_model::declare(&base, 0, "3", "5", "TAG", None).unwrap();
        for draft in [substituted, edited] {
            let (saved, entry) = save(&frame, "Literal design", &draft).unwrap();
            assert!(load(&saved, &entry.id).unwrap() == draft);
        }
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
    }
    #[test]
    fn library_choice_and_paired_bottom_sources_use_only_their_typed_reference_paths() {
        let (frame, base) = fixture();
        let library = crate::assembly_alternatives::enable(&base, "saved").unwrap();
        let (saved, entry) = save(&frame, "Source alternatives", &library).unwrap();
        let payload = stored(&saved.attachments[&entry.attachment_hash]).unwrap();
        let choices = payload.draft["alternatives"]["bins"].as_array().unwrap();
        assert!(choices.iter().all(|bin| {
            bin["choices"][0]["slot"]["captured"]
                .get("molecule")
                .is_none()
        }));
        assert_eq!(source_values(&payload.draft).unwrap().len(), 4);
        assert!(load(&saved, &entry.id).unwrap() == library);

        let mut pair = base;
        pair.methods.mode = Mode::Ligation;
        pair.methods.rows[1].preparation = OverlapPreparation::AnnealedOligos;
        pair.methods.rows[1].anneal_options.min_paired_bases = 4;
        pair.slots[1].anneal_bottom = Some(captured(&frame, 1));
        pair.sync_rows();
        pair.validate().unwrap();
        let (saved, entry) = save(&frame, "Paired strands", &pair).unwrap();
        assert_eq!(entry.source_count, 2);
        let payload = stored(&saved.attachments[&entry.attachment_hash]).unwrap();
        assert!(
            payload.draft["slots"][1]["anneal_bottom"]
                .get("molecule")
                .is_none()
        );
        assert!(load(&saved, &entry.id).unwrap() == pair);
    }
    #[test]
    fn closed_run_origin_capsule_survives_named_design_without_replacing_its_evidence() {
        use collagene_core::{
            analysis::Primer,
            cloning_outputs,
            run_lineage::InputEventBinding,
            run_material,
            verification::{DesignPlan, Operation, PlanInput, Source},
        };
        let source =
            Molecule::dna("Original PCR parent", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
        let target = Molecule::dna(
            "Independent PCR target",
            "ACGACCCCGGGGTAAC",
            Topology::Linear,
        )
        .unwrap();
        let frame =
            Workspace::from_documents("Retained run", vec![source.clone(), target.clone()], 0)
                .unwrap();
        let source_event = captured(&frame, 0).source.binding.clone();
        let target_event = captured(&frame, 1).source.binding.clone();
        let plan = DesignPlan {
            name: "Actual original-source PCR".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "original".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps: vec![Operation::Pcr {
                template_strand: collagene_core::chemistry::Strand::Top,
                template: Source::Input("original".into()),
                forward: Primer::new("Original forward", "ACGACCCC", "").unwrap(),
                reverse: Primer::new("Original reverse", "GTTACCCC", "").unwrap(),
            }],
        };
        let committed = cloning_outputs::finalize(
            &frame,
            &plan,
            target_event,
            vec![InputEventBinding {
                input_id: "original".into(),
                event: source_event,
            }],
            &cloning_outputs::Options {
                name: "Saved actual amplicon".into(),
                save_primers: false,
                show_prepared: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let selection = run_material::Selection {
            run_id: committed.receipt.run.run_id.clone(),
            step_index: 0,
            event: committed.receipt.construct.clone(),
        };
        let capsule = run_material::capsule(&committed.workspace, &selection).unwrap();
        let exact_capsule = serde_json::to_value(&capsule).unwrap();
        let index = committed
            .workspace
            .documents
            .iter()
            .position(|d| d.id == selection.event.document_id)
            .unwrap();
        let capture = captured(&committed.workspace, index);
        let mut draft = Draft::default();
        draft.circular = false;
        draft.methods.mode = Mode::Ligation;
        draft.slots.truncate(1);
        draft.methods.rows.truncate(1);
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = capture.molecule.sequence.len();
        }
        draft.slots[0].captured = Some(capture);
        draft.methods.rows[0].preparation = OverlapPreparation::WholeLinear;
        draft.methods.rows[0].run_origin = Some(
            crate::run_origin_model::PinnedOrigin::from_material(
                collagene_core::run_material::MaterialOrigin::Run(capsule),
            )
            .unwrap(),
        );
        draft.slots[0].material_origin = draft.methods.rows[0].run_origin.clone();
        draft.sync_rows();
        draft.validate().unwrap();
        assert_eq!(draft.version, assembly_draft_model::CURRENT_VERSION);
        let before = serde_json::to_value(&committed.workspace).unwrap();
        let (saved, entry) = save(&committed.workspace, "Prepared source design", &draft).unwrap();
        let restored = load(&saved, &entry.id).unwrap();
        assert!(restored == draft);
        let collagene_core::run_material::MaterialOrigin::Run(restored_origin) =
            restored.methods.rows[0]
                .run_origin
                .as_ref()
                .unwrap()
                .as_ref()
        else {
            panic!("Saved run material changed origin type");
        };
        assert_eq!(
            serde_json::to_value(restored_origin).unwrap(),
            exact_capsule
        );
        let (_, rows) = crate::cloning_preparations::prepared_rows(&restored.methods, 1).unwrap();
        let collagene_core::assembly_recipe::Preparation::FromRun { origin, .. } =
            &rows[0].preparation
        else {
            panic!("Saved material must retain explicit native reconstruction");
        };
        assert_eq!(serde_json::to_value(origin).unwrap(), exact_capsule);
        assert_eq!(serde_json::to_value(&committed.workspace).unwrap(), before);
        assert_eq!(saved.candidates.len(), committed.workspace.candidates.len());
        assert_eq!(saved.active_candidate, committed.workspace.active_candidate);
    }
    #[test]
    fn historical_events_are_restored_exactly_and_equal_dna_in_another_document_is_not_a_substitute()
     {
        let (frame, draft) = fixture();
        let (saved, entry) = save(&frame, "Pinned history", &draft).unwrap();
        let mut docs = DocumentHistory::restore_all(&saved).unwrap();
        let replacement =
            Molecule::dna("New current reference", "TTTTGGGG", Topology::Linear).unwrap();
        docs[0]
            .commit(replacement.clone(), "Reference edit")
            .unwrap();
        let mut advanced = saved.clone();
        advanced.documents.clear();
        advanced.revisions.clear();
        advanced.history.clear();
        for doc in &docs {
            doc.contribute(&mut advanced).unwrap();
        }
        advanced.validate().unwrap();
        let restored = load(&advanced, &entry.id).unwrap();
        assert!(restored == draft);
        assert_ne!(
            advanced.documents[0].head,
            restored.slots[0]
                .captured
                .as_ref()
                .unwrap()
                .source
                .binding
                .revision
        );

        let mut missing = saved.clone();
        missing.documents.clear();
        missing.revisions.clear();
        missing.history.clear();
        DocumentHistory::new(saved.documents[0].id.clone(), replacement)
            .unwrap()
            .contribute(&mut missing)
            .unwrap();
        DocumentHistory::new(
            saved.documents[1].id.clone(),
            draft.slots[0].captured.as_ref().unwrap().molecule.clone(),
        )
        .unwrap()
        .contribute(&mut missing)
        .unwrap();
        missing.validate().unwrap();
        let before = serde_json::to_value(&missing).unwrap();
        assert!(
            load(&missing, &entry.id)
                .err()
                .unwrap()
                .contains("source event is unavailable")
        );
        assert!(save(&missing, "No guessed source", &draft).is_err());
        assert_eq!(serde_json::to_value(&missing).unwrap(), before);
    }
    #[test]
    fn forged_references_descriptions_and_positive_fields_are_not_hydrated() {
        let (frame, draft) = fixture();
        let (saved, entry) = save(&frame, "Bound declaration", &draft).unwrap();
        for changed in [
            replace_payload(&saved, &entry.attachment_hash, |s| {
                s.draft["slots"][0]["captured"]["source"]["binding"]["event_id"] =
                    Value::from("0".repeat(64))
            }),
            replace_payload(&saved, &entry.attachment_hash, |s| {
                s.draft["slots"][0]["captured"]["source"]["name"] = Value::from("Forged label")
            }),
            replace_payload(&saved, &entry.attachment_hash, |s| {
                s.draft["verification"] = serde_json::json!({"status":"ConsistentUnderAssumptions"})
            }),
        ] {
            let before = serde_json::to_value(&changed).unwrap();
            assert!(load(&changed, &entry.id).is_err());
            assert_eq!(serde_json::to_value(&changed).unwrap(), before);
        }
        let mut changed = draft;
        let mut source = changed.slots[0].captured.as_ref().unwrap().as_ref().clone();
        source.molecule.sequence = collagene_model::Sequence::dna("AAAAAAAA").unwrap();
        changed.slots[0].captured = Some(Arc::new(source));
        assert!(save(&frame, "Changed snapshot", &changed).is_err());
    }
    #[test]
    fn renaming_and_deleting_are_atomic_and_preserve_active_autosave_and_originals() {
        let (frame, draft) = fixture();
        let active = frame.attachments.clone();
        let before = serde_json::to_value(&frame).unwrap();
        let (saved, entry) = save(&frame, "First saved design", &draft).unwrap();
        let renamed = rename(&saved, &entry.id, "Renamed design").unwrap();
        assert_eq!(list(&renamed).unwrap()[0].id, entry.id);
        assert_eq!(list(&renamed).unwrap()[0].name, "Renamed design");
        assert!(load(&renamed, &entry.id).unwrap() == draft);
        assert_ne!(
            list(&renamed).unwrap()[0].attachment_hash,
            entry.attachment_hash
        );
        let (second, second_entry) = save(&renamed, "Second saved design", &draft).unwrap();
        let unchanged = serde_json::to_value(&second).unwrap();
        assert!(rename(&second, &entry.id, "second saved design").is_err());
        assert!(save(&second, " RENAMED DESIGN ", &draft).is_err());
        assert!(rename(&second, &entry.id, "x".repeat(161).as_str()).is_err());
        assert_eq!(serde_json::to_value(&second).unwrap(), unchanged);
        let removed = delete(&second, &entry.id).unwrap();
        assert_eq!(list(&removed).unwrap()[0].id, second_entry.id);
        assert!(load(&removed, &entry.id).is_err());
        let mut different = draft.clone();
        different.name = "Another intended construct".into();
        let (recreated, new_entry) = save(&removed, "Renamed design", &different).unwrap();
        assert_ne!(new_entry.id, entry.id);
        for (hash, asset) in active {
            assert_eq!(
                serde_json::to_value(&recreated.attachments[&hash]).unwrap(),
                serde_json::to_value(asset).unwrap()
            );
        }
        assert_eq!(serde_json::to_value(&frame).unwrap(), before);
    }
    #[test]
    fn unreadable_named_attachment_can_be_recovered_without_sources_or_touching_other_data() {
        let (frame, draft) = fixture();
        let (saved, healthy) = save(&frame, "Healthy named design", &draft).unwrap();
        let mut damaged = saved.clone();
        let bytes = br#"{"schema":1,"not_a_saved_design":true}"#.to_vec();
        let bad_hash = blake3::hash(&bytes).to_hex().to_string();
        damaged.attachments.insert(
            bad_hash.clone(),
            WorkspaceAttachment {
                kind: KIND.into(),
                name: "Unreadable named design".into(),
                media_type: "application/json".into(),
                bytes,
            },
        );
        // Orphaned source references must not prevent attachment-only cleanup.
        damaged.documents.clear();
        damaged.revisions.clear();
        damaged.history.clear();
        DocumentHistory::new(
            "replacement-source",
            Molecule::dna("Other current document", "TTTTAAAA", Topology::Linear).unwrap(),
        )
        .unwrap()
        .contribute(&mut damaged)
        .unwrap();
        damaged.active_document = "replacement-source".into();
        damaged.validate().unwrap();
        let before = serde_json::to_value(&damaged).unwrap();
        assert!(list(&damaged).is_err());
        let issues = malformed(&damaged).unwrap();
        assert_eq!(issues.len(), 1);
        assert_eq!(issues[0].attachment_hash, bad_hash);
        assert_eq!(issues[0].name, "Unreadable named design");
        assert!(!issues[0].error.is_empty());
        assert!(discard_malformed(&damaged, &healthy.attachment_hash).is_err());
        let active_hash = damaged
            .attachments
            .iter()
            .find(|(_, asset)| asset.kind == crate::assembly_draft_model::ATTACHMENT_KIND)
            .unwrap()
            .0;
        assert!(discard_malformed(&damaged, active_hash).is_err());
        let cleaned = discard_malformed(&damaged, &bad_hash).unwrap();
        assert_eq!(list(&cleaned).unwrap(), vec![healthy.clone()]);
        assert!(load(&cleaned, &healthy.id).is_err());
        assert!(malformed(&cleaned).unwrap().is_empty());
        let mut expected = damaged.clone();
        expected.attachments.remove(&bad_hash);
        assert_eq!(
            serde_json::to_value(&cleaned).unwrap(),
            serde_json::to_value(&expected).unwrap()
        );
        assert_eq!(serde_json::to_value(&damaged).unwrap(), before);
    }
    #[test]
    fn repeated_source_references_are_bounded_before_rich_hydration_copies() {
        let mut molecule = Molecule::dna("Rich original", "AACCGGTT", Topology::Linear).unwrap();
        molecule.features.push(collagene_model::Feature {
            kind: "misc_feature".into(),
            location: collagene_model::Location::range(0, 8),
            qualifiers: (0..3)
                .map(|i| (format!("note{i}"), Some("x".repeat(1_000_000))))
                .collect(),
            needs_review: false,
        });
        let frame = Workspace::from_documents("Rich source", vec![molecule], 0).unwrap();
        let mut draft = Draft::default();
        for slot in &mut draft.slots {
            slot.captured = Some(captured(&frame, 0));
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    start: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = 0;
            }
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    end: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = 8;
            }
        }
        draft.sync_rows();
        draft.validate().unwrap();
        let (saved, entry) = save(&frame, "Source reference budget", &draft).unwrap();
        assert!(saved.attachments[&entry.attachment_hash].bytes.len() < 20_000);
        let excessive = replace_payload(&saved, &entry.attachment_hash, |saved| {
            let slot = saved.draft["slots"][0].clone();
            let row = saved.draft["methods"]["rows"][0].clone();
            saved.draft["slots"] = Value::Array(vec![slot; 12]);
            saved.draft["methods"]["rows"] = Value::Array(
                (0..12)
                    .map(|i| {
                        let mut row = row.clone();
                        row["document"] = Value::from(i);
                        row
                    })
                    .collect(),
            );
        });
        let error = load(&excessive, &entry.id).err().unwrap();
        assert!(error.contains("before source copies"), "{error}");
    }

    #[test]
    fn saved_policy_keeps_inactive_unfinished_settings_and_original_event_pins() {
        let (frame, mut draft) = fixture();
        let sources = draft
            .slots
            .iter()
            .map(|s| s.captured.as_ref().unwrap().source.clone())
            .collect::<Vec<_>>();
        let mut form = crate::cloning_methods_model::HomologyForm::default();
        form.enabled = false;
        form.exact_homology = false;
        form.cuts = "12:unfinished, :".into();
        form.outer_right = "2e".into();
        draft.methods.homology = Some(form);
        draft.sync_rows();
        let (saved, entry) = save(&frame, "Physical homology settings", &draft).unwrap();
        assert_eq!(entry.draft_version, assembly_draft_model::CURRENT_VERSION);
        let restored = load(&saved, &entry.id).unwrap();
        assert!(restored == draft);
        assert!(
            restored
                .slots
                .iter()
                .map(|s| s.captured.as_ref().unwrap().source.clone())
                .collect::<Vec<_>>()
                == sources
        );
        assert_eq!(
            serde_json::to_value(&saved.documents).unwrap(),
            serde_json::to_value(&frame.documents).unwrap()
        );
        assert!(saved.candidates.is_empty());
    }
    #[test]
    fn saved_digest_chain_retains_unfinished_six_oligos_pins_and_source_events() {
        let (frame, mut draft) = fixture();
        let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
        let id = catalog
            .edition()
            .enzymes
            .iter()
            .find(|e| e.name == "EcoRI")
            .unwrap()
            .id
            .clone();
        let pin = collagene_core::catalog_digest::select(&catalog, &[id]).unwrap();
        let row = &mut draft.methods.rows[0];
        row.preparation =
            crate::cloning_methods_model::OverlapPreparation::OverlapExtensionPcrCatalogDigest;
        row.catalog_selection = Some(pin.clone());
        row.fragment = "unfinished index".into();
        let mut raw = crate::cloning_methods_model::OverlapExtensionPcrInput::default();
        raw.left.forward.annealing = "unfinished PCR".into();
        row.overlap_extension_pcr = Some(raw);
        draft.sync_rows();
        let (saved, entry) = save(&frame, "Edited insert then digest", &draft).unwrap();
        assert_eq!(entry.draft_version, assembly_draft_model::CURRENT_VERSION);
        let restored = load(&saved, &entry.id).unwrap();
        assert!(restored == draft);
        assert_eq!(
            restored.methods.rows[0].catalog_selection.as_ref(),
            Some(&pin)
        );
        assert_eq!(
            serde_json::to_value(&saved.documents).unwrap(),
            serde_json::to_value(&frame.documents).unwrap()
        );
        assert!(saved.candidates.is_empty());
    }

    #[test]
    fn saved_design_retains_six_oligos_without_copying_or_replacing_source_events() {
        let (frame, mut draft) = fixture();
        let original = draft.slots[0].captured.as_ref().unwrap().source.clone();
        let mut raw = crate::cloning_methods_model::OverlapExtensionPcrInput::default();
        raw.left.forward.annealing = "unfinished left".into();
        raw.left.reverse.tail = "GGTCTGAAAA".into();
        raw.right.forward.tail = "GCTACGTATT".into();
        raw.right.reverse.annealing = "TGACGATC".into();
        raw.outer.forward.tail = "independent assembly tail".into();
        raw.outer.reverse.name = "Final separate outer primer".into();
        raw.overlap = "18?".into();
        draft.methods.rows[0].preparation =
            crate::cloning_methods_model::OverlapPreparation::OverlapExtensionPcr;
        draft.methods.rows[0].overlap_extension_pcr = Some(raw);
        draft.sync_rows();
        let (saved, entry) = save(&frame, "Edited insert preparation", &draft).unwrap();
        assert_eq!(entry.draft_version, assembly_draft_model::CURRENT_VERSION);
        let restored = load(&saved, &entry.id).unwrap();
        assert!(restored == draft);
        assert!(restored.slots[0].captured.as_ref().unwrap().source == original);
        assert_eq!(
            serde_json::to_value(&frame.documents).unwrap(),
            serde_json::to_value(&saved.documents).unwrap()
        );
        let text = std::str::from_utf8(&saved.attachments[&entry.attachment_hash].bytes).unwrap();
        assert!(!text.contains("\"molecule\""));
        assert!(saved.candidates.is_empty());
    }
}
