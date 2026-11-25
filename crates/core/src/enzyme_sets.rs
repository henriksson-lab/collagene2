//! Named, ordered restriction-enzyme preferences, retained as workspace attachments.
//! These preferences never validate a digest. Missing or changed geometry requires
//! explicit review; resolving a set returns only exact supported matches.
use crate::{
    Error,
    catalog_digest::catalog,
    restriction_inspection::EnzymeSpec,
    workspace::{self, Workspace, WorkspaceAttachment},
};
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeSet,
    io::{self, Write},
};

pub const KIND: &str = "restriction-enzyme-sets-v1";
pub const MAX_SETS: usize = 128;
pub const MAX_CHOICES: usize = 1000;
pub const MAX_BYTES: usize = 256 * 1024;
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CatalogPin {
    pub edition: String,
    pub content_hash: String,
    /// Stable canonical catalog IDs, in the same order as `choices`.
    pub enzyme_ids: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Set {
    pub schema_version: u16,
    pub name: String,
    pub choices: Vec<EnzymeSpec>,
    pub catalog_revision: String,
    pub catalog_pin: Option<CatalogPin>,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Resolution {
    pub available_names: Vec<String>,
    pub missing: Vec<String>,
    pub changed: Vec<String>,
    pub catalog_changed: bool,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("enzyme sets: {s}"))
}
fn text(s: &str, max: usize) -> bool {
    !s.is_empty() && s == s.trim() && s.len() <= max && !s.chars().any(char::is_control)
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn specs_guard(choices: &[EnzymeSpec]) -> Result<(), Error> {
    if choices.len() > MAX_CHOICES {
        return Err(Error::TooLarge);
    }
    let mut names = BTreeSet::new();
    for e in choices {
        if !text(&e.name, 64)
            || !names.insert(e.name.to_ascii_lowercase())
            || e.recognition.is_empty()
            || e.recognition.len() > 64
            || !e
                .recognition
                .bytes()
                .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
            || !(-1024..=1024).contains(&e.top_cut)
            || !(-1024..=1024).contains(&e.bottom_cut)
        {
            return Err(bad("invalid or duplicate enzyme identity/geometry"));
        }
    }
    Ok(())
}
impl Set {
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != 1 || !text(&self.name, 256) || !text(&self.catalog_revision, 256)
        {
            return Err(bad("invalid version, name or catalog revision"));
        }
        specs_guard(&self.choices)?;
        if let Some(pin) = &self.catalog_pin {
            if !text(&pin.edition, 128)
                || !hex(&pin.content_hash)
                || pin.enzyme_ids.len() != self.choices.len()
            {
                return Err(bad("invalid catalog pin"));
            }
            let mut ids = BTreeSet::new();
            if pin
                .enzyme_ids
                .iter()
                .any(|id| !text(id, 96) || !ids.insert(id))
            {
                return Err(bad("invalid or duplicate stable enzyme ID"));
            }
            if self.catalog_revision != format!("{}:{}", pin.edition, pin.content_hash) {
                return Err(bad("catalog revision/pin mismatch"));
            }
        }
        Ok(())
    }
}
/// Empty sets are valid: applying one explicitly clears the selected enzymes.
pub fn build(name: &str, choices: &[EnzymeSpec], catalog_revision: &str) -> Result<Set, Error> {
    if !text(name.trim(), 256) || !text(catalog_revision, 256) {
        return Err(bad("invalid name or catalog revision"));
    }
    specs_guard(choices)?;
    Ok(Set {
        schema_version: 1,
        name: name.trim().into(),
        choices: choices.to_vec(),
        catalog_revision: catalog_revision.into(),
        catalog_pin: None,
    })
}
pub fn catalog_revision(c: &catalog::Catalog) -> String {
    format!("{}:{}", c.edition().edition, c.content_hash())
}
fn spec(e: &catalog::Enzyme) -> EnzymeSpec {
    EnzymeSpec {
        name: e.name.clone(),
        recognition: e.geometry.recognition.clone(),
        top_cut: e.geometry.top_cut,
        bottom_cut: e.geometry.bottom_cut,
    }
}
/// Prefer this constructor for catalog-backed choices: retain canonical IDs as
/// well as the exact edition, content hash, name and cleavage geometry.
pub fn build_for_catalog(
    name: &str,
    choices: &[EnzymeSpec],
    c: &catalog::Catalog,
) -> Result<Set, Error> {
    let mut set = build(name, choices, &catalog_revision(c))?;
    let mut ids = Vec::new();
    for choice in choices {
        let record = c
            .edition()
            .enzymes
            .iter()
            .find(|e| spec(e) == *choice)
            .ok_or_else(|| bad("choice is not exact catalog identity/geometry"))?;
        ids.push(record.id.clone());
    }
    set.catalog_pin = Some(CatalogPin {
        edition: c.edition().edition.clone(),
        content_hash: c.content_hash().into(),
        enzyme_ids: ids,
    });
    set.validate()?;
    Ok(set)
}
/// Geometry/name matching for custom catalogs. Pinned sets require the catalog
/// resolver so an ID disappearance or reassignment cannot be hidden by a name.
pub fn resolve(set: &Set, current: &[EnzymeSpec], revision: &str) -> Result<Resolution, Error> {
    set.validate()?;
    specs_guard(current)?;
    if !text(revision, 256) {
        return Err(bad("invalid current catalog revision"));
    }
    if set.catalog_pin.is_some() {
        return Err(bad("use resolve_for_catalog for stable-ID pinned sets"));
    }
    Ok(resolve_matches(set, revision, |i| {
        current
            .iter()
            .find(|e| e.name == set.choices[i].name)
            .cloned()
    }))
}
fn resolve_matches(
    set: &Set,
    revision: &str,
    mut current: impl FnMut(usize) -> Option<EnzymeSpec>,
) -> Resolution {
    let mut result = Resolution {
        available_names: Vec::new(),
        missing: Vec::new(),
        changed: Vec::new(),
        catalog_changed: set.catalog_revision != revision,
    };
    for (i, choice) in set.choices.iter().enumerate() {
        match current(i) {
            Some(e) if e == *choice => result.available_names.push(choice.name.clone()),
            Some(_) => result.changed.push(choice.name.clone()),
            None => result.missing.push(choice.name.clone()),
        }
    }
    result
}
pub fn resolve_for_catalog(set: &Set, c: &catalog::Catalog) -> Result<Resolution, Error> {
    set.validate()?;
    Ok(resolve_matches(set, &catalog_revision(c), |i| {
        let record = if let Some(pin) = &set.catalog_pin {
            c.edition()
                .enzymes
                .iter()
                .find(|e| e.id == pin.enzyme_ids[i])
        } else {
            c.edition()
                .enzymes
                .iter()
                .find(|e| e.name == set.choices[i].name)
        };
        record.map(spec)
    }))
}
// Reject oversized native DTOs before workspace validation or cloning.
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
    for n in frame
        .revisions
        .values()
        .map(|m| m.sequence.len())
        .chain(frame.attachments.values().map(|a| a.bytes.len()))
    {
        if n > workspace::MAX_ENTRY_BYTES {
            return Err(Error::TooLarge);
        }
        bytes = bytes
            .checked_add(n)
            .filter(|n| *n <= workspace::MAX_EXPANDED_BYTES)
            .ok_or(Error::TooLarge)?;
    }
    for m in frame.revisions.values() {
        crate::model::bounds::model_preflight(&[m])?;
    }
    for c in &frame.candidates {
        if c.plan.inputs.len() > 1000 || c.plan.steps.len() > 1000 {
            return Err(Error::TooLarge);
        }
        for i in &c.plan.inputs {
            crate::model::bounds::model_preflight(&[&i.molecule])?;
        }
    }
    frame.validate()
}
struct Limited(Vec<u8>);
impl Write for Limited {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        if self.0.len().saturating_add(b.len()) > MAX_BYTES {
            return Err(io::Error::other("enzyme set byte limit"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn load_guarded(frame: &Workspace, hash: &str) -> Result<Set, Error> {
    if !hex(hash) {
        return Err(bad("invalid attachment identity"));
    }
    let a = frame
        .attachments
        .get(hash)
        .ok_or_else(|| bad("set is unavailable"))?;
    if a.kind != KIND
        || a.media_type != "application/json"
        || a.bytes.len() > MAX_BYTES
        || blake3::hash(&a.bytes).to_hex().as_str() != hash
    {
        return Err(bad("invalid set attachment type, size or hash"));
    }
    let set: Set = serde_json::from_slice(&a.bytes).map_err(|_| bad("invalid set JSON"))?;
    set.validate()?;
    if a.name != set.name {
        return Err(bad("attachment/set name mismatch"));
    }
    Ok(set)
}
fn list_guarded(frame: &Workspace) -> Result<Vec<(String, Set)>, Error> {
    let mut sets = Vec::new();
    let mut names = BTreeSet::new();
    for (hash, a) in &frame.attachments {
        if a.kind != KIND {
            continue;
        }
        if sets.len() >= MAX_SETS {
            return Err(Error::TooLarge);
        }
        let set = load_guarded(frame, hash)?;
        if !names.insert(set.name.to_ascii_lowercase()) {
            return Err(bad("duplicate named sets"));
        }
        sets.push((hash.clone(), set));
    }
    sets.sort_by(|a, b| a.1.name.cmp(&b.1.name));
    Ok(sets)
}
pub fn list(frame: &Workspace) -> Result<Vec<(String, Set)>, Error> {
    workspace_guard(frame)?;
    list_guarded(frame)
}
pub fn load(frame: &Workspace, hash: &str) -> Result<Set, Error> {
    workspace_guard(frame)?;
    list_guarded(frame)?;
    load_guarded(frame, hash)
}
/// Replace the same case-insensitive name atomically, retaining order and pinned
/// geometry. Only attachments change; documents, history and candidates do not.
pub fn store(frame: &Workspace, set: &Set) -> Result<(Workspace, String), Error> {
    set.validate()?;
    workspace_guard(frame)?;
    let sets = list_guarded(frame)?;
    let replaced = sets
        .iter()
        .find(|(_, s)| s.name.eq_ignore_ascii_case(&set.name))
        .map(|(h, _)| h.clone());
    if replaced.is_none() && sets.len() >= MAX_SETS {
        return Err(Error::TooLarge);
    }
    let mut encoded = Limited(Vec::new());
    serde_json::to_writer(&mut encoded, set).map_err(|_| Error::TooLarge)?;
    let hash = blake3::hash(&encoded.0).to_hex().to_string();
    if let Some(a) = frame.attachments.get(&hash) {
        if a.kind != KIND
            || a.name != set.name
            || a.media_type != "application/json"
            || a.bytes != encoded.0
        {
            return Err(bad(
                "attachment identity already has a different type or metadata",
            ));
        }
    }
    if !frame.attachments.contains_key(&hash)
        && replaced.is_none()
        && frame.attachments.len() >= workspace::MAX_ATTACHMENTS
    {
        return Err(Error::TooLarge);
    }
    let mut staged = frame.clone();
    if let Some(old) = replaced {
        staged.attachments.remove(&old);
    }
    staged.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: KIND.into(),
            name: set.name.clone(),
            media_type: "application/json".into(),
            bytes: encoded.0,
        },
    );
    staged.validate()?;
    Ok((staged, hash))
}
/// Delete only a validated named-set attachment, by exact content identity.
pub fn delete(frame: &Workspace, hash: &str) -> Result<Workspace, Error> {
    workspace_guard(frame)?;
    list_guarded(frame)?;
    load_guarded(frame, hash)?;
    let mut staged = frame.clone();
    staged.attachments.remove(hash);
    staged.validate()?;
    Ok(staged)
}
