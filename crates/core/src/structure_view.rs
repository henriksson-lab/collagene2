//! Documentary structure viewing; no source-to-sequence or material authority.
use crate::{Error, workspace::WorkspaceAttachment};
pub use collagene_structures as coordinates;
use coordinates::{
    Atom, Format, ImportNotes,
    scene::{self, Selection},
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE: &str = "structure-view-1";
fn err(e: impl std::fmt::Display) -> Error {
    Error::Format(e.to_string())
}
fn selection_bound(s: &Option<Selection>) -> Result<(), Error> {
    if let Some(s) = s {
        let mut strings = vec![];
        match &s.chains {
            scene::Chains::Label(v) | scene::Chains::Author(Some(v)) => strings.push(v),
            _ => {}
        }
        if let scene::Alternates::Named(v) = &s.alternates {
            strings.push(v);
        }
        if s.model == 0
            || strings
                .into_iter()
                .any(|v| v.is_empty() || v.len() > 256 || v.chars().any(char::is_control))
        {
            return Err(err("Invalid structure display selection"));
        }
    }
    Ok(())
}
fn binding(format: Format, hash: &str, selection: &Option<Selection>) -> Result<String, Error> {
    selection_bound(selection)?;
    if hash.len() != 64
        || !hash
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(err("Invalid structure source hash"));
    }
    Ok(blake3::hash(
        &serde_json::to_vec(&(ENGINE, coordinates::ENGINE, format, hash, selection))
            .map_err(err)?,
    )
    .to_hex()
    .to_string())
}
pub fn revision(
    bytes: &[u8],
    format: Format,
    selection: &Option<Selection>,
) -> Result<String, Error> {
    if bytes.is_empty() || bytes.len() > coordinates::MAX_SOURCE_BYTES {
        return Err(Error::TooLarge);
    }
    binding(format, blake3::hash(bytes).to_hex().as_str(), selection)
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Prepared {
    pub engine: String,
    pub revision: String,
    pub source_hash: String,
    pub format: Format,
    pub requested_selection: Option<Selection>,
    pub atoms: Vec<Atom>,
    pub notes: ImportNotes,
    pub scene: scene::Scene,
}
impl Prepared {
    /// Transport/geometry consistency only. Original source stays with the owner;
    /// only a fresh worker import resolves its deposited connection declarations.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let bad = || err("Structure report/source/scene mismatch");
        if self.engine != ENGINE
            || self.revision != expected
            || binding(self.format, &self.source_hash, &self.requested_selection)? != expected
            || self.scene.source_hash != self.source_hash
            || self.atoms.is_empty()
            || self.atoms.len() > coordinates::MAX_ATOMS
            || self.scene.points.is_empty()
            || self.scene.points.len() > self.atoms.len()
            || self.scene.edges.len() > coordinates::MAX_CONNECTIONS
            || self.notes.connection_declarations > coordinates::MAX_CONNECTIONS
            || self.scene.unresolved_connections > self.notes.connection_declarations
            || self.scene.withheld_symmetry_or_relation > self.notes.connection_declarations
        {
            return Err(bad());
        }
        if self.notes.ignored_pdb_records.len() > 4096
            || self.notes.unprojected_cif_categories.len() > 4096
        {
            return Err(bad());
        }
        let mut ids = BTreeSet::new();
        let mut models = BTreeSet::new();
        for a in &self.atoms {
            models.insert(a.model);
            let strings = [
                Some(&a.id),
                a.label.atom.as_ref(),
                a.label.component.as_ref(),
                a.label.chain.as_ref(),
                a.label.sequence.as_ref(),
                a.author.atom.as_ref(),
                a.author.component.as_ref(),
                a.author.chain.as_ref(),
                a.author.sequence.as_ref(),
                a.entity.as_ref(),
                a.alternate.as_ref(),
                a.insertion.as_ref(),
                a.element.as_ref(),
            ];
            if a.model == 0
                || a.id.is_empty()
                || !ids.insert((a.model, &a.id))
                || strings
                    .into_iter()
                    .flatten()
                    .any(|s| s.len() > 256 || s.chars().any(char::is_control))
                || !a.position.into_iter().all(f64::is_finite)
                || a.occupancy.is_some_and(|v| !v.is_finite())
                || a.b_factor.is_some_and(|v| !v.is_finite())
            {
                return Err(bad());
            }
        }
        if models.len() > coordinates::MAX_MODELS || !models.contains(&self.scene.selection.model) {
            return Err(bad());
        }
        let expected_selection = self
            .requested_selection
            .clone()
            .unwrap_or_else(|| Selection {
                model: *models.first().unwrap(),
                chains: scene::Chains::All,
                alternates: scene::Alternates::All,
                hydrogens: true,
            });
        if self.scene.selection != expected_selection
            || !self.scene.coordinate_magnitude.is_finite()
            || self.scene.coordinate_magnitude <= 0.0
            || !self.scene.scaled_radius.is_finite()
            || self.scene.scaled_radius <= 0.0
            || !self.scene.scaled_center.into_iter().all(f64::is_finite)
        {
            return Err(bad());
        }
        if self.scene.points.len()
            != self
                .atoms
                .iter()
                .filter(|a| scene::includes(a, &self.scene.selection))
                .count()
            || self.notes.atoms_without_element
                != self.atoms.iter().filter(|a| a.element.is_none()).count()
            || self.notes.atoms_with_alternates
                != self.atoms.iter().filter(|a| a.alternate.is_some()).count()
        {
            return Err(bad());
        }
        let mut seen = BTreeSet::new();
        for p in &self.scene.points {
            let a = self.atoms.get(p.atom_index as usize).ok_or_else(bad)?;
            if !seen.insert(p.atom_index)
                || a.model != self.scene.selection.model
                || p.color
                    .into_iter()
                    .any(|v| !v.is_finite() || !(0.0..=1.0).contains(&v))
            {
                return Err(bad());
            }
            for axis in 0..3 {
                let v = ((a.position[axis] / self.scene.coordinate_magnitude
                    - self.scene.scaled_center[axis])
                    / self.scene.scaled_radius) as f32;
                if !p.position[axis].is_finite()
                    || p.position[axis].abs() > 1.00001
                    || (p.position[axis] - v).abs() > 1e-5
                {
                    return Err(bad());
                }
            }
            match &self.scene.selection.chains {
                scene::Chains::Label(c) if a.label.chain.as_ref() != Some(c) => return Err(bad()),
                scene::Chains::Author(c) if &a.author.chain != c => return Err(bad()),
                _ => {}
            }
            if !self.scene.selection.hydrogens && matches!(a.element.as_deref(), Some("H" | "D")) {
                return Err(bad());
            }
            if let scene::Alternates::Named(n) = &self.scene.selection.alternates {
                if a.alternate.as_ref().is_some_and(|v| v != n) {
                    return Err(bad());
                }
            }
        }
        let mut edges = BTreeSet::new();
        for e in &self.scene.edges {
            if e.atoms[0] >= e.atoms[1]
                || !seen.contains(&e.atoms[0])
                || !seen.contains(&e.atoms[1])
                || e.declaration_index as usize >= self.notes.connection_declarations
                || !edges.insert(e.atoms)
            {
                return Err(bad());
            }
        }
        Ok(())
    }
}
pub fn prepare(
    bytes: &[u8],
    format: Format,
    selection: Option<Selection>,
) -> Result<Prepared, Error> {
    let revision = revision(bytes, format, &selection)?;
    let source = coordinates::import(bytes, format).map_err(err)?;
    let chosen = selection.clone().unwrap_or_else(|| Selection {
        model: *source.models().first().unwrap(),
        chains: scene::Chains::All,
        alternates: scene::Alternates::All,
        hydrogens: true,
    });
    let scene = scene::build(&source, chosen).map_err(err)?;
    let out = Prepared {
        engine: ENGINE.into(),
        revision,
        source_hash: source.source_hash().into(),
        format,
        requested_selection: selection,
        atoms: source.atoms().to_vec(),
        notes: source.notes().clone(),
        scene,
    };
    out.validate_binding(&out.revision)?;
    Ok(out)
}
pub fn source_kind(format: Format) -> &'static str {
    match format {
        Format::Pdb => "structure-source-pdb",
        Format::Mmcif => "structure-source-mmcif",
    }
}
pub fn source_format(kind: &str) -> Option<Format> {
    match kind {
        "structure-source-pdb" => Some(Format::Pdb),
        "structure-source-mmcif" => Some(Format::Mmcif),
        _ => None,
    }
}
pub fn retain(
    assets: &mut BTreeMap<String, WorkspaceAttachment>,
    name: &str,
    bytes: Vec<u8>,
    format: Format,
) -> Result<String, Error> {
    revision(&bytes, format, &None)?;
    if name.is_empty() || name.len() > 256 || name.chars().any(char::is_control) {
        return Err(err("Invalid structure filename"));
    }
    let hash = blake3::hash(&bytes).to_hex().to_string();
    if let Some(a) = assets.get(&hash) {
        if a.kind != source_kind(format) || a.bytes != bytes {
            return Err(err("Existing attachment has another role"));
        }
        return Ok(hash);
    }
    if assets.len() >= crate::workspace::MAX_ATTACHMENTS
        || assets
            .values()
            .try_fold(bytes.len(), |n, a| n.checked_add(a.bytes.len()))
            .is_none_or(|n| n > crate::workspace::MAX_ARCHIVE_BYTES)
    {
        return Err(Error::TooLarge);
    }
    assets.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: source_kind(format).into(),
            name: name.into(),
            media_type: if format == Format::Pdb {
                "chemical/x-pdb"
            } else {
                "chemical/x-mmcif"
            }
            .into(),
            bytes,
        },
    );
    Ok(hash)
}
