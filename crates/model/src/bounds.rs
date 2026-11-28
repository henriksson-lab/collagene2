//! Borrowed native-model guards. No serialization, cloning, recursive traversal or
//! unbounded child stack is used before these aggregate limits have been checked.
use crate::{Error, Feature, Location, Molecule};
pub const MAX_DEPTH: usize = 32;
#[derive(Clone, Copy)]
struct Limits {
    features: usize,
    nodes: usize,
    qualifiers: usize,
    metadata: usize,
    text_bytes: usize,
    original_bytes: usize,
    original_shares_text: bool,
    retained_metadata_original: bool,
}
impl Limits {
    const MODEL: Self = Self {
        features: 4096,
        nodes: 16384,
        qualifiers: 16384,
        metadata: 4096,
        text_bytes: 8 * 1024 * 1024,
        original_bytes: 4 * 1024 * 1024,
        original_shares_text: true,
        retained_metadata_original: false,
    };
    const ANNOTATION_EDIT: Self = Self {
        features: 100_000,
        nodes: 1_000_000,
        qualifiers: 1_000_000,
        metadata: 4096,
        text_bytes: 8 * 1024 * 1024,
        original_bytes: 32 * 1024 * 1024,
        original_shares_text: false,
        retained_metadata_original: true,
    };
}
pub struct ModelBudget {
    features: usize,
    nodes: usize,
    qualifiers: usize,
    metadata: usize,
    text: usize,
    original_text: usize,
    limits: Limits,
}
impl Default for ModelBudget {
    fn default() -> Self {
        Self::with_limits(Limits::MODEL)
    }
}
impl ModelBudget {
    fn with_limits(limits: Limits) -> Self {
        Self {
            features: 0,
            nodes: 0,
            qualifiers: 0,
            metadata: 0,
            text: 0,
            original_text: 0,
            limits,
        }
    }
    /// Native reference edits preserve import-scale annotation/original-record
    /// support. Annotation/metadata text retains its independent 8 MiB budget;
    /// exact original records have a separate aggregate 32 MiB allowance.
    pub fn annotation_edit() -> Self {
        Self::with_limits(Limits::ANNOTATION_EDIT)
    }
    /// Aggregate allocation limits for one admitted native envelope, distinct
    /// from the import-scale limits of each individual molecule.
    pub(crate) fn native_frame(maximum_bytes: usize) -> Self {
        Self::with_limits(Limits {
            features: 1_000_000,
            nodes: 1_000_000,
            qualifiers: 1_000_000,
            metadata: 1_000_000,
            text_bytes: maximum_bytes,
            original_bytes: maximum_bytes,
            original_shares_text: true,
            retained_metadata_original: true,
        })
    }
    pub(crate) fn retained_metadata(key: &str) -> bool {
        matches!(
            key,
            "trace_original_base64"
                | "source:snapgene:original_base64"
                | "legacy:original_base64"
                | "fastq:original_base64"
                | "fastq:quality_base64"
        ) || key
            .strip_prefix("annotation_import:")
            .and_then(|s| s.strip_suffix(":original_base64"))
            .is_some_and(|s| {
                s.len() == 64
                    && s.bytes()
                        .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
            })
    }
    pub(crate) fn metadata_value(&mut self, key: &str, value: &str) -> Result<(), Error> {
        if self.limits.retained_metadata_original && Self::retained_metadata(key) {
            self.original(value)
        } else {
            self.text(value, 1024 * 1024)
        }
    }
    fn models(&mut self, molecules: &[&Molecule]) -> Result<(), Error> {
        if molecules.len() > 1024 {
            return Err(Error::TooLarge);
        }
        for molecule in molecules {
            self.molecule(molecule)?;
        }
        Ok(())
    }

    fn add(value: &mut usize, n: usize, max: usize) -> Result<(), Error> {
        *value = value
            .checked_add(n)
            .filter(|v| *v <= max)
            .ok_or(Error::TooLarge)?;
        Ok(())
    }
    pub(crate) fn text(&mut self, s: &str, cap: usize) -> Result<(), Error> {
        if s.len() > cap {
            return Err(Error::TooLarge);
        }
        Self::add(&mut self.text, s.len(), self.limits.text_bytes)
    }
    pub(crate) fn feature_count(&mut self, count: usize) -> Result<(), Error> {
        Self::add(&mut self.features, count, self.limits.features)
    }
    pub(crate) fn qualifier_count(&mut self, count: usize) -> Result<(), Error> {
        Self::add(&mut self.qualifiers, count, self.limits.qualifiers)
    }
    pub(crate) fn node_count(&mut self, count: usize) -> Result<(), Error> {
        Self::add(&mut self.nodes, count, self.limits.nodes)
    }
    pub(crate) fn metadata_count(&mut self, count: usize) -> Result<(), Error> {
        Self::add(&mut self.metadata, count, self.limits.metadata)
    }
    pub(crate) fn remaining_features(&self) -> usize {
        self.limits.features - self.features
    }
    pub(crate) fn remaining_qualifiers(&self) -> usize {
        self.limits.qualifiers - self.qualifiers
    }
    pub(crate) fn remaining_nodes(&self) -> usize {
        self.limits.nodes - self.nodes
    }
    pub(crate) fn remaining_metadata(&self) -> usize {
        self.limits.metadata - self.metadata
    }
    pub(crate) fn original(&mut self, text: &str) -> Result<(), Error> {
        if self.limits.original_shares_text {
            self.text(text, self.limits.original_bytes)
        } else {
            Self::add(
                &mut self.original_text,
                text.len(),
                self.limits.original_bytes,
            )
        }
    }
    pub fn molecule(&mut self, m: &Molecule) -> Result<(), Error> {
        self.text(&m.name, 1024)?;
        self.metadata_count(m.metadata.len())?;
        for (k, v) in &m.metadata {
            self.text(k, 4096)?;
            self.metadata_value(k, v)?;
        }
        if let Some(s) = &m.original_genbank {
            self.original(s)?;
        }
        self.features(&m.features)
    }
    pub fn features(&mut self, features: &[Feature]) -> Result<(), Error> {
        self.feature_count(features.len())?;
        for feature in features {
            self.text(&feature.kind, 256)?;
            if feature.qualifiers.len() > 16384 {
                return Err(Error::TooLarge);
            }
            self.qualifier_count(feature.qualifiers.len())?;
            for (k, v) in &feature.qualifiers {
                self.text(k, 4096)?;
                if let Some(v) = v {
                    self.text(v, 1024 * 1024)?;
                }
            }
            self.location(&feature.location)?;
        }
        Ok(())
    }
    fn location(&mut self, root: &Location) -> Result<(), Error> {
        struct Frame<'a> {
            node: &'a Location,
            next: usize,
            entered: bool,
        }
        let mut stack = Vec::with_capacity(MAX_DEPTH + 1);
        stack.push(Frame {
            node: root,
            next: 0,
            entered: false,
        });
        while let Some(frame) = stack.last_mut() {
            if !frame.entered {
                frame.entered = true;
                self.node_count(1)?;
                match frame.node {
                    Location::Join(parts) | Location::Order(parts) => {
                        if parts.len() > 16384 {
                            return Err(Error::TooLarge);
                        }
                    }
                    Location::Unresolved(s) => self.text(s, 65536)?,
                    _ => {}
                }
            }
            let child = match frame.node {
                Location::Complement(child) if frame.next == 0 => {
                    frame.next = 1;
                    Some(child.as_ref())
                }
                Location::Join(parts) | Location::Order(parts) if frame.next < parts.len() => {
                    let child = &parts[frame.next];
                    frame.next += 1;
                    Some(child)
                }
                _ => None,
            };
            if let Some(node) = child {
                if stack.len() == MAX_DEPTH {
                    return Err(Error::TooLarge);
                }
                stack.push(Frame {
                    node,
                    next: 0,
                    entered: false,
                });
            } else {
                stack.pop();
            }
        }
        Ok(())
    }
}

/// Inspect borrowed annotation trees and strings before cloning, recursive validation,
/// serialization or hashing. This is a resource guard, not molecular validation.
/// All supplied models share one aggregate budget.
pub fn model_preflight(molecules: &[&Molecule]) -> Result<(), Error> {
    ModelBudget::default().models(molecules)
}

/// Borrowed admission for native annotation edits. Original source records use
/// a separate 32 MiB aggregate budget; feature/text/depth checks run through the
/// same implementation as ordinary native-model admission.
pub fn annotation_edit_preflight(molecules: &[&Molecule]) -> Result<(), Error> {
    ModelBudget::annotation_edit().models(molecules)
}

/// Guard retained original feature evidence without copying it into a dummy model.
/// Groups share one aggregate feature/tree/string budget. This grants no annotation
/// mapping or molecular validity to the evidence.
pub fn feature_preflight(groups: &[&[Feature]]) -> Result<(), Error> {
    if groups.len() > 1024 {
        return Err(Error::TooLarge);
    }
    let mut budget = ModelBudget::default();
    for features in groups {
        budget.features(features)?;
    }
    Ok(())
}

/// Borrowed admission matching native decoding: each molecule has its own import-
/// scale budget, while all models share one enclosing envelope allocation budget.
/// Checks precede revision hashing or serialization and do not clone model data.
pub fn native_frame_preflight<'a>(
    molecules: impl IntoIterator<Item = &'a Molecule>,
    maximum_bytes: usize,
) -> Result<(), Error> {
    let mut frame = ModelBudget::native_frame(maximum_bytes);
    for molecule in molecules {
        ModelBudget::annotation_edit().molecule(molecule)?;
        frame.molecule(molecule)?;
    }
    Ok(())
}
