//! Borrowed native-model guards. No serialization, cloning, recursive traversal or
//! unbounded child stack is used before these aggregate limits have been checked.
use crate::{Error, Feature, Location, Molecule};
pub const MAX_DEPTH: usize = 32;
#[derive(Clone, Copy)]
struct Limits {
    features: usize,
    nodes: usize,
    qualifiers: usize,
    original_bytes: usize,
    original_shares_text: bool,
}
impl Limits {
    const MODEL: Self = Self {
        features: 4096,
        nodes: 16384,
        qualifiers: 16384,
        original_bytes: 4 * 1024 * 1024,
        original_shares_text: true,
    };
    const ANNOTATION_EDIT: Self = Self {
        features: 100_000,
        nodes: 1_000_000,
        qualifiers: 1_000_000,
        original_bytes: 32 * 1024 * 1024,
        original_shares_text: false,
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
    fn text(&mut self, s: &str, cap: usize) -> Result<(), Error> {
        if s.len() > cap {
            return Err(Error::TooLarge);
        }
        Self::add(&mut self.text, s.len(), 8 * 1024 * 1024)
    }
    pub fn molecule(&mut self, m: &Molecule) -> Result<(), Error> {
        self.text(&m.name, 1024)?;
        Self::add(&mut self.metadata, m.metadata.len(), 4096)?;
        for (k, v) in &m.metadata {
            self.text(k, 4096)?;
            self.text(v, 1024 * 1024)?;
        }
        if let Some(s) = &m.original_genbank {
            if self.limits.original_shares_text {
                self.text(s, self.limits.original_bytes)?;
            } else {
                Self::add(&mut self.original_text, s.len(), self.limits.original_bytes)?;
            }
        }
        self.features(&m.features)
    }
    pub fn features(&mut self, features: &[Feature]) -> Result<(), Error> {
        Self::add(&mut self.features, features.len(), self.limits.features)?;
        for feature in features {
            self.text(&feature.kind, 256)?;
            if feature.qualifiers.len() > 16384 {
                return Err(Error::TooLarge);
            }
            Self::add(
                &mut self.qualifiers,
                feature.qualifiers.len(),
                self.limits.qualifiers,
            )?;
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
                Self::add(&mut self.nodes, 1, self.limits.nodes)?;
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
