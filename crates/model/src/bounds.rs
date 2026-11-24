//! Borrowed native-model guards. No serialization, cloning, recursive traversal or
//! unbounded child stack is used before these aggregate limits have been checked.
use crate::{Error, Feature, Location, Molecule};
pub const MAX_DEPTH: usize = 32;
#[derive(Default)]
pub struct ModelBudget {
    features: usize,
    nodes: usize,
    qualifiers: usize,
    metadata: usize,
    text: usize,
}
impl ModelBudget {
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
            self.text(s, 4 * 1024 * 1024)?;
        }
        self.features(&m.features)
    }
    pub fn features(&mut self, features: &[Feature]) -> Result<(), Error> {
        Self::add(&mut self.features, features.len(), 4096)?;
        for feature in features {
            self.text(&feature.kind, 256)?;
            Self::add(&mut self.qualifiers, feature.qualifiers.len(), 16384)?;
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
                Self::add(&mut self.nodes, 1, 16384)?;
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
    if molecules.len() > 1024 {
        return Err(Error::TooLarge);
    }
    let mut budget = ModelBudget::default();
    for molecule in molecules {
        budget.molecule(molecule)?;
    }
    Ok(())
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
