//! Finite local slot grammar. No inferred biological ordering or assembly search.
use crate::{Error, Result, concrete, require, text};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Slot {
    pub id: String,
    pub name: String,
    pub minimum: u16,
    pub maximum: Option<u16>,
    pub left_overhang: String,
    pub right_overhang: String,
    pub left_enzyme: String,
    pub right_enzyme: String,
    pub left_addition: String,
    pub right_addition: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 64>")]
    pub extra_avoided_sites: Vec<String>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, BodyRule, 16>")]
    pub body_rules: Vec<BodyRule>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 16>")]
    pub source_ids: Vec<String>,
}
// External enum representations work in binary Postcard; no internal tags.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum BodyRule {
    RequiredPrefix(String),
    PreferredInFrameSuffix(String),
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Transition {
    pub from: String,
    pub to: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Grammar {
    pub minimum_parts: u16,
    pub maximum_parts: Option<u16>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Slot, 32>")]
    pub slots: Vec<Slot>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 32>")]
    pub starts: Vec<String>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 32>")]
    pub ends: Vec<String>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Transition, 256>")]
    pub transitions: Vec<Transition>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum GrammarIssue {
    UnknownSlot {
        index: usize,
        slot: String,
    },
    InvalidStart {
        slot: String,
    },
    InvalidTransition {
        index: usize,
        from: String,
        to: String,
    },
    CardinalityExceeded {
        slot: Option<String>,
    },
    CardinalityUnsatisfied {
        slot: Option<String>,
    },
    InvalidEnd {
        slot: String,
    },
    OverhangMismatch {
        index: usize,
        from: String,
        to: String,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Successors {
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 32>")]
    pub allowed_slots: Vec<String>,
    pub may_finish: bool,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, GrammarIssue, 1024>")]
    pub prefix_issues: Vec<GrammarIssue>,
    pub scope: String,
}
impl Grammar {
    pub fn validate(&self) -> Result<()> {
        require(
            !self.slots.is_empty()
                && self.slots.len() <= 32
                && self.transitions.len() <= 256
                && self.starts.len() <= 32
                && self.ends.len() <= 32
                && self.minimum_parts > 0,
            "grammar inventory",
        )?;
        require(
            self.maximum_parts.is_none_or(|n| n >= self.minimum_parts),
            "total cardinality",
        )?;
        let mut ids = BTreeSet::new();
        for slot in &self.slots {
            require(
                text(&slot.id, 96)
                    && ids.insert(slot.id.as_str())
                    && text(&slot.name, 128)
                    && slot.maximum.is_none_or(|n| n >= slot.minimum && n > 0)
                    && concrete(&slot.left_overhang, 16)
                    && concrete(&slot.right_overhang, 16)
                    && text(&slot.left_enzyme, 96)
                    && text(&slot.right_enzyme, 96)
                    && slot.left_addition.len() <= 256
                    && slot.left_addition.bytes().all(|b| b"ACGT".contains(&b))
                    && slot.right_addition.len() <= 256
                    && slot.right_addition.bytes().all(|b| b"ACGT".contains(&b))
                    && slot.extra_avoided_sites.len() <= 64
                    && slot.extra_avoided_sites.iter().all(|s| text(s, 96))
                    && slot.body_rules.len() <= 16
                    && slot.source_ids.len() <= 16
                    && slot.source_ids.iter().all(|s| text(s, 96)),
                "slot fields/cardinality",
            )?;
            for rule in &slot.body_rules {
                let dna = match rule {
                    BodyRule::RequiredPrefix(s) | BodyRule::PreferredInFrameSuffix(s) => s,
                };
                require(concrete(dna, 64), "slot body rule DNA")?;
            }
        }
        for endpoints in [&self.starts, &self.ends] {
            require(!endpoints.is_empty(), "grammar endpoints")?;
            let mut seen = BTreeSet::new();
            for id in endpoints {
                require(
                    ids.contains(id.as_str()) && seen.insert(id),
                    "endpoint identity",
                )?;
            }
        }
        let mut edges = BTreeSet::new();
        for edge in &self.transitions {
            require(
                ids.contains(edge.from.as_str())
                    && ids.contains(edge.to.as_str())
                    && edges.insert((&edge.from, &edge.to)),
                "transition identity",
            )?;
        }
        Ok(())
    }
    pub fn slot(&self, id: &str) -> Option<&Slot> {
        self.slots.iter().find(|s| s.id == id)
    }

    fn inspect(&self, prefix: &[String], finished: bool) -> Vec<GrammarIssue> {
        let mut issues = Vec::new();
        let mut counts = BTreeMap::<&str, usize>::new();
        if self
            .maximum_parts
            .is_some_and(|n| prefix.len() > usize::from(n))
        {
            issues.push(GrammarIssue::CardinalityExceeded { slot: None });
        }
        for (index, id) in prefix.iter().enumerate() {
            let Some(slot) = self.slot(id) else {
                issues.push(GrammarIssue::UnknownSlot {
                    index,
                    slot: id.clone(),
                });
                continue;
            };
            let count = counts.entry(id.as_str()).or_default();
            *count += 1;
            if slot.maximum.is_some_and(|n| *count > usize::from(n)) {
                issues.push(GrammarIssue::CardinalityExceeded {
                    slot: Some(id.clone()),
                });
            }
            if index == 0 && !self.starts.contains(id) {
                issues.push(GrammarIssue::InvalidStart { slot: id.clone() });
            }
            if index > 0 {
                let previous = &prefix[index - 1];
                if !self
                    .transitions
                    .iter()
                    .any(|e| e.from == *previous && e.to == *id)
                {
                    issues.push(GrammarIssue::InvalidTransition {
                        index,
                        from: previous.clone(),
                        to: id.clone(),
                    });
                }
                if let Some(left) = self.slot(previous) {
                    if left.right_overhang != slot.left_overhang {
                        issues.push(GrammarIssue::OverhangMismatch {
                            index,
                            from: previous.clone(),
                            to: id.clone(),
                        });
                    }
                }
            }
        }
        if finished {
            if prefix.len() < usize::from(self.minimum_parts) {
                issues.push(GrammarIssue::CardinalityUnsatisfied { slot: None });
            }
            for slot in &self.slots {
                if counts.get(slot.id.as_str()).copied().unwrap_or(0) < usize::from(slot.minimum) {
                    issues.push(GrammarIssue::CardinalityUnsatisfied {
                        slot: Some(slot.id.clone()),
                    });
                }
            }
            if let Some(last) = prefix.last() {
                if !self.ends.contains(last) {
                    issues.push(GrammarIssue::InvalidEnd { slot: last.clone() });
                }
            }
        }
        issues
    }
    fn prefix_preflight(&self, prefix: &[String]) -> Result<()> {
        if prefix.len() > 128 || prefix.iter().any(|id| !text(id, 96)) {
            return Err(Error::Limit("grammar prefix"));
        }
        self.validate()
    }
    pub fn violations(&self, slots: &[String]) -> Result<Vec<GrammarIssue>> {
        self.prefix_preflight(slots)?;
        Ok(self.inspect(slots, true))
    }
    /// Locally legal one-step extensions only, not a proof of future reachability,
    /// physical ligation, reading frame, biological function or reaction success.
    pub fn legal_successors(&self, prefix: &[String]) -> Result<Successors> {
        self.prefix_preflight(prefix)?;
        let prefix_issues = self.inspect(prefix, false);
        let may_finish = prefix_issues.is_empty() && self.inspect(prefix, true).is_empty();
        let mut allowed_slots = Vec::new();
        if prefix_issues.is_empty() && prefix.len() < 128 {
            for slot in &self.slots {
                let mut candidate = prefix.to_vec();
                candidate.push(slot.id.clone());
                if self.inspect(&candidate, false).is_empty() {
                    allowed_slots.push(slot.id.clone());
                }
            }
        }
        allowed_slots.sort();
        Ok(Successors { allowed_slots, may_finish, prefix_issues,
            scope: "One-step declared slot/order/cardinality/overhang grammar only; no assembly simulation or future-path enumeration".into() })
    }
}
