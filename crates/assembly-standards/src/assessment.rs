//! Part-format examiner with explicit sequence, output and comparison budgets.
use crate::grammar::{BodyRule, GrammarIssue};
use crate::profile::{AvoidedSite, Profile, ProfileKind, Requirement, Source};
use crate::{Error, MAX_REPORT_BYTES, Result, bounded_json, mask, require, text};
use collagene_model::{Alphabet, Sequence, complement};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Part {
    pub id: String,
    pub slot: String,
    /// Oriented formatted linear part; not a whole vector. No origin is inferred.
    pub sequence: Sequence,
    /// Zero-based half-open body interval, excluding declared format additions.
    pub body_start: usize,
    pub body_end: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Limits {
    pub max_parts: usize,
    pub max_bases: usize,
    pub max_comparisons: usize,
    pub max_findings: usize,
    pub max_report_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_parts: 64,
            max_bases: 1_000_000,
            max_comparisons: 8_000_000,
            max_findings: 1024,
            max_report_bytes: 512 * 1024,
        }
    }
}
impl Limits {
    pub fn validate(&self) -> Result<()> {
        require(
            (1..=128).contains(&self.max_parts)
                && (1..=2_000_000).contains(&self.max_bases)
                && (1..=64_000_000).contains(&self.max_comparisons)
                && (1..=4096).contains(&self.max_findings)
                && (4096..=MAX_REPORT_BYTES).contains(&self.max_report_bytes),
            "assessment limits",
        )
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Severity {
    Violation,
    Indeterminate,
    Advisory,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Certainty {
    Certain,
    Possible,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Side {
    Left,
    Right,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum FindingKind {
    Grammar(GrammarIssue),
    Addition {
        side: Side,
        expected: String,
        certainty: Option<Certainty>,
    },
    BodyRule {
        rule: BodyRule,
        certainty: Option<Certainty>,
    },
    AvoidedSite {
        rule_id: String,
        start: usize,
        end: usize,
        reverse: bool,
        certainty: Certainty,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Finding {
    pub part_index: Option<usize>,
    pub severity: Severity,
    pub kind: FindingKind,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum BudgetHit {
    Comparisons,
    Findings,
    ReportBytes,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Outcome {
    Conforms,
    Violations,
    Indeterminate,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PartBinding {
    pub id: String,
    pub slot: String,
    pub sequence_blake3: String,
    pub sequence_length: usize,
    pub body_start: usize,
    pub body_end: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Assessment {
    pub engine: String,
    pub profile_identity: String,
    pub profile_id: String,
    pub profile_revision: String,
    pub profile_kind: ProfileKind,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, 16>")]
    pub source_evidence: Vec<Source>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 32>")]
    pub limitations: Vec<String>,
    pub catalog_edition: String,
    pub catalog_content_hash: String,
    pub scope: String,
    pub authority: String,
    pub request_revision: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, PartBinding, 128>")]
    pub parts: Vec<PartBinding>,
    pub limits: Limits,
    pub outcome: Outcome,
    pub complete: bool,
    pub comparisons: usize,
    pub violations_observed: usize,
    pub indeterminate_observed: usize,
    pub advisories_observed: usize,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Finding, 4096>")]
    pub findings: Vec<Finding>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, BudgetHit, 3>")]
    pub budgets_hit: Vec<BudgetHit>,
}
impl Assessment {
    pub fn identity(&self) -> Result<String> {
        Ok(blake3::hash(&bounded_json(
            self,
            self.limits.max_report_bytes.min(MAX_REPORT_BYTES),
        )?)
        .to_hex()
        .to_string())
    }
    /// Fresh replay is required; a self-reported checksum is not evidence.
    pub fn validate_for(&self, profile: &Profile, parts: &[Part], limits: &Limits) -> Result<()> {
        let fresh = profile.assess(parts, limits)?;
        if &fresh != self {
            return Err(Error::Invalid("assessment source/profile/options changed"));
        }
        Ok(())
    }
}
pub(crate) struct Examiner {
    report: Assessment,
}
impl Examiner {
    pub(crate) fn assess(profile: &Profile, parts: &[Part], limits: &Limits) -> Result<Assessment> {
        limits.validate()?;
        if parts.len() > limits.max_parts {
            return Err(Error::Limit("part count"));
        }
        let mut ids = BTreeSet::new();
        let mut bases = 0_usize;
        for part in parts {
            require(
                text(&part.id, 96)
                    && ids.insert(part.id.as_str())
                    && text(&part.slot, 96)
                    && part.sequence.alphabet() == Alphabet::Dna
                    && part.body_start < part.body_end
                    && part.body_end <= part.sequence.len(),
                "part identity/alphabet/body interval",
            )?;
            bases = bases
                .checked_add(part.sequence.len())
                .ok_or(Error::Limit("base count overflow"))?;
            if bases > limits.max_bases {
                return Err(Error::Limit("aggregate source bases"));
            }
        }
        let bindings = parts
            .iter()
            .map(|p| PartBinding {
                id: p.id.clone(),
                slot: p.slot.clone(),
                sequence_blake3: blake3::hash(p.sequence.as_str().as_bytes())
                    .to_hex()
                    .to_string(),
                sequence_length: p.sequence.len(),
                body_start: p.body_start,
                body_end: p.body_end,
            })
            .collect::<Vec<_>>();
        let request = bounded_json(&(profile.identity(), &bindings, limits), 128 * 1024)?;
        let mut examiner = Self { report: Assessment {
            engine: crate::ENGINE_VERSION.into(), profile_identity: profile.identity().into(),
            profile_id: profile.definition().id.clone(), profile_revision: profile.definition().revision.clone(),
            profile_kind: profile.kind(),
            source_evidence: profile.definition().sources.clone(), limitations: profile.definition().limitations.clone(),
            catalog_edition: profile.catalog_pin().edition.clone(), catalog_content_hash: profile.catalog_pin().content_hash.clone(),
            scope: profile.definition().scope.clone(), request_revision: blake3::hash(&request).to_hex().to_string(),
            authority: match profile.kind() {
                ProfileKind::Builtin => "Exact installed part-sequence profile only; citations and hashes are documentary provenance, not experimental verification or complete-standard compliance",
                ProfileKind::Custom => "Declared custom grammar; source URIs, labels and asserted standard names do not authenticate its rules as a supported built-in standard",
            }.into(),
            parts: bindings, limits: limits.clone(), outcome: Outcome::Indeterminate, complete: true,
            comparisons: 0, violations_observed: 0, indeterminate_observed: 0, advisories_observed: 0,
            findings: Vec::new(), budgets_hit: Vec::new(),
        }};
        let slots = parts.iter().map(|p| p.slot.clone()).collect::<Vec<_>>();
        for issue in profile.definition().grammar.violations(&slots)? {
            examiner.add(Finding {
                part_index: None,
                severity: Severity::Violation,
                kind: FindingKind::Grammar(issue),
            });
            if !examiner.report.complete {
                break;
            }
        }
        for (index, part) in parts.iter().enumerate() {
            if !examiner.report.complete {
                break;
            }
            let Some(slot) = profile.definition().grammar.slot(&part.slot) else {
                continue;
            };
            let dna = part.sequence.as_str().as_bytes();
            examiner.addition(
                index,
                Side::Left,
                &dna[..part.body_start],
                &slot.left_addition,
            );
            examiner.addition(
                index,
                Side::Right,
                &dna[part.body_end..],
                &slot.right_addition,
            );
            let body = &dna[part.body_start..part.body_end];
            for rule in &slot.body_rules {
                if !examiner.report.complete {
                    break;
                }
                examiner.body_rule(index, body, rule);
            }
            let mut seen = BTreeSet::new();
            for id in profile
                .definition()
                .default_avoided_sites
                .iter()
                .chain(&slot.extra_avoided_sites)
            {
                if !examiner.report.complete {
                    break;
                }
                if !seen.insert(id.as_str()) {
                    continue;
                }
                let site = profile
                    .definition()
                    .avoided_sites
                    .iter()
                    .find(|s| &s.id == id)
                    .ok_or(Error::Invalid("missing installed site rule"))?;
                examiner.site(index, body, part.body_start, site);
            }
        }
        examiner.finish()
    }
    fn limit(&mut self, reason: BudgetHit) {
        self.report.complete = false;
        if !self.report.budgets_hit.contains(&reason) {
            self.report.budgets_hit.push(reason);
        }
    }
    fn add(&mut self, finding: Finding) {
        match finding.severity {
            Severity::Violation => self.report.violations_observed += 1,
            Severity::Indeterminate => self.report.indeterminate_observed += 1,
            Severity::Advisory => self.report.advisories_observed += 1,
        }
        if self.report.findings.len() == self.report.limits.max_findings {
            self.limit(BudgetHit::Findings);
        } else {
            self.report.findings.push(finding);
        }
    }
    fn matched(&mut self, observed: &[u8], expected: &[u8]) -> Option<Certainty> {
        if observed.len() != expected.len() {
            return None;
        }
        let mut certainty = Certainty::Certain;
        for (&a, &b) in observed.iter().zip(expected) {
            if self.report.comparisons == self.report.limits.max_comparisons {
                self.limit(BudgetHit::Comparisons);
                return None;
            }
            self.report.comparisons += 1;
            let (a, b) = (mask(a), mask(b));
            if a & b == 0 {
                return None;
            }
            if a & !b != 0 {
                certainty = Certainty::Possible;
            }
        }
        Some(certainty)
    }
    fn addition(&mut self, index: usize, side: Side, observed: &[u8], expected: &str) {
        if !self.report.complete {
            return;
        }
        let certainty = self.matched(observed, expected.as_bytes());
        if !self.report.complete || certainty == Some(Certainty::Certain) {
            return;
        }
        self.add(Finding {
            part_index: Some(index),
            severity: if certainty.is_some() {
                Severity::Indeterminate
            } else {
                Severity::Violation
            },
            kind: FindingKind::Addition {
                side,
                expected: expected.into(),
                certainty,
            },
        });
    }
    fn body_rule(&mut self, index: usize, body: &[u8], rule: &BodyRule) {
        let (observed, expected, advisory) = match rule {
            BodyRule::RequiredPrefix(dna) => (body.get(..dna.len()), dna.as_bytes(), false),
            BodyRule::PreferredInFrameSuffix(dna) => (
                body.len()
                    .is_multiple_of(3)
                    .then(|| body.get(body.len().saturating_sub(dna.len())..))
                    .flatten(),
                dna.as_bytes(),
                true,
            ),
        };
        let certainty = observed.and_then(|s| self.matched(s, expected));
        if !self.report.complete || certainty == Some(Certainty::Certain) {
            return;
        }
        let severity = if advisory {
            Severity::Advisory
        } else if certainty.is_some() {
            Severity::Indeterminate
        } else {
            Severity::Violation
        };
        self.add(Finding {
            part_index: Some(index),
            severity,
            kind: FindingKind::BodyRule {
                rule: rule.clone(),
                certainty,
            },
        });
    }
    fn site(&mut self, index: usize, body: &[u8], offset: usize, site: &AvoidedSite) {
        if body.len() < site.recognition.len() {
            return;
        }
        let forward = site.recognition.as_bytes();
        let reverse = forward
            .iter()
            .rev()
            .map(|&b| complement(b, Alphabet::Dna))
            .collect::<Vec<_>>();
        for (reversed, motif) in [(false, forward), (true, reverse.as_slice())] {
            if reversed && reverse.as_slice() == forward {
                continue;
            }
            for start in 0..=body.len() - motif.len() {
                let certainty = self.matched(&body[start..start + motif.len()], motif);
                if !self.report.complete {
                    return;
                }
                if let Some(certainty) = certainty {
                    let severity = match (site.requirement, certainty) {
                        (Requirement::Required, Certainty::Certain) => Severity::Violation,
                        (Requirement::Required, Certainty::Possible) => Severity::Indeterminate,
                        (Requirement::Advisory, _) => Severity::Advisory,
                    };
                    self.add(Finding {
                        part_index: Some(index),
                        severity,
                        kind: FindingKind::AvoidedSite {
                            rule_id: site.id.clone(),
                            start: offset + start,
                            end: offset + start + motif.len(),
                            reverse: reversed,
                            certainty,
                        },
                    });
                    if !self.report.complete {
                        return;
                    }
                }
            }
        }
    }
    fn finish(mut self) -> Result<Assessment> {
        loop {
            self.report.outcome = if self.report.violations_observed > 0 {
                Outcome::Violations
            } else if !self.report.complete || self.report.indeterminate_observed > 0 {
                Outcome::Indeterminate
            } else {
                Outcome::Conforms
            };
            if bounded_json(&self.report, self.report.limits.max_report_bytes).is_ok() {
                return Ok(self.report);
            }
            self.limit(BudgetHit::ReportBytes);
            if self.report.findings.pop().is_none() {
                return Err(Error::Limit("report metadata exceeds byte budget"));
            }
        }
    }
}
