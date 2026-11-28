//! The single current typed contribution qualifier; imported claims are documentary.
use super::*;
use crate::{formats::annotations::AnnotationFormat, material::lineage::EventBinding};
use serde::{Deserialize, Serialize};
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum SourcePin {
    Library {
        library_revision: String,
        library_name: String,
        library_version: String,
        entry_id: String,
        entry_revision: String,
        source_revision: String,
        source_feature_index: usize,
        proposal_id: String,
        concrete_matches: usize,
        compatible_ambiguous: usize,
        mismatches: usize,
        template_length: usize,
    },
    Dataset {
        source_blake3: String,
        format: AnnotationFormat,
        seqid: String,
        proposal_id: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 10000>")]
        source_rows: Vec<usize>,
    },
    Donor {
        event: EventBinding,
        feature_index: usize,
        annotation_revision: String,
        transform: Transform,
    },
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Contribution {
    pub source: SourcePin,
    /// Proposed annotation content excluding this qualifier and external row evidence.
    pub feature_revision: String,
}
impl Contribution {
    /// Stable source evidence identity. Proposal IDs retain the first declared
    /// application context but do not create another source on target revision.
    /// Full PartialEq remains the exact documentary/report content comparison.
    pub fn same_evidence(&self, other: &Self) -> bool {
        if self.feature_revision != other.feature_revision {
            return false;
        }
        match (&self.source, &other.source) {
            (
                SourcePin::Library {
                    library_revision: a_library,
                    library_name: a_name,
                    library_version: a_version,
                    entry_id: a_entry,
                    entry_revision: a_entry_revision,
                    source_revision: a_source,
                    source_feature_index: a_feature,
                    proposal_id: _,
                    concrete_matches: a_concrete,
                    compatible_ambiguous: a_ambiguous,
                    mismatches: a_mismatches,
                    template_length: a_length,
                },
                SourcePin::Library {
                    library_revision: b_library,
                    library_name: b_name,
                    library_version: b_version,
                    entry_id: b_entry,
                    entry_revision: b_entry_revision,
                    source_revision: b_source,
                    source_feature_index: b_feature,
                    proposal_id: _,
                    concrete_matches: b_concrete,
                    compatible_ambiguous: b_ambiguous,
                    mismatches: b_mismatches,
                    template_length: b_length,
                },
            ) => {
                (
                    a_library,
                    a_name,
                    a_version,
                    a_entry,
                    a_entry_revision,
                    a_source,
                    a_feature,
                    a_concrete,
                    a_ambiguous,
                    a_mismatches,
                    a_length,
                ) == (
                    b_library,
                    b_name,
                    b_version,
                    b_entry,
                    b_entry_revision,
                    b_source,
                    b_feature,
                    b_concrete,
                    b_ambiguous,
                    b_mismatches,
                    b_length,
                )
            }
            (
                SourcePin::Dataset {
                    source_blake3: a_source,
                    format: a_format,
                    seqid: a_seqid,
                    proposal_id: _,
                    source_rows: a_rows,
                },
                SourcePin::Dataset {
                    source_blake3: b_source,
                    format: b_format,
                    seqid: b_seqid,
                    proposal_id: _,
                    source_rows: b_rows,
                },
            ) => (a_source, a_format, a_seqid, a_rows) == (b_source, b_format, b_seqid, b_rows),
            (SourcePin::Donor { .. }, SourcePin::Donor { .. }) => self.source == other.source,
            _ => false,
        }
    }
    pub(crate) fn new(source: SourcePin, feature: &Feature) -> Result<Self, Error> {
        let contribution = Self {
            source,
            feature_revision: semantic_revision(feature)?,
        };
        contribution.validate()?;
        Ok(contribution)
    }
    pub(super) fn validate(&self) -> Result<(), Error> {
        if !hex(&self.feature_revision) {
            return Err(bad("Invalid contribution feature pin"));
        }
        match &self.source {
            SourcePin::Library {
                library_revision,
                library_name,
                library_version,
                entry_id,
                entry_revision,
                source_revision,
                source_feature_index,
                proposal_id,
                concrete_matches,
                compatible_ambiguous,
                mismatches,
                template_length,
            } => {
                if !hex(library_revision)
                    || !hex(entry_revision)
                    || !hex(source_revision)
                    || entry_id != &format!("entry-{entry_revision}")
                    || !text(library_name, 1024)
                    || !text(library_version, 128)
                    || *source_feature_index >= 100_000
                    || !text(proposal_id, 128)
                    || *template_length == 0
                    || *template_length > library::MAX_TEMPLATE_BASES
                    || concrete_matches
                        .checked_add(*compatible_ambiguous)
                        .and_then(|n| n.checked_add(*mismatches))
                        != Some(*template_length)
                {
                    return Err(bad("Invalid documentary library contribution"));
                }
            }
            SourcePin::Dataset {
                source_blake3,
                seqid,
                proposal_id,
                source_rows,
                ..
            } => {
                if !hex(source_blake3)
                    || !text(seqid, 8192)
                    || !text(proposal_id, 128)
                    || source_rows.is_empty()
                    || source_rows.len() > 10_000
                    || source_rows.iter().any(|n| *n == 0)
                    || source_rows.windows(2).any(|pair| pair[0] >= pair[1])
                {
                    return Err(bad("Invalid documentary dataset contribution"));
                }
            }
            SourcePin::Donor {
                event,
                feature_index,
                annotation_revision,
                transform,
            } => {
                crate::material::run::event_shape(event)?;
                if *feature_index >= 100_000
                    || !hex(annotation_revision)
                    || transform.origin > collagene_model::MAX_SEQUENCE_LEN
                {
                    return Err(bad("Invalid documentary donor contribution"));
                }
            }
        }
        bounded_json::preflight(self, MAX_LEDGER_BYTES)
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Ledger {
    pub schema_version: u16,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Contribution, 32>")]
    pub contributions: Vec<Contribution>,
}
impl Ledger {
    /// Membership by stable evidence, preserving the first documentary proposal ID.
    pub fn contains_contribution(&self, contribution: &Contribution) -> bool {
        self.contributions
            .iter()
            .any(|old| old.same_evidence(contribution))
    }
    pub fn from_feature(feature: &Feature) -> Result<Self, Error> {
        collagene_model::bounds::feature_preflight(&[std::slice::from_ref(feature)])?;
        let mut values = feature
            .qualifiers
            .iter()
            .filter(|(key, _)| key == QUALIFIER);
        let Some((_, value)) = values.next() else {
            return Ok(Self {
                schema_version: 1,
                contributions: vec![],
            });
        };
        if values.next().is_some() {
            return Err(bad("Duplicate annotation contribution qualifier"));
        }
        let value = value
            .as_ref()
            .ok_or_else(|| bad("Contribution qualifier requires typed content"))?;
        if value.len() > MAX_LEDGER_BYTES {
            return Err(Error::TooLarge);
        }
        let ledger: Self = serde_json::from_str(value)
            .map_err(|_| bad("Invalid annotation contribution ledger"))?;
        ledger.validate()?;
        Ok(ledger)
    }
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != 1 {
            return Err(Error::Unsupported(
                "Unknown annotation contribution schema".into(),
            ));
        }
        if self.contributions.len() > MAX_CONTRIBUTIONS {
            return Err(Error::TooLarge);
        }
        for (i, contribution) in self.contributions.iter().enumerate() {
            contribution.validate()?;
            if self.contributions[..i]
                .iter()
                .any(|old| old.same_evidence(contribution))
            {
                return Err(bad("Duplicate annotation contribution"));
            }
        }
        bounded_json::preflight(self, MAX_LEDGER_BYTES)
    }
    pub(crate) fn add_to(
        feature: &mut Feature,
        contribution: &Contribution,
    ) -> Result<bool, Error> {
        let mut ledger = Self::from_feature(feature)?;
        contribution.validate()?;
        if ledger.contains_contribution(contribution) {
            return Ok(false);
        }
        if ledger.contributions.len() == MAX_CONTRIBUTIONS {
            return Err(Error::TooLarge);
        }
        ledger.contributions.push(contribution.clone());
        ledger.validate()?;
        let value = serde_json::to_string(&ledger).map_err(|_| Error::TooLarge)?;
        if let Some((_, previous)) = feature
            .qualifiers
            .iter_mut()
            .find(|(key, _)| key == QUALIFIER)
        {
            *previous = Some(value);
        } else {
            feature.qualifiers.push((QUALIFIER.into(), Some(value)));
        }
        Ok(true)
    }
}
pub(super) fn semantic_qualifiers(
    feature: &Feature,
) -> impl Iterator<Item = &(String, Option<String>)> {
    feature
        .qualifiers
        .iter()
        .filter(|(key, _)| key != QUALIFIER && key != "annotation_exchange:source_record")
}
pub(super) fn semantic_revision(feature: &Feature) -> Result<String, Error> {
    collagene_model::bounds::feature_preflight(&[std::slice::from_ref(feature)])?;
    bounded_json::hash(&(
        &feature.kind,
        &feature.location,
        semantic_qualifiers(feature).collect::<Vec<_>>(),
    ))
}
pub(super) fn same_meaning(a: &Feature, b: &Feature) -> bool {
    a.kind == b.kind && semantic_qualifiers(a).eq(semantic_qualifiers(b))
}
