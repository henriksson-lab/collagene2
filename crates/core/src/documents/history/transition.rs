//! Typed native event provenance. Descriptive labels never supply replay authority.
use crate::{
    documents::{reference_edit, workspace},
    material::lineage::EventBinding,
};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnalysisProducer {
    MaterialDerivation,
    RestrictionInspection,
    ReadReview,
    ReadContig,
    Translation,
    CodonDesign,
    CrisprDesign,
    MutationIntent,
    DomesticationIntent,
    AssemblyIntent,
    SelectedAllele,
    RnaChemistry,
    AnnotationCuration,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum AnalysisRecord {
    Attachment {
        hash: String,
    },
    /// The named producer must resolve and freshly validate this native context.
    /// A context identifier alone is not evidence that its record is available.
    Context {
        revision: String,
    },
    /// A real computational preview for which no durable recipe was retained.
    /// This is not a causal-edit receipt and cannot stand in for one.
    UnretainedPreview {
        revision: String,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnalysisOrigin {
    pub producer: AnalysisProducer,
    pub record: AnalysisRecord,
}
impl AnalysisOrigin {
    pub fn context(producer: AnalysisProducer, revision: impl Into<String>) -> Self {
        Self {
            producer,
            record: AnalysisRecord::Context {
                revision: revision.into(),
            },
        }
    }
    pub fn preview(producer: AnalysisProducer, revision: impl Into<String>) -> Self {
        Self {
            producer,
            record: AnalysisRecord::UnretainedPreview {
                revision: revision.into(),
            },
        }
    }
    pub fn validate(&self) -> Result<(), Error> {
        let digest = match &self.record {
            AnalysisRecord::Attachment { hash } => hash,
            AnalysisRecord::Context { revision }
            | AnalysisRecord::UnretainedPreview { revision } => revision,
        };
        if !workspace::digest(digest) {
            return Err(Error::Format(
                "invalid native analysis record locator".into(),
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum Transition {
    /// Opened or created as an independent reference, without an edit claim.
    IndependentReference,
    ReferenceEdit {
        source: EventBinding,
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, reference_edit::Operation, 1024>"
        )]
        actions: Vec<reference_edit::Operation>,
        /// Explicit ancillary ancestry, not inferred residue correspondence.
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 127>")]
        supplemental_sources: Vec<String>,
    },
    Fork {
        source: EventBinding,
    },
    CloningStep {
        run_id: String,
        step_index: usize,
    },
    Analysis {
        origin: AnalysisOrigin,
    },
}
impl Transition {
    pub(crate) fn check_fields(
        &self,
        document: &str,
        revision: &str,
        parents: &[String],
    ) -> Result<(), Error> {
        let source_matches = |source: &EventBinding| {
            workspace::identifier(&source.document_id)
                && workspace::digest(&source.revision)
                && parents.first() == Some(&source.event_id)
        };
        let valid = match self {
            Self::IndependentReference => parents.is_empty(),
            Self::ReferenceEdit {
                source,
                actions,
                supplemental_sources,
            } => {
                if actions.is_empty()
                    || actions.len() > reference_edit::MAX_ACTIONS
                    || supplemental_sources.len() > 127
                {
                    return Err(Error::TooLarge);
                }
                for action in actions {
                    action.preflight()?;
                }
                for action in actions {
                    if let Some(source_event) = action.source_event() {
                        if source_event != source
                            && !supplemental_sources.contains(&source_event.event_id)
                        {
                            return Err(Error::Format(
                                "reference edit source is not retained among event parents".into(),
                            ));
                        }
                    }
                }
                crate::encoding::bounded_json::preflight(&actions, reference_edit::MAX_EDIT_BYTES)?;
                source_matches(source)
                    && source.document_id == document
                    && parents.len() == 1 + supplemental_sources.len()
                    && parents.get(1..) == Some(supplemental_sources.as_slice())
            }
            Self::Fork { source } => {
                source_matches(source)
                    && source.document_id != document
                    && source.revision == revision
                    && parents.len() == 1
            }
            Self::CloningStep { run_id, step_index } => {
                workspace::digest(run_id) && *step_index < 1000 && !parents.is_empty()
            }
            Self::Analysis { origin } => {
                origin.validate()?;
                !parents.is_empty()
            }
        };
        if !valid {
            return Err(Error::Format(
                "native event transition and source bindings disagree".into(),
            ));
        }
        Ok(())
    }
    pub(crate) fn validate_sources(
        &self,
        events: &BTreeMap<&str, &workspace::RevisionEvent>,
        attachments: &BTreeMap<String, workspace::WorkspaceAttachment>,
    ) -> Result<(), Error> {
        if let Self::ReferenceEdit { source, .. } | Self::Fork { source } = self {
            let event = events
                .get(source.event_id.as_str())
                .ok_or_else(|| Error::Format("typed reference source event is absent".into()))?;
            if event.document_id != source.document_id || event.revision != source.revision {
                return Err(Error::Format(
                    "typed reference source ownership changed".into(),
                ));
            }
        }
        if let Self::ReferenceEdit { actions, .. } = self {
            for action in actions {
                if let Some(source_event) = action.source_event() {
                    let event = events.get(source_event.event_id.as_str()).ok_or_else(|| {
                        Error::Format("reference edit source event is absent".into())
                    })?;
                    if event.document_id != source_event.document_id
                        || event.revision != source_event.revision
                    {
                        return Err(Error::Format(
                            "reference edit source event binding changed".into(),
                        ));
                    }
                }
            }
        }
        if let Self::Analysis {
            origin:
                AnalysisOrigin {
                    record: AnalysisRecord::Attachment { hash },
                    ..
                },
        } = self
        {
            if !attachments.contains_key(hash) {
                return Err(Error::Format("native analysis attachment is absent".into()));
            }
        }
        Ok(())
    }
}
