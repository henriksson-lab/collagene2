//! Exact user roles and event/settings guards; scientific mapping stays native.
use collagene_core::{
    alignment::{ComparisonOptions, design_consensus as native},
    construct_qc,
    documents::workspace::Workspace,
    material::lineage::{self, EventBinding},
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct CodingDraft {
    pub enabled: bool,
    pub feature: String,
    pub genetic_code: String,
    pub rationale: String,
}
impl Default for CodingDraft {
    fn default() -> Self {
        Self {
            enabled: false,
            feature: String::new(),
            genetic_code: "1".into(),
            rationale: String::new(),
        }
    }
}
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct Draft {
    pub single_circular_turn: bool,
    pub axis_rationale: String,
    pub coding: CodingDraft,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Context {
    pub design: EventBinding,
    pub consensus: EventBinding,
    pub epoch: u64,
    pub comparison: ComparisonOptions,
    pub draft: Draft,
}
#[derive(Clone)]
pub(super) struct Capture {
    pub context: Context,
    pub request: native::Request,
    pub revision: String,
    pub fingerprint: String,
    pub design_name: String,
    pub consensus_name: String,
}
impl Capture {
    pub fn prepare(workspace: &Workspace, context: Context) -> Result<Self, String> {
        let consensus_axis = if context.draft.single_circular_turn {
            native::ConsensusAxis::DeclaredSingleCircularTurn {
                rationale: context.draft.axis_rationale.clone(),
            }
        } else {
            native::ConsensusAxis::StoredTopology
        };
        let mut request = native::Request {
            schema_version: 1,
            design: context.design.clone(),
            consensus: context.consensus.clone(),
            comparison: context.comparison,
            consensus_axis,
            cds: Vec::new(),
            limits: native::Limits::default(),
        };
        let (design, consensus) = request
            .sources(workspace)
            .map_err(|error| error.to_string())?;
        if context.draft.coding.enabled {
            let draft = &context.draft.coding;
            let index = draft
                .feature
                .parse::<usize>()
                .map_err(|_| "Choose an explicit design CDS annotation.")?;
            let feature = design
                .features
                .get(index)
                .filter(|feature| feature.kind.eq_ignore_ascii_case("CDS"))
                .ok_or("The selected design annotation is not a CDS.")?;
            if draft.rationale.trim().is_empty() || draft.rationale.len() > 1024 {
                return Err("Explain the declared design CDS role (1–1024 bytes).".into());
            }
            request.cds.push(native::CdsSelection {
                id: "selected-design-cds".into(),
                design_cds: construct_qc::AnnotationAssumption {
                    feature_index: index,
                    source_revision: context.design.revision.clone(),
                    annotation_revision: construct_qc::annotation_revision(feature)
                        .map_err(|error| error.to_string())?,
                    rationale: draft.rationale.clone(),
                },
                genetic_code: draft
                    .genetic_code
                    .parse::<u8>()
                    .map_err(|_| "Choose an explicit genetic code.")?,
            });
        }
        request
            .validate(workspace)
            .map_err(|error| error.to_string())?;
        let revision = native::revision(workspace, &request).map_err(|error| error.to_string())?;
        let fingerprint =
            lineage::workspace_fingerprint(workspace).map_err(|error| error.to_string())?;
        Ok(Self {
            context,
            request,
            revision,
            fingerprint,
            design_name: design.name.clone(),
            consensus_name: consensus.name.clone(),
        })
    }
    pub fn guard(&self, workspace: &Workspace, current: &Context) -> Result<(), String> {
        if self.context != *current
            || lineage::workspace_fingerprint(workspace).map_err(|error| error.to_string())?
                != self.fingerprint
        {
            return Err(
                "The design, supplied consensus, settings or workspace changed. Compare again."
                    .into(),
            );
        }
        Ok(())
    }
    pub fn validate_report(
        &self,
        workspace: &Workspace,
        report: &native::Report,
    ) -> Result<(), String> {
        if report.evidence.context_revision != self.revision {
            return Err("The semantic report belongs to another comparison.".into());
        }
        report
            .validate_request(workspace, &self.request)
            .map_err(|error| error.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_model::{Molecule, Topology};

    fn fixture() -> (Workspace, Context) {
        let sequence = Molecule::dna("Design", "ATGAAATAA", Topology::Circular).unwrap();
        let workspace =
            Workspace::from_documents("Roles", vec![sequence.clone(), sequence], 0).unwrap();
        let binding = |index: usize| {
            let document = &workspace.documents[index];
            EventBinding {
                document_id: document.id.clone(),
                event_id: document.head_event.clone().unwrap(),
                revision: document.head.clone(),
            }
        };
        let context = Context {
            design: binding(0),
            consensus: binding(1),
            epoch: 3,
            comparison: ComparisonOptions::default(),
            draft: Draft::default(),
        };
        (workspace, context)
    }
    #[test]
    fn identical_calls_do_not_merge_roles_or_retarget_captured_events() {
        let (workspace, context) = fixture();
        assert_eq!(context.design.revision, context.consensus.revision);
        assert_ne!(context.design, context.consensus);
        let capture = Capture::prepare(&workspace, context.clone()).unwrap();
        capture.guard(&workspace, &context).unwrap();
        let mut other = context.clone();
        std::mem::swap(&mut other.design, &mut other.consensus);
        assert!(capture.guard(&workspace, &other).is_err());
        let mut other = context.clone();
        other.epoch += 1;
        assert!(capture.guard(&workspace, &other).is_err());
        let mut other = context.clone();
        other.comparison.both_query_strands = !other.comparison.both_query_strands;
        assert!(capture.guard(&workspace, &other).is_err());
        let mut other = context.clone();
        other.draft.single_circular_turn = true;
        other.draft.axis_rationale = "Declared sequenced turn".into();
        assert!(capture.guard(&workspace, &other).is_err());
        let mut other = context;
        other.consensus = other.design.clone();
        assert!(Capture::prepare(&workspace, other).is_err());
    }
    #[test]
    fn annotation_role_is_explicit_and_workspace_changes_invalidate_export() {
        let (workspace, context) = fixture();
        let capture = Capture::prepare(&workspace, context.clone()).unwrap();
        let mut changed = workspace.clone();
        changed.name = "Different retained workspace".into();
        assert!(capture.guard(&changed, &context).is_err());
        let mut coding = context;
        coding.draft.coding.enabled = true;
        coding.draft.coding.feature = "0".into();
        coding.draft.coding.rationale = "Selected intended CDS".into();
        assert!(Capture::prepare(&workspace, coding).is_err());
    }
}
