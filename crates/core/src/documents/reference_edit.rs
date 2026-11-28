//! Canonical typed reference edits, shared by native history, editor and MCP.
//! These operations change reference data; they never grant physical authority.
use crate::sequence::{annotations, editing};
use collagene_model::{Error, Feature, Molecule, Topology};
use serde::{Deserialize, Serialize};

mod ancestry;

pub const MAX_ACTIONS: usize = 1024;
pub const MAX_EDIT_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_APPLY_BYTES: usize = 64 * 1024 * 1024;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum Operation {
    ReverseComplement,
    /// Literal native edit; annotation claims follow the shared model policy.
    ReplaceSequence {
        start: usize,
        end: usize,
        replacement: String,
    },
    /// Annotation meaning changes without changing the sequence calls.
    FlipAnnotationStrands {
        indices: Vec<usize>,
    },
    /// The actual source-bound annotation batch requested by the user.
    EditAnnotations {
        request: annotations::batch::Request,
    },
    MergeAnnotations {
        request: annotations::merge::Request,
    },
    PasteReference {
        source: Molecule,
        source_event: crate::material::lineage::EventBinding,
        selection: crate::sequence::reference_transfer::Selection,
        fragment_revision: String,
        start: usize,
        end: usize,
    },
    RotateOrigin {
        new_origin: usize,
    },
    ExtractRegion {
        start: usize,
        end: usize,
    },
    /// None appends a feature; Some replaces that exact annotation index.
    SetFeature {
        index: Option<usize>,
        feature: Feature,
    },
    DeleteFeature {
        index: usize,
    },
    Rename {
        name: String,
    },
    /// Reference topology only; this never asserts physical ligation or cleavage.
    SetTopology {
        topology: Topology,
    },
}
impl Operation {
    /// Apply the exact ordered actions as one reference transaction. Every
    /// intermediate is admitted before its next copy; callers receive nothing
    /// until the whole transaction succeeds.
    pub fn apply_all(source: &Molecule, actions: &[Self]) -> Result<Molecule, Error> {
        Self::apply_all_native(source, actions, None).map(|(molecule, _)| molecule)
    }
    pub(crate) fn apply_all_with_ancestry(
        source: &Molecule,
        source_event: &crate::material::lineage::EventBinding,
        actions: &[Self],
        budget: crate::verification::ancestry::TraceBudget,
    ) -> Result<
        (
            Molecule,
            Result<crate::verification::ancestry::ReconstructionTrace, Error>,
        ),
        Error,
    > {
        let (molecule, trace) =
            Self::apply_all_native(source, actions, Some((source_event, budget)))?;
        Ok((
            molecule,
            trace.ok_or_else(|| Error::Format("native edit trace was not requested".into()))?,
        ))
    }
    pub(crate) fn actions_revision(
        source: &crate::material::lineage::EventBinding,
        actions: &[Self],
    ) -> Result<String, Error> {
        if actions.is_empty()
            || actions.len() > MAX_ACTIONS
            || source.document_id.len() > 4096
            || source.event_id.len() > 4096
            || !crate::documents::workspace::digest(&source.revision)
        {
            return Err(Error::Format(
                "bounded exact reference edit context required".into(),
            ));
        }
        for action in actions {
            action.preflight()?;
        }
        crate::encoding::bounded_json::preflight(&actions, MAX_EDIT_BYTES)?;
        crate::encoding::bounded_json::hash(&("collagene-reference-edit-v1", source, actions))
    }
    fn apply_all_native(
        source: &Molecule,
        actions: &[Self],
        trace: Option<(
            &crate::material::lineage::EventBinding,
            crate::verification::ancestry::TraceBudget,
        )>,
    ) -> Result<
        (
            Molecule,
            Option<Result<crate::verification::ancestry::ReconstructionTrace, Error>>,
        ),
        Error,
    > {
        if actions.is_empty() || actions.len() > MAX_ACTIONS {
            return Err(Error::Format("reference edit needs 1–1024 actions".into()));
        }
        collagene_model::bounds::native_frame_preflight(std::iter::once(source), MAX_EDIT_BYTES)?;
        for action in actions {
            action.preflight()?;
        }
        crate::encoding::bounded_json::preflight(&actions, MAX_EDIT_BYTES)?;
        let mut recorder =
            trace.map(|(event, budget)| ancestry::Recorder::new(source, event, actions, budget));
        let mut work = 0usize;
        let mut scan_work = 0usize;
        let mut current = source.clone();
        for (index, action) in actions.iter().enumerate() {
            if let Self::MergeAnnotations { request } = action {
                scan_work = scan_work
                    .checked_add(request.limits.max_comparisons)
                    .filter(|n| *n <= annotations::library::MAX_SCAN_COMPARISONS)
                    .ok_or(Error::TooLarge)?;
            }
            let bytes = crate::verification::evidence::evidence_bytes(&current, MAX_EDIT_BYTES)?;
            work = work
                .checked_add(bytes)
                .filter(|n| *n <= MAX_APPLY_BYTES)
                .ok_or(Error::TooLarge)?;
            let capture = recorder
                .as_mut()
                .is_some_and(|recorder| recorder.capture(&current));
            let (product, receipt) = action.clone().execute(&current, capture)?;
            if let Some(recorder) = &mut recorder {
                recorder.record(index, action, &current, &product, receipt);
            }
            current = product;
        }
        Ok((current, recorder.map(ancestry::Recorder::finish)))
    }
    pub fn preflight(&self) -> Result<(), Error> {
        match self {
            Self::ReplaceSequence { replacement, .. }
                if replacement.len() > collagene_model::MAX_SEQUENCE_LEN.saturating_mul(2) =>
            {
                return Err(Error::TooLarge);
            }
            Self::FlipAnnotationStrands { indices } if indices.len() > 4096 => {
                return Err(Error::TooLarge);
            }
            Self::EditAnnotations { request } if request.indices.len() > 4096 => {
                return Err(Error::TooLarge);
            }
            _ => {}
        }
        if let Self::SetFeature { feature, .. } = self {
            collagene_model::bounds::ModelBudget::default()
                .features(std::slice::from_ref(feature))?;
        }
        match self {
            Self::MergeAnnotations { request } => request.preflight()?,
            Self::PasteReference {
                source,
                source_event,
                ..
            } => {
                collagene_model::bounds::native_frame_preflight(
                    std::iter::once(source),
                    MAX_EDIT_BYTES,
                )?;
                source.validate()?;
                if source.revision() != source_event.revision {
                    return Err(Error::Format(
                        "copied reference event and snapshot disagree".into(),
                    ));
                }
            }
            _ => {}
        }
        crate::encoding::bounded_json::preflight(self, MAX_EDIT_BYTES)?;
        Ok(())
    }
    /// Actual cross-document event claim retained by the workspace/history owner.
    pub fn source_event(&self) -> Option<&crate::material::lineage::EventBinding> {
        match self {
            Self::PasteReference { source_event, .. } => Some(source_event),
            Self::MergeAnnotations { request } => match &request.source {
                annotations::merge::Source::Donor { event, .. } => Some(event),
                _ => None,
            },
            _ => None,
        }
    }
    pub fn name(&self) -> &'static str {
        match self {
            Self::ReverseComplement => "ReverseComplement",
            Self::ReplaceSequence { .. } => "ReplaceSequence",
            Self::FlipAnnotationStrands { .. } => "FlipAnnotationStrands",
            Self::EditAnnotations { .. } => "EditAnnotations",
            Self::MergeAnnotations { .. } => "MergeAnnotations",
            Self::PasteReference { .. } => "PasteReference",
            Self::RotateOrigin { .. } => "RotateOrigin",
            Self::ExtractRegion { .. } => "ExtractRegion",
            Self::SetFeature { .. } => "SetFeature",
            Self::DeleteFeature { .. } => "DeleteFeature",
            Self::Rename { .. } => "Rename",
            Self::SetTopology { .. } => "SetTopology",
        }
    }
    pub fn apply(self, source: &Molecule) -> Result<Molecule, Error> {
        self.execute(source, false).map(|(molecule, _)| molecule)
    }
    fn execute(
        self,
        source: &Molecule,
        capture: bool,
    ) -> Result<(Molecule, Option<Receipt>), Error> {
        self.preflight()?;
        collagene_model::bounds::native_frame_preflight(std::iter::once(source), MAX_EDIT_BYTES)?;
        source.validate()?;
        macro_rules! observed {
            ($native:expr,$ordinary:expr,$kind:path) => {
                if capture {
                    $native.map(|(molecule, receipt)| (molecule, Some($kind(receipt))))
                } else {
                    $ordinary.map(|molecule| (molecule, None))
                }
            };
        }
        match self {
            Self::ReverseComplement => observed!(
                annotations::reverse_complement_with_projection(source),
                annotations::reverse_complement(source),
                Receipt::Projection
            ),
            Self::ReplaceSequence {
                start,
                end,
                replacement,
            } => observed!(
                source.replace_with_projection(start, end, &replacement),
                source.replace(start, end, &replacement),
                Receipt::Replacement
            ),
            Self::FlipAnnotationStrands { indices } => {
                let request = annotations::batch::Request {
                    source_revision: source.revision(),
                    indices,
                    kind: None,
                    flip_strand: true,
                    color: annotations::batch::ColorEdit::Keep,
                };
                observed!(
                    annotations::batch::apply_with_projection(source, &request),
                    annotations::batch::apply(source, &request),
                    Receipt::Annotations
                )
            }
            Self::EditAnnotations { request } => observed!(
                annotations::batch::apply_with_projection(source, &request),
                annotations::batch::apply(source, &request),
                Receipt::Annotations
            ),
            Self::MergeAnnotations { request } => {
                if capture {
                    annotations::merge::apply_with_projection(source, &request).map(
                        |(applied, identity)| {
                            (applied.molecule, Some(Receipt::Annotations(identity)))
                        },
                    )
                } else {
                    annotations::merge::apply(source, &request)
                        .map(|applied| (applied.molecule, None))
                }
            }
            Self::PasteReference {
                source: donor,
                source_event: _,
                selection,
                fragment_revision,
                start,
                end,
            } => {
                let request = crate::sequence::reference_transfer::RequestPaste {
                    source: donor,
                    destination: source.clone(),
                    selection,
                    fragment_revision,
                    start,
                    end,
                };
                if capture {
                    let (report, destination, donor) =
                        crate::sequence::reference_transfer::paste_with_projections(&request)?;
                    Ok((report.product, Some(Receipt::Paste { destination, donor })))
                } else {
                    crate::sequence::reference_transfer::paste(&request)
                        .map(|report| (report.product, None))
                }
            }
            Self::RotateOrigin { new_origin } => observed!(
                annotations::rotate_origin_with_projection(source, new_origin),
                annotations::rotate_origin(source, new_origin),
                Receipt::Projection
            ),
            Self::ExtractRegion { start, end } => observed!(
                annotations::extract_region_with_projection(source, start, end),
                annotations::extract_region(source, start, end),
                Receipt::Projection
            ),
            Self::SetFeature { index, feature } => observed!(
                editing::set_feature_with_projection(source, index, feature),
                editing::set_feature(source, index, feature),
                Receipt::Annotations
            ),
            Self::DeleteFeature { index } => observed!(
                editing::delete_feature_with_projection(source, index),
                editing::delete_feature(source, index),
                Receipt::Annotations
            ),
            Self::Rename { name } => editing::rename(source, &name)
                .map(|molecule| (molecule, capture.then_some(Receipt::Retained))),
            Self::SetTopology { topology } => observed!(
                annotations::change_reference_topology_with_projection(source, topology),
                editing::set_reference_topology(source, topology),
                Receipt::Projection
            ),
        }
    }
}
/// Fresh owner receipts only; no deserialization or coordinate-recovery path.
enum Receipt {
    Retained,
    Projection(annotations::ReferenceProjectionReceipt),
    Replacement(collagene_model::annotations::ReplacementProjection),
    Annotations(editing::AnnotationIdentity),
    Paste {
        destination: collagene_model::annotations::ReplacementProjection,
        donor: annotations::ReferenceProjectionReceipt,
    },
}

#[cfg(test)]
mod tests;
