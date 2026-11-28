//! Immutable reference edits use the same validated operations as the editor.
use crate::session::references::DocumentRef;
use collagene_core::sequence::{annotations, editing};
use collagene_model::{Error, Feature, Molecule, Topology};
use schemars::JsonSchema;
use serde::Deserialize;

#[derive(Debug, Clone, Deserialize, JsonSchema)]
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
    pub(super) fn preflight(&self) -> Result<(), Error> {
        match self {
            Self::ReplaceSequence { replacement, .. }
                if replacement.len() > collagene_model::MAX_SEQUENCE_LEN.saturating_mul(2) =>
            {
                return Err(Error::TooLarge);
            }
            Self::FlipAnnotationStrands { indices } if indices.len() > 4096 => {
                return Err(Error::TooLarge);
            }
            _ => {}
        }
        if let Self::SetFeature { feature, .. } = self {
            collagene_model::bounds::ModelBudget::default()
                .features(std::slice::from_ref(feature))?;
        }
        Ok(())
    }
    pub(super) fn name(&self) -> &'static str {
        match self {
            Self::ReverseComplement => "ReverseComplement",
            Self::ReplaceSequence { .. } => "ReplaceSequence",
            Self::FlipAnnotationStrands { .. } => "FlipAnnotationStrands",
            Self::RotateOrigin { .. } => "RotateOrigin",
            Self::ExtractRegion { .. } => "ExtractRegion",
            Self::SetFeature { .. } => "SetFeature",
            Self::DeleteFeature { .. } => "DeleteFeature",
            Self::Rename { .. } => "Rename",
            Self::SetTopology { .. } => "SetTopology",
        }
    }
    pub(super) fn apply(self, source: &Molecule) -> Result<Molecule, Error> {
        match self {
            Self::ReverseComplement => annotations::reverse_complement(source),
            Self::ReplaceSequence {
                start,
                end,
                replacement,
            } => source.replace(start, end, &replacement),
            Self::FlipAnnotationStrands { indices } => annotations::batch::apply(
                source,
                &annotations::batch::Request {
                    source_revision: source.revision(),
                    indices,
                    kind: None,
                    flip_strand: true,
                    color: annotations::batch::ColorEdit::Keep,
                },
            ),
            Self::RotateOrigin { new_origin } => annotations::rotate_origin(source, new_origin),
            Self::ExtractRegion { start, end } => annotations::extract_region(source, start, end),
            Self::SetFeature { index, feature } => editing::set_feature(source, index, feature),
            Self::DeleteFeature { index } => editing::delete_feature(source, index),
            Self::Rename { name } => editing::rename(source, &name),
            Self::SetTopology { topology } => editing::set_reference_topology(source, topology),
        }
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub document: DocumentRef,
    pub operation: Operation,
}
