use super::*;
use crate::{formats::annotations::AnnotationImportReport, material::lineage::EventBinding};
use serde::{Deserialize, Serialize};
macro_rules! dto {
    ($item:item) => {
        #[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        $item
    };
}
dto! {
    pub struct Transform {
        /// Rotate the optionally reverse-complemented donor by this origin.
        pub origin: usize,
        pub reverse: bool,
    }
}
dto! {
    pub enum Source {
        Library {
            library: library::FeatureLibrary,
            report: library::ScanReport,
        },
        Dataset {
            report: AnnotationImportReport,
        },
        /// Event ownership must be resolved by the surrounding Workspace owner.
        /// A supplied molecule/event tuple alone is only a documentary claim.
        Donor {
            molecule: Molecule,
            event: EventBinding,
            transform: Transform,
            #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 1000>")]
            indices: Vec<usize>,
        },
    }
}
dto! {
    pub enum Decision {
        Append,
        Skip,
        AddProvenance {
            target_index: usize,
            annotation_revision: String,
        },
    }
}
dto! {
    pub struct Selection {
        pub proposal_id: String,
        pub decision: Decision,
    }
}
dto! {
    pub struct Limits {
        pub max_comparisons: usize,
        pub max_relations: usize,
        pub max_output_bytes: usize,
    }
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_comparisons: 1_000_000,
            max_relations: 10_000,
            max_output_bytes: MAX_REPORT_BYTES,
        }
    }
}
dto! {
    pub struct Request {
        pub schema_version: u16,
        pub target_revision: String,
        pub source: Source,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Selection, 1000>")]
        pub selections: Vec<Selection>,
        pub limits: Limits,
    }
}
dto! {
    pub enum RelationKind {
        ExactDuplicate,
        EquivalentAnnotation,
        SameGeometryConflict,
        OppositeOrientation,
        Overlap,
        Unresolved,
    }
}
dto! {
    pub struct Relation {
        pub target_index: usize,
        pub annotation_revision: String,
        pub kind: RelationKind,
    }
}
dto! {
    pub struct Candidate {
        pub proposal_id: String,
        pub feature: Feature,
        pub contribution: Contribution,
        pub decision: Decision,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Relation, 100000>")]
        pub relations: Vec<Relation>,
        pub relations_complete: bool,
    }
}
dto! {
    pub struct Preview {
        pub engine_version: String,
        pub context_revision: String,
        pub target_revision: String,
        pub source_complete: bool,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Candidate, 1000>")]
        pub candidates: Vec<Candidate>,
        pub comparisons: usize,
        pub complete: bool,
        pub report_revision: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 256>")]
        pub warnings: Vec<String>,
    }
}
dto! {
    pub struct ActionReceipt {
        pub proposal_id: String,
        pub decision: Decision,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub destination_index: Option<usize>,
        pub contribution_added: bool,
    }
}
dto! {
    pub struct MergeReceipt {
        pub context_revision: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, ActionReceipt, 1000>")]
        pub actions: Vec<ActionReceipt>,
        pub retained_original_features: usize,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 1000>")]
        pub added_feature_indices: Vec<usize>,
        pub product_revision: String,
        /// Actual source/review/application comparison work, one request ledger.
        pub comparisons: usize,
    }
}
dto! {
    pub struct Applied {
        pub molecule: Molecule,
        pub review: Preview,
        pub receipt: MergeReceipt,
    }
}
