//! Current documentary contract for native-history CDS comparisons.
use crate::{
    construct_qc::{AnnotationAssumption, Citation, CodePin, Verdict},
    material::lineage::EventBinding,
};
use collagene_model::Location;
use serde::{Deserialize, Serialize};
macro_rules! dto {
    ($($item:item)*) => {$(
        #[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        $item
    )*};
}

dto! {
    pub struct Limits {
        pub max_history_events: usize,
        pub max_coding_calls: usize,
        pub max_work: usize,
        pub max_segments: usize,
        pub max_stops: usize,
    }
    pub struct CdsPair {
        pub id: String,
        pub source: EventBinding,
        pub product: EventBinding,
        pub source_cds: AnnotationAssumption,
        pub product_cds: AnnotationAssumption,
        pub genetic_code: u8,
    }
    pub struct Request {
        pub schema_version: u16,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, CdsPair, 32>")]
        pub pairs: Vec<CdsPair>,
        pub limits: Limits,
    }
    pub enum Rule {
        NativeHistory,
        AnnotationBasis,
        AnnotationCorrespondence,
        CodingCorrespondence,
        JunctionRegister,
        NewInternalStop,
        Scope,
    }
    pub struct Finding {
        pub rule: Rule,
        pub verdict: Verdict,
        pub reason: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Location, 16>")]
        pub product_coordinates: Vec<Location>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Citation, 8>")]
        pub citations: Vec<Citation>,
        pub under_declared_annotation_assumptions: bool,
    }
    pub enum TransitionPin {
        ReferenceEdit { event: EventBinding, source: EventBinding, actions_revision: String },
        Fork { event: EventBinding, source: EventBinding },
        CloningStep { event: EventBinding, run_id: String, step_index: usize, plan_revision: String },
        MaterialDerivation { event: EventBinding, context_revision: String, record_hash: String },
        TopoGraph { record_hash: String, state_revision: String, context_revision: String },
    }
    /// Half-open biological coding ordinals, after each selected codon_start.
    /// First/last parent ordinals are anchors, not an assertion of contiguous retention.
    pub struct RegisterSegment {
        pub product_coding_start: usize,
        pub product_coding_end: usize,
        pub first_parent_coding_ordinal: usize,
        pub last_parent_coding_ordinal: usize,
        pub offset_mod3: u8,
        pub native_junction_proven: bool,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, TransitionPin, 64>")]
        pub causes: Vec<TransitionPin>,
    }
    /// A definite product internal stop with an explicit comparison verdict.
    /// Indeterminate entries never claim that this stop was newly introduced.
    pub struct StopEvidence {
        pub verdict: Verdict,
        pub product_codon_index: usize,
        pub product_coordinates: Location,
        pub product_codon: String,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub parent_codon_index: Option<usize>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub parent_coordinates: Option<Location>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub parent_codon: Option<String>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub parent_amino_acid: Option<char>,
        pub reason: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, TransitionPin, 64>")]
        pub causes: Vec<TransitionPin>,
    }
    pub struct AnnotationCorrespondence {
        pub verdict: Verdict,
        pub reason: String,
    }
    pub struct PairReport {
        pub id: String,
        pub source: EventBinding,
        pub product: EventBinding,
        pub source_cds: AnnotationAssumption,
        pub product_cds: AnnotationAssumption,
        pub code: CodePin,
        pub annotation_correspondence: AnnotationCorrespondence,
        pub mapping_complete: bool,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Finding, 32>")]
        pub findings: Vec<Finding>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, RegisterSegment, 512>")]
        pub register_segments: Vec<RegisterSegment>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, StopEvidence, 256>")]
        pub new_stops: Vec<StopEvidence>,
    }
    pub struct Report {
        pub engine_version: String,
        pub context_revision: String,
        pub workspace_fingerprint: String,
        pub request_revision: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, PairReport, 32>")]
        pub pairs: Vec<PairReport>,
        pub work: usize,
        pub coding_calls: usize,
        pub history_events: usize,
        pub statement: String,
        pub report_revision: String,
    }
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_history_events: 256,
            max_coding_calls: 1_000_000,
            max_work: 20_000_000,
            max_segments: 512,
            max_stops: 256,
        }
    }
}
impl Default for Request {
    fn default() -> Self {
        Self {
            schema_version: 1,
            pairs: vec![],
            limits: Limits::default(),
        }
    }
}
