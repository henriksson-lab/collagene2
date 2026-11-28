//! One source/event-bound semantic record; the actual aligned rows live once.
use crate::{
    alignment::{ComparisonOptions, review},
    construct_qc::{AnnotationAssumption, CodePin, Verdict},
    material::lineage::EventBinding,
    sequence::translation::TranslationReport,
};
use collagene_model::{Error, Topology};
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
    pub enum ConsensusAxis {
        StoredTopology,
        DeclaredSingleCircularTurn { rationale: String },
    }
}
dto! {
    pub struct CdsSelection {
        pub id: String,
        pub design_cds: AnnotationAssumption,
        pub genetic_code: u8,
    }
}
dto! {
    pub struct Limits {
        pub max_coding_calls: usize,
        pub max_mapping_work: usize,
        pub max_segments: usize,
        pub max_consequences: usize,
        pub max_report_bytes: usize,
    }
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_coding_calls: 1_000_000,
            max_mapping_work: 8_000_000,
            max_segments: 512,
            max_consequences: 4096,
            max_report_bytes: super::MAX_REPORT_BYTES,
        }
    }
}
impl Limits {
    pub fn validate(&self) -> Result<(), Error> {
        if !(1..=2_000_000).contains(&self.max_coding_calls)
            || !(1..=64_000_000).contains(&self.max_mapping_work)
            || !(1..=512).contains(&self.max_segments)
            || !(1..=4096).contains(&self.max_consequences)
            || !(4096..=super::MAX_REPORT_BYTES).contains(&self.max_report_bytes)
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
dto! {
    pub struct Request {
        pub schema_version: u16,
        pub design: EventBinding,
        pub consensus: EventBinding,
        pub comparison: ComparisonOptions,
        pub consensus_axis: ConsensusAxis,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, CdsSelection, 32>")]
        pub cds: Vec<CdsSelection>,
        pub limits: Limits,
    }
}
dto! {
    pub enum SourceRole {
        Design,
        Consensus,
    }
}
dto! {
    pub struct ProvenanceField {
        pub key: String,
        pub value: String,
        /// False when this display prefix omits bytes or control characters.
        pub value_complete: bool,
    }
}
dto! {
    pub enum OriginalState {
        Unavailable,
        IntegrityCheckedMatchesCalls,
        IntegrityCheckedDifferentCalls,
        InvalidOrUnsupported,
    }
}
dto! {
    pub struct OriginalEvidence {
        pub state: OriginalState,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub format: Option<String>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub blake3: Option<String>,
        pub byte_count: usize,
        pub message: String,
    }
}
dto! {
    pub struct SourceProvenance {
        pub role: SourceRole,
        pub event: EventBinding,
        pub topology: Topology,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, ProvenanceField, 8>")]
        pub fields: Vec<ProvenanceField>,
        pub original: OriginalEvidence,
    }
}
dto! {
    pub struct CallSpan {
        pub start: usize,
        pub end: usize,
        pub reverse: bool,
    }
}
dto! {
    pub enum Rule {
        AnnotationBasis,
        CoordinateScope,
        Register,
        CodonConsequence,
        InternalStop,
        TerminalContext,
        Scope,
    }
}
dto! {
    pub struct Finding {
        pub rule: Rule,
        pub verdict: Verdict,
        pub message: String,
    }
}
dto! {
    pub struct RegisterSegment {
        /// Half-open coding ordinals after the design codon_start prefix.
        pub design_start: usize,
        pub design_end: usize,
        /// Half-open projected coding ordinals, not stored genomic coordinates.
        pub projected_start: usize,
        pub projected_end: usize,
        pub offset_mod3: u8,
    }
}
dto! {
    pub enum ConsequenceKind {
        Synonymous,
        Missense,
        DefiniteInternalStop,
        StopLost,
        InsertedCodon,
        Unresolved,
    }
}
dto! {
    pub struct CodonConsequence {
        pub projected_codon_index: usize,
        pub projected_codon: String,
        pub projected_amino_acid: char,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub design_codon_index: Option<usize>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub design_codon: Option<String>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub design_amino_acid: Option<char>,
        pub kind: ConsequenceKind,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, CallSpan, 3>")]
        pub design_coordinates: Vec<CallSpan>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, CallSpan, 3>")]
        pub consensus_coordinates: Vec<CallSpan>,
        pub message: String,
    }
}
dto! {
    pub struct CdsReport {
        pub id: String,
        pub design_cds: AnnotationAssumption,
        pub code: CodePin,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub design_translation: Option<TranslationReport>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub projected_translation: Option<TranslationReport>,
        pub mapping_complete: bool,
        pub register_complete: bool,
        pub consequences_complete: bool,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, RegisterSegment, 512>")]
        pub register_segments: Vec<RegisterSegment>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, CodonConsequence, 4096>")]
        pub consequences: Vec<CodonConsequence>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Finding, 16>")]
        pub findings: Vec<Finding>,
    }
}
dto! {
    pub struct WorkSummary {
        pub mapping_work: usize,
        pub coding_calls: usize,
        pub segments: usize,
        pub consequences: usize,
    }
}
dto! {
    pub struct Evidence {
        pub engine_version: String,
        pub context_revision: String,
        pub workspace_fingerprint: String,
        pub request_revision: String,
        pub pairwise_context_revision: String,
        pub pairwise_report_revision: String,
        pub request: Request,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, SourceProvenance, 2>")]
        pub provenance: Vec<SourceProvenance>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, CdsReport, 32>")]
        pub cds: Vec<CdsReport>,
        pub work: WorkSummary,
        pub coding_complete: bool,
        pub selected_path_conditional: bool,
        pub statement: String,
        pub evidence_revision: String,
    }
}
dto! {
    pub struct Report {
        pub pairwise: review::Report,
        pub evidence: Evidence,
    }
}
