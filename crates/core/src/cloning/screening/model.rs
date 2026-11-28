use crate::{
    chemistry::{
        Strand,
        tracked_duplex::{CheckedCleavageEvidence, GenerationAssumptions},
    },
    cloning::restriction::catalog::CatalogSelection,
    material::derivation::SourceOrigin,
    sequence::analysis::Primer,
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};

pub const ENGINE: &str = "named-product-screening-1";
pub const MAX_COMPONENTS: usize = 16;
pub const MAX_SUBJECTS: usize = 16;
pub const MAX_TESTS: usize = 64;
pub const MAX_REQUEST_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_REPORT_BYTES: usize = 8 * 1024 * 1024;
pub const SCOPE: &str = "Ideal modeled size signatures under explicit native exact-PCR or positional checked-cleavage hypotheses. All enumerated exact PCR geometries are possible modeled products, not a claim that they amplify experimentally. The conditional ideal union assumes every enumerated modeled size is available for comparison. Mixtures have no declared concentrations; equal-size products collapse to one observable size. No intensity, gel migration, detection sensitivity, yield, probability, experimental completion or automatic biological winner is inferred. A user-declared size tolerance is a comparison rule, not calibrated gel resolution. Saved-origin replay and screening share the actual exact-PCR and cleavage work ledgers; other native algorithms retain their separately documented caps. Incomplete source evidence withholds screening conclusions.";
macro_rules! data {
    ($item:item) => {
        #[derive(Debug, Clone, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        $item
    };
}
macro_rules! value {
    ($item:item) => {
        #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        $item
    };
}
data! {
    pub struct Component {
        pub id: String,
        pub label: String,
        pub revision: String,
        pub molecule: Molecule,
        pub origin: SourceOrigin,
        /// Only OriginalReference may declare generation. None retains unknown modifications.
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub original_generation: Option<GenerationAssumptions>,
    }
}
value! {
    pub enum DecoyKind {
        EmptyVector,
        WrongOrientation,
        TandemInsert,
        ParentCarryover,
        Other,
    }
}
value! {
    pub enum Role {
        Intended,
        Decoy(DecoyKind),
    }
}
value! {
    pub struct Subject {
        pub id: String,
        pub label: String,
        pub role: Role,
        /// Explicit retained components; neither labels nor roles construct molecules.
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 16>")]
        pub components: Vec<String>,
    }
}
data! {
    pub enum Method {
        Digest { selection: CatalogSelection },
        ExactPcr { forward: Primer, reverse: Primer, template_strand: Strand },
    }
}
data! {
    pub struct Test {
        pub id: String,
        pub label: String,
        pub method: Method,
    }
}
value! {
    pub struct IdealSizeReadout {
        /// Integer units of half a bp, using mean top/bottom strand length for digest fragments.
        /// A witness must have NO opposite-signature size at distance <= this tolerance.
        pub counterpart_tolerance_twice_bp: usize,
    }
}
value! {
    pub struct Limits {
        pub max_cleavage_work: usize,
        pub max_pcr_work: usize,
        pub max_sites_per_cell: usize,
        pub max_products_per_cell: usize,
        pub max_panel_checks: usize,
        pub max_report_bytes: usize,
    }
}
impl Limits {
    pub fn validate(&self) -> Result<(), Error> {
        if !(1..=100_000_000).contains(&self.max_cleavage_work)
            || !(1..=200_000_000).contains(&self.max_pcr_work)
            || !(1..=64).contains(&self.max_sites_per_cell)
            || !(1..=64).contains(&self.max_products_per_cell)
            || !(1..=100_000).contains(&self.max_panel_checks)
            || !(4096..=MAX_REPORT_BYTES).contains(&self.max_report_bytes)
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_cleavage_work: 100_000_000,
            max_pcr_work: 200_000_000,
            max_sites_per_cell: 64,
            max_products_per_cell: 64,
            max_panel_checks: 100_000,
            max_report_bytes: MAX_REPORT_BYTES,
        }
    }
}
data! {
    pub struct Request {
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Component, 16>")]
        pub components: Vec<Component>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Subject, 16>")]
        pub subjects: Vec<Subject>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Test, 64>")]
        pub tests: Vec<Test>,
        pub readout: IdealSizeReadout,
        pub limits: Limits,
    }
}
value! {
    pub enum CellState {
        Complete,
        Ambiguous,
        Limited,
        Unsearched,
        Unsupported,
    }
}
value! {
    pub struct PcrGeometry {
        pub forward_start: usize,
        pub reverse_start: usize,
        pub crosses_origin: bool,
        pub product_length: usize,
    }
}
data! {
    pub struct Cell {
        pub test_index: usize,
        pub component_index: usize,
        pub source_revision: String,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub material_revision: Option<String>,
        pub state: CellState,
        /// Sorted distinct sizes; incomplete cells never certify absence.
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 64>")]
        pub sizes_twice_bp: Vec<usize>,
        pub native_product_count: usize,
        /// Enumeration completeness is independent of uniqueness or experimental amplification.
        pub enumeration_complete: bool,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, PcrGeometry, 64>")]
        pub pcr_geometries: Vec<PcrGeometry>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub cleavage: Option<CheckedCleavageEvidence>,
        pub reason: String,
    }
}
value! {
    pub struct Signature {
        pub subject_index: usize,
        pub state: CellState,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 1024>")]
        pub sizes_twice_bp: Vec<usize>,
    }
}
value! {
    pub enum WitnessSide {
        Intended,
        Decoy,
    }
}
value! {
    pub struct Witness {
        pub side: WitnessSide,
        pub size_twice_bp: usize,
    }
}
value! {
    pub enum Distinction {
        Separated(Witness),
        IdenticalSignature,
        WithinDeclaredTolerance,
        Unresolved,
    }
}
value! {
    pub struct Comparison {
        pub decoy_subject_index: usize,
        pub distinction: Distinction,
    }
}
data! {
    pub struct TestResult {
        pub test_index: usize,
        pub definition: Test,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Signature, 16>")]
        pub signatures: Vec<Signature>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Comparison, 15>")]
        pub comparisons: Vec<Comparison>,
        pub separated_decoys: usize,
        pub complete_comparisons: usize,
    }
}
value! {
    pub struct PanelSearch {
        /// Best among the combinations actually searched, of size one through three.
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 3>")]
        pub test_indices: Vec<usize>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 15>")]
        pub separated_decoy_indices: Vec<usize>,
        pub evaluated_combinations: usize,
        pub comparison_work: usize,
        pub complete_through_size_three: bool,
    }
}
value! {
    pub enum ComponentOrigin {
        OriginalReference,
        SavedNative,
    }
}
data! {
    pub struct ComponentBinding {
        pub id: String,
        pub label: String,
        pub source_revision: String,
        pub origin: ComponentOrigin,
        pub origin_revision: String,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub original_generation: Option<GenerationAssumptions>,
    }
}
data! {
    pub struct Report {
        pub engine_version: String,
        pub context_revision: String,
        pub readout: IdealSizeReadout,
        pub limits: Limits,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, ComponentBinding, 16>")]
        pub components: Vec<ComponentBinding>,
        /// Original order and labels retained, even when identical DNA is nondiscriminable.
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Subject, 16>")]
        pub subjects: Vec<Subject>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Cell, 1024>")]
        pub cells: Vec<Cell>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, TestResult, 64>")]
        pub tests: Vec<TestResult>,
        /// Deterministic documentary distinction count, then complete comparisons, then input order.
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 64>")]
        pub ranked_test_indices: Vec<usize>,
        pub panel: PanelSearch,
        pub cleavage_work: usize,
        pub pcr_work: usize,
        pub statement: String,
        pub report_revision: String,
    }
}
