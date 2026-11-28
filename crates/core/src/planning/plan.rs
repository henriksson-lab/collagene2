//! The single typed cloning plan and its owned dependency/identity behavior.
use crate::planning::editing::{MAX_INPUTS, MAX_SOURCES_PER_STEP, MAX_STEPS};
use crate::{
    chemistry::{self, EndTreatment, tracked_duplex::GenerationAssumptions},
    sequence::analysis::Primer,
};
use collagene_model::Molecule;
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PlanInput {
    pub id: String,
    pub revision: String,
    pub molecule: Molecule,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Source {
    Input(String),
    Step(usize),
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct OrientedSource {
    pub source: Source,
    pub reverse: bool,
}
/// Explicit ordered role binding for a fresh TOPO preparation hypothesis.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct TopoSource {
    pub role: String,
    pub source: Source,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum Operation {
    Pcr {
        template: Source,
        forward: Primer,
        reverse: Primer,
        template_strand: chemistry::Strand,
    },
    /// An explicit ordered candidate, not an exhaustive assembly search.
    Overlap {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, MAX_SOURCES_PER_STEP>"
        )]
        fragments: Vec<Source>,
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, MAX_SOURCES_PER_STEP>"
        )]
        overlaps: Vec<usize>,
        circular: bool,
    },
    /// Complete digestion; select a zero-based fragment in ascending top-cut order.
    Digest {
        template: Source,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 128>")]
        enzymes: Vec<String>,
        fragment_index: usize,
    },
    TreatEnds {
        template: Source,
        treatment: EndTreatment,
    },
    /// Declared intended ligation and bounded one-copy alternatives, with explicit orientation and nick acceptance.
    Ligate {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, OrientedSource, MAX_SOURCES_PER_STEP>"
        )]
        fragments: Vec<OrientedSource>,
        circular: bool,
        allow_nicked: bool,
        options: crate::cloning::ligation::outcomes::Options,
    },
    PcrWithMismatches {
        template: Source,
        forward: Primer,
        reverse: Primer,
        options: crate::cloning::pcr::mutagenesis::PcrOptions,
    },
    GoldenGate {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, MAX_SOURCES_PER_STEP>"
        )]
        templates: Vec<Source>,
        options: crate::cloning::golden_gate::GoldenGateOptions,
        candidate_index: usize,
    },
    RepairEnds {
        template: Source,
        method: chemistry::EndRepairMethod,
    },
    /// Two explicitly declared linear 5′→3′ ssDNA oligo inputs. Terminal chemistry
    /// is an explicit plan assumption, never inferred from imported metadata.
    AnnealOligos {
        top: Source,
        bottom: Source,
        options: chemistry::AnnealOptions,
        candidate_index: usize,
    },
    Tail3Prime {
        template: Source,
        base: crate::cloning::ligation::tailing::TailingBase,
    },
    ProfiledOverlap {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, OrientedSource, MAX_SOURCES_PER_STEP>"
        )]
        fragments: Vec<OrientedSource>,
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, MAX_SOURCES_PER_STEP>"
        )]
        overlaps: Vec<usize>,
        circular: bool,
        options: crate::cloning::overlap::profiles::ProfileOptions,
    },
    /// Explicit closed Dam/Dcm/CpG generation model on an original immutable input.
    DeclareSourceModifications {
        template: Source,
        assumptions: GenerationAssumptions,
    },
    /// Positional modification-aware cleavage under an explicitly pinned profile.
    CheckedDigest {
        template: Source,
        selection: crate::cloning::restriction::catalog::CatalogSelection,
        options: crate::cloning::restriction::inspection::InspectionOptions,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        conditions: Option<crate::cloning::restriction::conditions::ReactionConditions>,
        fragment_index: usize,
    },
    /// Bounded reference-sequence exchange only. No chemical duplex state is established.
    RecombineReference {
        sources: [Source; 2],
        request: crate::cloning::recombination::Request,
        selection: crate::cloning::recombination::Selection,
    },
    /// Fresh full physical graph. A proven contour is only a restrictive reference
    /// view; no target acceptance or downstream chemical state is inferred.
    Topo {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, TopoSource, MAX_SOURCES_PER_STEP>"
        )]
        sources: Vec<TopoSource>,
        request: crate::cloning::topo::Request,
        selection: crate::cloning::topo::Selection,
    },
    /// Final-only independent-target acceptance under an explicit physical policy.
    /// Accepted reference coordinates never establish downstream duplex authority.
    CheckedTopo {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, TopoSource, MAX_SOURCES_PER_STEP>"
        )]
        sources: Vec<TopoSource>,
        request: crate::cloning::topo::Request,
        selection: crate::cloning::topo::Selection,
        policy: crate::cloning::topo::AcceptancePolicy,
    },
    /// Operator hypothesis for an original intact linear duplex only; never a
    /// measurement, synthesis, repair or restoration of imported physical state.
    DeclareSourceTermini {
        template: Source,
        left: chemistry::Terminus,
        right: chemistry::Terminus,
    },
    /// Complete cleavage under an exact embedded-catalog identity/geometry pin.
    /// Supplier documentary data never establishes reaction completion.
    CatalogDigest {
        template: Source,
        selection: crate::cloning::restriction::catalog::CatalogSelection,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        conditions: Option<crate::cloning::restriction::conditions::ReactionConditions>,
        fragment_index: usize,
    },
    /// Explicit path identity checked against the same fresh native enumeration.
    /// Neither an imported report nor correspondence to the target chooses a path.
    SelectedGoldenGate {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, MAX_SOURCES_PER_STEP>"
        )]
        templates: Vec<Source>,
        options: crate::cloning::golden_gate::GoldenGateOptions,
        candidate_index: usize,
        directed_sequence_key: String,
        physical_sequence_key: String,
    },
    /// Accounted explicit graph actions; includes every retained byproduct.
    TopoActions {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, TopoSource, MAX_SOURCES_PER_STEP>"
        )]
        sources: Vec<TopoSource>,
        request: crate::cloning::topo::actions::Request,
    },
    /// Complete declared multi-parent reference network; no physical duplex authority.
    MultiSiteRecombineReference {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, TopoSource, MAX_SOURCES_PER_STEP>"
        )]
        sources: Vec<TopoSource>,
        request: crate::cloning::recombination::multisite::Request,
        selection: crate::cloning::recombination::Selection,
    },
    /// Reexecute a saved product's closed origin scope. Metadata on the selected
    /// input never restores material; the private tracked state is freshly built.
    ReconstructRunProduct {
        template: Source,
        origin: crate::material::capsule::Capsule,
    },
    /// Explicit reciprocal polymerase extension of two intact overlapping duplexes.
    /// Products arise from the selected physical strands, never the intended target.
    OverlapExtension {
        left: OrientedSource,
        right: OrientedSource,
        options: crate::cloning::overlap::extension::Options,
    },
    /// Fresh physical resection/pairing/fill/seal under a versioned declared model.
    /// Historical ProfiledOverlap remains a sequence-reference operation.
    HomologyAssembly {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, OrientedSource, MAX_SOURCES_PER_STEP>"
        )]
        fragments: Vec<OrientedSource>,
        options: crate::cloning::overlap::profiles::ProfileOptions,
        policy: crate::cloning::overlap::assembly::Options,
    },
    /// Explicit hydrolytic3′ shortening of an actual intact tracked operand.
    ResectThreePrime {
        source: Source,
        cuts: crate::cloning::overlap::resection::Cuts,
    },
    /// Generic3′ processing/exposed5′tail annealing only. Every joined strand
    /// remains separate; no polymerase, ligase or cellular rescue is inferred.
    HomologyAnnealing {
        #[serde(
            deserialize_with = "collagene_model::decode::bounded_vec::<_, OrientedSource, MAX_SOURCES_PER_STEP>"
        )]
        fragments: Vec<OrientedSource>,
        options: crate::cloning::overlap::profiles::ProfileOptions,
        policy: crate::cloning::overlap::annealing::Options,
    },

    /// Replay a target-free typed derivation; serialized summaries confer no material.
    ReconstructMaterial {
        template: Source,
        origin: crate::material::derivation::Capsule,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DesignPlan {
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub standard: Option<crate::cloning::assembly::standards::declaration::Declaration>,
    pub name: String,
    pub target_revision: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, PlanInput, MAX_INPUTS>")]
    pub inputs: Vec<PlanInput>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Operation, MAX_STEPS>")]
    pub steps: Vec<Operation>,
}

impl Operation {
    /// Actual typed dependencies; producer metadata never supplies graph authority.
    pub fn sources(&self) -> Vec<&Source> {
        match self {
            Operation::RecombineReference { sources, .. } => sources.iter().collect(),
            Operation::Topo { sources, .. }
            | Operation::CheckedTopo { sources, .. }
            | Operation::TopoActions { sources, .. }
            | Operation::MultiSiteRecombineReference { sources, .. } => {
                sources.iter().map(|s| &s.source).collect()
            }
            Operation::Pcr { template, .. }
            | Operation::PcrWithMismatches { template, .. }
            | Operation::Digest { template, .. }
            | Operation::CatalogDigest { template, .. }
            | Operation::CheckedDigest { template, .. }
            | Operation::DeclareSourceModifications { template, .. }
            | Operation::DeclareSourceTermini { template, .. }
            | Operation::TreatEnds { template, .. }
            | Operation::RepairEnds { template, .. }
            | Operation::Tail3Prime { template, .. } => vec![template],
            Operation::ReconstructRunProduct { template, .. }
            | Operation::ReconstructMaterial { template, .. } => vec![template],
            Operation::ResectThreePrime { source, .. } => vec![source],
            Operation::Overlap { fragments, .. }
            | Operation::GoldenGate {
                templates: fragments,
                ..
            }
            | Operation::SelectedGoldenGate {
                templates: fragments,
                ..
            } => fragments.iter().collect(),
            Operation::Ligate { fragments, .. }
            | Operation::ProfiledOverlap { fragments, .. }
            | Operation::HomologyAssembly { fragments, .. }
            | Operation::HomologyAnnealing { fragments, .. } => {
                fragments.iter().map(|f| &f.source).collect()
            }
            Operation::OverlapExtension { left, right, .. } => vec![&left.source, &right.source],
            Operation::AnnealOligos { top, bottom, .. } => vec![top, bottom],
        }
    }

    /// Mutable typed dependency references only. Callers must revalidate ordering and
    /// source membership after edits; this visitor confers no scientific authority.
    pub fn sources_mut(&mut self) -> Vec<&mut Source> {
        match self {
            Operation::RecombineReference { sources, .. } => sources.iter_mut().collect(),
            Operation::Topo { sources, .. }
            | Operation::CheckedTopo { sources, .. }
            | Operation::TopoActions { sources, .. }
            | Operation::MultiSiteRecombineReference { sources, .. } => {
                sources.iter_mut().map(|s| &mut s.source).collect()
            }
            Operation::Pcr { template, .. }
            | Operation::PcrWithMismatches { template, .. }
            | Operation::Digest { template, .. }
            | Operation::CatalogDigest { template, .. }
            | Operation::CheckedDigest { template, .. }
            | Operation::DeclareSourceModifications { template, .. }
            | Operation::DeclareSourceTermini { template, .. }
            | Operation::TreatEnds { template, .. }
            | Operation::RepairEnds { template, .. }
            | Operation::Tail3Prime { template, .. } => vec![template],
            Operation::ReconstructRunProduct { template, .. }
            | Operation::ReconstructMaterial { template, .. } => vec![template],
            Operation::ResectThreePrime { source, .. } => vec![source],
            Operation::Overlap { fragments, .. }
            | Operation::GoldenGate {
                templates: fragments,
                ..
            }
            | Operation::SelectedGoldenGate {
                templates: fragments,
                ..
            } => fragments.iter_mut().collect(),
            Operation::Ligate { fragments, .. }
            | Operation::ProfiledOverlap { fragments, .. }
            | Operation::HomologyAssembly { fragments, .. }
            | Operation::HomologyAnnealing { fragments, .. } => {
                fragments.iter_mut().map(|f| &mut f.source).collect()
            }
            Operation::OverlapExtension { left, right, .. } => {
                vec![&mut left.source, &mut right.source]
            }
            Operation::AnnealOligos { top, bottom, .. } => vec![top, bottom],
        }
    }
}

impl DesignPlan {
    pub fn revision(&self) -> String {
        blake3::hash(&serde_json::to_vec(self).expect("serializable plan"))
            .to_hex()
            .to_string()
    }
}
