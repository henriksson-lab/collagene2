//! Typed plans are reconstructed from immutable inputs, never from claimed products.
use crate::{
    analysis::{self, Primer},
    annotations::{self, MappingSegment},
    chemistry::{self, Duplex, EndTreatment},
    model::*,
    modified_duplex::{self, CheckedDigestOptions, GenerationAssumptions, TrackedDuplex},
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PlanInput {
    pub id: String,
    pub revision: String,
    pub molecule: Molecule,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
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
    },
    /// An explicit ordered candidate, not an exhaustive assembly search.
    Overlap {
        fragments: Vec<Source>,
        overlaps: Vec<usize>,
        circular: bool,
    },
    /// Complete digestion; select a zero-based fragment in ascending top-cut order.
    Digest {
        template: Source,
        enzymes: Vec<String>,
        fragment_index: usize,
    },
    TreatEnds {
        template: Source,
        treatment: EndTreatment,
    },
    /// Declared intended ligation and bounded one-copy alternatives, with explicit orientation and nick acceptance.
    Ligate {
        fragments: Vec<OrientedSource>,
        circular: bool,
        allow_nicked: bool,
        options: crate::ligation_outcomes::Options,
    },
    PcrWithMismatches {
        template: Source,
        forward: Primer,
        reverse: Primer,
        options: crate::mutagenesis::PcrOptions,
    },
    GoldenGate {
        templates: Vec<Source>,
        options: crate::golden_gate::GoldenGateOptions,
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
        base: crate::tailing::TailingBase,
    },
    ProfiledOverlap {
        fragments: Vec<OrientedSource>,
        overlaps: Vec<usize>,
        circular: bool,
        options: crate::overlap_profiles::ProfileOptions,
    },
    /// Explicit closed Dam/Dcm/CpG generation model on an original immutable input.
    DeclareSourceModifications {
        template: Source,
        assumptions: GenerationAssumptions,
    },
    /// Positional modification-aware cleavage under an explicitly pinned profile.
    CheckedDigest {
        template: Source,
        selection: crate::catalog_digest::CatalogSelection,
        options: crate::restriction_inspection::InspectionOptions,
        fragment_index: usize,
    },
    /// Bounded reference-sequence exchange only. No chemical duplex state is established.
    RecombineReference {
        sources: [Source; 2],
        request: crate::recombination::Request,
        selection: crate::recombination::Selection,
    },
    /// Fresh full physical graph. A proven contour is only a restrictive reference
    /// view; no target acceptance or downstream chemical state is inferred.
    Topo {
        sources: Vec<TopoSource>,
        request: crate::topo::Request,
        selection: crate::topo::Selection,
    },
    /// Final-only independent-target acceptance under an explicit physical policy.
    /// Accepted reference coordinates never establish downstream duplex authority.
    CheckedTopo {
        sources: Vec<TopoSource>,
        request: crate::topo::Request,
        selection: crate::topo::Selection,
        policy: crate::topo::AcceptancePolicy,
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
        selection: crate::catalog_digest::CatalogSelection,
        fragment_index: usize,
    },
    /// Explicit path identity checked against the same fresh native enumeration.
    /// Neither an imported report nor correspondence to the target chooses a path.
    SelectedGoldenGate {
        templates: Vec<Source>,
        options: crate::golden_gate::GoldenGateOptions,
        candidate_index: usize,
        directed_sequence_key: String,
        physical_sequence_key: String,
    },
    /// Fresh declared original sources, transfer/sealing graph and selected paired
    /// contour. Material requires a faithful bridge; no imported graph is input.
    ProcessedTopo {
        sources: Vec<TopoSource>,
        request: crate::topo_processed_product::Request,
    },
    /// Accounted explicit graph actions; includes every retained byproduct.
    TopoActions {
        sources: Vec<TopoSource>,
        request: crate::topo_actions::Request,
    },
    /// Complete declared multi-parent reference network; no physical duplex authority.
    MultiSiteRecombineReference {
        sources: Vec<TopoSource>,
        request: crate::multisite_recombination::Request,
        selection: crate::recombination::Selection,
    },
    /// Reexecute a saved product's closed origin scope. Metadata on the selected
    /// input never restores material; the private tracked state is freshly built.
    ReconstructRunProduct {
        template: Source,
        origin: crate::run_material::Capsule,
    },
    /// Explicit reciprocal polymerase extension of two intact overlapping duplexes.
    /// Products arise from the selected physical strands, never the intended target.
    OverlapExtension {
        left: OrientedSource,
        right: OrientedSource,
        options: crate::overlap_extension::Options,
    },
    /// Fresh physical resection/pairing/fill/seal under a versioned declared model.
    /// Historical ProfiledOverlap remains a sequence-reference operation.
    HomologyAssembly {
        fragments: Vec<OrientedSource>,
        options: crate::overlap_profiles::ProfileOptions,
        policy: crate::homology_assembly::Options,
    },
    /// Explicit hydrolytic3′ shortening of an actual intact tracked operand.
    ResectThreePrime {
        source: Source,
        cuts: crate::three_prime_resection::Cuts,
    },
    /// Generic3′ processing/exposed5′tail annealing only. Every joined strand
    /// remains separate; no polymerase, ligase or cellular rescue is inferred.
    HomologyAnnealing {
        fragments: Vec<OrientedSource>,
        options: crate::overlap_profiles::ProfileOptions,
        policy: crate::homology_annealing::Options,
    },

    /// Replay a target-free typed derivation; serialized summaries confer no material.
    ReconstructMaterial {
        template: Source,
        origin: crate::material_derivation::Capsule,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DesignPlan {
    pub name: String,
    pub target_revision: String,
    pub inputs: Vec<PlanInput>,
    pub steps: Vec<Operation>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Status {
    Invalid,
    ConsistentUnderAssumptions,
    Incomplete,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct VerificationReport {
    pub plan_revision: String,
    pub target_revision: String,
    pub engine_version: String,
    pub status: Status,
    pub findings: Vec<String>,
    pub assumptions: Vec<String>,
    pub products: Vec<Molecule>,
    pub exhaustive_alternative_search: bool,
    /// Full per-step scientific evidence, including withheld TOPO contours.
    #[serde(default)]
    pub step_reports: Vec<StepReport>,
    /// Conservative annotation correspondence, separate from physical graph evidence.
    #[serde(default)]
    pub annotation_reports: Vec<AnnotationStepReport>,
    /// Full additive processed graph evidence, including withheld material.
    #[serde(default)]
    pub processed_step_reports: Vec<ProcessedTopoStepReport>,
    #[serde(default)]
    pub physical_processing_reports: Vec<PhysicalProcessingStepReport>,
    /// Typed reference-design evidence, distinct from actual material state.
    #[serde(default)]
    pub reference_step_reports: Vec<ReferenceStepReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct StepReport {
    pub step_index: usize,
    #[serde(default)]
    pub topo: Option<crate::topo::Product>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AnnotationStepReport {
    pub step_index: usize,
    pub evidence: crate::topo_annotations::AnnotationEvidence,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProcessedTopoStepReport {
    pub step_index: usize,
    pub evidence: crate::topo_processed_product::Report,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum PhysicalProcessingEvidence {
    TopoActions(crate::topo_actions::Report),
    RunProduct(crate::run_material::ReplayEvidence),
    OverlapExtension(crate::overlap_extension::Geometry),
    HomologyAssembly(crate::homology_assembly::Report),
    ThreePrimeResection(crate::three_prime_resection::Report),
    HomologyAnnealing(crate::homology_annealing::Report),
    DerivedOrigin(DerivedOriginEvidence),
    LigationOutcomes(LigationOutcomeEvidence),
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct LigationOutcomeBinding {
    pub key: String,
    pub source: Source,
    pub source_revision: String,
    pub material_revision: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct LigationOutcomeEvidence {
    pub bindings: Vec<LigationOutcomeBinding>,
    pub report: crate::ligation_outcomes::Report,
}
impl LigationOutcomeEvidence {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        self.report
            .validate_binding(&self.report.context_revision)?;
        if self.bindings.len() != self.report.sources.len() {
            return Err(Error::Format(
                "Ligation outcome source count changed".into(),
            ));
        }
        for (i, binding) in self.bindings.iter().enumerate() {
            if binding.key != format!("fragment-{}", i + 1)
                || !self
                    .report
                    .sources
                    .iter()
                    .any(|s| s.key == binding.key && s.state_revision == binding.material_revision)
                || [&binding.source_revision, &binding.material_revision]
                    .iter()
                    .any(|s| s.len() != 64 || !s.bytes().all(|b| b.is_ascii_hexdigit()))
                || matches!(&binding.source,Source::Input(id) if id.is_empty() || id.len()>128 || id.chars().any(char::is_control))
            {
                return Err(Error::Format(
                    "Ligation outcome exact source/material binding changed".into(),
                ));
            }
        }
        evidence_bytes(self, 32 * 1024 * 1024)?;
        Ok(())
    }
    pub fn validate_for_operation(
        &self,
        fragments: &[OrientedSource],
        circular: bool,
        allow_nicked: bool,
        options: &crate::ligation_outcomes::Options,
    ) -> Result<(), Error> {
        self.validate_documentary()?;
        if fragments.len() != self.bindings.len()
            || self.report.request.circular != circular
            || self.report.request.allow_nicked != allow_nicked
            || crate::assembly_recipe::hash(options)?
                != crate::assembly_recipe::hash(&self.report.request.options)?
        {
            return Err(Error::Format(
                "Ligation outcomes declaration changed".into(),
            ));
        }
        for ((fragment, binding), placement) in fragments
            .iter()
            .zip(&self.bindings)
            .zip(&self.report.request.intended)
        {
            if crate::assembly_recipe::hash(&fragment.source)?
                != crate::assembly_recipe::hash(&binding.source)?
                || placement.key != binding.key
                || placement.reverse != fragment.reverse
            {
                return Err(Error::Format(
                    "Ligation intended path/source correspondence changed".into(),
                ));
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct MaterialStateEvidence {
    pub material_revision: String,
    pub duplex: chemistry::DuplexSummary,
    pub modifications: modified_duplex::ModificationSummary,
    pub origin_incomplete: bool,
    pub origin_findings: Vec<String>,
    pub origin_findings_total: usize,
    pub origin_findings_truncated: bool,
}
impl MaterialStateEvidence {
    /// Documentary shape and exact recipe binding, never imported material authority.
    pub fn validate_against(
        &self,
        expected_state: &str,
        expected_incomplete: bool,
    ) -> Result<(), Error> {
        let hex = |value: &str| {
            value.len() == 64
                && value
                    .bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        };
        let count = self
            .modifications
            .known_modified
            .checked_add(self.modifications.unknown_eligible)
            .and_then(|n| n.checked_add(self.modifications.known_unmodified))
            .ok_or(Error::TooLarge)?;
        let total = self
            .duplex
            .top_length
            .checked_add(self.duplex.bottom_length)
            .ok_or(Error::TooLarge)?;
        if !hex(&self.material_revision)
            || self.material_revision != self.modifications.state_revision
            || self.material_revision != expected_state
            || self.origin_incomplete != expected_incomplete
            || self.modifications.engine_version != modified_duplex::ENGINE_VERSION
            || self.duplex.top_length == 0
            || self.duplex.bottom_length == 0
            || self.duplex.top_length > crate::model::MAX_SEQUENCE_LEN
            || self.duplex.bottom_length > crate::model::MAX_SEQUENCE_LEN
            || self.modifications.top_length != self.duplex.top_length
            || self.modifications.bottom_length != self.duplex.bottom_length
            || count > total
            || self.origin_findings.len() > 128
            || self
                .origin_findings
                .iter()
                .any(|s| s.len() > 8192 || s.chars().count() > 2048)
            || self.origin_findings_total < self.origin_findings.len()
            || (!self.origin_findings_truncated
                && self.origin_findings_total > self.origin_findings.len())
        {
            return Err(Error::Format(
                "Derived material evidence changed its bounded source/material binding".into(),
            ));
        }
        evidence_bytes(self, 32 * 1024 * 1024)?;
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DerivedOriginEvidence {
    pub origin: crate::material_derivation::Capsule,
    pub material: MaterialStateEvidence,
}
impl DerivedOriginEvidence {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        self.origin.validate_documentary()?;
        let record = crate::material_derivation::load_record_shallow(&self.origin.record_bytes)?;
        self.material
            .validate_against(&record.product_state_revision, record.incomplete)?;
        evidence_bytes(self, 32 * 1024 * 1024)?;
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PhysicalProcessingStepReport {
    pub step_index: usize,
    pub evidence: PhysicalProcessingEvidence,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ReferenceStepEvidence {
    MultiSiteRecombination {
        selection: crate::recombination::Selection,
        batch: crate::multisite_recombination::Batch,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ReferenceStepReport {
    pub step_index: usize,
    pub evidence: ReferenceStepEvidence,
}

pub const ENGINE_VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "/science-27");

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PcrDocument {
    pub molecule: Molecule,
    pub geometry: analysis::PcrProduct,
    pub findings: Vec<String>,
}

/// Enumerate sequence-reference previews. These are not target verification reports or
/// evidence of physical end chemistry; the immutable source remains the parent document.
pub fn simulate_pcr_documents(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
) -> Result<Vec<PcrDocument>, Error> {
    source.validate()?;
    if chemistry::has_unmapped_physical_evidence(source) {
        return Err(Error::Unsupported("PCR of imported physical-intermediate evidence requires reconstruction of its originating graph".into()));
    }
    build_pcr_documents(source, forward, reverse)
}

/// Only fresh native run reconstruction may bypass the imported-marker guard.
pub(crate) fn simulate_pcr_documents_from_fresh(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    tracked: &modified_duplex::TrackedDuplex,
) -> Result<Vec<PcrDocument>, Error> {
    source.validate()?;
    tracked.validate()?;
    pcr_material_preflight(source, tracked, forward, reverse)?;
    let mut documents = build_pcr_documents(source, forward, reverse)?;
    if tracked.duplex.bottom_start != 0 || tracked.duplex.bottom.len() != tracked.duplex.top.len() {
        for document in &mut documents {
            document.findings.push("PCR uses the actual intact top strand as a first-cycle template. Its initially staggered complement is not silently filled or repaired; both-strand coverage was checked for additional geometries outside the supported top reference.".into());
        }
    }
    Ok(documents)
}
/// Exact PCR may seed from an intact covalent top strand even when its initial
/// complement is staggered. Search the known strand union first: extra geometries
/// outside the top reference are explicit unsupported work, never false absence.
pub(crate) fn pcr_material_preflight(
    source: &Molecule,
    tracked: &modified_duplex::TrackedDuplex,
    forward: &Primer,
    reverse: &Primer,
) -> Result<(), Error> {
    source.validate()?;
    tracked.validate()?;
    let duplex = &tracked.duplex;
    if duplex.top != source.sequence
        || duplex.topology != source.topology
        || !duplex.nicks.is_empty()
    {
        return Err(Error::Unsupported("PCR requires reference-matching intact covalent template strands; nick propagation is not modeled".into()));
    }
    if duplex.bottom_start == 0 && duplex.bottom.len() == duplex.top.len() {
        return Ok(());
    }
    if duplex.topology != Topology::Linear {
        return Err(Error::Unsupported(
            "Staggered circular PCR templates require an explicit supported linearization".into(),
        ));
    }
    let bottom_end = duplex
        .bottom_start
        .checked_add(i64::try_from(duplex.bottom.len()).map_err(|_| Error::TooLarge)?)
        .ok_or(Error::TooLarge)?;
    let start = 0i64.min(duplex.bottom_start);
    let end = (duplex.top.len() as i64).max(bottom_end);
    let length = usize::try_from(end.checked_sub(start).ok_or(Error::TooLarge)?)
        .map_err(|_| Error::TooLarge)?;
    if length > crate::model::MAX_SEQUENCE_LEN {
        return Err(Error::TooLarge);
    }
    let bottom_plus = duplex.bottom.reverse_complement()?;
    let mut bases = Vec::with_capacity(length);
    for i in 0..length {
        let coordinate = start.checked_add(i as i64).ok_or(Error::TooLarge)?;
        let base = if coordinate >= 0 && (coordinate as usize) < duplex.top.len() {
            duplex.top.as_str().as_bytes()[coordinate as usize]
        } else {
            let p = coordinate
                .checked_sub(duplex.bottom_start)
                .and_then(|n| usize::try_from(n).ok())
                .filter(|p| *p < bottom_plus.len())
                .ok_or_else(|| {
                    Error::Unsupported("PCR strand coverage contains an unmodeled gap".into())
                })?;
            bottom_plus.as_str().as_bytes()[p]
        };
        bases.push(base);
    }
    let projection = Molecule::dna(
        "Actual PCR strand-coverage projection",
        std::str::from_utf8(&bases)
            .map_err(|_| Error::Format("PCR strand coverage is not ASCII".into()))?,
        Topology::Linear,
    )?;
    let offset = usize::try_from(start.checked_neg().ok_or(Error::TooLarge)?)
        .map_err(|_| Error::TooLarge)?;
    for geometry in analysis::simulate_pcr_geometry(&projection, forward, reverse)? {
        let begin = geometry.product.forward_start;
        let stop = begin
            .checked_add(geometry.template_length)
            .ok_or(Error::TooLarge)?;
        if begin < offset
            || stop
                > offset
                    .checked_add(duplex.top.len())
                    .ok_or(Error::TooLarge)?
        {
            return Err(Error::Unsupported("Actual staggered strands admit a PCR geometry outside the top-reference-contained profile. Inspect the strand extensions; no complete absence, unique product, or silent repair is asserted".into()));
        }
    }
    Ok(())
}

fn build_pcr_documents(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
) -> Result<Vec<PcrDocument>, Error> {
    let geometries = analysis::simulate_pcr_geometry(source, forward, reverse)?;
    let mut total = 0usize;
    let mut documents = Vec::new();
    for (index, geometry) in geometries.into_iter().enumerate() {
        let mut molecule = Molecule::dna(
            format!("{}_PCR_{}", source.name, index + 1),
            geometry.product.sequence.as_str(),
            Topology::Linear,
        )?;
        let mut findings = Vec::new();
        project_pcr_annotations(source, &geometry, &mut molecule, &mut findings)?;
        molecule
            .metadata
            .insert("pcr:template_revision".into(), source.revision());
        molecule.metadata.insert(
            "pcr:primers".into(),
            serde_json::to_string(&(forward, reverse)).map_err(|e| Error::Format(e.to_string()))?,
        );
        molecule.metadata.insert("pcr:model".into(), "Exact annealing reference preview; alternative geometries remain separate. No target verification or physical chemistry is asserted.".into());
        molecule
            .metadata
            .insert("pcr:engine_version".into(), ENGINE_VERSION.into());
        let size = serde_json::to_vec(&molecule)
            .map_err(|e| Error::Format(e.to_string()))?
            .len();
        total = total.checked_add(size).ok_or(Error::TooLarge)?;
        if total > crate::formats::MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
        documents.push(PcrDocument {
            molecule,
            geometry: geometry.product,
            findings,
        });
    }
    Ok(documents)
}

pub(crate) fn project_pcr_annotations(
    source: &Molecule,
    geometry: &analysis::PcrGeometryProduct,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<(), Error> {
    let map = region_map(
        source.sequence.len(),
        geometry.product.forward_start,
        geometry.template_length,
        geometry.forward_tail_length,
    )?;
    project_into(source, &map, product, findings)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MutagenesisPreview {
    pub report: crate::mutagenesis::PcrReport,
    pub documents: Vec<Molecule>,
    pub findings: Vec<String>,
}
pub fn simulate_mutagenesis_preview(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::mutagenesis::PcrOptions,
) -> Result<MutagenesisPreview, Error> {
    if chemistry::has_unmapped_physical_evidence(source) {
        return Err(Error::Unsupported(
            "PCR preview of physical-intermediate evidence requires its originating graph".into(),
        ));
    }
    let report = crate::mutagenesis::simulate_pcr(source, forward, reverse, options)?;
    project_mutagenesis_preview(source, forward, reverse, options, report)
}
/// Project the one freshly computed native report, preserving its full search evidence.
pub(crate) fn project_mutagenesis_preview(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::mutagenesis::PcrOptions,
    report: crate::mutagenesis::PcrReport,
) -> Result<MutagenesisPreview, Error> {
    let mut documents = Vec::new();
    let mut findings = Vec::new();
    let mut bytes = 0usize;
    for (i, geometry) in report.products.iter().enumerate() {
        let mut product = Molecule::dna(
            format!("{}_mutant_PCR_{}", source.name, i + 1),
            geometry.sequence.as_str(),
            Topology::Linear,
        )?;
        project_mutant_pcr(source, geometry, &mut product, &mut findings)?;
        product
            .metadata
            .insert("pcr:template_revision".into(), source.revision());
        product.metadata.insert(
            "pcr:mismatch_options".into(),
            serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?,
        );
        product.metadata.insert(
            "pcr:primers".into(),
            serde_json::to_string(&(forward, reverse)).map_err(|e| Error::Format(e.to_string()))?,
        );
        product.metadata.insert("pcr:model".into(),"Intentional-mismatch reference preview; no target verification, amplification or physical end chemistry is asserted. Retain the source document and report.".into());
        product
            .metadata
            .insert("pcr:engine_version".into(), ENGINE_VERSION.into());
        bytes = bytes.saturating_add(
            serde_json::to_vec(&product)
                .map_err(|e| Error::Format(e.to_string()))?
                .len(),
        );
        if bytes > 16 * 1024 * 1024 {
            return Err(Error::TooLarge);
        }
        documents.push(product);
    }
    Ok(MutagenesisPreview {
        report,
        documents,
        findings,
    })
}
pub(crate) fn project_mutant_pcr(
    source: &Molecule,
    geometry: &crate::mutagenesis::MutantPcrProduct,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<(), Error> {
    let map = region_map(
        source.sequence.len(),
        geometry.forward_start,
        geometry.template_length,
        geometry.forward_tail_length,
    )?;
    project_into(source, &map, product, findings)?;
    let mut changed = geometry
        .substitutions
        .iter()
        .map(|change| change.product_position)
        .collect::<Vec<_>>();
    changed.sort_unstable();
    changed.dedup();
    for feature in &mut product.features {
        if feature.location.spans().iter().any(|(a, b, _)| {
            changed
                .get(changed.partition_point(|position| *position < *a))
                .is_some_and(|position| *position < *b)
        }) {
            feature.needs_review = true;
            let mut evidence = Vec::new();
            feature.qualifiers.retain(|(key, value)| {
                if key == "translation" || key == "codon_start" {
                    evidence.push(format!(
                        "Source assertion /{key}={}",
                        value.as_deref().unwrap_or("")
                    ));
                    false
                } else {
                    true
                }
            });
            feature.qualifiers.push(("note".into(),Some("Primer-introduced changes overlap this annotation; review its sequence and function.".into())));
            for note in evidence {
                feature.qualifiers.push(("note".into(), Some(note)));
            }
        }
    }
    product.validate()
}

pub fn plan_revision(plan: &DesignPlan) -> String {
    blake3::hash(&serde_json::to_vec(plan).expect("serializable plan"))
        .to_hex()
        .to_string()
}

struct EvidenceCounter {
    bytes: usize,
    max: usize,
}
impl std::io::Write for EvidenceCounter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|&n| n <= self.max)
            .ok_or_else(|| std::io::Error::other("bounded plan step evidence"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn evidence_bytes(value: &impl Serialize, max: usize) -> Result<usize, Error> {
    let mut counter = EvidenceCounter { bytes: 0, max };
    serde_json::to_writer(&mut counter, value).map_err(|_| Error::TooLarge)?;
    Ok(counter.bytes)
}

/// Cheap borrowed snapshot/resource guard before plan hashing or source clones.
/// Actual sequence intervals, strand pairing and scalar engine policies are
/// independently checked by the source adapter against freshly bound molecules.
fn validate_topo_shape(
    sources: &[TopoSource],
    request: &crate::topo::Request,
) -> Result<(), Error> {
    if sources.is_empty()
        || sources.len() > 32
        || request.recipe.schema_version != 1
        || request.recipe.materials.is_empty()
        || request.recipe.materials.len() > 16
        || request.transfers.len() > 16
    {
        return Err(Error::TooLarge);
    }
    let mut roles = std::collections::BTreeSet::new();
    for source in sources {
        if source.role.len() > 256
            || source.role.trim().is_empty()
            || !roles.insert(source.role.as_str())
            || matches!(&source.source,Source::Input(id) if id.len()>128 || id.trim().is_empty())
        {
            return Err(Error::Format(
                "invalid or duplicate TOPO source role".into(),
            ));
        }
    }
    let mut strands = 0usize;
    let mut pairing = 0usize;
    let mut segments = 0usize;
    for material in &request.recipe.materials {
        match material {
            crate::topo::MaterialRecipe::Prepared {
                declaration_id,
                assumption,
                strands: items,
                pairing: spans,
            } => {
                strands = strands.checked_add(items.len()).ok_or(Error::TooLarge)?;
                pairing = pairing.checked_add(spans.len()).ok_or(Error::TooLarge)?;
                if strands > 32
                    || pairing > 1024
                    || declaration_id.len() > 256
                    || assumption.len() > 4096
                {
                    return Err(Error::TooLarge);
                }
                for strand in items {
                    if strand.id.len() > 256
                        || strand.segments.is_empty()
                        || strand.segments.len() > 1024
                    {
                        return Err(Error::TooLarge);
                    }
                    if let crate::topo::physical::ThreePrime::Phosphotyrosyl {
                        activation_id, ..
                    } = &strand.three_prime
                    {
                        if activation_id.len() > 256 {
                            return Err(Error::TooLarge);
                        }
                    }
                    segments = segments
                        .checked_add(strand.segments.len())
                        .ok_or(Error::TooLarge)?;
                    if segments > 2048 {
                        return Err(Error::TooLarge);
                    }
                    for segment in &strand.segments {
                        if segment.source_role.len() > 256
                            || !roles.contains(segment.source_role.as_str())
                        {
                            return Err(Error::Format(
                                "TOPO interval uses an unbound source role".into(),
                            ));
                        }
                    }
                }
            }
            crate::topo::MaterialRecipe::ExactPcr {
                material_id,
                template_role,
                forward,
                reverse,
                ..
            } => {
                strands = strands.checked_add(2).ok_or(Error::TooLarge)?;
                if strands > 32
                    || material_id.len() > 256
                    || template_role.len() > 256
                    || !roles.contains(template_role.as_str())
                {
                    return Err(Error::Format(
                        "TOPO PCR uses an unbound role or exceeds preparation bounds".into(),
                    ));
                }
                for primer in [forward, reverse] {
                    if primer.name.len() > 256
                        || primer.annealing.len() > 1000
                        || primer.tail.len() > 1000
                    {
                        return Err(Error::TooLarge);
                    }
                }
            }
        }
    }
    evidence_bytes(request, 1024 * 1024)?;
    Ok(())
}

/// Borrowed literal request bounds; actual source geometry is independently
/// reconstructed only after source models and aggregate resource caps pass.
fn validate_processed_topo_shape(
    sources: &[TopoSource],
    r: &crate::topo_processed_product::Request,
) -> Result<(), Error> {
    validate_topo_shape(sources, &r.treatment.preparation)?;
    let original = &r.treatment.preparation;
    crate::topo::physical::nick_sealing::request_preflight(
        &r.treatment.sealing,
        &original.engine_limits,
    )
    .map_err(|e| Error::Unsupported(e.to_string()))?;
    let l = &r.options.limits;
    let hash = &r.selection.sealing_context_revision;
    if original
        .transfers
        .len()
        .checked_add(r.treatment.sealing.nicks.len())
        .is_none_or(|n| n > 16 || n > original.engine_limits.max_steps)
        || hash.len() != 64
        || !hash
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        || r.selection.component >= 64
        || r.selection.reference_seed.strand >= 32
        || r.selection.reference_seed.index >= 100_000
        || l.max_pairs == 0
        || l.max_pairs > 500_000
        || l.max_work_cells == 0
        || l.max_work_cells > 128_000_000
        || l.max_work_bytes == 0
        || l.max_work_bytes > 128 * 1024 * 1024
        || l.max_mapping_segments == 0
        || l.max_mapping_segments > 16384
        || l.max_report_bytes == 0
        || l.max_report_bytes > 32 * 1024 * 1024
    {
        return Err(Error::TooLarge);
    }
    evidence_bytes(r, 2 * 1024 * 1024)?;
    Ok(())
}

/// Bind the exact ordered plan sources to explicit recipe roles. Display names
/// and imported producer metadata cannot resolve a role or establish chemistry.
pub fn topo_inputs(
    sources: &[TopoSource],
    molecules: &[&Molecule],
    request: &crate::topo::Request,
) -> Result<Vec<crate::topo::Input>, Error> {
    validate_topo_shape(sources, request)?;
    if sources.len() != molecules.len() {
        return Err(Error::Format("TOPO source count mismatch".into()));
    }
    collagene_model::bounds::model_preflight(molecules)?;
    let mut bases = 0usize;
    for molecule in molecules {
        bases = bases
            .checked_add(molecule.sequence.len())
            .filter(|&n| n <= 100_000)
            .ok_or(Error::TooLarge)?;
        if molecule.sequence.alphabet() != Alphabet::Dna {
            return Err(Error::Format("TOPO sources must be DNA".into()));
        }
    }
    let inputs = sources
        .iter()
        .zip(molecules)
        .map(|(source, molecule)| crate::topo::Input {
            role: source.role.clone(),
            molecule: (*molecule).clone(),
        })
        .collect::<Vec<_>>();
    crate::topo::source_preflight(&inputs, request)
        .map_err(|error| Error::Format(format!("TOPO source preparation rejected: {error}")))?;
    Ok(inputs)
}

/// A declaration may constrain an original intact reference hypothesis only.
/// This guard never interprets producer summaries as chemical authority.
pub fn declared_source_termini_preflight(source: &Molecule) -> Result<(), Error> {
    crate::model::bounds::model_preflight(&[source])?;
    source.validate()?;
    if source.topology != Topology::Linear
        || source.sequence.alphabet() != Alphabet::Dna
        || source
            .sequence
            .as_str()
            .bytes()
            .any(|b| !b"ACGT".contains(&b))
    {
        return Err(Error::Unsupported("Declared source termini require an original intact linear unambiguous DNA reference; topology or missing calls are not repaired.".into()));
    }
    let raw_producer_marker = source.original_genbank.as_ref().is_some_and(|raw| {
        raw.lines()
            .filter_map(|line| {
                line.strip_prefix("COMMENT     ")
                    .or_else(|| line.strip_prefix("CC   "))
            })
            .any(|value| {
                matches!(
                    value.trim_end(),
                    "Collagene physical product v1: unverified sequence-only evidence"
                        | "Collagene reference-only v1: unverified sequence-only evidence"
                )
            })
    });
    if raw_producer_marker
        || chemistry::has_unmapped_physical_evidence(source)
        || source.metadata.keys().any(|key| {
            matches!(
                key.as_str(),
                "tracked_modifications"
                    | "source:declared_termini"
                    | "source:legacy:chemistry_unmapped"
                    | "source:snapgene:chemistry_unmapped"
                    | "topo:physical_graph_required"
                    | "source:recombination:reference_only"
                    | "source:recombination:unverified_recipe"
                    | "source:restriction_inspection:recipe"
            ) || key == crate::recombination::REFERENCE_ONLY_KEY
                || key == crate::recombination::RECIPE_KEY
                || key == crate::restriction_inspection::RECIPE_KEY
                || key.starts_with("variant:")
                || key.starts_with("topo:")
                || key.starts_with("source:restriction_inspection:")
                || key.starts_with("restriction_catalog:")
                || key.starts_with("source:restriction_catalog:")
        })
    {
        return Err(Error::Unsupported("Source terminal declarations cannot promote imported physical, mixed/single-strand, variant or reference-only evidence. Reconstruct the actual originating graph; this declaration is not materialization, nick repair or chemical recovery.".into()));
    }
    Ok(())
}

pub(crate) fn validate_reference_steps(plan: &DesignPlan) -> Result<(), Error> {
    crate::run_material::plan_preflight(plan)?;
    if plan.steps.len() > 1000 || plan.inputs.len() > 1000 {
        return Err(Error::TooLarge);
    }
    let pinned_catalog = if plan.steps.iter().any(|s| {
        matches!(
            s,
            Operation::CatalogDigest { .. } | Operation::CheckedDigest { .. }
        )
    }) {
        Some(crate::catalog_digest::builtin_catalog()?)
    } else {
        None
    };
    let mut count = 0usize;
    for step in &plan.steps {
        if let Operation::CheckedDigest {
            selection,
            options,
            fragment_index,
            ..
        } = step
        {
            selection.validate_for(
                pinned_catalog
                    .as_ref()
                    .ok_or_else(|| Error::Format("missing checked catalog".into()))?,
            )?;
            options.validate()?;
            if *fragment_index >= options.max_fragments {
                return Err(Error::TooLarge);
            }
        }
        if let Operation::ResectThreePrime { cuts, .. } = step {
            count += 1;
            if count > 64
                || cuts.left > crate::three_prime_resection::MAX_SOURCE_BASES / 2
                || cuts.right > crate::three_prime_resection::MAX_SOURCE_BASES / 2
            {
                return Err(Error::TooLarge);
            }
        }
        if let Operation::HomologyAnnealing {
            fragments,
            options,
            policy,
        } = step
        {
            count += 1;
            if count > 64 {
                return Err(Error::TooLarge);
            }
            policy.validate_shape(fragments.len())?;
            options.validate()?;
            if fragments.len() == 1 && policy.topology != Topology::Circular {
                return Err(Error::Unsupported("Single linear3′ processing uses ResectThreePrime, not a joining assembly operation".into()));
            }
            if fragments
                .iter()
                .any(|f| matches!(&f.source,Source::Input(id)if id.trim().is_empty()||id.len()>128))
            {
                return Err(Error::Format("Invalid annealing source identity".into()));
            }
        }
        if let Operation::HomologyAssembly {
            fragments,
            options,
            policy,
        } = step
        {
            count += 1;
            if count > 64 {
                return Err(Error::TooLarge);
            }
            policy.validate_shape(fragments.len())?;
            options.validate()?;
            if options.method == crate::overlap_profiles::OverlapMethod::InFusion {
                return Err(Error::Unsupported("The declared 5′ resection/fill/seal model cannot represent In-Fusion 3′ processing and host rescue; retain its reference-only interpretation.".into()));
            }
            if fragments.iter().any(
                |f| matches!(&f.source, Source::Input(id) if id.trim().is_empty() || id.len()>128),
            ) {
                return Err(Error::Format(
                    "Invalid physical homology source identity".into(),
                ));
            }
        }
        if let Operation::SelectedGoldenGate {
            templates,
            options,
            candidate_index,
            directed_sequence_key,
            physical_sequence_key,
        } = step
        {
            if templates.is_empty()
                || templates.len() > 32
                || options.min_fragments == 0
                || options.min_fragments > options.max_fragments
                || options.max_fragments > 32
                || options.max_products == 0
                || options.max_products > 1000
                || options.max_states == 0
                || options.max_states > 1_000_000
                || *candidate_index >= options.max_products
                || [directed_sequence_key, physical_sequence_key]
                    .iter()
                    .any(|key| {
                        key.len() != 64
                            || !key
                                .bytes()
                                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
                    })
                || templates.iter().any(
                    |source| matches!(source, Source::Input(id) if id.is_empty() || id.len() > 128),
                )
            {
                return Err(Error::Format(
                    "Invalid bounded key-bound Golden Gate selection".into(),
                ));
            }
        }
        if let Operation::CatalogDigest {
            template,
            selection,
            fragment_index,
        } = step
        {
            selection
                .validate_for(pinned_catalog.as_ref().ok_or(Error::TooLarge)?)
                .map_err(|error| Error::Format(format!("Catalog digest pin rejected: {error}")))?;
            if *fragment_index > 100_000
                || matches!(template,Source::Input(id) if id.is_empty()||id.len()>128)
            {
                return Err(Error::TooLarge);
            }
        }
        if let Operation::DeclareSourceTermini {
            template: Source::Input(id),
            ..
        } = step
        {
            if id.is_empty() || id.len() > 128 {
                return Err(Error::TooLarge);
            }
        }
        if let Operation::Topo {
            sources,
            request,
            selection,
        }
        | Operation::CheckedTopo {
            sources,
            request,
            selection,
            ..
        } = step
        {
            count += 1;
            if count > 64 {
                return Err(Error::TooLarge);
            }
            validate_topo_shape(sources, request)?;
            if selection.batch_revision.len() != 64
                || !selection
                    .batch_revision
                    .bytes()
                    .all(|b| b.is_ascii_hexdigit())
                || selection.component >= 64
            {
                return Err(Error::Format("invalid TOPO source selection".into()));
            }
        }
        if let Operation::MultiSiteRecombineReference {
            sources,
            request,
            selection,
        } = step
        {
            count += 1;
            if count > 64 || !(2..=5).contains(&sources.len()) {
                return Err(Error::TooLarge);
            }
            let mut roles = BTreeSet::new();
            for source in sources {
                if source.role.trim().is_empty()
                    || source.role.len() > 256
                    || !roles.insert(source.role.as_str())
                {
                    return Err(Error::Format(
                        "Multi-site roles must be unique and bounded".into(),
                    ));
                }
            }
            crate::multisite_recombination::request_preflight(request)?;
            if selection.batch_revision.len() != 64
                || !selection
                    .batch_revision
                    .bytes()
                    .all(|b| b.is_ascii_hexdigit())
                || selection.candidate_id >= request.reaction.limits.max_candidates
                || selection.component_id >= request.reaction.limits.max_components_per_candidate
            {
                return Err(Error::Format("Invalid multi-site selection".into()));
            }
        }
        if let Operation::TopoActions { sources, request } = step {
            count += 1;
            if count > 64 {
                return Err(Error::TooLarge);
            }
            validate_topo_shape(sources, &request.preparation)?;
            crate::topo_actions::request_preflight(request)?;
        }
        if let Operation::ProcessedTopo { sources, request } = step {
            count += 1;
            if count > 64 {
                return Err(Error::TooLarge);
            }
            validate_processed_topo_shape(sources, request)?;
        }
        if let Operation::RecombineReference {
            sources,
            request,
            selection,
        } = step
        {
            count += 1;
            if count > 64
                || selection.batch_revision.len() != 64
                || !selection
                    .batch_revision
                    .bytes()
                    .all(|b| b.is_ascii_hexdigit())
            {
                return Err(Error::TooLarge);
            }
            request.validate()?;
            if selection.candidate_id >= request.limits.max_candidates
                || selection.component_id >= request.limits.max_components_per_candidate
                || sources.iter().any(
                    |s| matches!(s, Source::Input(id) if id.len() > 128 || id.trim().is_empty()),
                )
            {
                return Err(Error::Format(
                    "invalid reference-exchange selection or source shape".into(),
                ));
            }
        }
    }
    if count > 0
        || plan.steps.iter().any(|s| {
            matches!(
                s,
                Operation::DeclareSourceTermini { .. }
                    | Operation::CatalogDigest { .. }
                    | Operation::SelectedGoldenGate { .. }
                    | Operation::HomologyAssembly { .. }
                    | Operation::ResectThreePrime { .. }
                    | Operation::HomologyAnnealing { .. }
            )
        })
    {
        if plan.name.len() > 1024 || plan.target_revision.len() != 64 {
            return Err(Error::TooLarge);
        }
        let mut bases = 0usize;
        for input in &plan.inputs {
            if input.id.len() > 128 || input.revision.len() != 64 {
                return Err(Error::TooLarge);
            }
            bases = bases
                .checked_add(input.molecule.sequence.len())
                .filter(|n| *n <= 32_000_000)
                .ok_or(Error::TooLarge)?;
        }
        crate::recombination::reference::model_preflight(
            &plan.inputs.iter().map(|i| &i.molecule).collect::<Vec<_>>(),
        )?;
    }
    Ok(())
}

/// Actual typed dependencies; producer metadata never supplies graph authority.
pub fn operation_sources(op: &Operation) -> Vec<&Source> {
    match op {
        Operation::RecombineReference { sources, .. } => sources.iter().collect(),
        Operation::Topo { sources, .. }
        | Operation::CheckedTopo { sources, .. }
        | Operation::ProcessedTopo { sources, .. }
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
pub fn operation_sources_mut(op: &mut Operation) -> Vec<&mut Source> {
    match op {
        Operation::RecombineReference { sources, .. } => sources.iter_mut().collect(),
        Operation::Topo { sources, .. }
        | Operation::CheckedTopo { sources, .. }
        | Operation::ProcessedTopo { sources, .. }
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

/// Roles bind ordered slots and typed ancestry, independently of molecule display names.
pub fn recombination_inputs(
    sources: &[Source; 2],
    molecules: [&Molecule; 2],
    request: &crate::recombination::Request,
) -> Result<Vec<crate::recombination::reference::Input>, Error> {
    crate::recombination::source_preflight(molecules, request)?;
    if sources
        .iter()
        .any(|s| matches!(s, Source::Input(id) if id.len() > 128 || id.trim().is_empty()))
    {
        return Err(Error::TooLarge);
    }
    sources
        .iter()
        .zip(molecules)
        .enumerate()
        .map(|(slot, (source, molecule))| {
            let bytes = serde_json::to_vec(&("plan-reference-exchange-source-v1", slot, source))
                .map_err(|e| Error::Format(e.to_string()))?;
            Ok(crate::recombination::reference::Input {
                role: blake3::hash(&bytes).to_hex().to_string(),
                molecule: molecule.clone(),
            })
        })
        .collect()
}

pub fn verify(plan: &DesignPlan, target: &Molecule) -> VerificationReport {
    evaluate(plan, target).0
}

/// Fresh private material state from the same reconstruction as the report.
/// Only in-process recipe adapters can access it; serialized product metadata or
/// imported reports cannot supply this vector. A failed step leaves only the
/// actually reconstructed leading prefix, aligned with report.products.
pub(crate) fn verify_with_physical(
    plan: &DesignPlan,
    target: &Molecule,
) -> (VerificationReport, Vec<Option<Duplex>>) {
    let (report, physical) = evaluate(plan, target);
    (
        report,
        physical
            .into_iter()
            .map(|state| state.map(|tracked| tracked.duplex))
            .collect(),
    )
}

/// Keep the full freshly reconstructed strand-local state for serial cloning.
/// This is deliberately private to the native core: a saved report, metadata or
/// deserialized tracked state cannot substitute for this reconstruction pass.
pub(crate) fn verify_with_tracked(
    plan: &DesignPlan,
    target: &Molecule,
) -> (VerificationReport, Vec<Option<TrackedDuplex>>) {
    evaluate(plan, target)
}

/// Capture the actual final Golden Gate enumeration from the same native pass.
/// Serialized reports cannot provide this private reconstruction receipt.
pub(crate) fn verify_with_golden_gate_trace(
    plan: &DesignPlan,
    target: &Molecule,
) -> (
    VerificationReport,
    Vec<Option<Duplex>>,
    Option<crate::golden_gate::GoldenGateReport>,
) {
    let mut captured = None;
    let (report, physical) = evaluate_captured(
        plan,
        target,
        &mut captured,
        &mut modified_duplex::CleavageBudget::new(100_000_000)
            .expect("fixed native cleavage bound"),
    );
    (
        report,
        physical.into_iter().map(|s| s.map(|t| t.duplex)).collect(),
        captured,
    )
}

fn evaluate(
    plan: &DesignPlan,
    target: &Molecule,
) -> (VerificationReport, Vec<Option<TrackedDuplex>>) {
    evaluate_captured(
        plan,
        target,
        &mut None,
        &mut modified_duplex::CleavageBudget::new(100_000_000)
            .expect("fixed native cleavage bound"),
    )
}
pub(crate) fn verify_with_tracked_budget(
    plan: &DesignPlan,
    target: &Molecule,
    budget: &mut modified_duplex::CleavageBudget,
) -> (VerificationReport, Vec<Option<TrackedDuplex>>) {
    evaluate_captured(plan, target, &mut None, budget)
}
fn evaluate_captured(
    plan: &DesignPlan,
    target: &Molecule,
    final_golden_gate: &mut Option<crate::golden_gate::GoldenGateReport>,
    cleavage_budget: &mut modified_duplex::CleavageBudget,
) -> (VerificationReport, Vec<Option<TrackedDuplex>>) {
    // New embedded profile snapshots are bounded before hashing their public native DTOs.
    let preflight = validate_reference_steps(plan).and_then(|()| {
        if plan.steps.iter().any(|s| {
            matches!(
                s,
                Operation::RecombineReference { .. }
                    | Operation::MultiSiteRecombineReference { .. }
                    | Operation::Topo { .. }
                    | Operation::CheckedTopo { .. }
                    | Operation::ProcessedTopo { .. }
                    | Operation::TopoActions { .. }
                    | Operation::DeclareSourceTermini { .. }
                    | Operation::CatalogDigest { .. }
                    | Operation::SelectedGoldenGate { .. }
                    | Operation::HomologyAssembly { .. }
                    | Operation::ResectThreePrime { .. }
                    | Operation::HomologyAnnealing { .. }
            )
        }) {
            let mut models = plan.inputs.iter().map(|i| &i.molecule).collect::<Vec<_>>();
            models.push(target);
            crate::recombination::reference::model_preflight(&models)?;
        }
        Ok(())
    });
    if let Err(error) = preflight {
        return (VerificationReport { plan_revision: String::new(), target_revision: String::new(), engine_version: ENGINE_VERSION.into(), status: Status::Invalid,
            findings: vec![error.to_string()], assumptions: vec!["Reference/TOPO request rejected before plan hashing; no current result binding established.".into()], products: vec![], exhaustive_alternative_search: false, step_reports: vec![], annotation_reports: vec![], processed_step_reports: vec![], physical_processing_reports: vec![], reference_step_reports: vec![] }, Vec::new());
    }
    let mut report = VerificationReport {
        plan_revision: plan_revision(plan), target_revision: target.revision(), engine_version: ENGINE_VERSION.into(),
        status: Status::Invalid, findings: vec![], products: vec![], exhaustive_alternative_search: false, step_reports: vec![], annotation_reports: vec![], processed_step_reports: vec![], physical_processing_reports: vec![], reference_step_reports: vec![],
        assumptions: vec!["Exact annealing/overlap model; no experimental outcome is guaranteed.".into(),
            "RecombineReference is a bounded full-site-profile reciprocal reference exchange. All components are computed before an explicit selection; neither its selected DNA nor target agreement establishes chemical strand state, biological repair, efficacy or a selected experimental winner. Reference-only ancestry cannot supply physical PCR or annealing material.".into(),
            "TOPO reconstructs explicitly declared strands, pairing and elementary transfers. Its full physical graph is retained; a selected proven paired contour supplies a restrictive reference view only. No implicit target acceptance policy, nick/flap repair, supplier-kit state or experimental success is asserted. Downstream physical graph consumption is not yet modeled.".into(),
            "ProcessedTopo freshly composes original declared preparations, transfer and selected nick-sealing steps before proving a source-call-selected paired contour. A faithful bridge can feed private downstream physical operations under explicit projection options. Unknown modifications remain unknown; imported reports/sequence summaries cannot restore authority. Other components, unsealed nicks and unapplied branches remain evidence; component selection is not isolation or experimental success.".into(),
            "Exact-PCR operations consider exact non-overlapping binding sites within one template turn, including competing oligo orientations and single-oligo products. Intentional-mismatch PCR uses its separately recorded mismatch and 3′ policy; ambiguous and overlapping-domain products remain incomplete previews.".into(),
            "Golden Gate enumerates a bounded one-copy fragment pool under explicit options; recut paths, carryover and alternative source histories remain evidence, not a claim of exhaustive kinetic mixture uniqueness.".into(),
            "Ordinary digest is an explicit complete-cleavage assumption. Golden Gate and CheckedDigest evaluate actual tracked positions under the current pinned sensitivity profile; reaction conditions, secondary structure and genome-wide specificity are not established.".into(),
            "Ordered overlap/ligation candidates are checked for consistency; alternative mixture products are not exhaustively searched.".into(),
            "Legacy Digest assumes complete cleavage with 5′-phosphate/3′-OH ends. DeclareSourceModifications applies only to original inputs under a closed Dam/Dcm/CpG generation model; it is not a methylation measurement or evidence about other modification systems. Tracked descendants preserve actual strand states; unknown/impaired checked sites withhold fully determined products.".into(),
            "Exact PCR physical intermediates are ideal blunt duplexes with unmodified 5′-OH/3′-OH termini. Polymerase-specific terminal additions and modified primers are not inferred; single-nucleotide additions require an explicit Tail3Prime operation.".into(),
            "End repair performs only explicit 3′ fill-in or removal under its accessible-OH policy; 5′ chemistry and internal nicks are not repaired. Oligo annealing assumes independently declared linear ssDNA inputs and explicitly supplied terminal chemistry, searches exact antiparallel placements, and forms no covalent bonds.".into(),
            "Single-nucleotide 3′ tailing is an explicit ideal blunt/OH endpoint transform, not predicted polymerase activity. Method-specific overlap recommendations are advisory sequence checks, not proprietary reagent or experimental compatibility evidence.".into(),
            "Product annotations are projected from source coordinates; clipped, interrupted and uncertain features require review. Unmapped annotations remain source-graph evidence. Imported metadata cannot assert end chemistry.".into()],
    };
    let plan_hash = report.plan_revision.clone();
    let mut physical = Vec::new();
    match reconstruct(
        plan,
        target,
        &plan_hash,
        &mut report.products,
        &mut report.findings,
        &mut report.step_reports,
        &mut report.annotation_reports,
        &mut report.processed_step_reports,
        &mut report.physical_processing_reports,
        &mut report.reference_step_reports,
        &mut physical,
        final_golden_gate,
        cleavage_budget,
    ) {
        Ok(incomplete) => {
            let product = report.products.last().expect("nonempty reconstructed plan");
            let matches = plan_target_matches(plan, product, target);
            if matches {
                report.status = if incomplete.is_empty() {
                    Status::ConsistentUnderAssumptions
                } else {
                    Status::Incomplete
                };
                report.findings.push(if matches!(plan.steps.last(), Some(Operation::CheckedTopo { .. })) {
                    "Fresh final TOPO graph satisfies the explicitly declared physical and independent-target comparison policy. Coordinate agreement and declared graph assumptions do not establish experimental success or repaired chemistry.".into()
                } else {
                    "Reconstructed product matches the separately bound target sequence and topology. Orientation is preserved.".into()
                });
            } else {
                report.findings.push(
                    "Reconstructed product does not match the target sequence and topology.".into(),
                );
            }
            report.findings.extend(incomplete);
        }
        Err(Error::Unsupported(reason)) => {
            report.status = Status::Incomplete;
            report.findings.push(reason);
        }
        Err(error) => report.findings.push(error.to_string()),
    }
    (report, physical)
}

/// Coordinate comparison only; this helper confers no physical or report authority.
/// Callers bind and bound models separately. Legacy plans retain their oriented
/// comparison; final CheckedTopo uses only its explicitly declared allowances.
pub fn plan_target_matches(plan: &DesignPlan, product: &Molecule, target: &Molecule) -> bool {
    if product.topology != target.topology
        || product.sequence.alphabet() != target.sequence.alphabet()
    {
        return false;
    }
    let policy = match plan.steps.last() {
        Some(Operation::CheckedTopo { policy, .. }) => Some(policy),
        _ => None,
    };
    let same = |candidate: &str| {
        candidate == target.sequence.as_str()
            || target.topology == Topology::Circular
                && policy.is_none_or(|p| p.allow_circular_rotation)
                && analysis::rotation_offset(candidate, target.sequence.as_str()).is_some()
    };
    same(product.sequence.as_str())
        || policy.is_some_and(|p| p.allow_reverse_complement)
            && product
                .sequence
                .reverse_complement()
                .is_ok_and(|rc| same(rc.as_str()))
}

fn region_map(
    length: usize,
    start: usize,
    count: usize,
    destination: usize,
) -> Result<Vec<MappingSegment>, Error> {
    if start >= length || count > length || count == 0 {
        return Err(Error::Coordinates(
            "annotation source interval is outside one source turn".into(),
        ));
    }
    let first = count.min(length - start);
    let mut map = vec![MappingSegment {
        source_start: start,
        source_end: start + first,
        destination_start: destination,
        reverse: false,
    }];
    if first < count {
        map.push(MappingSegment {
            source_start: 0,
            source_end: count - first,
            destination_start: destination + first,
            reverse: false,
        });
    }
    Ok(map)
}
/// Old top base x is on the reversed bottom strand at boundary E-x.
/// Once ligated, neighboring top segments can restore its pairing, so clip
/// against the final reference rather than the original paired intersection.
fn reverse_reference_map(
    source_length: usize,
    bottom_end: i64,
    destination: usize,
    product_length: usize,
    topology: Topology,
) -> Result<Vec<MappingSegment>, Error> {
    let end = destination as i64 + bottom_end;
    let n = product_length as i64;
    if topology == Topology::Linear {
        let start = (end - n).max(0).min(source_length as i64) as usize;
        let stop = end.max(0).min(source_length as i64) as usize;
        return Ok(if start < stop {
            vec![MappingSegment {
                source_start: start,
                source_end: stop,
                destination_start: (end - stop as i64) as usize,
                reverse: true,
            }]
        } else {
            vec![]
        });
    }
    if source_length > product_length {
        return Err(Error::Unsupported(
            "reversed strand exceeds one final circular reference turn".into(),
        ));
    }
    let mut map = Vec::new();
    let mut start = 0;
    while start < source_length {
        let boundary = (end - start as i64).rem_euclid(n) as usize;
        let upper = if boundary == 0 {
            product_length
        } else {
            boundary
        };
        let count = (source_length - start).min(upper);
        map.push(MappingSegment {
            source_start: start,
            source_end: start + count,
            destination_start: upper - count,
            reverse: true,
        });
        start += count;
    }
    Ok(map)
}

fn project_into(
    source: &Molecule,
    map: &[MappingSegment],
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<(), Error> {
    project_into_with_collapse(source, map, product, findings, false)
}
/// Shared source-top coordinate projection for a freshly reconstructed digest fragment.
/// The caller owns cleavage eligibility; this helper creates no physical authority.
pub(crate) fn project_digest_annotations(
    source: &Molecule,
    fragment: &chemistry::DigestFragment,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<(), Error> {
    let map = region_map(
        source.sequence.len(),
        fragment.top_start,
        fragment.duplex.top.len(),
        0,
    )?;
    project_into(source, &map, product, findings)
}
fn project_into_with_collapse(
    source: &Molecule,
    map: &[MappingSegment],
    product: &mut Molecule,
    findings: &mut Vec<String>,
    allow_homology_collapse: bool,
) -> Result<(), Error> {
    if source.features.is_empty() {
        return Ok(());
    }
    let projection = if allow_homology_collapse {
        annotations::project_collapsing_features(
            source,
            map,
            product.sequence.len(),
            product.topology,
        )?
    } else {
        annotations::project_features_with_topology(
            source,
            map,
            product.sequence.len(),
            product.topology,
        )?
    };
    product.features.extend(projection.features);
    if !projection.unmapped_features.is_empty() {
        product.metadata.insert(
            format!("annotation_mapping:unmapped:{}", source.revision()),
            serde_json::to_string(&projection.unmapped_features)
                .map_err(|e| Error::Format(e.to_string()))?,
        );
    }
    for warning in projection.warnings.iter().take(100) {
        findings.push(warning.chars().take(2048).collect());
    }
    if projection.warnings.len() > 100 {
        findings.push("Additional annotation-projection findings omitted; source evidence and review flags remain in the plan.".into());
    }
    // Identical annotations shared by homology arms are represented once, with sources in the graph.
    let mut seen = std::collections::BTreeSet::new();
    product.features.retain(|f| {
        seen.insert(*blake3::hash(&serde_json::to_vec(f).expect("serializable feature")).as_bytes())
    });
    product.validate()?;
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RepairPreview {
    pub source_revision: String,
    pub enzymes: Vec<String>,
    pub fragment_index: usize,
    pub method: chemistry::EndRepairMethod,
    pub molecule: Molecule,
    pub repair: chemistry::EndRepairReport,
    pub findings: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnnealPreview {
    pub top_revision: String,
    pub bottom_revision: String,
    pub annealing: chemistry::AnnealReport,
}
fn bounded_physical_preview(value: &impl Serialize) -> Result<(), Error> {
    struct Count(usize);
    impl std::io::Write for Count {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0 = self
                .0
                .checked_add(bytes.len())
                .filter(|n| *n <= chemistry::MAX_CHEMISTRY_REPORT_BYTES)
                .ok_or_else(|| {
                    std::io::Error::other("physical preview exceeds normal message budget")
                })?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(&mut Count(0), value).map_err(|_| Error::TooLarge)
}
fn carry_unmapped_evidence(source: &Molecule, product: &mut Molecule) {
    for (key, value) in &source.metadata {
        if key.starts_with("annotation_mapping:unmapped") {
            product
                .metadata
                .entry(key.clone())
                .or_insert_with(|| value.clone());
        }
    }
}
pub(crate) fn repaired_document(
    source: &Molecule,
    report: &chemistry::EndRepairReport,
    name: String,
    findings: &mut Vec<String>,
) -> Result<Molecule, Error> {
    let mut product = report.duplex.molecule(name)?;
    let map = report
        .mappings
        .iter()
        .filter(|m| {
            m.source_strand == chemistry::Strand::Top
                && m.product_strand == chemistry::Strand::Top
                && !m.synthesized
        })
        .map(|m| MappingSegment {
            source_start: m.source_start,
            source_end: m.source_end,
            destination_start: m.product_start,
            reverse: m.reverse,
        })
        .collect::<Vec<_>>();
    project_into(source, &map, &mut product, findings)?;
    carry_unmapped_evidence(source, &mut product);
    product
        .metadata
        .insert("end_repair:source_revision".into(), source.revision());
    product
        .metadata
        .insert("end_repair:method".into(), format!("{:?}", report.method));
    product
        .metadata
        .insert("end_repair:engine_version".into(), ENGINE_VERSION.into());
    product.validate()?;
    Ok(product)
}
/// Standalone preview reconstructs the digest rather than accepting an imported
/// duplex summary. The returned molecule is reference/evidence, not physical authority.
pub fn simulate_end_repair(
    source: &Molecule,
    enzymes: &[String],
    fragment_index: usize,
    method: chemistry::EndRepairMethod,
) -> Result<RepairPreview, Error> {
    source.validate()?;
    if enzymes.is_empty() || enzymes.len() > 128 {
        return Err(Error::Format(
            "end-repair preview requires 1–128 explicit restriction enzymes".into(),
        ));
    }
    let catalog = analysis::starter_enzymes();
    let selected = enzymes
        .iter()
        .map(|name| {
            catalog
                .iter()
                .find(|enzyme| enzyme.name.eq_ignore_ascii_case(name))
                .cloned()
                .ok_or_else(|| {
                    Error::Unsupported(format!("enzyme {name} is not in the embedded catalog"))
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let initial = Duplex::from_molecule(source)?;
    let fragments = chemistry::digest(&initial, &selected)?;
    let fragment = fragments.get(fragment_index).ok_or_else(|| {
        Error::Format(format!(
            "digest has {} fragments; selected fragment index {fragment_index} is unavailable",
            fragments.len()
        ))
    })?;
    let mut findings = Vec::new();
    let mut intermediate = fragment.duplex.molecule("Reconstructed digest fragment")?;
    project_into(
        source,
        &region_map(
            source.sequence.len(),
            fragment.top_start,
            fragment.duplex.top.len(),
            0,
        )?,
        &mut intermediate,
        &mut findings,
    )?;
    carry_unmapped_evidence(source, &mut intermediate);
    let repair = chemistry::repair_ends(&fragment.duplex, method)?;
    let mut molecule = repaired_document(
        &intermediate,
        &repair,
        format!("{}_repaired_{}", source.name, fragment_index + 1),
        &mut findings,
    )?;
    molecule.metadata.insert(
        "end_repair:original_source_revision".into(),
        source.revision(),
    );
    molecule.metadata.insert(
        "end_repair:enzymes".into(),
        serde_json::to_string(enzymes).map_err(|e| Error::Format(e.to_string()))?,
    );
    molecule.metadata.insert(
        "end_repair:fragment_index".into(),
        fragment_index.to_string(),
    );
    findings.extend(repair.warnings.iter().cloned());
    crate::protocol::save_document(&molecule)?;
    let preview = RepairPreview {
        source_revision: source.revision(),
        enzymes: enzymes.to_vec(),
        fragment_index,
        method,
        molecule,
        repair,
        findings,
    };
    bounded_physical_preview(&preview)?;
    Ok(preview)
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum TailSubstrate {
    Digest {
        enzymes: Vec<String>,
        fragment_index: usize,
        blunt_before: bool,
    },
    ExactPcr {
        forward: Primer,
        reverse: Primer,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TailPreview {
    pub source_revision: String,
    pub recipe: TailSubstrate,
    pub base: crate::tailing::TailingBase,
    pub molecule: Molecule,
    pub tailing: crate::tailing::TailingReport,
    pub findings: Vec<String>,
}
fn tailed_document(
    source: &Molecule,
    report: &crate::tailing::TailingReport,
    name: String,
    findings: &mut Vec<String>,
) -> Result<Molecule, Error> {
    let mut product = report.duplex.molecule(name)?;
    let maps = report
        .mappings
        .iter()
        .filter(|m| {
            m.source_strand == chemistry::Strand::Top
                && m.product_strand == chemistry::Strand::Top
                && !m.synthesized
        })
        .map(|m| MappingSegment {
            source_start: m.source_start,
            source_end: m.source_end,
            destination_start: m.product_start,
            reverse: m.reverse,
        })
        .collect::<Vec<_>>();
    project_into(source, &maps, &mut product, findings)?;
    carry_unmapped_evidence(source, &mut product);
    product
        .metadata
        .insert("tailing:source_revision".into(), source.revision());
    product
        .metadata
        .insert("tailing:base".into(), report.base.as_str().into());
    product
        .metadata
        .insert("tailing:engine_version".into(), ENGINE_VERSION.into());
    product.validate()?;
    Ok(product)
}
/// Readonly evidence reconstructed from an explicit digest or unique exact PCR.
/// The returned molecule metadata never establishes physical authority in a later plan.
pub fn simulate_tailing(
    source: &Molecule,
    recipe: &TailSubstrate,
    base: crate::tailing::TailingBase,
) -> Result<TailPreview, Error> {
    source.validate()?;
    let mut findings = Vec::new();
    let (intermediate, duplex) = match recipe {
        TailSubstrate::Digest {
            enzymes,
            fragment_index,
            blunt_before,
        } => {
            if enzymes.is_empty() || enzymes.len() > 128 {
                return Err(Error::Format(
                    "tailing preview needs1–128 explicit digest enzymes".into(),
                ));
            }
            let catalog = analysis::starter_enzymes();
            let selected = enzymes
                .iter()
                .map(|name| {
                    catalog
                        .iter()
                        .find(|e| e.name.eq_ignore_ascii_case(name))
                        .cloned()
                        .ok_or_else(|| {
                            Error::Unsupported(format!(
                                "enzyme {name} is absent from the embedded catalog"
                            ))
                        })
                })
                .collect::<Result<Vec<_>, _>>()?;
            let initial = Duplex::from_molecule(source)?;
            let fragments = chemistry::digest(&initial, &selected)?;
            let fragment = fragments
                .get(*fragment_index)
                .ok_or_else(|| Error::Format("selected digest fragment is unavailable".into()))?;
            let mut molecule = fragment.duplex.molecule("Tailing digest substrate")?;
            project_into(
                source,
                &region_map(
                    source.sequence.len(),
                    fragment.top_start,
                    fragment.duplex.top.len(),
                    0,
                )?,
                &mut molecule,
                &mut findings,
            )?;
            carry_unmapped_evidence(source, &mut molecule);
            if *blunt_before {
                let report =
                    chemistry::repair_ends(&fragment.duplex, chemistry::EndRepairMethod::Blunt)?;
                let repaired = repaired_document(
                    &molecule,
                    &report,
                    "Blunted tailing substrate".into(),
                    &mut findings,
                )?;
                (repaired, report.duplex)
            } else {
                (molecule, fragment.duplex.clone())
            }
        }
        TailSubstrate::ExactPcr { forward, reverse } => {
            if chemistry::has_unmapped_physical_evidence(source) {
                return Err(Error::Unsupported("PCR tailing source has imported physical evidence; reconstruct its originating graph".into()));
            }
            let matches = analysis::simulate_pcr_geometry(source, forward, reverse)?;
            if matches.len() != 1 {
                return Err(Error::Unsupported(format!(
                    "PCR tailing requires exactly one exact product; found {}",
                    matches.len()
                )));
            }
            let geometry = &matches[0];
            let duplex = Duplex::exact_pcr(&geometry.product.sequence)?;
            let mut product = duplex.molecule("Exact PCR tailing substrate")?;
            project_pcr_annotations(source, geometry, &mut product, &mut findings)?;
            (product, duplex)
        }
    };
    let tailing = crate::tailing::tail_blunt_duplex(&duplex, base)?;
    let mut molecule = tailed_document(
        &intermediate,
        &tailing,
        format!("{}_{}-tailed", source.name, base.as_str()),
        &mut findings,
    )?;
    molecule
        .metadata
        .insert("tailing:original_source_revision".into(), source.revision());
    molecule.metadata.insert(
        "tailing:recipe".into(),
        serde_json::to_string(recipe).map_err(|e| Error::Format(e.to_string()))?,
    );
    findings.extend(tailing.warnings.iter().cloned());
    let preview = TailPreview {
        source_revision: source.revision(),
        recipe: recipe.clone(),
        base,
        molecule,
        tailing,
        findings,
    };
    bounded_physical_preview(&preview)?;
    Ok(preview)
}
pub(crate) fn declared_oligo(source: &Molecule) -> Result<(), Error> {
    source.validate()?;
    if source.sequence.alphabet() != Alphabet::Dna || source.topology != Topology::Linear {
        return Err(Error::Unsupported(
            "annealing requires declared linear DNA oligo inputs".into(),
        ));
    }
    if source.metadata.contains_key("duplex_model")
        || source
            .metadata
            .get("source:legacy:chemistry_unmapped")
            .is_some_and(|value| value == "true")
        || source
            .metadata
            .get("source:snapgene:chemistry_unmapped")
            .is_some_and(|value| value == "true")
    {
        return Err(Error::Unsupported("imported physical intermediates cannot be flattened into ssDNA oligos; use independently declared oligo inputs".into()));
    }
    if source
        .metadata
        .get("genbank.molecule_type")
        .is_some_and(|value| {
            let value = value.to_ascii_uppercase();
            value.starts_with("DS-") || value.starts_with("MS-")
        })
        || source
            .original_genbank
            .as_ref()
            .and_then(|text| text.lines().find(|line| line.starts_with("LOCUS")))
            .is_some_and(|line| {
                line.split_ascii_whitespace().any(|word| {
                    let word = word.to_ascii_uppercase();
                    word.starts_with("DS-") || word.starts_with("MS-")
                })
            })
        || source.metadata.contains_key("source:snapgene:dna_flags")
    {
        return Err(Error::Unsupported(
            "an explicitly double/mixed-stranded input is not a declared ssDNA oligo".into(),
        ));
    }
    Ok(())
}
pub fn simulate_oligo_annealing(
    top: &Molecule,
    bottom: &Molecule,
    options: &chemistry::AnnealOptions,
) -> Result<AnnealPreview, Error> {
    declared_oligo(top)?;
    declared_oligo(bottom)?;
    let annealing = chemistry::anneal_oligos(&top.sequence, &bottom.sequence, options)?;
    let preview = AnnealPreview {
        top_revision: top.revision(),
        bottom_revision: bottom.revision(),
        annealing,
    };
    bounded_physical_preview(&preview)?;
    Ok(preview)
}
fn annealed_document_inner(
    top: &Molecule,
    bottom: &Molecule,
    report: &chemistry::AnnealReport,
    candidate_index: usize,
    findings: &mut Vec<String>,
) -> Result<Molecule, Error> {
    let candidate = report
        .candidates
        .get(candidate_index)
        .ok_or_else(|| Error::Format("selected annealing candidate is unavailable".into()))?;
    let mut product = candidate.duplex.molecule("Annealed oligos")?;
    project_into(
        top,
        &[MappingSegment {
            source_start: 0,
            source_end: top.sequence.len(),
            destination_start: 0,
            reverse: false,
        }],
        &mut product,
        findings,
    )?;
    let end = candidate.duplex.bottom_start + bottom.sequence.len() as i64;
    project_into(
        bottom,
        &reverse_reference_map(
            bottom.sequence.len(),
            end,
            0,
            product.sequence.len(),
            Topology::Linear,
        )?,
        &mut product,
        findings,
    )?;
    carry_unmapped_evidence(top, &mut product);
    carry_unmapped_evidence(bottom, &mut product);
    product
        .metadata
        .insert("anneal:top_revision".into(), top.revision());
    product
        .metadata
        .insert("anneal:bottom_revision".into(), bottom.revision());
    product
        .metadata
        .insert("anneal:candidate_id".into(), candidate.id.clone());
    product
        .metadata
        .insert("anneal:engine_version".into(), ENGINE_VERSION.into());
    product.validate()?;
    crate::protocol::save_document(&product)?;
    Ok(product)
}
/// Reconstruct standalone read-only preview evidence; no imported report can
/// provide physical authority or edited annotation coordinates.
pub fn annealed_document(
    top: &Molecule,
    bottom: &Molecule,
    preview: &AnnealPreview,
    candidate_index: usize,
) -> Result<Molecule, Error> {
    if preview.top_revision != top.revision() || preview.bottom_revision != bottom.revision() {
        return Err(Error::Format(
            "annealing preview source revision changed".into(),
        ));
    }
    let rebuilt = simulate_oligo_annealing(top, bottom, &preview.annealing.options)?;
    bounded_physical_preview(preview)?;
    if rebuilt != *preview {
        return Err(Error::Format(
            "annealing preview differs from independently reconstructed evidence".into(),
        ));
    }
    annealed_document_inner(
        top,
        bottom,
        &rebuilt.annealing,
        candidate_index,
        &mut Vec::new(),
    )
}

fn assemble_overlap(
    fragments: &[&Molecule],
    overlaps: &[usize],
    circular: bool,
    index: usize,
    findings: &mut Vec<String>,
) -> Result<Molecule, Error> {
    if fragments.is_empty()
        || fragments.len() > 1000
        || overlaps.len() != fragments.len() - 1 + usize::from(circular)
    {
        return Err(Error::Format(
            "ordered assembly requires one overlap per junction".into(),
        ));
    }
    for f in fragments {
        if f.topology != Topology::Linear || f.sequence.alphabet() != Alphabet::Dna {
            return Err(Error::Unsupported(
                "overlap inputs must be linear DNA fragments".into(),
            ));
        }
        if f.sequence
            .as_str()
            .bytes()
            .any(|b| analysis::base_mask(b).count_ones() != 1)
        {
            return Err(Error::Unsupported("ambiguous assembly sequence".into()));
        }
    }
    for (i, fragment) in fragments.iter().enumerate() {
        let incoming = if i > 0 {
            overlaps[i - 1]
        } else if circular {
            *overlaps.last().unwrap()
        } else {
            0
        };
        let outgoing = if i + 1 < fragments.len() {
            overlaps[i]
        } else if circular {
            *overlaps.last().unwrap()
        } else {
            0
        };
        if incoming.saturating_add(outgoing) >= fragment.sequence.len() {
            return Err(Error::Unsupported(format!(
                "fragment {} has overlapping or exhaustive overlap arms; a more detailed assembly model is required",
                i + 1
            )));
        }
    }
    let mut sequence = fragments[0].sequence.as_str().to_owned();
    for (i, fragment) in fragments.iter().enumerate().skip(1) {
        let overlap = overlaps[i - 1];
        let next = fragment.sequence.as_str();
        let previous = fragments[i - 1].sequence.as_str();
        if overlap == 0
            || overlap >= previous.len()
            || overlap >= next.len()
            || previous[previous.len() - overlap..] != next[..overlap]
        {
            return Err(Error::Format(format!(
                "step {}: invalid overlap at junction {i}",
                index + 1
            )));
        }
        if sequence.len() + next.len() - overlap > MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        sequence.push_str(&next[overlap..]);
    }
    if circular {
        let overlap = overlaps[overlaps.len() - 1];
        if overlap == 0
            || overlap >= fragments[0].sequence.len()
            || overlap >= fragments.last().unwrap().sequence.len()
            || overlap >= sequence.len()
            || sequence[sequence.len() - overlap..] != sequence[..overlap]
        {
            return Err(Error::Format(
                "circular closure overlap does not match".into(),
            ));
        }
        sequence.truncate(sequence.len() - overlap);
    }
    let mut product = Molecule::dna(
        format!("Assembly {}", index + 1),
        &sequence,
        if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )?;
    let mut offset = 0usize;
    for (i, fragment) in fragments.iter().enumerate() {
        let length = fragment.sequence.len();
        let n = product.sequence.len();
        let destination = offset % n;
        let first = length.min(n - destination);
        let mut map = vec![MappingSegment {
            source_start: 0,
            source_end: first,
            destination_start: destination,
            reverse: false,
        }];
        if first < length {
            map.push(MappingSegment {
                source_start: first,
                source_end: length,
                destination_start: 0,
                reverse: false,
            });
        }
        project_into_with_collapse(
            fragment,
            &map,
            &mut product,
            findings,
            circular && fragments.len() == 1,
        )?;
        offset += length;
        if i < overlaps.len() {
            offset -= overlaps[i];
        }
    }
    Ok(product)
}

fn reconstruct(
    plan: &DesignPlan,
    target: &Molecule,
    plan_hash: &str,
    products: &mut Vec<Molecule>,
    findings: &mut Vec<String>,
    step_reports: &mut Vec<StepReport>,
    annotation_reports: &mut Vec<AnnotationStepReport>,
    processed_step_reports: &mut Vec<ProcessedTopoStepReport>,
    physical_processing_reports: &mut Vec<PhysicalProcessingStepReport>,
    reference_step_reports: &mut Vec<ReferenceStepReport>,
    physical: &mut Vec<Option<TrackedDuplex>>,
    final_golden_gate: &mut Option<crate::golden_gate::GoldenGateReport>,
    cleavage_budget: &mut modified_duplex::CleavageBudget,
) -> Result<Vec<String>, Error> {
    target.validate()?;
    if plan.target_revision != target.revision() {
        return Err(Error::Format(
            "target revision changed or was not bound separately".into(),
        ));
    }
    if plan.steps.is_empty() || plan.steps.len() > 1000 || plan.inputs.len() > 1000 {
        return Err(Error::Format(
            "plan needs 1–1000 steps and at most 1000 inputs".into(),
        ));
    }
    let mut inputs = BTreeMap::new();
    for input in &plan.inputs {
        input.molecule.validate()?;
        if input.id.trim().is_empty() || inputs.insert(input.id.clone(), &input.molecule).is_some()
        {
            return Err(Error::Format("empty or duplicate plan input ID".into()));
        }
        if input.revision != input.molecule.revision() {
            return Err(Error::Format(format!("stale input revision: {}", input.id)));
        }
    }
    let mut reference_only: Vec<bool> = vec![];
    let mut product_bases = 0usize;
    let mut product_bytes = 0usize;
    let mut step_evidence_bytes = 0usize;
    let mut incomplete = Vec::new();
    for (index, step) in plan.steps.iter().enumerate() {
        if matches!(step, Operation::CheckedTopo { .. }) && index + 1 != plan.steps.len() {
            return Err(Error::Unsupported("CheckedTopo independent-target acceptance is supported only as the final plan step; intermediate target acceptance and downstream physical graph consumption are not modeled.".into()));
        }
        let inherited_reference_only = operation_sources(step).iter().any(|source| match source {
            Source::Step(i) => reference_only.get(*i).copied().unwrap_or(false),
            Source::Input(_) => false,
        });
        if inherited_reference_only
            && matches!(
                step,
                Operation::Pcr { .. }
                    | Operation::PcrWithMismatches { .. }
                    | Operation::AnnealOligos { .. }
                    | Operation::Topo { .. }
                    | Operation::CheckedTopo { .. }
                    | Operation::ProcessedTopo { .. }
                    | Operation::TopoActions { .. }
            )
        {
            return Err(Error::Unsupported("Reference-only ancestry cannot establish a physical PCR/annealing/TOPO substrate. Its full originating graph must be consumed by an explicit supported physical operation; no materialization or graph repair is inferred.".into()));
        }
        if inherited_reference_only
            && matches!(
                step,
                Operation::HomologyAssembly { .. }
                    | Operation::HomologyAnnealing { .. }
                    | Operation::ResectThreePrime { .. }
            )
        {
            return Err(Error::Unsupported("Reference-only ancestry cannot provide a physical homology substrate; reconstruct actual original strands first.".into()));
        }
        let resolve = |source: &Source| -> Result<&Molecule, Error> {
            match source {
                Source::Input(id) => inputs
                    .get(id)
                    .copied()
                    .ok_or_else(|| Error::Format(format!("unknown input {id}"))),
                Source::Step(i) => products.get(*i).ok_or_else(|| {
                    Error::Format(format!(
                        "step {index} has a cyclic or forward dependency on step {i}"
                    ))
                }),
            }
        };
        let resolve_tracked = |source: &Source| -> Result<TrackedDuplex, Error> {
            let m = resolve(source)?;
            match source {
                Source::Input(_) => TrackedDuplex::from_unknown_duplex(&Duplex::from_molecule(m)?),
                Source::Step(i) => physical.get(*i).and_then(Option::as_ref).cloned().ok_or_else(|| Error::Unsupported(format!("step {i} has no established physical strand/modification state; reconstruct a supported physical graph"))),
            }
        };
        let resolve_physical =
            |source: &Source| -> Result<Duplex, Error> { Ok(resolve_tracked(source)?.duplex) };
        let mut next_physical = None;
        let mut diagnostic_reference_only = false;
        let mut product = match step {
            Operation::ResectThreePrime { source, cuts } => {
                if physical_processing_reports.len() >= 64 {
                    return Err(Error::TooLarge);
                }
                let source_model = resolve(source)?;
                let tracked = resolve_tracked(source)?;
                let (processed, geometry) = tracked.resect_three_prime(cuts)?;
                geometry.validate_documentary()?;
                let mut product = processed
                    .duplex
                    .molecule(format!("3′-resected operand {}", index + 1))?;
                let top_end = product.sequence.len();
                if source_model.sequence.as_str().get(..top_end) != Some(product.sequence.as_str())
                {
                    return Err(Error::Format(
                        "Resected top reference does not match the actual original source domain"
                            .into(),
                    ));
                }
                project_into(
                    source_model,
                    &[MappingSegment {
                        source_start: 0,
                        source_end: top_end,
                        destination_start: 0,
                        reverse: false,
                    }],
                    &mut product,
                    findings,
                )?;
                carry_unmapped_evidence(source_model, &mut product);
                findings.extend(geometry.warnings.iter().take(100).cloned());
                findings.push(format!("Step {} removes bottom3′ suffix {} nt and top3′ suffix {} nt. Original5′ end groups and every surviving strand/modification remain; any bottom5′ extension beyond the displayed top reference stays in the actual duplex. No filling, nick sealing or host repair is inferred.",index+1,cuts.left,cuts.right));
                step_evidence_bytes = step_evidence_bytes
                    .checked_add(evidence_bytes(&geometry, 16 * 1024 * 1024)?)
                    .filter(|n| *n <= 32 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
                physical_processing_reports.push(PhysicalProcessingStepReport {
                    step_index: index,
                    evidence: PhysicalProcessingEvidence::ThreePrimeResection(geometry),
                });
                next_physical = Some(processed);
                product
            }
            Operation::HomologyAnnealing {
                fragments,
                options,
                policy,
            } => {
                if physical_processing_reports.len() >= 64 {
                    return Err(Error::TooLarge);
                }
                policy.validate_shape(fragments.len())?;
                options.validate()?;
                if fragments.len() == 1 && policy.topology != Topology::Circular {
                    return Err(Error::Unsupported("Single linear3′ processing uses ResectThreePrime, not a joining assembly operation".into()));
                }
                let mut tracked = Vec::with_capacity(fragments.len());
                let mut oriented = Vec::with_capacity(fragments.len());
                for fragment in fragments {
                    let material = resolve_tracked(&fragment.source)?;
                    tracked.push(if fragment.reverse {
                        material.reverse()?
                    } else {
                        material
                    });
                    let source = resolve(&fragment.source)?;
                    if fragment.reverse {
                        let mut model = Molecule::dna(
                            format!("Oriented {}", source.name),
                            source.sequence.reverse_complement()?.as_str(),
                            source.topology,
                        )?;
                        project_into(
                            source,
                            &[MappingSegment {
                                source_start: 0,
                                source_end: source.sequence.len(),
                                destination_start: 0,
                                reverse: true,
                            }],
                            &mut model,
                            findings,
                        )?;
                        carry_unmapped_evidence(source, &mut model);
                        oriented.push(model);
                    } else {
                        oriented.push(source.clone());
                    }
                }
                let (_separate_processed_operands, geometry) =
                    modified_duplex::homology_anneal(&tracked.iter().collect::<Vec<_>>(), policy)?;
                geometry.validate_documentary()?;
                let mut product = assemble_overlap(
                    &oriented.iter().collect::<Vec<_>>(),
                    &policy.overlaps,
                    policy.topology == Topology::Circular,
                    index,
                    findings,
                )?;
                if product.sequence != geometry.reference || product.topology != policy.topology {
                    return Err(Error::Format("Annotated source-derived contour differs from independently annealed strand reference".into()));
                }
                product.name = format!("Annealed reference {}", index + 1);
                product.metadata.insert(
                    "homology_annealing:policy".into(),
                    serde_json::to_string(policy).map_err(|e| Error::Format(e.to_string()))?,
                );
                product.metadata.insert(
                    "homology_annealing:profile".into(),
                    serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?,
                );
                findings.extend(
                    geometry
                        .warnings
                        .iter()
                        .chain(&geometry.limits)
                        .take(100)
                        .cloned(),
                );
                incomplete.push(format!("Step {} retains {} separate source strands and {} unsealed joining discontinuities. This is a diagnostic sequence contour only: no synthesis, ligation or cellular rescue establishes a reusable assembled product.",index+1,geometry.segments.len(),2*geometry.junctions.len()));
                if !geometry.terminal_search_complete {
                    incomplete.push(format!(
                        "Step {} bounded exposed-terminal alternative search is incomplete.",
                        index + 1
                    ));
                }
                step_evidence_bytes = step_evidence_bytes
                    .checked_add(evidence_bytes(&geometry, 16 * 1024 * 1024)?)
                    .filter(|n| *n <= 32 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
                physical_processing_reports.push(PhysicalProcessingStepReport {
                    step_index: index,
                    evidence: PhysicalProcessingEvidence::HomologyAnnealing(geometry),
                });
                diagnostic_reference_only = true;
                product
            }
            Operation::HomologyAssembly {
                fragments,
                options,
                policy,
            } => {
                if physical_processing_reports.len() >= 64 {
                    return Err(Error::TooLarge);
                }
                policy.validate_shape(fragments.len())?;
                options.validate()?;
                // Resolve actual native states first. Reference-only steps and
                // imported producer metadata cannot supply this operation.
                let mut tracked = Vec::with_capacity(fragments.len());
                for fragment in fragments {
                    let material = resolve_tracked(&fragment.source)?;
                    tracked.push(if fragment.reverse {
                        material.reverse()?
                    } else {
                        material
                    });
                }
                let refs = tracked.iter().collect::<Vec<_>>();
                let (material, geometry) = modified_duplex::homology_assemble(&refs, policy)?;
                step_evidence_bytes = step_evidence_bytes
                    .checked_add(evidence_bytes(&geometry, 32 * 1024 * 1024)?)
                    .filter(|n| *n <= 32 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
                findings.extend(geometry.warnings.iter().take(100).cloned());
                findings.extend(geometry.limits.iter().take(100).cloned());
                if geometry.incomplete {
                    incomplete.push(format!("Step {} physical homology assembly retains unresolved nick chemistry, terminal alternatives or bounded-search limits.",index+1));
                }
                let reference_top_start = geometry.reference_top_start;
                let reference_top_end = geometry.reference_top_end;
                let reference_length = geometry.reference_length;
                let diagnostic = if material.is_none() {
                    geometry
                        .product
                        .as_ref()
                        .map(|d| d.molecule(format!("Homology diagnostic {}", index + 1)))
                        .transpose()?
                } else {
                    None
                };
                physical_processing_reports.push(PhysicalProcessingStepReport {
                    step_index: index,
                    evidence: PhysicalProcessingEvidence::HomologyAssembly(geometry),
                });
                let mut product = if let Some(material) = &material {
                    material
                        .duplex
                        .molecule(format!("Homology assembly {}", index + 1))?
                } else {
                    diagnostic_reference_only = true;
                    incomplete.push(format!("Step {} retains paired diagnostic strands but withholds reusable material: internal nick chemistry cannot be faithfully represented by the current private duplex model.",index+1));
                    diagnostic.ok_or_else(|| Error::Unsupported("Physical homology assembly withheld its product: no unsupported gap fill or reference-derived material is invented.".into()))?
                };
                // Advisory sequence guidance is checked on the exact source
                // molecules, after actual strand reconstruction. It cannot
                // choose a physical path or create the private material above.
                let sources = fragments
                    .iter()
                    .map(|f| resolve(&f.source).cloned())
                    .collect::<Result<Vec<_>, _>>()?;
                let reverse = fragments.iter().map(|f| f.reverse).collect::<Vec<_>>();
                let circular = policy.topology == Topology::Circular;
                let profile = crate::overlap_profiles::check_oriented_fragments(
                    &sources,
                    &reverse,
                    &policy.overlaps,
                    circular,
                    options,
                )?;
                if !profile.geometry_valid {
                    return Err(Error::Format("Source homology guidance disagrees with freshly reconstructed material geometry".into()));
                }
                if !profile.complete {
                    incomplete.push(format!("Step {} advisory homology repeat scan reached its bound; alternatives remain unsearched.",index+1));
                }
                findings.extend(profile.findings.iter().take(100).cloned());
                findings.extend(profile.warnings.iter().take(100).cloned());
                let mut oriented = Vec::with_capacity(sources.len());
                for (source, reverse) in sources.iter().zip(&reverse) {
                    if *reverse {
                        let mut model = Molecule::dna(
                            format!("Oriented {}", source.name),
                            source.sequence.reverse_complement()?.as_str(),
                            source.topology,
                        )?;
                        project_into(
                            source,
                            &[MappingSegment {
                                source_start: 0,
                                source_end: source.sequence.len(),
                                destination_start: 0,
                                reverse: true,
                            }],
                            &mut model,
                            findings,
                        )?;
                        carry_unmapped_evidence(source, &mut model);
                        oriented.push(model);
                    } else {
                        oriented.push(source.clone());
                    }
                }
                let refs = oriented.iter().collect::<Vec<_>>();
                let reference =
                    assemble_overlap(&refs, &policy.overlaps, circular, index, findings)?;
                if reference.sequence.len() != reference_length
                    || reference_top_start >= reference_top_end
                    || reference
                        .sequence
                        .as_str()
                        .get(reference_top_start..reference_top_end)
                        != Some(product.sequence.as_str())
                {
                    return Err(Error::Format("Source-only annotation reference differs from independently constructed physical strands".into()));
                }
                project_into(
                    &reference,
                    &[MappingSegment {
                        source_start: reference_top_start,
                        source_end: reference_top_end,
                        destination_start: 0,
                        reverse: false,
                    }],
                    &mut product,
                    findings,
                )?;
                carry_unmapped_evidence(&reference, &mut product);
                product.metadata.insert(
                    "homology_assembly:policy".into(),
                    serde_json::to_string(policy).map_err(|e| Error::Format(e.to_string()))?,
                );
                product.metadata.insert(
                    "homology_assembly:profile".into(),
                    serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?,
                );
                findings.push(format!("Step {} reconstructs retained source strands, declared 5′ cuts, template-supported synthesis and eligible nick seals before target comparison. This explicit ideal model is not experimental success, host rescue or proprietary kit chemistry.",index+1));
                next_physical = material;
                product
            }
            Operation::OverlapExtension {
                left,
                right,
                options,
            } => {
                if physical_processing_reports.len() >= 64 {
                    return Err(Error::TooLarge);
                }
                let mut a = resolve_tracked(&left.source)?;
                let mut b = resolve_tracked(&right.source)?;
                if left.reverse {
                    a = a.reverse()?;
                }
                if right.reverse {
                    b = b.reverse()?;
                }
                let (extended, geometry) = a.overlap_extend(&b, options)?;
                let mut product = extended
                    .duplex
                    .molecule(format!("Overlap extension {}", index + 1))?;
                let offset = a.duplex.top.len() - options.overlap;
                for (oriented, destination_start) in [(left, 0), (right, offset)] {
                    let source = resolve(&oriented.source)?;
                    project_into(
                        source,
                        &[MappingSegment {
                            source_start: 0,
                            source_end: source.sequence.len(),
                            destination_start,
                            reverse: oriented.reverse,
                        }],
                        &mut product,
                        findings,
                    )?;
                }
                findings.extend(geometry.warnings.iter().take(100).cloned());
                findings.push(format!("Step {} denatures two intact duplexes, pairs the left top strand with the right bottom strand over {} bp and extends their paired 3′ ends. Retained strands preserve modifications; synthesized complements are assumed unmodified. Selection does not prove isolation, amplification efficiency or reaction-mixture uniqueness.", index+1, options.overlap));
                if !geometry.terminal_pairings.complete || !geometry.terminal_pairings.unique {
                    incomplete.push(format!("Step {} selected terminal pairing is not uniquely established within the bounded terminal-pair model; alternative products remain possible.", index+1));
                }
                step_evidence_bytes = step_evidence_bytes
                    .checked_add(evidence_bytes(&geometry, 32 * 1024 * 1024)?)
                    .filter(|n| *n <= 32 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
                physical_processing_reports.push(PhysicalProcessingStepReport {
                    step_index: index,
                    evidence: PhysicalProcessingEvidence::OverlapExtension(geometry),
                });
                next_physical = Some(extended);
                product
            }
            Operation::ReconstructMaterial { template, origin } => {
                if !matches!(template, Source::Input(_)) {
                    return Err(Error::Unsupported(
                        "Derived material reconstruction requires an exact original input event"
                            .into(),
                    ));
                }
                let fresh = crate::material_derivation::replay(origin, resolve(template)?)?;
                let modification_summary = modified_duplex::summary(fresh.tracked())?;
                let duplex_summary = fresh.tracked().duplex.summary()?;
                let bytes = evidence_bytes(origin, crate::material_derivation::MAX_CAPSULE_BYTES)?
                    .checked_add(evidence_bytes(&modification_summary, 32 * 1024 * 1024)?)
                    .ok_or(Error::TooLarge)?;
                step_evidence_bytes = step_evidence_bytes
                    .checked_add(bytes)
                    .filter(|n| *n <= 32 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
                let (
                    reference,
                    tracked,
                    origin_incomplete,
                    origin_findings,
                    origin_findings_total,
                    origin_findings_truncated,
                ) = fresh.into_parts();
                let evidence = DerivedOriginEvidence {
                    origin: origin.clone(),
                    material: MaterialStateEvidence {
                        material_revision: modification_summary.state_revision.clone(),
                        duplex: duplex_summary,
                        modifications: modification_summary,
                        origin_incomplete,
                        origin_findings,
                        origin_findings_total,
                        origin_findings_truncated,
                    },
                };
                evidence.validate_documentary()?;
                findings.push(format!("Step {} freshly reconstructs saved material derivation event {} from its exact original source and typed recipe. Actual strand termini and positional modifications are retained; this is material retrieval, not a new target check or proof of experimental isolation.",index+1,origin.product.event_id));
                if origin_incomplete {
                    incomplete.push(format!("Step {} retains incomplete origin findings; saved material reuse does not upgrade them.",index+1));
                }
                physical_processing_reports.push(PhysicalProcessingStepReport {
                    step_index: index,
                    evidence: PhysicalProcessingEvidence::DerivedOrigin(evidence),
                });
                next_physical = Some(tracked);
                reference
            }
            Operation::ReconstructRunProduct { template, origin } => {
                if !matches!(template, Source::Input(_)) {
                    return Err(Error::Unsupported("Saved material reconstruction requires an explicit saved input event, not an enclosing step projection".into()));
                }
                let fresh = crate::run_material_capsule::replay(origin, resolve(template)?)?;
                let bytes = evidence_bytes(
                    &fresh.evidence,
                    crate::run_material_capsule::MAX_REPORT_BYTES,
                )?;
                step_evidence_bytes = step_evidence_bytes
                    .checked_add(bytes)
                    .filter(|n| *n <= 32 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
                findings.push(format!("Step {} freshly reconstructed saved run {} product {} against its original independent target; strands, nicks, termini, positional modifications and generation provenance remain actual tracked state. The selected event is a declared source, not proof of isolation or experimental success.", index+1, origin.selection.run_id, origin.selection.step_index+1));
                if fresh.incomplete {
                    incomplete.push(format!("Step {} originating run retains incomplete findings; material reuse does not upgrade that assessment. Complete originating evidence is retained.", index+1));
                }
                physical_processing_reports.push(PhysicalProcessingStepReport {
                    step_index: index,
                    evidence: PhysicalProcessingEvidence::RunProduct(fresh.evidence),
                });
                next_physical = Some(fresh.tracked);
                fresh.reference
            }
            Operation::MultiSiteRecombineReference {
                sources,
                request,
                selection,
            } => {
                let resolved = sources
                    .iter()
                    .map(|source| resolve(&source.source))
                    .collect::<Result<Vec<_>, _>>()?;
                crate::multisite_recombination::source_preflight(&resolved, request)?;
                let inputs = sources
                    .iter()
                    .zip(&resolved)
                    .map(
                        |(source, molecule)| crate::recombination::reference::Input {
                            role: source.role.clone(),
                            molecule: (*molecule).clone(),
                        },
                    )
                    .collect::<Vec<_>>();
                let fresh = crate::multisite_recombination::simulate(&inputs, request)?;
                let bytes = evidence_bytes(&fresh, 32 * 1024 * 1024)?;
                step_evidence_bytes = step_evidence_bytes
                    .checked_add(bytes)
                    .filter(|n| *n <= 32 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
                let selected = fresh.selected(selection);
                findings.push(format!("Step {} retains a whole {}-parent multi-site reference exchange: {} candidate networks and every reciprocal component. No target selects pairings; reference provenance is not physical strand/modification authority.",index+1,inputs.len(),fresh.reference.candidates.len()));
                reference_step_reports.push(ReferenceStepReport {
                    step_index: index,
                    evidence: ReferenceStepEvidence::MultiSiteRecombination {
                        selection: selection.clone(),
                        batch: fresh,
                    },
                });
                let selected = selected?;
                incomplete.push(format!("Step {} is a declared multi-site reference design. Recombination chemistry, source-strand inheritance, reaction efficacy and mixture selection remain unestablished.",index+1));
                selected.molecule
            }
            Operation::TopoActions { sources, request } => {
                let molecules = sources
                    .iter()
                    .map(|s| resolve(&s.source))
                    .collect::<Result<Vec<_>, _>>()?;
                validate_topo_shape(sources, &request.preparation)?;
                collagene_model::bounds::model_preflight(&molecules)?;
                molecules.iter().try_fold(0usize, |sum, molecule| {
                    sum.checked_add(molecule.sequence.len())
                        .filter(|n| *n <= 100_000)
                        .ok_or(Error::TooLarge)
                })?;
                let mut templates = Vec::new();
                for (source, molecule) in sources.iter().zip(&molecules) {
                    match &source.source {
                        Source::Step(_) => {
                            // Only actual fresh private material can enter the
                            // PCR-only synthesis route. No imported summary or
                            // matching sequence establishes this permission.
                            templates.push((source.role.clone(), resolve_tracked(&source.source)?));
                        }
                        Source::Input(_) => {
                            if chemistry::has_unmapped_physical_evidence(molecule)
                                || molecule.metadata.contains_key("tracked_modifications")
                                || molecule
                                    .metadata
                                    .contains_key(crate::recombination::REFERENCE_ONLY_KEY)
                                || molecule.metadata.keys().any(|key| key.starts_with("topo:"))
                            {
                                return Err(Error::Unsupported("TopoActions requires original reference sources with explicit preparation declarations; imported physical/reference-only summaries cannot restore their lost graph or tracked modification state.".into()));
                            }
                        }
                    }
                }
                // Native input cloning remains bounded before preparation. The
                // private adapter checks each Step role is PCR-only and agrees
                // exactly with its tracked template; source metadata is retained.
                let inputs = sources
                    .iter()
                    .zip(&molecules)
                    .map(|(source, molecule)| crate::topo::Input {
                        role: source.role.clone(),
                        molecule: (*molecule).clone(),
                    })
                    .collect::<Vec<_>>();
                let tracked_templates = templates
                    .iter()
                    .map(|(role, material)| (role.clone(), material))
                    .collect::<Vec<_>>();
                let fresh = crate::topo_actions::inspect_reconstructed(
                    &inputs,
                    request,
                    &tracked_templates,
                )?;
                if !templates.is_empty() {
                    findings.push(format!("Step {} uses {} fresh upstream material role(s) only as templates for explicitly declared new PCR synthesis. Template atoms, internal modifications and generation history are not inherited into the new strands; arbitrary Prepared correspondence from prior products remains unsupported.", index + 1, templates.len()));
                }
                let fresh_bytes = evidence_bytes(&fresh, crate::topo_actions::MAX_REPORT_BYTES)?;
                step_evidence_bytes = step_evidence_bytes
                    .checked_add(fresh_bytes)
                    .filter(|n| *n <= 32 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
                let selected_bases = fresh
                    .product
                    .as_ref()
                    .map_or(0, |m| m.sequence.len())
                    .checked_add(
                        fresh
                            .tracked
                            .as_ref()
                            .map_or(0, |t| t.top.len() + t.bottom.len()),
                    )
                    .ok_or(Error::TooLarge)?;
                let tracked_bytes = fresh
                    .tracked
                    .as_ref()
                    .map(|t| evidence_bytes(t, 32 * 1024 * 1024))
                    .transpose()?
                    .unwrap_or(0);
                let selected_document_bytes = fresh
                    .product
                    .as_ref()
                    .map(|m| evidence_bytes(m, 32 * 1024 * 1024))
                    .transpose()?
                    .unwrap_or(0);
                if product_bases
                    .checked_add(selected_bases)
                    .is_none_or(|n| n > 32_000_000)
                    || product_bytes
                        .checked_add(step_evidence_bytes)
                        .and_then(|n| n.checked_add(tracked_bytes))
                        .and_then(|n| n.checked_add(selected_document_bytes))
                        .is_none_or(|n| n > 64 * 1024 * 1024)
                {
                    return Err(Error::Unsupported("TopoActions selected material plus full retained graph exceeds the plan material/evidence budget before duplicating its selected product.".into()));
                }
                let graph = &fresh.projection.processing.product;
                findings.push(format!("Step {} fresh TopoActions retains {} components, {} original source calls, {} remaining graph nicks and {} unpaired calls. Explicit component selection is not physical isolation or a mixture winner; all other components and released activations remain typed evidence.",index+1,graph.components.len(),fresh.projection.processing.source_bases,graph.nicks.len(),graph.unpaired.iter().map(|r|r.end-r.start).sum::<usize>()));
                if !fresh
                    .projection
                    .processing
                    .parent
                    .all_requested_transfers_applied
                    || !fresh.projection.processing.all_requested_actions_applied
                {
                    findings.push(format!("Step {} retains an unapplied/unknown requested branch under its explicit projection options; target agreement cannot establish that operation occurred.",index+1));
                }
                findings.extend(fresh.bridge_notes.iter().cloned());
                let actual = fresh.product.as_ref().cloned();
                let tracked = fresh.tracked.as_ref().cloned();
                physical_processing_reports.push(PhysicalProcessingStepReport {
                    step_index: index,
                    evidence: PhysicalProcessingEvidence::TopoActions(fresh),
                });
                let (Some(actual), Some(tracked)) = (actual, tracked) else {
                    return Err(Error::Unsupported(format!(
                        "Step {} TopoActions material is withheld: no faithfully representable paired duplex. Full preparation, parent/sealed graphs, flaps, components and terminal evidence are retained; no target-derived reference product was substituted.",
                        index + 1
                    )));
                };
                product_bytes = product_bytes
                    .checked_add(tracked_bytes)
                    .ok_or(Error::TooLarge)?;
                next_physical = Some(tracked);
                actual
            }
            Operation::ProcessedTopo { sources, request } => {
                // Converting a private tracked descendant to a plain sequence
                // would erase known/unknown modification and graph ancestry.
                if sources.iter().any(|s| matches!(s.source, Source::Step(_))) {
                    return Err(Error::Unsupported("ProcessedTopo composition currently requires original source inputs. Tracked Step ancestry is not mapped through its declared preparation; use the existing private tracked state directly for downstream digest/repair instead of resetting chemical state.".into()));
                }
                let molecules = sources
                    .iter()
                    .map(|s| resolve(&s.source))
                    .collect::<Result<Vec<_>, _>>()?;
                for molecule in &molecules {
                    if chemistry::has_unmapped_physical_evidence(molecule)
                        || molecule.metadata.contains_key("tracked_modifications")
                        || molecule
                            .metadata
                            .contains_key(crate::recombination::REFERENCE_ONLY_KEY)
                        || molecule.metadata.keys().any(|key| key.starts_with("topo:"))
                    {
                        return Err(Error::Unsupported("ProcessedTopo requires original reference sources with explicit preparation declarations; imported physical/reference-only summaries cannot restore their lost graph or tracked modification state.".into()));
                    }
                }
                let inputs = topo_inputs(sources, &molecules, &request.treatment.preparation)?;
                let fresh = crate::topo_processed_product::inspect(&inputs, request)?;
                let fresh_bytes =
                    evidence_bytes(&fresh, crate::topo_processed_product::MAX_REPORT_BYTES)?;
                step_evidence_bytes = step_evidence_bytes
                    .checked_add(fresh_bytes)
                    .filter(|n| *n <= 32 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
                let selected_bases = fresh
                    .product
                    .as_ref()
                    .map_or(0, |m| m.sequence.len())
                    .checked_add(
                        fresh
                            .tracked
                            .as_ref()
                            .map_or(0, |t| t.top.len() + t.bottom.len()),
                    )
                    .ok_or(Error::TooLarge)?;
                let tracked_bytes = fresh
                    .tracked
                    .as_ref()
                    .map(|t| evidence_bytes(t, 32 * 1024 * 1024))
                    .transpose()?
                    .unwrap_or(0);
                let selected_document_bytes = fresh
                    .product
                    .as_ref()
                    .map(|m| evidence_bytes(m, 32 * 1024 * 1024))
                    .transpose()?
                    .unwrap_or(0);
                if product_bases
                    .checked_add(selected_bases)
                    .is_none_or(|n| n > 32_000_000)
                    || product_bytes
                        .checked_add(step_evidence_bytes)
                        .and_then(|n| n.checked_add(tracked_bytes))
                        .and_then(|n| n.checked_add(selected_document_bytes))
                        .is_none_or(|n| n > 64 * 1024 * 1024)
                {
                    return Err(Error::Unsupported("ProcessedTopo selected material plus full retained graph exceeds the plan material/evidence budget before duplicating its selected product.".into()));
                }
                let graph = &fresh.projection.sealing.product;
                findings.push(format!("Step {} fresh ProcessedTopo retains {} components, {} original source calls, {} remaining graph nicks and {} unpaired calls. Explicit component selection is not physical isolation or a mixture winner; all other components and released activations remain typed evidence.",index+1,graph.components.len(),fresh.projection.sealing.source_bases,graph.nicks.len(),graph.unpaired.iter().map(|r|r.end-r.start).sum::<usize>()));
                if !fresh
                    .projection
                    .sealing
                    .parent
                    .all_requested_transfers_applied
                    || fresh
                        .projection
                        .sealing
                        .steps
                        .iter()
                        .any(|s| s.outcome != crate::topo::physical::Outcome::Applied)
                {
                    findings.push(format!("Step {} retains an unapplied/unknown requested branch under its explicit projection options; target agreement cannot establish that operation occurred.",index+1));
                }
                findings.extend(fresh.bridge_notes.iter().cloned());
                let actual = fresh.product.as_ref().cloned();
                let tracked = fresh.tracked.as_ref().cloned();
                processed_step_reports.push(ProcessedTopoStepReport {
                    step_index: index,
                    evidence: fresh,
                });
                let (Some(actual), Some(tracked)) = (actual, tracked) else {
                    return Err(Error::Unsupported(format!(
                        "Step {} ProcessedTopo material is withheld: no faithfully representable paired duplex. Full preparation, parent/sealed graphs, flaps, components and terminal evidence are retained; no target-derived reference product was substituted.",
                        index + 1
                    )));
                };
                product_bytes = product_bytes
                    .checked_add(tracked_bytes)
                    .ok_or(Error::TooLarge)?;
                next_physical = Some(tracked);
                actual
            }
            Operation::Topo {
                sources,
                request,
                selection,
            }
            | Operation::CheckedTopo {
                sources,
                request,
                selection,
                ..
            } => {
                let molecules = sources
                    .iter()
                    .map(|source| resolve(&source.source))
                    .collect::<Result<Vec<_>, _>>()?;
                let inputs = topo_inputs(sources, &molecules, request)?;
                // Only a fresh source-backed reconstruction may populate this
                // evidence. Imported Product/Selection data grants no authority.
                let (fresh, outcome) = if let Operation::CheckedTopo { policy, .. } = step {
                    let assessment = crate::topo::check_target_annotated(
                        &inputs, request, selection, target, policy,
                    )
                    .map_err(|error| {
                        Error::Format(format!("TOPO target assessment rejected: {error}"))
                    })?;
                    (assessment.product, Some(assessment.outcome))
                } else {
                    (
                        crate::topo::annotated_product(&inputs, request, selection).map_err(
                            |error| Error::Format(format!("TOPO reconstruction rejected: {error}")),
                        )?,
                        None,
                    )
                };
                let fully_specified = fresh.physical.projection.physical.fully_specified_dna;
                let physical_bytes =
                    evidence_bytes(&fresh.physical, crate::topo::MAX_REPORT_BYTES)?;
                let annotation_bytes = fresh
                    .annotations
                    .as_ref()
                    .map(|evidence| evidence_bytes(evidence, crate::topo::MAX_REPORT_BYTES))
                    .transpose()?
                    .unwrap_or(0);
                step_evidence_bytes = step_evidence_bytes
                    .checked_add(physical_bytes)
                    .and_then(|bytes| bytes.checked_add(annotation_bytes))
                    .ok_or(Error::TooLarge)?;
                if step_evidence_bytes > 32 * 1024 * 1024 {
                    return Err(Error::Unsupported("TOPO full-graph evidence exceeds the plan's 32 MiB evidence budget; no reference-only substitute is published".into()));
                }
                let crate::topo::AnnotatedProduct {
                    physical,
                    annotations,
                    reference: preview,
                    ..
                } = fresh;
                step_reports.push(StepReport {
                    step_index: index,
                    topo: Some(physical),
                });
                if let Some(evidence) = annotations {
                    findings.push(format!("Step {} retains {} projected TOPO annotations and per-role unmapped/review evidence. These are source-coordinate correspondences, not inherited chemical ancestry or biological validity.",index+1,evidence.features.len()));
                    annotation_reports.push(AnnotationStepReport {
                        step_index: index,
                        evidence,
                    });
                }
                findings.push(format!("Step {} retains its complete freshly reconstructed TOPO preparation, transfer graph, all physical components and paired-contour proof. A selected reference view does not seal nicks, remove flaps or establish experimental success.",index+1));
                match &outcome {
                    Some(crate::topo::TargetOutcome::Rejected { reasons }) => return Err(Error::Format(format!("CheckedTopo independent-target policy rejected: {}", reasons.join("; ")))),
                    Some(crate::topo::TargetOutcome::Withheld { reason }) => return Err(Error::Unsupported(format!("CheckedTopo target assessment withheld: {reason}"))),
                    Some(crate::topo::TargetOutcome::AcceptedUnderDeclaredGraph) if !fully_specified => return Err(Error::Unsupported("CheckedTopo graph retains ambiguous DNA; target agreement cannot establish fully determined sequence.".into())),
                    _ => {}
                }
                let Some(preview) = preview else {
                    return Err(Error::Unsupported(format!(
                        "Step {} has no uniquely proven paired TOPO contour. Its full physical graph is retained; no product DNA or target-derived path is invented.",
                        index + 1
                    )));
                };
                if outcome.is_none() {
                    incomplete.push(format!("Step {} TOPO contour matches only an explicitly selected reference path. No physical target-acceptance policy was supplied, and downstream graph consumption is not modeled; nicks, other components and transfer outcomes remain in full step evidence.",index+1));
                }
                preview
            }
            Operation::RecombineReference {
                sources,
                request,
                selection,
            } => {
                let molecules = [resolve(&sources[0])?, resolve(&sources[1])?];
                // Borrowed preflight happens before cloning the immutable source snapshots.
                let inputs = recombination_inputs(sources, molecules, request)?;
                let fresh = crate::recombination::simulate(&inputs, request)?;
                if fresh.batch_revision != selection.batch_revision {
                    return Err(Error::Format("Reference component selection belongs to a different source/profile/options request".into()));
                }
                findings.push(format!("Step {} reference search: {} candidate(s), {} site(s), complete under declared profile: {}. Profile snapshot {} retains all reciprocal components; selected candidate {} component {} only.", index+1, fresh.reference.candidates.len(), fresh.reference.sites.len(), fresh.reference.complete, request.profile.revision, selection.candidate_id, selection.component_id));
                findings.extend(fresh.reference.diagnostics);
                findings.extend(fresh.warnings);
                let selected = fresh
                    .documents
                    .into_iter()
                    .find(|p| &p.selection == selection)
                    .ok_or_else(|| {
                        Error::Format(
                            "Selected reference candidate/component was not reconstructed".into(),
                        )
                    })?;
                findings.extend(selected.warnings.iter().cloned());
                incomplete.push(format!("Step {} selects candidate {} component {} under pinned reference profile {} / {} ({}). All reciprocal components were reconstructed; this establishes reference sequence only, not recombinase efficacy, physical strand/modification states, cellular repair, selection or chemical uniqueness.", index + 1, selection.candidate_id, selection.component_id, request.profile.profile.id, request.profile.profile.edition, request.profile.revision));
                selected.molecule
            }
            Operation::Pcr {
                template,
                forward,
                reverse,
            } => {
                let substrate = resolve(template)?;
                match template {
                    Source::Input(_) => {
                        if chemistry::has_unmapped_physical_evidence(substrate) {
                            return Err(Error::Unsupported("PCR of imported physical-intermediate evidence requires reconstruction of its originating graph".into()));
                        }
                    }
                    Source::Step(i) => {
                        let duplex = physical.get(*i).and_then(Option::as_ref).ok_or_else(|| {
                            Error::Unsupported("PCR requires an established physical template from the originating step; a reference-only assembly cannot be materialized by PCR".into())
                        })?;
                        pcr_material_preflight(substrate, duplex, forward, reverse)?;
                        if duplex.bottom_start != 0 || duplex.bottom.len() != duplex.top.len() {
                            findings.push(format!("Step {} exact PCR seeds from an actual intact top strand; initially staggered complementary ends remain on the template and are not silently repaired. Known strand extensions were checked for unsupported additional geometries.",index+1));
                        }
                    }
                }
                let matches = analysis::simulate_pcr_geometry(substrate, forward, reverse)?;
                if matches.is_empty() {
                    return Err(Error::Format(format!("step {}: no PCR product", index + 1)));
                }
                if matches.len() != 1 {
                    return Err(Error::Unsupported(format!(
                        "step {}: {} PCR products require disambiguation",
                        index + 1,
                        matches.len()
                    )));
                }
                let geometry = &matches[0];
                let duplex = Duplex::exact_pcr(&geometry.product.sequence)?;
                let mut product = duplex.molecule(format!("PCR {}", index + 1))?;
                project_pcr_annotations(substrate, geometry, &mut product, findings)?;
                next_physical = Some(TrackedDuplex::exact_pcr(&duplex.top)?);
                product
            }
            Operation::PcrWithMismatches {
                template,
                forward,
                reverse,
                options,
            } => {
                let substrate = resolve(template)?;
                let physical_template = match template {
                    Source::Input(_) => {
                        if chemistry::has_unmapped_physical_evidence(substrate) {
                            return Err(Error::Unsupported("Mismatched PCR of imported physical evidence requires its originating graph".into()));
                        }
                        None
                    }
                    Source::Step(i) => {
                        let duplex = physical.get(*i).and_then(Option::as_ref).ok_or_else(|| {
                            Error::Unsupported("Mismatched PCR requires an established physical template from the originating step; a reference-only assembly cannot be materialized by PCR".into())
                        })?;
                        if duplex.bottom_start != 0
                            || duplex.bottom.len() != duplex.top.len()
                            || !duplex.nicks.is_empty()
                        {
                            return Err(Error::Unsupported(
                                "Mismatched PCR of staggered/nicked templates is not modeled"
                                    .into(),
                            ));
                        }
                        Some(duplex)
                    }
                };
                let report = match physical_template {
                    Some(duplex) => crate::mutagenesis::simulate_pcr_reconstructed(
                        substrate, duplex, forward, reverse, options,
                    )?,
                    None => crate::mutagenesis::simulate_pcr(substrate, forward, reverse, options)?,
                };
                findings.extend(
                    report
                        .warnings
                        .iter()
                        .take(20)
                        .map(|s| s.chars().take(2048).collect()),
                );
                if report.products.len() != 1 {
                    return Err(
                        if report.enumeration_complete && report.products.is_empty() {
                            Error::Format(
                                "No mismatched PCR product under the recorded binding policy"
                                    .into(),
                            )
                        } else {
                            Error::Unsupported(format!(
                                "{} mismatched PCR geometries require disambiguation; enumeration complete: {}",
                                report.products.len(),
                                report.enumeration_complete
                            ))
                        },
                    );
                }
                let geometry = &report.products[0];
                if !geometry.certain || geometry.overlapping_domains {
                    return Err(Error::Unsupported("Ambiguous or overlapping-primer products remain previews, not established amplified duplexes".into()));
                }
                if !report.enumeration_complete {
                    incomplete.push(
                        "Mismatched PCR alternative-product search or geometry is incomplete"
                            .to_string(),
                    );
                }
                let duplex = Duplex::exact_pcr(&geometry.sequence)?;
                let mut product = duplex.molecule(format!("Mismatched PCR {}", index + 1))?;
                project_mutant_pcr(substrate, geometry, &mut product, findings)?;
                findings.push(format!("Mismatched PCR introduced {} template-base changes under explicit 3′/mismatch policy; amplification efficiency is not established.",geometry.substitutions.len()));
                next_physical = Some(TrackedDuplex::exact_pcr(&duplex.top)?);
                product
            }
            Operation::GoldenGate {
                templates,
                options,
                candidate_index,
            }
            | Operation::SelectedGoldenGate {
                templates,
                options,
                candidate_index,
                ..
            } => {
                let sources = templates
                    .iter()
                    .map(|s| resolve(s).cloned())
                    .collect::<Result<Vec<_>, _>>()?;
                let substrates = templates
                    .iter()
                    .map(&resolve_tracked)
                    .collect::<Result<Vec<_>, _>>()?;
                let fresh_gate = crate::golden_gate::simulate_tracked(
                    &sources,
                    &substrates,
                    options,
                    Some(*candidate_index),
                    cleavage_budget,
                )?;
                let selected_tracked = fresh_gate.selected;
                let mut local_report = Some(fresh_gate.report);
                let report = if index + 1 == plan.steps.len() {
                    *final_golden_gate = local_report.take();
                    final_golden_gate.as_ref().expect("fresh final enumeration")
                } else {
                    local_report
                        .as_ref()
                        .expect("fresh intermediate enumeration")
                };
                let selected=report.candidates.get(*candidate_index).ok_or_else(||if report.enumeration_complete{Error::Format("Golden Gate candidate index does not exist".into())}else{Error::Unsupported("Golden Gate search is incomplete and the selected candidate is unavailable".into())})?;
                if let Operation::SelectedGoldenGate {
                    directed_sequence_key,
                    physical_sequence_key,
                    ..
                } = step
                {
                    if directed_sequence_key.len() != 64
                        || physical_sequence_key.len() != 64
                        || selected.directed_sequence_key != *directed_sequence_key
                        || selected.physical_sequence_key != *physical_sequence_key
                    {
                        return Err(Error::Format("Golden Gate explicit candidate keys disagree with the fresh native enumeration".into()));
                    }
                }
                if !selected.cleavage.complete {
                    return Err(Error::Unsupported("Golden Gate candidate site absence is unproven: bounded recut evidence is incomplete; no site-free cycling-product claim is established".into()));
                }
                if !selected.restriction_site_free {
                    return Err(Error::Format("Golden Gate candidate retains recognition sites; it is not a site-free cycling product".into()));
                }
                if !report.enumeration_complete {
                    incomplete.push(format!(
                        "Golden Gate enumeration incomplete: {}",
                        report.truncation_reasons.join(" ")
                    ));
                }
                findings.push(format!("Golden Gate checked candidate {} of {} enumerated paths; {} stable physical sequence groups, {} uncut carryover sources. This finite one-copy model does not establish exhaustive reaction-mixture uniqueness.",candidate_index+1,report.candidates.len(),report.stable_physical_product_count,report.carryover_sources.len()));
                findings.extend(
                    report
                        .warnings
                        .iter()
                        .take(20)
                        .map(|s| s.chars().take(2048).collect()),
                );
                let duplex = selected.duplex.clone();
                let mut product = duplex.molecule(format!("Golden Gate {}", index + 1))?;
                let mut offset = 0usize;
                for oriented in &selected.path {
                    let fragment = &report.fragments[oriented.fragment];
                    let source = &sources[fragment.source_index];
                    let mut intermediate =
                        fragment.duplex.molecule("Golden Gate source fragment")?;
                    let map = region_map(
                        source.sequence.len(),
                        fragment.top_start,
                        fragment.duplex.top.len(),
                        0,
                    )?;
                    project_into(source, &map, &mut intermediate, findings)?;
                    let oriented_duplex = if oriented.reverse {
                        fragment.duplex.reverse()?
                    } else {
                        fragment.duplex.clone()
                    };
                    let map = if oriented.reverse {
                        reverse_reference_map(
                            intermediate.sequence.len(),
                            fragment.duplex.bottom_start + fragment.duplex.bottom.len() as i64,
                            offset,
                            product.sequence.len(),
                            product.topology,
                        )?
                    } else {
                        vec![MappingSegment {
                            source_start: 0,
                            source_end: intermediate.sequence.len(),
                            destination_start: offset,
                            reverse: false,
                        }]
                    };
                    project_into(&intermediate, &map, &mut product, findings)?;
                    offset += oriented_duplex.top.len();
                }
                next_physical=Some(selected_tracked.ok_or_else(||Error::Unsupported("Selected Golden Gate material was not freshly reconstructed under checked cleavage".into()))?);
                product
            }
            Operation::Overlap {
                fragments,
                overlaps,
                circular,
            } => {
                if fragments.is_empty()
                    || fragments.len() > 1000
                    || overlaps.len() != fragments.len() - 1 + usize::from(*circular)
                {
                    return Err(Error::Format(
                        "ordered assembly requires one overlap per junction".into(),
                    ));
                }
                for source in fragments {
                    if let Source::Input(_) = source {
                        let input = resolve(source)?;
                        if chemistry::has_unmapped_physical_evidence(input) {
                            return Err(Error::Unsupported("overlap assembly of imported physical-intermediate evidence requires reconstruction of its originating graph".into()));
                        }
                    }
                    if let Source::Step(i) = source {
                        if let Some(Some(duplex)) = physical.get(*i) {
                            let ends = duplex.ends()?;
                            if duplex.left.three_prime == chemistry::ThreePrime::Blocked
                                || duplex.right.three_prime == chemistry::ThreePrime::Blocked
                            {
                                return Err(Error::Unsupported("overlap assembly of a known blocked3′ terminus requires an explicitly modeled unblocking/resection mechanism".into()));
                            }
                            if !ends.0.overhang.is_empty()
                                || !ends.1.overhang.is_empty()
                                || !duplex.nicks.is_empty()
                            {
                                return Err(Error::Unsupported("overlap assembly of staggered or nicked physical substrates is not yet modeled".into()));
                            }
                        }
                    }
                }
                let fragments = fragments
                    .iter()
                    .map(resolve)
                    .collect::<Result<Vec<_>, _>>()?;
                assemble_overlap(&fragments, overlaps, *circular, index, findings)?
            }
            Operation::DeclareSourceTermini {
                template,
                left,
                right,
            } => {
                if !matches!(template, Source::Input(_)) {
                    return Err(Error::Unsupported("Source termini may only be declared on an original input; descendants require their actual physical graph, not a new intact-source hypothesis.".into()));
                }
                let substrate = resolve(template)?;
                declared_source_termini_preflight(substrate)?;
                // Start from the same guarded intact reference assumption as
                // other input chemistry; summaries never supply strand state.
                let mut duplex = Duplex::from_molecule(substrate)?;
                duplex.left = *left;
                duplex.right = *right;
                duplex.validate()?;
                let tracked = TrackedDuplex::from_unknown_duplex(&duplex)?;
                // DNA and all source feature/citation evidence are unchanged.
                let mut product = substrate.clone();
                product.metadata.insert(
                    "duplex_model".into(),
                    serde_json::to_string(&duplex.summary()?)
                        .map_err(|e| Error::Format(e.to_string()))?,
                );
                product.metadata.insert("source:declared_termini".into(),serde_json::to_string(&(substrate.revision(),left,right,"Operator hypothesis for intact original linear DNA; not measured chemistry, synthesis, repair or proof of material state.")).map_err(|e|Error::Format(e.to_string()))?);
                findings.push(format!("Step {} explicitly assumes original intact linear duplex terminal groups: left {:?}, right {:?}. This is an operator hypothesis, not measured chemistry, synthesis, nick repair or recovery of imported physical state. Original DNA/features/citations are retained; positional modification states remain unknown.",index+1,left,right));
                next_physical = Some(tracked);
                product
            }
            Operation::DeclareSourceModifications {
                template,
                assumptions,
            } => {
                if !matches!(template, Source::Input(_)) {
                    return Err(Error::Unsupported("Source-generation assumptions may only bind an original input, never rescan a reconstructed descendant".into()));
                }
                let substrate = resolve(template)?;
                let tracked = TrackedDuplex::from_input(substrate, assumptions)?;
                let mut product = tracked
                    .duplex
                    .molecule(format!("Declared source generation {}", index + 1))?;
                project_into(
                    substrate,
                    &[MappingSegment {
                        source_start: 0,
                        source_end: substrate.sequence.len(),
                        destination_start: 0,
                        reverse: false,
                    }],
                    &mut product,
                    findings,
                )?;
                findings.push(format!("Step {} records a closed Dam/Dcm/CpG source-generation assumption on exact input revision {}: {}. This does not measure modifications or infer any other methylase activity.",index+1,substrate.revision(),serde_json::to_string(assumptions).map_err(|e|Error::Format(e.to_string()))?));
                next_physical = Some(tracked);
                product
            }
            Operation::CheckedDigest {
                template,
                selection,
                options,
                fragment_index,
            } => {
                options.validate()?;
                let tracked = resolve_tracked(template)?;
                let report = modified_duplex::checked_digest_with_budget(
                    &tracked,
                    &CheckedDigestOptions {
                        selection: selection.clone(),
                        max_comparisons: usize::try_from(options.max_comparisons)
                            .map_err(|_| Error::TooLarge)?,
                        max_sites: options.max_sites,
                    },
                    cleavage_budget,
                )?;
                evidence_bytes(&report, options.max_report_bytes)?;
                if report.fragments.len() > options.max_fragments {
                    return Err(Error::TooLarge);
                }
                for site in report.sites.iter().take(100) {
                    findings.push(format!(
                        "Step {} checked {} site at source base {}: {:?}. {}",
                        index + 1,
                        site.enzyme,
                        site.recognition_start + 1,
                        site.outcome,
                        site.reason
                    ));
                }
                findings.extend(report.warnings.iter().take(100).cloned());
                if !report.complete || !report.products_determined {
                    return Err(Error::Unsupported(format!(
                        "Checked digest cannot establish complete products under the pinned positional profile; unknown/impaired or bounded evidence remains. {}",
                        report.limits.join(" ")
                    )));
                }
                let fragment=report.fragments.get(*fragment_index).ok_or_else(||Error::Format(format!("checked digest has {} fragments; selected zero-based index {fragment_index} is unavailable",report.fragments.len())))?;
                let mut product = fragment.duplex.duplex.molecule(format!(
                    "Checked digest {} · fragment {}",
                    index + 1,
                    fragment_index + 1
                ))?;
                let substrate = resolve(template)?;
                project_into(
                    substrate,
                    &region_map(
                        substrate.sequence.len(),
                        fragment.top_start,
                        fragment.duplex.duplex.top.len(),
                        0,
                    )?,
                    &mut product,
                    findings,
                )?;
                findings.push(format!("Step {} selected fragment {} of {} fully determined positional-profile digest products; supplier conditions and experimental completion remain assumptions.",index+1,fragment_index+1,report.fragments.len()));
                next_physical = Some(fragment.duplex.clone());
                product
            }
            Operation::Digest {
                template,
                fragment_index,
                ..
            }
            | Operation::CatalogDigest {
                template,
                fragment_index,
                ..
            } => {
                let source = resolve_tracked(template)?;
                let enzymes = match step {
                    Operation::CatalogDigest { selection, .. } => {
                        let catalog = crate::catalog_digest::builtin_catalog()?;
                        let resolved = selection.resolved_enzymes(&catalog).map_err(|error| {
                            Error::Format(format!("Catalog digest pin rejected: {error}"))
                        })?;
                        findings.push(format!("Step {} uses explicit embedded catalog edition {} / content {} and exact ordered enzyme IDs {}. Selected identity/geometry is validated against the catalog; documentary supplier conditions do not establish complete cleavage.",index+1,selection.catalog_edition,selection.catalog_content_hash,selection.enzymes.iter().map(|e|e.id.as_str()).collect::<Vec<_>>().join(", ")));
                        resolved
                    }
                    Operation::Digest { enzymes, .. } => {
                        let catalog = analysis::starter_enzymes();
                        enzymes
                            .iter()
                            .map(|name| {
                                catalog
                                    .iter()
                                    .find(|e| e.name.eq_ignore_ascii_case(name))
                                    .cloned()
                                    .ok_or_else(|| {
                                        Error::Unsupported(format!(
                                            "enzyme {name} is not in the embedded catalog"
                                        ))
                                    })
                            })
                            .collect::<Result<Vec<_>, _>>()?
                    }
                    _ => unreachable!("matched digest operation"),
                };
                if !source.generation.is_empty() {
                    let state = modified_duplex::summary(&source)?;
                    if state.known_modified > 0 || state.unknown_eligible > 0 {
                        incomplete.push(format!("Step {} uses legacy assumed-complete cleavage on explicitly tracked modified/unknown material. Use CheckedDigest to establish cleavage under the positional profile.",index+1));
                    }
                }
                let fragments = modified_duplex::digest(&source, &enzymes)?;
                let fragment=fragments.get(*fragment_index).ok_or_else(||Error::Format(format!("digest has {} fragments; selected zero-based index {fragment_index} does not exist",fragments.len())))?;
                let mut product = fragment.duplex.molecule(format!(
                    "Digest {} · fragment {}",
                    index + 1,
                    fragment_index + 1
                ))?;
                let template_molecule = resolve(template)?;
                let map = region_map(
                    template_molecule.sequence.len(),
                    fragment.top_start,
                    fragment.duplex.top.len(),
                    0,
                )?;
                project_into(template_molecule, &map, &mut product, findings)?;
                findings.push(format!("Step {} selects fragment {} of {} from the complete digest (top reference boundary {}); top/bottom lengths {}/{}.",index+1,fragment_index+1,fragments.len(),fragment.top_start,fragment.duplex.top.len(),fragment.duplex.bottom.len()));
                next_physical = Some(fragment.duplex.clone());
                product
            }
            Operation::TreatEnds {
                template,
                treatment,
            } => {
                let tracked = resolve_tracked(template)?.treat(*treatment)?;
                let duplex = &tracked.duplex;
                let mut product = duplex.molecule(format!("End treatment {}", index + 1))?;
                let source = resolve(template)?;
                project_into(
                    source,
                    &[MappingSegment {
                        source_start: 0,
                        source_end: source.sequence.len(),
                        destination_start: 0,
                        reverse: false,
                    }],
                    &mut product,
                    findings,
                )?;
                next_physical = Some(tracked);
                product
            }
            Operation::Tail3Prime { template, base } => {
                let source = resolve_physical(template)?;
                let report = crate::tailing::tail_blunt_duplex(&source, *base)?;
                let product = tailed_document(
                    resolve(template)?,
                    &report,
                    format!("3′ tailing {}", index + 1),
                    findings,
                )?;
                findings.push(format!("Step {} explicitly appends one {} nucleotide to both strand 3′ termini. Added nucleotides have no source-template ancestry; five-prime chemistry and internal nicks are retained.",index+1,base.as_str()));
                findings.extend(report.warnings.iter().take(100).cloned());
                next_physical = Some(resolve_tracked(template)?.tail(*base)?);
                product
            }
            Operation::ProfiledOverlap {
                fragments,
                overlaps,
                circular,
                options,
            } => {
                let mut sources = Vec::new();
                let mut reverse = Vec::new();
                for oriented in fragments {
                    let source = resolve(&oriented.source)?;
                    if let Source::Input(_) = &oriented.source {
                        if chemistry::has_unmapped_physical_evidence(source) {
                            return Err(Error::Unsupported("profiled overlap assembly of imported physical evidence requires reconstruction of its originating graph".into()));
                        }
                    }
                    if let Source::Step(i) = &oriented.source {
                        if let Some(Some(d)) = physical.get(*i) {
                            let ends = d.ends()?;
                            if d.left.three_prime == chemistry::ThreePrime::Blocked
                                || d.right.three_prime == chemistry::ThreePrime::Blocked
                            {
                                return Err(Error::Unsupported("profiled overlap assembly of a known blocked3′ terminus requires an explicitly modeled unblocking/resection mechanism".into()));
                            }
                            if !ends.0.overhang.is_empty()
                                || !ends.1.overhang.is_empty()
                                || !d.nicks.is_empty()
                            {
                                return Err(Error::Unsupported("profiled overlap assembly of staggered or nicked physical substrates is not modeled".into()));
                            }
                        }
                    }
                    sources.push(source.clone());
                    reverse.push(oriented.reverse);
                }
                let report = crate::overlap_profiles::check_oriented_fragments(
                    &sources, &reverse, overlaps, *circular, options,
                )?;
                if !report.geometry_valid {
                    return Err(Error::Format(format!(
                        "profiled overlap geometry is invalid: {}",
                        report.findings.join(" ")
                    )));
                }
                if !report.complete {
                    incomplete.push(format!("Step {} profile repeat scan reached its configured bound; alternative placements remain unsearched.",index+1));
                }
                findings.extend(report.findings.iter().take(100).cloned());
                findings.extend(report.warnings.iter().take(100).cloned());
                findings.push(format!("Step {} uses {:?} sequence recommendations. Recommendation compliance {:?} is advisory and does not establish reagent activity, experimental outcome or uniqueness.",index+1,options.method,report.recommendation_compliant));
                let mut oriented_sources = Vec::new();
                for (source, reverse) in sources.iter().zip(&reverse) {
                    if *reverse {
                        let mut oriented = Molecule::dna(
                            format!("Oriented {}", source.name),
                            source.sequence.reverse_complement()?.as_str(),
                            source.topology,
                        )?;
                        project_into(
                            source,
                            &[MappingSegment {
                                source_start: 0,
                                source_end: source.sequence.len(),
                                destination_start: 0,
                                reverse: true,
                            }],
                            &mut oriented,
                            findings,
                        )?;
                        carry_unmapped_evidence(source, &mut oriented);
                        oriented_sources.push(oriented);
                    } else {
                        oriented_sources.push(source.clone());
                    }
                }
                let refs = oriented_sources.iter().collect::<Vec<_>>();
                let mut product = assemble_overlap(&refs, overlaps, *circular, index, findings)?;
                product.metadata.insert(
                    "overlap_profile:options".into(),
                    serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?,
                );
                product.metadata.insert(
                    "overlap_profile:source_revisions".into(),
                    serde_json::to_string(&report.source_revisions)
                        .map_err(|e| Error::Format(e.to_string()))?,
                );
                product.metadata.insert(
                    "overlap_profile:orientations".into(),
                    serde_json::to_string(&reverse).map_err(|e| Error::Format(e.to_string()))?,
                );
                product
            }
            Operation::RepairEnds { template, method } => {
                let source = resolve_physical(template)?;
                let repair = chemistry::repair_ends(&source, *method)?;
                let product = repaired_document(
                    resolve(template)?,
                    &repair,
                    format!("End repair {}", index + 1),
                    findings,
                )?;
                findings.push(format!("Step {} end repair {:?}: top added/removed {}/{} bases; bottom added/removed {}/{} bases. Five-prime chemistry and internal nicks are retained.",index+1,method,repair.top_added_3prime.len(),repair.top_removed_3prime.len(),repair.bottom_added_3prime.len(),repair.bottom_removed_3prime.len()));
                findings.extend(repair.warnings.iter().take(100).cloned());
                next_physical = Some(resolve_tracked(template)?.repair(*method)?);
                product
            }
            Operation::AnnealOligos {
                top,
                bottom,
                options,
                candidate_index,
            } => {
                if !matches!(top, Source::Input(_)) || !matches!(bottom, Source::Input(_)) {
                    return Err(Error::Unsupported("annealing requires explicitly declared input oligos; duplex-producing steps do not implicitly denature into ssDNA".into()));
                }
                let top = resolve(top)?;
                let bottom = resolve(bottom)?;
                let preview = simulate_oligo_annealing(top, bottom, options)?;
                let candidate = preview
                    .annealing
                    .candidates
                    .get(*candidate_index)
                    .ok_or_else(|| {
                        Error::Format("selected annealing candidate is unavailable".into())
                    })?;
                let product = annealed_document_inner(
                    top,
                    bottom,
                    &preview.annealing,
                    *candidate_index,
                    findings,
                )?;
                findings.push(format!("Step {} anneals declared 5′→3′ ssDNA oligos at bottom offset {}, with {} paired bases. Explicit terminal chemistry is a plan assumption; annealing creates no covalent bonds.",index+1,candidate.duplex.bottom_start,candidate.paired_bases));
                if !preview.annealing.complete {
                    incomplete.push(format!("Step {} annealing placement search was incomplete; unsearched alignments cannot establish uniqueness.",index+1));
                }
                if !preview.annealing.unique {
                    incomplete.push(format!("Step {} annealing has competing or unproven placements; an explicit selected candidate does not establish a unique annealed product.",index+1));
                }
                findings.extend(preview.annealing.warnings.iter().take(100).cloned());
                next_physical = Some(TrackedDuplex::annealed_unmodified(&candidate.duplex)?);
                product
            }
            Operation::Ligate {
                fragments,
                circular,
                allow_nicked,
                options,
            } => {
                if fragments.is_empty() || fragments.len() > 1000 {
                    return Err(Error::Format(
                        "ligation requires1–1000 explicitly oriented fragments".into(),
                    ));
                }
                options.validate()?;
                if fragments.len() > crate::ligation_outcomes::MAX_FRAGMENTS {
                    return Err(Error::TooLarge);
                }
                let mut raw = Vec::with_capacity(fragments.len());
                let mut oriented_lengths = Vec::with_capacity(fragments.len());
                let mut bases = 0usize;
                for source in fragments {
                    let duplex = resolve_tracked(&source.source)?;
                    bases = bases
                        .saturating_add(duplex.top.len())
                        .saturating_add(duplex.bottom.len());
                    if bases > 2 * MAX_SEQUENCE_LEN {
                        return Err(Error::TooLarge);
                    }
                    oriented_lengths.push(if source.reverse {
                        duplex.bottom.len()
                    } else {
                        duplex.top.len()
                    });
                    raw.push(duplex);
                }
                let (duplex, outcomes_receipt) = {
                    if physical_processing_reports.len() >= 64 {
                        return Err(Error::TooLarge);
                    }
                    let keys = (0..fragments.len())
                        .map(|i| format!("fragment-{}", i + 1))
                        .collect::<Vec<_>>();
                    let operands = keys
                        .iter()
                        .zip(&raw)
                        .map(|(key, material)| crate::ligation_outcomes::Operand { key, material })
                        .collect::<Vec<_>>();
                    let request = crate::ligation_outcomes::Request {
                        intended: keys
                            .iter()
                            .zip(fragments)
                            .map(|(key, source)| crate::ligation_outcomes::Placement {
                                key: key.clone(),
                                reverse: source.reverse,
                            })
                            .collect(),
                        circular: *circular,
                        allow_nicked: *allow_nicked,
                        options: options.clone(),
                    };
                    let (outcomes, actual) =
                        crate::ligation_outcomes::analyze_fresh(&operands, &request)?.into_parts();
                    let bindings = keys
                        .iter()
                        .zip(fragments)
                        .map(|(key, source)| {
                            Ok(LigationOutcomeBinding {
                                key: key.clone(),
                                source: source.source.clone(),
                                source_revision: resolve(&source.source)?.revision(),
                                material_revision: outcomes
                                    .sources
                                    .iter()
                                    .find(|s| &s.key == key)
                                    .ok_or_else(|| {
                                        Error::Format(
                                            "Missing native ligation source identity".into(),
                                        )
                                    })?
                                    .state_revision
                                    .clone(),
                            })
                        })
                        .collect::<Result<Vec<_>, Error>>()?;
                    let evidence = LigationOutcomeEvidence {
                        bindings,
                        report: outcomes,
                    };
                    evidence.validate_for_operation(
                        fragments,
                        *circular,
                        *allow_nicked,
                        options,
                    )?;
                    step_evidence_bytes = step_evidence_bytes.checked_add(evidence_bytes(&evidence,32*1024*1024)?).filter(|n|*n<=32*1024*1024).ok_or_else(||Error::Unsupported("Complete ligation outcomes plus prior processing evidence exceed the32MiB plan evidence bound; reduce explicit outcome bounds before capturing this plan. The intended chemistry was not rejected and no alternatives were silently dropped.".into()))?;
                    let outcomes_receipt = (
                        crate::assembly_recipe::hash(&evidence)?,
                        evidence.report.context_revision.clone(),
                    );
                    findings.push(format!("Step {} retains the declared intended ligation and {} alternative one-copy paths; enumeration complete {}, exhaustive within that scope {}. Alternative existence does not invalidate the intended path or establish isolation; other copy numbers and incomplete assemblies remain outside this scope.",index+1,evidence.report.alternatives.len(),evidence.report.enumeration_complete,evidence.report.exhaustive_within_scope));
                    physical_processing_reports.push(PhysicalProcessingStepReport {
                        step_index: index,
                        evidence: PhysicalProcessingEvidence::LigationOutcomes(evidence),
                    });
                    (actual, outcomes_receipt)
                };
                findings.push(format!("Step {} ligation reconstructs both strands with {} retained nicks; covalently closed circular duplex in this model: {}. No nick repair is assumed.",index+1,duplex.nicks.len(),duplex.topology==Topology::Circular && duplex.nicks.is_empty()));
                let mut product = duplex.molecule(format!("Ligation {}", index + 1))?;
                {
                    let (hash, context) = outcomes_receipt;
                    product
                        .metadata
                        .insert("source:ligation_outcomes:evidence_hash".into(), hash);
                    product
                        .metadata
                        .insert("source:ligation_outcomes:context".into(), context);
                }
                let mut destination = 0usize;
                for (oriented_source, oriented_top_length) in
                    fragments.iter().zip(&oriented_lengths)
                {
                    let source = resolve(&oriented_source.source)?;
                    let map = if oriented_source.reverse {
                        let original = resolve_physical(&oriented_source.source)?;
                        let end = original.bottom_start + original.bottom.len() as i64;
                        reverse_reference_map(
                            source.sequence.len(),
                            end,
                            destination,
                            product.sequence.len(),
                            product.topology,
                        )?
                    } else {
                        vec![MappingSegment {
                            source_start: 0,
                            source_end: source.sequence.len(),
                            destination_start: destination,
                            reverse: false,
                        }]
                    };
                    project_into(source, &map, &mut product, findings)?;
                    destination += *oriented_top_length;
                }
                next_physical = Some(duplex);
                product
            }
        };
        let is_reference_only = diagnostic_reference_only
            || inherited_reference_only
            || matches!(
                step,
                Operation::RecombineReference { .. }
                    | Operation::MultiSiteRecombineReference { .. }
                    | Operation::Topo { .. }
                    | Operation::CheckedTopo { .. }
            );
        if is_reference_only {
            if next_physical.is_some() {
                return Err(Error::Unsupported("Reference-only ancestry cannot acquire chemical authority without a modeled materialization assumption".into()));
            }
            // Restrictive portable evidence only. Actual graph authority lives in the private flag.
            product.metadata.insert(
                crate::recombination::REFERENCE_ONLY_KEY.into(),
                "true".into(),
            );
            product.metadata.insert(
                crate::recombination::CHEMISTRY_UNMAPPED_KEY.into(),
                "true".into(),
            );
        }
        if let Some(tracked) = &next_physical {
            let summary = modified_duplex::summary(tracked)?;
            product.metadata.insert(
                "tracked_modifications".into(),
                serde_json::to_string(&summary).map_err(|e| Error::Format(e.to_string()))?,
            );
        }
        product_bases = product_bases
            .saturating_add(product.sequence.len())
            .saturating_add(
                next_physical
                    .as_ref()
                    .map_or(0, |d| d.top.len() + d.bottom.len()),
            );
        if product_bases > 32_000_000 {
            return Err(Error::Unsupported(
                "plan products exceed the 32 million base materialization budget".into(),
            ));
        }
        product
            .metadata
            .insert("plan_revision".into(), plan_hash.into());
        product
            .metadata
            .insert("engine_version".into(), ENGINE_VERSION.into());
        product_bytes = product_bytes.saturating_add(
            serde_json::to_vec(&product)
                .map_err(|e| Error::Format(e.to_string()))?
                .len(),
        );
        if product_bytes.saturating_add(step_evidence_bytes) > 64 * 1024 * 1024 {
            return Err(Error::Unsupported(
                "plan products exceed the 64 MiB serialized annotation/evidence budget".into(),
            ));
        }
        products.push(product);
        physical.push(next_physical);
        reference_only.push(is_reference_only);
    }
    Ok(incomplete)
}

pub fn escape_html(value: &str) -> String {
    value
        .chars()
        .map(|c| match c {
            '&' => "&amp;".into(),
            '<' => "&lt;".into(),
            '>' => "&gt;".into(),
            '"' => "&quot;".into(),
            '\'' => "&#39;".into(),
            _ => c.to_string(),
        })
        .collect()
}

/// Resolve only a validated bound input or an actually reconstructed preceding
/// step. Invalid/forward/missing references must not invent ordered oligo bases.
pub(crate) fn protocol_oligo_sequence<'a>(
    plan: &'a DesignPlan,
    report: &'a VerificationReport,
    source: &Source,
    step_index: usize,
) -> Option<&'a str> {
    let molecule = match source {
        Source::Input(id) => {
            let mut matches = plan.inputs.iter().filter(|input| &input.id == id);
            let input = matches.next()?;
            if matches.next().is_some() || input.revision != input.molecule.revision() {
                return None;
            }
            &input.molecule
        }
        Source::Step(index) if *index < step_index => report.products.get(*index)?,
        Source::Step(_) => return None,
    };
    if molecule.sequence.alphabet() != Alphabet::Dna
        || molecule.sequence.len() > chemistry::MAX_ANNEAL_OLIGO_BASES
    {
        return None;
    }
    Some(molecule.sequence.as_str())
}

/// Native-generated documentary appendix. It cannot supply molecular status.
pub(crate) struct ProtocolAppendix {
    pub paragraphs: Vec<String>,
    pub tables: Vec<ProtocolTable>,
}
pub(crate) struct ProtocolTable {
    pub heading: String,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
}
impl ProtocolAppendix {
    pub(crate) fn validate(&self) -> Result<(), Error> {
        let mut bytes = 0usize;
        let mut rows = 0usize;
        if self.paragraphs.len() > 8192 || self.tables.len() > 128 {
            return Err(Error::TooLarge);
        }
        let mut charge = |text: &str| -> Result<(), Error> {
            bytes = bytes
                .checked_add(text.len())
                .filter(|n| *n <= 512 * 1024)
                .ok_or(Error::TooLarge)?;
            Ok(())
        };
        for p in &self.paragraphs {
            charge(p)?;
        }
        for table in &self.tables {
            rows = rows
                .checked_add(table.rows.len())
                .filter(|n| *n <= 8192)
                .ok_or(Error::TooLarge)?;
            if table.columns.is_empty()
                || table.columns.len() > 8
                || table.rows.iter().any(|r| r.len() != table.columns.len())
            {
                return Err(Error::Format("invalid protocol appendix table".into()));
            }
            charge(&table.heading)?;
            for cell in table.columns.iter().chain(table.rows.iter().flatten()) {
                charge(cell)?;
            }
        }
        Ok(())
    }
}
/// Fresh in-process protocol context. It cannot be deserialized or constructed
/// from an imported report; immutable borrows keep plan/target paired with this
/// exact native reconstruction while several output formats are rendered.
pub(crate) struct CheckedProtocol<'a> {
    plan: &'a DesignPlan,
    target: &'a Molecule,
    report: VerificationReport,
    appendix: Option<&'a ProtocolAppendix>,
}
impl<'a> CheckedProtocol<'a> {
    pub(crate) fn with_appendix(mut self, appendix: &'a ProtocolAppendix) -> Result<Self, Error> {
        appendix.validate()?;
        self.appendix = Some(appendix);
        Ok(self)
    }
    pub(crate) fn appendix(&self) -> Option<&ProtocolAppendix> {
        self.appendix
    }

    pub(crate) fn plan(&self) -> &'a DesignPlan {
        self.plan
    }
    pub(crate) fn target(&self) -> &'a Molecule {
        self.target
    }
    pub(crate) fn report(&self) -> &VerificationReport {
        &self.report
    }
}
pub(crate) fn checked_protocol<'a>(
    plan: &'a DesignPlan,
    target: &'a Molecule,
) -> CheckedProtocol<'a> {
    // Native verify owns the source/operation guards and independently rebuilds
    // the product. Invalid inputs produce fresh draft findings, never authority.
    CheckedProtocol {
        plan,
        target,
        report: verify(plan, target),
        appendix: None,
    }
}
/// Rechecks the plan at export time; the caller cannot attach an unrelated report.
pub fn export_protocol_html(plan: &DesignPlan, target: &Molecule) -> String {
    export_protocol_html_checked(&checked_protocol(plan, target))
}
pub(crate) fn export_protocol_html_checked(checked: &CheckedProtocol<'_>) -> String {
    let plan = checked.plan();
    let target = checked.target();
    let report = checked.report();
    let e = escape_html;
    let mut body = format!(
        "<h1>{}</h1><p><strong>{:?}</strong> · computed sequence checks; laboratory conditions require review.</p><p>Target: {} · {} bases · {:?}</p><h2>Starting materials</h2><ul>",
        e(&plan.name),
        report.status,
        e(&target.name),
        target.sequence.len(),
        target.topology
    );
    for input in &plan.inputs {
        body.push_str(&format!(
            "<li>Input <strong>{}</strong>: {} · {} bases · {:?}<br>Revision <code>{}</code></li>",
            e(&input.id),
            e(&input.molecule.name),
            input.molecule.sequence.len(),
            input.molecule.topology,
            e(&input.revision)
        ));
    }
    body.push_str("</ul><h2>Cloning steps</h2><ol>");
    let source = |s: &Source| {
        e(&match s {
            Source::Input(id) => format!("input {id}"),
            Source::Step(i) => i
                .checked_add(1)
                .map(|number| format!("product of step {number}"))
                .unwrap_or_else(|| format!("invalid step index {i}")),
        })
    };
    for (step_index, step) in plan.steps.iter().enumerate() {
        match step {
            Operation::ResectThreePrime { source: operand, cuts } => {
                body.push_str(&format!("<li>Explicit3′→5′ resection of {}: remove {} nt from bottom3′ left end and {} nt from top3′ right end. Positive cuts require accessible3′OH; original5′ groups and surviving modifications stay unchanged. No fill, ligase, unblocking or host repair.</li>",source(operand),cuts.left,cuts.right));
            },
            Operation::HomologyAnnealing {fragments,options,policy} => {
                body.push_str(&format!("<li>Generic3′ processing and exposed5′tail annealing of ordered substrates {}. Explicit policy <code>{}</code>; advisory profile <code>{}</code>. Every nick/gap retains separate covalent strands. The reference contour does not establish a reusable assembled molecule or proprietary In-Fusion/cellular rescue.</li>",fragments.iter().map(|f|format!("{} ({})",source(&f.source),if f.reverse{"reverse complement"}else{"forward"})).collect::<Vec<_>>().join(", "),e(&serde_json::to_string(policy).unwrap_or_default()),e(&serde_json::to_string(options).unwrap_or_default())));
            },
            Operation::HomologyAssembly { fragments, options, policy } => {
                body.push_str(&format!("<li>Physical homology assembly of ordered substrates {}. Declared geometry and 5′ processing/fill/sealing policy: <code>{}</code>. Advisory sequence profile: <code>{}</code>. Surviving source strands retain modifications; only actual template-supported synthesis is unmodified. Eligible paired 3′OH/5′P nicks may seal; unresolved chemistry stays explicit. Linear outer-end processing is retained, not silently blunted. This model does not establish experimental success or a proprietary formulation.</li>",fragments.iter().map(|f|format!("{} ({})",source(&f.source),if f.reverse{"reverse complement"}else{"forward"})).collect::<Vec<_>>().join(", "),e(&serde_json::to_string(policy).unwrap_or_default()),e(&serde_json::to_string(options).unwrap_or_default())));
            },
            Operation::OverlapExtension { left, right, options } => {
                body.push_str(&format!("<li>Denature overlapping duplexes {} ({}) and {} ({}), reanneal left top/right bottom strands across {} bp, and extend both paired 3′ ends. Minimum terminal pairing {} bp. Complementary synthesis is assumed unmodified; selected geometry does not establish isolation or experimental success.</li>", source(&left.source), if left.reverse {"reverse complement"} else {"forward"}, source(&right.source), if right.reverse {"reverse complement"} else {"forward"}, options.overlap, options.min_paired_bases));
            },
            Operation::ReconstructMaterial { template, origin } => {
                body.push_str(&format!("<li>Retrieve saved material {} by replaying its exact source and native derivation for event <code>{}</code>. Preserve actual strands, termini, modifications and inherited uncertainty. This retrieves material; it does not rerun historical reactions, create a target assessment or prove isolation.</li>",source(template),e(&origin.product.event_id)));
            },
            Operation::ReconstructRunProduct { template, origin } => {
                body.push_str(&format!("<li>Use saved material {} by reconstructing run <code>{}</code>, product {} against its unchanged original target <code>{}</code>. Preserve actual nicks, termini and positional modifications. Complete origin evidence is retained; this is not a newly performed reaction, inferred repair or proof of isolation.</li>", source(template), e(&origin.selection.run_id), origin.selection.step_index+1, e(&origin.target_event.revision)));
            },
            Operation::TopoActions{sources,request}=>{
                body.push_str(&format!("<li>Accounted TOPO chemical actions from original roles {}. Every original base and excised fragment remains in the full graph. Explicit selected component {} and source seed strand {}, index {} (zero-based); component selection is not isolation. Ideal hydrolytic profiles are declared model assumptions, not a named nuclease, kit or host-repair guarantee.<ul>",sources.iter().map(|r|format!("{} = {}",e(&r.role),source(&r.source))).collect::<Vec<_>>().join(", "),request.selection.component,request.selection.reference_seed.strand,request.selection.reference_seed.index));
                for action in &request.processing.actions {body.push_str(&format!("<li>{}</li>",e(&format!("{action:?}"))));}
                body.push_str("</ul></li>");
            },
            Operation::ProcessedTopo {sources,request} => {
                body.push_str(&format!("<li>Fresh processed TOPO graph from original roles {}. Explicit ligase assumption: {}. Requested nick seals: {}. Select component {} using source call strand {}, index {} (zero-based), bound to <code>{}</code>. Require all requested transfers/seals applied: {}; require known terminal chemistry: {}. Only a faithfully representable paired contour becomes private duplex material; all other components and graph evidence are retained. Component choice is not isolation or a reaction winner; independent final-target verification remains separate.<ul>",sources.iter().map(|s|format!("{} = {}",e(&s.role),source(&s.source))).collect::<Vec<_>>().join(", "),e(&request.treatment.sealing.ligase_assumption),request.treatment.sealing.nicks.len(),request.selection.component,request.selection.reference_seed.strand,request.selection.reference_seed.index,e(&request.selection.sealing_context_revision),request.options.require_all_requested_transfers_applied,request.options.require_known_terminal_chemistry));
                for material in &request.treatment.preparation.recipe.materials {
                    match material {
                        crate::topo::MaterialRecipe::Prepared{declaration_id,assumption,strands,pairing}=>{
                            body.push_str(&format!("<li>Declared preparation {}: {}<ul>",e(declaration_id),e(assumption)));
                            for strand in strands {
                                body.push_str(&format!("<li>Prepared strand {}: declared 5′ {:?}, 3′ {}. Original source intervals (zero-based start; length; orientation): {}.</li>",e(&strand.id),strand.five_prime,e(&format!("{:?}",strand.three_prime)),strand.segments.iter().map(|segment|format!("{}: {}; {}; {}",e(&segment.source_role),segment.start,segment.length,if segment.reverse{"reverse"}else{"forward"})).collect::<Vec<_>>().join("; ")));
                            }
                            body.push_str(&format!("<li>{} explicitly declared antiparallel pairing spans; complete source-coordinate graph remains in matching JSON evidence.</li></ul></li>",pairing.len()));
                        },
                        crate::topo::MaterialRecipe::ExactPcr{material_id,template_role,forward,reverse,tail}=>body.push_str(&format!("<li>Actual exact PCR {} from role {}: ordered oligos 5′→3′ <code>{}</code> and <code>{}</code>; explicit tail {:?}.</li>",e(material_id),e(template_role),e(&forward.ordered_sequence()),e(&reverse.ordered_sequence()),tail)),
                    }
                }
                for transfer in &request.treatment.preparation.transfers {
                    body.push_str(&format!("<li>TOPO transfer: donor prepared strand {} to acceptor prepared strand {}, geometry {:?} (zero-based strand indexes).</li>",transfer.donor,transfer.acceptor,transfer.geometry));
                }
                for nick in &request.treatment.sealing.nicks {
                    body.push_str(&format!("<li>Ligase nick: strand {} call {} (3′) to strand {} call {} (5′), zero-based original source coordinates. Actual exposed chemistry/template eligibility is checked during reconstruction.</li>",nick.three.strand,nick.three.index,nick.five.strand,nick.five.index));
                }
                body.push_str("</ul></li>");
            },
            Operation::Topo { sources, request, selection } | Operation::CheckedTopo { sources, request, selection, .. } => {
                body.push_str(&format!("<li>Source-backed TOPO elementary-transfer graph with explicit ordered roles: {}. Preparation/transfers/limits <code>{}</code>; selection <code>{}</code>. Full physical evidence is retained; a paired-contour reference view does not establish nick repair, flap removal, supplier state, target acceptance or experimental success.<ul>",sources.iter().map(|s|format!("{} = {}",e(&s.role),source(&s.source))).collect::<Vec<_>>().join(", "),e(&serde_json::to_string(request).unwrap_or_default()),e(&serde_json::to_string(selection).unwrap_or_default())));
                if let Operation::CheckedTopo { policy, .. } = step {
                    body.push_str(&format!("<li>Final-only independent-target physical policy: <code>{}</code>. Acceptance is conditional on this declared graph and fully specified DNA. Rotation/reverse comparison is allowed only as explicitly listed; source-coordinate annotation correspondence does not establish chemical ancestry, graph repair or experimental success.</li>", e(&serde_json::to_string(policy).unwrap_or_default())));
                }
                for material in &request.recipe.materials {
                    match material {
                        crate::topo::MaterialRecipe::Prepared {declaration_id,assumption,..} => body.push_str(&format!("<li>Preparation {}: {}</li>",e(declaration_id),e(assumption))),
                        crate::topo::MaterialRecipe::ExactPcr {material_id,template_role,forward,reverse,tail} => body.push_str(&format!("<li>Fresh exact PCR {} from role {}. Actual ordered oligos 5′→3′: <code>{}</code> and <code>{}</code>. Optional separately declared 3′ tail: {:?}. Unique current geometry and ideal newly synthesized OH-ended strands are assumptions, not measured amplification.</li>",e(material_id),e(template_role),e(&forward.ordered_sequence()),e(&reverse.ordered_sequence()),tail)),
                    }
                }
                body.push_str("</ul></li>");
            },
            Operation::RecombineReference { sources, request, selection } => body.push_str(&format!("<li>Reference-only {:?} exchange of {} and {}. Pinned full profile/data, work/annotation limits: <code>{}</code>. Explicit component selection: <code>{}</code>. All reciprocal components are computed; selected reference DNA does not establish chemical state, efficacy, repair or experimental selection.</li>",request.reaction,source(&sources[0]),source(&sources[1]),e(&serde_json::to_string(request).unwrap_or_default()),e(&serde_json::to_string(selection).unwrap_or_default()))),
            Operation::PcrWithMismatches{template,forward,reverse,options}=>body.push_str(&format!("<li>Intentional-mismatch PCR of {}. Oligos 5′→3′: <code>{}</code> and <code>{}</code>. Binding/search policy: <code>{}</code>. Primer-derived changes are reconstructed; amplification is not established.</li>",source(template),e(&forward.ordered_sequence()),e(&reverse.ordered_sequence()),e(&serde_json::to_string(options).unwrap_or_default()))),
            Operation::GoldenGate{templates,options,candidate_index}=>body.push_str(&format!("<li>Golden Gate with {:?}, inputs {}. Explicit zero-based candidate {} under policy <code>{}</code>. Cleavage, orientation, cohesive ends and recognition-site removal are checked. Finite one-copy fragment search does not establish exhaustive mixture uniqueness; conditions and quantities are not specified by the molecular model.</li>",options.enzyme,templates.iter().map(source).collect::<Vec<_>>().join(", "),candidate_index,e(&serde_json::to_string(options).unwrap_or_default()))),
            Operation::SelectedGoldenGate{templates,options,candidate_index,directed_sequence_key,physical_sequence_key}=>body.push_str(&format!("<li>Key-bound Golden Gate with {:?}, inputs {}. Explicit zero-based candidate {} under policy <code>{}</code>. Cleavage, orientation, cohesive ends and recognition-site removal are checked. Finite one-copy fragment search does not establish exhaustive mixture uniqueness; conditions and quantities are not specified by the molecular model. Fresh enumeration must match directed key <code>{}</code> and physical key <code>{}</code>.</li>",options.enzyme,templates.iter().map(source).collect::<Vec<_>>().join(", "),candidate_index,e(&serde_json::to_string(options).unwrap_or_default()),e(directed_sequence_key),e(physical_sequence_key))),
            Operation::MultiSiteRecombineReference{sources,request,selection}=>body.push_str(&format!("<li>Declared multi-site {:?} reference exchange using sources {}. Profile <code>{}</code>; specificity order {:?}; select candidate {}, component {} from bound batch <code>{}</code>. Retain all reciprocal components. This is a reference design, not established chemical material, experimental selection or kit success.</li>",request.reaction.reaction,e(&sources.iter().map(|s|format!("{}: {}",s.role,source(&s.source))).collect::<Vec<_>>().join(", ")),e(&request.reaction.profile.id),request.reaction.specificities,selection.candidate_id,selection.component_id,e(&selection.batch_revision))),
            Operation::Pcr { template, forward, reverse } => body.push_str(&format!("<li>PCR of {} with {}: <code>{}</code> and {}: <code>{}</code>. Primers shown 5′→3′. The molecular model does not specify reaction conditions or quantities.</li>",source(template),e(&forward.name),e(&forward.ordered_sequence()),e(&reverse.name),e(&reverse.ordered_sequence()))),
            Operation::DeclareSourceTermini{template,left,right}=>body.push_str(&format!("<li>Declare original intact linear duplex source termini for {}: left <code>{}</code>, right <code>{}</code>. Explicit operator hypothesis only; not measured chemistry, synthesis, repair, recovered imported physical state or experimental success. Original sequence/features/citations are retained; methylation remains unknown.</li>",source(template),e(&serde_json::to_string(left).unwrap_or_default()),e(&serde_json::to_string(right).unwrap_or_default()))),
            Operation::DeclareSourceModifications{template,assumptions}=>body.push_str(&format!("<li>Declare closed Dam/Dcm/CpG source generation for {} under explicit assumptions <code>{}</code>. Bind the original input revision; these are not measured modifications. Descendant strand states are propagated without resynthesizing methylation.</li>",source(template),e(&serde_json::to_string(assumptions).unwrap_or_default()))),
            Operation::CheckedDigest{template,selection,options,fragment_index}=>body.push_str(&format!("<li>Checked positional-modification digest of {} with {}. Profile <code>{}</code>; zero-based fragment {}. Analysis limits: {} comparisons, {} sites, {} fragments, {} report bytes. Blocked sites remain intact; unknown/impaired or incomplete evidence withholds fully determined products. Reaction conditions and experimental success remain unspecified.</li>",source(template),e(&selection.enzymes.iter().map(|e|e.enzyme.name.clone()).collect::<Vec<_>>().join(", ")),e(&modified_duplex::profile_revision()),fragment_index,options.max_comparisons,options.max_sites,options.max_fragments,options.max_report_bytes)),
            Operation::CatalogDigest{template,selection,fragment_index}=>body.push_str(&format!("<li>Catalog-pinned complete digest of {} under exact selection <code>{}</code>; select zero-based fragment {} in top-reference cut ordering. Identity/name/geometry are checked against the embedded catalog. Complete cleavage and compatible conditions remain assumptions; supplier observations do not establish kinetics, methylation compatibility or experimental completion.</li>",source(template),e(&serde_json::to_string(selection).unwrap_or_default()),fragment_index)),
            Operation::Digest{template,enzymes,fragment_index}=>body.push_str(&format!("<li>Complete restriction digest of {} with {}. Select zero-based fragment index {} in ascending top-reference cut order. Cleavage assumes compatible conditions and no blocking methylation.</li>",source(template),e(&enzymes.join(", ")),fragment_index)),
            Operation::TreatEnds{template,treatment}=>body.push_str(&format!("<li>Terminal {:?} treatment of {}. Internal nick repair is not assumed; the molecular model does not specify reaction conditions or quantities.</li>",treatment,source(template))),
            Operation::Tail3Prime{template,base}=>body.push_str(&format!("<li>Append one {} nucleotide to both strand 3′ ends of {}. Requires reconstructed blunt linear DNA with accessible 3′-OH. Five-prime chemistry and internal nicks remain unchanged; no template ancestry or polymerase/kit activity is inferred.</li>",base.as_str(),source(template))),
            Operation::ProfiledOverlap{fragments,overlaps,circular,options}=>body.push_str(&format!("<li>Profiled {:?} ordered overlap assembly: {}. Overlaps {:?} bp; {} product. Explicit source orientation is retained. Method recommendations are advisory, not evidence of reagent compatibility, experimental completion or uniqueness. Policy <code>{}</code>.</li>",options.method,fragments.iter().map(|f|format!("{} ({})",source(&f.source),if f.reverse{"reverse complement"}else{"forward"})).collect::<Vec<_>>().join(" → "),overlaps,if *circular{"circular"}else{"linear"},e(&serde_json::to_string(options).unwrap_or_default()))),
            Operation::RepairEnds{template,method}=>body.push_str(&format!("<li>End repair {:?} of {}. Only explicitly modeled 3′ fill-in/exonuclease changes occur. Five-prime chemistry and internal nicks are retained; polymerase activity, reagents and experimental success are assumptions.</li>",method,source(template))),
            Operation::AnnealOligos{top,bottom,options,candidate_index}=>body.push_str(&format!("<li>Anneal declared 5′→3′ ssDNA inputs {} and {}.<br>Top ordered oligo 5′→3′: <code>{}</code> (5′ {:?}; 3′ {:?}).<br>Bottom ordered oligo 5′→3′: <code>{}</code> (5′ {:?}; 3′ {:?}).<br>Select zero-based candidate {} under policy <code>{}</code>. Exact antiparallel pairing and placement alternatives are reconstructed; no covalent ligation or implicit phosphorylation occurs.</li>",source(top),source(bottom),e(protocol_oligo_sequence(plan,&report,top,step_index).unwrap_or("Unresolved or unsupported bound oligo sequence")),options.top_chemistry.five_prime,options.top_chemistry.three_prime,e(protocol_oligo_sequence(plan,&report,bottom,step_index).unwrap_or("Unresolved or unsupported bound oligo sequence")),options.bottom_chemistry.five_prime,options.bottom_chemistry.three_prime,candidate_index,e(&serde_json::to_string(options).unwrap_or_default()))),
            Operation::Ligate{fragments,circular,allow_nicked,options}=>body.push_str(&format!("<li>Intended ordered ligation: {}. Product {}; retained nicks {}. Enumerate one-copy alternatives under explicit bounds <code>{}</code>; retain the intended path independently. No probabilities or experimental isolation are inferred. Full extractable outcome paths/sequences are retained below.</li>",fragments.iter().map(|f|format!("{} ({})",source(&f.source),if f.reverse{"reverse"}else{"forward"})).collect::<Vec<_>>().join(", "),if *circular{"circular"}else{"linear"},if *allow_nicked{"allowed"}else{"not allowed"},e(&serde_json::to_string(options).unwrap_or_default()))),
            Operation::Overlap { fragments, overlaps, circular } => body.push_str(&format!("<li>Ordered overlap assembly, in this order: {}. Junction overlaps {:?} bp; {} product. The molecular model does not specify reagent conditions or quantities.</li>",fragments.iter().map(source).collect::<Vec<_>>().join(", "),overlaps,if *circular {"circular"} else {"linear"})),
        }
    }
    body.push_str("</ol><h2>Reconstructed products</h2><ul>");
    for (index, product) in report.products.iter().enumerate() {
        body.push_str(&format!(
            "<li>Step {}: {} · {} bases · {:?}</li>",
            index + 1,
            e(&product.name),
            product.sequence.len(),
            product.topology
        ));
        if let Some(summary) = product
            .metadata
            .get("duplex_model")
            .and_then(|s| serde_json::from_str::<chemistry::DuplexSummary>(s).ok())
        {
            body.push_str(&format!(
                "<li><pre style=\"white-space:pre-wrap\">{}</pre></li>",
                e(&summary.describe())
            ));
        }
    }
    body.push_str("</ul><h2>Intended target map</h2>");
    match crate::svg::render_map(target) {
        Ok(svg) => body.push_str(&format!("<div class=\"protocol-map\">{svg}</div>")),
        Err(error) => body.push_str(&format!(
            "<p>Map unavailable: {}</p>",
            e(&error.to_string())
        )),
    }
    body.push_str("<h2>Verification findings</h2><ul>");
    for finding in &report.findings {
        body.push_str(&format!("<li>{}</li>", e(finding)));
    }
    body.push_str("</ul>");
    if !report.step_reports.is_empty() {
        body.push_str("<h2>TOPO physical step evidence</h2><ul>");
        for step in &report.step_reports {
            if let Some(topo) = &step.topo {
                let physical = &topo.projection.physical;
                let graph = &physical.product;
                let unpaired = graph
                    .unpaired
                    .iter()
                    .map(|run| run.end - run.start)
                    .sum::<usize>();
                let contour = match &topo.projection.outcome {
                    crate::topo::physical::projection::Outcome::Proven(view) => format!(
                        "Proven {:?} paired contour with {} retained nicks; known terminal chemistry {}; all requested transfers applied {}. This is no target-acceptance or repair claim.",
                        view.topology,
                        view.nicks.len(),
                        view.fully_known_terminal_chemistry,
                        view.all_requested_transfers_applied
                    ),
                    crate::topo::physical::projection::Outcome::Withheld(reason) => format!(
                        "Contour withheld: {reason:?}. No target-derived path, flap removal or nick repair was attempted."
                    ),
                };
                body.push_str(&format!("<li>Step {}: {} components, {} actual chains, {} graph nicks and {} unpaired nucleotides retained. {} The full bounded graph and source-call maps remain in the downloadable verification JSON.<ul>",step.step_index+1,graph.components.len(),graph.chains.len(),graph.nicks.len(),unpaired,e(&contour)));
                for assumption in &topo.preparation.assumptions {
                    body.push_str(&format!("<li>{}</li>", e(assumption)));
                }
                body.push_str("</ul></li>");
            }
        }
        body.push_str("</ul>");
    }
    if !report.reference_step_reports.is_empty() {
        body.push_str("<h2>Whole reference-network evidence</h2><ul>");
        for step in &report.reference_step_reports {
            let ReferenceStepEvidence::MultiSiteRecombination { selection, batch } = &step.evidence;
            body.push_str(&format!("<li>Step {}: {} original parents; {} candidate networks; {} reciprocal components retained. Complete within declared model: {}. Selected candidate {}, component {}. No target determined a pairing; reference agreement does not establish physical chemistry.<ul>",step.step_index+1,batch.sources.len(),batch.reference.candidates.len(),batch.documents.len(),batch.reference.complete,selection.candidate_id+1,selection.component_id+1));
            for note in batch
                .reference
                .warnings
                .iter()
                .chain(&batch.reference.limits)
                .chain(&batch.reference.unsupported)
                .chain(&batch.warnings)
            {
                body.push_str(&format!("<li>{}</li>", e(note)));
            }
            body.push_str("</ul></li>");
        }
        body.push_str("</ul>");
    }
    if !report.physical_processing_reports.is_empty() {
        body.push_str("<h2>Accounted physical processing evidence</h2><ul>");
        for step in &report.physical_processing_reports {
            match &step.evidence {
                PhysicalProcessingEvidence::LigationOutcomes(evidence) => {
                    body.push_str(&format!("<li>Step {} ligation outcomes: intended path retained; {} alternatives, {} checked paths, {} incompatible and {} unresolved. Enumeration complete {}, exhaustive within one-copy scope {}. Limits: {}. These are possible native paths, not probabilities, a complete mixture model or proof of isolation.<ol>",step.step_index+1,evidence.report.alternatives.len(),evidence.report.paths_checked,evidence.report.incompatible_paths,evidence.report.unresolved_paths,evidence.report.enumeration_complete,evidence.report.exhaustive_within_scope,e(&format!("{:?}",evidence.report.limits))));
                    for (i, product) in std::iter::once(&evidence.report.intended)
                        .chain(&evidence.report.alternatives)
                        .enumerate()
                    {
                        body.push_str(&format!("<li>{}: {} · {:?}, {} bases · material <code>{}</code>; {} retained nicks.</li>",if i==0{"Intended".into()}else{format!("Alternative {i}")},e(&product.path.iter().map(|p|format!("{} ({})",p.key,if p.reverse{"reverse"}else{"forward"})).collect::<Vec<_>>().join(" → ")),product.topology,product.sequence.len(),e(&product.modifications.state_revision),product.duplex.nicks.len()));
                    }
                    body.push_str(
                        "</ol><details><summary>Full extractable outcome FASTA</summary><pre>",
                    );
                    body.push_str(&e(&evidence.report.fasta().unwrap_or_else(|error| {
                        format!("Outcome extraction refused: {error}")
                    })));
                    body.push_str("</pre></details><details><summary>Full typed outcome paths, strands, modifications and bounds (JSON)</summary><pre>");
                    body.push_str(&e(&serde_json::to_string_pretty(evidence)
                        .unwrap_or_else(|error| format!("Outcome JSON refused: {error}"))));
                    body.push_str("</pre></details></li>");
                }
                PhysicalProcessingEvidence::ThreePrimeResection(geometry) => {
                    body.push_str(&format!("<li>Step {} actual3′ resection: bottom left {} nt / top right {} nt; retained top {} nt, bottom {} nt at offset {}. Exact copied/removed source spans and unchanged5′ groups remain in JSON evidence.</li>",step.step_index+1,geometry.cuts.left,geometry.cuts.right,geometry.product.top.len(),geometry.product.bottom.len(),geometry.product.bottom_start));
                }
                PhysicalProcessingEvidence::HomologyAnnealing(geometry) => {
                    body.push_str(&format!("<li>Step {} annealed diagnostic: {} actual surviving strands, {} junctions; paired gapless diagnostic {}; terminal search complete {}. No joining bond is sealed and no host repair is assumed.<ul>",step.step_index+1,geometry.segments.len(),geometry.junctions.len(),geometry.paired_diagnostic.is_some(),geometry.terminal_search_complete));
                    for j in &geometry.junctions {
                        body.push_str(&format!("<li>Sources {} → {}: {} bp homology; top discontinuity [{}, {}) 3′ {:?}/5′ {:?}; bottom [{}, {}) 5′ {:?}/3′ {:?}. Equal boundaries mean an unsealed nick; nonempty intervals remain unfilled gaps.</li>",j.left_source+1,j.right_source+1,j.overlap,j.top.start,j.top.end,j.top.three_prime,j.top.five_prime,j.bottom.start,j.bottom.end,j.bottom.five_prime,j.bottom.three_prime));
                    }
                    body.push_str("</ul></li>");
                }
                PhysicalProcessingEvidence::HomologyAssembly(geometry) => {
                    body.push_str(&format!("<li>Step {} physical homology assembly: {} junctions, {} resected strand spans; paired diagnostic product available {}; incomplete {}; terminal search complete {}; {} comparisons.<ul>",step.step_index+1,geometry.junctions.len(),geometry.resected.len(),geometry.product.is_some(),geometry.incomplete,geometry.terminal_search_complete,geometry.comparisons));
                    for junction in &geometry.junctions {
                        body.push_str(&format!("<li>Sources {} → {}: {} bp homology; gaps top/bottom {}/{}; seals top/bottom {}/{}; top nick at {} and bottom nick at {} (zero-based product boundaries).</li>",junction.left_source+1,junction.right_source+1,junction.overlap,junction.top_gap,junction.bottom_gap,junction.top_sealed,junction.bottom_sealed,junction.top_nick,junction.bottom_nick));
                    }
                    for note in geometry.warnings.iter().chain(&geometry.limits) {
                        body.push_str(&format!("<li>{}</li>", e(note)));
                    }
                    body.push_str("</ul></li>");
                }
                PhysicalProcessingEvidence::OverlapExtension(geometry) => {
                    body.push_str(&format!("<li>Step {} reciprocal extension: selected homology {} bp; retained terminal pairings {:?}; complete search {}; unique within terminal model {}; {} comparisons. Original source geometry revisions: {} / {}. Exact copied/synthesized strand mappings remain in the reconstruction report.</li>", step.step_index+1, geometry.options.overlap, geometry.terminal_pairings.overlaps, geometry.terminal_pairings.complete, geometry.terminal_pairings.unique, geometry.terminal_pairings.comparisons, e(&geometry.source_geometry_revisions[0]), e(&geometry.source_geometry_revisions[1])));
                }
                PhysicalProcessingEvidence::DerivedOrigin(evidence) => {
                    let material = &evidence.material;
                    body.push_str(&format!("<li>Step {}: replayed derivation event <code>{}</code>; material <code>{}</code>. {} Origin incomplete: {}; showing {} of {} findings{}.<ul>",step.step_index+1,e(&evidence.origin.product.event_id),e(&material.material_revision),e(&material.duplex.describe()),material.origin_incomplete,material.origin_findings.len(),material.origin_findings_total,if material.origin_findings_truncated {" (truncated)"} else {""}));
                    for finding in &material.origin_findings {
                        body.push_str(&format!("<li>{}</li>", e(finding)));
                    }
                    body.push_str("</ul></li>");
                }
                PhysicalProcessingEvidence::RunProduct(evidence) => {
                    body.push_str(&format!("<li>Step {}: saved run <code>{}</code>, product {}, actual material summary <code>{}</code>. Full originating report is retained with the reconstruction recipe.<ul>", step.step_index+1, e(&evidence.origin.selection.run_id), evidence.origin.selection.step_index+1, e(&evidence.material_revision)));
                    if let Ok(origin_report) = evidence.decode_report() {
                        body.push_str(&format!("<li>Origin assessment: {:?}; engine {}. Reuse does not upgrade incomplete findings.</li>", origin_report.status,e(&origin_report.engine_version)));
                        for finding in origin_report
                            .findings
                            .iter()
                            .chain(&origin_report.assumptions)
                        {
                            body.push_str(&format!("<li>{}</li>", e(finding)));
                        }
                    }
                    body.push_str("</ul></li>");
                }
                PhysicalProcessingEvidence::TopoActions(report) => {
                    let graph = &report.projection.processing.product;
                    body.push_str(&format!("<li>Step {}: {} original source calls remain across {} components, {} chains and {} unpaired calls. Faithful selected duplex available: {}. Excised material is retained, not discarded from evidence.<ul>",step.step_index+1,report.projection.processing.source_bases,graph.components.len(),graph.chains.len(),graph.unpaired.iter().map(|r|r.end-r.start).sum::<usize>(),report.tracked.is_some()));
                    for action in &report.projection.processing.steps {
                        body.push_str(&format!("<li>{}</li>", e(&format!("{action:?}"))));
                    }
                    for note in &report.bridge_notes {
                        body.push_str(&format!("<li>{}</li>", e(note)));
                    }
                    body.push_str("</ul></li>");
                }
            }
        }
        body.push_str("</ul>");
    }
    if !report.processed_step_reports.is_empty() {
        body.push_str("<h2>Processed TOPO physical step evidence</h2><ul>");
        for step in &report.processed_step_reports {
            let evidence = &step.evidence;
            let sealing = &evidence.projection.sealing;
            let graph = &sealing.product;
            body.push_str(&format!("<li>Step {}: original graph {} nicks; processed graph {} nicks, {} components, {} source calls, {} released activations. Faithful duplex bridge available: {}. Full parent/child graphs, per-seal chemical outcomes and source maps remain in typed verification evidence.<ul>",step.step_index+1,sealing.parent.product.nicks.len(),graph.nicks.len(),graph.components.len(),sealing.source_bases,graph.released_activations.len(),evidence.tracked.is_some()));
            for outcome in &sealing.steps {
                body.push_str(&format!(
                    "<li>Nick seal {}: {:?}; {}. Declared chemistry 3′ {}, 5′ {:?}.</li>",
                    outcome.step + 1,
                    outcome.outcome,
                    e(&outcome.reason),
                    e(&format!("{:?}", outcome.three_prime)),
                    outcome.five_prime
                ));
            }
            for note in &evidence.bridge_notes {
                body.push_str(&format!("<li>{}</li>", e(note)));
            }
            if let Some(annotation) = &evidence.annotations {
                body.push_str(&format!("<li>{} source-correspondence annotations retained; {} warnings. Annotation correspondence establishes neither inherited modifications nor biological validity.</li>",annotation.features.len(),annotation.warnings.len()));
            }
            body.push_str("</ul></li>");
        }
        body.push_str("</ul>");
    }
    if !report.annotation_reports.is_empty() {
        body.push_str("<h2>TOPO annotation correspondence</h2><p>Projected source coordinates establish neither inherited chemical ancestry nor biological feature validity. Unmapped annotations remain original source evidence.</p><ul>");
        for annotation in &report.annotation_reports {
            let mut offset = 0usize;
            for role in &annotation.evidence.sources {
                let end = offset.saturating_add(role.projected);
                // The projector emits each role's features in source-role order.
                let reviewed = annotation
                    .evidence
                    .features
                    .get(offset..end)
                    .unwrap_or(&[])
                    .iter()
                    .filter(|feature| feature.needs_review)
                    .count();
                body.push_str(&format!("<li>Step {}, role {}: {} projected ({} require review), {} original annotations retained unmapped, {} outside selected correspondence. Projection supported: {}.</li>",annotation.step_index+1,e(&role.role),role.projected,reviewed,role.unmapped.len(),role.dropped,role.projection_supported));
                offset = end;
            }
        }
        body.push_str("</ul>");
    }
    body.push_str("<h2>Assumptions and limits</h2><ul>");
    for assumption in &report.assumptions {
        body.push_str(&format!("<li>{}</li>", e(assumption)));
    }
    body.push_str(&format!("</ul><h2>Provenance</h2><p>Engine {}<br>Plan <code>{}</code><br>Target <code>{}</code></p>",e(&report.engine_version),e(&report.plan_revision),e(&report.target_revision)));
    if let Some(appendix) = checked.appendix() {
        body.push_str("<section><h2>Declared reaction quantities</h2><p>Arithmetic for the declared subset is separate from the molecular verification above. Unlisted conditions and reagents remain unspecified.</p>");
        for line in &appendix.paragraphs {
            body.push_str(&format!("<p>{}</p>", e(line)));
        }
        for table in &appendix.tables {
            body.push_str(&format!("<h3>{}</h3><table><thead><tr>", e(&table.heading)));
            for column in &table.columns {
                body.push_str(&format!("<th>{}</th>", e(column)));
            }
            body.push_str("</tr></thead><tbody>");
            for row in &table.rows {
                body.push_str("<tr>");
                for cell in row {
                    body.push_str(&format!("<td>{}</td>", e(cell)));
                }
                body.push_str("</tr>");
            }
            body.push_str("</tbody></table>");
        }
        body.push_str("</section>");
    }
    format!(
        "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width\"><title>{}</title><style>body{{font:16px/1.6 system-ui;max-width:850px;margin:40px auto;padding:0 24px;color:#18302d}}code{{overflow-wrap:anywhere}}.protocol-map>svg{{display:block;max-width:100%;height:auto}}h1,h2{{line-height:1.2}}li{{margin:12px 0}}@media print{{body{{margin:0;max-width:none}}h2{{break-after:avoid}}li{{break-inside:avoid}}}}</style>{body}</html>",
        e(&plan.name)
    )
}
