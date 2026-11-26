//! Versioned bounded binary messages. Schema changes require a new version.
use crate::{
    analysis::{CutSite, Enzyme, Match},
    formats::MAX_IMPORT_BYTES,
    model::*,
};
use serde::{Deserialize, Serialize};

pub const WORKER_SCHEMA_VERSION: u16 = 82;
pub const CANDIDATE_SCHEMA_VERSION: u16 = 18;

fn current_version(magic: &[u8; 4]) -> u16 {
    match magic {
        b"CGQ\0" | b"CGR\0" => WORKER_SCHEMA_VERSION,
        b"CGP\0" => CANDIDATE_SCHEMA_VERSION,
        _ => 1,
    }
}
const HEADER: usize = 42;

pub const MAX_WORKSPACE_MESSAGE_BYTES: usize = 128 * 1024 * 1024;
struct BoundedPayload {
    bytes: Vec<u8>,
    limit: usize,
}
impl postcard::ser_flavors::Flavor for BoundedPayload {
    type Output = Vec<u8>;
    fn try_push(&mut self, byte: u8) -> Result<(), postcard::Error> {
        if self.bytes.len() >= self.limit {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.bytes.push(byte);
        Ok(())
    }
    fn try_extend(&mut self, bytes: &[u8]) -> Result<(), postcard::Error> {
        if bytes.len() > self.limit.saturating_sub(self.bytes.len()) {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.bytes.extend_from_slice(bytes);
        Ok(())
    }
    fn finalize(self) -> Result<Self::Output, postcard::Error> {
        Ok(self.bytes)
    }
}
pub(crate) fn encode<T: Serialize>(magic: &[u8; 4], value: &T) -> Result<Vec<u8>, Error> {
    encode_with_limit(magic, value, MAX_IMPORT_BYTES)
}
fn encode_with_limit<T: Serialize>(
    magic: &[u8; 4],
    value: &T,
    maximum: usize,
) -> Result<Vec<u8>, Error> {
    let payload = postcard::serialize_with_flavor(
        value,
        BoundedPayload {
            bytes: Vec::new(),
            limit: maximum,
        },
    )
    .map_err(|error| {
        if matches!(error, postcard::Error::SerializeBufferFull) {
            Error::TooLarge
        } else {
            Error::Format(error.to_string())
        }
    })?;
    let mut bytes = Vec::with_capacity(HEADER + payload.len());
    bytes.extend_from_slice(magic);
    bytes.extend_from_slice(&current_version(magic).to_le_bytes());
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
    bytes.extend_from_slice(&payload);
    Ok(bytes)
}
fn envelope_info(magic: &[u8; 4], bytes: &[u8], maximum: usize) -> Result<(u16, usize), Error> {
    if bytes.len() < HEADER || &bytes[..4] != magic {
        return Err(Error::Format(
            "not a supported Collagene binary file/message".into(),
        ));
    }
    let version = u16::from_le_bytes([bytes[4], bytes[5]]);
    if version != current_version(magic) && !(version >= 1 && version < current_version(magic)) {
        return Err(Error::Unsupported(format!(
            "binary schema version {version}"
        )));
    }
    let len = u32::from_le_bytes(bytes[6..10].try_into().unwrap()) as usize;
    if len > maximum {
        return Err(Error::TooLarge);
    }
    if HEADER.checked_add(len) != Some(bytes.len()) {
        return Err(Error::Format("binary length mismatch".into()));
    }
    Ok((version, len))
}
pub(crate) fn decode<T: for<'de> Deserialize<'de>>(
    magic: &[u8; 4],
    bytes: &[u8],
) -> Result<T, Error> {
    decode_with_limit(magic, bytes, MAX_IMPORT_BYTES)
}
fn decode_with_limit<T: for<'de> Deserialize<'de>>(
    magic: &[u8; 4],
    bytes: &[u8],
    maximum: usize,
) -> Result<T, Error> {
    envelope_info(magic, bytes, maximum)?;
    if blake3::hash(&bytes[HEADER..]).as_bytes() != &bytes[10..42] {
        return Err(Error::Format("binary checksum mismatch".into()));
    }
    let (result, rest) =
        postcard::take_from_bytes(&bytes[HEADER..]).map_err(|e| Error::Format(e.to_string()))?;
    if !rest.is_empty() {
        return Err(Error::Format("unexpected trailing payload".into()));
    }
    Ok(result)
}
// Postcard structs encode their fields consecutively; the borrowed prefix probes only
// id/revision/variant and cannot allocate the untrusted large variant's payload.
fn worker_limit(
    magic: &[u8; 4],
    bytes: &[u8],
    allowed_large_variants: &[(u32, u16)],
) -> Result<usize, Error> {
    let (version, len) = envelope_info(magic, bytes, MAX_WORKSPACE_MESSAGE_BYTES)?;
    if len <= MAX_IMPORT_BYTES {
        return Ok(MAX_IMPORT_BYTES);
    }
    if version < 6 {
        return Err(Error::TooLarge);
    }
    let ((_, revision, variant), _) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
            .map_err(|e| Error::Format(e.to_string()))?;
    if revision.len() > 128
        || !allowed_large_variants
            .iter()
            .any(|&(tag, minimum)| tag == variant && version >= minimum)
    {
        return Err(Error::TooLarge);
    }
    Ok(MAX_WORKSPACE_MESSAGE_BYTES)
}

pub fn save_document(mol: &Molecule) -> Result<Vec<u8>, Error> {
    mol.validate()?;
    encode(b"CG2\0", mol)
}
pub fn load_document(bytes: &[u8]) -> Result<Molecule, Error> {
    let mol: Molecule = decode(b"CG2\0", bytes)?;
    mol.validate()?;
    Ok(mol)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AnalysisJob {
    Search {
        molecule: Molecule,
        query: String,
    },
    Restriction {
        molecule: Molecule,
        enzyme: Enzyme,
    },
    Pcr {
        molecule: Molecule,
        forward: crate::analysis::Primer,
        reverse: crate::analysis::Primer,
    },
    Verify {
        target: Molecule,
        plan: crate::verification::DesignPlan,
    },
    DuplexDigest {
        molecule: Molecule,
        enzymes: Vec<Enzyme>,
    },
    Translate {
        molecule: Molecule,
        options: crate::translation::TranslationOptions,
    },
    TranslateCds {
        molecule: Molecule,
        feature_index: usize,
    },
    FindOrfs {
        molecule: Molecule,
        options: crate::translation::OrfOptions,
    },
    Align {
        molecule: Molecule,
        query: Molecule,
        options: crate::alignment::AlignmentOptions,
    },
    PcrDocuments {
        molecule: Molecule,
        forward: crate::analysis::Primer,
        reverse: crate::analysis::Primer,
    },
    AnalyzePrimer {
        molecule: Molecule,
        primer: crate::analysis::Primer,
        conditions: crate::primers::ThermoConditions,
        binding_options: crate::primers::BindingOptions,
    },
    DesignPrimers {
        molecule: Molecule,
        options: crate::primers::DesignOptions,
    },
    ExportProtocolPdf {
        target: Molecule,
        plan: crate::verification::DesignPlan,
    },
    GoldenGate {
        molecule: Molecule,
        templates: Vec<Molecule>,
        options: crate::golden_gate::GoldenGateOptions,
    },
    MismatchPcr {
        molecule: Molecule,
        forward: crate::analysis::Primer,
        reverse: crate::analysis::Primer,
        options: crate::mutagenesis::PcrOptions,
    },
    InversePcr {
        molecule: Molecule,
        target: Molecule,
        options: crate::mutagenesis::InversePcrOptions,
    },
    PredictGel {
        molecule: Molecule,
        lanes: Vec<crate::gel::LaneInput>,
        calibration: crate::gel::Calibration,
    },
    DiagnosticDigests {
        molecule: Molecule,
        candidates: Vec<Molecule>,
        options: crate::gel::DiagnosticOptions,
    },
    MismatchPcrPreview {
        molecule: Molecule,
        forward: crate::analysis::Primer,
        reverse: crate::analysis::Primer,
        options: crate::mutagenesis::PcrOptions,
    },
    SaveWorkspace {
        workspace: crate::workspace::Workspace,
    },
    LoadWorkspace {
        bytes: Vec<u8>,
    },
    ValidateTraces {
        molecule: Molecule,
        reads: Vec<Molecule>,
        options: crate::trace_validation::TraceValidationOptions,
    },
    ScanFeatureLibrary {
        molecule: Molecule,
        library: crate::feature_library::FeatureLibrary,
        options: crate::feature_library::ScanOptions,
    },
    ApplyFeatureLibrary {
        molecule: Molecule,
        library: crate::feature_library::FeatureLibrary,
        report: crate::feature_library::ScanReport,
        selected: Vec<String>,
    },
    ImportAnnotations {
        molecule: Molecule,
        text: String,
        format: crate::annotation_formats::AnnotationFormat,
        seqid: String,
    },
    ApplyAnnotations {
        molecule: Molecule,
        report: crate::annotation_formats::AnnotationImportReport,
        selected: Vec<usize>,
    },
    ExportAnnotations {
        molecule: Molecule,
        format: crate::annotation_formats::AnnotationFormat,
        seqid: String,
    },
    RepairPreview {
        molecule: Molecule,
        enzymes: Vec<String>,
        fragment_index: usize,
        method: crate::chemistry::EndRepairMethod,
    },
    AnnealPreview {
        molecule: Molecule,
        bottom: Molecule,
        options: crate::chemistry::AnnealOptions,
    },
    AnnealedDocument {
        molecule: Molecule,
        bottom: Molecule,
        preview: crate::verification::AnnealPreview,
        candidate_index: usize,
    },
    ViewTracks {
        molecule: Molecule,
        window: crate::view_tracks::TrackWindow,
        options: crate::view_tracks::TrackOptions,
    },
    SearchDocument {
        molecule: Molecule,
        options: crate::search::SearchOptions,
    },
    PrepareTracks {
        molecule: Molecule,
    },
    CachedViewTracks {
        window: crate::view_tracks::TrackWindow,
        options: crate::view_tracks::TrackOptions,
    },
    ImportAlignment {
        bytes: Vec<u8>,
        format: crate::alignment_formats::AlignmentFormat,
        alphabet: Alphabet,
        name: String,
    },
    ExportAlignment {
        document: crate::alignment_formats::AlignmentDocument,
        format: crate::alignment_formats::AlignmentFormat,
    },
    ExtractAlignmentRows {
        document: crate::alignment_formats::AlignmentDocument,
        row_ids: Vec<String>,
    },
    SaveAlignment {
        document: crate::alignment_formats::AlignmentDocument,
    },
    LoadAlignment {
        bytes: Vec<u8>,
    },
    ProteinProperties {
        molecule: Molecule,
        options: crate::protein_properties::ProteinOptions,
    },
    TailPreview {
        molecule: Molecule,
        recipe: crate::verification::TailSubstrate,
        base: crate::tailing::TailingBase,
    },
    CheckOverlapProfile {
        molecules: Vec<Molecule>,
        reverse: Vec<bool>,
        overlaps: Vec<usize>,
        circular: bool,
        options: crate::overlap_profiles::ProfileOptions,
    },
    RecordCloningRun {
        workspace: crate::workspace::Workspace,
        plan: crate::verification::DesignPlan,
        target_binding: crate::run_lineage::EventBinding,
        input_bindings: Vec<crate::run_lineage::InputEventBinding>,
        allow_incomplete: bool,
    },
    ReplayCloningRun {
        workspace: crate::workspace::Workspace,
        run_id: String,
    },
    DesignCodons {
        molecule: Molecule,
        options: crate::codon_design::DesignOptions,
    },
    CheckCodonDesign {
        source: Molecule,
        candidate: Molecule,
        options: crate::codon_design::DesignOptions,
    },
    DotPlot {
        molecule: Molecule,
        query: Molecule,
        options: crate::dotplot::DotPlotOptions,
    },
    AssembleContigs {
        reads: Vec<crate::contig::ReadInput>,
        options: crate::contig::ContigOptions,
    },
    AssessMethylation {
        molecule: Molecule,
        options: crate::methylation::MethylationOptions,
    },
    DesignCrisprGuides {
        molecule: Molecule,
        options: crate::crispr::GuideOptions,
    },
    SearchCrisprOffTargets {
        query: crate::crispr::GuideQuery,
        references: Vec<crate::crispr::ReferenceInput>,
        options: crate::crispr::OffTargetOptions,
    },
    DesignEditTemplate {
        source: Molecule,
        target: Molecule,
        options: crate::edit_template::EditTemplateOptions,
        guide: Option<crate::edit_template::GuideBinding>,
    },
    CheckEditTemplate {
        source: Molecule,
        target: Molecule,
        options: crate::edit_template::EditTemplateOptions,
        guide: Option<crate::edit_template::GuideBinding>,
        donor: String,
    },
    BeginGenomeFasta {
        build: crate::genome::BuildInfo,
        options: crate::genome::IngestOptions,
        expected_source_bytes: u64,
    },
    PushGenomeFasta {
        import_id: String,
        offset: u64,
        bytes: Vec<u8>,
    },
    FinishGenomeFasta {
        import_id: String,
    },
    BeginGenomeNative {
        header: Vec<u8>,
        tail: Vec<u8>,
        file_length: u64,
        manifest: Vec<u8>,
    },
    PushGenomeNativeBlock {
        import_id: String,
        block: crate::genome::EncodedBlock,
    },
    FinishGenomeNative {
        import_id: String,
    },
    AbortGenomeImport {
        import_id: String,
    },
    RegisterGenome {
        manifest: crate::genome::GenomeManifest,
    },
    GenomeIntervalGuides {
        genome: crate::genome_session::GenomeRef,
        options: crate::genome_crispr::GenomeGuideOptions,
        blocks: Vec<crate::genome::EncodedBlock>,
    },
    StartGenomeSearch {
        genome: crate::genome_session::GenomeRef,
        query: crate::crispr::GuideQuery,
        options: crate::genome_crispr::GenomeSearchOptions,
    },
    AdvanceGenomeSearch {
        job_id: String,
        blocks: Vec<crate::genome::EncodedBlock>,
    },
    PageGenomeSearch {
        job_id: String,
        offset: usize,
        limit: usize,
    },
    CancelGenomeJob {
        id: String,
    },
    ReadGenomeInterval {
        genome: crate::genome_session::GenomeRef,
        options: crate::genome_crispr::GenomeIntervalReadOptions,
        blocks: Vec<crate::genome::EncodedBlock>,
    },
    InspectRestriction {
        source_role: String,
        molecule: Molecule,
        enzymes: Vec<Enzyme>,
        options: crate::restriction_inspection::InspectionOptions,
    },
    MaterializeDigestFragment {
        source_role: String,
        molecule: Molecule,
        enzymes: Vec<Enzyme>,
        options: crate::restriction_inspection::InspectionOptions,
        fragment_index: usize,
    },
    InstallDigestFragment {
        workspace: crate::workspace::Workspace,
        source: crate::run_lineage::EventBinding,
        enzymes: Vec<Enzyme>,
        options: crate::restriction_inspection::InspectionOptions,
        fragment_index: usize,
    },
    InspectCatalogRestriction {
        selection: crate::catalog_digest::CatalogSelection,
        source_role: String,
        molecule: Molecule,
        options: crate::restriction_inspection::InspectionOptions,
    },
    MaterializeCatalogDigestFragment {
        selection: crate::catalog_digest::CatalogSelection,
        source_role: String,
        molecule: Molecule,
        options: crate::restriction_inspection::InspectionOptions,
        fragment_index: usize,
    },
    InstallCatalogDigestFragment {
        selection: crate::catalog_digest::CatalogSelection,
        workspace: crate::workspace::Workspace,
        source: crate::run_lineage::EventBinding,
        options: crate::restriction_inspection::InspectionOptions,
        fragment_index: usize,
    },
    InspectRecombination {
        inputs: Vec<crate::recombination::reference::Input>,
        request: crate::recombination::Request,
    },
    CheckRecombinationTarget {
        inputs: Vec<crate::recombination::reference::Input>,
        request: crate::recombination::Request,
        target: Molecule,
        policy: crate::recombination::reference::TargetPolicy,
    },
    SelectRecombinationComponent {
        inputs: Vec<crate::recombination::reference::Input>,
        request: crate::recombination::Request,
        selection: crate::recombination::Selection,
        target: Molecule,
        policy: crate::recombination::reference::TargetPolicy,
    },
    /// Recognition scan only; appended to preserve every older wire index.
    ScanRestriction {
        source_role: String,
        molecule: Molecule,
        enzymes: Vec<Enzyme>,
        options: crate::restriction_inspection::InspectionOptions,
    },
    PrepareReadAlignments {
        bytes: Vec<u8>,
        format: crate::read_alignment_session::FileFormat,
    },
    ReadAlignmentPage {
        source_hash: String,
        start: usize,
        count: usize,
    },
    CheckReadAlignment {
        source_hash: String,
        index: usize,
        reference: Molecule,
        options: crate::read_alignment_reference::Options,
    },
    ExportReadAlignments {
        source_hash: String,
        format: crate::read_alignment_session::FileFormat,
        max_output_bytes: usize,
    },
    DisposeReadAlignments {
        source_hash: String,
    },
    SimulateTopo {
        inputs: Vec<crate::topo::Input>,
        request: crate::topo::Request,
    },
    InspectTopoOrientations {
        inputs: Vec<crate::topo::Input>,
        request: crate::topo::Request,
        vector: [usize; 2],
        insert: [usize; 2],
        directional: bool,
    },
    ProjectTopoProduct {
        inputs: Vec<crate::topo::Input>,
        request: crate::topo::Request,
        selection: crate::topo::Selection,
    },
    CheckTopoTarget {
        inputs: Vec<crate::topo::Input>,
        request: crate::topo::Request,
        selection: crate::topo::Selection,
        target: Molecule,
        policy: crate::topo::AcceptancePolicy,
    },
    ProjectAnnotatedTopo {
        inputs: Vec<crate::topo::Input>,
        request: crate::topo::Request,
        selection: crate::topo::Selection,
    },
    CheckAnnotatedTopoTarget {
        inputs: Vec<crate::topo::Input>,
        request: crate::topo::Request,
        selection: crate::topo::Selection,
        target: Molecule,
        policy: crate::topo::AcceptancePolicy,
    },
    PrepareVariants {
        bytes: Vec<u8>,
        format: crate::variants_session::FileFormat,
    },
    VariantPage {
        source_hash: String,
        start: usize,
        count: usize,
    },
    CheckVariantReference {
        source_hash: String,
        reference: Molecule,
        options: crate::variant_reference::Options,
    },
    PreviewVariantAllele {
        source_hash: String,
        reference: Molecule,
        options: crate::variant_reference::Options,
    },
    ExportVariants {
        source_hash: String,
        format: crate::variants_session::ExportFormat,
        max_output_bytes: usize,
    },
    DisposeVariants {
        source_hash: String,
    },
    PreviewOverlapRecipe {
        inputs: Vec<crate::verification::PlanInput>,
        target: Molecule,
        recipe: crate::assembly_recipe::Recipe,
    },
    PreviewLigationRecipe {
        inputs: Vec<crate::verification::PlanInput>,
        target: Molecule,
        recipe: crate::ligation_recipe::Recipe,
    },
    InspectAssemblyJunctions {
        plan: crate::verification::DesignPlan,
        target: Molecule,
        options: crate::assembly_junctions::Options,
    },
    DesignOverlapPrimers {
        inputs: Vec<crate::verification::PlanInput>,
        target: Molecule,
        request: crate::overlap_primer_design::Request,
    },
    BuildAssemblyIntent {
        inputs: Vec<crate::verification::PlanInput>,
        request: crate::assembly_intent::Request,
    },
    SelectDigestFragment {
        selection: crate::catalog_digest::CatalogSelection,
        source_role: String,
        molecule: Molecule,
        options: crate::restriction_inspection::InspectionOptions,
        expected_batch_revision: String,
        choice: crate::digest_cut_selection::Choice,
    },
    ExportPlanOligos {
        plan: crate::verification::DesignPlan,
        target: Molecule,
        options: crate::plan_oligos::Options,
        delimiter: crate::oligo_formats::Delimiter,
    },
    ExportProtocolHtml {
        plan: crate::verification::DesignPlan,
        target: Molecule,
    },
    FinalizeCloningOutputs {
        workspace: crate::workspace::Workspace,
        plan: crate::verification::DesignPlan,
        target: crate::run_lineage::EventBinding,
        inputs: Vec<crate::run_lineage::InputEventBinding>,
        options: crate::cloning_outputs::Options,
    },
    DesignGoldenGatePrimers {
        inputs: Vec<crate::verification::PlanInput>,
        target: Molecule,
        request: crate::golden_gate_primer_design::Request,
    },
    BuildMutationIntent {
        molecule: Molecule,
        request: crate::mutation_intent::Request,
    },
    CaptureReadCheck {
        workspace: crate::workspace::Workspace,
        reference: crate::run_lineage::EventBinding,
        reads: Vec<crate::run_lineage::EventBinding>,
        options: crate::trace_validation::TraceValidationOptions,
    },
    ReplayReadCheck {
        workspace: crate::workspace::Workspace,
        attachment_hash: String,
    },
    PreviewReadReview {
        workspace: crate::workspace::Workspace,
        request: crate::read_review::Request,
    },
    CommitReadReview {
        workspace: crate::workspace::Workspace,
        request: crate::read_review::Request,
    },
    BuildEditedAssemblyIntent {
        inputs: Vec<crate::verification::PlanInput>,
        base: crate::assembly_intent::Request,
        edits: Vec<crate::assembly_junction_edit::Edit>,
    },
    DesignEditedAssemblyPrimers {
        inputs: Vec<crate::verification::PlanInput>,
        request: crate::assembly_junction_edit::Request,
    },
    DesignSelectiveAssemblyPrimers {
        inputs: Vec<crate::verification::PlanInput>,
        request: crate::selective_junction_edit::Request,
    },
    DesignGoldenGateJunctionPrimers {
        inputs: Vec<crate::verification::PlanInput>,
        request: crate::golden_gate_junction_edit::Request,
    },
    InspectPreparedGoldenGatePool {
        inputs: Vec<crate::verification::PlanInput>,
        target: Molecule,
        request: crate::prepared_golden_gate::Request,
    },
    SelectPreparedGoldenGateCandidate {
        report: Box<crate::prepared_golden_gate::Report>,
        selection: crate::prepared_golden_gate::Selection,
    },
    InspectGoldenGateReaction {
        inputs: Vec<crate::verification::PlanInput>,
        target: Molecule,
        request: crate::golden_gate_reaction::Request,
    },
    SelectGoldenGateReaction {
        report: Box<crate::golden_gate_reaction::Report>,
        selection: crate::prepared_golden_gate::Selection,
    },
    PreparePcrOutputs {
        request: crate::pcr_outputs::Request,
    },
    PrepareMismatchPcrOutputs {
        request: crate::mismatch_pcr_outputs::Request,
    },
    CompareSequences {
        reference: Molecule,
        query: Molecule,
        options: crate::pairwise_review::Options,
    },
    CopyReferenceFragment {
        request: crate::reference_transfer::RequestCopy,
    },
    PasteReferenceFragment {
        request: crate::reference_transfer::RequestPaste,
    },
    TranslateReference {
        request: crate::translation_reference::Request,
    },
    PrepareTranslatedReference {
        request: crate::translation_reference::Request,
        name: String,
    },
    PredictPinnedGel {
        molecule: Molecule,
        lanes: Vec<crate::gel::LaneInput>,
        calibration: crate::gel::Calibration,
    },
    ExportPairwiseComparison {
        record: crate::pairwise_history::Record,
        reference_name: String,
        query_name: String,
        fresh: bool,
        format: crate::pairwise_export::Format,
    },
    SearchPairwiseComparison {
        report: crate::pairwise_review::Report,
        pattern: String,
        options: crate::pairwise_search::Options,
    },
    PreparePairwiseContext {
        report: crate::pairwise_review::Report,
        reference: Molecule,
        query: Molecule,
        translate: bool,
    },
    BeginCombinatorialBatch {
        workspace: crate::workspace::Workspace,
        matrix: crate::combinatorial_batch::Request,
    },
    AdvanceCombinatorialBatch {
        context: String,
        index: usize,
    },
    CaptureCombinatorialBatch {
        workspace: crate::workspace::Workspace,
        current_context: String,
    },
    BeginCombinatorialRelease {
        workspace: crate::workspace::Workspace,
        request: crate::combinatorial_release::Request,
    },
    AdvanceCombinatorialRelease {
        context: String,
        index: usize,
    },
    InspectCombinatorialRelease {
        context: String,
        index: usize,
    },
    FinishCombinatorialRelease {
        workspace: crate::workspace::Workspace,
        current_context: String,
    },
    ExportCombinatorialOutputs {
        workspace: crate::workspace::Workspace,
        manifest_hash: String,
        options: crate::combinatorial_export::Options,
    },
    SearchGoldenGateFusions {
        inputs: Vec<crate::verification::PlanInput>,
        request: crate::golden_gate_fusion_design::Request,
    },
    MaterializeGoldenGateFusions {
        inputs: Vec<crate::verification::PlanInput>,
        request: crate::golden_gate_fusion_design::Request,
        selection: crate::golden_gate_fusion_design::Selection,
    },
    ScoreLigationFidelity {
        dataset: crate::ligation_fidelity::Dataset,
        ends: Vec<String>,
    },
    InspectDomestication {
        source: Molecule,
        request: crate::domestication::Request,
    },
    ProposeDomestication {
        source: Molecule,
        request: crate::domestication::Request,
    },
    CheckDomestication {
        source: Molecule,
        request: crate::domestication::Request,
        selection: crate::domestication::Selection,
    },
    PrepareDomestication {
        input: crate::verification::PlanInput,
        request: crate::domestication_preparation::Request,
    },
    ListLigationProfiles,
    ScorePublishedLigation {
        profile: crate::ligation_profiles::ProfileRef,
        ends: Vec<String>,
        physical_competitors: Vec<String>,
    },
    RankGoldenGateFusions {
        inputs: Vec<crate::verification::PlanInput>,
        request: crate::golden_gate_fusion_rank::Request,
    },
    ScoreLigationCompetition {
        dataset: crate::ligation_fidelity::Dataset,
        ends: Vec<String>,
        physical_competitors: Vec<String>,
    },
    SealTopoNicks {
        inputs: Vec<crate::topo::Input>,
        request: crate::topo_nick_sealing::Request,
    },
    ProjectProcessedTopo {
        inputs: Vec<crate::topo::Input>,
        request: crate::topo_processed_product::Request,
    },
    LoadProcessedRunEvidence {
        workspace: crate::workspace::Workspace,
        run_id: String,
    },
    ProcessTopoActions {
        inputs: Vec<crate::topo::Input>,
        request: crate::topo_actions::Request,
    },
    InspectTopoActions {
        inputs: Vec<crate::topo::Input>,
        request: crate::topo_actions::ActionsRequest,
    },
    SimulateMultiSiteRecombination {
        inputs: Vec<crate::recombination::reference::Input>,
        request: crate::multisite_recombination::Request,
    },
    ProjectMultiSiteRecombination {
        inputs: Vec<crate::recombination::reference::Input>,
        request: crate::multisite_recombination::Request,
        selection: crate::multisite_recombination::Selection,
    },
    CheckMultiSiteRecombination {
        inputs: Vec<crate::recombination::reference::Input>,
        request: crate::multisite_recombination::Request,
        selection: crate::multisite_recombination::Selection,
        target: Molecule,
    },
    InspectCombinatorialMultiSite {
        workspace: crate::workspace::Workspace,
        design: crate::combinatorial_design::Design,
        construct: crate::combinatorial_design::SelectedConstruct,
    },
    PreviewTailingRecipe {
        inputs: Vec<crate::verification::PlanInput>,
        target: Molecule,
        recipe: crate::tailing_recipe::Recipe,
    },
    InspectPcrDigest {
        molecule: Molecule,
        request: crate::pcr_digest_inspection::Request,
    },
    InspectRunMaterial {
        workspace: crate::workspace::Workspace,
        event: crate::run_lineage::EventBinding,
        selection: Option<crate::run_material::Selection>,
    },
    InspectRunPcrDigest {
        molecule: Molecule,
        origin: crate::run_material::Capsule,
        request: crate::pcr_digest_inspection::Request,
    },
    InspectRunCatalogRestriction {
        source_role: String,
        molecule: Molecule,
        origin: crate::run_material::Capsule,
        selection: crate::catalog_digest::CatalogSelection,
        options: crate::restriction_inspection::InspectionOptions,
    },
    BuildSubstitutedAssemblyIntent {
        inputs: Vec<crate::verification::PlanInput>,
        request: crate::assembly_substitution_intent::Request,
    },
    InspectMismatchPcrDigest {
        molecule: Molecule,
        origin: Option<crate::run_material::Capsule>,
        request: crate::mismatch_pcr_digest_inspection::Request,
    },
    BuildContributionEditIntent {
        inputs: Vec<crate::verification::PlanInput>,
        request: crate::assembly_edit_intent::Request,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Request {
    pub id: u32,
    pub revision: String,
    pub job: AnalysisJob,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AnalysisResult {
    Search(Vec<Match>),
    Restriction(Vec<CutSite>),
    Error(String),
    Pcr(Vec<crate::analysis::PcrProduct>),
    Verification(crate::verification::VerificationReport),
    DuplexDigest(Vec<crate::chemistry::DigestFragment>),
    Translation(crate::translation::TranslationReport),
    CdsTranslation(crate::translation::CdsTranslation),
    Orfs(crate::translation::OrfReport),
    Alignment(crate::alignment::Alignment),
    PcrDocuments(Vec<crate::verification::PcrDocument>),
    PrimerEvidence(PrimerEvidence),
    PrimerDesign(crate::primers::DesignReport),
    ProtocolPdf(crate::pdf::PdfReport),
    GoldenGate(crate::golden_gate::GoldenGateReport),
    MismatchPcr(crate::mutagenesis::PcrReport),
    InversePcr(crate::mutagenesis::InversePcrDesign),
    Gel(crate::gel::GelReport),
    DiagnosticDigest(crate::gel::DiagnosticReport),
    MismatchPcrPreview(crate::verification::MutagenesisPreview),
    WorkspaceSaved(Vec<u8>),
    WorkspaceLoaded(crate::workspace::Workspace),
    TraceValidation(crate::trace_validation::TraceValidationReport),
    FeatureLibraryScan(crate::feature_library::ScanReport),
    AnnotatedDocument(Molecule),
    AnnotationImport(crate::annotation_formats::AnnotationImportReport),
    AnnotationExport(AnnotationExportReport),
    RepairPreview(crate::verification::RepairPreview),
    AnnealPreview(crate::verification::AnnealPreview),
    PhysicalPreviewDocument(Molecule),
    ViewTracks(crate::view_tracks::TrackReport),
    DocumentSearch(crate::search::SearchReport),
    TracksPrepared(String),
    AlignmentDocument(crate::alignment_formats::AlignmentDocument),
    AlignmentExport(crate::alignment_formats::AlignmentExport),
    AlignmentRows(Vec<Molecule>),
    AlignmentSaved(Vec<u8>),
    ProteinProperties(crate::protein_properties::ProteinReport),
    TailPreview(crate::verification::TailPreview),
    OverlapProfile(crate::overlap_profiles::ProfileReport),
    RunCommit(crate::run_lineage::RunCommit),
    RunReplay(crate::run_lineage::ReplayReport),
    CodonDesign(crate::codon_design::DesignReport),
    CodonCheck(crate::codon_design::CheckReport),
    DotPlot(crate::dotplot::DotPlotReport),
    Contigs(crate::contig::ContigReport),
    Methylation(crate::methylation::MethylationReport),
    CrisprGuides(crate::crispr::GuideReport),
    CrisprOffTargets(crate::crispr::OffTargetReport),
    EditTemplate(crate::edit_template::EditTemplateReport),
    EditTemplateCheck(crate::edit_template::EditTemplateCheck),
    GenomeImportStarted(crate::genome_session::ImportStarted),
    GenomeImportProgress(crate::genome_session::ImportProgress),
    GenomeFastaFinished(crate::genome_session::FastaFinished),
    GenomeRegistered(crate::genome_session::RegisteredGenome),
    GenomeIntervalGuides(crate::genome_session::IntervalGuides),
    GenomeSearchState(crate::genome_session::SearchState),
    GenomeSearchPage(crate::genome_session::SearchPage),
    GenomeDisposed(crate::genome_session::Disposed),
    GenomeIntervalRead(crate::genome_session::IntervalRead),
    RestrictionInspection(crate::restriction_inspection::RestrictionInspection),
    DigestProduct(crate::restriction_inspection::DigestProduct),
    DigestInstalled(crate::restriction_history::InstalledDigest),
    CatalogRestrictionInspection(crate::catalog_digest::CatalogRestrictionInspection),
    CatalogDigestProduct(crate::catalog_digest::CatalogDigestProduct),
    CatalogDigestInstalled(crate::catalog_digest::CatalogDigestInstalled),
    RecombinationInspection(crate::recombination::RecombinationReport),
    RecombinationTargetCheck(crate::recombination::TargetAssessment),
    RecombinationProduct(crate::recombination::CheckedProduct),
    ReadAlignmentsPrepared(crate::read_alignment_session::Prepared),
    ReadAlignmentPage(crate::read_alignment_session::Page),
    ReadAlignmentChecked(crate::read_alignment_session::Checked),
    ReadAlignmentsExported(crate::read_alignment_session::Exported),
    ReadAlignmentsDisposed(crate::read_alignment_session::Disposed),
    TopoSimulation(crate::topo::TopoReport),
    TopoOrientations(crate::topo::OrientationInspection),
    TopoProjection(crate::topo::Product),
    TopoTarget(crate::topo::TargetAssessment),
    AnnotatedTopoProjection(crate::topo::AnnotatedProduct),
    AnnotatedTopoTarget(crate::topo::AnnotatedAssessment),
    VariantsPrepared(crate::variants_session::Prepared),
    VariantPage(crate::variants_session::Page),
    VariantReferenceChecked(crate::variants_session::Checked),
    VariantAllelePreview(crate::variants_session::Previewed),
    VariantsExported(crate::variants_session::Exported),
    VariantsDisposed(crate::variants_session::Disposed),
    OverlapRecipePreview(crate::assembly_recipe::Preview),
    LigationRecipePreview(crate::ligation_recipe::Preview),
    AssemblyJunctionInspection(crate::assembly_junctions::Report),
    OverlapPrimersDesigned(crate::overlap_primer_design::Report),
    AssemblyIntent(crate::assembly_intent::Report),
    DigestFragmentSelected(crate::digest_cut_selection::Report),
    PlanOligosExported(crate::plan_oligos::Exported),
    ProtocolHtml(String),
    CloningOutputsFinalized(crate::cloning_outputs::Commit),
    GoldenGatePrimersDesigned(crate::golden_gate_primer_design::Report),
    MutationIntent(crate::mutation_intent::Report),
    ReadCheckCaptured(crate::read_checks::Commit),
    ReadCheckReplayed(crate::read_checks::Checked),
    ReadReviewPreview(crate::read_review::Preview),
    ReadReviewCommitted(crate::read_review::Commit),
    EditedAssemblyIntent(crate::assembly_junction_edit::Intent),
    EditedAssemblyPrimers(crate::assembly_junction_edit::Report),
    SelectiveAssemblyPrimers(crate::selective_junction_edit::Report),
    GoldenGateJunctionPrimers(crate::golden_gate_junction_edit::Report),
    PreparedGoldenGatePool(crate::prepared_golden_gate::Report),
    PreparedGoldenGateCandidate(crate::prepared_golden_gate::SelectionReport),
    GoldenGateReaction(crate::golden_gate_reaction::Report),
    GoldenGateReactionCandidate(crate::prepared_golden_gate::SelectionReport),
    PcrOutputs(crate::pcr_outputs::Prepared),
    MismatchPcrOutputs(crate::mismatch_pcr_outputs::Prepared),
    ReviewedComparison(crate::pairwise_review::Report),
    ReferenceFragment(crate::reference_transfer::Fragment),
    ReferencePaste(crate::reference_transfer::Report),
    TranslatedReference(crate::translation_reference::Report),
    TranslatedOutputs(crate::translation_reference::Outputs),
    PinnedGel(crate::gel::GelReport),
    ExportedPairwiseComparison(crate::pairwise_export::Export),
    SearchedPairwiseComparison(crate::pairwise_search::Matches),
    PreparedPairwiseContext(crate::pairwise_context::Context),
    CombinatorialBatchStarted(crate::combinatorial_batch::Progress),
    CombinatorialBatchAdvanced(crate::combinatorial_batch::Progress),
    CombinatorialBatchCaptured(crate::combinatorial_batch::Commit),
    CombinatorialReleaseStarted(crate::combinatorial_release::Progress),
    CombinatorialReleaseAdvanced(crate::combinatorial_release::Progress),
    CombinatorialReleaseInspected(crate::combinatorial_release::Inspection),
    CombinatorialReleaseFinished(crate::combinatorial_release::Commit),
    CombinatorialOutputsArchive(crate::combinatorial_export::Archive),
    GoldenGateFusions(crate::golden_gate_fusion_design::Report),
    GoldenGateFusionMaterialized(crate::golden_gate_fusion_design::Materialized),
    LigationFidelity(crate::ligation_fidelity::Report),
    DomesticationInspected(crate::domestication::Inspection),
    DomesticationProposed(crate::domestication::Proposal),
    DomesticationChecked {
        report: crate::domestication::CheckReport,
        target: Molecule,
    },
    DomesticationPrepared(crate::domestication_preparation::Prepared),
    LigationProfiles(Vec<crate::ligation_profiles::Summary>),
    PublishedLigationScore(crate::ligation_profiles::Report),
    RankedGoldenGateFusions(crate::golden_gate_fusion_rank::Report),
    LigationCompetition(crate::ligation_fidelity::CompetitionReport),
    TopoNicksSealed(crate::topo_nick_sealing::Report),
    ProcessedTopoProduct(crate::topo_processed_product::Report),
    HistoricalProcessedRunEvidence(crate::processed_run_evidence::Download),
    TopoActionsProduct(crate::topo_actions::Report),
    TopoActionsInspected(crate::topo_actions::Inspection),
    MultiSiteRecombination(crate::multisite_recombination::Batch),
    MultiSiteRecombinationProduct(crate::multisite_recombination::ProductDocument),
    MultiSiteRecombinationChecked(crate::multisite_recombination::CheckedProduct),
    CombinatorialMultiSiteInspection(crate::combinatorial_design::MultiSiteInspection),
    TailingRecipePreview(crate::tailing_recipe::Preview),
    PcrDigestInspection(crate::pcr_digest_inspection::Preview),
    RunMaterialInspection(crate::run_material::Inspection),
    RunPcrDigestInspection(crate::pcr_digest_inspection::RunPreview),
    RunCatalogRestrictionInspection(crate::run_catalog_digest::Inspection),
    SubstitutedAssemblyIntent(crate::assembly_substitution_intent::Report),
    MismatchPcrDigestInspection(crate::mismatch_pcr_digest_inspection::Preview),
    ContributionEditIntent(crate::assembly_edit_intent::Report),
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationExportReport {
    pub text: String,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerEvidence {
    pub analysis: crate::primers::PrimerAnalysis,
    pub bindings: crate::primers::BindingReport,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Response {
    pub id: u32,
    pub revision: String,
    pub result: AnalysisResult,
}

/// Versions 1/2 did not record the scientific engine in verification reports.
/// Migration must never attribute those historical findings to today's engine.
pub const LEGACY_ENGINE_VERSION: &str = "unrecorded/legacy-worker-v1-v2";

#[derive(Deserialize)]
struct LegacyVerificationReport {
    plan_revision: String,
    target_revision: String,
    status: crate::verification::Status,
    findings: Vec<String>,
    assumptions: Vec<String>,
    products: Vec<Molecule>,
    exhaustive_alternative_search: bool,
}
impl From<LegacyVerificationReport> for crate::verification::VerificationReport {
    fn from(old: LegacyVerificationReport) -> Self {
        Self {
            plan_revision: old.plan_revision,
            target_revision: old.target_revision,
            engine_version: LEGACY_ENGINE_VERSION.into(),
            status: old.status,
            findings: old.findings,
            assumptions: old.assumptions,
            products: old.products,
            exhaustive_alternative_search: old.exhaustive_alternative_search,
            step_reports: Vec::new(),
            annotation_reports: Vec::new(),
            processed_step_reports: Vec::new(),
            physical_processing_reports: Vec::new(),
            reference_step_reports: Vec::new(),
        }
    }
}
// Versions 3–24 include engine provenance but no per-step physical graph.
// Keep this positional layout separate: serde defaults do not extend postcard
// structs at EOF. Historical reports gain no invented physical evidence.
#[derive(Deserialize)]
struct LegacyVerificationReportV24 {
    plan_revision: String,
    target_revision: String,
    engine_version: String,
    status: crate::verification::Status,
    findings: Vec<String>,
    assumptions: Vec<String>,
    products: Vec<Molecule>,
    exhaustive_alternative_search: bool,
}
#[derive(Deserialize)]
#[allow(dead_code)] // Positional selectors 0–3 precede the historical verification layout.
enum LegacyVerificationResultV24 {
    Search(Vec<Match>),
    Restriction(Vec<CutSite>),
    Error(String),
    Pcr(Vec<crate::analysis::PcrProduct>),
    Verification(LegacyVerificationReportV24),
}
#[derive(Deserialize)]
struct LegacyVerificationResponseV24 {
    id: u32,
    revision: String,
    result: LegacyVerificationResultV24,
}
impl From<LegacyVerificationReportV24> for crate::verification::VerificationReport {
    fn from(old: LegacyVerificationReportV24) -> Self {
        Self {
            plan_revision: old.plan_revision,
            target_revision: old.target_revision,
            engine_version: old.engine_version,
            status: old.status,
            findings: old.findings,
            assumptions: old.assumptions,
            products: old.products,
            exhaustive_alternative_search: old.exhaustive_alternative_search,
            step_reports: Vec::new(),
            annotation_reports: Vec::new(),
            processed_step_reports: Vec::new(),
            physical_processing_reports: Vec::new(),
            reference_step_reports: Vec::new(),
        }
    }
}
// Schemas25/26 retain physical step reports but predate annotation evidence.
// Their positional payload is migrated explicitly, without inventing mappings.
#[derive(Deserialize)]
struct LegacyVerificationReportV26 {
    plan_revision: String,
    target_revision: String,
    engine_version: String,
    status: crate::verification::Status,
    findings: Vec<String>,
    assumptions: Vec<String>,
    products: Vec<Molecule>,
    exhaustive_alternative_search: bool,
    step_reports: Vec<crate::verification::StepReport>,
}
#[derive(Deserialize)]
#[allow(dead_code)]
enum LegacyVerificationResultV26 {
    Search(Vec<Match>),
    Restriction(Vec<CutSite>),
    Error(String),
    Pcr(Vec<crate::analysis::PcrProduct>),
    Verification(LegacyVerificationReportV26),
}
#[derive(Deserialize)]
struct LegacyVerificationResponseV26 {
    id: u32,
    revision: String,
    result: LegacyVerificationResultV26,
}
impl From<LegacyVerificationReportV26> for crate::verification::VerificationReport {
    fn from(old: LegacyVerificationReportV26) -> Self {
        Self {
            plan_revision: old.plan_revision,
            target_revision: old.target_revision,
            engine_version: old.engine_version,
            status: old.status,
            findings: old.findings,
            assumptions: old.assumptions,
            products: old.products,
            exhaustive_alternative_search: old.exhaustive_alternative_search,
            step_reports: old.step_reports,
            annotation_reports: Vec::new(),
            processed_step_reports: Vec::new(),
            physical_processing_reports: Vec::new(),
            reference_step_reports: Vec::new(),
        }
    }
}
// Separate historical enums encode the actual old layouts, not the present DTO
// under an older header. In particular, old Verification lacks engine_version.
#[derive(Deserialize)]
enum LegacyResultV1 {
    Search(Vec<Match>),
    Restriction(Vec<CutSite>),
    Error(String),
    Pcr(Vec<crate::analysis::PcrProduct>),
    Verification(LegacyVerificationReport),
}
#[derive(Deserialize)]
struct LegacyResponseV1 {
    id: u32,
    revision: String,
    result: LegacyResultV1,
}
#[derive(Deserialize)]
enum LegacyResultV2 {
    Search(Vec<Match>),
    Restriction(Vec<CutSite>),
    Error(String),
    Pcr(Vec<crate::analysis::PcrProduct>),
    Verification(LegacyVerificationReport),
    DuplexDigest(Vec<crate::chemistry::DigestFragment>),
}
#[derive(Deserialize)]
struct LegacyResponseV2 {
    id: u32,
    revision: String,
    result: LegacyResultV2,
}
impl From<LegacyResultV1> for AnalysisResult {
    fn from(old: LegacyResultV1) -> Self {
        match old {
            LegacyResultV1::Search(value) => Self::Search(value),
            LegacyResultV1::Restriction(value) => Self::Restriction(value),
            LegacyResultV1::Error(value) => Self::Error(value),
            LegacyResultV1::Pcr(value) => Self::Pcr(value),
            LegacyResultV1::Verification(value) => Self::Verification(value.into()),
        }
    }
}
impl From<LegacyResultV2> for AnalysisResult {
    fn from(old: LegacyResultV2) -> Self {
        match old {
            LegacyResultV2::Search(value) => Self::Search(value),
            LegacyResultV2::Restriction(value) => Self::Restriction(value),
            LegacyResultV2::Error(value) => Self::Error(value),
            LegacyResultV2::Pcr(value) => Self::Pcr(value),
            LegacyResultV2::Verification(value) => Self::Verification(value.into()),
            LegacyResultV2::DuplexDigest(value) => Self::DuplexDigest(value),
        }
    }
}
pub fn encode_request(r: &Request) -> Result<Vec<u8>, Error> {
    if let AnalysisJob::ExportPairwiseComparison {
        record,
        reference_name,
        query_name,
        ..
    } = &r.job
    {
        crate::pairwise_export::preflight(record, reference_name, query_name)?;
    }
    if let AnalysisJob::SearchPairwiseComparison {
        report,
        pattern,
        options,
    } = &r.job
    {
        crate::pairwise_search::preflight(report, pattern, options)?;
    }
    if let AnalysisJob::PreparePairwiseContext {
        report,
        reference,
        query,
        translate,
    } = &r.job
    {
        crate::pairwise_context::preflight(report, reference, query, *translate)?;
    }
    match &r.job {
        AnalysisJob::InspectMismatchPcrDigest {
            molecule,
            origin,
            request,
        } => {
            crate::mismatch_pcr_digest_inspection::source_preflight(
                molecule,
                origin.as_ref(),
                request,
            )?;
        }
        AnalysisJob::BuildContributionEditIntent { inputs, request } => {
            if inputs.len() > 32 {
                return Err(Error::TooLarge);
            }
            let sources = inputs
                .iter()
                .map(|i| (i.id.as_str(), &i.molecule))
                .collect::<Vec<_>>();
            crate::assembly_edit_intent::source_preflight(&sources, request)?;
        }
        AnalysisJob::BuildSubstitutedAssemblyIntent { inputs, request } => {
            if inputs.len() > 32 {
                return Err(Error::TooLarge);
            }
            let sources = inputs
                .iter()
                .map(|i| (i.id.as_str(), &i.molecule))
                .collect::<Vec<_>>();
            crate::assembly_substitution_intent::source_preflight(&sources, request)?;
        }
        AnalysisJob::InspectRunCatalogRestriction {
            source_role,
            molecule,
            origin,
            selection,
            options,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::run_catalog_digest::source_preflight(
                &catalog,
                selection,
                source_role,
                molecule,
                origin,
                options,
            )?;
        }
        AnalysisJob::InspectRunPcrDigest {
            molecule,
            origin,
            request,
        } => {
            crate::pcr_digest_inspection::preflight_from_run(molecule, origin, request)?;
        }
        AnalysisJob::InspectRunMaterial {
            workspace,
            event,
            selection,
        } => {
            crate::run_material::inspection_revision(workspace, event, selection)?;
        }
        AnalysisJob::InspectPcrDigest { molecule, request } => {
            crate::pcr_digest_inspection::source_preflight(molecule, request)?;
        }
        AnalysisJob::PreviewTailingRecipe {
            inputs,
            target,
            recipe,
        } => {
            crate::tailing_recipe::preflight(inputs, target, recipe)?;
        }
        AnalysisJob::InspectCombinatorialMultiSite {
            workspace,
            design,
            construct,
        } => {
            crate::combinatorial_design::inspect_multisite_revision(design, construct, workspace)?;
        }
        AnalysisJob::SimulateMultiSiteRecombination { inputs, request } => {
            crate::multisite_recombination::preflight(inputs, request)?;
        }
        AnalysisJob::ProjectMultiSiteRecombination {
            inputs,
            request,
            selection,
        } => {
            crate::multisite_recombination::product_revision(inputs, request, selection)?;
        }
        AnalysisJob::CheckMultiSiteRecombination {
            inputs,
            request,
            selection,
            target,
        } => {
            crate::multisite_recombination::check_revision(inputs, request, selection, target)?;
        }
        AnalysisJob::ProcessTopoActions { inputs, request } => {
            crate::topo_actions::preflight(inputs, request)?
        }
        AnalysisJob::InspectTopoActions { inputs, request } => {
            crate::topo_actions::preflight_actions(inputs, request)?
        }
        AnalysisJob::ProjectProcessedTopo { inputs, request } => {
            crate::topo_processed_product::preflight(inputs, request)?
        }
        AnalysisJob::SealTopoNicks { inputs, request } => {
            crate::topo_nick_sealing::preflight(inputs, request)?;
        }
        AnalysisJob::ListLigationProfiles => {
            crate::ligation_profiles::catalog_revision()?;
        }
        AnalysisJob::ScorePublishedLigation {
            profile,
            ends,
            physical_competitors,
        } => {
            crate::ligation_profiles::score_revision_with_competitors(
                profile,
                ends,
                physical_competitors,
            )?;
        }
        AnalysisJob::ScoreLigationCompetition {
            dataset,
            ends,
            physical_competitors,
        } => {
            crate::ligation_fidelity::competition_revision(dataset, ends, physical_competitors)?;
        }
        AnalysisJob::RankGoldenGateFusions { inputs, request } => {
            crate::golden_gate_fusion_rank::preflight(
                &inputs
                    .iter()
                    .map(|i| (i.id.as_str(), &i.molecule))
                    .collect::<Vec<_>>(),
                request,
            )?;
        }
        AnalysisJob::InspectDomestication { source, request }
        | AnalysisJob::ProposeDomestication { source, request }
        | AnalysisJob::CheckDomestication {
            source, request, ..
        } => {
            crate::domestication::preflight(source, request)?;
        }
        AnalysisJob::PrepareDomestication { input, request } => {
            crate::domestication_preparation::preflight(input, request)?;
        }
        AnalysisJob::SearchGoldenGateFusions { inputs, request }
        | AnalysisJob::MaterializeGoldenGateFusions {
            inputs, request, ..
        } => {
            crate::golden_gate_fusion_design::preflight(
                &inputs
                    .iter()
                    .map(|input| (input.id.as_str(), &input.molecule))
                    .collect::<Vec<_>>(),
                request,
            )?;
        }
        AnalysisJob::ScoreLigationFidelity { dataset, ends } => {
            crate::ligation_fidelity::revision(dataset, ends)?;
        }
        AnalysisJob::BeginCombinatorialBatch { workspace, matrix } => {
            crate::combinatorial_batch::begin_revision(workspace, matrix)?;
        }
        AnalysisJob::CaptureCombinatorialBatch {
            workspace,
            current_context,
        } => {
            crate::combinatorial_batch::capture_revision(workspace, current_context)?;
        }
        AnalysisJob::BeginCombinatorialRelease { workspace, request } => {
            crate::combinatorial_release::preflight(workspace, request)?;
        }
        AnalysisJob::ExportCombinatorialOutputs {
            workspace,
            manifest_hash,
            options,
        } => {
            crate::combinatorial_export::preflight(workspace, manifest_hash, options)?;
        }
        _ => {}
    }
    let maximum = if matches!(
        r.job,
        AnalysisJob::SaveWorkspace { .. }
            | AnalysisJob::LoadWorkspace { .. }
            | AnalysisJob::RecordCloningRun { .. }
            | AnalysisJob::ReplayCloningRun { .. }
            | AnalysisJob::LoadProcessedRunEvidence { .. }
            | AnalysisJob::FinalizeCloningOutputs { .. }
            | AnalysisJob::CaptureReadCheck { .. }
            | AnalysisJob::ReplayReadCheck { .. }
            | AnalysisJob::PreviewReadReview { .. }
            | AnalysisJob::CommitReadReview { .. }
            | AnalysisJob::InstallDigestFragment { .. }
            | AnalysisJob::InstallCatalogDigestFragment { .. }
            | AnalysisJob::BeginCombinatorialBatch { .. }
            | AnalysisJob::CaptureCombinatorialBatch { .. }
            | AnalysisJob::BeginCombinatorialRelease { .. }
            | AnalysisJob::FinishCombinatorialRelease { .. }
            | AnalysisJob::ExportCombinatorialOutputs { .. }
            | AnalysisJob::InspectCombinatorialMultiSite { .. }
            | AnalysisJob::InspectRunMaterial { .. }
    ) {
        MAX_WORKSPACE_MESSAGE_BYTES
    } else {
        MAX_IMPORT_BYTES
    };
    encode_with_limit(b"CGQ\0", r, maximum)
}
pub fn decode_response(bytes: &[u8]) -> Result<Response, Error> {
    // Read the selector only after checking its bounds; decode still validates
    // magic/version/length/checksum and the complete selected payload layout.
    let maximum = worker_limit(
        b"CGR\0",
        bytes,
        &[
            (20, 6),
            (21, 6),
            (40, 12),
            (41, 12),
            (62, 20),
            (65, 21),
            (94, 39),
            (97, 42),
            (100, 43),
            (122, 57),
            (126, 58),
            (127, 58),
        ],
    )?;
    preflight_codon_tag(bytes, &[42, 43])?;
    preflight_dotplot_tag(bytes, &[44])?;
    preflight_contig_tag(bytes, &[45])?;
    preflight_methylation_tag(bytes, &[46])?;
    preflight_crispr_tag(bytes, &[47, 48, 49, 50])?;
    preflight_genome_tag(bytes, &[51, 52, 53, 54, 55, 56, 57, 58, 59])?;
    preflight_restriction_tag(bytes, &[60, 61, 62])?;
    preflight_catalog_tag(bytes, &[63, 64, 65])?;
    preflight_recombination_tag(bytes, &[66, 67, 68])?;
    preflight_read_alignment_tag(bytes, &[69, 70, 71, 72, 73])?;
    preflight_topo_tag(bytes, &[74, 75, 76, 77])?;
    preflight_annotated_topo_tag(bytes, &[78, 79])?;
    preflight_variants_tag(bytes, &[80, 81, 82, 83, 84, 85])?;
    preflight_overlap_recipe_tag(bytes, &[86])?;
    preflight_ligation_recipe_tag(bytes, &[87])?;
    preflight_assembly_junctions_tag(bytes, &[88])?;
    preflight_overlap_primer_design_tag(bytes, &[89])?;
    preflight_assembly_intent_tag(bytes, &[90])?;
    preflight_cloning_output_tag(bytes, &[91, 92, 93])?;
    preflight_finalize_tag(bytes, &[94])?;
    preflight_golden_gate_primer_tag(bytes, &[95])?;
    preflight_mutation_intent_tag(bytes, &[96])?;
    preflight_read_check_tag(bytes, &[97, 98])?;
    preflight_read_review_tag(bytes, &[99, 100])?;
    preflight_edited_assembly_tag(bytes, &[101, 102])?;
    preflight_selective_assembly_tag(bytes, &[103])?;
    preflight_golden_gate_junction_tag(bytes, &[104])?;
    preflight_prepared_golden_gate_tag(bytes, &[105, 106])?;
    preflight_golden_gate_reaction_tag(bytes, &[107, 108])?;
    preflight_pcr_outputs_tag(bytes, &[109])?;
    preflight_mismatch_pcr_outputs_tag(bytes, &[110])?;
    preflight_pairwise_review_tag(bytes, &[111])?;
    preflight_reference_transfer_tag(bytes, &[112, 113])?;
    preflight_translation_reference_tag(bytes, &[114, 115])?;
    preflight_pinned_gel_tag(bytes, &[116])?;
    preflight_pairwise_export_tag(bytes, &[117, 118, 119])?;
    preflight_combinatorial_batch_tag(bytes, &[120, 121, 122])?;
    preflight_combinatorial_release_tag(bytes, &[123, 124, 125, 126, 127])?;
    preflight_fusion_design_tag(bytes, &[128, 129, 130])?;
    preflight_domestication_tag(bytes, &[131, 132, 133, 134])?;
    preflight_published_ligation_tag(bytes, &[135, 136, 137, 138])?;
    preflight_nick_sealing_tag(bytes, &[139])?;
    preflight_processed_topo_tag(bytes, &[140])?;
    preflight_processed_run_tag(bytes, &[141])?;
    preflight_topo_actions_tag(bytes, &[142, 143])?;
    preflight_multisite_recombination_tag(bytes, &[144, 145, 146])?;
    preflight_combinatorial_multisite_tag(bytes, &[147])?;
    preflight_tailing_recipe_tag(bytes, &[148])?;
    preflight_pcr_digest_tag(bytes, &[149])?;
    preflight_run_material_tag(bytes, &[150, 151])?;
    preflight_run_catalog_tag(bytes, &[152])?;
    preflight_substitution_intent_tag(bytes, &[153])?;
    preflight_mismatch_digest_tag(bytes, &[154])?;
    preflight_contribution_edit_tag(bytes, &[155])?;
    let version = bytes.get(4..6).map(|v| u16::from_le_bytes([v[0], v[1]]));
    let response: Response = match version {
        Some(1) => {
            let old: LegacyResponseV1 = decode(b"CGR\0", bytes)?;
            Response {
                id: old.id,
                revision: old.revision,
                result: old.result.into(),
            }
        }
        Some(2) => {
            let old: LegacyResponseV2 = decode(b"CGR\0", bytes)?;
            Response {
                id: old.id,
                revision: old.revision,
                result: old.result.into(),
            }
        }
        Some(3..=24)
            if postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
                .is_ok_and(|((_, _, tag), _)| tag == 4) =>
        {
            let old: LegacyVerificationResponseV24 = decode_with_limit(b"CGR\0", bytes, maximum)?;
            let LegacyVerificationResultV24::Verification(report) = old.result else {
                return Err(Error::Format(
                    "historical verification selector mismatch".into(),
                ));
            };
            Response {
                id: old.id,
                revision: old.revision,
                result: AnalysisResult::Verification(report.into()),
            }
        }
        Some(25 | 26)
            if postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
                .is_ok_and(|((_, _, tag), _)| tag == 4) =>
        {
            let old: LegacyVerificationResponseV26 = decode_with_limit(b"CGR\0", bytes, maximum)?;
            let LegacyVerificationResultV26::Verification(report) = old.result else {
                return Err(Error::Format(
                    "historical physical verification selector mismatch".into(),
                ));
            };
            Response {
                id: old.id,
                revision: old.revision,
                result: AnalysisResult::Verification(report.into()),
            }
        }
        _ => match crate::protocol_legacy_v64::decode(bytes, maximum)? {
            Some(old) => old,
            None => decode_with_limit(b"CGR\0", bytes, maximum)?,
        },
    };
    let required = match &response.result {
        AnalysisResult::ContributionEditIntent(_) => 79,
        AnalysisResult::MismatchPcrDigestInspection(_) => 78,
        AnalysisResult::SubstitutedAssemblyIntent(_) => 77,
        AnalysisResult::RunCatalogRestrictionInspection(_) => 75,
        AnalysisResult::RunMaterialInspection(_) | AnalysisResult::RunPcrDigestInspection(_) => 74,
        AnalysisResult::PcrDigestInspection(_) => 72,
        AnalysisResult::TailingRecipePreview(report) => {
            check_variant_version(bytes, plan_worker_version(&report.compiled.plan))?;
            70
        }
        AnalysisResult::CombinatorialMultiSiteInspection(_) => 69,
        AnalysisResult::MultiSiteRecombination(_)
        | AnalysisResult::MultiSiteRecombinationProduct(_)
        | AnalysisResult::MultiSiteRecombinationChecked(_) => 68,
        AnalysisResult::TopoActionsProduct(_) | AnalysisResult::TopoActionsInspected(_) => 66,
        AnalysisResult::HistoricalProcessedRunEvidence(_) => 65,
        AnalysisResult::ProcessedTopoProduct(_) => 64,
        AnalysisResult::TopoNicksSealed(_) => 63,
        AnalysisResult::LigationProfiles(_)
        | AnalysisResult::PublishedLigationScore(_)
        | AnalysisResult::RankedGoldenGateFusions(_)
        | AnalysisResult::LigationCompetition(_) => 62,
        AnalysisResult::DomesticationInspected(_)
        | AnalysisResult::DomesticationProposed(_)
        | AnalysisResult::DomesticationChecked { .. }
        | AnalysisResult::DomesticationPrepared(_) => 61,
        AnalysisResult::GoldenGateFusions(_)
        | AnalysisResult::GoldenGateFusionMaterialized(_)
        | AnalysisResult::LigationFidelity(_) => 60,
        AnalysisResult::PcrOutputs(_) => 50,
        AnalysisResult::MismatchPcrOutputs(_) => 51,
        AnalysisResult::ReviewedComparison(_) => 52,
        AnalysisResult::ReferenceFragment(_) | AnalysisResult::ReferencePaste(_) => 53,
        AnalysisResult::TranslatedReference(_) | AnalysisResult::TranslatedOutputs(_) => 54,
        AnalysisResult::PinnedGel(_) => 55,
        AnalysisResult::ExportedPairwiseComparison(_) => 56,
        AnalysisResult::SearchedPairwiseComparison(_) => 56,
        AnalysisResult::PreparedPairwiseContext(_) => 56,
        AnalysisResult::CombinatorialBatchStarted(_)
        | AnalysisResult::CombinatorialBatchAdvanced(_)
        | AnalysisResult::CombinatorialBatchCaptured(_) => 57,
        AnalysisResult::CombinatorialReleaseStarted(_)
        | AnalysisResult::CombinatorialReleaseAdvanced(_)
        | AnalysisResult::CombinatorialReleaseInspected(_)
        | AnalysisResult::CombinatorialReleaseFinished(_)
        | AnalysisResult::CombinatorialOutputsArchive(_) => 58,
        AnalysisResult::GoldenGateReaction(_) | AnalysisResult::GoldenGateReactionCandidate(_) => {
            49
        }
        AnalysisResult::EditedAssemblyIntent(_) | AnalysisResult::EditedAssemblyPrimers(_) => 44,
        AnalysisResult::SelectiveAssemblyPrimers(_) => 45,
        AnalysisResult::GoldenGateJunctionPrimers(_) => 46,
        AnalysisResult::PreparedGoldenGatePool(_)
        | AnalysisResult::PreparedGoldenGateCandidate(_) => 47,
        AnalysisResult::ReadReviewPreview(_) | AnalysisResult::ReadReviewCommitted(_) => 43,
        AnalysisResult::ReadCheckCaptured(_) | AnalysisResult::ReadCheckReplayed(_) => 42,
        AnalysisResult::MutationIntent(_) => 41,
        AnalysisResult::GoldenGatePrimersDesigned(_) => 40,
        AnalysisResult::CloningOutputsFinalized(_) => 39,
        AnalysisResult::DigestFragmentSelected(_)
        | AnalysisResult::PlanOligosExported(_)
        | AnalysisResult::ProtocolHtml(_) => 38,
        AnalysisResult::AssemblyIntent(_) => 36,
        AnalysisResult::OverlapPrimersDesigned(_) => 35,
        AnalysisResult::AssemblyJunctionInspection(_) => 33,
        AnalysisResult::LigationRecipePreview(report) => {
            check_variant_version(bytes, plan_worker_version(&report.compiled.plan))?;
            plan_worker_version(&report.compiled.plan).max(32)
        }
        AnalysisResult::OverlapRecipePreview(report) => {
            check_variant_version(bytes, plan_worker_version(&report.compiled.plan))?;
            plan_worker_version(&report.compiled.plan).max(31)
        }
        AnalysisResult::VariantsPrepared(_)
        | AnalysisResult::VariantPage(_)
        | AnalysisResult::VariantReferenceChecked(_)
        | AnalysisResult::VariantAllelePreview(_)
        | AnalysisResult::VariantsExported(_)
        | AnalysisResult::VariantsDisposed(_) => 29,
        AnalysisResult::AnnotatedTopoProjection(_) | AnalysisResult::AnnotatedTopoTarget(_) => 27,
        AnalysisResult::TopoSimulation(_)
        | AnalysisResult::TopoOrientations(_)
        | AnalysisResult::TopoProjection(_)
        | AnalysisResult::TopoTarget(_) => 26,
        AnalysisResult::ReadAlignmentsPrepared(_)
        | AnalysisResult::ReadAlignmentPage(_)
        | AnalysisResult::ReadAlignmentChecked(_)
        | AnalysisResult::ReadAlignmentsExported(_)
        | AnalysisResult::ReadAlignmentsDisposed(_) => 24,
        AnalysisResult::RecombinationInspection(_)
        | AnalysisResult::RecombinationTargetCheck(_)
        | AnalysisResult::RecombinationProduct(_) => 22,
        AnalysisResult::CatalogRestrictionInspection(_)
        | AnalysisResult::CatalogDigestProduct(_)
        | AnalysisResult::CatalogDigestInstalled(_) => 21,
        AnalysisResult::RestrictionInspection(_)
        | AnalysisResult::DigestProduct(_)
        | AnalysisResult::DigestInstalled(_) => 20,
        AnalysisResult::GenomeImportStarted(_)
        | AnalysisResult::GenomeImportProgress(_)
        | AnalysisResult::GenomeFastaFinished(_)
        | AnalysisResult::GenomeRegistered(_)
        | AnalysisResult::GenomeIntervalGuides(_)
        | AnalysisResult::GenomeSearchState(_)
        | AnalysisResult::GenomeSearchPage(_)
        | AnalysisResult::GenomeDisposed(_)
        | AnalysisResult::GenomeIntervalRead(_) => 19,
        AnalysisResult::Search(_)
        | AnalysisResult::Restriction(_)
        | AnalysisResult::Error(_)
        | AnalysisResult::Pcr(_)
        | AnalysisResult::Verification(_) => 1,
        AnalysisResult::DuplexDigest(_) => 2,
        AnalysisResult::Translation(_)
        | AnalysisResult::CdsTranslation(_)
        | AnalysisResult::Orfs(_)
        | AnalysisResult::Alignment(_) => 3,
        AnalysisResult::PcrDocuments(_)
        | AnalysisResult::PrimerEvidence(_)
        | AnalysisResult::PrimerDesign(_)
        | AnalysisResult::ProtocolPdf(_) => 4,
        AnalysisResult::WorkspaceSaved(_)
        | AnalysisResult::WorkspaceLoaded(_)
        | AnalysisResult::TraceValidation(_) => 6,
        AnalysisResult::CrisprGuides(_)
        | AnalysisResult::CrisprOffTargets(_)
        | AnalysisResult::EditTemplate(_)
        | AnalysisResult::EditTemplateCheck(_) => 18,
        AnalysisResult::Methylation(_) => 16,
        AnalysisResult::Contigs(_) => 15,
        AnalysisResult::DotPlot(_) => 14,
        AnalysisResult::CodonDesign(_) | AnalysisResult::CodonCheck(_) => 13,
        AnalysisResult::RunCommit(_) | AnalysisResult::RunReplay(_) => 12,
        AnalysisResult::TailPreview(_) | AnalysisResult::OverlapProfile(_) => 11,
        AnalysisResult::ProteinProperties(_) => 10,
        AnalysisResult::AlignmentDocument(_)
        | AnalysisResult::AlignmentExport(_)
        | AnalysisResult::AlignmentRows(_)
        | AnalysisResult::AlignmentSaved(_) => 10,
        AnalysisResult::ViewTracks(_)
        | AnalysisResult::DocumentSearch(_)
        | AnalysisResult::TracksPrepared(_) => 9,
        AnalysisResult::RepairPreview(_)
        | AnalysisResult::AnnealPreview(_)
        | AnalysisResult::PhysicalPreviewDocument(_) => 8,
        AnalysisResult::FeatureLibraryScan(_)
        | AnalysisResult::AnnotatedDocument(_)
        | AnalysisResult::AnnotationImport(_)
        | AnalysisResult::AnnotationExport(_) => 7,
        _ => 5,
    };
    check_variant_version(bytes, required)?;
    check_preparation_version(bytes, &response.result)?;
    // Existing generic AnnealOligos plans keep their historical version. These
    // recipe previews newly admit annealing as a row preparation in schema71.
    match &response.result {
        AnalysisResult::OverlapRecipePreview(p)
            if p.compiled
                .rows
                .iter()
                .any(|r| annealed_preparation(&p.compiled.plan, r)) =>
        {
            check_variant_version(bytes, 71)?
        }
        AnalysisResult::LigationRecipePreview(p)
            if p.compiled
                .rows
                .iter()
                .any(|r| annealed_preparation(&p.compiled.plan, &r.preparation)) =>
        {
            check_variant_version(bytes, 71)?
        }
        AnalysisResult::TailingRecipePreview(p)
            if p.compiled
                .rows
                .iter()
                .any(|r| annealed_preparation(&p.compiled.plan, &r.preparation)) =>
        {
            check_variant_version(bytes, 71)?
        }
        _ => {}
    }
    let recipe_plan = match &response.result {
        AnalysisResult::OverlapRecipePreview(p) => Some(&p.compiled.plan),
        AnalysisResult::LigationRecipePreview(p) => Some(&p.compiled.plan),
        AnalysisResult::TailingRecipePreview(p) => Some(&p.compiled.plan),
        _ => None,
    };
    // The existing Operation stays worker81-compatible, but interpreting it as
    // a forked recipe-row preparation requires the ownership model of schema82.
    if recipe_plan.is_some_and(|p| {
        p.steps
            .iter()
            .any(|op| matches!(op, crate::verification::Operation::OverlapExtension { .. }))
    }) {
        check_variant_version(bytes, 82)?;
    }
    if recipe_plan.is_some_and(|p| {
        p.steps.iter().any(|op| {
            matches!(
                op,
                crate::verification::Operation::ReconstructRunProduct { .. }
            )
        })
    }) {
        check_variant_version(bytes, 74)?;
    }
    if recipe_plan.is_some_and(|p| {
        p.steps
            .iter()
            .any(|op| matches!(op, crate::verification::Operation::PcrWithMismatches { .. }))
    }) {
        check_variant_version(bytes, 76)?;
    }
    if recipe_plan.is_some_and(|p| {
        p.steps.iter().any(|op| match op {
            crate::verification::Operation::CatalogDigest {
                template: crate::verification::Source::Step(index),
                ..
            } => matches!(
                p.steps.get(*index),
                Some(crate::verification::Operation::PcrWithMismatches { .. })
            ),
            _ => false,
        })
    }) {
        check_variant_version(bytes, 78)?;
    }
    if let AnalysisResult::AnnotationImport(report) = &response.result {
        check_annotation_format_version(bytes, report.format)?;
    }
    if let AnalysisResult::WorkspaceLoaded(workspace) = &response.result {
        check_workspace_version(bytes, workspace)?;
        if workspace
            .attachments
            .values()
            .any(|a| a.kind == crate::run_lineage::WORKSPACE_ATTACHMENT_KIND)
        {
            check_variant_version(bytes, 12)?;
        }
        if workspace
            .attachments
            .values()
            .any(|a| a.kind == "read-assembly-v1")
        {
            check_variant_version(bytes, 15)?;
        }
        for c in &workspace.candidates {
            check_variant_version(bytes, plan_worker_version(&c.plan))?;
        }
    }
    match &response.result {
        AnalysisResult::TopoActionsProduct(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::TopoActionsInspected(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::HistoricalProcessedRunEvidence(value) => {
            if value.run_id.len() != 64
                || value.recorded_engine.len() > 128
                || value.step_count > 64
                || value
                    .bytes
                    .as_ref()
                    .is_some_and(|b| b.len() > crate::processed_run_evidence::MAX_BYTES)
                || value.bytes.is_none() != (value.step_count == 0)
            {
                return Err(Error::Format(
                    "Invalid historical processed graph download".into(),
                ));
            }
        }
        AnalysisResult::ProcessedTopoProduct(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::TopoNicksSealed(report) => report.validate_binding(&response.revision)?,
        AnalysisResult::LigationProfiles(items) => {
            let expected = crate::ligation_profiles::list()?;
            if items.len() != expected.len()
                || response.revision != crate::ligation_profiles::catalog_revision()?
                || crate::assembly_recipe::hash(items)? != crate::assembly_recipe::hash(&expected)?
            {
                return Err(Error::Format("Published ligation catalog changed".into()));
            }
        }
        AnalysisResult::PublishedLigationScore(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::RankedGoldenGateFusions(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::LigationCompetition(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::DomesticationInspected(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::DomesticationProposed(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::DomesticationChecked { report, target } => {
            report.validate_binding(&response.revision)?;
            target.validate()?;
            if report.target_revision != target.revision() {
                return Err(Error::Format("Domestication reference changed".into()));
            }
        }
        AnalysisResult::DomesticationPrepared(report) => {
            check_verification_report_version(bytes, &report.verification)?;
            check_variant_version(bytes, plan_worker_version(&report.preparation.plan))?;
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::GoldenGateFusions(report) => report.validate_binding(&response.revision)?,
        AnalysisResult::GoldenGateFusionMaterialized(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::LigationFidelity(report) => report.validate_binding(&response.revision)?,
        AnalysisResult::EditedAssemblyIntent(intent) => {
            intent.validate_binding(&response.revision)?
        }
        AnalysisResult::EditedAssemblyPrimers(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::SelectiveAssemblyPrimers(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::GoldenGateJunctionPrimers(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::PreparedGoldenGatePool(report) => {
            check_verification_report_version(bytes, &report.preview.verification)?;
            check_variant_version(bytes, plan_worker_version(&report.preview.compiled.plan))?;
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::PreparedGoldenGateCandidate(report) => {
            check_variant_version(bytes, plan_worker_version(&report.plan))?;
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::PcrOutputs(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::MismatchPcrOutputs(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::ReviewedComparison(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::ExportedPairwiseComparison(export) => export.validate()?,
        AnalysisResult::SearchedPairwiseComparison(matches) => matches.validate()?,
        AnalysisResult::PreparedPairwiseContext(context) => {
            context.validate()?;
            if context.context_revision != response.revision {
                return Err(Error::Format("Pairwise context response changed.".into()));
            }
        }
        AnalysisResult::CombinatorialBatchStarted(progress)
        | AnalysisResult::CombinatorialBatchAdvanced(progress) => {
            progress.validate()?;
            if progress.context_revision != response.revision {
                return Err(Error::Format(
                    "Combinatorial batch response context changed.".into(),
                ));
            }
        }
        AnalysisResult::CombinatorialBatchCaptured(commit) => {
            check_workspace_version(bytes, &commit.workspace)?;
            commit.validate_binding(&response.revision)?;
        }
        AnalysisResult::CombinatorialReleaseStarted(progress)
        | AnalysisResult::CombinatorialReleaseAdvanced(progress) => {
            progress.validate()?;
            if progress.context_revision != response.revision {
                return Err(Error::Format(
                    "Combinatorial release response context changed.".into(),
                ));
            }
        }
        AnalysisResult::CombinatorialReleaseInspected(inspection) => {
            inspection.validate()?;
            if inspection.context_revision != response.revision {
                return Err(Error::Format(
                    "Combinatorial inspection response context changed.".into(),
                ));
            }
        }
        AnalysisResult::CombinatorialReleaseFinished(commit) => {
            check_workspace_version(bytes, &commit.workspace)?;
            commit.validate_binding(&response.revision, &commit.source_workspace_fingerprint)?;
        }
        AnalysisResult::CombinatorialOutputsArchive(archive) => {
            archive.validate()?;
            if archive.context_revision != response.revision {
                return Err(Error::Format(
                    "Combinatorial archive response context changed.".into(),
                ));
            }
        }
        AnalysisResult::ReferenceFragment(fragment) => {
            fragment.validate_binding(&response.revision)?;
        }
        AnalysisResult::ReferencePaste(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::TranslatedReference(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::TranslatedOutputs(outputs) => {
            outputs.validate_binding(&response.revision)?;
        }
        AnalysisResult::GoldenGateReaction(report) => {
            check_variant_version(bytes, plan_worker_version(&report.probe_plan))?;
            check_verification_report_version(bytes, &report.verification)?;
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::GoldenGateReactionCandidate(report) => {
            check_variant_version(bytes, plan_worker_version(&report.plan))?;
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::ReadReviewPreview(preview) => {
            preview.validate_binding(&response.revision)?
        }
        AnalysisResult::ReadReviewCommitted(commit) => {
            check_workspace_version(bytes, &commit.workspace)?;
            commit.validate_binding(&response.revision)?;
        }
        AnalysisResult::ReadCheckCaptured(commit) => {
            check_workspace_version(bytes, &commit.workspace)?;
            commit.validate_binding(&response.revision)?;
        }
        AnalysisResult::ReadCheckReplayed(checked) => {
            checked.validate_binding(&response.revision)?
        }
        AnalysisResult::MutationIntent(report) => report.validate_binding(&response.revision)?,
        AnalysisResult::GoldenGatePrimersDesigned(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::CloningOutputsFinalized(commit) => {
            check_workspace_version(bytes, &commit.workspace)?;
            commit.validate_binding(&response.revision)?;
        }
        AnalysisResult::DigestFragmentSelected(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::PlanOligosExported(report) => {
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::AssemblyIntent(report) => report.validate_binding(&response.revision)?,
        AnalysisResult::OverlapPrimersDesigned(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::AssemblyJunctionInspection(report) => {
            check_verification_report_version(bytes, &report.verification)?;
            if matches!(
                report.method,
                Some(crate::assembly_junctions::Method::GoldenGate { .. })
            ) {
                check_variant_version(bytes, 48)?;
            }
            report.validate_binding(&response.revision)?
        }
        AnalysisResult::MismatchPcrDigestInspection(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::ContributionEditIntent(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::SubstitutedAssemblyIntent(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::RunCatalogRestrictionInspection(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::RunPcrDigestInspection(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::RunMaterialInspection(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::PcrDigestInspection(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::TailingRecipePreview(report) => {
            check_verification_report_version(bytes, &report.verification)?;
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::LigationRecipePreview(report) => {
            check_verification_report_version(bytes, &report.verification)?;
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::OverlapRecipePreview(report) => {
            check_verification_report_version(bytes, &report.verification)?;
            report.validate()?;
            if report.context_revision != response.revision {
                return Err(Error::Format(
                    "overlap preparation response context changed".into(),
                ));
            }
        }
        AnalysisResult::VariantsPrepared(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format("variant response context changed".into()));
            }
        }
        AnalysisResult::VariantPage(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format("variant response context changed".into()));
            }
        }
        AnalysisResult::VariantReferenceChecked(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format("variant response context changed".into()));
            }
        }
        AnalysisResult::VariantAllelePreview(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format("variant response context changed".into()));
            }
        }
        AnalysisResult::VariantsExported(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format("variant response context changed".into()));
            }
        }
        AnalysisResult::VariantsDisposed(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format("variant response context changed".into()));
            }
        }
        AnalysisResult::AnnotatedTopoProjection(report) => {
            crate::topo_worker::validate_annotated_projection(report, &response.revision)?
        }
        AnalysisResult::AnnotatedTopoTarget(report) => {
            crate::topo_worker::validate_annotated_target(report, &response.revision)?
        }
        AnalysisResult::TopoSimulation(report) => {
            crate::topo_worker::validate_simulation(report, &response.revision)?
        }
        AnalysisResult::TopoOrientations(report) => {
            crate::topo_worker::validate_orientations(report, &response.revision)?
        }
        AnalysisResult::TopoProjection(report) => {
            crate::topo_worker::validate_projection(report, &response.revision)?
        }
        AnalysisResult::TopoTarget(report) => {
            crate::topo_worker::validate_target(report, &response.revision)?
        }
        AnalysisResult::ReadAlignmentsPrepared(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format(
                    "read alignment preparation response mismatch".into(),
                ));
            }
        }
        AnalysisResult::ReadAlignmentPage(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format(
                    "read alignment page response mismatch".into(),
                ));
            }
        }
        AnalysisResult::ReadAlignmentChecked(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format(
                    "read alignment check response mismatch".into(),
                ));
            }
        }
        AnalysisResult::ReadAlignmentsExported(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format(
                    "read alignment export response mismatch".into(),
                ));
            }
        }
        AnalysisResult::ReadAlignmentsDisposed(report) => {
            report.validate()?;
            if report.revision != response.revision {
                return Err(Error::Format(
                    "read alignment disposal response mismatch".into(),
                ));
            }
        }
        AnalysisResult::CombinatorialMultiSiteInspection(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::MultiSiteRecombination(report) => {
            report.validate_binding(&response.revision)?;
        }
        AnalysisResult::MultiSiteRecombinationProduct(report) => {
            report.validate_binding(&report.selection.batch_revision)?;
            if report.context_revision()? != response.revision {
                return Err(Error::Format(
                    "multisite selected-component reply context mismatch".into(),
                ));
            }
        }
        AnalysisResult::MultiSiteRecombinationChecked(report) => {
            report.validate_binding(&report.product.selection.batch_revision)?;
            if report.context_revision()? != response.revision {
                return Err(Error::Format(
                    "multisite independent-target reply context mismatch".into(),
                ));
            }
        }
        AnalysisResult::RecombinationInspection(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "recombination response binding mismatch".into(),
                ));
            }
        }
        AnalysisResult::RecombinationTargetCheck(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "recombination response binding mismatch".into(),
                ));
            }
        }
        AnalysisResult::RecombinationProduct(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "recombination response binding mismatch".into(),
                ));
            }
        }
        AnalysisResult::CatalogRestrictionInspection(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "catalog restriction batch revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::CatalogDigestProduct(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "catalog product batch revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::CatalogDigestInstalled(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "catalog installation batch revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::DigestInstalled(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "installed digest batch revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::RestrictionInspection(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "restriction inspection batch revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::DigestProduct(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "digest product batch revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::GenomeImportStarted(_)
        | AnalysisResult::GenomeImportProgress(_)
        | AnalysisResult::GenomeFastaFinished(_)
        | AnalysisResult::GenomeRegistered(_)
        | AnalysisResult::GenomeIntervalGuides(_)
        | AnalysisResult::GenomeSearchState(_)
        | AnalysisResult::GenomeSearchPage(_)
        | AnalysisResult::GenomeDisposed(_)
        | AnalysisResult::GenomeIntervalRead(_) => {
            crate::genome_session::validate_result(&response.result, &response.revision)?
        }
        AnalysisResult::InversePcr(report) => {
            check_verification_report_version(bytes, &report.verification)?;
            check_variant_version(bytes, plan_worker_version(&report.plan))?;
        }
        AnalysisResult::Verification(report) => {
            check_verification_report_version(bytes, report)?;
        }
        AnalysisResult::ProtocolPdf(report) => {
            check_scientific_report_version(bytes, &report.engine_version)?;
        }
        AnalysisResult::EditTemplate(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "edit-template batch revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::EditTemplateCheck(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "edit-template check batch revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::CrisprGuides(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format("CRISPR guide batch revision mismatch".into()));
            }
        }
        AnalysisResult::CrisprOffTargets(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format(
                    "CRISPR off-target batch revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::Methylation(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format("methylation batch revision mismatch".into()));
            }
        }
        AnalysisResult::Contigs(report) => {
            report.validate()?;
            if report.batch_revision != response.revision {
                return Err(Error::Format("contig batch revision mismatch".into()));
            }
        }
        AnalysisResult::DotPlot(report) => {
            report.validate()?;
            if report.reference_revision != response.revision {
                return Err(Error::Format("dot plot source revision mismatch".into()));
            }
        }
        AnalysisResult::CodonDesign(report) => {
            report.validate()?;
            if report.source_revision != response.revision {
                return Err(Error::Format(
                    "codon design source revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::CodonCheck(report) => {
            report.validate()?;
            if report.source_revision != response.revision {
                return Err(Error::Format("codon check source revision mismatch".into()));
            }
        }
        AnalysisResult::RunCommit(commit) => {
            check_workspace_version(bytes, &commit.workspace)?;
            check_scientific_report_version(bytes, &commit.receipt.recorded_engine)?;
            commit.validate()?;
            if commit.receipt.source_workspace_fingerprint != response.revision {
                return Err(Error::Format("run commit source workspace mismatch".into()));
            }
        }
        AnalysisResult::RunReplay(report) => {
            check_scientific_report_version(bytes, &report.recorded_engine)?;
            check_verification_report_version(bytes, &report.current_report)?;
            report.validate()?;
            if report.source_workspace_fingerprint != response.revision {
                return Err(Error::Format("run replay source workspace mismatch".into()));
            }
        }
        AnalysisResult::TailPreview(preview) => {
            preview.molecule.validate()?;
            preview.tailing.validate()?;
            if preview.source_revision != response.revision
                || preview.base != preview.tailing.base
                || preview.molecule.sequence != preview.tailing.duplex.top
            {
                return Err(Error::Format(
                    "tailing preview response binding mismatch".into(),
                ));
            }
        }
        AnalysisResult::OverlapProfile(report) => {
            report.validate()?;
            if profile_binding(
                &report.source_revisions,
                &report.orientations,
                &report.overlaps,
                report.circular,
                &report.options,
            )? != response.revision
            {
                return Err(Error::Format(
                    "overlap profile response binding mismatch".into(),
                ));
            }
        }
        AnalysisResult::ProteinProperties(report) => {
            report.validate()?;
            if report.source_revision != response.revision {
                return Err(Error::Format(
                    "protein-property response revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::AlignmentDocument(document) => document.validate()?,
        AnalysisResult::AlignmentRows(rows) => {
            if rows.len() > crate::alignment_formats::MAX_ALIGNMENT_ROWS {
                return Err(Error::TooLarge);
            }
            for row in rows {
                row.validate()?;
            }
        }
        AnalysisResult::AlignmentSaved(bytes) => {
            crate::alignment_formats::load_alignment(bytes)?;
            if blake3::hash(bytes).to_hex().as_str() != response.revision {
                return Err(Error::Format(
                    "saved alignment response revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::AlignmentExport(report) => {
            if report.text.len() > crate::alignment_formats::MAX_ALIGNMENT_BYTES
                || report.losses.len() > crate::alignment_formats::MAX_ALIGNMENT_ANNOTATIONS
                || report.renamings.len() > crate::alignment_formats::MAX_ALIGNMENT_ROWS
            {
                return Err(Error::TooLarge);
            }
        }
        AnalysisResult::TracksPrepared(revision) => {
            if revision != &response.revision
                || revision.len() != 64
                || !revision.bytes().all(|b| b.is_ascii_hexdigit())
            {
                return Err(Error::Format(
                    "prepared-track response revision mismatch".into(),
                ));
            }
        }
        AnalysisResult::DocumentSearch(report) => {
            report.validate()?;
            if report.source_revision != response.revision {
                return Err(Error::Format("search response revision mismatch".into()));
            }
        }
        AnalysisResult::ViewTracks(report) => {
            report.validate()?;
            if report.source_revision != response.revision {
                return Err(Error::Format("track response revision mismatch".into()));
            }
        }
        AnalysisResult::RepairPreview(preview) => {
            preview.molecule.validate()?;
            preview.repair.source.validate()?;
            preview.repair.duplex.validate()?;
            if preview.molecule.sequence != preview.repair.duplex.top {
                return Err(Error::Format(
                    "repair preview sequence/duplex mismatch".into(),
                ));
            }
        }
        AnalysisResult::AnnealPreview(preview) => {
            for candidate in &preview.annealing.candidates {
                candidate.duplex.validate()?;
            }
        }
        AnalysisResult::AnnotatedDocument(molecule)
        | AnalysisResult::PhysicalPreviewDocument(molecule) => molecule.validate()?,
        AnalysisResult::WorkspaceLoaded(workspace) => workspace.validate()?,
        AnalysisResult::WorkspaceSaved(data) => {
            if data.len() > crate::workspace::MAX_ARCHIVE_BYTES {
                return Err(Error::TooLarge);
            }
            if u16::from_le_bytes([bytes[4], bytes[5]]) < WORKER_SCHEMA_VERSION {
                let workspace = crate::workspace::load_workspace(data)?;
                check_workspace_version(bytes, &workspace)?;
            }
        }
        _ => {}
    }
    Ok(response)
}
// Probe the borrowed selector before decoding a schema-13 scientific payload.
// These new jobs/results are never eligible for the larger workspace budget.
fn preflight_codon_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 13)?;
    }
    Ok(())
}
fn preflight_fusion_design_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 60)?;
    }
    Ok(())
}
fn preflight_domestication_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 61)?;
    }
    Ok(())
}
fn preflight_topo_actions_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 66)?;
    }
    Ok(())
}
fn preflight_processed_run_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 65)?;
    }
    Ok(())
}
fn preflight_processed_topo_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 64)?;
    }
    Ok(())
}
fn preflight_nick_sealing_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 63)?;
    }
    Ok(())
}
fn preflight_published_ligation_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 62)?;
    }
    Ok(())
}
// Dot plots retain the ordinary scientific budget; gate the new selector before
// deserializing any schema-14 source/query body under a historical envelope.
fn preflight_dotplot_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 14)?;
    }
    Ok(())
}
// Assembly retains ordinary scientific limits and requires its explicit schema
// before any read batch is deserialized under a historical envelope.
fn preflight_contig_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 15)?;
    }
    Ok(())
}
fn preflight_methylation_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 16)?;
    }
    Ok(())
}
fn preflight_topo_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 26)?;
    }
    Ok(())
}
fn preflight_pinned_gel_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 55)?;
    }
    Ok(())
}
fn preflight_pairwise_export_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 56)?;
    }
    Ok(())
}
fn preflight_combinatorial_batch_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 57)?;
    }
    Ok(())
}
fn preflight_combinatorial_release_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 58)?;
    }
    Ok(())
}
fn preflight_translation_reference_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 54)?;
    }
    Ok(())
}
fn preflight_reference_transfer_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 53)?;
    }
    Ok(())
}
fn preflight_pairwise_review_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 52)?;
    }
    Ok(())
}
fn preflight_mismatch_pcr_outputs_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 51)?;
    }
    Ok(())
}
fn preflight_pcr_outputs_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 50)?;
    }
    Ok(())
}
fn preflight_golden_gate_reaction_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 49)?;
    }
    Ok(())
}
fn preflight_prepared_golden_gate_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 47)?;
    }
    Ok(())
}
fn preflight_golden_gate_junction_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 46)?;
    }
    Ok(())
}
fn preflight_selective_assembly_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 45)?;
    }
    Ok(())
}
fn preflight_edited_assembly_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 44)?;
    }
    Ok(())
}
fn preflight_read_review_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 43)?;
    }
    Ok(())
}
fn preflight_read_check_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 42)?;
    }
    Ok(())
}
fn preflight_mutation_intent_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 41)?;
    }
    Ok(())
}
fn preflight_golden_gate_primer_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 40)?;
    }
    Ok(())
}
fn preflight_finalize_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 39)?;
    }
    Ok(())
}
fn preflight_cloning_output_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 38)?;
    }
    Ok(())
}
fn preflight_assembly_intent_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 36)?;
    }
    Ok(())
}
fn preflight_overlap_primer_design_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 35)?;
    }
    Ok(())
}
fn preflight_assembly_junctions_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 33)?;
    }
    Ok(())
}
fn preflight_mismatch_digest_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 78)?;
    }
    Ok(())
}
fn preflight_contribution_edit_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 79)?;
    }
    Ok(())
}
fn preflight_substitution_intent_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 77)?;
    }
    Ok(())
}
fn preflight_run_catalog_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 75)?;
    }
    Ok(())
}
fn preflight_run_material_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 74)?;
    }
    Ok(())
}
fn preflight_pcr_digest_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 72)?;
    }
    Ok(())
}
fn preflight_tailing_recipe_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 70)?;
    }
    Ok(())
}
fn preflight_ligation_recipe_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 32)?;
    }
    Ok(())
}
fn preflight_overlap_recipe_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 31)?;
    }
    Ok(())
}
fn preflight_variants_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 29)?;
    }
    Ok(())
}
fn preflight_annotated_topo_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 27)?;
    }
    Ok(())
}
fn preflight_read_alignment_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 24)?;
    }
    Ok(())
}
fn preflight_combinatorial_multisite_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 69)?;
    }
    Ok(())
}
fn preflight_multisite_recombination_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 68)?;
    }
    Ok(())
}
fn preflight_recombination_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 22)?;
    }
    Ok(())
}
fn preflight_catalog_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 21)?;
    }
    Ok(())
}
fn preflight_restriction_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 20)?;
    }
    Ok(())
}
fn preflight_genome_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 19)?;
    }
    Ok(())
}
fn preflight_crispr_tag(bytes: &[u8], tags: &[u32]) -> Result<(), Error> {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tags.contains(&tag) {
        check_variant_version(bytes, 18)?;
    }
    Ok(())
}
fn check_verification_report_version(
    bytes: &[u8],
    report: &crate::verification::VerificationReport,
) -> Result<(), Error> {
    check_scientific_report_version(bytes, &report.engine_version)?;
    if !report.reference_step_reports.is_empty() {
        check_variant_version(bytes, 68)?;
        if report.reference_step_reports.len() > 64 {
            return Err(Error::TooLarge);
        }
        let mut seen = std::collections::BTreeSet::new();
        for step in &report.reference_step_reports {
            if step.step_index >= 1000 || !seen.insert(step.step_index) {
                return Err(Error::Format("Invalid reference step index".into()));
            }
            match &step.evidence {
                crate::verification::ReferenceStepEvidence::MultiSiteRecombination {
                    selection,
                    batch,
                } => {
                    batch.validate_binding(&selection.batch_revision)?;
                    if !batch.documents.iter().any(|p| p.selection == *selection) {
                        return Err(Error::Format(
                            "Reference selection is missing from its batch".into(),
                        ));
                    }
                }
            }
        }
    }
    if !report.physical_processing_reports.is_empty() {
        check_variant_version(bytes, 66)?;
        if report.physical_processing_reports.len() > 64 {
            return Err(Error::TooLarge);
        }
        let mut seen = std::collections::BTreeSet::new();
        for step in &report.physical_processing_reports {
            if step.step_index >= 1000 || !seen.insert(step.step_index) {
                return Err(Error::Format(
                    "Invalid physical processing step index".into(),
                ));
            }
            match &step.evidence {
                crate::verification::PhysicalProcessingEvidence::RunProduct(r) => {
                    check_variant_version(bytes, 73)?;
                    r.validate_documentary()?;
                    check_capsule_candidate_version(
                        u16::from_le_bytes([bytes[4], bytes[5]]),
                        &r.origin.candidate,
                    )?;
                    // Inspect the retained scientific stamp without replaying or
                    // conferring freshness on this documentary report.
                    check_scientific_report_version(bytes, &r.decode_report()?.engine_version)?;
                }
                crate::verification::PhysicalProcessingEvidence::TopoActions(r) => {
                    r.validate_binding(&r.context_revision)?
                }
                crate::verification::PhysicalProcessingEvidence::OverlapExtension(r) => {
                    check_variant_version(bytes, 81)?;
                    r.validate_documentary()?;
                }
            }
        }
    }
    if !report.processed_step_reports.is_empty() {
        check_variant_version(bytes, 65)?;
        if report.processed_step_reports.len() > 64 {
            return Err(Error::TooLarge);
        }
        let mut seen = std::collections::BTreeSet::new();
        for step in &report.processed_step_reports {
            if step.step_index >= 1000 || !seen.insert(step.step_index) {
                return Err(Error::Format(
                    "Invalid processed step evidence index".into(),
                ));
            }
            step.evidence
                .validate_binding(&step.evidence.context_revision)?;
        }
    }
    Ok(())
}
fn check_scientific_report_version(bytes: &[u8], engine: &str) -> Result<(), Error> {
    if engine.ends_with("/science-22") {
        check_variant_version(bytes, 81)?;
    }
    if engine.ends_with("/science-21") {
        check_variant_version(bytes, 80)?;
    }
    if engine.ends_with("/science-20") {
        check_variant_version(bytes, 73)?;
    }
    if engine.ends_with("/science-19") {
        check_variant_version(bytes, 68)?;
    }
    if engine.ends_with("/science-18") {
        check_variant_version(bytes, 67)?;
    }
    if engine.ends_with("/science-17") {
        check_variant_version(bytes, 66)?;
    }
    if engine.ends_with("/science-16") {
        check_variant_version(bytes, 65)?;
    }
    if engine.ends_with("/science-15") {
        check_variant_version(bytes, 59)?;
    }
    if engine.ends_with("/science-14") {
        check_variant_version(bytes, 37)?;
    }
    if engine.ends_with("/science-13") {
        check_variant_version(bytes, 34)?;
    }
    if engine.ends_with("/science-12") {
        check_variant_version(bytes, 28)?;
    }
    if engine.ends_with("/science-11") {
        check_variant_version(bytes, 27)?;
    }
    if engine.ends_with("/science-10") {
        check_variant_version(bytes, 25)?;
    }
    if engine.ends_with("/science-9") {
        check_variant_version(bytes, 22)?;
    }
    if engine.ends_with("/science-8") {
        check_variant_version(bytes, 17)?;
    }
    Ok(())
}
fn annealed_preparation(
    plan: &crate::verification::DesignPlan,
    row: &crate::assembly_recipe::PreparedRow,
) -> bool {
    row.preparation_step
        .and_then(|i| plan.steps.get(i))
        .is_some_and(|step| matches!(step, crate::verification::Operation::AnnealOligos { .. }))
}
// A post-decode semantic gate for an appended nested enum. This visits borrowed
// typed shapes without allocating an intermediate JSON value or copying DNA.
// Byte payloads remain opaque; recognized workspace records have their own
// decoded typed request checks below.
mod preparation_shape {
    use serde::Serialize;
    use serde::ser::{
        self, SerializeMap, SerializeSeq, SerializeStruct, SerializeStructVariant, SerializeTuple,
        SerializeTupleStruct, SerializeTupleVariant,
    };
    #[derive(Debug)]
    pub struct Found(u16);
    impl std::fmt::Display for Found {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(
                f,
                "recipe preparation capability requires worker schema {}",
                self.0
            )
        }
    }
    impl std::error::Error for Found {}
    impl ser::Error for Found {
        fn custom<T: std::fmt::Display>(_: T) -> Self {
            Self(74)
        }
    }
    pub struct Visitor {
        pub worker: u16,
        pub candidate: bool,
    }
    fn variant(worker: u16, name: &str, variant: &str) -> Result<(), Found> {
        let minimum = match (name, variant) {
            ("Preparation", "OverlapExtensionPcr")
            | ("RunPreparation", "OverlapExtensionPcr")
            | ("Stage", "OverlapExtension") => 82,
            ("Operation", "OverlapExtension")
            | ("PhysicalProcessingEvidence", "OverlapExtension") => 81,
            ("Preparation", "AnnealedOligos") => 71,
            ("Preparation", "FromRun") => 74,
            ("Preparation", "MismatchPcr") | ("RunPreparation", "MismatchPcr") => 76,
            ("Preparation", "MismatchPcrCatalogDigest")
            | ("RunPreparation", "MismatchPcrCatalogDigest") => 78,
            _ => 0,
        };
        if worker < minimum {
            Err(Found(minimum))
        } else {
            Ok(())
        }
    }
    macro_rules! scalar {
        ($($method:ident($ty:ty)),* $(,)?) => { $(fn $method(self, _: $ty) -> Result<(),Found> { Ok(()) })* };
    }
    impl ser::Serializer for Visitor {
        fn is_human_readable(&self) -> bool {
            false
        }
        type Ok = ();
        type Error = Found;
        type SerializeSeq = Self;
        type SerializeTuple = Self;
        type SerializeTupleStruct = Self;
        type SerializeTupleVariant = Self;
        type SerializeMap = Self;
        type SerializeStruct = Self;
        type SerializeStructVariant = Self;
        scalar!(
            serialize_bool(bool),
            serialize_i8(i8),
            serialize_i16(i16),
            serialize_i32(i32),
            serialize_i64(i64),
            serialize_i128(i128),
            serialize_u8(u8),
            serialize_u16(u16),
            serialize_u32(u32),
            serialize_u64(u64),
            serialize_u128(u128),
            serialize_f32(f32),
            serialize_f64(f64),
            serialize_char(char),
            serialize_str(&str)
        );
        fn serialize_bytes(self, bytes: &[u8]) -> Result<(), Found> {
            if self.candidate {
                super::check_capsule_candidate_version(self.worker, bytes)
                    .map_err(|_| Found(81))?;
            }
            Ok(())
        }
        fn serialize_none(self) -> Result<(), Found> {
            Ok(())
        }
        fn serialize_some<T: ?Sized + Serialize>(self, v: &T) -> Result<(), Found> {
            v.serialize(Self {
                worker: self.worker,
                candidate: false,
            })
        }
        fn serialize_unit(self) -> Result<(), Found> {
            Ok(())
        }
        fn serialize_unit_struct(self, _: &'static str) -> Result<(), Found> {
            Ok(())
        }
        fn serialize_unit_variant(
            self,
            n: &'static str,
            _: u32,
            v: &'static str,
        ) -> Result<(), Found> {
            variant(self.worker, n, v)
        }
        fn serialize_newtype_struct<T: ?Sized + Serialize>(
            self,
            _: &'static str,
            v: &T,
        ) -> Result<(), Found> {
            v.serialize(Self {
                worker: self.worker,
                candidate: false,
            })
        }
        fn serialize_newtype_variant<T: ?Sized + Serialize>(
            self,
            n: &'static str,
            _: u32,
            v: &'static str,
            value: &T,
        ) -> Result<(), Found> {
            variant(self.worker, n, v)?;
            // Native archive bytes are inspected through check_workspace_version,
            // not traversed one scalar byte at a time by this shape-only visitor.
            if n == "AnalysisResult" && v == "WorkspaceSaved" {
                return Ok(());
            }
            value.serialize(Self {
                worker: self.worker,
                candidate: false,
            })
        }
        fn serialize_seq(self, _: Option<usize>) -> Result<Self, Found> {
            Ok(self)
        }
        fn serialize_tuple(self, _: usize) -> Result<Self, Found> {
            Ok(self)
        }
        fn serialize_tuple_struct(self, _: &'static str, _: usize) -> Result<Self, Found> {
            Ok(self)
        }
        fn serialize_tuple_variant(
            self,
            n: &'static str,
            _: u32,
            v: &'static str,
            _: usize,
        ) -> Result<Self, Found> {
            variant(self.worker, n, v)?;
            Ok(self)
        }
        fn serialize_map(self, _: Option<usize>) -> Result<Self, Found> {
            Ok(self)
        }
        fn serialize_struct(self, _: &'static str, _: usize) -> Result<Self, Found> {
            Ok(self)
        }
        fn serialize_struct_variant(
            self,
            n: &'static str,
            _: u32,
            v: &'static str,
            _: usize,
        ) -> Result<Self, Found> {
            variant(self.worker, n, v)?;
            Ok(self)
        }
        fn collect_str<T: ?Sized + std::fmt::Display>(self, _: &T) -> Result<(), Found> {
            Ok(())
        }
    }
    macro_rules! sequence {
        ($trait:ident,$method:ident) => {
            impl $trait for Visitor {
                type Ok = ();
                type Error = Found;
                fn $method<T: ?Sized + Serialize>(&mut self, v: &T) -> Result<(), Found> {
                    v.serialize(Visitor {
                        worker: self.worker,
                        candidate: false,
                    })
                }
                fn end(self) -> Result<(), Found> {
                    Ok(())
                }
            }
        };
    }
    sequence!(SerializeSeq, serialize_element);
    sequence!(SerializeTuple, serialize_element);
    sequence!(SerializeTupleStruct, serialize_field);
    sequence!(SerializeTupleVariant, serialize_field);
    impl SerializeMap for Visitor {
        type Ok = ();
        type Error = Found;
        fn serialize_key<T: ?Sized + Serialize>(&mut self, v: &T) -> Result<(), Found> {
            v.serialize(Visitor {
                worker: self.worker,
                candidate: false,
            })
        }
        fn serialize_value<T: ?Sized + Serialize>(&mut self, v: &T) -> Result<(), Found> {
            v.serialize(Visitor {
                worker: self.worker,
                candidate: false,
            })
        }
        fn end(self) -> Result<(), Found> {
            Ok(())
        }
    }
    macro_rules! structure {
        ($trait:ident) => {
            impl $trait for Visitor {
                type Ok = ();
                type Error = Found;
                fn serialize_field<T: ?Sized + Serialize>(
                    &mut self,
                    key: &'static str,
                    v: &T,
                ) -> Result<(), Found> {
                    if key == "bytes" {
                        Ok(())
                    } else {
                        v.serialize(Visitor {
                            worker: self.worker,
                            candidate: key == "candidate",
                        })
                    }
                }
                fn end(self) -> Result<(), Found> {
                    Ok(())
                }
            }
        };
    }
    structure!(SerializeStruct);
    structure!(SerializeStructVariant);
}
fn check_preparation_version(bytes: &[u8], value: &impl Serialize) -> Result<(), Error> {
    let worker = u16::from_le_bytes([bytes[4], bytes[5]]);
    if worker < 82 {
        value
            .serialize(preparation_shape::Visitor {
                worker,
                candidate: false,
            })
            .map_err(|e| Error::Format(e.to_string()))?;
    }
    Ok(())
}

fn check_combinatorial_request_version(
    bytes: &[u8],
    request: &crate::combinatorial_batch::Request,
) -> Result<(), Error> {
    check_preparation_version(bytes, request)?;
    use crate::combinatorial_design::Method;
    let required = |method: &Method| match method {
        Method::MultiSiteReference { .. } => 69,
        Method::GoldenGatePool { .. } => 59,
        _ => 57,
    };
    if let Some(design) = &request.design {
        check_variant_version(bytes, required(&design.method))?;
    }
    for row in &request.rows {
        if let crate::combinatorial_batch::RowRequest::Ready(row) = row
            && let Some(method) = row.overrides.as_ref().and_then(|o| o.method.as_ref())
        {
            check_variant_version(bytes, required(method))?;
        }
    }
    Ok(())
}
fn check_workspace_version(
    bytes: &[u8],
    workspace: &crate::workspace::Workspace,
) -> Result<(), Error> {
    if u16::from_le_bytes([bytes[4], bytes[5]]) >= WORKER_SCHEMA_VERSION {
        return Ok(());
    }
    for candidate in &workspace.candidates {
        check_variant_version(bytes, plan_worker_version(&candidate.plan))?;
        check_preparation_version(bytes, &candidate.plan)?;
    }
    for attachment in workspace.attachments.values() {
        if attachment.kind == "assembly-draft-v8" {
            check_variant_version(bytes, 82)?;
        } else if attachment.kind == "assembly-draft-v7" {
            check_variant_version(bytes, 79)?;
        } else if attachment.kind == "assembly-draft-v6" {
            check_variant_version(bytes, 77)?;
        } else if attachment.kind == "assembly-draft-v5" {
            check_variant_version(bytes, 74)?;
        } else if attachment.kind == "assembly-draft-v4" {
            check_variant_version(bytes, 71)?;
        } else if attachment.kind == crate::run_lineage::WORKSPACE_ATTACHMENT_KIND {
            check_variant_version(bytes, 12)?;
            let record = crate::run_lineage::load_run_record(&attachment.bytes)?;
            check_variant_version(bytes, plan_worker_version(&record.plan))?;
            check_preparation_version(bytes, &record.plan)?;
            if record.plan_schema >= 16 {
                check_variant_version(bytes, 81)?;
            }
            if record.plan_schema >= 15 {
                check_variant_version(bytes, 73)?;
            }
            if record.plan_schema >= 14 {
                check_variant_version(bytes, 68)?;
            }
            if record.plan_schema >= 13 {
                check_variant_version(bytes, 67)?;
            }
            if record.plan_schema >= 12 {
                check_variant_version(bytes, 66)?;
            }
            if record.plan_schema >= 11 {
                check_variant_version(bytes, 65)?;
            }
            if record.plan_schema >= 10 {
                check_variant_version(bytes, 59)?;
            }
            if record.plan_schema >= 4 {
                check_variant_version(bytes, 17)?;
            }
            if record.plan_schema >= 9 {
                check_variant_version(bytes, 37)?;
            }
            if record.plan_schema >= 8 {
                check_variant_version(bytes, 34)?;
            }
            if record.plan_schema >= 7 {
                check_variant_version(bytes, 28)?;
            }
            if record.plan_schema >= 6 {
                check_variant_version(bytes, 25)?;
            }
            if record.plan_schema >= 5 {
                check_variant_version(bytes, 22)?;
            }
            check_scientific_report_version(bytes, &record.assessment.engine_version)?;
        } else if attachment.kind == crate::combinatorial_batch::KIND {
            check_variant_version(bytes, 57)?;
            let record = crate::combinatorial_batch::load(
                workspace,
                &blake3::hash(&attachment.bytes).to_hex().to_string(),
            )?;
            check_combinatorial_request_version(bytes, &record.request)?;
            check_scientific_report_version(bytes, &record.engine_version)?;
        } else if attachment.kind == crate::combinatorial_release::KIND {
            check_variant_version(bytes, 58)?;
            let manifest = crate::combinatorial_release::load(
                workspace,
                &blake3::hash(&attachment.bytes).to_hex().to_string(),
            )?;
            check_scientific_report_version(bytes, &manifest.engine_version)?;
        } else if attachment.kind == crate::genome::WORKSPACE_ATTACHMENT_KIND {
            check_variant_version(bytes, 19)?;
        } else if attachment.kind == "read-assembly-v1" {
            check_variant_version(bytes, 15)?;
        }
    }
    Ok(())
}
// Closed run origins retain typed plans in compact candidate envelopes. Walk
// those packages iteratively so a new operation cannot hide under an old worker
// envelope. This is capability validation only, never scientific replay.
fn check_capsule_candidate_version(worker: u16, first: &[u8]) -> Result<(), Error> {
    if first.len() > crate::run_material_capsule::MAX_CAPSULE_BYTES {
        return Err(Error::TooLarge);
    }
    let mut pending = vec![(first.to_vec(), 0usize)];
    let mut total = 0usize;
    while let Some((bytes, depth)) = pending.pop() {
        total = total.checked_add(bytes.len()).ok_or(Error::TooLarge)?;
        if total > crate::run_material_capsule::MAX_CAPSULE_BYTES
            || depth > crate::run_material_capsule::MAX_NESTING
        {
            return Err(Error::TooLarge);
        }
        let candidate = crate::transfer::decode_candidate_shallow(&bytes)?;
        if worker < plan_worker_version(&candidate.plan) {
            return Err(Error::Unsupported(
                "nested origin requires a newer worker schema".into(),
            ));
        }
        for operation in candidate.plan.steps {
            if let crate::verification::Operation::ReconstructRunProduct { origin, .. } = operation
            {
                pending.push((origin.candidate, depth + 1));
            }
        }
    }
    Ok(())
}
// AnnotationFormat is append-only inside existing jobs/reports. A new selector
// must not acquire semantics under an older envelope, even when the outer tag
// predates it. Historical GFF3/BED payloads retain their schema-7 contract.
fn check_annotation_format_version(
    bytes: &[u8],
    format: crate::annotation_formats::AnnotationFormat,
) -> Result<(), Error> {
    if format == crate::annotation_formats::AnnotationFormat::Gtf {
        check_variant_version(bytes, 30)?;
    }
    Ok(())
}
fn check_variant_version(bytes: &[u8], required: u16) -> Result<(), Error> {
    if u16::from_le_bytes([bytes[4], bytes[5]]) < required {
        return Err(Error::Unsupported(format!(
            "message variant requires binary schema {required} or later"
        )));
    }
    Ok(())
}
fn profile_binding(
    revisions: &[String],
    reverse: &[bool],
    overlaps: &[usize],
    circular: bool,
    options: &crate::overlap_profiles::ProfileOptions,
) -> Result<String, Error> {
    if revisions.is_empty()
        || revisions.len() > 1000
        || reverse.len() != revisions.len()
        || overlaps.len() > 1000
    {
        return Err(Error::Format(
            "profile binding needs explicit bounded source/orientation lists".into(),
        ));
    }
    let bytes = serde_json::to_vec(&(revisions, reverse, overlaps, circular, options))
        .map_err(|e| Error::Format(e.to_string()))?;
    Ok(blake3::hash(&bytes).to_hex().to_string())
}
pub fn overlap_profile_revision(
    molecules: &[Molecule],
    reverse: &[bool],
    overlaps: &[usize],
    circular: bool,
    options: &crate::overlap_profiles::ProfileOptions,
) -> Result<String, Error> {
    if molecules.is_empty() || molecules.len() > 1000 {
        return Err(Error::TooLarge);
    }
    profile_binding(
        &molecules.iter().map(Molecule::revision).collect::<Vec<_>>(),
        reverse,
        overlaps,
        circular,
        options,
    )
}
fn plan_worker_version(plan: &crate::verification::DesignPlan) -> u16 {
    plan.steps
        .iter()
        .map(|step| match step {
            crate::verification::Operation::OverlapExtension { .. } => 81,
            crate::verification::Operation::ReconstructRunProduct { .. } => 73,
            crate::verification::Operation::MultiSiteRecombineReference { .. } => 68,
            step if crate::verification::topo_uses_step_source(step) => 67,
            crate::verification::Operation::TopoActions { .. } => 66,
            crate::verification::Operation::ProcessedTopo { .. } => 65,
            crate::verification::Operation::SelectedGoldenGate { .. } => 59,
            crate::verification::Operation::CatalogDigest { .. } => 37,
            crate::verification::Operation::DeclareSourceTermini { .. } => 34,
            crate::verification::Operation::CheckedTopo { .. } => 28,
            crate::verification::Operation::Topo { .. } => 25,
            crate::verification::Operation::RecombineReference { .. } => 22,
            crate::verification::Operation::DeclareSourceModifications { .. }
            | crate::verification::Operation::CheckedDigest { .. } => 17,
            crate::verification::Operation::Tail3Prime { .. }
            | crate::verification::Operation::ProfiledOverlap { .. } => 11,
            crate::verification::Operation::RepairEnds { .. }
            | crate::verification::Operation::AnnealOligos { .. } => 8,
            crate::verification::Operation::PcrWithMismatches { .. }
            | crate::verification::Operation::GoldenGate { .. } => 5,
            crate::verification::Operation::Digest { .. }
            | crate::verification::Operation::TreatEnds { .. }
            | crate::verification::Operation::Ligate { .. } => 2,
            _ => 1,
        })
        .max()
        .unwrap_or(1)
}
pub fn process_request(bytes: &[u8]) -> Result<Vec<u8>, Error> {
    let maximum = worker_limit(
        b"CGQ\0",
        bytes,
        &[
            (19, 6),
            (20, 6),
            (42, 12),
            (43, 12),
            (69, 20),
            (72, 21),
            (102, 39),
            (105, 42),
            (106, 42),
            (107, 43),
            (108, 43),
            (128, 57),
            (130, 57),
            (131, 58),
            (134, 58),
            (135, 58),
            (149, 65),
            (155, 69),
            (158, 74),
        ],
    )?;
    preflight_codon_tag(bytes, &[44, 45])?;
    preflight_dotplot_tag(bytes, &[46])?;
    preflight_contig_tag(bytes, &[47])?;
    preflight_methylation_tag(bytes, &[48])?;
    preflight_crispr_tag(bytes, &[49, 50, 51, 52])?;
    preflight_genome_tag(
        bytes,
        &[53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66],
    )?;
    preflight_restriction_tag(bytes, &[67, 68, 69])?;
    preflight_catalog_tag(bytes, &[70, 71, 72])?;
    preflight_recombination_tag(bytes, &[73, 74, 75])?;
    preflight_read_alignment_tag(bytes, &[77, 78, 79, 80, 81])?;
    preflight_topo_tag(bytes, &[82, 83, 84, 85])?;
    preflight_annotated_topo_tag(bytes, &[86, 87])?;
    preflight_variants_tag(bytes, &[88, 89, 90, 91, 92, 93])?;
    preflight_overlap_recipe_tag(bytes, &[94])?;
    preflight_ligation_recipe_tag(bytes, &[95])?;
    preflight_assembly_junctions_tag(bytes, &[96])?;
    preflight_overlap_primer_design_tag(bytes, &[97])?;
    preflight_assembly_intent_tag(bytes, &[98])?;
    preflight_cloning_output_tag(bytes, &[99, 100, 101])?;
    preflight_finalize_tag(bytes, &[102])?;
    preflight_golden_gate_primer_tag(bytes, &[103])?;
    preflight_mutation_intent_tag(bytes, &[104])?;
    preflight_read_check_tag(bytes, &[105, 106])?;
    preflight_read_review_tag(bytes, &[107, 108])?;
    preflight_edited_assembly_tag(bytes, &[109, 110])?;
    preflight_selective_assembly_tag(bytes, &[111])?;
    preflight_golden_gate_junction_tag(bytes, &[112])?;
    preflight_prepared_golden_gate_tag(bytes, &[113, 114])?;
    preflight_golden_gate_reaction_tag(bytes, &[115, 116])?;
    preflight_pcr_outputs_tag(bytes, &[117])?;
    preflight_mismatch_pcr_outputs_tag(bytes, &[118])?;
    preflight_pairwise_review_tag(bytes, &[119])?;
    preflight_reference_transfer_tag(bytes, &[120, 121])?;
    preflight_translation_reference_tag(bytes, &[122, 123])?;
    preflight_pinned_gel_tag(bytes, &[124])?;
    preflight_pairwise_export_tag(bytes, &[125, 126, 127])?;
    preflight_combinatorial_batch_tag(bytes, &[128, 129, 130])?;
    preflight_combinatorial_release_tag(bytes, &[131, 132, 133, 134, 135])?;
    preflight_fusion_design_tag(bytes, &[136, 137, 138])?;
    preflight_domestication_tag(bytes, &[139, 140, 141, 142])?;
    preflight_published_ligation_tag(bytes, &[143, 144, 145, 146])?;
    preflight_nick_sealing_tag(bytes, &[147])?;
    preflight_processed_topo_tag(bytes, &[148])?;
    preflight_processed_run_tag(bytes, &[149])?;
    preflight_topo_actions_tag(bytes, &[150, 151])?;
    preflight_multisite_recombination_tag(bytes, &[152, 153, 154])?;
    preflight_combinatorial_multisite_tag(bytes, &[155])?;
    preflight_tailing_recipe_tag(bytes, &[156])?;
    preflight_pcr_digest_tag(bytes, &[157])?;
    preflight_run_material_tag(bytes, &[158, 159])?;
    preflight_run_catalog_tag(bytes, &[160])?;
    preflight_substitution_intent_tag(bytes, &[161])?;
    preflight_mismatch_digest_tag(bytes, &[162])?;
    preflight_contribution_edit_tag(bytes, &[163])?;
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tag == 76 {
        check_variant_version(bytes, 23)?;
    }
    let r: Request = decode_with_limit(b"CGQ\0", bytes, maximum)?;
    let required = match &r.job {
        AnalysisJob::BuildContributionEditIntent { .. } => 79,
        AnalysisJob::InspectMismatchPcrDigest { .. } => 78,
        AnalysisJob::BuildSubstitutedAssemblyIntent { .. } => 77,
        AnalysisJob::InspectRunCatalogRestriction { .. } => 75,
        AnalysisJob::InspectRunMaterial { .. } | AnalysisJob::InspectRunPcrDigest { .. } => 74,
        AnalysisJob::InspectPcrDigest { .. } => 72,
        AnalysisJob::PreviewTailingRecipe { .. } => 70,
        AnalysisJob::InspectCombinatorialMultiSite { .. } => 69,
        AnalysisJob::SimulateMultiSiteRecombination { .. }
        | AnalysisJob::ProjectMultiSiteRecombination { .. }
        | AnalysisJob::CheckMultiSiteRecombination { .. } => 68,
        AnalysisJob::ProcessTopoActions { .. } | AnalysisJob::InspectTopoActions { .. } => 66,
        AnalysisJob::LoadProcessedRunEvidence { .. } => 65,
        AnalysisJob::ProjectProcessedTopo { .. } => 64,
        AnalysisJob::SealTopoNicks { .. } => 63,
        AnalysisJob::ListLigationProfiles
        | AnalysisJob::ScorePublishedLigation { .. }
        | AnalysisJob::RankGoldenGateFusions { .. }
        | AnalysisJob::ScoreLigationCompetition { .. } => 62,
        AnalysisJob::InspectDomestication { .. }
        | AnalysisJob::ProposeDomestication { .. }
        | AnalysisJob::CheckDomestication { .. }
        | AnalysisJob::PrepareDomestication { .. } => 61,
        AnalysisJob::SearchGoldenGateFusions { .. }
        | AnalysisJob::MaterializeGoldenGateFusions { .. }
        | AnalysisJob::ScoreLigationFidelity { .. } => 60,
        AnalysisJob::PreparePcrOutputs { .. } => 50,
        AnalysisJob::PrepareMismatchPcrOutputs { .. } => 51,
        AnalysisJob::CompareSequences { .. } => 52,
        AnalysisJob::CopyReferenceFragment { .. } => 53,
        AnalysisJob::PasteReferenceFragment { .. } => 53,
        AnalysisJob::TranslateReference { .. } => 54,
        AnalysisJob::PrepareTranslatedReference { .. } => 54,
        AnalysisJob::PredictPinnedGel { .. } => 55,
        AnalysisJob::ExportPairwiseComparison { .. } => 56,
        AnalysisJob::SearchPairwiseComparison { .. } => 56,
        AnalysisJob::PreparePairwiseContext { .. } => 56,
        AnalysisJob::BeginCombinatorialBatch { .. }
        | AnalysisJob::AdvanceCombinatorialBatch { .. }
        | AnalysisJob::CaptureCombinatorialBatch { .. } => 57,
        AnalysisJob::BeginCombinatorialRelease { .. }
        | AnalysisJob::AdvanceCombinatorialRelease { .. }
        | AnalysisJob::InspectCombinatorialRelease { .. }
        | AnalysisJob::FinishCombinatorialRelease { .. }
        | AnalysisJob::ExportCombinatorialOutputs { .. } => 58,
        AnalysisJob::PredictGel { lanes, .. }
            if lanes.iter().any(|lane| {
                matches!(
                    lane.preparation,
                    crate::gel::Preparation::CatalogDigest { .. }
                )
            }) =>
        {
            55
        }
        AnalysisJob::InspectGoldenGateReaction { .. }
        | AnalysisJob::SelectGoldenGateReaction { .. } => 49,
        AnalysisJob::BuildEditedAssemblyIntent { .. }
        | AnalysisJob::DesignEditedAssemblyPrimers { .. } => 44,
        AnalysisJob::DesignSelectiveAssemblyPrimers { .. } => 45,
        AnalysisJob::DesignGoldenGateJunctionPrimers { .. } => 46,
        AnalysisJob::InspectPreparedGoldenGatePool { .. }
        | AnalysisJob::SelectPreparedGoldenGateCandidate { .. } => 47,
        AnalysisJob::PreviewReadReview { .. } | AnalysisJob::CommitReadReview { .. } => 43,
        AnalysisJob::CaptureReadCheck { .. } | AnalysisJob::ReplayReadCheck { .. } => 42,
        AnalysisJob::BuildMutationIntent { .. } => 41,
        AnalysisJob::DesignGoldenGatePrimers { .. } => 40,
        AnalysisJob::PreviewOverlapRecipe { recipe, .. }
            if recipe.fragments.iter().any(|row| {
                matches!(
                    row.preparation,
                    crate::assembly_recipe::Preparation::PcrCatalogDigest { .. }
                )
            }) =>
        {
            40
        }
        AnalysisJob::PreviewLigationRecipe { recipe, .. }
            if recipe.fragments.iter().any(|row| {
                matches!(
                    row.source.preparation,
                    crate::assembly_recipe::Preparation::PcrCatalogDigest { .. }
                )
            }) =>
        {
            40
        }
        AnalysisJob::FinalizeCloningOutputs { .. } => 39,
        AnalysisJob::SelectDigestFragment { .. }
        | AnalysisJob::ExportPlanOligos { .. }
        | AnalysisJob::ExportProtocolHtml { .. } => 38,
        AnalysisJob::BuildAssemblyIntent { .. } => 36,
        AnalysisJob::DesignOverlapPrimers { .. } => 35,
        AnalysisJob::InspectAssemblyJunctions { plan, .. } => plan_worker_version(plan).max(33),
        AnalysisJob::PreviewOverlapRecipe { recipe, .. }
            if recipe.fragments.iter().any(|row| {
                matches!(
                    row.preparation,
                    crate::assembly_recipe::Preparation::CatalogDigest { .. }
                )
            }) =>
        {
            37
        }
        AnalysisJob::PreviewLigationRecipe { recipe, .. }
            if recipe.fragments.iter().any(|row| {
                matches!(
                    row.source.preparation,
                    crate::assembly_recipe::Preparation::CatalogDigest { .. }
                )
            }) =>
        {
            37
        }
        AnalysisJob::PreviewOverlapRecipe { recipe, .. }
            if recipe.fragments.iter().any(|row| {
                matches!(
                    row.preparation,
                    crate::assembly_recipe::Preparation::DeclaredDuplex { .. }
                )
            }) =>
        {
            34
        }
        AnalysisJob::PreviewLigationRecipe { recipe, .. }
            if recipe.fragments.iter().any(|row| {
                matches!(
                    row.source.preparation,
                    crate::assembly_recipe::Preparation::DeclaredDuplex { .. }
                )
            }) =>
        {
            34
        }
        AnalysisJob::PreviewLigationRecipe { .. } => 32,
        AnalysisJob::PreviewOverlapRecipe { .. } => 31,
        AnalysisJob::PrepareVariants { .. }
        | AnalysisJob::VariantPage { .. }
        | AnalysisJob::CheckVariantReference { .. }
        | AnalysisJob::PreviewVariantAllele { .. }
        | AnalysisJob::ExportVariants { .. }
        | AnalysisJob::DisposeVariants { .. } => 29,
        AnalysisJob::ProjectAnnotatedTopo { .. } | AnalysisJob::CheckAnnotatedTopoTarget { .. } => {
            27
        }
        AnalysisJob::SimulateTopo { .. }
        | AnalysisJob::InspectTopoOrientations { .. }
        | AnalysisJob::ProjectTopoProduct { .. }
        | AnalysisJob::CheckTopoTarget { .. } => 26,
        AnalysisJob::PrepareReadAlignments { .. }
        | AnalysisJob::ReadAlignmentPage { .. }
        | AnalysisJob::CheckReadAlignment { .. }
        | AnalysisJob::ExportReadAlignments { .. }
        | AnalysisJob::DisposeReadAlignments { .. } => 24,
        AnalysisJob::ScanRestriction { .. } => 23,
        AnalysisJob::InspectRecombination { .. }
        | AnalysisJob::CheckRecombinationTarget { .. }
        | AnalysisJob::SelectRecombinationComponent { .. } => 22,
        AnalysisJob::InspectCatalogRestriction { .. }
        | AnalysisJob::MaterializeCatalogDigestFragment { .. }
        | AnalysisJob::InstallCatalogDigestFragment { .. } => 21,
        AnalysisJob::InspectRestriction { .. }
        | AnalysisJob::MaterializeDigestFragment { .. }
        | AnalysisJob::InstallDigestFragment { .. } => 20,
        AnalysisJob::BeginGenomeFasta { .. }
        | AnalysisJob::PushGenomeFasta { .. }
        | AnalysisJob::FinishGenomeFasta { .. }
        | AnalysisJob::BeginGenomeNative { .. }
        | AnalysisJob::PushGenomeNativeBlock { .. }
        | AnalysisJob::FinishGenomeNative { .. }
        | AnalysisJob::AbortGenomeImport { .. }
        | AnalysisJob::RegisterGenome { .. }
        | AnalysisJob::GenomeIntervalGuides { .. }
        | AnalysisJob::StartGenomeSearch { .. }
        | AnalysisJob::AdvanceGenomeSearch { .. }
        | AnalysisJob::PageGenomeSearch { .. }
        | AnalysisJob::CancelGenomeJob { .. }
        | AnalysisJob::ReadGenomeInterval { .. } => 19,
        AnalysisJob::Search { .. }
        | AnalysisJob::Restriction { .. }
        | AnalysisJob::Pcr { .. }
        | AnalysisJob::Verify { .. } => 1,
        AnalysisJob::DuplexDigest { .. } => 2,
        AnalysisJob::Translate { .. }
        | AnalysisJob::TranslateCds { .. }
        | AnalysisJob::FindOrfs { .. }
        | AnalysisJob::Align { .. } => 3,
        AnalysisJob::PcrDocuments { .. }
        | AnalysisJob::AnalyzePrimer { .. }
        | AnalysisJob::DesignPrimers { .. }
        | AnalysisJob::ExportProtocolPdf { .. } => 4,
        AnalysisJob::SaveWorkspace { .. }
        | AnalysisJob::LoadWorkspace { .. }
        | AnalysisJob::ValidateTraces { .. } => 6,
        AnalysisJob::DesignCrisprGuides { .. }
        | AnalysisJob::SearchCrisprOffTargets { .. }
        | AnalysisJob::DesignEditTemplate { .. }
        | AnalysisJob::CheckEditTemplate { .. } => 18,
        AnalysisJob::AssessMethylation { .. } => 16,
        AnalysisJob::AssembleContigs { .. } => 15,
        AnalysisJob::DotPlot { .. } => 14,
        AnalysisJob::DesignCodons { .. } | AnalysisJob::CheckCodonDesign { .. } => 13,
        AnalysisJob::RecordCloningRun { .. } | AnalysisJob::ReplayCloningRun { .. } => 12,
        AnalysisJob::TailPreview { .. } | AnalysisJob::CheckOverlapProfile { .. } => 11,
        AnalysisJob::ProteinProperties { .. } => 10,
        AnalysisJob::ImportAlignment { .. }
        | AnalysisJob::ExportAlignment { .. }
        | AnalysisJob::ExtractAlignmentRows { .. }
        | AnalysisJob::SaveAlignment { .. }
        | AnalysisJob::LoadAlignment { .. } => 10,
        AnalysisJob::ViewTracks { .. }
        | AnalysisJob::SearchDocument { .. }
        | AnalysisJob::PrepareTracks { .. }
        | AnalysisJob::CachedViewTracks { .. } => 9,
        AnalysisJob::RepairPreview { .. }
        | AnalysisJob::AnnealPreview { .. }
        | AnalysisJob::AnnealedDocument { .. } => 8,
        AnalysisJob::ScanFeatureLibrary { .. }
        | AnalysisJob::ApplyFeatureLibrary { .. }
        | AnalysisJob::ImportAnnotations { .. }
        | AnalysisJob::ApplyAnnotations { .. }
        | AnalysisJob::ExportAnnotations { .. } => 7,
        _ => 5,
    };
    check_variant_version(bytes, required)?;
    check_preparation_version(bytes, &r.job)?;
    if let AnalysisJob::BeginCombinatorialBatch { workspace, matrix } = &r.job {
        check_workspace_version(bytes, workspace)?;
        check_combinatorial_request_version(bytes, matrix)?;
    }
    match &r.job {
        AnalysisJob::SelectGoldenGateReaction { report, .. } => {
            check_verification_report_version(bytes, &report.verification)?;
            check_variant_version(bytes, plan_worker_version(&report.probe_plan))?;
        }
        AnalysisJob::SelectPreparedGoldenGateCandidate { report, .. } => {
            check_verification_report_version(bytes, &report.preview.verification)?;
            check_variant_version(bytes, plan_worker_version(&report.preview.compiled.plan))?;
            check_variant_version(bytes, plan_worker_version(&report.preview.compiled.plan))?;
        }
        AnalysisJob::ImportAnnotations { format, .. }
        | AnalysisJob::ExportAnnotations { format, .. } => {
            check_annotation_format_version(bytes, *format)?;
        }
        AnalysisJob::ApplyAnnotations { report, .. } => {
            check_annotation_format_version(bytes, report.format)?;
        }
        _ => {}
    }
    if let AnalysisJob::SaveWorkspace { workspace } = &r.job {
        check_workspace_version(bytes, workspace)?;
        if workspace
            .attachments
            .values()
            .any(|a| a.kind == crate::run_lineage::WORKSPACE_ATTACHMENT_KIND)
        {
            check_variant_version(bytes, 12)?;
        }
        if workspace
            .attachments
            .values()
            .any(|a| a.kind == "read-assembly-v1")
        {
            check_variant_version(bytes, 15)?;
        }
        for c in &workspace.candidates {
            check_variant_version(bytes, plan_worker_version(&c.plan))?;
        }
    }
    if let AnalysisJob::RecordCloningRun {
        workspace, plan, ..
    } = &r.job
    {
        check_workspace_version(bytes, workspace)?;
        check_variant_version(bytes, plan_worker_version(plan))?;
    }
    if let AnalysisJob::LoadProcessedRunEvidence { workspace, run_id } = &r.job {
        if run_id.len() != 64 {
            return Err(Error::Format("Invalid run identity".into()));
        }
        check_workspace_version(bytes, workspace)?;
    }
    if let AnalysisJob::ReplayCloningRun { workspace, .. } = &r.job {
        check_workspace_version(bytes, workspace)?;
    }
    if let AnalysisJob::PreviewReadReview { workspace, .. }
    | AnalysisJob::CommitReadReview { workspace, .. } = &r.job
    {
        check_workspace_version(bytes, workspace)?;
    }
    if let AnalysisJob::CaptureReadCheck { workspace, .. }
    | AnalysisJob::ReplayReadCheck { workspace, .. } = &r.job
    {
        check_workspace_version(bytes, workspace)?;
    }
    if let AnalysisJob::FinalizeCloningOutputs {
        workspace, plan, ..
    } = &r.job
    {
        check_workspace_version(bytes, workspace)?;
        check_variant_version(bytes, plan_worker_version(plan))?;
    }
    if let AnalysisJob::Verify { plan, .. }
    | AnalysisJob::ExportProtocolPdf { plan, .. }
    | AnalysisJob::ExportProtocolHtml { plan, .. }
    | AnalysisJob::ExportPlanOligos { plan, .. } = &r.job
    {
        check_variant_version(bytes, plan_worker_version(plan))?;
        if plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::Tail3Prime { .. }
                    | crate::verification::Operation::ProfiledOverlap { .. }
            )
        }) {
            check_variant_version(bytes, 11)?;
        }
        if plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::RepairEnds { .. }
                    | crate::verification::Operation::AnnealOligos { .. }
            )
        }) {
            check_variant_version(bytes, 8)?;
        }
        if plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::PcrWithMismatches { .. }
                    | crate::verification::Operation::GoldenGate { .. }
            )
        }) {
            check_variant_version(bytes, 5)?;
        }
        if plan.steps.iter().any(|step| {
            matches!(
                step,
                crate::verification::Operation::Digest { .. }
                    | crate::verification::Operation::TreatEnds { .. }
                    | crate::verification::Operation::Ligate { .. }
            )
        }) {
            check_variant_version(bytes, 2)?;
        }
    }
    if let AnalysisJob::LoadWorkspace { bytes: archive } = &r.job {
        // A valid archive cannot hide new operation semantics in an old envelope.
        // Malformed archives retain the ordinary typed worker-error response below.
        if u16::from_le_bytes([bytes[4], bytes[5]]) < WORKER_SCHEMA_VERSION {
            if let Ok(workspace) = crate::workspace::load_workspace(archive) {
                check_workspace_version(bytes, &workspace)?;
            }
        }
    }
    if let AnalysisJob::PredictGel { lanes, .. } = &r.job {
        if lanes.iter().any(|lane| {
            matches!(
                lane.preparation,
                crate::gel::Preparation::CatalogDigest { .. }
            )
        }) {
            return Err(Error::Format(
                "Catalog-pinned lanes require the pinned-gel worker operation.".into(),
            ));
        }
    }
    let mol = match &r.job {
        AnalysisJob::FinalizeCloningOutputs {
            workspace,
            plan,
            target,
            inputs,
            options,
        } => {
            if crate::cloning_outputs::finalize_revision(
                workspace,
                plan,
                target.clone(),
                inputs.clone(),
                options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Cloning output source/options context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SelectDigestFragment {
            selection,
            source_role,
            molecule,
            options,
            expected_batch_revision,
            choice,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            if crate::digest_cut_selection::selection_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                options,
                expected_batch_revision,
                choice,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Digest cut selection source/boundaries changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ExportPlanOligos {
            plan,
            target,
            options,
            delimiter,
        } => {
            if crate::plan_oligos::export_revision(plan, target, options, *delimiter)? != r.revision
            {
                return Err(Error::Format("Plan oligo export context changed".into()));
            }
            None
        }
        AnalysisJob::InspectMismatchPcrDigest {
            molecule,
            origin,
            request,
        } => {
            if crate::mismatch_pcr_digest_inspection::inspection_revision(
                molecule,
                origin.as_ref(),
                request,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Mismatch PCR digest original source/origin/policy/catalog context changed"
                        .into(),
                ));
            }
            None
        }
        AnalysisJob::BuildContributionEditIntent { inputs, request } => {
            if crate::assembly_edit_intent::intent_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Contribution edit original source/declaration changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BuildSubstitutedAssemblyIntent { inputs, request } => {
            if crate::assembly_substitution_intent::intent_revision(inputs, request)? != r.revision
            {
                return Err(Error::Format(
                    "Substituted assembly original source/declaration changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BuildAssemblyIntent { inputs, request } => {
            if crate::assembly_intent::intent_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "assembly intent source/region/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BuildEditedAssemblyIntent {
            inputs,
            base,
            edits,
        } => {
            if crate::assembly_junction_edit::intent_revision(inputs, base, edits)? != r.revision {
                return Err(Error::Format(
                    "Edited assembly source/intent changed".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignEditedAssemblyPrimers { inputs, request } => {
            if crate::assembly_junction_edit::design_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Edited assembly source/primers changed".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignSelectiveAssemblyPrimers { inputs, request } => {
            if crate::selective_junction_edit::design_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Selective assembly source/primers changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectGoldenGateReaction {
            inputs,
            target,
            request,
        } => {
            if crate::golden_gate_reaction::inspect_revision(inputs, target, request)? != r.revision
            {
                return Err(Error::Format("Golden Gate reaction context changed".into()));
            }
            None
        }
        AnalysisJob::PreparePcrOutputs { request } => {
            if crate::pcr_outputs::revision(request)? != r.revision {
                return Err(Error::Format("PCR output selection changed".into()));
            }
            None
        }
        AnalysisJob::CopyReferenceFragment { request } => {
            if crate::reference_transfer::copy_revision(request)? != r.revision {
                return Err(Error::Format(
                    "Reference copy source or selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PasteReferenceFragment { request } => {
            if crate::reference_transfer::paste_revision(request)? != r.revision {
                return Err(Error::Format(
                    "Reference paste sources or selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreparePairwiseContext {
            report,
            reference,
            query,
            translate,
        } => {
            if crate::pairwise_context::revision(report, reference, query, *translate)?
                != r.revision
            {
                return Err(Error::Format("Pairwise annotation context changed.".into()));
            }
            None
        }
        AnalysisJob::SearchPairwiseComparison {
            report,
            pattern,
            options,
        } => {
            if crate::pairwise_search::revision(report, pattern, options)? != r.revision {
                return Err(Error::Format("Pairwise search context changed.".into()));
            }
            None
        }
        AnalysisJob::ExportPairwiseComparison {
            record,
            reference_name,
            query_name,
            fresh,
            format,
        } => {
            if crate::pairwise_export::revision(
                record,
                reference_name,
                query_name,
                *fresh,
                *format,
            )? != r.revision
            {
                return Err(Error::Format("Pairwise export context changed.".into()));
            }
            None
        }
        AnalysisJob::PredictPinnedGel {
            molecule,
            lanes,
            calibration,
        } => {
            if lanes.is_empty() || lanes.len() > 32 {
                return Err(Error::TooLarge);
            }
            calibration.validate()?;
            molecule.validate()?;
            let catalog = crate::catalog_digest::builtin_catalog()?;
            let mut bases = 0usize;
            for lane in lanes {
                lane.molecule.validate()?;
                bases = bases
                    .checked_add(lane.molecule.sequence.len())
                    .filter(|n| *n <= MAX_IMPORT_BYTES)
                    .ok_or(Error::TooLarge)?;
                if let crate::gel::Preparation::CatalogDigest { recipe } = &lane.preparation {
                    crate::model::bounds::model_preflight(&[&lane.molecule])?;
                    recipe.validate()?;
                    recipe.selection.validate_for(&catalog)?;
                    if lane.molecule.revision() != recipe.source.revision {
                        return Err(Error::Format("Pinned gel source revision changed.".into()));
                    }
                }
            }
            Some(molecule)
        }
        AnalysisJob::TranslateReference { request } => {
            if crate::translation_reference::inspect_revision(request)? != r.revision {
                return Err(Error::Format(
                    "Translation source or selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PrepareTranslatedReference { request, name } => {
            if crate::translation_reference::output_revision(request, name)? != r.revision {
                return Err(Error::Format(
                    "Translation source, selection or output name changed".into(),
                ));
            }
            None
        }
        AnalysisJob::CompareSequences {
            reference,
            query,
            options,
        } => {
            if crate::pairwise_review::revision(reference, query, options)? != r.revision {
                return Err(Error::Format(
                    "Comparison sources or options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PrepareMismatchPcrOutputs { request } => {
            if crate::mismatch_pcr_outputs::revision(request)? != r.revision {
                return Err(Error::Format(
                    "Mismatch-PCR output selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SelectGoldenGateReaction { report, selection } => {
            if crate::golden_gate_reaction::selection_revision(report, selection)? != r.revision {
                return Err(Error::Format(
                    "Golden Gate reaction selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ProcessTopoActions { inputs, request } => {
            if crate::topo_actions::revision(inputs, request)? != r.revision {
                return Err(Error::Format("TOPO actions product context changed".into()));
            }
            None
        }
        AnalysisJob::InspectTopoActions { inputs, request } => {
            if crate::topo_actions::revision_actions(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "TOPO action inspection context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ProjectProcessedTopo { inputs, request } => {
            if crate::topo_processed_product::revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Processed TOPO source/treatment/contour context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SealTopoNicks { inputs, request } => {
            if crate::topo_nick_sealing::revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "TOPO nick-sealing source/request context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ListLigationProfiles => {
            if crate::ligation_profiles::catalog_revision()? != r.revision {
                return Err(Error::Format(
                    "Published ligation catalog context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ScorePublishedLigation {
            profile,
            ends,
            physical_competitors,
        } => {
            if crate::ligation_profiles::score_revision_with_competitors(
                profile,
                ends,
                physical_competitors,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Published ligation profile/end context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ScoreLigationCompetition {
            dataset,
            ends,
            physical_competitors,
        } => {
            if crate::ligation_fidelity::competition_revision(dataset, ends, physical_competitors)?
                != r.revision
            {
                return Err(Error::Format("Ligation competition context changed".into()));
            }
            None
        }
        AnalysisJob::RankGoldenGateFusions { inputs, request } => {
            if crate::golden_gate_fusion_rank::inspect_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Fusion ranking source/profile context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectDomestication { source, request } => {
            if crate::domestication::inspect_revision(source, request)? != r.revision {
                return Err(Error::Format(
                    "Domestication inspection context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ProposeDomestication { source, request } => {
            if crate::domestication::propose_revision(source, request)? != r.revision {
                return Err(Error::Format(
                    "Domestication proposal context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::CheckDomestication {
            source,
            request,
            selection,
        } => {
            if crate::domestication::check_revision(source, request, selection)? != r.revision {
                return Err(Error::Format(
                    "Domestication substitution context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PrepareDomestication { input, request } => {
            if crate::domestication_preparation::revision(input, request)? != r.revision {
                return Err(Error::Format(
                    "Domestication preparation context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SearchGoldenGateFusions { inputs, request } => {
            if crate::golden_gate_fusion_design::inspect_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Golden Gate fusion search context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::MaterializeGoldenGateFusions {
            inputs,
            request,
            selection,
        } => {
            if crate::golden_gate_fusion_design::materialize_revision(inputs, request, selection)?
                != r.revision
            {
                return Err(Error::Format("Golden Gate fusion selection changed".into()));
            }
            None
        }
        AnalysisJob::ScoreLigationFidelity { dataset, ends } => {
            if crate::ligation_fidelity::revision(dataset, ends)? != r.revision {
                return Err(Error::Format(
                    "Ligation fidelity data or ends changed".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignGoldenGateJunctionPrimers { inputs, request } => {
            if crate::golden_gate_junction_edit::design_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "Golden Gate junction source/primers changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectPreparedGoldenGatePool {
            inputs,
            target,
            request,
        } => {
            if crate::prepared_golden_gate::inspect_revision(inputs, target, request)? != r.revision
            {
                return Err(Error::Format(
                    "Prepared Golden Gate pool context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SelectPreparedGoldenGateCandidate { report, selection } => {
            if crate::prepared_golden_gate::selection_revision(report, selection)? != r.revision {
                return Err(Error::Format(
                    "Prepared Golden Gate selection context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreviewReadReview { workspace, request } => {
            if crate::read_review::preview_revision(workspace, request)? != r.revision {
                return Err(Error::Format("Read review preview context changed".into()));
            }
            None
        }
        AnalysisJob::CommitReadReview { workspace, request } => {
            if crate::read_review::commit_revision(workspace, request)? != r.revision {
                return Err(Error::Format("Read review commit context changed".into()));
            }
            None
        }
        AnalysisJob::CaptureReadCheck {
            workspace,
            reference,
            reads,
            options,
        } => {
            if crate::read_checks::capture_revision(workspace, reference, reads, options)?
                != r.revision
            {
                return Err(Error::Format("Read-check capture context changed".into()));
            }
            None
        }
        AnalysisJob::ReplayReadCheck {
            workspace,
            attachment_hash,
        } => {
            if crate::read_checks::replay_revision(workspace, attachment_hash)? != r.revision {
                return Err(Error::Format("Read-check replay context changed".into()));
            }
            None
        }
        AnalysisJob::BuildMutationIntent { molecule, request } => {
            if crate::mutation_intent::intent_revision(molecule, request)? != r.revision {
                return Err(Error::Format("Mutation intent source/edit changed".into()));
            }
            None
        }
        AnalysisJob::DesignGoldenGatePrimers {
            inputs,
            target,
            request,
        } => {
            if crate::golden_gate_primer_design::preview_revision(inputs, target, request)?
                != r.revision
            {
                return Err(Error::Format(
                    "Golden Gate primer sources/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignOverlapPrimers {
            inputs,
            target,
            request,
        } => {
            if crate::overlap_primer_design::design_revision(inputs, target, request)? != r.revision
            {
                return Err(Error::Format(
                    "overlap primer sources/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectAssemblyJunctions {
            plan,
            target,
            options,
        } => {
            if crate::assembly_junctions::inspection_revision(plan, target, options)? != r.revision
            {
                return Err(Error::Format(
                    "assembly junction plan/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectRunMaterial {
            workspace,
            event,
            selection,
        } => {
            check_workspace_version(bytes, workspace)?;
            if crate::run_material::inspection_revision(workspace, event, selection)? != r.revision
            {
                return Err(Error::Format(
                    "Run material workspace/event/origin selection changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectRunCatalogRestriction {
            source_role,
            molecule,
            origin,
            selection,
            options,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            if crate::run_catalog_digest::inspection_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                origin,
                options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Saved-run catalog inspection origin/source/pin/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectRunPcrDigest {
            molecule,
            origin,
            request,
        } => {
            if crate::pcr_digest_inspection::inspection_revision_from_run(
                molecule, origin, request,
            )? != r.revision
            {
                return Err(Error::Format(
                    "Saved-run PCR digest origin/source/primer/catalog context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectPcrDigest { molecule, request } => {
            if crate::pcr_digest_inspection::inspection_revision(molecule, request)? != r.revision {
                return Err(Error::Format(
                    "PCR digest original source, primers or catalog options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreviewTailingRecipe {
            inputs,
            target,
            recipe,
        } => {
            if crate::tailing_recipe::preview_revision(inputs, target, recipe)? != r.revision {
                return Err(Error::Format(
                    "tailing preparation source/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreviewLigationRecipe {
            inputs,
            target,
            recipe,
        } => {
            if crate::ligation_recipe::preview_revision(inputs, target, recipe)? != r.revision {
                return Err(Error::Format(
                    "ligation preparation source/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::PreviewOverlapRecipe {
            inputs,
            target,
            recipe,
        } => {
            if crate::assembly_recipe::preview_revision(inputs, target, recipe)? != r.revision {
                return Err(Error::Format(
                    "overlap preparation source/target/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BeginCombinatorialBatch { .. }
        | AnalysisJob::AdvanceCombinatorialBatch { .. }
        | AnalysisJob::CaptureCombinatorialBatch { .. } => {
            return Err(Error::Unsupported(
                "combinatorial batches require a persistent private worker session".into(),
            ));
        }
        AnalysisJob::BeginCombinatorialRelease { .. }
        | AnalysisJob::AdvanceCombinatorialRelease { .. }
        | AnalysisJob::InspectCombinatorialRelease { .. }
        | AnalysisJob::FinishCombinatorialRelease { .. } => {
            return Err(Error::Unsupported(
                "combinatorial release requires a persistent private worker session".into(),
            ));
        }
        AnalysisJob::ExportCombinatorialOutputs {
            workspace,
            manifest_hash,
            options,
        } => {
            if crate::combinatorial_export::revision(workspace, manifest_hash, options)?
                != r.revision
            {
                return Err(Error::Format("Combinatorial export inputs changed.".into()));
            }
            None
        }
        AnalysisJob::PrepareVariants { .. }
        | AnalysisJob::VariantPage { .. }
        | AnalysisJob::CheckVariantReference { .. }
        | AnalysisJob::PreviewVariantAllele { .. }
        | AnalysisJob::ExportVariants { .. }
        | AnalysisJob::DisposeVariants { .. } => {
            return Err(Error::Unsupported(
                "variant operations require a persistent worker session".into(),
            ));
        }
        AnalysisJob::ProjectAnnotatedTopo {
            inputs,
            request,
            selection,
        } => {
            if crate::topo_worker::annotated_projection_revision(inputs, request, selection)?
                != r.revision
            {
                return Err(Error::Format(
                    "annotated TOPO projection context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::CheckAnnotatedTopoTarget {
            inputs,
            request,
            selection,
            target,
            policy,
        } => {
            if crate::topo_worker::annotated_target_revision(
                inputs, request, selection, target, policy,
            )? != r.revision
            {
                return Err(Error::Format(
                    "annotated TOPO target/policy context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SimulateTopo { inputs, request } => {
            if crate::topo_worker::simulation_revision(inputs, request)? != r.revision {
                return Err(Error::Format("TOPO simulation context changed".into()));
            }
            None
        }
        AnalysisJob::InspectTopoOrientations {
            inputs,
            request,
            vector,
            insert,
            directional,
        } => {
            if crate::topo_worker::orientations_revision(
                inputs,
                request,
                *vector,
                *insert,
                *directional,
            )? != r.revision
            {
                return Err(Error::Format("TOPO orientation context changed".into()));
            }
            None
        }
        AnalysisJob::ProjectTopoProduct {
            inputs,
            request,
            selection,
        } => {
            if crate::topo_worker::projection_revision(inputs, request, selection)? != r.revision {
                return Err(Error::Format("TOPO projection context changed".into()));
            }
            None
        }
        AnalysisJob::CheckTopoTarget {
            inputs,
            request,
            selection,
            target,
            policy,
        } => {
            if crate::topo_worker::target_revision(inputs, request, selection, target, policy)?
                != r.revision
            {
                return Err(Error::Format("TOPO target/policy context changed".into()));
            }
            None
        }
        AnalysisJob::PrepareReadAlignments { .. }
        | AnalysisJob::ReadAlignmentPage { .. }
        | AnalysisJob::CheckReadAlignment { .. }
        | AnalysisJob::ExportReadAlignments { .. }
        | AnalysisJob::DisposeReadAlignments { .. } => {
            return Err(Error::Unsupported(
                "read alignment operations require a persistent worker session".into(),
            ));
        }
        AnalysisJob::InspectCombinatorialMultiSite {
            workspace,
            design,
            construct,
        } => {
            check_workspace_version(bytes, workspace)?;
            if crate::combinatorial_design::inspect_multisite_revision(
                design, construct, workspace,
            )? != r.revision
            {
                return Err(Error::Format(
                    "combinatorial multisite inspection source/declaration context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::SimulateMultiSiteRecombination { inputs, request } => {
            if crate::multisite_recombination::revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "multisite sources/profile/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ProjectMultiSiteRecombination {
            inputs,
            request,
            selection,
        } => {
            if crate::multisite_recombination::product_revision(inputs, request, selection)?
                != r.revision
            {
                return Err(Error::Format("multisite selected component changed".into()));
            }
            None
        }
        AnalysisJob::CheckMultiSiteRecombination {
            inputs,
            request,
            selection,
            target,
        } => {
            if crate::multisite_recombination::check_revision(inputs, request, selection, target)?
                != r.revision
            {
                return Err(Error::Format(
                    "multisite independent target context changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectRecombination { inputs, request } => {
            if crate::recombination::batch_revision(inputs, request)? != r.revision {
                return Err(Error::Format(
                    "recombination source/profile/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::CheckRecombinationTarget {
            inputs,
            request,
            target,
            policy,
        } => {
            if crate::recombination::check_revision(inputs, request, target, policy)? != r.revision
            {
                return Err(Error::Format("recombination target binding changed".into()));
            }
            None
        }
        AnalysisJob::SelectRecombinationComponent {
            inputs,
            request,
            selection,
            target,
            policy,
        } => {
            if crate::recombination::checked_product_revision(
                inputs, request, selection, target, policy,
            )? != r.revision
            {
                return Err(Error::Format(
                    "recombination selected component binding changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectCatalogRestriction {
            selection,
            source_role,
            molecule,
            options,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            if crate::catalog_digest::batch_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "catalog restriction source/pin/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::MaterializeCatalogDigestFragment {
            selection,
            source_role,
            molecule,
            options,
            fragment_index,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            if crate::catalog_digest::product_revision(
                &catalog,
                selection,
                source_role,
                molecule,
                options,
                *fragment_index,
            )? != r.revision
            {
                return Err(Error::Format(
                    "catalog product source/pin/options/index changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InstallCatalogDigestFragment {
            selection,
            workspace,
            source,
            options,
            fragment_index,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            if crate::catalog_digest::installation_revision(
                &catalog,
                selection,
                workspace,
                source,
                options,
                *fragment_index,
            )? != r.revision
            {
                return Err(Error::Format(
                    "catalog installation workspace/source/pin/options/index changed".into(),
                ));
            }
            None
        }

        AnalysisJob::InstallDigestFragment {
            workspace,
            source,
            enzymes,
            options,
            fragment_index,
        } => {
            if crate::restriction_history::installation_revision(
                workspace,
                source,
                enzymes,
                options,
                *fragment_index,
            )? != r.revision
            {
                return Err(Error::Format(
                    "digest installation workspace/source/geometry/options/index changed".into(),
                ));
            }
            None
        }
        AnalysisJob::InspectRestriction {
            source_role,
            molecule,
            enzymes,
            options,
        }
        | AnalysisJob::ScanRestriction {
            source_role,
            molecule,
            enzymes,
            options,
        } => {
            if crate::restriction_inspection::batch_revision(
                source_role,
                molecule,
                enzymes,
                options,
            )? != r.revision
            {
                return Err(Error::Format(
                    "restriction inspection source/role/geometry/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::MaterializeDigestFragment {
            source_role,
            molecule,
            enzymes,
            options,
            fragment_index,
        } => {
            if crate::restriction_inspection::product_revision(
                source_role,
                molecule,
                enzymes,
                options,
                *fragment_index,
            )? != r.revision
            {
                return Err(Error::Format(
                    "digest product source/role/geometry/options/index changed".into(),
                ));
            }
            None
        }
        AnalysisJob::BeginGenomeFasta { .. }
        | AnalysisJob::PushGenomeFasta { .. }
        | AnalysisJob::FinishGenomeFasta { .. }
        | AnalysisJob::BeginGenomeNative { .. }
        | AnalysisJob::PushGenomeNativeBlock { .. }
        | AnalysisJob::FinishGenomeNative { .. }
        | AnalysisJob::AbortGenomeImport { .. }
        | AnalysisJob::RegisterGenome { .. }
        | AnalysisJob::GenomeIntervalGuides { .. }
        | AnalysisJob::StartGenomeSearch { .. }
        | AnalysisJob::AdvanceGenomeSearch { .. }
        | AnalysisJob::PageGenomeSearch { .. }
        | AnalysisJob::CancelGenomeJob { .. }
        | AnalysisJob::ReadGenomeInterval { .. } => None,
        AnalysisJob::DesignEditTemplate {
            source,
            target,
            options,
            guide,
        } => {
            if crate::edit_template::batch_revision(source, target, options, guide.as_ref())?
                != r.revision
            {
                return Err(Error::Format(
                    "edit-template source/target/options/guide changed".into(),
                ));
            }
            None
        }
        AnalysisJob::CheckEditTemplate {
            source,
            target,
            options,
            guide,
            donor,
        } => {
            if crate::edit_template::check_batch_revision(
                source,
                target,
                options,
                guide.as_ref(),
                donor,
            )? != r.revision
            {
                return Err(Error::Format(
                    "edit-template candidate/source/target/options/guide changed".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignCrisprGuides { molecule, options } => {
            if crate::crispr::guide_revision(molecule, options)? != r.revision {
                return Err(Error::Format("CRISPR guide source/options changed".into()));
            }
            None
        }
        AnalysisJob::SearchCrisprOffTargets {
            query,
            references,
            options,
        } => {
            if crate::crispr::off_target_revision(query, references, options)? != r.revision {
                return Err(Error::Format(
                    "CRISPR off-target query/references/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::AssessMethylation { molecule, options } => {
            if crate::methylation::batch_revision(molecule, options)? != r.revision {
                return Err(Error::Format(
                    "methylation source/options/profile changed".into(),
                ));
            }
            None
        }
        AnalysisJob::AssembleContigs { reads, options } => {
            if crate::contig::batch_revision(reads, options)? != r.revision {
                return Err(Error::Format("contig batch source/options changed".into()));
            }
            None
        }
        AnalysisJob::RecordCloningRun { workspace, .. }
        | AnalysisJob::ReplayCloningRun { workspace, .. }
        | AnalysisJob::LoadProcessedRunEvidence { workspace, .. } => {
            if crate::run_lineage::workspace_fingerprint(workspace)? != r.revision {
                return Err(Error::Format("run request source workspace changed".into()));
            }
            None
        }
        AnalysisJob::CachedViewTracks { .. } => None,
        AnalysisJob::CheckOverlapProfile {
            molecules,
            reverse,
            overlaps,
            circular,
            options,
        } => {
            if overlap_profile_revision(molecules, reverse, overlaps, *circular, options)?
                != r.revision
            {
                return Err(Error::Format(
                    "overlap profile request source/options changed".into(),
                ));
            }
            None
        }
        AnalysisJob::ImportAlignment { bytes, .. } | AnalysisJob::LoadAlignment { bytes } => {
            if blake3::hash(bytes).to_hex().as_str() != r.revision {
                return Err(Error::Format("alignment input digest mismatch".into()));
            }
            None
        }
        AnalysisJob::ExportAlignment { document, .. }
        | AnalysisJob::ExtractAlignmentRows { document, .. }
        | AnalysisJob::SaveAlignment { document } => {
            if alignment_revision(document)? != r.revision {
                return Err(Error::Format("alignment document revision mismatch".into()));
            }
            None
        }
        AnalysisJob::DotPlot {
            molecule,
            query,
            options,
        } => {
            if crate::dotplot::source_revisions(molecule, query, options)?.0 != r.revision {
                return Err(Error::Format(
                    "dot plot request source revision mismatch".into(),
                ));
            }
            None
        }
        AnalysisJob::DesignCodons { molecule, options }
        | AnalysisJob::CheckCodonDesign {
            source: molecule,
            options,
            ..
        } => {
            // Bound source evidence and derive the required protein before hashing.
            // These jobs must not take the generic molecule-revision path first.
            if crate::codon_design::source_revision(molecule, options)? != r.revision {
                return Err(Error::Format(
                    "codon request source revision mismatch".into(),
                ));
            }
            None
        }
        AnalysisJob::TailPreview { molecule, .. }
        | AnalysisJob::ProteinProperties { molecule, .. }
        | AnalysisJob::Search { molecule, .. }
        | AnalysisJob::PrepareTracks { molecule }
        | AnalysisJob::ViewTracks { molecule, .. }
        | AnalysisJob::SearchDocument { molecule, .. }
        | AnalysisJob::Restriction { molecule, .. }
        | AnalysisJob::Pcr { molecule, .. }
        | AnalysisJob::PcrDocuments { molecule, .. }
        | AnalysisJob::AnalyzePrimer { molecule, .. }
        | AnalysisJob::DesignPrimers { molecule, .. }
        | AnalysisJob::GoldenGate { molecule, .. }
        | AnalysisJob::MismatchPcr { molecule, .. }
        | AnalysisJob::InversePcr { molecule, .. }
        | AnalysisJob::PredictGel { molecule, .. }
        | AnalysisJob::DiagnosticDigests { molecule, .. }
        | AnalysisJob::MismatchPcrPreview { molecule, .. }
        | AnalysisJob::DuplexDigest { molecule, .. }
        | AnalysisJob::Translate { molecule, .. }
        | AnalysisJob::TranslateCds { molecule, .. }
        | AnalysisJob::FindOrfs { molecule, .. }
        | AnalysisJob::Align { molecule, .. }
        | AnalysisJob::ValidateTraces { molecule, .. }
        | AnalysisJob::ScanFeatureLibrary { molecule, .. }
        | AnalysisJob::ApplyFeatureLibrary { molecule, .. }
        | AnalysisJob::ImportAnnotations { molecule, .. }
        | AnalysisJob::ApplyAnnotations { molecule, .. }
        | AnalysisJob::ExportAnnotations { molecule, .. }
        | AnalysisJob::RepairPreview { molecule, .. }
        | AnalysisJob::AnnealPreview { molecule, .. }
        | AnalysisJob::AnnealedDocument { molecule, .. } => Some(molecule),
        AnalysisJob::Verify { target, .. }
        | AnalysisJob::ExportProtocolPdf { target, .. }
        | AnalysisJob::ExportProtocolHtml { target, .. } => Some(target),
        AnalysisJob::SaveWorkspace { workspace } => {
            workspace.validate()?;
            let head = workspace.head(&workspace.active_document)?;
            if head.revision() != r.revision {
                return Err(Error::Format(
                    "workspace worker input revision mismatch".into(),
                ));
            }
            None
        }
        AnalysisJob::LoadWorkspace { bytes } => {
            if bytes.len() > crate::workspace::MAX_ARCHIVE_BYTES {
                return Err(Error::TooLarge);
            }
            if blake3::hash(bytes).to_hex().as_str() != r.revision {
                return Err(Error::Format(
                    "workspace archive worker digest mismatch".into(),
                ));
            }
            None
        }
    };
    if let Some(mol) = mol {
        mol.validate()?;
        if mol.revision() != r.revision {
            return Err(Error::Format("worker input revision mismatch".into()));
        }
    }
    let result = match r.job {
        AnalysisJob::ProcessTopoActions { inputs, request } => {
            crate::topo_actions::inspect(&inputs, &request).map(AnalysisResult::TopoActionsProduct)
        }
        AnalysisJob::InspectTopoActions { inputs, request } => {
            crate::topo_actions::inspect_actions(&inputs, &request)
                .map(AnalysisResult::TopoActionsInspected)
        }
        AnalysisJob::LoadProcessedRunEvidence { workspace, run_id } => {
            crate::processed_run_evidence::download(&workspace, &run_id)
                .map(AnalysisResult::HistoricalProcessedRunEvidence)
        }
        AnalysisJob::ProjectProcessedTopo { inputs, request } => {
            crate::topo_processed_product::inspect(&inputs, &request)
                .map(AnalysisResult::ProcessedTopoProduct)
        }
        AnalysisJob::SealTopoNicks { inputs, request } => {
            crate::topo_nick_sealing::inspect(&inputs, &request)
                .map(AnalysisResult::TopoNicksSealed)
        }
        AnalysisJob::ListLigationProfiles => {
            crate::ligation_profiles::list().map(AnalysisResult::LigationProfiles)
        }
        AnalysisJob::ScorePublishedLigation {
            profile,
            ends,
            physical_competitors,
        } => {
            crate::ligation_profiles::score_with_competitors(&profile, &ends, &physical_competitors)
                .map(AnalysisResult::PublishedLigationScore)
        }
        AnalysisJob::ScoreLigationCompetition {
            dataset,
            ends,
            physical_competitors,
        } => crate::ligation_fidelity::score_competition(&dataset, &ends, &physical_competitors)
            .map(AnalysisResult::LigationCompetition),
        AnalysisJob::RankGoldenGateFusions { inputs, request } => {
            crate::golden_gate_fusion_rank::inspect(&inputs, &request)
                .map(AnalysisResult::RankedGoldenGateFusions)
        }
        AnalysisJob::InspectDomestication { source, request } => {
            crate::domestication::inspect(&source, &request)
                .map(AnalysisResult::DomesticationInspected)
        }
        AnalysisJob::ProposeDomestication { source, request } => {
            crate::domestication::propose(&source, &request)
                .map(AnalysisResult::DomesticationProposed)
        }
        AnalysisJob::CheckDomestication {
            source,
            request,
            selection,
        } => crate::domestication::check(&source, &request, &selection).map(|checked| {
            AnalysisResult::DomesticationChecked {
                report: checked.report().clone(),
                target: checked.target().clone(),
            }
        }),
        AnalysisJob::PrepareDomestication { input, request } => {
            crate::domestication_preparation::materialize(&input, &request)
                .map(AnalysisResult::DomesticationPrepared)
        }
        AnalysisJob::SearchGoldenGateFusions { inputs, request } => {
            crate::golden_gate_fusion_design::inspect(&inputs, &request)
                .map(AnalysisResult::GoldenGateFusions)
        }
        AnalysisJob::MaterializeGoldenGateFusions {
            inputs,
            request,
            selection,
        } => crate::golden_gate_fusion_design::materialize(&inputs, &request, &selection)
            .map(AnalysisResult::GoldenGateFusionMaterialized),
        AnalysisJob::ScoreLigationFidelity { dataset, ends } => {
            crate::ligation_fidelity::score(&dataset, &ends).map(AnalysisResult::LigationFidelity)
        }
        AnalysisJob::FinalizeCloningOutputs {
            workspace,
            plan,
            target,
            inputs,
            options,
        } => crate::cloning_outputs::finalize(&workspace, &plan, target, inputs, &options)
            .map(AnalysisResult::CloningOutputsFinalized),
        AnalysisJob::SelectDigestFragment {
            selection,
            source_role,
            molecule,
            options,
            expected_batch_revision,
            choice,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::digest_cut_selection::preview(
                &catalog,
                &selection,
                &source_role,
                &molecule,
                &options,
                &expected_batch_revision,
                &choice,
            )
            .map(AnalysisResult::DigestFragmentSelected)
        }
        AnalysisJob::ExportPlanOligos {
            plan,
            target,
            options,
            delimiter,
        } => crate::plan_oligos::export(&plan, &target, &options, delimiter)
            .map(AnalysisResult::PlanOligosExported),
        AnalysisJob::ExportProtocolHtml { plan, target } => Ok(AnalysisResult::ProtocolHtml(
            crate::verification::export_protocol_html(&plan, &target),
        )),
        AnalysisJob::InspectMismatchPcrDigest {
            molecule,
            origin,
            request,
        } => crate::mismatch_pcr_digest_inspection::inspect(&molecule, origin.as_ref(), &request)
            .map(AnalysisResult::MismatchPcrDigestInspection),
        AnalysisJob::BuildContributionEditIntent { inputs, request } => {
            crate::assembly_edit_intent::build(&inputs, &request)
                .map(AnalysisResult::ContributionEditIntent)
        }
        AnalysisJob::BuildSubstitutedAssemblyIntent { inputs, request } => {
            crate::assembly_substitution_intent::build(&inputs, &request)
                .map(AnalysisResult::SubstitutedAssemblyIntent)
        }
        AnalysisJob::BuildAssemblyIntent { inputs, request } => {
            crate::assembly_intent::build(&inputs, &request).map(AnalysisResult::AssemblyIntent)
        }
        AnalysisJob::BuildEditedAssemblyIntent {
            inputs,
            base,
            edits,
        } => crate::assembly_junction_edit::build_intent(&inputs, &base, &edits)
            .map(AnalysisResult::EditedAssemblyIntent),
        AnalysisJob::DesignEditedAssemblyPrimers { inputs, request } => {
            crate::assembly_junction_edit::design(&inputs, &request)
                .map(AnalysisResult::EditedAssemblyPrimers)
        }
        AnalysisJob::DesignSelectiveAssemblyPrimers { inputs, request } => {
            crate::selective_junction_edit::design(&inputs, &request)
                .map(AnalysisResult::SelectiveAssemblyPrimers)
        }
        AnalysisJob::InspectGoldenGateReaction {
            inputs,
            target,
            request,
        } => crate::golden_gate_reaction::inspect(&inputs, &target, &request)
            .map(AnalysisResult::GoldenGateReaction),
        AnalysisJob::PreparePcrOutputs { request } => {
            crate::pcr_outputs::prepare(&request).map(AnalysisResult::PcrOutputs)
        }
        AnalysisJob::PrepareMismatchPcrOutputs { request } => {
            crate::mismatch_pcr_outputs::prepare(&request).map(AnalysisResult::MismatchPcrOutputs)
        }
        AnalysisJob::CompareSequences {
            reference,
            query,
            options,
        } => crate::pairwise_review::compare(&reference, &query, &options)
            .map(AnalysisResult::ReviewedComparison),
        AnalysisJob::ExportPairwiseComparison {
            record,
            reference_name,
            query_name,
            fresh,
            format,
        } => crate::pairwise_export::export(&record, &reference_name, &query_name, fresh, format)
            .map(AnalysisResult::ExportedPairwiseComparison),
        AnalysisJob::SearchPairwiseComparison {
            report,
            pattern,
            options,
        } => crate::pairwise_search::search(&report, &pattern, options)
            .map(AnalysisResult::SearchedPairwiseComparison),
        AnalysisJob::PreparePairwiseContext {
            report,
            reference,
            query,
            translate,
        } => crate::pairwise_context::prepare(&report, &reference, &query, translate)
            .map(AnalysisResult::PreparedPairwiseContext),
        AnalysisJob::CopyReferenceFragment { request } => {
            crate::reference_transfer::copy(&request).map(AnalysisResult::ReferenceFragment)
        }
        AnalysisJob::PasteReferenceFragment { request } => {
            crate::reference_transfer::paste(&request).map(AnalysisResult::ReferencePaste)
        }
        AnalysisJob::TranslateReference { request } => {
            crate::translation_reference::inspect(&request).map(AnalysisResult::TranslatedReference)
        }
        AnalysisJob::PrepareTranslatedReference { request, name } => {
            crate::translation_reference::prepare(&request, &name)
                .map(AnalysisResult::TranslatedOutputs)
        }
        AnalysisJob::SelectGoldenGateReaction { report, selection } => {
            crate::golden_gate_reaction::select(&report, &selection)
                .map(AnalysisResult::GoldenGateReactionCandidate)
        }
        AnalysisJob::DesignGoldenGateJunctionPrimers { inputs, request } => {
            crate::golden_gate_junction_edit::design(&inputs, &request)
                .map(AnalysisResult::GoldenGateJunctionPrimers)
        }
        AnalysisJob::InspectPreparedGoldenGatePool {
            inputs,
            target,
            request,
        } => crate::prepared_golden_gate::inspect(&inputs, &target, &request)
            .map(AnalysisResult::PreparedGoldenGatePool),
        AnalysisJob::SelectPreparedGoldenGateCandidate { report, selection } => {
            crate::prepared_golden_gate::select(&report, &selection)
                .map(AnalysisResult::PreparedGoldenGateCandidate)
        }
        AnalysisJob::PreviewReadReview { workspace, request } => {
            crate::read_review::preview(&workspace, &request).map(AnalysisResult::ReadReviewPreview)
        }
        AnalysisJob::CommitReadReview { workspace, request } => {
            crate::read_review::commit(&workspace, &request)
                .map(AnalysisResult::ReadReviewCommitted)
        }
        AnalysisJob::CaptureReadCheck {
            workspace,
            reference,
            reads,
            options,
        } => crate::read_checks::capture(&workspace, reference, reads, &options)
            .map(AnalysisResult::ReadCheckCaptured),
        AnalysisJob::ReplayReadCheck {
            workspace,
            attachment_hash,
        } => crate::read_checks::replay(&workspace, &attachment_hash)
            .map(AnalysisResult::ReadCheckReplayed),
        AnalysisJob::BuildMutationIntent { molecule, request } => {
            crate::mutation_intent::build(&molecule, &request).map(AnalysisResult::MutationIntent)
        }
        AnalysisJob::DesignGoldenGatePrimers {
            inputs,
            target,
            request,
        } => crate::golden_gate_primer_design::design(&inputs, &target, &request)
            .map(AnalysisResult::GoldenGatePrimersDesigned),
        AnalysisJob::DesignOverlapPrimers {
            inputs,
            target,
            request,
        } => crate::overlap_primer_design::design(&inputs, &target, &request)
            .map(AnalysisResult::OverlapPrimersDesigned),
        AnalysisJob::InspectAssemblyJunctions {
            plan,
            target,
            options,
        } => crate::assembly_junctions::inspect_plan(&plan, &target, &options)
            .map(AnalysisResult::AssemblyJunctionInspection),
        AnalysisJob::InspectRunMaterial {
            workspace,
            event,
            selection,
        } => crate::run_material::inspect(&workspace, &event, &selection)
            .map(AnalysisResult::RunMaterialInspection),
        AnalysisJob::InspectRunCatalogRestriction {
            source_role,
            molecule,
            origin,
            selection,
            options,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::run_catalog_digest::inspect(
                &catalog,
                &selection,
                &source_role,
                &molecule,
                &origin,
                &options,
            )
            .map(AnalysisResult::RunCatalogRestrictionInspection)
        }
        AnalysisJob::InspectRunPcrDigest {
            molecule,
            origin,
            request,
        } => crate::pcr_digest_inspection::inspect_from_run(&molecule, &origin, &request)
            .map(AnalysisResult::RunPcrDigestInspection),
        AnalysisJob::InspectPcrDigest { molecule, request } => {
            crate::pcr_digest_inspection::inspect(&molecule, &request)
                .map(AnalysisResult::PcrDigestInspection)
        }
        AnalysisJob::PreviewTailingRecipe {
            inputs,
            target,
            recipe,
        } => crate::tailing_recipe::preview(&inputs, &target, &recipe)
            .map(AnalysisResult::TailingRecipePreview),
        AnalysisJob::PreviewLigationRecipe {
            inputs,
            target,
            recipe,
        } => crate::ligation_recipe::preview(&inputs, &target, &recipe)
            .map(AnalysisResult::LigationRecipePreview),
        AnalysisJob::PreviewOverlapRecipe {
            inputs,
            target,
            recipe,
        } => crate::assembly_recipe::preview(&inputs, &target, &recipe)
            .map(AnalysisResult::OverlapRecipePreview),
        AnalysisJob::PrepareVariants { .. }
        | AnalysisJob::VariantPage { .. }
        | AnalysisJob::CheckVariantReference { .. }
        | AnalysisJob::PreviewVariantAllele { .. }
        | AnalysisJob::ExportVariants { .. }
        | AnalysisJob::DisposeVariants { .. } => Err(Error::Unsupported(
            "variant operations require a persistent worker session".into(),
        )),
        AnalysisJob::ProjectAnnotatedTopo {
            inputs,
            request,
            selection,
        } => crate::topo::annotated_product(&inputs, &request, &selection)
            .map(AnalysisResult::AnnotatedTopoProjection),
        AnalysisJob::CheckAnnotatedTopoTarget {
            inputs,
            request,
            selection,
            target,
            policy,
        } => crate::topo::check_target_annotated(&inputs, &request, &selection, &target, &policy)
            .map(AnalysisResult::AnnotatedTopoTarget),
        AnalysisJob::SimulateTopo { inputs, request } => {
            crate::topo::simulate(&inputs, &request).map(AnalysisResult::TopoSimulation)
        }
        AnalysisJob::InspectTopoOrientations {
            inputs,
            request,
            vector,
            insert,
            directional,
        } => crate::topo::inspect_orientations(&inputs, &request, vector, insert, directional)
            .map(AnalysisResult::TopoOrientations),
        AnalysisJob::ProjectTopoProduct {
            inputs,
            request,
            selection,
        } => {
            crate::topo::product(&inputs, &request, &selection).map(AnalysisResult::TopoProjection)
        }
        AnalysisJob::CheckTopoTarget {
            inputs,
            request,
            selection,
            target,
            policy,
        } => crate::topo::check_target(&inputs, &request, &selection, &target, &policy)
            .map(AnalysisResult::TopoTarget),
        AnalysisJob::PrepareReadAlignments { .. }
        | AnalysisJob::ReadAlignmentPage { .. }
        | AnalysisJob::CheckReadAlignment { .. }
        | AnalysisJob::ExportReadAlignments { .. }
        | AnalysisJob::DisposeReadAlignments { .. } => Err(Error::Unsupported(
            "read alignment operations require a persistent worker session".into(),
        )),
        AnalysisJob::InspectCombinatorialMultiSite {
            workspace,
            design,
            construct,
        } => crate::combinatorial_design::inspect_multisite(&design, &construct, &workspace)
            .map(AnalysisResult::CombinatorialMultiSiteInspection),
        AnalysisJob::SimulateMultiSiteRecombination { inputs, request } => {
            crate::multisite_recombination::simulate(&inputs, &request)
                .map(AnalysisResult::MultiSiteRecombination)
        }
        AnalysisJob::ProjectMultiSiteRecombination {
            inputs,
            request,
            selection,
        } => crate::multisite_recombination::product(&inputs, &request, &selection)
            .map(AnalysisResult::MultiSiteRecombinationProduct),
        AnalysisJob::CheckMultiSiteRecombination {
            inputs,
            request,
            selection,
            target,
        } => crate::multisite_recombination::check_target(&inputs, &request, &selection, &target)
            .map(AnalysisResult::MultiSiteRecombinationChecked),
        AnalysisJob::InspectRecombination { inputs, request } => {
            crate::recombination::simulate(&inputs, &request)
                .map(AnalysisResult::RecombinationInspection)
        }
        AnalysisJob::CheckRecombinationTarget {
            inputs,
            request,
            target,
            policy,
        } => crate::recombination::check_target(&inputs, &request, &target, &policy)
            .map(AnalysisResult::RecombinationTargetCheck),
        AnalysisJob::SelectRecombinationComponent {
            inputs,
            request,
            selection,
            target,
            policy,
        } => crate::recombination::checked_product(&inputs, &request, &selection, &target, &policy)
            .map(AnalysisResult::RecombinationProduct),
        AnalysisJob::InspectCatalogRestriction {
            selection,
            source_role,
            molecule,
            options,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::catalog_digest::inspect(&catalog, &selection, &source_role, &molecule, &options)
                .map(AnalysisResult::CatalogRestrictionInspection)
        }
        AnalysisJob::MaterializeCatalogDigestFragment {
            selection,
            source_role,
            molecule,
            options,
            fragment_index,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::catalog_digest::product(
                &catalog,
                &selection,
                &source_role,
                &molecule,
                &options,
                fragment_index,
            )
            .map(AnalysisResult::CatalogDigestProduct)
        }
        AnalysisJob::InstallCatalogDigestFragment {
            selection,
            workspace,
            source,
            options,
            fragment_index,
        } => {
            let catalog = crate::catalog_digest::builtin_catalog()?;
            crate::catalog_digest::install(
                &catalog,
                &selection,
                &workspace,
                &source,
                &options,
                fragment_index,
            )
            .map(AnalysisResult::CatalogDigestInstalled)
        }

        AnalysisJob::InstallDigestFragment {
            workspace,
            source,
            enzymes,
            options,
            fragment_index,
        } => crate::restriction_history::install(
            &workspace,
            &source,
            &enzymes,
            &options,
            fragment_index,
        )
        .map(AnalysisResult::DigestInstalled),
        AnalysisJob::ScanRestriction {
            source_role,
            molecule,
            enzymes,
            options,
        } => crate::restriction_inspection::scan(&source_role, &molecule, &enzymes, &options)
            .map(AnalysisResult::RestrictionInspection),
        AnalysisJob::InspectRestriction {
            source_role,
            molecule,
            enzymes,
            options,
        } => crate::restriction_inspection::inspect(&source_role, &molecule, &enzymes, &options)
            .map(AnalysisResult::RestrictionInspection),
        AnalysisJob::MaterializeDigestFragment {
            source_role,
            molecule,
            enzymes,
            options,
            fragment_index,
        } => crate::restriction_inspection::product(
            &source_role,
            &molecule,
            &enzymes,
            &options,
            fragment_index,
        )
        .map(AnalysisResult::DigestProduct),
        AnalysisJob::BeginGenomeFasta { .. }
        | AnalysisJob::PushGenomeFasta { .. }
        | AnalysisJob::FinishGenomeFasta { .. }
        | AnalysisJob::BeginGenomeNative { .. }
        | AnalysisJob::PushGenomeNativeBlock { .. }
        | AnalysisJob::FinishGenomeNative { .. }
        | AnalysisJob::AbortGenomeImport { .. }
        | AnalysisJob::RegisterGenome { .. }
        | AnalysisJob::GenomeIntervalGuides { .. }
        | AnalysisJob::StartGenomeSearch { .. }
        | AnalysisJob::AdvanceGenomeSearch { .. }
        | AnalysisJob::PageGenomeSearch { .. }
        | AnalysisJob::CancelGenomeJob { .. }
        | AnalysisJob::ReadGenomeInterval { .. } => Err(Error::Unsupported(
            "genome operations require a persistent private worker session".into(),
        )),
        AnalysisJob::DesignEditTemplate {
            source,
            target,
            options,
            guide,
        } => crate::edit_template::design(&source, &target, &options, guide.as_ref())
            .map(AnalysisResult::EditTemplate),
        AnalysisJob::CheckEditTemplate {
            source,
            target,
            options,
            guide,
            donor,
        } => crate::edit_template::check(&source, &target, &options, guide.as_ref(), &donor)
            .map(AnalysisResult::EditTemplateCheck),
        AnalysisJob::DesignCrisprGuides { molecule, options } => {
            crate::crispr::design_guides(&molecule, &options).map(AnalysisResult::CrisprGuides)
        }
        AnalysisJob::SearchCrisprOffTargets {
            query,
            references,
            options,
        } => crate::crispr::search_off_targets(&query, &references, &options)
            .map(AnalysisResult::CrisprOffTargets),
        AnalysisJob::AssessMethylation { molecule, options } => {
            crate::methylation::assess(&molecule, &options).map(AnalysisResult::Methylation)
        }
        AnalysisJob::AssembleContigs { reads, options } => {
            crate::contig::assemble(&reads, &options).map(AnalysisResult::Contigs)
        }
        AnalysisJob::DotPlot {
            molecule,
            query,
            options,
        } => crate::dotplot::analyze(&molecule, &query, &options).map(AnalysisResult::DotPlot),
        AnalysisJob::DesignCodons { molecule, options } => {
            crate::codon_design::design(&molecule, &options).map(AnalysisResult::CodonDesign)
        }
        AnalysisJob::CheckCodonDesign {
            source,
            candidate,
            options,
        } => crate::codon_design::check(&source, &candidate, &options)
            .map(AnalysisResult::CodonCheck),
        AnalysisJob::RecordCloningRun {
            workspace,
            plan,
            target_binding,
            input_bindings,
            allow_incomplete,
        } => crate::run_lineage::record_run(
            &workspace,
            &plan,
            target_binding,
            input_bindings,
            allow_incomplete,
        )
        .map(AnalysisResult::RunCommit),
        AnalysisJob::ReplayCloningRun { workspace, run_id } => {
            crate::run_lineage::replay_run(&workspace, &run_id).map(AnalysisResult::RunReplay)
        }
        AnalysisJob::TailPreview {
            molecule,
            recipe,
            base,
        } => crate::verification::simulate_tailing(&molecule, &recipe, base)
            .map(AnalysisResult::TailPreview),
        AnalysisJob::CheckOverlapProfile {
            molecules,
            reverse,
            overlaps,
            circular,
            options,
        } => crate::overlap_profiles::check_oriented_fragments(
            &molecules, &reverse, &overlaps, circular, &options,
        )
        .map(AnalysisResult::OverlapProfile),
        AnalysisJob::ProteinProperties { molecule, options } => {
            crate::protein_properties::analyze_protein(&molecule, &options)
                .map(AnalysisResult::ProteinProperties)
        }
        AnalysisJob::ImportAlignment {
            bytes,
            format,
            alphabet,
            name,
        } => crate::alignment_formats::import_alignment(&bytes, format, alphabet, &name)
            .map(AnalysisResult::AlignmentDocument),
        AnalysisJob::LoadAlignment { bytes } => {
            crate::alignment_formats::load_alignment(&bytes).map(AnalysisResult::AlignmentDocument)
        }
        AnalysisJob::ExportAlignment { document, format } => {
            crate::alignment_formats::export_alignment(&document, format)
                .map(AnalysisResult::AlignmentExport)
        }
        AnalysisJob::ExtractAlignmentRows { document, row_ids } => document
            .extract_rows(&row_ids)
            .map(AnalysisResult::AlignmentRows),
        AnalysisJob::SaveAlignment { document } => {
            crate::alignment_formats::save_alignment(&document).map(AnalysisResult::AlignmentSaved)
        }
        AnalysisJob::BeginCombinatorialBatch { .. }
        | AnalysisJob::AdvanceCombinatorialBatch { .. }
        | AnalysisJob::CaptureCombinatorialBatch { .. } => Err(Error::Unsupported(
            "combinatorial batches require a persistent private worker session".into(),
        )),
        AnalysisJob::BeginCombinatorialRelease { .. }
        | AnalysisJob::AdvanceCombinatorialRelease { .. }
        | AnalysisJob::InspectCombinatorialRelease { .. }
        | AnalysisJob::FinishCombinatorialRelease { .. } => Err(Error::Unsupported(
            "combinatorial release requires a persistent private worker session".into(),
        )),
        AnalysisJob::ExportCombinatorialOutputs {
            workspace,
            manifest_hash,
            options,
        } => crate::combinatorial_export::export(&workspace, &manifest_hash, &options)
            .map(AnalysisResult::CombinatorialOutputsArchive),
        AnalysisJob::PrepareTracks { .. } | AnalysisJob::CachedViewTracks { .. } => Err(
            Error::Unsupported("cached track jobs require an AnalysisSession".into()),
        ),
        AnalysisJob::ViewTracks {
            molecule,
            window,
            options,
        } => crate::view_tracks::view_tracks(&molecule, window, &options)
            .map(AnalysisResult::ViewTracks),
        AnalysisJob::SearchDocument { molecule, options } => {
            crate::search::search_document(&molecule, &options).map(AnalysisResult::DocumentSearch)
        }
        AnalysisJob::Search { molecule, query } => {
            crate::analysis::find_matches(&molecule, &query, true).map(AnalysisResult::Search)
        }
        AnalysisJob::Restriction { molecule, enzyme } => {
            crate::analysis::restriction_sites(&molecule, &enzyme).map(AnalysisResult::Restriction)
        }
        AnalysisJob::Pcr {
            molecule,
            forward,
            reverse,
        } => crate::analysis::simulate_pcr(&molecule, &forward, &reverse).map(AnalysisResult::Pcr),
        AnalysisJob::RepairPreview {
            molecule,
            enzymes,
            fragment_index,
            method,
        } => crate::verification::simulate_end_repair(&molecule, &enzymes, fragment_index, method)
            .map(AnalysisResult::RepairPreview),
        AnalysisJob::AnnealPreview {
            molecule,
            bottom,
            options,
        } => crate::verification::simulate_oligo_annealing(&molecule, &bottom, &options)
            .map(AnalysisResult::AnnealPreview),
        AnalysisJob::AnnealedDocument {
            molecule,
            bottom,
            preview,
            candidate_index,
        } => crate::verification::annealed_document(&molecule, &bottom, &preview, candidate_index)
            .map(AnalysisResult::PhysicalPreviewDocument),
        AnalysisJob::PcrDocuments {
            molecule,
            forward,
            reverse,
        } => crate::verification::simulate_pcr_documents(&molecule, &forward, &reverse)
            .map(AnalysisResult::PcrDocuments),
        AnalysisJob::AnalyzePrimer {
            molecule,
            primer,
            conditions,
            binding_options,
        } => crate::primers::analyze_primer(&primer, &conditions).and_then(|analysis| {
            crate::primers::find_binding_sites(&molecule, &primer, &binding_options).map(
                |bindings| AnalysisResult::PrimerEvidence(PrimerEvidence { analysis, bindings }),
            )
        }),
        AnalysisJob::DesignPrimers { molecule, options } => {
            crate::primers::design_primers(&molecule, &options).map(AnalysisResult::PrimerDesign)
        }
        AnalysisJob::ExportProtocolPdf { target, plan } => {
            crate::pdf::export_protocol_pdf(&plan, &target).map(AnalysisResult::ProtocolPdf)
        }
        AnalysisJob::GoldenGate {
            templates, options, ..
        } => crate::golden_gate::simulate(&templates, &options).map(AnalysisResult::GoldenGate),
        AnalysisJob::MismatchPcr {
            molecule,
            forward,
            reverse,
            options,
        } => crate::mutagenesis::simulate_pcr(&molecule, &forward, &reverse, &options)
            .map(AnalysisResult::MismatchPcr),
        AnalysisJob::InversePcr {
            molecule,
            target,
            options,
        } => crate::mutagenesis::design_inverse_pcr(&molecule, &target, &options)
            .map(AnalysisResult::InversePcr),
        AnalysisJob::PredictGel {
            lanes, calibration, ..
        } => crate::gel::predict_gel(&lanes, &calibration).map(AnalysisResult::Gel),
        AnalysisJob::PredictPinnedGel {
            lanes, calibration, ..
        } => crate::gel::predict_gel(&lanes, &calibration).map(AnalysisResult::PinnedGel),
        AnalysisJob::DiagnosticDigests {
            candidates,
            options,
            ..
        } => crate::gel::diagnostic_digests(&candidates, &options)
            .map(AnalysisResult::DiagnosticDigest),
        AnalysisJob::MismatchPcrPreview {
            molecule,
            forward,
            reverse,
            options,
        } => crate::verification::simulate_mutagenesis_preview(
            &molecule, &forward, &reverse, &options,
        )
        .map(AnalysisResult::MismatchPcrPreview),
        AnalysisJob::DuplexDigest { molecule, enzymes } => {
            crate::chemistry::Duplex::from_molecule(&molecule)
                .and_then(|d| crate::chemistry::digest(&d, &enzymes))
                .map(AnalysisResult::DuplexDigest)
        }
        AnalysisJob::Translate { molecule, options } => {
            crate::translation::translate(&molecule.sequence, &options)
                .map(AnalysisResult::Translation)
        }
        AnalysisJob::TranslateCds {
            molecule,
            feature_index,
        } => molecule
            .features
            .get(feature_index)
            .ok_or_else(|| Error::Coordinates("CDS annotation no longer exists".into()))
            .and_then(|f| crate::translation::translate_cds(&molecule, f))
            .map(AnalysisResult::CdsTranslation),
        AnalysisJob::FindOrfs { molecule, options } => {
            crate::translation::find_orfs(&molecule, &options).map(AnalysisResult::Orfs)
        }
        AnalysisJob::Align {
            molecule,
            query,
            options,
        } => query
            .validate()
            .and_then(|_| crate::alignment::align(&molecule.sequence, &query.sequence, &options))
            .map(AnalysisResult::Alignment),
        AnalysisJob::Verify { target, plan } => Ok(AnalysisResult::Verification(
            crate::verification::verify(&plan, &target),
        )),
        AnalysisJob::SaveWorkspace { workspace } => {
            crate::workspace::save_workspace(&workspace).map(AnalysisResult::WorkspaceSaved)
        }
        AnalysisJob::LoadWorkspace { bytes } => {
            crate::workspace::load_workspace(&bytes).map(AnalysisResult::WorkspaceLoaded)
        }
        AnalysisJob::ValidateTraces {
            molecule,
            reads,
            options,
        } => crate::trace_validation::validate_reads(&molecule, &reads, &options)
            .map(AnalysisResult::TraceValidation),
        AnalysisJob::ScanFeatureLibrary {
            molecule,
            library,
            options,
        } => crate::feature_library::scan(&molecule, &library, &options)
            .map(AnalysisResult::FeatureLibraryScan),
        AnalysisJob::ApplyFeatureLibrary {
            molecule,
            library,
            report,
            selected,
        } => crate::feature_library::apply_selected(&molecule, &library, &report, &selected)
            .map(AnalysisResult::AnnotatedDocument),
        AnalysisJob::ImportAnnotations {
            molecule,
            text,
            format,
            seqid,
        } => crate::annotation_formats::import_annotations(&molecule, &text, format, &seqid)
            .map(AnalysisResult::AnnotationImport),
        AnalysisJob::ApplyAnnotations {
            molecule,
            report,
            selected,
        } => crate::annotation_formats::apply_annotations(&molecule, &report, &selected)
            .map(AnalysisResult::AnnotatedDocument),
        AnalysisJob::ExportAnnotations {
            molecule,
            format,
            seqid,
        } => crate::annotation_formats::export_annotations(&molecule, format, &seqid).map(|r| {
            AnalysisResult::AnnotationExport(AnnotationExportReport {
                text: r.text,
                warnings: r.warnings,
            })
        }),
    }
    .unwrap_or_else(|e| AnalysisResult::Error(e.to_string()));
    let maximum = if matches!(
        result,
        AnalysisResult::WorkspaceSaved(_)
            | AnalysisResult::WorkspaceLoaded(_)
            | AnalysisResult::RunCommit(_)
            | AnalysisResult::RunReplay(_)
            | AnalysisResult::CloningOutputsFinalized(_)
            | AnalysisResult::ReadCheckCaptured(_)
            | AnalysisResult::ReadReviewCommitted(_)
            | AnalysisResult::DigestInstalled(_)
            | AnalysisResult::CatalogDigestInstalled(_)
            | AnalysisResult::CombinatorialOutputsArchive(_)
    ) {
        MAX_WORKSPACE_MESSAGE_BYTES
    } else {
        MAX_IMPORT_BYTES
    };
    encode_with_limit(
        b"CGR\0",
        &Response {
            id: r.id,
            revision: r.revision,
            result,
        },
        maximum,
    )
}

/// A persistent worker retains a validated immutable track source/index. Replacing
/// it is atomic; viewport jobs carry only options/window and the expected revision.
#[derive(Default)]
pub struct AnalysisSession {
    tracks: Option<crate::view_tracks::PreparedTracks>,
    genomes: crate::genome_session::GenomeSession,
    read_alignments: crate::read_alignment_session::Session,
    variants: crate::variants_session::Session,
    combinatorial: Option<crate::combinatorial_batch::Session>,
    release: Option<crate::combinatorial_release::ReleaseSession>,
}
impl AnalysisSession {
    pub fn process_request(&mut self, bytes: &[u8]) -> Result<Vec<u8>, Error> {
        let maximum = worker_limit(
            b"CGQ\0",
            bytes,
            &[
                (19, 6),
                (20, 6),
                (42, 12),
                (43, 12),
                (69, 20),
                (72, 21),
                (102, 39),
                (105, 42),
                (106, 42),
                (107, 43),
                (108, 43),
                (128, 57),
                (130, 57),
                (131, 58),
                (134, 58),
                (135, 58),
                (149, 65),
                (155, 69),
                (158, 74),
            ],
        )?;
        // Route by a borrowed prefix so ordinary jobs are not fully decoded twice.
        // Stable append-only track, genome and read-alignment variants retain
        // private worker state; all other jobs use ordinary stateless dispatch.
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
            .map_err(|e| Error::Format(e.to_string()))?;
        if !matches!(tag, 32 | 33 | 53..=66 | 77..=81 | 88..=93 | 128..=134) {
            return process_request(bytes);
        }
        preflight_combinatorial_batch_tag(bytes, &[128, 129, 130])?;
        preflight_combinatorial_release_tag(bytes, &[131, 132, 133, 134])?;
        if (53..=66).contains(&tag) {
            check_variant_version(bytes, 19)?;
        }
        if (88..=93).contains(&tag) {
            check_variant_version(bytes, 29)?;
        }
        if (77..=81).contains(&tag) {
            check_variant_version(bytes, 24)?;
        }
        let r: Request = decode_with_limit(b"CGQ\0", bytes, maximum)?;
        let result = match r.job {
            AnalysisJob::BeginCombinatorialRelease { workspace, request } => {
                check_workspace_version(bytes, &workspace)?;
                let pending =
                    crate::combinatorial_release::ReleaseSession::begin(workspace, request)?;
                let progress = pending.summary();
                if progress.context_revision != r.revision {
                    return Err(Error::Format(
                        "Combinatorial release source or options changed.".into(),
                    ));
                }
                self.release = Some(pending);
                AnalysisResult::CombinatorialReleaseStarted(progress)
            }
            AnalysisJob::AdvanceCombinatorialRelease { context, index } => {
                if context != r.revision {
                    return Err(Error::Format(
                        "Combinatorial release advance context changed.".into(),
                    ));
                }
                let session = self.release.as_mut().ok_or_else(|| {
                    Error::Unsupported("No private combinatorial release has been started.".into())
                })?;
                let progress = match session.advance(&context, index) {
                    Ok(progress) => progress,
                    Err(error) => {
                        // A failed reconstruction poisons the complete stage.
                        // A rejected out-of-order request has not executed a row.
                        if session.summary().failed.is_some() {
                            self.release = None;
                        }
                        return Err(error);
                    }
                };
                if progress.failed.is_some() {
                    self.release = None;
                }
                AnalysisResult::CombinatorialReleaseAdvanced(progress)
            }
            AnalysisJob::InspectCombinatorialRelease { context, index } => {
                if context != r.revision {
                    return Err(Error::Format(
                        "Combinatorial release inspection context changed.".into(),
                    ));
                }
                let session = self.release.as_ref().ok_or_else(|| {
                    Error::Unsupported("No private combinatorial release has been started.".into())
                })?;
                if session.summary().context_revision != context {
                    return Err(Error::Format(
                        "Combinatorial release inspection session changed.".into(),
                    ));
                }
                AnalysisResult::CombinatorialReleaseInspected(session.inspect(index)?)
            }
            AnalysisJob::FinishCombinatorialRelease {
                workspace,
                current_context,
            } => {
                check_workspace_version(bytes, &workspace)?;
                if current_context != r.revision
                    || self
                        .release
                        .as_ref()
                        .is_none_or(|session| session.summary().context_revision != current_context)
                {
                    return Err(Error::Format(
                        "Combinatorial release finish context changed.".into(),
                    ));
                }
                // Consuming the private stage makes every failed finish atomic:
                // no decoded prefix or subsequent retry can publish its products.
                let session = self.release.take().ok_or_else(|| {
                    Error::Unsupported("No private combinatorial release has been started.".into())
                })?;
                AnalysisResult::CombinatorialReleaseFinished(
                    session.finish_to(&workspace, &current_context)?,
                )
            }
            AnalysisJob::BeginCombinatorialBatch { workspace, matrix } => {
                check_workspace_version(bytes, &workspace)?;
                check_combinatorial_request_version(bytes, &matrix)?;
                let pending = crate::combinatorial_batch::Session::begin(workspace, matrix)?;
                let progress = pending.summary();
                if progress.context_revision != r.revision {
                    return Err(Error::Format(
                        "Combinatorial batch source or matrix changed.".into(),
                    ));
                }
                self.combinatorial = Some(pending);
                AnalysisResult::CombinatorialBatchStarted(progress)
            }
            AnalysisJob::AdvanceCombinatorialBatch { context, index } => {
                if context != r.revision {
                    return Err(Error::Format(
                        "Combinatorial batch advance context changed.".into(),
                    ));
                }
                let session = self.combinatorial.as_mut().ok_or_else(|| {
                    Error::Unsupported("No private combinatorial batch has been started.".into())
                })?;
                AnalysisResult::CombinatorialBatchAdvanced(session.advance(&context, index)?)
            }
            AnalysisJob::CaptureCombinatorialBatch {
                workspace,
                current_context,
            } => {
                check_workspace_version(bytes, &workspace)?;
                if crate::combinatorial_batch::capture_revision(&workspace, &current_context)?
                    != r.revision
                {
                    return Err(Error::Format(
                        "Combinatorial batch capture workspace changed.".into(),
                    ));
                }
                let session = self.combinatorial.as_ref().ok_or_else(|| {
                    Error::Unsupported("No private combinatorial batch has been started.".into())
                })?;
                AnalysisResult::CombinatorialBatchCaptured(
                    session.capture_to(&workspace, &current_context)?,
                )
            }
            job @ (AnalysisJob::PrepareVariants { .. }
            | AnalysisJob::VariantPage { .. }
            | AnalysisJob::CheckVariantReference { .. }
            | AnalysisJob::PreviewVariantAllele { .. }
            | AnalysisJob::ExportVariants { .. }
            | AnalysisJob::DisposeVariants { .. }) => self.variants.process(job, &r.revision)?,
            job @ (AnalysisJob::PrepareReadAlignments { .. }
            | AnalysisJob::ReadAlignmentPage { .. }
            | AnalysisJob::CheckReadAlignment { .. }
            | AnalysisJob::ExportReadAlignments { .. }
            | AnalysisJob::DisposeReadAlignments { .. }) => {
                self.read_alignments.process(job, &r.revision)?
            }
            job @ (AnalysisJob::BeginGenomeFasta { .. }
            | AnalysisJob::PushGenomeFasta { .. }
            | AnalysisJob::FinishGenomeFasta { .. }
            | AnalysisJob::BeginGenomeNative { .. }
            | AnalysisJob::PushGenomeNativeBlock { .. }
            | AnalysisJob::FinishGenomeNative { .. }
            | AnalysisJob::AbortGenomeImport { .. }
            | AnalysisJob::RegisterGenome { .. }
            | AnalysisJob::GenomeIntervalGuides { .. }
            | AnalysisJob::StartGenomeSearch { .. }
            | AnalysisJob::AdvanceGenomeSearch { .. }
            | AnalysisJob::PageGenomeSearch { .. }
            | AnalysisJob::CancelGenomeJob { .. }
            | AnalysisJob::ReadGenomeInterval { .. }) => self.genomes.process(job, &r.revision)?,
            AnalysisJob::PrepareTracks { molecule } => {
                check_variant_version(bytes, 9)?;
                let prepared =
                    crate::view_tracks::PreparedTracks::new(std::sync::Arc::new(molecule))?;
                if prepared.source_revision() != r.revision {
                    return Err(Error::Format("track preparation revision mismatch".into()));
                }
                self.tracks = Some(prepared);
                AnalysisResult::TracksPrepared(r.revision.clone())
            }
            AnalysisJob::CachedViewTracks { window, options } => {
                check_variant_version(bytes, 9)?;
                let prepared = self.tracks.as_ref().ok_or_else(|| {
                    Error::Unsupported("cached track source has not been prepared".into())
                })?;
                AnalysisResult::ViewTracks(prepared.report_for_revision(
                    &r.revision,
                    window,
                    &options,
                )?)
            }
            _ => return process_request(bytes),
        };
        let maximum = if matches!(
            result,
            AnalysisResult::CombinatorialBatchCaptured(_)
                | AnalysisResult::CombinatorialReleaseFinished(_)
        ) {
            MAX_WORKSPACE_MESSAGE_BYTES
        } else {
            MAX_IMPORT_BYTES
        };
        encode_with_limit(
            b"CGR\0",
            &Response {
                id: r.id,
                revision: r.revision,
                result,
            },
            maximum,
        )
    }
}

/// Exact immutable alignment content identity, including original source evidence.
pub fn alignment_revision(
    document: &crate::alignment_formats::AlignmentDocument,
) -> Result<String, Error> {
    Ok(
        blake3::hash(&crate::alignment_formats::save_alignment(document)?)
            .to_hex()
            .to_string(),
    )
}
