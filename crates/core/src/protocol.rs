//! Versioned bounded binary messages. Schema changes require a new version.
use crate::{
    analysis::{CutSite, Enzyme, Match},
    formats::MAX_IMPORT_BYTES,
    model::*,
};
use serde::{Deserialize, Serialize};

fn current_version(magic: &[u8; 4]) -> u16 {
    match magic {
        b"CGQ\0" | b"CGR\0" => 23,
        b"CGP\0" => 7,
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
    let maximum = if matches!(
        r.job,
        AnalysisJob::SaveWorkspace { .. }
            | AnalysisJob::LoadWorkspace { .. }
            | AnalysisJob::RecordCloningRun { .. }
            | AnalysisJob::ReplayCloningRun { .. }
            | AnalysisJob::InstallDigestFragment { .. }
            | AnalysisJob::InstallCatalogDigestFragment { .. }
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
        &[(20, 6), (21, 6), (40, 12), (41, 12), (62, 20), (65, 21)],
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
        _ => decode_with_limit(b"CGR\0", bytes, maximum)?,
    };
    let required = match &response.result {
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
        AnalysisResult::Verification(report) => {
            check_scientific_report_version(bytes, &report.engine_version)?;
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
            check_scientific_report_version(bytes, &report.current_report.engine_version)?;
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
            if u16::from_le_bytes([bytes[4], bytes[5]]) < 19 {
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
fn check_scientific_report_version(bytes: &[u8], engine: &str) -> Result<(), Error> {
    if engine.ends_with("/science-9") {
        check_variant_version(bytes, 22)?;
    }
    if engine.ends_with("/science-8") {
        check_variant_version(bytes, 17)?;
    }
    Ok(())
}
fn check_workspace_version(
    bytes: &[u8],
    workspace: &crate::workspace::Workspace,
) -> Result<(), Error> {
    if u16::from_le_bytes([bytes[4], bytes[5]]) >= 22 {
        return Ok(());
    }
    for candidate in &workspace.candidates {
        check_variant_version(bytes, plan_worker_version(&candidate.plan))?;
    }
    for attachment in workspace.attachments.values() {
        if attachment.kind == crate::run_lineage::WORKSPACE_ATTACHMENT_KIND {
            check_variant_version(bytes, 12)?;
            let record = crate::run_lineage::load_run_record(&attachment.bytes)?;
            check_variant_version(bytes, plan_worker_version(&record.plan))?;
            if record.plan_schema >= 4 {
                check_variant_version(bytes, 17)?;
            }
            if record.plan_schema >= 5 {
                check_variant_version(bytes, 22)?;
            }
            check_scientific_report_version(bytes, &record.assessment.engine_version)?;
        } else if attachment.kind == crate::genome::WORKSPACE_ATTACHMENT_KIND {
            check_variant_version(bytes, 19)?;
        } else if attachment.kind == "read-assembly-v1" {
            check_variant_version(bytes, 15)?;
        }
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
        &[(19, 6), (20, 6), (42, 12), (43, 12), (69, 20), (72, 21)],
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
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
        .map_err(|e| Error::Format(e.to_string()))?;
    if tag == 76 {
        check_variant_version(bytes, 23)?;
    }
    let r: Request = decode_with_limit(b"CGQ\0", bytes, maximum)?;
    let required = match &r.job {
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
    if let AnalysisJob::ReplayCloningRun { workspace, .. } = &r.job {
        check_workspace_version(bytes, workspace)?;
    }
    if let AnalysisJob::Verify { plan, .. } | AnalysisJob::ExportProtocolPdf { plan, .. } = &r.job {
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
        if u16::from_le_bytes([bytes[4], bytes[5]]) < 22 {
            if let Ok(workspace) = crate::workspace::load_workspace(archive) {
                check_workspace_version(bytes, &workspace)?;
            }
        }
    }
    let mol = match &r.job {
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
        | AnalysisJob::ReplayCloningRun { workspace, .. } => {
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
        AnalysisJob::Verify { target, .. } | AnalysisJob::ExportProtocolPdf { target, .. } => {
            Some(target)
        }
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
            | AnalysisResult::DigestInstalled(_)
            | AnalysisResult::CatalogDigestInstalled(_)
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
}
impl AnalysisSession {
    pub fn process_request(&mut self, bytes: &[u8]) -> Result<Vec<u8>, Error> {
        let maximum = worker_limit(
            b"CGQ\0",
            bytes,
            &[(19, 6), (20, 6), (42, 12), (43, 12), (69, 20), (72, 21)],
        )?;
        // Route by a borrowed prefix so ordinary jobs are not fully decoded twice.
        // Stable append-only tags 32/33 are the only session-specific variants.
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
            .map_err(|e| Error::Format(e.to_string()))?;
        if !matches!(tag, 32 | 33 | 53..=66) {
            return process_request(bytes);
        }
        if (53..=66).contains(&tag) {
            check_variant_version(bytes, 19)?;
        }
        let r: Request = decode_with_limit(b"CGQ\0", bytes, maximum)?;
        let result = match r.job {
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
        encode(
            b"CGR\0",
            &Response {
                id: r.id,
                revision: r.revision,
                result,
            },
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
