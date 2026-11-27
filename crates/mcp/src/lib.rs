//! Local, session-scoped MCP access to the same deterministic core as the browser.
//! Client text and feature notes are data; no tool accepts paths, shell commands or permissions.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{
    Molecule, alignment, alignment_formats,
    analysis::{self, Primer},
    annotation_formats, chemistry, feature_library, formats, oligo_formats, primers,
    trace_validation, transfer, translation,
    verification::{self, DesignPlan, Operation, PlanInput},
    workspace,
};
use rmcp::{
    ErrorData, RoleServer, ServerHandler,
    handler::server::{router::tool::ToolRouter, wrapper::Parameters},
    model::{
        CallToolResult, ListResourcesResult, PaginatedRequestParams, ReadResourceRequestParams,
        ReadResourceResponse, ReadResourceResult, Resource, ResourceContents, ServerCapabilities,
        ServerConfig,
    },
    service::RequestContext,
    tool, tool_handler, tool_router,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::{Arc, Mutex},
};

/// Number of tools in the current native registry (historical manifests remain literal).
pub const CURRENT_TOOL_COUNT: usize = 189;

const MAX_STORE_BYTES: usize = 64 * 1024 * 1024;
const MAX_ITEMS: usize = 128;
const MAX_DOCUMENT_REVISIONS: usize = 10_000;
const MAX_UPLOADS: usize = 4;
const MAX_UPLOAD_CHUNK: usize = 256 * 1024;
const MAX_READ_ALIGNMENT_SOURCE: usize = 4 * 1024 * 1024;

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DocumentRef {
    /// Immutable session document ID, returned by import_document.
    pub id: String,
    /// Exact expected revision, returned by import_document or list_documents.
    pub revision: String,
}
/// Explicit writer profiles; unrelated options cannot be silently ignored.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum AnnotatedSequenceFormat {
    GenBank,
    Embl,
    GenPept {
        /// Actual export date DD-MMM-YYYY; absent uses retained source date only.
        locus_date: Option<String>,
        entry_name: Option<String>,
    },
    UniProt {
        entry_name: Option<String>,
        /// Explicit valid UniProt accession; never allocated or invented.
        accession: Option<String>,
    },
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnnotatedSequenceExportQuery {
    pub document: DocumentRef,
    pub format: AnnotatedSequenceFormat,
    /// Default false returns warning review without an artifact when losses exist.
    #[serde(default)]
    pub accept_omissions: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoSourceRef {
    pub role: String,
    pub document: DocumentRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoSimulateQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Preparation,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoOrientationQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Preparation,
    /// Exact prepared-strand indexes; not source-document indexes or a kit guess.
    pub vector: [usize; 2],
    pub insert: [usize; 2],
    pub directional: bool,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoProjectQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Preparation,
    pub selection: collagene_core::topo::Selection,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TopoTargetQuery {
    pub sources: Vec<TopoSourceRef>,
    pub request: topo_graph_tools::Preparation,
    pub selection: collagene_core::topo::Selection,
    pub target: DocumentRef,
    /// Every allowance is explicit; the target never chooses a component/seed.
    pub policy: collagene_core::topo::AcceptancePolicy,
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
/// Saved template context is attested against an exact session-owned archive;
/// the task can never supply an alleged reconstructed duplex.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum SavedTopoTask {
    Orientations {
        request: topo_graph_tools::Preparation,
        vector: [usize; 2],
        insert: [usize; 2],
        directional: bool,
    },
    Project {
        request: topo_graph_tools::Preparation,
        selection: collagene_core::topo::Selection,
    },
    Check {
        request: topo_graph_tools::Preparation,
        selection: collagene_core::topo::Selection,
        target: DocumentRef,
        policy: collagene_core::topo::AcceptancePolicy,
    },
    Simulate {
        request: topo_graph_tools::Preparation,
    },
    ProcessActions {
        request: topo_graph_tools::ActionProduct,
    },
}
impl SavedTopoTask {
    fn preparation(&self) -> &topo_graph_tools::Preparation {
        match self {
            Self::Orientations { request, .. }
            | Self::Project { request, .. }
            | Self::Check { request, .. }
            | Self::Simulate { request } => request,
            Self::ProcessActions { request } => &request.preparation,
        }
    }
    fn native(
        &self,
        target: Option<Molecule>,
        preparation: collagene_core::topo::Request,
    ) -> Result<collagene_core::topo_run_preparation::Task, ErrorData> {
        use collagene_core::topo_run_preparation::Task;
        Ok(match self {
            Self::Orientations {
                request: _,
                vector,
                insert,
                directional,
            } => Task::Orientations {
                request: preparation,
                vector: *vector,
                insert: *insert,
                directional: *directional,
            },
            Self::Project {
                request: _,
                selection,
            } => Task::Project {
                request: preparation,
                selection: selection.clone(),
            },
            Self::Check {
                request: _,
                selection,
                policy,
                ..
            } => Task::Check {
                request: preparation,
                selection: selection.clone(),
                target: target.ok_or_else(|| bad("Independent target reference is required"))?,
                policy: policy.clone(),
            },
            Self::Simulate { request: _ } => Task::Simulate {
                request: preparation,
            },
            Self::ProcessActions { request } => Task::ProcessActions {
                request: collagene_core::topo_actions::Request {
                    preparation,
                    selection: request.selection.clone(),
                    options: request.options.clone(),
                },
            },
        })
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SavedTopoInspectionQuery {
    pub inputs: Vec<BoundInput>,
    #[serde(deserialize_with = "Option::deserialize")]
    pub template: Option<topo_graph_tools::BoundTemplate>,
    pub origin_workspace: WorkspaceRef,
    pub task: SavedTopoTask,
}
#[derive(Serialize)]
struct TopoEvidence<'a, T: Serialize> {
    origin_workspace: &'a Option<WorkspaceRef>,
    source_bindings: &'a [TopoSourceRef],
    request: &'a collagene_core::topo::Request,
    target_requirement: Option<&'a DocumentRef>,
    report: &'a T,
    binding_scope: &'static str,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportText {
    /// GenBank, FASTA or plain DNA contents. Data only, never a file path.
    pub text: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportUrl {
    /// One nucleotide record: bare/NCBI accession.version or NCBI URL, ena:ACCESSION.VERSION or an official ENA record URL, ddbj:ACCESSION.VERSION or an official DDBJ record URL, Addgene plasmid page/ID, or direct public HTTPS sequence-file URL. Bare accessions still select NCBI; numeric IDs still select Addgene. Collections, ranges and search URLs are not complete-record imports.
    pub url: String,
    /// Addgene-only catalog sequence ID selected using list_source_sequences; omit for NCBI, ENA, DDBJ and direct HTTPS.
    #[serde(default)]
    pub sequence_id: Option<u64>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SourceQuery {
    pub url: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SearchDocumentQuery {
    pub document: DocumentRef,
    pub options: collagene_core::search::SearchOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ViewTracksQuery {
    pub document: DocumentRef,
    pub window: collagene_core::view_tracks::TrackWindow,
    pub options: collagene_core::view_tracks::TrackOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AlignmentImport {
    /// Exact file bytes as base64, at most8MiB decoded. No paths, URLs or aligner invocation.
    pub data_base64: String,
    pub format: alignment_formats::AlignmentFormat,
    pub alphabet: collagene_core::Alphabet,
    pub name: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AlignmentDocumentImport {
    pub data_base64: String,
}

#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum VariantFormat {
    Vcf,
    Bcf,
}
impl VariantFormat {
    fn core(self) -> collagene_core::variant_reference::exchange::Format {
        match self {
            Self::Vcf => collagene_core::variant_reference::exchange::Format::Vcf,
            Self::Bcf => collagene_core::variant_reference::exchange::Format::Bcf,
        }
    }
}
#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum VariantExportFormat {
    Original,
    Vcf,
    Bcf,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VariantImport {
    pub data_base64: String,
    pub format: VariantFormat,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VariantQuery {
    pub collection: EvidenceRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VariantPage {
    pub collection: EvidenceRef,
    #[serde(default)]
    pub start: usize,
    #[serde(default = "read_preview")]
    pub limit: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VariantCheck {
    pub collection: EvidenceRef,
    pub record_index: usize,
    pub reference: DocumentRef,
    /// Explicit file contig to actual immutable reference DNA mapping.
    pub reference_name: String,
    /// One-based alternative. Never selected from sample genotype columns.
    pub alternate_index: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VariantReplay {
    /// Exact source hash of a private completed collection; no decoded receipt authority.
    pub source_id: String,
    pub reference: DocumentRef,
    pub design: DocumentRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VariantExport {
    pub collection: EvidenceRef,
    pub format: VariantExportFormat,
    #[serde(default = "read_export_cap")]
    pub max_bytes: usize,
}

#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum ReadAlignmentFormat {
    Sam,
    Bam,
}
impl ReadAlignmentFormat {
    fn core(self) -> collagene_core::read_alignment_store::Format {
        match self {
            Self::Sam => collagene_core::read_alignment_store::Format::Sam,
            Self::Bam => collagene_core::read_alignment_store::Format::Bam,
        }
    }
}
#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum ReadAlignmentExportFormat {
    Original,
    Sam,
    Bam,
}
fn read_source_cap() -> usize {
    MAX_READ_ALIGNMENT_SOURCE
}
fn read_retained_cap() -> usize {
    32 * 1024 * 1024
}
fn read_record_cap() -> usize {
    10_000
}
fn read_preview() -> usize {
    20
}
fn read_export_cap() -> usize {
    8 * 1024 * 1024
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReadAlignmentImport {
    /// Complete SAM or BAM/BGZF bytes, base64; at most4MiB decoded. No paths or URLs.
    pub data_base64: String,
    pub format: ReadAlignmentFormat,
    #[serde(default = "read_source_cap")]
    pub max_source_bytes: usize,
    /// Managed retention charge, capped at64MiB and remaining shared session budget.
    #[serde(default = "read_retained_cap")]
    pub max_retained_bytes: usize,
    #[serde(default = "read_record_cap")]
    pub max_records: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReadAlignmentPage {
    /// Source BLAKE3 returned only after complete validated file parsing.
    pub collection: EvidenceRef,
    #[serde(default)]
    pub start: usize,
    #[serde(default = "read_preview")]
    pub limit: usize,
    #[serde(default)]
    pub reference_start: usize,
    #[serde(default = "read_preview")]
    pub reference_limit: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReadAlignmentCheck {
    pub collection: EvidenceRef,
    pub record_index: usize,
    /// Explicit actual DNA reference; its revision supplies expected_reference_revision.
    pub reference: DocumentRef,
    /// Explicit SAM/BAM dictionary mapping, never inferred from document name/M5/UR.
    pub reference_name: String,
    #[serde(default = "default_read_check_columns")]
    pub max_columns: usize,
    #[serde(default = "default_read_check_work")]
    pub max_work: u64,
    #[serde(default = "default_read_check_report")]
    pub max_report_bytes: usize,
    #[serde(default = "read_preview")]
    pub preview_columns: usize,
}
fn default_read_check_columns() -> usize {
    100_000
}
fn default_read_check_work() -> u64 {
    128_000_000
}
fn default_read_check_report() -> usize {
    16 * 1024 * 1024
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReadAlignmentExport {
    pub collection: EvidenceRef,
    /// original is byte-exact; BAM conversion may reframe BGZF and returns notices.
    pub format: ReadAlignmentExportFormat,
    #[serde(default = "read_export_cap")]
    pub max_bytes: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AlignmentQuery {
    pub alignment: EvidenceRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AlignmentRegionQuery {
    pub alignment: EvidenceRef,
    pub row_start: usize,
    pub max_rows: usize,
    pub column_start: usize,
    pub max_columns: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
pub enum AlignmentExportFormat {
    AlignedFasta,
    Clustal,
    Stockholm,
    Native,
    Original,
    PhylipSequential,
    PhylipInterleaved,
    PhylipRelaxedSequential,
    PhylipRelaxedInterleaved,
    Msf,
    Nexus,
    Pir,
    Gde,
    Selex,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AlignmentExportQuery {
    pub alignment: EvidenceRef,
    pub format: AlignmentExportFormat,
    /// Additional exchange profiles only; omission retains lossless defaults.
    #[serde(default)]
    pub options: Option<collagene_core::alignment_exchange::ExportOptions>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AlignmentRowsQuery {
    pub alignment: EvidenceRef,
    pub row_ids: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ContigReadRef {
    /// Stable, unique role identifier for this ordered read input.
    pub id: String,
    pub document: DocumentRef,
}
fn genome_chunk_default() -> usize {
    65536
}
fn genome_finish_default() -> usize {
    8
}
fn genome_page_default() -> usize {
    25
}
fn genome_advance_default() -> usize {
    1
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub(crate) struct GenomeBeginQuery {
    pub format: genome_store::ImportFormat,
    pub expected_bytes: u64,
    pub expected_hash: Option<String>,
    pub build: Option<collagene_core::genome::BuildInfo>,
    pub options: Option<collagene_core::genome::IngestOptions>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeAppendQuery {
    pub import_id: String,
    pub offset: u64,
    pub data_base64: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeFinishQuery {
    pub import_id: String,
    #[serde(default = "genome_finish_default")]
    pub max_blocks: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeInspectQuery {
    pub genome: Option<collagene_core::genome_session::GenomeRef>,
    #[serde(default)]
    pub contig_start: usize,
    #[serde(default = "genome_page_default")]
    pub max_contigs: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeExportQuery {
    pub genome: collagene_core::genome_session::GenomeRef,
    #[serde(default)]
    pub offset: u64,
    #[serde(default = "genome_chunk_default")]
    pub max_bytes: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeGuidesQuery {
    pub genome: collagene_core::genome_session::GenomeRef,
    pub options: collagene_core::genome_crispr::GenomeGuideOptions,
    #[serde(default)]
    pub candidate_start: usize,
    #[serde(default = "genome_page_default")]
    pub max_candidates: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeReadQuery {
    pub genome: collagene_core::genome_session::GenomeRef,
    pub options: collagene_core::genome_crispr::GenomeIntervalReadOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeSearchQuery {
    pub genome: collagene_core::genome_session::GenomeRef,
    pub query: collagene_core::crispr::GuideQuery,
    #[serde(default)]
    pub options: collagene_core::genome_crispr::GenomeSearchOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeAdvanceQuery {
    pub job_id: String,
    #[serde(default = "genome_advance_default")]
    pub max_blocks: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomePageQuery {
    pub job_id: String,
    #[serde(default)]
    pub offset: usize,
    #[serde(default = "genome_page_default")]
    pub limit: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenomeCancelQuery {
    pub id: String,
    pub genome: Option<collagene_core::genome_session::GenomeRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditTemplateQuery {
    pub source: DocumentRef,
    pub target: DocumentRef,
    pub options: collagene_core::edit_template::EditTemplateOptions,
    #[serde(default)]
    pub guide: Option<collagene_core::edit_template::GuideBinding>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditTemplateCheckQuery {
    pub source: DocumentRef,
    pub target: DocumentRef,
    pub donor: DocumentRef,
    pub options: collagene_core::edit_template::EditTemplateOptions,
    #[serde(default)]
    pub guide: Option<collagene_core::edit_template::GuideBinding>,
}
fn default_crispr_page() -> usize {
    20
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CrisprGuidesQuery {
    pub document: DocumentRef,
    #[serde(default)]
    pub options: collagene_core::crispr::GuideOptions,
    #[serde(default)]
    pub candidate_start: usize,
    #[serde(default = "default_crispr_page")]
    pub max_candidates: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CrisprReferenceRef {
    pub id: String,
    pub document: DocumentRef,
    pub build: Option<String>,
    pub accession: Option<String>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CrisprOffTargetsQuery {
    pub query: collagene_core::crispr::GuideQuery,
    pub references: Vec<CrisprReferenceRef>,
    #[serde(default)]
    pub options: collagene_core::crispr::OffTargetOptions,
    #[serde(default)]
    pub hit_start: usize,
    #[serde(default = "default_crispr_page")]
    pub max_hits: usize,
}
fn default_methylation_page() -> usize {
    20
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MethylationQuery {
    pub document: DocumentRef,
    pub options: collagene_core::methylation::MethylationOptions,
    #[serde(default)]
    pub site_start: usize,
    #[serde(default = "default_methylation_page")]
    pub max_sites: usize,
}
fn default_contig_page() -> usize {
    5
}
fn default_overlap_page() -> usize {
    10
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AssembleContigsQuery {
    pub reads: Vec<ContigReadRef>,
    pub options: collagene_core::contig::ContigOptions,
    #[serde(default)]
    pub contig_start: usize,
    #[serde(default = "default_contig_page")]
    pub max_contigs: usize,
    #[serde(default)]
    pub overlap_start: usize,
    #[serde(default = "default_overlap_page")]
    pub max_overlaps: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProteinPropertiesQuery {
    pub document: DocumentRef,
    pub options: collagene_core::protein_properties::ProteinOptions,
}
fn default_dotplot_page() -> usize {
    50
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DotPlotQuery {
    pub reference: DocumentRef,
    pub query: DocumentRef,
    pub options: collagene_core::dotplot::DotPlotOptions,
    #[serde(default)]
    pub hit_start: usize,
    #[serde(default = "default_dotplot_page")]
    pub max_hits: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CodonDesignQuery {
    pub document: DocumentRef,
    pub options: collagene_core::codon_design::DesignOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CodonCheckQuery {
    pub source: DocumentRef,
    pub candidate: DocumentRef,
    pub options: collagene_core::codon_design::DesignOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SavedPrimerQuery {
    pub document: DocumentRef,
    /// Exact zero-based primer_bind feature index in this document revision.
    pub feature_index: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PrimerQuery {
    pub document: DocumentRef,
    pub primer: Primer,
    pub conditions: primers::ThermoConditions,
    pub binding_options: primers::BindingOptions,
    /// Optional second oligo for a contiguous complementarity screen.
    #[serde(default)]
    pub partner: Option<Primer>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PrimerDesignQuery {
    pub document: DocumentRef,
    pub options: primers::DesignOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TailingQuery {
    pub document: DocumentRef,
    pub recipe: verification::TailSubstrate,
    pub base: collagene_core::tailing::TailingBase,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OverlapProfileQuery {
    pub documents: Vec<DocumentRef>,
    pub reverse: Vec<bool>,
    pub overlaps: Vec<usize>,
    pub circular: bool,
    pub options: collagene_core::overlap_profiles::ProfileOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EndRepairQuery {
    pub document: DocumentRef,
    pub enzymes: Vec<String>,
    pub fragment_index: usize,
    pub method: chemistry::EndRepairMethod,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnnealQuery {
    pub top: DocumentRef,
    pub bottom: DocumentRef,
    pub options: chemistry::AnnealOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GoldenGateQuery {
    pub documents: Vec<DocumentRef>,
    pub options: collagene_core::golden_gate::GoldenGateOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MismatchPcrQuery {
    pub document: DocumentRef,
    pub forward: Primer,
    pub reverse: Primer,
    pub options: collagene_core::mutagenesis::PcrOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InversePcrQuery {
    pub document: DocumentRef,
    pub target: DocumentRef,
    #[serde(deserialize_with = "Option::deserialize")]
    #[schemars(required)]
    pub origin: Option<collagene_core::run_material::MaterialSelection>,
    #[serde(deserialize_with = "Option::deserialize")]
    #[schemars(required)]
    pub origin_workspace: Option<WorkspaceRef>,
    pub options: collagene_core::mutagenesis::InversePcrOptions,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BoundGelLane {
    pub label: String,
    pub document: DocumentRef,
    pub preparation: collagene_core::gel::Preparation,
    pub mass_ng: Option<f64>,
    /// Required for CatalogDigest: native retained owner/event/revision proof.
    #[serde(default)]
    pub workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GelQuery {
    pub lanes: Vec<BoundGelLane>,
    pub calibration: collagene_core::gel::Calibration,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DiagnosticQuery {
    pub documents: Vec<DocumentRef>,
    pub options: collagene_core::gel::DiagnosticOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportWorkspace {
    /// Base64 of a .cg2design snapshot exported by the web workspace; at most 8 MiB decoded.
    pub data_base64: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReadValidationQuery {
    pub reference: DocumentRef,
    pub reads: Vec<DocumentRef>,
    pub options: trace_validation::TraceValidationOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FastqImport {
    pub text: String,
    /// Explicit caller encoding; never inferred from overlapping ASCII ranges.
    pub encoding: oligo_formats::QualityEncoding,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FastqExport {
    pub document: DocumentRef,
    pub encoding: oligo_formats::QualityEncoding,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OligoImport {
    pub text: String,
    pub options: oligo_formats::OligoImportOptions,
    /// Explicit source data-row indices to create as DNA references. Unsupported
    /// modified rows remain collection evidence and cannot become DNA.
    #[serde(default)]
    pub selection: Vec<usize>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CollectionRef {
    pub id: String,
    pub revision: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OligoExport {
    pub collection: CollectionRef,
    pub delimiter: oligo_formats::Delimiter,
    pub mode: oligo_formats::TableExportMode,
    /// None exports all rows; an explicit list exports those source rows in order.
    #[serde(default)]
    pub selection: Option<Vec<usize>>,
}
/// A reference to immutable session evidence; exact revisions are always required.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EvidenceRef {
    pub id: String,
    pub revision: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LibraryFeatureSource {
    pub document: DocumentRef,
    pub feature_index: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CreateFeatureLibrary {
    pub name: String,
    pub version: String,
    pub features: Vec<LibraryFeatureSource>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportFeatureLibrary {
    pub data_base64: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FeatureLibraryQuery {
    pub library: EvidenceRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FeatureLibraryScan {
    pub document: DocumentRef,
    pub library: EvidenceRef,
    pub options: feature_library::ScanOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FeatureLibraryApply {
    pub document: DocumentRef,
    pub library: EvidenceRef,
    pub scan: EvidenceRef,
    pub selected: Vec<String>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnnotationImport {
    pub document: DocumentRef,
    pub text: String,
    /// Explicit annotation syntax: Gff3, Gtf or Bed. No name-based reference mapping.
    pub format: annotation_formats::AnnotationFormat,
    pub seqid: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnnotationApply {
    pub document: DocumentRef,
    pub dataset: EvidenceRef,
    pub selected: Vec<usize>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnnotationExport {
    pub document: DocumentRef,
    /// Explicit annotation syntax: Gff3, Gtf or Bed. No name-based reference mapping.
    pub format: annotation_formats::AnnotationFormat,
    pub seqid: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnnotationDatasetQuery {
    pub dataset: EvidenceRef,
}
#[derive(Debug, Clone, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EmptyRequest {}
#[derive(Debug, Clone, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DocumentPage {
    #[serde(default)]
    pub start: usize,
    /// Default 128; at most 128 summaries per response.
    #[serde(default)]
    pub limit: Option<usize>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BeginWorkspaceUpload {
    /// Exact decoded native archive byte length, at most 64 MiB.
    pub expected_bytes: usize,
    /// Lower-case BLAKE3 digest of the complete decoded archive.
    pub blake3: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WorkspaceUploadRef {
    pub upload_id: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AppendWorkspaceUpload {
    pub upload_id: String,
    /// Exact next decoded byte offset. No overlap, gaps or implicit replay.
    pub offset: usize,
    /// Independently base64-encoded chunk, at most 256 KiB decoded.
    pub data_base64: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WorkspaceRef {
    pub id: String,
    pub revision: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunMaterialOriginsQuery {
    pub workspace: WorkspaceRef,
    pub event: collagene_core::run_lineage::EventBinding,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialOriginQuery {
    pub workspace: WorkspaceRef,
    pub event: collagene_core::run_lineage::EventBinding,
    /// None discovers choices; Some explicitly replays that exact retained origin.
    pub selection: Option<collagene_core::run_material::MaterialSelection>,
}
/// MCP callers pass small immutable event selections, never closed capsule bytes.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum MaterialDerivationOrigin {
    OriginalReference,
    Saved(collagene_core::run_material::MaterialSelection),
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialPcrChoice {
    pub selected_index: usize,
    pub expected_preview_revision: String,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum MaterialDerivationOperation {
    ExactPcr {
        forward: Primer,
        reverse: Primer,
        selection: MaterialPcrChoice,
        template_strand: collagene_core::chemistry::Strand,
    },
    MismatchPcr {
        forward: Primer,
        reverse: Primer,
        options: collagene_core::mutagenesis::PcrOptions,
        selection: MaterialPcrChoice,
    },
    EndRepair {
        method: collagene_core::chemistry::EndRepairMethod,
    },
    CatalogDigest {
        selection: collagene_core::catalog_digest::CatalogSelection,
        options: collagene_core::restriction_inspection::InspectionOptions,
        fragment_index: usize,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialDerivationInstallRequest {
    pub source: collagene_core::run_lineage::EventBinding,
    pub origin: MaterialDerivationOrigin,
    pub derivation: MaterialDerivationOperation,
    pub name: String,
    pub save_primers: bool,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum MaterialDerivationTask {
    PreviewPcr {
        source: collagene_core::run_lineage::EventBinding,
        origin: MaterialDerivationOrigin,
        forward: Primer,
        reverse: Primer,
        template_strand: collagene_core::chemistry::Strand,
    },
    PreviewMismatchPcr {
        source: collagene_core::run_lineage::EventBinding,
        origin: MaterialDerivationOrigin,
        forward: Primer,
        reverse: Primer,
        options: collagene_core::mutagenesis::PcrOptions,
    },
    PreviewRepair {
        source: collagene_core::run_lineage::EventBinding,
        origin: MaterialDerivationOrigin,
        method: collagene_core::chemistry::EndRepairMethod,
    },
    PreviewDigest {
        source: collagene_core::run_lineage::EventBinding,
        origin: MaterialDerivationOrigin,
        selection: collagene_core::catalog_digest::CatalogSelection,
        options: collagene_core::restriction_inspection::InspectionOptions,
    },
    Install {
        request: MaterialDerivationInstallRequest,
    },
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialDerivationQuery {
    pub workspace: WorkspaceRef,
    pub task: MaterialDerivationTask,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunMaterialPreparationQuery {
    pub workspace: WorkspaceRef,
    pub selection: collagene_core::run_material::Selection,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RecordWorkspaceRun {
    pub workspace: WorkspaceRef,
    pub plan_revision: String,
    /// Explicit target event/document/revision binding; never inferred by sequence hash.
    pub target_binding: collagene_core::run_lineage::EventBinding,
    /// Explicit plan input IDs and their exact historical workspace events.
    pub input_bindings: Vec<collagene_core::run_lineage::InputEventBinding>,
    #[serde(default)]
    pub allow_incomplete: bool,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FinalizeCloningOutputs {
    pub workspace: WorkspaceRef,
    pub plan_revision: String,
    /// Exact independent target event, never selected by DNA equality.
    pub target_binding: collagene_core::run_lineage::EventBinding,
    /// Exact original event binding for each plan-local input ID.
    pub input_bindings: Vec<collagene_core::run_lineage::InputEventBinding>,
    #[serde(default)]
    pub options: collagene_core::cloning_outputs::Options,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReplayWorkspaceRun {
    pub workspace: WorkspaceRef,
    pub run_id: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportArchive {
    /// Base64 native .cg2workspace bytes; <=8 MiB decoded, else use upload tools.
    pub data_base64: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Region {
    pub document: DocumentRef,
    /// Zero-based inclusive start; circular regions must be split at the origin.
    pub start: usize,
    /// Zero-based exclusive end, at most 10,000 bases after start.
    pub end: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Restriction {
    pub document: DocumentRef,
    /// One name from the embedded starter catalog, e.g. EcoRI or BsaI.
    pub enzyme: String,
}
/// Inspect actual PCR-product digest geometry while keeping the original source immutable.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PcrDigestInspectionQuery {
    pub document: DocumentRef,
    pub request: collagene_core::pcr_digest_inspection::Request,
    /// Explicit mismatch-PCR policy; omitted keeps the exact-PCR inspector unchanged.
    #[serde(default)]
    pub pcr_options: Option<collagene_core::mutagenesis::PcrOptions>,
    /// Explicit documentary origin for a saved-run source; never inferred from metadata.
    #[serde(default)]
    pub origin: Option<collagene_core::run_material::MaterialOrigin>,
    /// Required exact retained native workspace when a saved-run origin is supplied.
    #[serde(default)]
    pub workspace: Option<WorkspaceRef>,
    #[serde(default)]
    pub fragment_start: usize,
    #[serde(default = "default_methylation_page")]
    pub max_fragments: usize,
}
/// Inspect the actual final outer-PCR product of a complete physical OE preparation.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OverlapExtensionDigestInspectionQuery {
    pub document: DocumentRef,
    pub request: collagene_core::overlap_extension_digest_inspection::Request,
    #[serde(default)]
    pub origin: Option<collagene_core::run_material::MaterialOrigin>,
    /// Required exact retained native workspace when an origin is supplied.
    #[serde(default)]
    pub workspace: Option<WorkspaceRef>,
    #[serde(default)]
    pub fragment_start: usize,
    #[serde(default = "default_methylation_page")]
    pub max_fragments: usize,
}
/// Complete catalog-bound restriction evidence from one immutable source.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "source", deny_unknown_fields)]
pub enum RestrictionInspectionEvidence {
    Original {
        inspection: collagene_core::catalog_digest::CatalogRestrictionInspection,
    },
    SavedRun {
        inspection: collagene_core::run_catalog_digest::Inspection,
    },
}
/// Bounded recognition evidence and assumed-complete digest summaries for one immutable source.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RestrictionInspectionQuery {
    pub document: DocumentRef,
    /// Stable source role; included in every report binding independently of document identity.
    pub source_role: String,
    /// Explicit saved material source; requires its exact native workspace.
    #[serde(default)]
    pub origin: Option<collagene_core::run_material::Capsule>,
    #[serde(default)]
    pub workspace: Option<WorkspaceRef>,
    /// One to 64 exact catalog enzyme IDs or uniquely resolved names/aliases. Formulations are separate.
    pub enzymes: Vec<String>,
    /// Optional exact current catalog edition/hash/ordered-ID/geometry pin; stale or changed pins fail.
    pub catalog_pin: Option<collagene_core::catalog_digest::CatalogSelection>,
    #[serde(default)]
    pub options: collagene_core::restriction_inspection::InspectionOptions,
    #[serde(default)]
    pub site_start: usize,
    #[serde(default = "default_methylation_page")]
    pub max_sites: usize,
    #[serde(default)]
    pub fragment_start: usize,
    #[serde(default = "default_methylation_page")]
    pub max_fragments: usize,
}
/// Recognition/cut-coordinate evidence only; never constructs a physical digest.
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RestrictionScanQuery {
    pub document: DocumentRef,
    /// Stable explicit source role, independently bound to source revision.
    pub source_role: String,
    /// One to 64 exact catalog IDs or uniquely resolved names/aliases, in order.
    pub enzymes: Vec<String>,
    /// Optional exact edition/hash/ordered-ID/geometry pin; changed pins fail.
    pub catalog_pin: Option<collagene_core::catalog_digest::CatalogSelection>,
    #[serde(default)]
    pub options: collagene_core::restriction_inspection::InspectionOptions,
    #[serde(default)]
    pub site_start: usize,
    #[serde(default = "default_methylation_page")]
    pub max_sites: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DigestQuery {
    pub document: DocumentRef,
    pub enzymes: Vec<String>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Pcr {
    pub document: DocumentRef,
    pub forward: Primer,
    pub reverse: Primer,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TranslateQuery {
    pub document: DocumentRef,
    pub options: translation::TranslationOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CdsQuery {
    pub document: DocumentRef,
    pub feature_index: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OrfQuery {
    pub document: DocumentRef,
    pub options: translation::OrfOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AlignQuery {
    pub reference: DocumentRef,
    pub query: DocumentRef,
    pub options: alignment::AlignmentOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CompareQuery {
    pub reference: DocumentRef,
    pub query: DocumentRef,
    pub options: alignment::ComparisonOptions,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BoundInput {
    /// Plan-local source ID used by steps, e.g. left or template.
    pub id: String,
    pub document: DocumentRef,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Proposal {
    pub name: String,
    /// Independent desired target. Never inferred from a candidate product.
    pub target: DocumentRef,
    pub inputs: Vec<BoundInput>,
    /// Typed operations refer to input IDs or earlier steps. OverlapExtension explicitly denatures/reanneals two intact linear tracked duplexes and extends their selected terminal pairing; use actual PCR oligos before and after it. Reference overlap assemblies cannot substitute for physical PCR templates.
    pub steps: Vec<Operation>,
}
/// MCP accepts small exact native selections, never a caller's closed capsule.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct IntentMaterialBinding {
    pub input_id: String,
    pub event: collagene_core::run_lineage::EventBinding,
    pub selection: collagene_core::run_material::MaterialSelection,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AssemblyIntentRequest {
    pub name: String,
    pub fragments: Vec<collagene_core::assembly_intent::Contribution>,
    pub circular: bool,
    pub material_bindings: Vec<IntentMaterialBinding>,
}
impl AssemblyIntentRequest {
    fn native(
        &self,
        workspace: Option<&collagene_core::workspace::Workspace>,
    ) -> Result<collagene_core::assembly_intent::Request, ErrorData> {
        bounded_json_limit(self, 2 * 1024 * 1024)?;
        for f in &self.fragments {
            if f.role.trim().is_empty()
                || f.role.len() > 128
                || f.input_id.trim().is_empty()
                || f.input_id.len() > 128
            {
                return Err(bad("bounded source contribution role/input required"));
            }
        }
        if self.fragments.is_empty()
            || self.fragments.len() > 32
            || self.material_bindings.len() > 32
            || self.name.len() > 1024
        {
            return Err(bad(
                "bounded explicit contributions and material bindings required",
            ));
        }
        let mut ids = BTreeSet::new();
        for b in &self.material_bindings {
            if b.input_id.trim().is_empty()
                || b.input_id.len() > 128
                || !ids.insert(b.input_id.as_str())
                || b.selection.event() != &b.event
            {
                return Err(bad(
                    "material selection must name its exact unique input/event",
                ));
            }
            b.selection.validate().map_err(bad)?;
        }
        let origins = if self.material_bindings.is_empty() {
            vec![]
        } else {
            let workspace = workspace
                .ok_or_else(|| bad("native material selection requires origin_workspace"))?;
            let selections = self
                .material_bindings
                .iter()
                .map(|b| b.selection.clone())
                .collect::<Vec<_>>();
            collagene_core::run_material::documentary_origins(workspace, &selections, &[])
                .map_err(bad)?
                .materials
        };
        Ok(collagene_core::assembly_intent::Request {
            name: self.name.clone(),
            fragments: self.fragments.clone(),
            circular: self.circular,
            material_bindings: self
                .material_bindings
                .iter()
                .zip(origins)
                .map(
                    |(binding, origin)| collagene_core::assembly_intent::MaterialBinding {
                        input_id: binding.input_id.clone(),
                        event: binding.event.clone(),
                        origin,
                    },
                )
                .collect(),
        })
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SubstitutedIntentRequest {
    pub base: AssemblyIntentRequest,
    pub junction_edits: Vec<collagene_core::assembly_junction_edit::Edit>,
    pub substitutions: Vec<collagene_core::assembly_substitution_intent::Substitution>,
}
impl SubstitutedIntentRequest {
    fn native(
        &self,
        w: Option<&collagene_core::workspace::Workspace>,
    ) -> Result<collagene_core::assembly_substitution_intent::Request, ErrorData> {
        bounded_json_limit(self, 2 * 1024 * 1024)?;
        if self.junction_edits.len() > 32 || self.substitutions.len() > 256 {
            return Err(bad("intent edit bound"));
        }
        Ok(collagene_core::assembly_substitution_intent::Request {
            base: self.base.native(w)?,
            junction_edits: self.junction_edits.clone(),
            substitutions: self.substitutions.clone(),
        })
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ContributionEditIntentRequest {
    pub base: SubstitutedIntentRequest,
    pub edits: Vec<collagene_core::assembly_edit_intent::Edit>,
}
impl ContributionEditIntentRequest {
    fn native(
        &self,
        w: Option<&collagene_core::workspace::Workspace>,
    ) -> Result<collagene_core::assembly_edit_intent::Request, ErrorData> {
        bounded_json_limit(self, 2 * 1024 * 1024)?;
        if self.edits.len() > 256 {
            return Err(bad("intent edit bound"));
        }
        Ok(collagene_core::assembly_edit_intent::Request {
            base: self.base.native(w)?,
            edits: self.edits.clone(),
        })
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AssemblyIntentQuery {
    /// Exact original source snapshots; local IDs match request fragments.
    pub inputs: Vec<BoundInput>,
    /// Declared source regions, order and orientation only. No simulator-derived target.
    pub request: AssemblyIntentRequest,
    /// Pinned native archive owning all exact material selections; null for references.
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OverlapPrimerDesignQuery {
    /// Actual original templates with exact immutable source revisions.
    pub inputs: Vec<BoundInput>,
    /// Independently supplied intended DNA; never manufactured from chosen intervals.
    pub target: DocumentRef,
    /// Explicit ordered template intervals, orientation and bounded heuristic controls.
    pub request: collagene_core::overlap_primer_design::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OverlapExtensionPrimerDesignQuery {
    /// Immutable original or exact saved native source snapshot; saved sources require their owning workspace.
    pub document: DocumentRef,
    /// Construct-oriented contribution in the explicit stored-strand frame, independently declared edits and bounded primer policy.
    pub request: collagene_core::overlap_extension_primer_design::Request,
    #[serde(deserialize_with = "Option::deserialize")]
    #[schemars(required)]
    pub origin_workspace: Option<WorkspaceRef>,
}

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProtocolLayoutQuery {
    pub workspace: WorkspaceRef,
    pub task: ProtocolLayoutTask,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum ProtocolLayoutTask {
    Preview {
        request: collagene_core::protocol_layout::Request,
    },
    Capture {
        request: collagene_core::protocol_layout::Request,
        expected_context_revision: String,
    },
    List {
        offset: usize,
        limit: usize,
    },
    Open {
        layout_hash: String,
    },
    Export {
        layout_hash: String,
        format: collagene_core::protocol_layout::ExportFormat,
    },
    ReviewPublication {
        layout_hash: String,
    },
    ProjectPublication {
        layout_hash: String,
        expected_scope_revision: String,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MaterialOriginWorkspace {
    pub role: String,
    pub workspace: WorkspaceRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct HomologyRecipePreviewQuery {
    pub inputs: Vec<BoundInput>,
    pub target: DocumentRef,
    pub recipe: collagene_core::assembly_recipe::Recipe,
    /// Explicit native strand resection/fill/sealing policy; not inferred from a vendor label.
    pub policy: collagene_core::homology_assembly::Options,
    #[serde(default)]
    pub origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AnnealingRecipePreviewQuery {
    pub inputs: Vec<BoundInput>,
    pub target: DocumentRef,
    pub recipe: collagene_core::assembly_recipe::Recipe,
    /// Explicit native 3′ resection/annealing policy, with no fill or sealing; not inferred from a vendor label.
    pub policy: collagene_core::homology_annealing::Options,
    #[serde(default)]
    pub origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditedOverlapPrimerRequest {
    pub intent: ContributionEditIntentRequest,
    pub overlap: usize,
    pub options: collagene_core::overlap_extension_primer_design::Options,
    pub profile_options: collagene_core::overlap_profiles::ProfileOptions,
    pub finishing: collagene_core::assembly_recipe::OverlapFinishing,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditedOverlapPrimerDesignQuery {
    pub inputs: Vec<BoundInput>,
    pub request: EditedOverlapPrimerRequest,
    /// One immutable native archive owning every selected saved-source role; null for originals.
    #[serde(deserialize_with = "Option::deserialize")]
    #[schemars(required)]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SubstitutedAssemblyIntentQuery {
    /// Exact immutable MCP source snapshots, including retained historical revisions.
    pub inputs: Vec<BoundInput>,
    /// Independent regions, junction edits and expected-base substitutions; never primer-derived.
    pub request: SubstitutedIntentRequest,
    /// Pinned native archive owning all exact material selections; null for references.
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ContributionEditIntentQuery {
    /// Exact immutable original source snapshots, never products or guessed native parents.
    pub inputs: Vec<BoundInput>,
    /// Independent construct-oriented original contribution intervals and literal edits.
    pub request: ContributionEditIntentRequest,
    /// Pinned native archive owning all exact material selections; null for references.
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditedAssemblyIntentQuery {
    pub inputs: Vec<BoundInput>,
    pub base: AssemblyIntentRequest,
    pub edits: Vec<collagene_core::assembly_junction_edit::Edit>,
    /// Pinned native archive owning all exact material selections; null for references.
    #[serde(deserialize_with = "Option::deserialize")]
    pub origin_workspace: Option<WorkspaceRef>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditedAssemblyPrimersQuery {
    pub inputs: Vec<BoundInput>,
    pub request: collagene_core::assembly_junction_edit::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SelectiveAssemblyPrimersQuery {
    pub inputs: Vec<BoundInput>,
    pub request: collagene_core::selective_junction_edit::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GoldenGateJunctionPrimersQuery {
    pub inputs: Vec<BoundInput>,
    pub request: collagene_core::golden_gate_junction_edit::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReadReviewQuery {
    pub workspace: WorkspaceRef,
    pub request: collagene_core::read_review::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RecordReadCheck {
    pub workspace: WorkspaceRef,
    pub reference: collagene_core::run_lineage::EventBinding,
    pub reads: Vec<collagene_core::run_lineage::EventBinding>,
    pub options: collagene_core::trace_validation::TraceValidationOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ListReadChecks {
    pub workspace: WorkspaceRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReplayReadCheck {
    pub workspace: WorkspaceRef,
    pub attachment_hash: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MutationIntentQuery {
    pub document: DocumentRef,
    pub request: collagene_core::mutation_intent::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GoldenGatePrimerDesignQuery {
    pub inputs: Vec<BoundInput>,
    pub target: DocumentRef,
    pub request: collagene_core::golden_gate_primer_design::Request,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OverlapAssemblyProposal {
    /// Exact immutable source documents, identified by recipe-local input IDs.
    pub inputs: Vec<BoundInput>,
    /// Independently supplied intended DNA; never inferred from prepared fragments.
    pub target: DocumentRef,
    /// Ordered preparations, orientation, explicit homology lengths and method profile.
    pub recipe: collagene_core::assembly_recipe::Recipe,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PreparedGoldenGatePoolQuery {
    pub inputs: Vec<BoundInput>,
    pub target: DocumentRef,
    pub request: collagene_core::prepared_golden_gate::Request,
    pub origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PreparedGoldenGateSelectionQuery {
    pub artifact_id: String,
    pub selection: collagene_core::prepared_golden_gate::Selection,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GoldenGateReactionQuery {
    pub inputs: Vec<BoundInput>,
    pub target: DocumentRef,
    pub request: collagene_core::golden_gate_reaction::Request,
    pub origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct GoldenGateReactionEvidence {
    kind: String,
    inputs: Vec<BoundInput>,
    target: DocumentRef,
    report: collagene_core::golden_gate_reaction::Report,
    origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct PreparedGoldenGateEvidence {
    kind: String,
    inputs: Vec<BoundInput>,
    target: DocumentRef,
    report: collagene_core::prepared_golden_gate::Report,
    origin_workspaces: Vec<MaterialOriginWorkspace>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LigationAssemblyProposal {
    /// Original immutable DNA sources; preparations never replace these inputs.
    pub inputs: Vec<BoundInput>,
    /// Independent desired reference, not a sequence inferred from a ligated product.
    pub target: DocumentRef,
    /// Explicit ordered physical preparations, end treatments and ligation policy.
    pub recipe: collagene_core::ligation_recipe::Recipe,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TailingAssemblyProposal {
    /// Original immutable DNA sources; preparations never replace these inputs.
    pub inputs: Vec<BoundInput>,
    /// Independent desired reference, not a sequence inferred from a ligated product.
    pub target: DocumentRef,
    /// Explicit TA/GC profile, original-source preparations, per-row keep/add-tail actions, terminal treatments and ligation policy.
    pub recipe: collagene_core::tailing_recipe::Recipe,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PlanRef {
    /// Exact immutable plan revision returned by propose_plan.
    pub plan_revision: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AssemblyJunctionQuery {
    /// Exact immutable saved candidate; its separately bound target is reused.
    pub plan_revision: String,
    /// Explicit bounded sequence flanks and genetic-code inspection policy.
    pub options: collagene_core::assembly_junctions::Options,
}
#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
pub enum PlanOligoFormat {
    Csv,
    Tsv,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PlanOligoExport {
    pub plan_revision: String,
    pub format: PlanOligoFormat,
    #[serde(default)]
    pub options: collagene_core::plan_oligos::Options,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ArtifactRegion {
    pub artifact_id: String,
    /// UTF-8 byte offset, must be on a character boundary.
    pub start: usize,
    /// Maximum 16,384 UTF-8 bytes; next_start is returned for pagination.
    pub max_bytes: usize,
}
#[derive(Clone)]
struct StoredPlan {
    plan: DesignPlan,
    target: Molecule,
}
enum ArtifactData {
    Text(String),
    Binary(Vec<u8>),
}
struct Artifact {
    data: ArtifactData,
    mime_type: &'static str,
    encoding: &'static str,
    filename: &'static str,
}
#[derive(Clone, Serialize)]
struct StoredCandidate {
    plan_revision: String,
    name: String,
    target_revision: String,
    inputs: Vec<(String, String)>,
    steps: Vec<Operation>,
}
#[derive(Clone, Serialize)]
struct StoredWorkspace {
    name: String,
    active_document: String,
    documents: Vec<workspace::WorkspaceDocument>,
    revision_ids: Vec<String>,
    history: Vec<workspace::RevisionEvent>,
    candidates: Vec<StoredCandidate>,
    active_candidate: Option<String>,
    attachments: BTreeMap<String, workspace::WorkspaceAttachment>,
}
impl StoredWorkspace {
    fn split(value: workspace::Workspace) -> (Self, BTreeMap<String, Molecule>) {
        let revisions = value.revisions;
        (
            Self {
                name: value.name,
                active_document: value.active_document,
                documents: value.documents,
                revision_ids: revisions.keys().cloned().collect(),
                history: value.history,
                candidates: value
                    .candidates
                    .into_iter()
                    .map(|c| StoredCandidate {
                        plan_revision: c.plan_revision,
                        name: c.plan.name,
                        target_revision: c.plan.target_revision,
                        inputs: c
                            .plan
                            .inputs
                            .into_iter()
                            .map(|i| (i.id, i.revision))
                            .collect(),
                        steps: c.plan.steps,
                    })
                    .collect(),
                active_candidate: value.active_candidate,
                attachments: value.attachments,
            },
            revisions,
        )
    }
    fn restore(&self, store: &Store) -> Result<workspace::Workspace, ErrorData> {
        let molecule = |revision: &String| {
            store
                .documents
                .get(revision)
                .cloned()
                .ok_or_else(|| bad("workspace revision unavailable"))
        };
        Ok(workspace::Workspace {
            name: self.name.clone(),
            active_document: self.active_document.clone(),
            documents: self.documents.clone(),
            revisions: self
                .revision_ids
                .iter()
                .map(|r| Ok((r.clone(), molecule(r)?)))
                .collect::<Result<_, ErrorData>>()?,
            history: self.history.clone(),
            active_candidate: self.active_candidate.clone(),
            attachments: self.attachments.clone(),
            candidates: self
                .candidates
                .iter()
                .map(|c| {
                    Ok(transfer::Candidate {
                        plan_revision: c.plan_revision.clone(),
                        plan: DesignPlan {
                            name: c.name.clone(),
                            target_revision: c.target_revision.clone(),
                            steps: c.steps.clone(),
                            inputs: c
                                .inputs
                                .iter()
                                .map(|(id, r)| {
                                    Ok(PlanInput {
                                        id: id.clone(),
                                        revision: r.clone(),
                                        molecule: molecule(r)?,
                                    })
                                })
                                .collect::<Result<_, ErrorData>>()?,
                        },
                    })
                })
                .collect::<Result<_, ErrorData>>()?,
        })
    }
    fn summary(&self, id: &str) -> Value {
        json!({"id":id,"revision":id,"name":self.name,"active_document":self.active_document,
            "document_count":self.documents.len(),"revision_count":self.revision_ids.len(),"history_event_count":self.history.len(),
            "attachment_count":self.attachments.len(),"active_candidate":self.active_candidate,
            "candidate_count":self.candidates.len(),"candidate_summaries_truncated":self.candidates.len()>20,"candidates":self.candidates.iter().take(20).map(|c|json!({"plan_revision":c.plan_revision,"name":c.name.chars().take(128).collect::<String>(),"target_revision":c.target_revision,"state":"candidate","verified":false})).collect::<Vec<_>>(),
            "authority":"Archived candidates are preserved data. To compute a plan, explicitly rebind its inputs and independent target with propose_plan; import does not assert verification."})
    }
}
#[derive(Clone, Serialize)]
struct StoredSelection {
    collection_id: String,
    source_rows: Vec<usize>,
}
struct WorkspaceUpload {
    expected: usize,
    hash: String,
    bytes: Vec<u8>,
}
#[derive(Default)]
struct Store {
    documents: BTreeMap<String, Molecule>,
    plans: BTreeMap<String, StoredPlan>,
    artifacts: BTreeMap<String, Artifact>,
    used: usize,
    workspaces: BTreeMap<String, StoredWorkspace>,
    oligos: BTreeMap<String, oligo_formats::OligoTable>,
    uploads: BTreeMap<String, WorkspaceUpload>,
    oligo_selections: BTreeMap<String, StoredSelection>,
    upload_serial: u64,
    feature_libraries: BTreeMap<String, feature_library::FeatureLibrary>,
    feature_scans: BTreeMap<String, feature_library::ScanReport>,
    annotation_datasets: BTreeMap<String, annotation_formats::AnnotationImportReport>,
    alignments: BTreeMap<String, Arc<alignment_formats::AlignmentDocument>>,
    read_alignments: BTreeMap<String, Arc<collagene_core::read_alignment_store::Collection>>,
    variants: BTreeMap<String, Arc<collagene_core::variant_reference::exchange::Collection>>,
    contig_reports: BTreeMap<String, String>,
}
impl Store {
    fn variants(
        &self,
        reference: &EvidenceRef,
    ) -> Result<Arc<collagene_core::variant_reference::exchange::Collection>, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale variant collection revision"));
        }
        self.variants
            .get(&reference.id)
            .cloned()
            .ok_or_else(|| bad("unknown variant collection"))
    }
    fn read_alignments(
        &self,
        reference: &EvidenceRef,
    ) -> Result<Arc<collagene_core::read_alignment_store::Collection>, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale read alignment collection revision"));
        }
        self.read_alignments
            .get(&reference.id)
            .cloned()
            .ok_or_else(|| bad("unknown read alignment collection"))
    }
    fn alignment(
        &self,
        reference: &EvidenceRef,
    ) -> Result<Arc<alignment_formats::AlignmentDocument>, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale alignment revision"));
        }
        self.alignments
            .get(&reference.id)
            .cloned()
            .ok_or_else(|| bad("unknown alignment"))
    }
    fn document(&self, reference: &DocumentRef) -> Result<&Molecule, ErrorData> {
        let mol = self
            .documents
            .get(&reference.id)
            .ok_or_else(|| bad("unknown document ID"))?;
        if reference.id != reference.revision {
            return Err(bad("stale document revision"));
        }
        Ok(mol)
    }
    fn feature_library(
        &self,
        reference: &EvidenceRef,
    ) -> Result<&feature_library::FeatureLibrary, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale feature library revision"));
        }
        self.feature_libraries
            .get(&reference.id)
            .ok_or_else(|| bad("unknown feature library"))
    }
    fn feature_scan(
        &self,
        reference: &EvidenceRef,
    ) -> Result<&feature_library::ScanReport, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale feature scan revision"));
        }
        self.feature_scans
            .get(&reference.id)
            .ok_or_else(|| bad("unknown feature scan"))
    }
    fn annotation_dataset(
        &self,
        reference: &EvidenceRef,
    ) -> Result<&annotation_formats::AnnotationImportReport, ErrorData> {
        if reference.id != reference.revision {
            return Err(bad("stale annotation dataset revision"));
        }
        self.annotation_datasets
            .get(&reference.id)
            .ok_or_else(|| bad("unknown annotation dataset"))
    }
    fn reserve(&mut self, bytes: usize) -> Result<(), ErrorData> {
        if bytes > MAX_STORE_BYTES.saturating_sub(self.used) {
            return Err(bad(
                "session storage limit reached; restart or use fewer documents/plans",
            ));
        }
        self.used += bytes;
        Ok(())
    }
}
fn bad(message: impl ToString) -> ErrorData {
    ErrorData::invalid_params(message.to_string(), None)
}
// Bound temporary JSON materialization before session accounting, including binary arrays.
fn bounded_json<T: Serialize>(value: &T) -> Result<String, ErrorData> {
    bounded_json_limit(value, MAX_STORE_BYTES)
}
fn bounded_json_limit<T: Serialize>(value: &T, cap: usize) -> Result<String, ErrorData> {
    struct Sink {
        bytes: Vec<u8>,
        cap: usize,
    }
    impl std::io::Write for Sink {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if bytes.len() > self.cap.saturating_sub(self.bytes.len()) {
                return Err(std::io::Error::other(
                    "JSON artifact exceeds session budget",
                ));
            }
            self.bytes.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut sink = Sink {
        bytes: Vec::new(),
        cap,
    };
    serde_json::to_writer(&mut sink, value).map_err(bad)?;
    String::from_utf8(sink.bytes).map_err(bad)
}
fn topo_summary(report: &collagene_core::topo::physical::Report) -> Value {
    json!({"engine":report.engine,"model":report.model,"physical_batch_revision":report.batch_revision,
        "work_complete":report.work_complete,"all_requested_steps_applied":report.all_requested_steps_applied,
        "unsupported_or_unknown":report.unsupported_or_unknown,"fully_specified_dna":report.fully_specified_dna,
        "step_count":report.steps.len(),"component_count":report.product.components.len(),"chain_count":report.product.chains.len(),
        "nick_count":report.product.nicks.len(),"unpaired_run_count":report.product.unpaired.len(),"released_activation_count":report.product.released_activations.len(),
        "components":report.product.components.iter().take(64).enumerate().map(|(i,c)|json!({"component":i,"chains":c.chains,"paired_nucleotides":c.paired_nucleotides,"unpaired_nucleotides":c.unpaired_nucleotides,"covalently_closed":c.covalently_closed})).collect::<Vec<_>>(),
        "components_truncated":report.product.components.len()>64,
        "chain_lengths":report.product.chains.iter().take(64).map(|c|c.sequence.len()).collect::<Vec<_>>(),
        "chains_truncated":report.product.chains.len()>64})
}
fn topo_annotation_summary(
    evidence: Option<&collagene_core::topo_annotations::AnnotationEvidence>,
) -> Value {
    match evidence {
        None => {
            json!({"state":"withheld_without_paired_reference","authority":"No annotation coordinates or chemical state inferred."})
        }
        Some(evidence) => {
            json!({"state":"source_correspondence_only","features":evidence.features.len(),"sources":evidence.sources.iter().map(|source|json!({"role":source.role,"revision":source.source_revision,"projection_supported":source.projection_supported,"projected":source.projected,"dropped_from_selected_reference":source.dropped,"unmapped_original_features":source.unmapped.len(),"warnings":source.warnings.iter().take(8).collect::<Vec<_>>() })).collect::<Vec<_>>(),"authority":"Review disrupted/PCR-derived annotations; retained source features are not chemical ancestry or restored physical authority."})
        }
    }
}
fn topo_projection_summary(product: &collagene_core::topo::Product) -> Value {
    use collagene_core::topo::physical::projection::Outcome;
    match &product.projection.outcome {
        Outcome::Proven(view) => {
            json!({"state":"paired_contour","topology":view.topology,"length":view.top.len(),
            "nick_count":view.nicks.len(),"fully_known_terminal_chemistry":view.fully_known_terminal_chemistry,
            "all_requested_steps_applied":view.all_requested_steps_applied,
            "reference_available":product.reference.is_some(),"source_map_count":view.source_maps.len()})
        }
        Outcome::Withheld(reason) => {
            json!({"state":"withheld","reason":reason,"reference_available":false})
        }
    }
}
fn output(value: Value) -> Result<CallToolResult, ErrorData> {
    Ok(CallToolResult::structured(value))
}
fn capture_variant(
    store: &Store,
    args: &VariantCheck,
) -> Result<
    (
        Arc<collagene_core::variant_reference::exchange::Collection>,
        Molecule,
        collagene_core::variant_reference::Options,
    ),
    ErrorData,
> {
    if args.reference_name.trim().is_empty()
        || args.reference_name.len() > 1024
        || args.reference_name.chars().any(char::is_control)
        || !(1..=64).contains(&args.alternate_index)
    {
        return Err(bad(
            "explicit bounded contig name and 1-based ALT 1..64 required",
        ));
    }
    let collection = store.variants(&args.collection)?;
    let reference = store.document(&args.reference)?;
    collagene_core::model::bounds::model_preflight(&[reference]).map_err(bad)?;
    if collection.record(args.record_index).is_none() {
        return Err(bad("variant record index unavailable"));
    }
    let options = collagene_core::variant_reference::Options {
        source_blake3: collection.source_blake3().into(),
        record_index: args.record_index,
        reference_name: args.reference_name.clone(),
        expected_reference_revision: args.reference.revision.clone(),
        alternate_index: args.alternate_index,
    };
    Ok((collection, reference.clone(), options))
}
fn variant_summary(collection: &collagene_core::variant_reference::exchange::Collection) -> Value {
    json!({"collection":{"id":collection.source_blake3(),"revision":collection.source_blake3()},
        "source_blake3":collection.source_blake3(),"source_bytes":collection.original_bytes().len(),
        "format":match collection.format(){collagene_core::variant_reference::exchange::Format::Vcf=>"vcf",collagene_core::variant_reference::exchange::Format::Bcf=>"bcf"},
        "record_count":collection.len(),"sample_column_count":collection.header().sample_names().len(),
        "sample_columns":collection.header().sample_names().iter().take(32).map(|s|read_name(Some(s))).collect::<Vec<_>>(),
        "sample_columns_truncated":collection.header().sample_names().len()>32,
        "authority":"Declared file evidence only; header/build/sample/genotype/INFO assertions are unverified. No implicit DNA reference mapping, variant calling or clinical inference."})
}
fn read_alignment_summary(collection: &collagene_core::read_alignment_store::Collection) -> Value {
    let summary = collection.summary();
    json!({"collection":{"id":summary.source_blake3,"revision":summary.source_blake3},
        "source_blake3":summary.source_blake3,"source_bytes":summary.source_bytes,"format":summary.format,
        "record_count":summary.records,"reference_count":summary.references.len(),
        "flag_unmapped":summary.flag_unmapped,"flag_secondary":summary.flag_secondary,"flag_supplementary":summary.flag_supplementary,
        "retained_charge":collection.retained_charge(),"file_complete":true,
        "warnings":summary.warnings.iter().take(20).map(|s|s.chars().take(512).collect::<String>()).collect::<Vec<_>>(),
        "warnings_omitted":summary.warnings.len().saturating_sub(20),"interpretation":summary.interpretation})
}
fn read_name(value: Option<&str>) -> Option<String> {
    value.map(|s| s.chars().take(256).collect())
}
fn summary(id: &str, m: &Molecule) -> Value {
    let source = m
        .metadata
        .iter()
        .filter(|(k, _)| {
            [
                "source:provider",
                "source:canonical_url",
                "source:identifier",
                "source:accession_version",
                "source:retrieved_at",
                "source:content_blake3",
                "source:format",
                "source:sequence_id",
                "source:addgene_origin",
                "source:addgene_completeness",
            ]
            .contains(&k.as_str())
        })
        .collect::<Vec<_>>();
    json!({"id":id,"revision":id,"name":m.name.chars().take(256).collect::<String>(),"length":m.sequence.len(),"alphabet":m.sequence.alphabet(),"topology":m.topology,"feature_count":m.features.len(),"source":source.iter().map(|(k,v)|(k.as_str(),v.chars().take(256).collect::<String>())).collect::<BTreeMap<_,_>>(),"source_values_truncated":source.iter().any(|(_,v)|v.chars().count()>256)})
}

fn alignment_summary(id: &str, document: &alignment_formats::AlignmentDocument) -> Value {
    json!({"alignment":{"id":id,"revision":id},"name":document.name.chars().take(256).collect::<String>(),"alphabet":document.alphabet,"row_count":document.rows.len(),"column_count":document.column_count(),"source_format":document.source.format,"source_blake3":document.source.blake3,"rows":document.rows.iter().take(100).map(|r|json!({"id":r.id,"identifier":r.identifier.chars().take(256).collect::<String>(),"description":r.description.chars().take(256).collect::<String>(),"text_preview_truncated":r.identifier.chars().count()>256||r.description.chars().count()>256})).collect::<Vec<_>>(),"rows_preview_truncated":document.rows.len()>100,"persistence":"session only","semantics":"Imported alignment evidence, not a newly computed alignment or verified external reference coordinate map."})
}

struct NativeSearch {
    genome: collagene_core::genome_session::GenomeRef,
    state: collagene_core::genome_crispr::GenomeSearchSession,
}
#[derive(Default)]
struct GenomeResources {
    native: Option<genome_store::NativeGenomes>,
    searches: BTreeMap<String, NativeSearch>,
    serial: u64,
}
impl GenomeResources {
    fn native(&mut self) -> Result<&mut genome_store::NativeGenomes, ErrorData> {
        if self.native.is_none() {
            self.native = Some(genome_store::NativeGenomes::new().map_err(bad)?);
        }
        Ok(self.native.as_mut().unwrap())
    }
    fn existing(&mut self) -> Result<&mut genome_store::NativeGenomes, ErrorData> {
        self.native
            .as_mut()
            .ok_or_else(|| bad("no live large-reference session"))
    }
}
fn native_ref(d: &genome_store::NativeDataset) -> collagene_core::genome_session::GenomeRef {
    collagene_core::genome_session::GenomeRef {
        id: d.id.clone(),
        manifest_hash: d.reference.manifest_hash.clone(),
        content_revision: d.reference.manifest.content_revision.clone(),
    }
}
fn native_summary(d: &genome_store::NativeDataset) -> Value {
    let m = &d.reference.manifest;
    json!({"genome":native_ref(d),"name":m.build.name.chars().take(256).collect::<String>(),"assembly_accession":m.build.assembly_accession,"version":m.build.version,"contig_count":m.contigs.len(),"total_bases":m.total_bases,"input_bytes":d.input_bytes,"file_length":d.reference.file_length,"semantic_verified_at_import":d.semantic_verified,"expected_components_complete":m.expected_complete,"missing_expected_count":m.missing_expected.len(),"raw_source_hash":m.raw_source_hash,"storage":"Session-owned managed scratchfile; no arbitrary filesystem path or persistent workspace write."})
}
#[derive(Clone)]
pub struct CollageneServer {
    store: Arc<Mutex<Store>>,
    tool_router: ToolRouter<Self>,
    compute_gate: Arc<tokio::sync::Semaphore>,
    source_client: Arc<tokio::sync::OnceCell<collagene_sources::SourceClient>>,
    genomes: Arc<Mutex<GenomeResources>>,
}
impl CollageneServer {
    fn topo_sources(
        &self,
        sources: &[TopoSourceRef],
        request: &collagene_core::topo::Request,
        target: Option<&DocumentRef>,
    ) -> Result<(Vec<collagene_core::topo::Input>, Option<Molecule>), ErrorData> {
        if sources.len() > 32
            || request.recipe.materials.is_empty()
            || request.recipe.materials.len() > 16
            || request.steps.len() > 16
        {
            return Err(bad(
                "TOPO requires bounded explicit sources and current native recipe/materials/steps",
            ));
        }
        let mut roles = BTreeSet::new();
        for source in sources {
            if source.role.trim().is_empty()
                || source.role.len() > 256
                || source.role.chars().any(char::is_control)
                || !roles.insert(source.role.as_str())
            {
                return Err(bad(
                    "TOPO roles must be unique nonblank bounded identifiers",
                ));
            }
        }
        let store = self.store.lock().unwrap();
        let borrowed = sources
            .iter()
            .map(|s| store.document(&s.document))
            .collect::<Result<Vec<_>, _>>()?;
        let target = target.map(|t| store.document(t)).transpose()?;
        let mut aggregate = borrowed.clone();
        aggregate.extend(target);
        collagene_core::model::bounds::model_preflight(&aggregate).map_err(bad)?;
        let mut bases = 0usize;
        for molecule in &borrowed {
            if molecule.sequence.alphabet() != collagene_core::Alphabet::Dna {
                return Err(bad("TOPO requires explicit DNA sources"));
            }
            bases = bases
                .checked_add(molecule.sequence.len())
                .filter(|n| *n <= 100_000)
                .ok_or_else(|| bad("TOPO source aggregate exceeds100000 bases"))?;
        }
        if target.is_some_and(|m| {
            m.sequence.len() > 100_000 || m.sequence.alphabet() != collagene_core::Alphabet::Dna
        }) {
            return Err(bad("TOPO independent target requires bounded DNA"));
        }
        let inputs = sources
            .iter()
            .zip(borrowed)
            .map(|(s, m)| collagene_core::topo::Input {
                role: s.role.clone(),
                molecule: m.clone(),
            })
            .collect::<Vec<_>>();
        collagene_core::topo::source_preflight(&inputs, request).map_err(bad)?;
        Ok((inputs, target.cloned()))
    }
    fn text_artifact(
        &self,
        text: String,
        mime_type: &'static str,
        filename: &'static str,
    ) -> Result<String, ErrorData> {
        let id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let mut store = self.store.lock().unwrap();
        if !store.artifacts.contains_key(&id) {
            if store.artifacts.len() >= MAX_ITEMS {
                return Err(bad("at most 128 artifacts"));
            }
            store.reserve(text.len())?;
            store.artifacts.insert(
                id.clone(),
                Artifact {
                    data: ArtifactData::Text(text),
                    mime_type,
                    encoding: "utf-8",
                    filename,
                },
            );
        }
        Ok(id)
    }
    // Install an immutable codon report and optional candidate document together.
    // Capacity checks precede every mutation, including artifact reservation.
    fn codon_artifact<T: Serialize>(
        &self,
        report: &T,
        candidate: Option<Molecule>,
    ) -> Result<(String, Option<Value>), ErrorData> {
        let text = bounded_json(report)?;
        let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let candidate = candidate.map(|m| (m.revision(), m));
        let mut store = self.store.lock().unwrap();
        let new_artifact = !store.artifacts.contains_key(&artifact_id);
        if new_artifact && store.artifacts.len() >= MAX_ITEMS {
            return Err(bad("at most 128 artifacts"));
        }
        let mut bytes = if new_artifact { text.len() } else { 0 };
        let candidate_summary = if let Some((id, m)) = &candidate {
            m.validate().map_err(bad)?;
            if !store.documents.contains_key(id) {
                if store.documents.len() >= MAX_DOCUMENT_REVISIONS {
                    return Err(bad("at most 10000 session document revisions"));
                }
                bytes = bytes
                    .checked_add(bounded_json(m)?.len())
                    .ok_or_else(|| bad("session capacity"))?;
            }
            Some(summary(id, m))
        } else {
            None
        };
        store.reserve(bytes)?;
        if new_artifact {
            store.artifacts.insert(
                artifact_id.clone(),
                Artifact {
                    data: ArtifactData::Text(text),
                    mime_type: "application/json",
                    encoding: "utf-8",
                    filename: "codon-design.json",
                },
            );
        }
        if let Some((id, m)) = candidate {
            store.documents.entry(id).or_insert(m);
        }
        Ok((artifact_id, candidate_summary))
    }
    // Retain a fresh assembly report and every source-derived consensus together.
    // Normal capacity errors cannot leave an orphan document or report behind.
    fn contig_artifacts<T: Serialize>(
        &self,
        report: &T,
        candidates: Vec<Molecule>,
        batch_revision: &str,
    ) -> Result<(String, Vec<Value>), ErrorData> {
        let text = bounded_json(report)?;
        let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let mut rows = Vec::with_capacity(candidates.len());
        let mut unique = BTreeMap::new();
        for molecule in candidates {
            molecule.validate().map_err(bad)?;
            let encoded_size = bounded_json(&molecule)?.len();
            let id = molecule.revision();
            rows.push(summary(&id, &molecule));
            unique.entry(id).or_insert((molecule, encoded_size));
        }
        let mut store = self.store.lock().unwrap();
        let new_artifact = match store.artifacts.get(&artifact_id) {
            Some(existing) => {
                if !matches!(existing.data, ArtifactData::Text(_))
                    || existing.mime_type != "application/json"
                    || existing.filename != "contig-assembly.json"
                {
                    return Err(bad(
                        "assembly artifact already has a different representation",
                    ));
                }
                false
            }
            None => true,
        };
        if new_artifact && store.artifacts.len() >= MAX_ITEMS {
            return Err(bad("at most 128 artifacts"));
        }
        let new_documents = unique
            .keys()
            .filter(|id| !store.documents.contains_key(*id))
            .count();
        if store.documents.len().saturating_add(new_documents) > MAX_DOCUMENT_REVISIONS {
            return Err(bad("at most 10000 session document revisions"));
        }
        let new_cache = !store.contig_reports.contains_key(batch_revision);
        let mut bytes = if new_artifact { text.len() } else { 0 };
        if new_cache {
            bytes = bytes
                .checked_add(batch_revision.len() + artifact_id.len() + 64)
                .ok_or_else(|| bad("session capacity"))?;
        }
        for (id, (_, size)) in &unique {
            if !store.documents.contains_key(id) {
                bytes = bytes
                    .checked_add(*size)
                    .ok_or_else(|| bad("session capacity"))?;
            }
        }
        store.reserve(bytes)?;
        if new_artifact {
            store.artifacts.insert(
                artifact_id.clone(),
                Artifact {
                    data: ArtifactData::Text(text),
                    mime_type: "application/json",
                    encoding: "utf-8",
                    filename: "contig-assembly.json",
                },
            );
        }
        for (id, (molecule, _)) in unique {
            store.documents.entry(id).or_insert(molecule);
        }
        store
            .contig_reports
            .insert(batch_revision.to_owned(), artifact_id.clone());
        Ok((artifact_id, rows))
    }
    // Candidate donor installation is atomic with its source-bound computed report.
    fn edit_template_artifact<T: Serialize>(
        &self,
        report: &T,
        donor: Option<Molecule>,
    ) -> Result<(String, Option<Value>), ErrorData> {
        self.candidate_document_artifact(report, donor, "edit-template.json")
    }
    fn candidate_document_artifact<T: Serialize>(
        &self,
        report: &T,
        donor: Option<Molecule>,
        filename: &'static str,
    ) -> Result<(String, Option<Value>), ErrorData> {
        self.candidate_document_artifact_bounded(report, donor, filename, MAX_STORE_BYTES)
    }
    fn candidate_document_artifact_bounded<T: Serialize>(
        &self,
        report: &T,
        donor: Option<Molecule>,
        filename: &'static str,
        max_report_bytes: usize,
    ) -> Result<(String, Option<Value>), ErrorData> {
        let text = bounded_json_limit(report, max_report_bytes)?;
        let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let donor = donor
            .map(|m| {
                m.validate().map_err(bad)?;
                let bytes = bounded_json(&m)?.len();
                Ok::<_, ErrorData>((m.revision(), m, bytes))
            })
            .transpose()?;
        let mut store = self.store.lock().unwrap();
        let new_artifact = match store.artifacts.get(&artifact_id) {
            Some(a) => {
                if !matches!(a.data, ArtifactData::Text(_))
                    || a.mime_type != "application/json"
                    || a.filename != filename
                {
                    return Err(bad("candidate artifact already has another representation"));
                }
                false
            }
            None => true,
        };
        if new_artifact && store.artifacts.len() >= MAX_ITEMS {
            return Err(bad("at most 128 artifacts"));
        }
        let mut bytes = if new_artifact { text.len() } else { 0 };
        let donor_summary = if let Some((id, m, size)) = &donor {
            if !store.documents.contains_key(id) {
                if store.documents.len() >= MAX_DOCUMENT_REVISIONS {
                    return Err(bad("at most 10000 session document revisions"));
                }
                bytes = bytes
                    .checked_add(*size)
                    .ok_or_else(|| bad("session capacity"))?;
            }
            Some(summary(id, m))
        } else {
            None
        };
        store.reserve(bytes)?;
        if new_artifact {
            store.artifacts.insert(
                artifact_id.clone(),
                Artifact {
                    data: ArtifactData::Text(text),
                    mime_type: "application/json",
                    encoding: "utf-8",
                    filename,
                },
            );
        }
        if let Some((id, m, _)) = donor {
            store.documents.entry(id).or_insert(m);
        }
        Ok((artifact_id, donor_summary))
    }
    // Full dot-plot JSON and SVG are installed together after all capacity checks.
    fn dotplot_artifacts(
        &self,
        report: &collagene_core::dotplot::DotPlotReport,
        svg: String,
    ) -> Result<(String, String), ErrorData> {
        let text = bounded_json(report)?;
        let entries = [
            (text, "application/json", "dotplot.json"),
            (svg, "image/svg+xml", "dotplot.svg"),
        ];
        self.paired_text_artifacts(entries)
    }
    // Publish the documentary receipt, binary candidate and unverified session
    // plan atomically. Neither workspace heads nor historical status are adopted.
    fn saved_review_artifacts(
        &self,
        review: collagene_core::run_review::PreparedReview,
        reference: WorkspaceRef,
    ) -> Result<CallToolResult, ErrorData> {
        let text = bounded_json_limit(&review, collagene_core::run_review::MAX_REVIEW_BYTES)?;
        let bytes = transfer::save_candidate(&review.candidate).map_err(bad)?;
        let report_id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let candidate_id = blake3::hash(&bytes).to_hex().to_string();
        let plan_revision = review.candidate.plan_revision.clone();
        let plan_bytes =
            bounded_json_limit(&(&review.candidate.plan, &review.target), MAX_STORE_BYTES)?.len();
        let inputs = review.inputs.iter().take(100).collect::<Vec<_>>();
        let assessment = &review.recorded_assessment;
        let recorded_assessment = json!({"engine_version":assessment.engine_version,"status":assessment.status,"findings":assessment.findings.iter().take(20).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),"finding_count":assessment.findings.len(),"findings_truncated":assessment.findings.len()>20||assessment.findings.iter().take(20).any(|s|s.chars().count()>1024),"assumptions":assessment.assumptions.iter().take(20).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),"assumption_count":assessment.assumptions.len(),"assumptions_truncated":assessment.assumptions.len()>20||assessment.assumptions.iter().take(20).any(|s|s.chars().count()>1024),"exhaustive_alternative_search":assessment.exhaustive_alternative_search,"authority":"historical only"});
        let result = output(
            json!({"workspace":reference,"run_id":review.run_id,"context_revision":review.context_revision,"source_workspace_fingerprint":review.source_workspace_fingerprint,"run_attachment_hash":review.run_attachment_hash,"state":"candidate","verified":false,"plan_revision":plan_revision,"target_event":review.target_event,"target":{"id":review.target.revision(),"revision":review.target.revision()},"input_count":review.inputs.len(),"inputs":inputs,"inputs_truncated":review.inputs.len()>100,"recorded_assessment":recorded_assessment,"recorded_engine_matches":review.recorded_engine_matches,"current_verification_engine":review.current_verification_engine,"has_saved_material_origins":review.has_saved_material_origins,"warnings":review.warnings,"artifact_id":report_id,"candidate_artifact_id":candidate_id,"candidate_filename":"saved-run-review.cg2plan","candidate_encoding":"base64","authority":"Documentary preparation of exact retained run inputs and independent target. The session candidate is unverified; call verify_plan for fresh checks. No historical status/material authority was adopted and no workspace head or native catalog was changed."}),
        )?;
        let mut store = self.store.lock().unwrap();
        let new_plan = !store.plans.contains_key(&plan_revision);
        if let Some(existing) = store.plans.get(&plan_revision) {
            if existing.target != review.target
                || bounded_json(&existing.plan)? != bounded_json(&review.candidate.plan)?
            {
                return Err(bad(
                    "saved plan revision already has different exact source/target data",
                ));
            }
        }
        let representations = [
            (
                &report_id,
                false,
                "application/json",
                "saved-run-review.json",
            ),
            (
                &candidate_id,
                true,
                "application/octet-stream",
                "saved-run-review.cg2plan",
            ),
        ];
        let mut added = 0usize;
        let mut charge = if new_plan { plan_bytes } else { 0 };
        for (id, binary, mime, filename) in representations {
            if let Some(existing) = store.artifacts.get(id) {
                if matches!(&existing.data, ArtifactData::Binary(_)) != binary
                    || existing.mime_type != mime
                    || existing.filename != filename
                {
                    return Err(bad("saved review artifact has another representation"));
                }
            } else {
                added += 1;
                charge = charge
                    .checked_add(if binary { bytes.len() } else { text.len() })
                    .ok_or_else(|| bad("session capacity"))?;
            }
        }
        if (new_plan && store.plans.len() >= MAX_ITEMS)
            || store.artifacts.len().saturating_add(added) > MAX_ITEMS
        {
            return Err(bad(
                "session item limit; no saved-review artifacts or plan were published",
            ));
        }
        store.reserve(charge)?;
        store.artifacts.entry(report_id).or_insert(Artifact {
            data: ArtifactData::Text(text),
            mime_type: "application/json",
            encoding: "utf-8",
            filename: "saved-run-review.json",
        });
        store.artifacts.entry(candidate_id).or_insert(Artifact {
            data: ArtifactData::Binary(bytes),
            mime_type: "application/octet-stream",
            encoding: "base64",
            filename: "saved-run-review.cg2plan",
        });
        if new_plan {
            store.plans.insert(
                plan_revision,
                StoredPlan {
                    plan: review.candidate.plan,
                    target: review.target,
                },
            );
        }
        Ok(result)
    }
    fn paired_text_artifacts(
        &self,
        entries: [(String, &'static str, &'static str); 2],
    ) -> Result<(String, String), ErrorData> {
        let ids = entries
            .iter()
            .map(|(text, _, _)| blake3::hash(text.as_bytes()).to_hex().to_string())
            .collect::<Vec<_>>();
        let mut store = self.store.lock().unwrap();
        let mut bytes = 0usize;
        let mut added = 0usize;
        for ((text, mime_type, filename), id) in entries.iter().zip(&ids) {
            if let Some(existing) = store.artifacts.get(id) {
                if !matches!(existing.data, ArtifactData::Text(_))
                    || existing.mime_type != *mime_type
                    || existing.filename != *filename
                {
                    return Err(bad(
                        "artifact content already exists with a different representation",
                    ));
                }
            } else {
                bytes = bytes
                    .checked_add(text.len())
                    .ok_or_else(|| bad("session capacity"))?;
                added += 1;
            }
        }
        if store.artifacts.len().saturating_add(added) > MAX_ITEMS {
            return Err(bad("at most 128 artifacts"));
        }
        store.reserve(bytes)?;
        for ((text, mime_type, filename), id) in entries.into_iter().zip(&ids) {
            store.artifacts.entry(id.clone()).or_insert(Artifact {
                data: ArtifactData::Text(text),
                mime_type,
                encoding: "utf-8",
                filename,
            });
        }
        Ok((ids[0].clone(), ids[1].clone()))
    }
    fn analysis_artifact(&self, value: &Value) -> Result<String, ErrorData> {
        let text = serde_json::to_string(value).map_err(bad)?;
        let id = blake3::hash(text.as_bytes()).to_hex().to_string();
        let mut store = self.store.lock().unwrap();
        if !store.artifacts.contains_key(&id) {
            if store.artifacts.len() >= MAX_ITEMS {
                return Err(bad("at most 128 artifacts"));
            }
            store.reserve(text.len())?;
            store.artifacts.insert(
                id.clone(),
                Artifact {
                    data: ArtifactData::Text(text),
                    mime_type: "application/json",
                    encoding: "utf-8",
                    filename: "sequence-analysis.json",
                },
            );
        }
        Ok(id)
    }
    fn binary_artifact(
        &self,
        bytes: Vec<u8>,
        mime_type: &'static str,
        filename: &'static str,
    ) -> Result<String, ErrorData> {
        // Existing clients use binary artifact IDs as the decoded-content digest.
        let id = blake3::hash(&bytes).to_hex().to_string();
        let mut store = self.store.lock().unwrap();
        if let Some(existing) = store.artifacts.get(&id) {
            if !matches!(existing.data, ArtifactData::Binary(_))
                || existing.mime_type != mime_type
                || existing.filename != filename
            {
                return Err(bad(
                    "artifact content already exists with a different representation",
                ));
            }
        }
        if !store.artifacts.contains_key(&id) {
            if store.artifacts.len() >= MAX_ITEMS {
                return Err(bad("at most 128 artifacts"));
            }
            store.reserve(bytes.len())?;
            store.artifacts.insert(
                id.clone(),
                Artifact {
                    data: ArtifactData::Binary(bytes),
                    mime_type,
                    encoding: "base64",
                    filename,
                },
            );
        }
        Ok(id)
    }
    fn insert_alignment(
        &self,
        document: alignment_formats::AlignmentDocument,
    ) -> Result<CallToolResult, ErrorData> {
        let bytes = alignment_formats::save_alignment(&document).map_err(bad)?;
        let id = blake3::hash(&bytes).to_hex().to_string();
        let mut result = alignment_summary(&id, &document);
        result["source_is_current"] = json!(document.source_is_current().map_err(bad)?);
        let mut store = self.store.lock().unwrap();
        if !store.alignments.contains_key(&id) {
            if store.alignments.len() >= MAX_ITEMS {
                return Err(bad("at most128 alignments"));
            }
            store.reserve(bytes.len())?;
            store.alignments.insert(id, Arc::new(document));
        }
        output(result)
    }
    fn insert_feature_library(
        &self,
        library: feature_library::FeatureLibrary,
    ) -> Result<CallToolResult, ErrorData> {
        let bytes = feature_library::save_library(&library).map_err(bad)?;
        let id = library.revision();
        let mut store = self.store.lock().unwrap();
        let result = json!({"library":{"id":id,"revision":id},"name":library.name,"version":library.version,"entry_count":library.entries.len(),"entries":library.entries.iter().take(100).map(|entry|json!({"id":entry.id,"revision":entry.revision,"name":entry.name.chars().take(256).collect::<String>(),"length":entry.sequence.len(),"kind":entry.kind,"source":{"document_revision":entry.source.document_revision,"feature_index":entry.source.feature_index,"document_name":entry.source.document_name.chars().take(256).collect::<String>(),"source_needs_review":entry.source.source_needs_review}})).collect::<Vec<_>>(),"entries_truncated":library.entries.len()>100,"persistence":"session only"});
        if !store.feature_libraries.contains_key(&id) {
            if store.feature_libraries.len() >= MAX_ITEMS {
                return Err(bad("at most 128 feature libraries"));
            }
            store.reserve(bytes.len())?;
            store.feature_libraries.insert(id, library);
        }
        output(result)
    }
    /// Attest a documentary saved-source capsule against one exact retained
    /// workspace without synthesizing a plan or running scientific reconstruction.
    fn attest_run_source(
        &self,
        reference: &WorkspaceRef,
        origin: &collagene_core::run_material::Capsule,
        source: &Molecule,
    ) -> Result<(), ErrorData> {
        let frame = self.restore_workspace(reference)?;
        let recorded =
            collagene_core::run_lineage::load_run(&frame, &origin.selection.run_id).map_err(bad)?;
        let candidate = transfer::load_candidate(&origin.candidate).map_err(bad)?;
        if candidate.plan_revision != recorded.plan_revision
            || origin.target_event != recorded.target
            || origin.inputs != recorded.inputs
            || origin.products != recorded.products
            || frame.revisions.get(&recorded.target.revision) != Some(&origin.target)
            || frame.revisions.get(&origin.selection.event.revision) != Some(source)
        {
            return Err(bad(
                "saved-run capsule differs from the exact retained run/source/target events",
            ));
        }
        let origins =
            collagene_core::run_material::originating_products(&frame, &origin.selection.event)
                .map_err(bad)?;
        if !origins.contains(&origin.selection) {
            return Err(bad(
                "saved-run origin selection is not retained for this exact source event",
            ));
        }
        Ok(())
    }
    fn resolve_material_derivation_origin(
        &self,
        workspace: &workspace::Workspace,
        event: &collagene_core::run_lineage::EventBinding,
        origin: &MaterialDerivationOrigin,
    ) -> Result<collagene_core::material_derivation::SourceOrigin, ErrorData> {
        use collagene_core::material_derivation::SourceOrigin;
        match origin {
            MaterialDerivationOrigin::OriginalReference => {
                SourceOrigin::OriginalReference
                    .validate_shallow_for(event)
                    .map_err(bad)?;
                Ok(SourceOrigin::OriginalReference)
            }
            MaterialDerivationOrigin::Saved(selection) => self
                .resolve_material_origin(workspace, event, selection)
                .map(SourceOrigin::Saved),
        }
    }
    fn resolve_material_origin(
        &self,
        workspace: &workspace::Workspace,
        event: &collagene_core::run_lineage::EventBinding,
        selection: &collagene_core::run_material::MaterialSelection,
    ) -> Result<collagene_core::run_material::MaterialOrigin, ErrorData> {
        use collagene_core::run_material;

        selection.validate().map_err(bad)?;
        if selection.event() != event {
            return Err(bad(
                "Selected material handle names a different exact source event.",
            ));
        }
        let requested = Some(selection.clone());
        let expected = run_material::material_inspection_revision(workspace, event, &requested)
            .map_err(bad)?;
        let inspection =
            run_material::inspect_material(workspace, event, &requested).map_err(bad)?;
        inspection.validate_binding(&expected).map_err(bad)?;
        let selected=inspection.selected.ok_or_else(||bad("Selected retained material cannot currently be reconstructed; inspect its origin findings."))?;
        if selected.selection() != *selection {
            return Err(bad("Fresh origin differs from explicit selected handle."));
        }
        selected
            .validate_for_workspace(workspace, event)
            .map_err(bad)?;
        Ok(selected)
    }
    fn attest_material_origins(
        &self,
        inputs: &[PlanInput],
        bindings: &[MaterialOriginWorkspace],
        origins: &[(
            &str,
            &str,
            collagene_core::run_material::MaterialOriginRef<'_>,
        )],
    ) -> Result<(), ErrorData> {
        use collagene_core::run_material::{MaterialOrigin, MaterialOriginRef, MaterialSelection};
        if bindings.len() != origins.len() || bindings.len() > 128 {
            return Err(bad(
                "Every saved reaction role requires exactly one retained workspace binding.",
            ));
        }
        let mut seen = BTreeSet::new();
        for binding in bindings {
            if binding.role.is_empty()
                || binding.role.len() > 128
                || !seen.insert(binding.role.as_str())
                || !origins.iter().any(|(role, _, _)| *role == binding.role)
            {
                return Err(bad(
                    "Invalid, duplicate or unrelated saved reaction workspace role.",
                ));
            }
        }
        for (role, input_id, origin) in origins {
            let binding = bindings
                .iter()
                .find(|b| b.role == *role)
                .ok_or_else(|| bad("Saved reaction workspace is unavailable."))?;
            let selection = match origin {
                MaterialOriginRef::Run(c) => MaterialSelection::Run(c.selection.clone()),
                MaterialOriginRef::Derived(c) => MaterialSelection::Derived(c.product.clone()),
            };
            let frame = self.restore_workspace(&binding.workspace)?;
            let actual = self.resolve_material_origin(&frame, selection.event(), &selection)?;
            let source = inputs
                .iter()
                .find(|i| i.id == *input_id)
                .ok_or_else(|| bad("Saved reaction source role is unavailable."))?;
            if frame.revisions.get(&selection.event().revision) != Some(&source.molecule) {
                return Err(bad(
                    "Saved reaction source differs from its exact retained event.",
                ));
            }
            let same = match (&actual, origin) {
                (MaterialOrigin::Run(a), MaterialOriginRef::Run(b)) => {
                    bounded_json_limit(a, 32 * 1024 * 1024)?
                        == bounded_json_limit(b, 32 * 1024 * 1024)?
                }
                (MaterialOrigin::Derived(a), MaterialOriginRef::Derived(b)) => {
                    bounded_json_limit(a, 32 * 1024 * 1024)?
                        == bounded_json_limit(b, 32 * 1024 * 1024)?
                }
                _ => false,
            };
            if !same {
                return Err(bad(
                    "Supplied reaction origin differs from the freshly resolved retained native origin.",
                ));
            }
        }
        Ok(())
    }
    fn material_derivation_source(
        &self,
        workspace: &workspace::Workspace,
        event: &collagene_core::run_lineage::EventBinding,
        origin: &collagene_core::material_derivation::SourceOrigin,
    ) -> Result<Molecule, ErrorData> {
        origin.validate_shallow_for(event).map_err(bad)?;
        workspace.validate().map_err(bad)?;
        let exact = workspace
            .history
            .iter()
            .find(|e| e.id == event.event_id)
            .ok_or_else(|| bad("material source event is absent"))?;
        if exact.document_id != event.document_id || exact.revision != event.revision {
            return Err(bad("material source event/document/revision disagree"));
        }
        let document = workspace
            .documents
            .iter()
            .find(|d| d.id == event.document_id)
            .ok_or_else(|| bad("material source document is absent"))?;
        if document.head != event.revision
            || document.head_event.as_deref() != Some(event.event_id.as_str())
        {
            return Err(bad(
                "material source must be the exact retained current document head",
            ));
        }
        if let collagene_core::material_derivation::SourceOrigin::Saved(selected) = origin {
            selected
                .validate_for_workspace(workspace, event)
                .map_err(bad)?;
        }
        let source = workspace
            .revisions
            .get(&event.revision)
            .ok_or_else(|| bad("material source revision is absent"))?;
        collagene_core::model::bounds::model_preflight(&[source]).map_err(bad)?;
        Ok(source.clone())
    }
    fn restore_workspace(
        &self,
        reference: &WorkspaceRef,
    ) -> Result<workspace::Workspace, ErrorData> {
        let store = self.store.lock().unwrap();
        if reference.id != reference.revision {
            return Err(bad("stale workspace revision"));
        }
        store
            .workspaces
            .get(&reference.id)
            .ok_or_else(|| bad("unknown workspace"))?
            .restore(&store)
    }
    /// Trusted native embedding only: import already selected archive bytes into
    /// this memory-only MCP session. This grants no path or database authority.
    pub fn preload_archive_bytes(&self, bytes: &[u8]) -> Result<CallToolResult, ErrorData> {
        self.import_archive_bytes(bytes)
    }
    /// Load a literal native .cg2 document without inventing workspace events,
    /// run receipts, or reusable physical-material authority.
    pub fn preload_document_bytes(&self, bytes: &[u8]) -> Result<CallToolResult, ErrorData> {
        let molecule = collagene_core::protocol::load_document(bytes).map_err(bad)?;
        self.insert_documents(vec![molecule], vec![])
    }
    fn import_archive_bytes(&self, bytes: &[u8]) -> Result<CallToolResult, ErrorData> {
        let value = workspace::load_workspace(bytes).map_err(bad)?;
        self.commit_workspace(value)
    }
    fn commit_workspace(
        &self,
        workspace: workspace::Workspace,
    ) -> Result<CallToolResult, ErrorData> {
        self.commit_workspace_artifact(workspace, None)
    }
    fn commit_workspace_artifact(
        &self,
        workspace: workspace::Workspace,
        run_artifact: Option<(String, Value)>,
    ) -> Result<CallToolResult, ErrorData> {
        self.commit_workspace_outputs(workspace, run_artifact, None)
    }
    // Stage every workspace/revision, run receipt and optional native output
    // before reserving capacity or publishing any session entry.
    fn commit_workspace_outputs(
        &self,
        workspace: workspace::Workspace,
        run_artifact: Option<(String, Value)>,
        native_artifact: Option<Vec<u8>>,
    ) -> Result<CallToolResult, ErrorData> {
        let finalized = native_artifact.is_some();
        let receipt_filename = if finalized {
            "cloning-outputs.json"
        } else {
            "cloning-run.json"
        };
        self.commit_workspace_outputs_named(
            workspace,
            run_artifact,
            native_artifact,
            receipt_filename,
            "historical run data; replay recomputes current evidence",
        )
    }
    fn commit_workspace_outputs_named(
        &self,
        workspace: workspace::Workspace,
        run_artifact: Option<(String, Value)>,
        native_artifact: Option<Vec<u8>>,
        receipt_filename: &'static str,
        receipt_authority: &'static str,
    ) -> Result<CallToolResult, ErrorData> {
        let mut pending_alignments = BTreeMap::new();
        for attachment in workspace
            .attachments
            .values()
            .filter(|a| a.kind == alignment_formats::WORKSPACE_ATTACHMENT_KIND)
        {
            let document = alignment_formats::load_alignment(&attachment.bytes).map_err(bad)?;
            let bytes = alignment_formats::save_alignment(&document).map_err(bad)?;
            let id = blake3::hash(&bytes).to_hex().to_string();
            pending_alignments.insert(id, (Arc::new(document), bytes.len()));
        }
        let (structure, pending) = StoredWorkspace::split(workspace);
        let metadata = serde_json::to_vec(&structure).map_err(bad)?;
        let id = blake3::hash(&metadata).to_hex().to_string();
        let mut store = self.store.lock().unwrap();
        let result = structure.summary(&id);
        let artifact_id = run_artifact
            .as_ref()
            .map(|(text, _)| blake3::hash(text.as_bytes()).to_hex().to_string());
        let native_id = native_artifact
            .as_ref()
            .map(|bytes| blake3::hash(bytes).to_hex().to_string());
        let mut added_artifacts = 0usize;
        let mut artifact_bytes = 0usize;
        if let Some((text, _)) = &run_artifact {
            let key = artifact_id
                .as_ref()
                .ok_or_else(|| bad("missing staged receipt ID"))?;
            if let Some(existing) = store.artifacts.get(key) {
                if !matches!(existing.data, ArtifactData::Text(_))
                    || existing.mime_type != "application/json"
                    || existing.filename != receipt_filename
                {
                    return Err(bad(
                        "receipt artifact already has another representation; no workspace outputs published",
                    ));
                }
            } else {
                added_artifacts += 1;
                artifact_bytes = text.len();
            }
        }
        if let Some(bytes) = &native_artifact {
            let key = native_id
                .as_ref()
                .ok_or_else(|| bad("missing staged native archive ID"))?;
            if let Some(existing) = store.artifacts.get(key) {
                if !matches!(existing.data, ArtifactData::Binary(_))
                    || existing.mime_type != "application/octet-stream"
                    || existing.filename != "workspace.cg2workspace"
                {
                    return Err(bad(
                        "native workspace artifact already has another representation; no outputs published",
                    ));
                }
            } else {
                added_artifacts += 1;
                artifact_bytes = artifact_bytes
                    .checked_add(bytes.len())
                    .ok_or_else(|| bad("artifact accounting overflow"))?;
            }
        }
        if store.artifacts.len().saturating_add(added_artifacts) > MAX_ITEMS {
            return Err(bad(
                "final workspace artifacts exceed item capacity; no outputs published",
            ));
        }
        if store.workspaces.contains_key(&id) && run_artifact.is_none() && native_artifact.is_none()
        {
            return output(
                json!({"workspace":result,"deduplicated":true,"persistence":"session only"}),
            );
        }
        let deduplicated = store.workspaces.contains_key(&id);
        if !deduplicated && store.workspaces.len() >= MAX_ITEMS {
            return Err(bad("at most 128 workspaces"));
        }
        let (count, bytes) = pending
            .iter()
            .filter(|(r, _)| !store.documents.contains_key(*r))
            .try_fold(
                (0usize, if deduplicated { 0 } else { metadata.len() }),
                |(count, bytes), (_, m)| {
                    Ok::<_, ErrorData>((
                        count + 1,
                        bytes
                            .checked_add(serde_json::to_vec(m).map_err(bad)?.len())
                            .ok_or_else(|| bad("session accounting overflow"))?,
                    ))
                },
            )?;
        if store.documents.len() + count > MAX_DOCUMENT_REVISIONS {
            return Err(bad(
                "at most 10000 session document revisions; workspace import is atomic",
            ));
        }
        let (alignment_count, alignment_bytes) = pending_alignments
            .iter()
            .filter(|(id, _)| !store.alignments.contains_key(*id))
            .try_fold((0usize, 0usize), |(count, bytes), (_, (_, size))| {
                Ok::<_, ErrorData>((
                    count + 1,
                    bytes
                        .checked_add(*size)
                        .ok_or_else(|| bad("alignment accounting overflow"))?,
                ))
            })?;
        if store.alignments.len() + alignment_count > MAX_ITEMS {
            return Err(bad(
                "workspace import exceeds128 alignment collections; no state installed",
            ));
        }
        let total = bytes
            .checked_add(alignment_bytes)
            .and_then(|n| n.checked_add(artifact_bytes))
            .ok_or_else(|| bad("session accounting overflow"))?;
        store.reserve(total)?;
        store.alignments.extend(
            pending_alignments
                .into_iter()
                .map(|(id, (document, _))| (id, document)),
        );
        store.documents.extend(pending);
        store.workspaces.insert(id, structure);
        let mut response =
            json!({"workspace":result,"deduplicated":deduplicated,"persistence":"session only"});
        if let Some((text, receipt)) = run_artifact {
            let artifact_id = artifact_id.unwrap();
            store
                .artifacts
                .entry(artifact_id.clone())
                .or_insert(Artifact {
                    data: ArtifactData::Text(text),
                    mime_type: "application/json",
                    encoding: "utf-8",
                    filename: receipt_filename,
                });
            response["receipt"] = receipt;
            response["artifact_id"] = json!(artifact_id);
            response["authority"] = json!(receipt_authority);
        }
        if let Some(bytes) = native_artifact {
            let native_id = native_id.ok_or_else(|| bad("missing native artifact ID"))?;
            let length = bytes.len();
            store
                .artifacts
                .entry(native_id.clone())
                .or_insert(Artifact {
                    data: ArtifactData::Binary(bytes),
                    mime_type: "application/octet-stream",
                    encoding: "base64",
                    filename: "workspace.cg2workspace",
                });
            response["native_artifact_id"] = json!(native_id);
            response["native_bytes"] = json!(length);
            response["encoding"] = json!("base64");
            response["filename"] = json!("workspace.cg2workspace");
            response["authority"] = json!(
                "Finalized declared-model outputs in a new immutable session workspace only. All native preparation history is retained; presentation visibility and saved order data do not establish experimental material or inherited verification authority. Replay reconstructs again."
            );
        }
        output(response)
    }
    async fn compute<F>(&self, job: F) -> Result<CallToolResult, ErrorData>
    where
        F: FnOnce() -> Result<CallToolResult, ErrorData> + Send + 'static,
    {
        let permit = self
            .compute_gate
            .clone()
            .acquire_owned()
            .await
            .map_err(|e| bad(e.to_string()))?;
        tokio::task::spawn_blocking(move || {
            let _permit = permit;
            job()
        })
        .await
        .map_err(|e| ErrorData::internal_error(e.to_string(), None))?
    }
}
impl CollageneServer {
    /// Configure managed file-backed session storage. This grants no filesystem
    /// tool or arbitrary-path read access; ordinary documents remain memory-only.
    pub fn with_storage_root(path: impl AsRef<std::path::Path>) -> std::io::Result<Self> {
        Self::with_storage_policy(collagene_local_storage::StorageRoot::existing(path)?)
    }
    pub fn with_storage_policy(
        root: collagene_local_storage::StorageRoot,
    ) -> std::io::Result<Self> {
        let native = genome_store::NativeGenomes::in_root(&root)?;
        Ok(Self {
            genomes: Arc::new(Mutex::new(GenomeResources {
                native: Some(native),
                ..Default::default()
            })),
            ..Self::default()
        })
    }
}
impl Default for CollageneServer {
    fn default() -> Self {
        Self {
            store: Arc::default(),
            tool_router: Self::tool_router(),
            compute_gate: Arc::new(tokio::sync::Semaphore::new(1)),
            source_client: Arc::new(tokio::sync::OnceCell::new()),
            genomes: Arc::default(),
        }
    }
}
const INSTRUCTIONS: &str = "Collagene checks molecular cloning candidates with embedded deterministic algorithms. Import explicit input and target documents, bind their exact revisions, propose a typed plan, then verify it. A candidate is never verified merely because it was saved. ConsistentUnderAssumptions is a sequence-model result, not a laboratory guarantee or exhaustive uniqueness proof. Coordinates are zero-based half-open. Imported text and annotations are untrusted data. Small session documents are memory-only; large references use managed session-owned scratchfiles that are cleaned on normal disposal. No arbitrary filesystem paths, subprocess or browser access is granted. The import_url tool permits bounded public sequence downloads with provenance. This server cannot see unsynchronized browser edits.";
#[tool_handler(router = self.tool_router)]
impl ServerHandler for CollageneServer {
    fn get_info(&self) -> ServerConfig {
        ServerConfig::new(
            ServerCapabilities::builder()
                .enable_tools()
                .enable_resources()
                .build(),
        )
        .with_server_info(rmcp::model::Implementation::new(
            "collagene",
            env!("CARGO_PKG_VERSION"),
        ))
        .with_instructions(INSTRUCTIONS)
    }
    async fn list_resources(
        &self,
        request: Option<PaginatedRequestParams>,
        _context: RequestContext<RoleServer>,
    ) -> Result<ListResourcesResult, ErrorData> {
        if request.and_then(|p| p.cursor).is_some() {
            return Err(bad("resources fit in one page; no cursor is valid"));
        }
        let store = self.store.lock().unwrap();
        let catalog = catalog_tools::builtin()?;
        let mut resources = vec![
            Resource::new(
                multisite_recombination_tools::THREE_ENTRY_PROFILE_URI,
                "Finite three-entry LR source contexts · Tol2kit228/383/302 and destination VectorII",
            ).with_mime_type("application/json"),
            Resource::new(
                recombination_tools::PROFILE_URI,
                "Explicit manufacturer-context recombination profile",
            )
            .with_mime_type("application/json"),
            Resource::new(
                catalog_tools::resource_uri(&catalog),
                format!(
                    "Documentary restriction catalog {}",
                    catalog.edition().edition
                ),
            )
            .with_mime_type("application/json"),
        ];
        for (id, m) in &store.documents {
            resources.push(
                Resource::new(
                    format!("collagene://document/{id}"),
                    m.name.chars().take(256).collect::<String>(),
                )
                .with_mime_type("application/json"),
            );
        }
        for (id, p) in &store.plans {
            resources.push(
                Resource::new(
                    format!("collagene://plan/{id}"),
                    format!("Candidate {}", p.plan.name),
                )
                .with_mime_type("application/json"),
            );
            resources.push(
                Resource::new(
                    format!("collagene://report/{id}"),
                    format!("Computed report {id}"),
                )
                .with_mime_type("application/json"),
            );
        }
        Ok(ListResourcesResult::with_all_items(resources))
    }
    async fn read_resource(
        &self,
        request: ReadResourceRequestParams,
        _context: RequestContext<RoleServer>,
    ) -> Result<ReadResourceResponse, ErrorData> {
        let value = if request.uri == recombination_tools::PROFILE_URI {
            recombination_tools::profile()?
        } else if request.uri == multisite_recombination_tools::THREE_ENTRY_PROFILE_URI {
            multisite_recombination_tools::three_entry_profile()?
        } else if request.uri.starts_with("collagene://enzyme-catalog/") {
            let catalog = catalog_tools::builtin()?;
            if request.uri != catalog_tools::resource_uri(&catalog) {
                return Err(bad(
                    "unknown catalog edition/content hash; use the listed exact versioned resource",
                ));
            }
            catalog_tools::resource(&catalog)
        } else if let Some(id) = request.uri.strip_prefix("collagene://document/") {
            let store = self.store.lock().unwrap();
            let m = store
                .documents
                .get(id)
                .ok_or_else(|| bad("unknown document resource"))?;
            let mut value = summary(id, m);
            value["sequence_access"] = json!(
                "Use read_sequence_region with this ID and revision; zero-based half-open, at most 10,000 bases"
            );
            value
        } else if let Some(id) = request.uri.strip_prefix("collagene://plan/") {
            let store = self.store.lock().unwrap();
            let saved = store
                .plans
                .get(id)
                .ok_or_else(|| bad("unknown plan resource"))?;
            json!({"plan_revision":id,"state":"candidate","verified":false,"name":saved.plan.name,"target_revision":saved.plan.target_revision,"step_count":saved.plan.steps.len(),"inputs":saved.plan.inputs.iter().map(|i|json!({"id":i.id,"revision":i.revision})).collect::<Vec<_>>()})
        } else if let Some(id) = request.uri.strip_prefix("collagene://report/") {
            self.verify_plan(Parameters(PlanRef {
                plan_revision: id.into(),
            }))
            .await?
            .structured_content
            .unwrap()
        } else {
            return Err(bad(
                "unknown resource URI; filesystem and network resources are not supported",
            ));
        };
        Ok(ReadResourceResult::new(vec![ResourceContents::text(
            serde_json::to_string(&value).map_err(bad)?,
            request.uri,
        )])
        .into())
    }
}

#[tool_router]
impl CollageneServer {
    #[tool(
        description = r#"Plan explicit source-bound cloning reaction quantities using exact retained candidate/saved-run plan, independent target and ordered input event pins. Preview returns exact arithmetic plus missing-field diagnostics, not chemistry, stock availability, recommended recipes or experimental readiness. Capture requires reviewed context and atomically returns a NEW immutable workspace, full quantity artifact and archive without changing biological heads/status. List/Open inspect saved plans. Export formats Csv/Html/Pdf are quantity companions; ProtocolHtml/ProtocolPdf combine complete cloning instructions, fresh native molecular consistency and declared quantity tables. Molecular status and arithmetic completeness remain separate; neither establishes experimental success. Ordinary plasmid sharing omits this sidecar.

Exact quantity JSON: Volume={"unit":"nL","value":{"numerator":"20000","denominator":"1"}} is20uL. Rational numerator/denominator are decimal unsigned64-bit STRINGS; denominator positive. Never use float, exponent or NaN. Canonical mass concentration unit is"ng/nL":100ng/uL=value1/10. Canonical molar concentration unit is"pmol/nL":10uM=value1/100. Never interchange dimensional tags. Blank direct volume is{"Direct":{"volume":null}}, not zero. MassDilution has {"recipe":{"schema_version":1,"stock":<mass concentration or null>,"target":<mass concentration or null>,"final_volume":<Volume or null>},"quantization":null}; MolarDilution uses molar concentrations. Explicit quantization={"increment":<Volume>,"mode":"ExactOnly"|"Floor"|"Ceiling"|"NearestTiesUp"}; inspect achieved concentration/deviation. Nullable request fields must be present.

MINIMAL_PREVIEW={"workspace":{"id":"<workspace id>","revision":"<workspace revision>"},"task":{"Preview":{"request":{"name":"Explicit PCR quantities","plan_source":{"SavedRun":{"run_id":"<retained run id>"}},"plan_revision":"<plan revision>","target":{"document_id":"<target document>","event_id":"<original target event>","revision":"<target revision>"},"inputs":[{"input_id":"template","event":{"document_id":"<source document>","event_id":"<original source event>","revision":"<source revision>"}}],"reactions":[{"id":"pcr","step_index":0,"final_volume":{"unit":"nL","value":{"numerator":"20000","denominator":"1"}},"fill":{"ToFinal":{"component_id":"water","label":"Operator-declared water"}},"components":[{"id":"template-volume","label":"User-declared template volume","source":{"Input":"template"},"volume":{"Direct":{"volume":{"unit":"nL","value":{"numerator":"2000","denominator":"1"}}}}}]}]}}}} END_PREVIEW

Replace ALL angle-bracket pins with exact native values, never names or equal-DNA substitutions. This arithmetic example leaves all other PCR reagents unspecified;18uL fill is an operator declaration, NOT a recommended reaction recipe. Arithmetic completeness is not scientific readiness. Use plan_source="Candidate" only for an exactly retained candidate and explicitly declared events. Capture adds expected_context_revision from Preview. Open/Export use quantity_hash from capture; Export format Csv,Html,Pdf,ProtocolHtml,ProtocolPdf. Combined outputs include a fresh current-engine molecular receipt, actual findings/assumptions in the instructions and exact original target/input events. They never reuse historical scientific status as fresh acceptance. List requires offset and limit1..200."#
    )]
    async fn protocol_quantities(
        &self,
        Parameters(args): Parameters<quantity_tools::Query>,
    ) -> Result<CallToolResult, ErrorData> {
        self.reaction_quantities(args).await
    }
    #[tool(
        description = "Plan explicit96/384-well documentary cloning assignments using exact retained batch declarations or released run/product events. Preview never records; Capture requires the reviewed context and creates a NEW immutable workspace with a bounded layout attachment, preserving all biological heads and historical status. List/Open inspect saved layouts. Export emits native CSV/HTML/PDF planning documents with exact identifiers and missing-quantity notices. No automatic placement, inventory, schedules, instrument commands, fabricated amounts or experimental success. Ordinary plasmid publication omits plate layouts; complete plate publication requires a separate explicit closure review."
    )]
    async fn protocol_layout(
        &self,
        Parameters(args): Parameters<ProtocolLayoutQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            use collagene_core::protocol_layout as layout;
            let frame=server.restore_workspace(&args.workspace)?;
            match args.task {
                ProtocolLayoutTask::Preview{request}=>{
                    let expected=layout::preview_revision(&frame,&request).map_err(bad)?;
                    let report=layout::preview(&frame,&request).map_err(bad)?;
                    report.validate_binding(&expected).map_err(bad)?;
                    let text=bounded_json_limit(&report,layout::MAX_BYTES)?;
                    let artifact=server.text_artifact(text,"application/json","protocol-layout.json")?;
                    output(json!({"context_revision":report.context_revision,"assignment_count":report.resolved.len(),"warnings":report.warnings,"artifact_id":artifact,"state":"documentary preview","verified":false,"authority":"Well assignments do not verify reactions or laboratory execution. Quantities remain missing unless explicitly provided by a separate typed planning workflow."}))
                }
                ProtocolLayoutTask::Capture{request,expected_context_revision}=>{
                    let expected=layout::preview_revision(&frame,&request).map_err(bad)?;
                    if expected!=expected_context_revision{return Err(bad("Layout source or assignments changed; preview again before capture."));}
                    let commit=layout::capture(&frame,&request).map_err(bad)?;
                    commit.validate_binding(&expected).map_err(bad)?;
                    let text=bounded_json_limit(&commit.layout,layout::MAX_BYTES)?;
                    let receipt=json!({"attachment_hash":commit.attachment_hash,"context_revision":commit.context_revision,"source_workspace_fingerprint":commit.source_workspace_fingerprint,"assignment_count":commit.layout.resolved.len(),"state":"documentary","verified":false});
                    let bytes=workspace::save_workspace(&commit.workspace).map_err(bad)?;
                    server.commit_workspace_outputs_named(commit.workspace,Some((text,receipt)),Some(bytes),"protocol-layout.json","Documentary well assignments only; no laboratory execution or scientific verification authority.")
                }
                ProtocolLayoutTask::List{offset,limit}=>{
                    if limit==0||limit>200||offset>layout::MAX_LAYOUTS{return Err(bad("Invalid layout page bounds."));}
                    let mut entries=Vec::new();
                    for (hash,asset) in frame.attachments.iter().filter(|(_,a)|a.kind==layout::KIND){
                        if entries.len()>=layout::MAX_LAYOUTS{return Err(bad("Too many layout attachments."));}
                        let report=layout::load(&frame,hash).map_err(bad)?;
                        entries.push(json!({"layout_hash":hash,"name":asset.name,"context_revision":report.context_revision,"plate":report.request.plate,"assignment_count":report.request.assignments.len()}));
                    }
                    let total=entries.len();let page=entries.into_iter().skip(offset).take(limit).collect::<Vec<_>>();
                    output(json!({"workspace":args.workspace,"offset":offset,"total":total,"next_offset":(offset+page.len()<total).then_some(offset+page.len()),"layouts":page,"authority":"Historical documentary layouts only."}))
                }
                ProtocolLayoutTask::Open{layout_hash}=>{
                    let report=layout::load(&frame,&layout_hash).map_err(bad)?;
                    let artifact=server.text_artifact(bounded_json_limit(&report,layout::MAX_BYTES)?,"application/json","protocol-layout.json")?;
                    output(json!({"layout_hash":layout_hash,"context_revision":report.context_revision,"assignment_count":report.resolved.len(),"artifact_id":artifact,"state":"documentary","verified":false}))
                }
                ProtocolLayoutTask::ReviewPublication{layout_hash}=>{
                    let review=collagene_library_sharing::inspect_protocol_layout_closure(&frame,&layout_hash).map_err(bad)?;
                    let text=bounded_json_limit(&review,32*1024*1024)?;
                    let scope_revision=blake3::hash(text.as_bytes()).to_hex().to_string();
                    let artifact=server.text_artifact(text,"application/json","protocol-layout-publication.json")?;
                    output(json!({"layout_hash":layout_hash,"scope_revision":scope_revision,"documents":review.inventory.documents.len(),"events":review.inventory.events.len(),"attachments":review.inventory.attachments.len(),"omitted":review.inventory.omitted,"notices":review.notices,"artifact_id":artifact,"authority":"Explicit layout scope review; no archive or publication has been made."}))
                }
                ProtocolLayoutTask::ProjectPublication{layout_hash,expected_scope_revision}=>{
                    let review=collagene_library_sharing::inspect_protocol_layout_closure(&frame,&layout_hash).map_err(bad)?;
                    let text=bounded_json_limit(&review,32*1024*1024)?;
                    if blake3::hash(text.as_bytes()).to_hex().as_str()!=expected_scope_revision{return Err(bad("Layout publication scope changed; review the full dependency counts again."));}
                    let projected=collagene_library_sharing::project_protocol_layout_closure(&frame,&layout_hash).map_err(bad)?;
                    let bytes=workspace::save_workspace(&projected.workspace).map_err(bad)?;
                    server.commit_workspace_outputs_named(projected.workspace,Some((text,json!({"layout_hash":layout_hash,"scope_revision":expected_scope_revision,"authority":"Explicit complete layout projection, including indivisible batch/release dependencies. No server upload or experimental authority."}))),Some(bytes),"protocol-layout-publication.json","Reviewed explicit layout dependency projection only; no server upload, experiment or scientific verification.")
                }
                ProtocolLayoutTask::Export{layout_hash,format}=>{
                    let expected=layout::export_revision(&frame,&layout_hash,format).map_err(bad)?;
                    let exported=layout::export(&frame,&layout_hash,format).map_err(bad)?;
                    exported.validate_binding(&expected).map_err(bad)?;
                    let (mime,filename)=match format{layout::ExportFormat::Csv=>("text/csv","protocol-layout.csv"),layout::ExportFormat::Html=>("text/html","protocol-layout.html"),layout::ExportFormat::Pdf=>("application/pdf","protocol-layout.pdf")};
                    let length=exported.bytes.len();
                    let artifact=server.binary_artifact(exported.bytes,mime,filename)?;
                    output(json!({"layout_hash":layout_hash,"context_revision":expected,"artifact_id":artifact,"filename":filename,"mime_type":mime,"encoding":"base64","bytes":length,"warnings":exported.warnings,"authority":"Documentary planning export, not experimental execution evidence."}))
                }
            }
        }).await
    }
    #[tool(
        description = "Freshly reconstruct every explicitly selected saved combinatorial row against its independent target and atomically release named native constructs, retained preparation/run events and optional oligo collections in a NEW immutable workspace. Per-row output options are explicit. A late failure discards the entire private prefix and publishes no workspace/artifact. Returns full historical release manifest plus native archive; imported reports never authorize release, and experimental success is not asserted."
    )]
    async fn release_combinatorial_outputs(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ReleaseQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_release(args).await
    }
    #[tool(
        description = "Export a saved native combinatorial release as a bounded ZIP of selected GenBank/FASTA, declared-primer CSV/TSV and HTML/PDF protocols. Optional construct IDs select a subset. Stable ordinal/hash entry paths never use names as paths. Protocol files require fresh independent-declaration and native-product checks; sequence/primer-only exchange remains clearly historical. All requested files or an explicit failure, with complete format losses and exact source/run/output events in manifest.json. No subprocess, arbitrary filesystem path or experimental-success claim."
    )]
    async fn export_combinatorial_outputs(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_archive(args).await
    }
    #[tool(
        description = "List retained native combinatorial release manifests with bounded paging. The records identify actual native output events and historical lineage; listing does not recompute science or restore current authority."
    )]
    async fn list_combinatorial_releases(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ListQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_release_list(args).await
    }
    #[tool(
        description = "Open one native combinatorial release manifest as historical JSON, including exact input/target/output events, per-row options and retained run identities. Leaves heads unchanged and never treats imported manifest status as fresh scientific authority."
    )]
    async fn open_combinatorial_release(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ReleaseRef>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_release_open(args).await
    }
    #[tool(
        description = "Page a declared alternative-source matrix in an imported native workspace, with exact retained source-event ownership and original interval checks. Uses stable bin/choice identities and bounded lazy enumeration (at most64 rows), not reaction-product enumeration. Does not test preparation, junctions, independent target agreement or experimental success."
    )]
    async fn page_combinatorial_design(
        &self,
        Parameters(args): Parameters<combinatorial_tools::MatrixQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_matrix(args).await
    }
    #[tool(
        description = "Resolve one explicitly selected construct's exact retained source events, original contributions, per-construct overrides, independent target declaration and typed native preparation recipe. Returns complete source-bound JSON evidence. This inspection does not reconstruct the causal assembly or establish verification; use fresh batch evaluation for that."
    )]
    async fn inspect_combinatorial_construct(
        &self,
        Parameters(args): Parameters<combinatorial_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_inspect(args).await
    }
    #[tool(
        description = "List immutable assembly-library edit recovery records in an imported native workspace, with bounded paging and stable historical hashes/allocation counters. Declared source pins are retained native event references, not complete host-draft attestation or scientific evidence. Does not apply recovery, rewind identities or change source heads."
    )]
    async fn list_assembly_library_history(
        &self,
        Parameters(args): Parameters<assembly_library_history_tools::ListQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.assembly_library_history_list(args).await
    }
    #[tool(
        description = "Open an exact immutable assembly-library recovery hash. Returns complete metadata JSON and the exact raw before-draft JSON in separate artifacts, preserving unsaved text without integer-array expansion. Historical recovery only: no apply/restore action, allocation rewind, complete host-draft validation or scientific authority is granted."
    )]
    async fn open_assembly_library_history(
        &self,
        Parameters(args): Parameters<assembly_library_history_tools::OpenQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.assembly_library_history_open(args).await
    }
    #[tool(
        description = "Discover the finite native Golden Gate reaction pool for one source-bound combinatorial construct, including whole sources, complete PCR amplicons, discarded fragments, carryover and explicit preparation substrates. Returns exact source/target event bindings and source_revision, candidate indices and both sequence keys; complete native report is retained in the JSON artifact. Missing selection is allowed and never defaults to candidate zero. Independent target correspondence is advisory, never a selection rule. Explicit selection followed by fresh batch evaluation and release is required; imported reports confer no authority."
    )]
    async fn inspect_combinatorial_pool(
        &self,
        Parameters(args): Parameters<combinatorial_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_pool_inspect(args).await
    }
    #[tool(
        description = "Inspect whole-parent multi-site BP/LR reference networks for one exact-source library construct. Resolve its independently declared target or freeze original contribution regions/junctions before reaction reconstruction. Returns a source_revision for explicit component selection, bounded candidate summaries and complete typed batch evidence. Missing selection never means candidate zero; matching the target never selects a pairing. Reference-only scope remains incomplete, with no physical recombination or experimental success claim."
    )]
    async fn inspect_combinatorial_multisite(
        &self,
        Parameters(args): Parameters<combinatorial_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_multisite_inspect(args).await
    }
    #[tool(
        description = "Freshly evaluate explicitly selected constructs with the existing deterministic preparation compilers/verifier and retain a named immutable batch record in a NEW native workspace. Stable design/construct/bin/choice identities and exact original events are preserved. Per-construct preparation/method/topology overrides are explicit; unresolved draft rows remain errors. Imported reports never restore fresh authority. Original source heads remain unchanged; native archive and complete historical JSON are returned atomically. This does not release bulk outputs or assert experimental success."
    )]
    async fn evaluate_combinatorial_batch(
        &self,
        Parameters(args): Parameters<combinatorial_tools::BatchQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_evaluate(args).await
    }
    #[tool(
        description = "List retained combinatorial batch records in an imported native workspace with bounded paging. Recorded statuses are historical evidence; listing never recomputes preparations or restores current verification authority."
    )]
    async fn list_combinatorial_batches(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ListQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_list(args).await
    }
    #[tool(
        description = "Open one immutable combinatorial batch record as complete historical JSON. Exact source/target bindings, original choices, per-construct settings and recorded outcomes remain inspectable without changing current heads. Opening does not rerun science or establish fresh authority."
    )]
    async fn open_combinatorial_batch(
        &self,
        Parameters(args): Parameters<combinatorial_tools::OpenQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_open(args).await
    }
    #[tool(
        description = "Freshly replay one explicitly selected row from a saved combinatorial batch using its exact original source events, independent target and preparation settings. Returns recorded and current row evidence separately. Leaves the workspace, batch record and source heads unchanged; current computational consistency does not establish experimental success or release bulk products."
    )]
    async fn replay_combinatorial_batch(
        &self,
        Parameters(args): Parameters<combinatorial_tools::ReplayQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.combinatorial_replay(args).await
    }
    #[tool(
        description = "Search or exactly look up the embedded versioned restriction catalog. Returns edition/hash, stable enzyme geometry identities, separate formulation records, aliases and documentary source citations. Search/exact ambiguity pages are bounded to 100 records. Optional edition/hash pin rejects stale editions. Enzyme identity and formulation assay facts are distinct; lookup does not predict cleavage or transfer supplier conditions between related enzymes."
    )]
    fn lookup_enzyme_catalog(
        &self,
        Parameters(args): Parameters<catalog_tools::LookupQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        output(catalog_tools::lookup(args)?)
    }
    #[tool(
        description = "Resolve one explicit formulation/condition/assay context from the pinned embedded documentary catalog. Keeps Known, Unknown, Conflict and ContextRequired distinct, with every original qualifier, citation, source URL and context preserved. No evidence or conflicting evidence never becomes compatibility or activity. Buffer activity is a contextual vendor assay observation, not digestion efficiency; methylation categories are not site-specific permission. Optional edition/hash pin rejects stale catalog claims."
    )]
    fn resolve_enzyme_evidence(
        &self,
        Parameters(args): Parameters<catalog_tools::EvidenceQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        output(catalog_tools::evidence(args)?)
    }
    #[tool(
        description = "Calculate exact single-stock dilution from decimal-string quantities and explicit mass or molar units. Missing values and mixed dimensions are errors. Optional explicit stock-volume quantization reports absolute rounding deviation, refill diluent and achieved concentration; omitted quantization never rounds. Assumes C1*V1=C2*V2, additive volumes and zero-concentration diluent. No stock inventory lookup, solvent compatibility or reaction verification."
    )]
    fn calculate_dilution(
        &self,
        Parameters(args): Parameters<unit_tools::DilutionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        output(unit_tools::calculate_dilution(args)?)
    }
    #[tool(
        description = "Balance up to 128 uniquely named declared component volumes against an explicit final volume, returning exact specified and remaining fill-diluent volumes. Uses decimal strings and nL/uL/mL unit enums; missing component volumes, overfill and bounded arithmetic overflow remain errors. Additive volume arithmetic only, without concentration, inventory or chemical-compatibility inference."
    )]
    fn balance_mixture(
        &self,
        Parameters(args): Parameters<unit_tools::MixtureQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        output(unit_tools::balance_mixture(args)?)
    }
    #[tool(
        description = "Calculate duplex-DNA molecule-ratio masses, amounts and optional stock volumes for one weighed anchor and up to 127 explicitly named targets. Requires a caller-declared positive g/mol/base-pair mass coefficient, model identity and assumptions; no implicit 660 coefficient or recommended ratio. Exact decimal-string inputs and rational-string outputs, with explicit missing-stock/invalid/overflow errors. Amounts count molecules, not ends; no sequence chemistry, stock inventory or reaction success is verified."
    )]
    fn calculate_dna_molar_ratios(
        &self,
        Parameters(args): Parameters<unit_tools::DnaRatioQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        output(unit_tools::calculate_dna_molar_ratios(args)?)
    }
    #[tool(
        description = "Begin a session-owned large-reference import from bounded FASTA or native .cgref bytes. Pin the exact declared input byte length and optional BLAKE3; build/options apply only to FASTA. Returns an opaque import handle, never a filesystem path. Data are streamed into managed scratchfiles outside the small-document store; no giant Molecule or whole-file RAM buffer. Finish requires full private semantic readback before publication."
    )]
    async fn begin_genome_import(
        &self,
        Parameters(args): Parameters<GenomeBeginQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();let id=g.native()?.begin(args.format,args.expected_bytes,args.expected_hash,args.build,args.options).map_err(bad)?;output(json!({"import_id":id,"expected_bytes":args.expected_bytes,"received":0,"max_chunk_bytes":1048576,"ready":false,"authority":"Unpublished private streaming import; labels or supplied hashes do not establish genome completeness."}))}).await
    }
    #[tool(
        description = "Append one exact-offset base64 chunk (at most 1 MiB decoded) to a live owned large-reference import. Offset/size/base64 checks precede writes. Invalid stream/parser/hash evidence aborts the private import; no partial dataset is published. No caller-selected filesystem paths."
    )]
    async fn append_genome_import(
        &self,
        Parameters(args): Parameters<GenomeAppendQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > 1048576usize.div_ceil(3) * 4 {
            return Err(bad("genome chunk exceeds 1 MiB decoded"));
        }
        let bytes = STANDARD.decode(&args.data_base64).map_err(bad)?;
        if bytes.is_empty() || bytes.len() > 1048576 {
            return Err(bad("genome chunk requires 1..1048576 bytes"));
        }
        let server = self.clone();
        self.compute(move || {
            let mut g = server.genomes.lock().unwrap();
            let received = g
                .existing()?
                .append(&args.import_id, args.offset, &bytes)
                .map_err(bad)?;
            output(json!({"import_id":args.import_id,"received":received,"ready":false}))
        })
        .await
    }
    #[tool(
        description = "Finalize framing and cooperatively verify at most max_blocks (1..64) canonical packed blocks per call. Repeat while done=false. Every block and full canonical contig hash must pass a private readback cursor before the dataset handle is published. Serialized receipts/progress cannot restore authority. No one-shot whole-genome RAM allocation or unbounded verification job."
    )]
    async fn finish_genome_import(
        &self,
        Parameters(args): Parameters<GenomeFinishQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_blocks == 0 || args.max_blocks > 64 {
            return Err(bad("finish max_blocks requires 1..64"));
        }
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();let p=g.existing()?.finish(&args.import_id,args.max_blocks).map_err(bad)?;output(json!({"done":p.done,"import_id":p.import_id,"source_bytes":p.source_bytes,"blocks_verified":p.blocks_verified,"blocks_total":p.blocks_total,"dataset":p.dataset.as_ref().map(native_summary),"authority":"Private complete readback establishes stored DNA integrity only; declared build labels do not attest a complete public reference genome."}))}).await
    }
    #[tool(
        description = "List live large-reference dataset summaries or inspect a bounded contig metadata page for an exact opaque dataset/manifest/content binding. Build/accession/expected-component declarations remain visible and are not independently verified whole-genome coverage. Returns no source DNA or arbitrary paths."
    )]
    async fn inspect_genome(
        &self,
        Parameters(args): Parameters<GenomeInspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_contigs == 0 || args.max_contigs > 100 {
            return Err(bad("contig pages require 1..100 rows"));
        }
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();if let Some(r)=args.genome{r.validate().map_err(bad)?;let d=g.existing()?.inspect(&r.id,&r.manifest_hash,&r.content_revision).map_err(bad)?;let m=&d.reference.manifest;if args.contig_start>m.contigs.len(){return Err(bad("contig page offset exceeds count"));}let end=args.contig_start.saturating_add(args.max_contigs).min(m.contigs.len());let rows=m.contigs[args.contig_start..end].iter().map(|c|json!({"id":c.id,"source_id":c.source_id,"accession":c.accession,"description":c.description.chars().take(256).collect::<String>(),"length":c.length,"topology":c.topology,"sequence_hash":c.sequence_hash,"block_count":c.blocks.len()})).collect::<Vec<_>>();output(json!({"dataset":native_summary(&d),"contig_start":args.contig_start,"next_contig_start":end,"done":end==m.contigs.len(),"contigs":rows,"available_profiles":[collagene_core::crispr::spcas9(),collagene_core::crispr::cas12a()]}))}else{let rows=g.native.as_ref().map(|n|n.list().iter().map(native_summary).collect::<Vec<_>>()).unwrap_or_default();output(json!({"datasets":rows,"storage":"Managed session scratchfiles, separate from the 64 MiB small-document budget."}))}}).await
    }
    #[tool(
        description = "Read a bounded raw .cgref export byte chunk (at most 1 MiB) from an exact live owned dataset. Concatenate decoded chunks in offset order; no whole-file/base64 artifact is materialized. The native framing and hash manifest are retained. Reimport validates all encoded/canonical sequence hashes before publication; no caller-chosen path is accepted."
    )]
    async fn read_genome_export(
        &self,
        Parameters(args): Parameters<GenomeExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_bytes == 0 || args.max_bytes > 1048576 {
            return Err(bad("export chunks require 1..1048576 bytes"));
        }
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();let bytes=g.existing()?.read_export(&args.genome.id,&args.genome.manifest_hash,&args.genome.content_revision,args.offset,args.max_bytes).map_err(bad)?;let next=args.offset.checked_add(bytes.len()as u64).ok_or_else(||bad("export offset overflow"))?;output(json!({"genome":args.genome,"offset":args.offset,"next_offset":next,"done":bytes.len()<args.max_bytes,"data_base64":STANDARD.encode(&bytes),"encoding":"base64","mime_type":"application/vnd.collagene.cgref","bytes":bytes.len()}))}).await
    }
    #[tool(
        description = "Read at most 16384 verified local reference bases under exact global contig coordinates, role and dataset hashes, without converting a large genome into Molecule. Required packed blocks are read and verified locally. This interval does not establish fresh whole-contig or whole-manifest canonical hash coverage."
    )]
    async fn read_genome_interval(
        &self,
        Parameters(args): Parameters<GenomeReadQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args
            .options
            .interval
            .end
            .checked_sub(args.options.interval.start)
            .is_none_or(|n| n > 16384)
        {
            return Err(bad("MCP sequence windows are limited to 16384 bases"));
        }
        let server = self.clone();
        self.compute(move || {
            let mut g = server.genomes.lock().unwrap();
            let n = g.existing()?;
            let d = n
                .inspect(
                    &args.genome.id,
                    &args.genome.manifest_hash,
                    &args.genome.content_revision,
                )
                .map_err(bad)?;
            let req = collagene_core::genome_crispr::read_interval_required_blocks(
                &d.reference.manifest,
                &args.options,
            )
            .map_err(bad)?;
            let blocks = n
                .read_blocks(
                    &args.genome.id,
                    &args.genome.manifest_hash,
                    &args.genome.content_revision,
                    &req,
                )
                .map_err(bad)?;
            let report = collagene_core::genome_crispr::read_interval(
                &d.reference.manifest,
                &args.options,
                &blocks,
            )
            .map_err(bad)?;
            report
                .validate_bindings_for(&d.reference.manifest, &args.options)
                .map_err(bad)?;
            output(json!({"genome":args.genome,"read":report}))
        })
        .await
    }
    #[tool(
        description = "Design a bounded page of source-bound CRISPR guide/PAM candidates over global interval owner coordinates, including verified right halo/circular seam context. The exact source role, manifest/content hashes and nuclease/filter options remain bound. Returns full immutable JSON evidence plus at most100 candidates per page. Local integrity is not whole-contig or build completeness, editing activity or biological safety."
    )]
    async fn design_genome_interval_guides(
        &self,
        Parameters(args): Parameters<GenomeGuidesQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_candidates == 0 || args.max_candidates > 100 {
            return Err(bad("guide pages require1..100 rows"));
        }
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();let n=g.existing()?;let d=n.inspect(&args.genome.id,&args.genome.manifest_hash,&args.genome.content_revision).map_err(bad)?;let req=collagene_core::genome_crispr::interval_required_blocks(&d.reference.manifest,&args.options).map_err(bad)?;let blocks=n.read_blocks(&args.genome.id,&args.genome.manifest_hash,&args.genome.content_revision,&req).map_err(bad)?;let report=collagene_core::genome_crispr::design_interval(&d.reference.manifest,&args.options,&blocks).map_err(bad)?;report.validate_bindings_for(&d.reference.manifest,&args.options).map_err(bad)?;if args.candidate_start>report.candidates.len(){return Err(bad("candidate offset exceeds count"));}let end=args.candidate_start.saturating_add(args.max_candidates).min(report.candidates.len());drop(g);let artifact_id=server.text_artifact(bounded_json(&report)?,"application/json","genome-guides.json")?;output(json!({"genome":args.genome,"manifest_hash":report.manifest_hash,"content_revision":report.content_revision,"engine_version":report.engine_version,"options":report.options,"options_revision":report.options_revision,"complete":report.complete,"comparisons":report.comparisons,"limits":report.limits,"warnings":report.warnings,"candidate_count":report.candidates.len(),"candidate_start":args.candidate_start,"next_candidate_start":end,"done":end==report.candidates.len(),"candidates":&report.candidates[args.candidate_start..end],"artifact_id":artifact_id}))}).await
    }
    #[tool(
        description = "Start a private cooperative whole-manifest CRISPR Hamming/PAM search against one exact live file-backed dataset and explicit reference role. Each next owner block plus verified right halo is processed only by advance. All canonical contig hashes must be verified before exhaustive-for-manifest/absence evidence; missing declared components prevent complete-build flags. Labels/checkpoints cannot grant authority. No external aligner, physical activity or whole-public-build attestation."
    )]
    async fn start_genome_off_target_search(
        &self,
        Parameters(args): Parameters<GenomeSearchQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{let mut g=server.genomes.lock().unwrap();if g.searches.len()>=4{return Err(bad("at most four live genome searches"));}let d=g.existing()?.inspect(&args.genome.id,&args.genome.manifest_hash,&args.genome.content_revision).map_err(bad)?;let state=collagene_core::genome_crispr::GenomeSearchSession::begin(&d.reference.manifest,&args.query,&args.options).map_err(bad)?;g.serial=g.serial.checked_add(1).ok_or_else(||bad("search serial overflow"))?;let id=blake3::hash(format!("{}:{}:{}",args.genome.id,args.genome.manifest_hash,g.serial).as_bytes()).to_hex().to_string();let progress=state.progress();g.searches.insert(id.clone(),NativeSearch{genome:args.genome.clone(),state});output(json!({"job_id":id,"genome":args.genome,"progress":progress,"authority":"Private live cursor only; serialized progress cannot resume or prove fresh coverage."}))}).await
    }
    #[tool(
        description = "Advance a live whole-manifest search by at most max_blocks1..8 owner steps. The backend reads exact required owner/halo blocks from its managed file; clients cannot submit a forged continuation or skip components. Each call is bounded and cooperative. Coverage and source-integrity flags remain incomplete until private ordered canonical-hash replay reaches EOF. Errors invalidate the private search instead of returning false absence."
    )]
    async fn advance_genome_off_target_search(
        &self,
        Parameters(args): Parameters<GenomeAdvanceQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_blocks == 0 || args.max_blocks > 8 {
            return Err(bad("search advance requires1..8 owner blocks"));
        }
        let server = self.clone();
        self.compute(move || {
            let mut g = server.genomes.lock().unwrap();
            let mut search = g
                .searches
                .remove(&args.job_id)
                .ok_or_else(|| bad("unknown live genome search"))?;
            for _ in 0..args.max_blocks {
                let p = search.state.progress();
                if p.finished {
                    break;
                }
                let blocks = g
                    .existing()?
                    .read_blocks(
                        &search.genome.id,
                        &search.genome.manifest_hash,
                        &search.genome.content_revision,
                        &p.next_blocks,
                    )
                    .map_err(bad)?;
                search.state.advance(&blocks).map_err(bad)?;
            }
            let progress = search.state.progress();
            let genome = search.genome.clone();
            g.searches.insert(args.job_id.clone(), search);
            output(json!({"job_id":args.job_id,"genome":genome,"progress":progress}))
        })
        .await
    }
    #[tool(
        description = "Read a bounded exact hit page from a live privately computed whole-manifest search. Every hit retains global contig coordinates, orientation, ambiguity bounds and exact intended-role/manifest/coordinate classification. Partial/cancelled/missing-component evidence cannot imply whole-genome absence. Imported JSON or checkpoint fields cannot populate this private state."
    )]
    async fn read_genome_off_target_hits(
        &self,
        Parameters(args): Parameters<GenomePageQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.limit == 0 || args.limit > 100 {
            return Err(bad("hit pages require1..100 rows"));
        }
        let server = self.clone();
        self.compute(move||{let g=server.genomes.lock().unwrap();let s=g.searches.get(&args.job_id).ok_or_else(||bad("unknown live genome search"))?;let hits=s.state.page(args.offset,args.limit).map_err(bad)?;let progress=s.state.progress();let next=args.offset.checked_add(hits.len()).ok_or_else(||bad("hit offset overflow"))?;output(json!({"job_id":args.job_id,"genome":s.genome,"progress":progress,"offset":args.offset,"next_offset":next,"done":next==s.state.progress().hits_retained,"hits":hits}))}).await
    }
    #[tool(
        description = "Dispose a session-owned genome import/search handle, or release a dataset when its exact genome binding is supplied. Releasing a dataset cancels related searches and deletes its managed scratchfile; arbitrary paths are never accepted. Private parser/search cursors cannot be restored from serialized progress. No persistent application workspace changes."
    )]
    async fn cancel_genome_job(
        &self,
        Parameters(args): Parameters<GenomeCancelQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let mut g = server.genomes.lock().unwrap();
            if let Some(r) = args.genome {
                if r.id != args.id {
                    return Err(bad("dataset disposal binding mismatch"));
                }
                g.existing()?
                    .delete_dataset(&r.id, &r.manifest_hash, &r.content_revision)
                    .map_err(bad)?;
                g.searches.retain(|_, s| s.genome.id != r.id);
                return output(json!({"id":args.id,"removed":true,"kind":"dataset"}));
            }
            if let Some(mut search) = g.searches.remove(&args.id) {
                let progress = search.state.cancel();
                return output(
                    json!({"id":args.id,"removed":true,"kind":"search","progress":progress}),
                );
            }
            let removed = g
                .native
                .as_mut()
                .is_some_and(|n| n.cancel_import(&args.id).is_ok());
            output(json!({"id":args.id,"removed":removed,"kind":"import"}))
        })
        .await
    }
    #[tool(
        description = "Freshly reconstruct TOPO strand preparation from explicit immutable source roles, source intervals or exact PCR/tailing, declared end chemistry/activation, pairing and requested transfers. No imported molecule metadata or foreign physical product grants activation/state. Returns bounded component/nick/unpaired/transfer summaries and a full <=32MiB source/request-bound JSON artifact, retaining every component. A declared graph computation is not kit efficacy, kinetics, manufacturing conservation or experimental success."
    )]
    async fn simulate_topo(
        &self,
        Parameters(args): Parameters<TopoSimulateQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (native,_)=server.hydrate_topo(&args.request,args.origin_workspace.as_ref(),&[])?;

            let(inputs,_)=server.topo_sources(&args.sources,&native,None)?;
            let report=collagene_core::topo::simulate(&inputs,&native).map_err(bad)?;
            let evidence=TopoEvidence{origin_workspace:&args.origin_workspace,source_bindings:&args.sources,request:&native,target_requirement:None,report:&report,
                binding_scope:"Exact immutable MCP document roles/revisions; native event parents require the typed plan/run workflow."};
            let artifact=server.text_artifact(bounded_json_limit(&evidence,collagene_core::topo::MAX_REPORT_BYTES)?,"application/json","topo-simulation.json")?;
            output(json!({"sources":args.sources,"engine_version":report.engine_version,"batch_revision":report.batch_revision,
                "source_identities":report.sources,"preparation_assumptions":report.preparation.assumptions,
                "unrepresented_reference_interval_count":report.preparation.not_represented.len(),"physical":topo_summary(&report.physical),
                "artifact_id":artifact,"authority":"Fresh computation under explicit preparation/activation declarations; no inferred kit state, winner or experimental validation."}))
        }).await
    }
    #[tool(
        description = "Inspect both declared TOPO vector/insert orientations from freshly reconstructed sources. Explicit prepared-strand indexes and directional flag; input transfers must be empty. Retains all physical components, failed transfers and alternative orientation graphs in full <=32MiB JSON. Returns exact orientation transfer requests and fresh source batch pins for explicit later projection, without choosing a winner or using target DNA; two orientations are not an exhaustive biochemical mixture search."
    )]
    async fn inspect_topo_orientations(
        &self,
        Parameters(args): Parameters<TopoOrientationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (native,_)=server.hydrate_topo(&args.request,args.origin_workspace.as_ref(),&[])?;

            let(inputs,_)=server.topo_sources(&args.sources,&native,None)?;
            let report=collagene_core::topo::inspect_orientations(&inputs,&native,args.vector,args.insert,args.directional).map_err(bad)?;
            let candidates=report.reference.candidates.iter().zip(&report.candidate_requests).enumerate().map(|(index,(physical,request))|{
                Ok::<_,ErrorData>(json!({"index":index,"batch_revision":collagene_core::topo::batch_revision(&inputs,request).map_err(bad)?,"steps":request.steps,"physical":topo_summary(physical)}))
            }).collect::<Result<Vec<_>,_>>()?;
            let evidence=TopoEvidence{origin_workspace:&args.origin_workspace,source_bindings:&args.sources,request:&native,target_requirement:None,report:&report,
                binding_scope:"Exact immutable MCP roles/revisions; complete candidate requests are retained, never target-selected."};
            let artifact=server.text_artifact(bounded_json_limit(&evidence,collagene_core::topo::MAX_REPORT_BYTES)?,"application/json","topo-orientations.json")?;
            output(json!({"sources":args.sources,"engine_version":report.engine_version,"batch_revision":report.batch_revision,
                "orientation_search_complete":report.reference.orientation_search_complete,"biochemical_mixture_exhaustive":report.reference.biochemical_mixture_exhaustive,
                "candidate_count":candidates.len(),"candidates":candidates,"artifact_id":artifact,
                "authority":"All declared orientation alternatives retained; no preferred winner, target-derived selection or exhaustive mixture claim."}))
        }).await
    }
    #[tool(
        description = "Freshly reconstruct one explicitly selected TOPO component and source-call seed pinned to the source/recipe batch. Projection is withheld unless a supported paired contour can be proven; retains full physical graph, nicks, unpaired DNA, residual components and unapplied transfers. A projected immutable sequence document carries conservative source-coordinate annotations and remains a chemistry-unmapped view, atomically installed with <=32MiB full graph/annotation evidence; it cannot substitute for graph replay or independent target comparison. Imported reports/products are never accepted as physical state."
    )]
    async fn project_topo_product(
        &self,
        Parameters(args): Parameters<TopoProjectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (native,_)=server.hydrate_topo(&args.request,args.origin_workspace.as_ref(),&[])?;

            let(inputs,_)=server.topo_sources(&args.sources,&native,None)?;
            let product=collagene_core::topo::product(&inputs,&native,&args.selection).map_err(bad)?;
            let evidence=TopoEvidence{origin_workspace:&args.origin_workspace,source_bindings:&args.sources,request:&native,target_requirement:None,report:&product,
                binding_scope:"Source-backed physical reconstruction; sequence view is not standalone chemistry authority."};
            let (artifact,document)=server.candidate_document_artifact_bounded(&evidence,product.reference.clone(),"topo-projection.json",collagene_core::topo::MAX_REPORT_BYTES)?;
            output(json!({"sources":args.sources,"selection":args.selection,"engine_version":product.engine_version,
                "batch_revision":product.batch_revision,"projection":topo_projection_summary(&product),
                "physical":topo_summary(&product.projection.physical),"annotations":topo_annotation_summary(product.annotations.as_ref()),"document":document,"artifact_id":artifact,
                "authority":"Explicit component/seed only; no target-guided path selection, nick repair or standalone physical authority from the projected reference."}))
        }).await
    }
    #[tool(
        description = "Freshly reconstruct the explicitly selected TOPO component/seed and compare an independently supplied immutable target under a fully explicit AcceptancePolicy. Target DNA is a noncausal requirement and never chooses components, seeds, orientation or product sequence. Exact equality cannot silently seal nicks, specify unknown terminal chemistry, discard other components or apply failed transfers. Returns conditional accepted/rejected/withheld outcome and complete source/request/target-bound <=32MiB JSON; no experimental success guarantee."
    )]
    async fn check_topo_target(
        &self,
        Parameters(args): Parameters<TopoTargetQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (native,_)=server.hydrate_topo(&args.request,args.origin_workspace.as_ref(),&[])?;

            let(inputs,target)=server.topo_sources(&args.sources,&native,Some(&args.target))?;
            let target=target.ok_or_else(||bad("explicit independent TOPO target required"))?;
            let report=collagene_core::topo::check_target(&inputs,&native,&args.selection,&target,&args.policy).map_err(bad)?;
            let evidence=TopoEvidence{origin_workspace:&args.origin_workspace,source_bindings:&args.sources,request:&native,target_requirement:Some(&args.target),report:&report,
                binding_scope:"Independent target is a noncausal requirement, not an ancestry parent or source of generated DNA."};
            let artifact=server.text_artifact(bounded_json_limit(&evidence,collagene_core::topo::MAX_REPORT_BYTES)?,"application/json","topo-target-check.json")?;
            output(json!({"sources":args.sources,"target":args.target,"selection":args.selection,"policy":args.policy,
                "engine_version":report.engine_version,"batch_revision":report.batch_revision,"target_revision":report.target_revision,
                "outcome":report.outcome,"projection":topo_projection_summary(&report.product),
                "physical":topo_summary(&report.product.projection.physical),"artifact_id":artifact,"annotations":topo_annotation_summary(report.product.annotations.as_ref()),
                "authority":"Conditional comparison under declared physical graph and explicit allowances, not a repaired molecule or experimental verification."}))
        }).await
    }
    #[tool(
        description = "Inspect explicitly profiled BP/LR reciprocal reference-sequence alternatives from exactly two ordered immutable source roles. Mandatory full profile snapshot and resource limits; all reciprocal components and uncertainty retained in immutable JSON. Bounded candidate page only, no automatic winner or product installation, no chemical activity claims."
    )]
    async fn inspect_recombination(
        &self,
        Parameters(args): Parameters<recombination_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.recombination_inspect(args).await
    }
    #[tool(
        description = "Freshly compare independently supplied target DNA against all explicitly profiled BP/LR reference components. Preserves possible sites, unsupported geometry, ambiguity and incomplete-search flags; no-match claims are scoped to the declared full-site profile. The target is noncausal and never supplies product DNA."
    )]
    async fn check_recombination_target(
        &self,
        Parameters(args): Parameters<recombination_tools::CheckQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.recombination_check(args).await
    }
    #[tool(
        description = "Freshly reconstruct one exact candidate/component pinned to its source/profile batch and require a fully specified independent target match. Atomically installs its immutable reference-only document and full checked JSON; all reciprocal component identities retained. No physical reaction, strand state, efficacy or experimental winner is inferred. Native history uses the typed plan/run workflow."
    )]
    async fn select_recombination_component(
        &self,
        Parameters(args): Parameters<recombination_tools::SelectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.recombination_select(args).await
    }
    #[tool(
        description = "Inspect simultaneous multi-parent reciprocal BP/LR reference alternatives from 2..5 exact immutable source roles using a complete caller-supplied finite profile and explicit specificity set. Optional native source pins are attested by owner/event/revision and literal molecule equality. Returns bounded candidate/component summaries and full typed batch/annotation JSON; no curated biological profile, automatic winner or physical activity claim."
    )]
    async fn inspect_multisite_recombination(
        &self,
        Parameters(args): Parameters<multisite_recombination_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.multisite_inspect(args).await
    }
    #[tool(
        description = "Freshly reconstruct one explicitly selected multisite reference component pinned to its exact source/profile/options batch. Atomically retains the projected reference document and full JSON with exact source bindings, reciprocal identities and conservative annotations. No independent target match, repaired state or physical reaction is inferred; this is reference-only material."
    )]
    async fn select_multisite_recombination_component(
        &self,
        Parameters(args): Parameters<multisite_recombination_tools::ProductQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.multisite_select(args).await
    }
    #[tool(
        description = "Freshly compare one explicitly selected multisite reference component against an independent immutable target. Target never determines site matching, component selection or product DNA. Retains full typed evidence, mismatch and incomplete/unsupported search diagnostics; does not install documents or establish chemical success."
    )]
    async fn check_multisite_recombination_target(
        &self,
        Parameters(args): Parameters<multisite_recombination_tools::CheckQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.multisite_check(args).await
    }
    #[tool(
        description = "Design a literal source-derived homology editing donor against an independently supplied immutable target. Explicit edit/deletion/replacement/arm/orientation options and an optional exact source-bound guide are independently checked base by base. The target is a noncausal requirement, never donor source DNA. Returns bounded donor/edited-reference previews and full immutable JSON; installs an immutable donor document only when the independent check is consistent, atomically with the artifact. Sequence consistency is not cellular HDR, nuclease activity or editing success. No physical operation or persistent application writes."
    )]
    async fn design_edit_template(
        &self,
        Parameters(args): Parameters<EditTemplateQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let(source,target)={let store=server.store.lock().unwrap();let source=store.document(&args.source)?;let target=store.document(&args.target)?;
                collagene_core::edit_template::source_preflight(source,target,&args.options,args.guide.as_ref()).map_err(bad)?;(source.clone(),target.clone())};
            let report=collagene_core::edit_template::design(&source,&target,&args.options,args.guide.as_ref()).map_err(bad)?;
            report.validate_bindings_for(&source,&target,&args.options,args.guide.as_ref()).map_err(bad)?;
            let donor=report.check.consistent.then(||report.donor.clone());
            let(artifact_id,donor_document)=server.edit_template_artifact(&report,donor)?;
            output(json!({"source":args.source,"target":args.target,"engine_version":report.engine_version,"source_revision":report.source_revision,
                "target_revision":report.target_revision,"batch_revision":report.batch_revision,"options_revision":report.options_revision,"options":report.options,"guide":report.guide,
                "left_arm_length":report.left.sequence.len(),"right_arm_length":report.right.sequence.len(),"deleted_count":report.deleted_positions.len(),
                "donor_length":report.donor.sequence.len(),"donor_preview":report.donor.sequence.as_str().chars().take(200).collect::<String>(),"donor_preview_truncated":report.donor.sequence.len()>200,
                "edited_reference_length":report.edited_reference.sequence.len(),"edited_reference_preview":report.edited_reference.sequence.as_str().chars().take(200).collect::<String>(),
                "edited_reference_preview_truncated":report.edited_reference.sequence.len()>200,"consistent":report.check.consistent,"check":report.check,"warnings":report.warnings,
                "artifact_id":artifact_id,"donor_document":donor_document,
                "authority":"Literal source/edit/donor/independent-target sequence consistency only; the target is a noncausal requirement and this does not prove cellular HDR, nuclease activity or editing success."}))
        }).await
    }
    #[tool(
        description = "Independently check an immutable supplied donor DNA document against separately bound source/target documents, exact editing options and optional source-bound guide. Borrows/preflights source, target and linear bounded donor before copying DNA. Reconstructs expected arms/edit and checks the candidate base by base; imported metadata and foreign positive reports grant no authority. Returns exact donor/source/target/options hashes and a complete immutable check artifact; no donor search, product edits, repair prediction or persistent writes."
    )]
    async fn check_edit_template(
        &self,
        Parameters(args): Parameters<EditTemplateCheckQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let(source,target,donor)={let store=server.store.lock().unwrap();let source=store.document(&args.source)?;let target=store.document(&args.target)?;
                collagene_core::edit_template::source_preflight(source,target,&args.options,args.guide.as_ref()).map_err(bad)?;
                let donor=store.document(&args.donor)?;
                if donor.sequence.alphabet()!=collagene_core::Alphabet::Dna || donor.topology!=collagene_core::Topology::Linear || donor.sequence.len()>collagene_core::edit_template::MAX_DONOR{return Err(bad("donor candidate must be bounded linear DNA"));}
                (source.clone(),target.clone(),donor.sequence.as_str().to_owned())};
            let report=collagene_core::edit_template::check(&source,&target,&args.options,args.guide.as_ref(),&donor).map_err(bad)?;
            report.validate_bindings_for(&source,&target,&args.options,args.guide.as_ref(),&donor).map_err(bad)?;
            let artifact_id=server.text_artifact(bounded_json(&report)?,"application/json","edit-template-check.json")?;
            output(json!({"source":args.source,"target":args.target,"donor_document":args.donor,"check":report,"artifact_id":artifact_id,
                "authority":"Independent literal source/edit/candidate/target sequence check; no cellular repair, empirical activity, physical product or inherited verification authority."}))
        }).await
    }
    #[tool(
        description = "Design source-bound CRISPR guide/PAM candidates under an explicit canonical SpCas9, Cas12a or custom nuclease geometry and bounded GC/homopolymer/search filters. Default options use SpCas9; returns available builtin profiles, paged exact positions/strands, nominal cuts, ambiguity/filter evidence and complete immutable JSON. Model geometry and nominal cuts do not establish activity, editing outcome or genome-wide specificity. No sequence edits or persistent application writes."
    )]
    async fn design_crispr_guides(
        &self,
        Parameters(args): Parameters<CrisprGuidesQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_candidates == 0 || args.max_candidates > 100 {
            return Err(bad("guide pages require max_candidates 1..100"));
        }
        let server = self.clone();
        self.compute(move||{
            let molecule={let store=server.store.lock().unwrap();let molecule=store.document(&args.document)?;collagene_core::crispr::guide_preflight(molecule,&args.options).map_err(bad)?;molecule.clone()};
            let report=collagene_core::crispr::design_guides(&molecule,&args.options).map_err(bad)?;
            report.validate_bindings_for(&molecule,&args.options).map_err(bad)?;
            if args.candidate_start>report.candidates.len(){return Err(bad("guide page start exceeds retained candidate count"));}
            let end=args.candidate_start.saturating_add(args.max_candidates).min(report.candidates.len());
            let artifact_id=server.text_artifact(bounded_json(&report)?,"application/json","crispr-guides.json")?;
            output(json!({"document":args.document,"engine_version":report.engine_version,"source_revision":report.source_revision,
                "batch_revision":report.batch_revision,"options_revision":report.options_revision,"profile_revision":report.profile_revision,
                "options":report.options,"source_length":report.source_length,"circular":report.circular,
                "complete":report.complete,"comparisons":report.comparisons,"limits":report.limits,"warnings":report.warnings,
                "candidate_count":report.candidates.len(),"candidate_start":args.candidate_start,"next_candidate_start":end,"done":end==report.candidates.len(),
                "candidates":&report.candidates[args.candidate_start..end],"artifact_id":artifact_id,
                "available_profiles":[collagene_core::crispr::spcas9(),collagene_core::crispr::cas12a()],
                "authority":"Source-bound PAM/guide geometry under the declared model; no empirical activity, editing outcome or whole-genome specificity."}))
        }).await
    }
    #[tool(
        description = "Search an explicit immutable selected reference set for bounded CRISPR guide Hamming off-target evidence with exact PAM/strand/coordinate binding. Borrowed reference preflight precedes owned copies. Build/accession/scope labels are user declarations, not independently verified genome coverage. Every retained hit remains visible; is_intended requires exact role, revision, guide/PAM positions and orientation. Returns bounded hit pages and full JSON evidence. No bulge model, genome-wide safety claim, sequence editing or persistent writes."
    )]
    async fn search_crispr_off_targets(
        &self,
        Parameters(args): Parameters<CrisprOffTargetsQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_hits == 0 || args.max_hits > 100 {
            return Err(bad("off-target pages require max_hits 1..100"));
        }
        if args.references.is_empty() || args.references.len() > 64 {
            return Err(bad(
                "off-target search requires 1..64 explicit reference records",
            ));
        }
        let server = self.clone();
        self.compute(move||{
            let references={let store=server.store.lock().unwrap();
                let borrowed=args.references.iter().map(|r|Ok(collagene_core::crispr::ReferenceBorrow{id:&r.id,molecule:store.document(&r.document)?,build:r.build.as_deref(),accession:r.accession.as_deref()})).collect::<Result<Vec<_>,ErrorData>>()?;
                collagene_core::crispr::reference_batch_preflight(&args.query,&borrowed,&args.options).map_err(bad)?;
                borrowed.into_iter().map(|r|collagene_core::crispr::ReferenceInput{id:r.id.to_owned(),molecule:r.molecule.clone(),build:r.build.map(str::to_owned),accession:r.accession.map(str::to_owned)}).collect::<Vec<_>>()};
            let report=collagene_core::crispr::search_off_targets(&args.query,&references,&args.options).map_err(bad)?;
            report.validate_bindings_for(&args.query,&references,&args.options).map_err(bad)?;
            if args.hit_start>report.hits.len(){return Err(bad("off-target page start exceeds retained hit count"));}
            let end=args.hit_start.saturating_add(args.max_hits).min(report.hits.len());
            let artifact_id=server.text_artifact(bounded_json(&report)?,"application/json","crispr-off-targets.json")?;
            output(json!({"reference_documents":args.references,"engine_version":report.engine_version,"batch_revision":report.batch_revision,
                "options_revision":report.options_revision,"profile_revision":report.profile_revision,"query":report.query,"options":report.options,
                "references":report.references,"complete":report.complete,"exhaustive_for_selected_records":report.exhaustive_for_selected_records,
                "comparisons":report.comparisons,"index_entries":report.index_entries,"index_used":report.index_used,"limits":report.limits,"warnings":report.warnings,
                "hit_count":report.hits.len(),"hit_start":args.hit_start,"next_hit_start":end,"done":end==report.hits.len(),"hits":&report.hits[args.hit_start..end],"artifact_id":artifact_id,
                "authority":"Hamming/PAM evidence for the explicit immutable selected reference set only; declared build labels do not prove complete genome coverage, activity or empirical safety."}))
        }).await
    }
    #[tool(
        description = "Assess starter-enzyme recognition sites under explicitly declared Dam, Dcm and CpG modification assumptions using versioned supplier profiles and bounded sequence-context matching. Returns source/options/profile-bound site evidence, explicit ambiguity/unsupported/partial-search warnings, paged site summaries and a complete immutable JSON artifact. Hemi/mixed and unknown states cannot silently become unmodified. This assessment does not modify DNA, predict digest products, prove empirical activity or grant cloning verification authority."
    )]
    async fn assess_methylation(
        &self,
        Parameters(args): Parameters<MethylationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_sites == 0 || args.max_sites > 100 {
            return Err(bad("methylation site pages require max_sites 1..100"));
        }
        let server = self.clone();
        self.compute(move || {
            let molecule = {
                let store = server.store.lock().unwrap();
                let molecule = store.document(&args.document)?;
                collagene_core::methylation::source_preflight(molecule, &args.options).map_err(bad)?;
                molecule.clone()
            };
            let report = collagene_core::methylation::assess(&molecule, &args.options).map_err(bad)?;
            report.validate_bindings_for(&molecule, &args.options).map_err(bad)?;
            if args.site_start > report.sites.len() {
                return Err(bad("methylation page start exceeds retained site count"));
            }
            let end = args.site_start.saturating_add(args.max_sites).min(report.sites.len());
            let text = bounded_json(&report)?;
            let artifact_id = server.text_artifact(text, "application/json", "methylation-assessment.json")?;
            output(json!({
                "document":args.document,"engine_version":report.engine_version,"profile_version":report.profile_version,
                "profile_revision":report.profile_revision,"source_revision":report.source_revision,
                "batch_revision":report.batch_revision,"options_revision":report.options_revision,"options":report.options,
                "source_length":report.source_length,"circular":report.circular,"profiles":report.profiles,
                "complete":report.complete,"comparisons":report.comparisons,"limits":report.limits,"warnings":report.warnings,
                "site_count":report.sites.len(),"site_start":args.site_start,"next_site_start":end,"done":end==report.sites.len(),
                "sites":&report.sites[args.site_start..end],"artifact_id":artifact_id,
                "authority":"Sequence-context assessment under explicit source modification assumptions and declared profiles; no empirical activity, digest product or cloning verification authority."
            }))
        }).await
    }
    #[tool(
        description = "Assemble immutable linear DNA reads without an external reference using embedded affine-gap overlap, conservative layout and source-call consensus models. Explicit trim/raw-score/missing-quality policies, branches, tied overlaps, conflicts, gaps and computational limits remain visible; no invented calibrated probabilities or reference-filled consensus. Returns bounded contig/overlap previews, complete JSON evidence and immutable source-derived consensus document references installed atomically. Repeat the exact batch/options for cached deterministic pages; no imported report authority, empirical uniqueness claim or persistent application writes."
    )]
    async fn assemble_contigs(
        &self,
        Parameters(args): Parameters<AssembleContigsQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_contigs == 0
            || args.max_contigs > 10
            || args.max_overlaps == 0
            || args.max_overlaps > 50
        {
            return Err(bad(
                "assembly pages require max_contigs1..10 and max_overlaps1..50",
            ));
        }
        if args.reads.is_empty() || args.reads.len() > 32 {
            return Err(bad("assembly requires 1..32 read records"));
        }
        args.options.validate().map_err(bad)?;
        let mut ids = std::collections::BTreeSet::new();
        for read in &args.reads {
            if read.id.is_empty()
                || read.id.len() > 256
                || read.id.chars().any(char::is_control)
                || !ids.insert(&read.id)
            {
                return Err(bad(
                    "read role IDs must be distinct nonempty bounded strings",
                ));
            }
        }
        let server = self.clone();
        self.compute(move || {
            let reads = {
                let store = server.store.lock().unwrap();
                // Bound repeated source references before cloning owned worker inputs.
                let mut bytes = 0usize;
                for read in &args.reads {
                    let molecule = store.document(&read.document)?;
                    if molecule.sequence.len() > 10000 || molecule.features.len() > 4096 || molecule.metadata.len() > 4096 {
                        return Err(bad("assembly source exceeds the supported input bounds"));
                    }
                    bytes = bytes.checked_add(bounded_json(molecule)?.len()).ok_or_else(||bad("assembly source size"))?;
                    if bytes > 16 * 1024 * 1024 {
                        return Err(bad("assembly sources exceed 16 MiB"));
                    }
                }
                args.reads.iter().map(|r| Ok(collagene_core::contig::ReadInput {
                    id: r.id.clone(), molecule: store.document(&r.document)?.clone(),
                })).collect::<Result<Vec<_>, ErrorData>>()?
            };
            let batch_revision = collagene_core::contig::batch_revision(&reads, &args.options).map_err(bad)?;
            let cached = {
                let store = server.store.lock().unwrap();
                store.contig_reports.get(&batch_revision)
                    .and_then(|id| store.artifacts.get(id))
                    .and_then(|a| match &a.data { ArtifactData::Text(text) if a.mime_type=="application/json" && a.filename=="contig-assembly.json" => Some(text.clone()), _=>None })
            };
            let reused_result = cached.is_some();
            let report = match cached {
                Some(text) => serde_json::from_str::<collagene_core::contig::ContigReport>(&text).map_err(bad)?,
                None => collagene_core::contig::assemble(&reads, &args.options).map_err(bad)?,
            };
            report.validate_bindings_for(&reads, &args.options).map_err(bad)?;
            if args.contig_start > report.contigs.len() || args.overlap_start > report.overlaps.len() {
                return Err(bad("assembly page start exceeds the retained result count"));
            }
            let end = args.contig_start.saturating_add(args.max_contigs).min(report.contigs.len());
            let overlap_end = args.overlap_start.saturating_add(args.max_overlaps).min(report.overlaps.len());
            let molecules = report.contigs.iter().filter_map(|c| c.molecule.clone()).collect();
            let (artifact_id, documents) = server.contig_artifacts(&report, molecules, &batch_revision)?;
            let warning_preview = |values: &[String]| values.iter().take(8).map(|s|s.chars().take(256).collect::<String>()).collect::<Vec<_>>();
            let contigs = report.contigs[args.contig_start..end].iter().map(|c|json!({
                "id":c.id,"name":c.name.chars().take(128).collect::<String>(),"length":c.sequence.len(),
                "sequence_preview":c.sequence.chars().take(200).collect::<String>(),"sequence_preview_truncated":c.sequence.len()>200,
                "complete":c.complete,"layout_unique":c.layout_unique,"column_count":c.columns.len(),"placement_count":c.placements.len(),
                "consensus_document_revision":c.molecule.as_ref().map(Molecule::revision),"warnings":warning_preview(&c.warnings)
            })).collect::<Vec<_>>();
            let overlaps = report.overlaps[args.overlap_start..overlap_end].iter().map(|e|json!({
                "id":e.id,"from":e.from,"to":e.to,"from_reverse":e.from_reverse,"to_reverse":e.to_reverse,
                "score":e.score,"from_start":e.from_start,"from_end":e.from_end,"to_start":e.to_start,"to_end":e.to_end,
                "paired_bases":e.paired_bases,"certain_matches":e.certain_matches,"error_columns":e.error_columns,"columns":e.columns,
                "optimal_paths_saturated":e.optimal_paths_saturated,"accepted":e.accepted,"selected":e.selected,"kind":e.kind,
                "reason":e.reason.as_ref().map(|s|s.chars().take(256).collect::<String>())
            })).collect::<Vec<_>>();
            let read_preview = report.reads.iter().take(32).map(|r|json!({
                "id":r.id,"name":r.name.chars().take(128).collect::<String>(),"source_revision":r.source_revision,
                "source_length":r.source_length,"source_call_length":r.source_call_length,"source_reverse":r.source_reverse,
                "trim_start":r.trim_start,"trim_end":r.trim_end,"retained_start":r.retained_start,"retained_end":r.retained_end,
                "eligible_calls":r.eligible_calls,"quality_source":r.quality_source,"duplicate_group":r.duplicate_group,
                "error":r.error.as_ref().map(|s|s.chars().take(256).collect::<String>())
            })).collect::<Vec<_>>();
            output(json!({
                "reads":args.reads,"engine_version":report.engine_version,"batch_revision":report.batch_revision,
                "options_revision":report.options_revision,"options":report.options,"complete":report.complete,"cells_used":report.cells_used,
                "reused_result":reused_result,"read_count":report.reads.len(),"read_preview":read_preview,"read_preview_truncated":report.reads.len()>32,
                "contig_count":report.contigs.len(),"contig_start":args.contig_start,"next_contig_start":end,"contigs_done":end==report.contigs.len(),"contigs":contigs,
                "overlap_count":report.overlaps.len(),"overlap_start":args.overlap_start,"next_overlap_start":overlap_end,"overlaps_done":overlap_end==report.overlaps.len(),"overlaps":overlaps,
                "unassembled_count":report.unassembled.len(),"unassembled_preview":report.unassembled.iter().take(32).map(|r|json!({"read_index":r.read_index,"reason":r.reason.chars().take(256).collect::<String>()})).collect::<Vec<_>>(),
                "limits":warning_preview(&report.limits),"warnings":warning_preview(&report.warnings),"artifact_id":artifact_id,"consensus_documents":documents,
                "authority":"Session-computed read evidence under the declared assembly model; not physical independence, empirical sequence identity or imported verification authority."
            }))
        }).await
    }
    #[tool(
        description = "Analyze a revision-bound pair of DNA, RNA or protein sequences as a bounded window dot plot. Returns exact original-reference/query coordinates, forward/reverse orientation, compatible versus certain matching evidence, explicit sampling and work/hit/index completeness, a bounded hit page, and complete JSON plus vector SVG artifacts. No external aligner, inferred homology/ancestry, sequence edits or persistent writes. Repeat with hit_start/max_hits for deterministic bounded inspection pages; artifacts remain content-addressed and session only."
    )]
    async fn analyze_dotplot(
        &self,
        Parameters(args): Parameters<DotPlotQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_hits == 0 || args.max_hits > 200 {
            return Err(bad("dot plot hit pages require max_hits1..200"));
        }
        let server = self.clone();
        self.compute(move || {
            let (reference, query) = { let store=server.store.lock().unwrap(); (store.document(&args.reference)?.clone(), store.document(&args.query)?.clone()) };
            let report = collagene_core::dotplot::analyze(&reference,&query,&args.options).map_err(bad)?;
            report.validate_bindings_for(&reference,&query,&args.options).map_err(bad)?;
            if args.hit_start > report.hits.len() { return Err(bad("dot plot hit_start exceeds retained hit count")); }
            let end = args.hit_start.saturating_add(args.max_hits).min(report.hits.len());
            let svg=collagene_core::dotplot::export_svg(&report).map_err(bad)?;
            let (artifact_id,svg_artifact_id)=server.dotplot_artifacts(&report,svg)?;
            output(json!({"reference":args.reference,"query":args.query,"engine_version":report.engine_version,"reference_revision":report.reference_revision,"query_revision":report.query_revision,"options_revision":report.options_revision,"reference_length":report.reference_length,"query_length":report.query_length,"alphabet":report.alphabet,"reference_circular":report.reference_circular,"query_circular":report.query_circular,"strategy":report.strategy,"work_used":report.work_used,"examined_windows":report.examined_windows,"total_windows":report.total_windows,"grid_complete":report.grid_complete,"complete":report.complete,"sampled":report.sampled,"hit_count":report.hits.len(),"hit_start":args.hit_start,"next_hit_start":end,"done":end==report.hits.len(),"hits":&report.hits[args.hit_start..end],"limits":report.limits,"warnings":report.warnings,"artifact_id":artifact_id,"svg_artifact_id":svg_artifact_id,"authority":"Window-matching evidence only; not an alignment, homology proof or physical ancestry."}))
        }).await
    }
    #[tool(
        description = "Design synonymous DNA for an immutable concrete protein or exact DNA coding source using explicit genetic code, integer codon preferences and synthesis constraints. Deterministic bounded search independently checks every returned candidate against the source protein; incomplete search never claims unsatisfiability or a global optimum. Contexts constrain junctions but are excluded from the candidate CDS. Returns a bounded preview, full JSON artifact and immutable candidate document reference when successful. No organism expression, physical synthesis or cloning verification claim; no persistent application writes."
    )]
    async fn design_codons(
        &self,
        Parameters(args): Parameters<CodonDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let source = server.store.lock().unwrap().document(&args.document)?.clone();
            let report = collagene_core::codon_design::design(&source, &args.options).map_err(bad)?;
            report.validate_for_source(&source, &args.options).map_err(bad)?;
            let (artifact_id, candidate_document) = server.codon_artifact(&report, report.candidate.clone())?;
            output(json!({"document":args.document,"engine_version":report.engine_version,"source_revision":report.source_revision,"options_revision":report.options_revision,"target_protein_preview":report.target_protein.chars().take(200).collect::<String>(),"target_protein_length":report.target_protein.len(),"target_preview_truncated":report.target_protein.len()>200,"status":report.status,"search_complete":report.search_complete,"optimality_proven":report.optimality_proven,"work_used":report.work_used,"pruned_states":report.pruned_states,"score":report.score,"candidate":candidate_document,"candidate_sequence_preview":report.candidate.as_ref().map(|m|m.sequence.as_str().chars().take(200).collect::<String>()),"candidate_preview_truncated":report.candidate.as_ref().is_some_and(|m|m.sequence.len()>200),"candidate_check":report.candidate_check.as_ref().map(|check|json!({"complete":check.complete,"satisfies":check.satisfies,"candidate_revision":check.candidate_revision,"gc_count":check.gc_count,"coding_length":check.coding_length,"score":check.score})),"limits":report.limits.iter().take(100).collect::<Vec<_>>(),"limits_truncated":report.limits.len()>100,"warnings":report.warnings.iter().take(100).collect::<Vec<_>>(),"warnings_truncated":report.warnings.len()>100,"artifact_id":artifact_id,"persistence":"session only","authority":"Candidate sequence constraints are checked against the immutable required protein. No organism expression, physical synthesis or cloning verification is established."}))
        }).await
    }
    #[tool(
        description = "Independently check a supplied immutable DNA candidate against a separately bound concrete protein or exact DNA coding source and explicit synthesis constraints. Re-derives the required protein, translation, codon preference score, GC, forbidden motifs on both strands and supplied context junctions. Search is not used. Incomplete bounded checks never grant a positive pass. Returns bounded evidence and a full JSON artifact; no edits or persistent writes."
    )]
    async fn check_codon_design(
        &self,
        Parameters(args): Parameters<CodonCheckQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let (source, candidate) = { let store = server.store.lock().unwrap();
                (store.document(&args.source)?.clone(), store.document(&args.candidate)?.clone()) };
            let report = collagene_core::codon_design::check(&source, &candidate, &args.options).map_err(bad)?;
            report.validate_for(&source, &candidate, &args.options).map_err(bad)?;
            let (artifact_id, _) = server.codon_artifact(&report, None)?;
            output(json!({"source":args.source,"candidate":args.candidate,"engine_version":report.engine_version,"source_revision":report.source_revision,"candidate_revision":report.candidate_revision,"options_revision":report.options_revision,"target_protein_preview":report.target_protein.chars().take(200).collect::<String>(),"translated_protein_preview":report.translated_protein.chars().take(200).collect::<String>(),"target_protein_length":report.target_protein.len(),"translated_protein_length":report.translated_protein.len(),"protein_previews_truncated":report.target_protein.len()>200||report.translated_protein.len()>200,"complete":report.complete,"satisfies":report.satisfies,"gc_count":report.gc_count,"coding_length":report.coding_length,"score":report.score,"violations":report.violations.iter().take(100).collect::<Vec<_>>(),"violations_truncated":report.violations.len()>100,"work_used":report.work_used,"warnings":report.warnings.iter().take(100).collect::<Vec<_>>(),"warnings_truncated":report.warnings.len()>100,"artifact_id":artifact_id}))
        }).await
    }
    #[tool(
        description = "Analyze an immutable protein document with explicit terminal-stop, ambiguity, disulfide and modification assumptions. Embedded composition, average/monoisotopic mass, charge/pI and extinction models expose unavailable or bounded values rather than inventing exact ambiguity/modification chemistry. Full revision-bound report is a JSON artifact; no sequence edits or experimental claims."
    )]
    async fn analyze_protein(
        &self,
        Parameters(args): Parameters<ProteinPropertiesQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let molecule = server
                .store
                .lock()
                .unwrap()
                .document(&args.document)?
                .clone();
            let report =
                collagene_core::protein_properties::analyze_protein(&molecule, &args.options)
                    .map_err(bad)?;
            let value = serde_json::to_value(&report).map_err(bad)?;
            let artifact = server.analysis_artifact(&value)?;
            let mut summary = value;
            summary["document"] = json!(args.document);
            summary["artifact_id"] = json!(artifact);
            output(summary)
        })
        .await
    }
    #[tool(
        description = "Import complete explicit VCF or BCF/BGZF bytes as base64, <=4 MiB decoded under a fixed bounded MCP parser profile. Publishes a private immutable source-hash collection only after complete parsing. Original lexical bytes, header/sample columns, typed INFO/FORMAT and genotype phasing/missingness remain evidence; no build identity, measured genotype, variant calling or clinical inference. No filesystem paths, reference fetching or external executable."
    )]
    async fn import_variants(
        &self,
        Parameters(args): Parameters<VariantImport>,
    ) -> Result<CallToolResult, ErrorData> {
        const CAP: usize = 4 * 1024 * 1024;
        if args.data_base64.len() > CAP.div_ceil(3) * 4 {
            return Err(bad("variant MCP source exceeds 4 MiB decoded"));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            if bytes.len() > CAP {
                return Err(bad("variant MCP source exceeds 4 MiB decoded"));
            }
            let collection = collagene_core::variant_reference::exchange::Collection::import(
                &bytes,
                args.format.core(),
                &collagene_core::variants_session::mcp_limits(),
            )
            .map_err(bad)?;
            let id = collection.source_blake3().to_owned();
            let response = variant_summary(&collection);
            let mut store = server.store.lock().unwrap();
            if let Some(existing) = store.variants.get(&id) {
                if existing.original_bytes() != collection.original_bytes()
                    || existing.format() != collection.format()
                {
                    return Err(bad("variant collection source identity conflict"));
                }
            } else {
                if store.variants.len() >= MAX_ITEMS {
                    return Err(bad("at most 128 variant collections"));
                }
                store.reserve(collection.retained_charge())?;
                store.variants.insert(id, Arc::new(collection));
            }
            output(response)
        })
        .await
    }
    #[tool(
        description = "Inspect at most 100 records from a private complete VCF/BCF collection. Full bounded normalized VCF row evidence is in a JSON artifact, retaining sample columns/missingness/phasing and symbolic/breakend declarations; stdio previews explicitly truncate long lines. Positions are declared 1-based, absent stays absent. Header/sample/build assertions are not validated DNA, inventory, measured genotype or clinical evidence."
    )]
    async fn inspect_variants(
        &self,
        Parameters(args): Parameters<VariantPage>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.limit) {
            return Err(bad("variant pages must be 1..100"));
        }
        let server = self.clone();
        self.compute(move||{
            let collection=server.store.lock().unwrap().variants(&args.collection)?;
            let rows=collection.page(args.start,args.limit).map_err(bad)?;
            let full=rows.iter().map(|r|json!({"index":r.index,"reference_name":r.reference_name,"position_1based":r.position_1based,"vcf":r.vcf})).collect::<Vec<_>>();
            let evidence=json!({"collection":args.collection,"source_blake3":collection.source_blake3(),"start":args.start,"records":full,"interpretation":"Normalized declared variant rows; exact original lexical bytes remain separately available."});
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","variant-record-page.json")?;
            let mut response=variant_summary(&collection);response["start"]=json!(args.start);response["next_start"]=json!(args.start+rows.len());response["done"]=json!(args.start+rows.len()==collection.len());
            response["records"]=json!(rows.iter().map(|r|json!({"index":r.index,"reference_name":read_name(Some(&r.reference_name)),"reference_name_truncated":r.reference_name.chars().count()>256,"position_1based":r.position_1based,"vcf_preview":r.vcf.chars().take(1024).collect::<String>(),"vcf_preview_truncated":r.vcf.chars().count()>1024})).collect::<Vec<_>>());response["artifact_id"]=json!(artifact);output(response)
        }).await
    }
    #[tool(
        description = "Freshly compare one explicitly selected 1-based ALT and declared REF against actual immutable DNA DocumentRef with exact file contig mapping. Requires private collection source hash and current document revision; wrong contig/stale identity rejected. Returns bounded source/record/reference-bound JSON evidence. Match only establishes REF correspondence, never genotype, build identity, sample status or variant validity; no name-based automatic reference lookup."
    )]
    async fn check_variant_reference(
        &self,
        Parameters(args): Parameters<VariantCheck>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (collection,reference,options)=capture_variant(&server.store.lock().unwrap(),&args)?;
            let report=collagene_core::variant_reference::check_record(&collection,&reference,&options).map_err(bad)?;
            let evidence=json!({"collection":args.collection,"reference":args.reference,"report":report});
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","variant-reference-check.json")?;
            output(json!({"collection":args.collection,"reference":args.reference,"record_index":args.record_index,"alternate_index":args.alternate_index,"revision":report.revision,"outcome":report.outcome,"interpretation":report.interpretation,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Freshly check a literal selected variant REF/ALT against an explicitly mapped immutable DNA reference, then create an edited reference-only preview. Symbolic/SV/breakend/ambiguous or mismatching alleles cannot produce a preview. Atomically installs reference-only document plus bounded source/reference-bound JSON artifact; never physical material, confirmed genotype or target-guided variant selection. Original sources remain unchanged."
    )]
    async fn preview_variant_allele(
        &self,
        Parameters(args): Parameters<VariantCheck>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (collection,reference,options)=capture_variant(&server.store.lock().unwrap(),&args)?;
            let preview=collagene_core::variant_reference::preview(&collection,&reference,&options).map_err(bad)?;
            let evidence=json!({"collection":args.collection,"reference":args.reference,"report":preview.report,"interpretation":"Reference-only explicit selected-allele design; no measured genotype or physical-material authority."});
            let (artifact,document)=server.candidate_document_artifact_bounded(&evidence,Some(preview.molecule),"variant-allele-preview.json",32*1024*1024)?;
            output(json!({"collection":args.collection,"reference":args.reference,"record_index":args.record_index,"alternate_index":args.alternate_index,"revision":preview.report.revision,"outcome":preview.report.outcome,"reference_only":true,"document":document,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Replay a retained variant reference-only preview from its bounded unverified recipe against the original private collection source_id and actual immutable reference/design DocumentRefs. Freshly reconstructs and compares literal DNA plus topology only. Returns full bounded source/reference/design-bound JSON evidence and claimed/fresh revisions; installs nothing and restores no genotype, physical chemistry, biological function or historical verification authority. Designs without retained recipes require explicit original record/ALT selection instead."
    )]
    async fn replay_variant_preview(
        &self,
        Parameters(args): Parameters<VariantReplay>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.source_id.len() != 64
            || !args
                .source_id
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(bad("source_id must be the exact private collection BLAKE3"));
        }
        let server = self.clone();
        self.compute(move||{
            let (collection,reference,design)={
                let store=server.store.lock().unwrap();
                let collection=store.variants(&EvidenceRef{id:args.source_id.clone(),revision:args.source_id.clone()})?;
                let reference=store.document(&args.reference)?;let design=store.document(&args.design)?;
                collagene_core::model::bounds::model_preflight(&[reference,design]).map_err(bad)?;
                (collection,reference.clone(),design.clone())
            };
            let replay=collagene_core::variant_reference::replay_preview(&collection,&reference,&design).map_err(bad)?;
            let fresh_revision=replay.fresh.molecule.revision();
            let evidence=json!({"source_id":args.source_id,"reference":args.reference,"design":args.design,"fresh_revision":fresh_revision,"replay":replay});
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","variant-preview-replay.json")?;
            output(json!({"source_id":args.source_id,"reference":args.reference,"design":args.design,
                "literal_dna_and_topology_match":replay.literal_dna_and_topology_match,"claimed_revision":replay.claimed_revision,"fresh_revision":fresh_revision,"interpretation":replay.interpretation,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Export a private completed variant collection as exact original bytes or normalized VCF/BCF, with explicit conversion notices/errors for unsupported representations. Returns bounded binary artifact read via read_artifact; any codec/size/storage failure publishes no partial artifact. No paths or external bcftools. Original is the byte-exact evidence mode."
    )]
    async fn export_variants(
        &self,
        Parameters(args): Parameters<VariantExport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_bytes == 0 || args.max_bytes > MAX_STORE_BYTES {
            return Err(bad("variant export byte cap must be 1..64 MiB"));
        }
        let server = self.clone();
        self.compute(move||{
            use collagene_core::variant_reference::exchange::{Format,ExportFormat};
            let collection=server.store.lock().unwrap().variants(&args.collection)?;
            let mode=match args.format{VariantExportFormat::Original=>ExportFormat::Original,VariantExportFormat::Vcf=>ExportFormat::Vcf,VariantExportFormat::Bcf=>ExportFormat::Bcf};
            let format=match args.format{VariantExportFormat::Original=>collection.format(),VariantExportFormat::Vcf=>Format::Vcf,VariantExportFormat::Bcf=>Format::Bcf};
            let exported=collection.export(mode,args.max_bytes).map_err(bad)?;let len=exported.bytes.len();
            let artifact=server.binary_artifact(exported.bytes,"application/octet-stream",if format==Format::Vcf{"variants.vcf"}else{"variants.bcf"})?;
            output(json!({"collection":args.collection,"source_blake3":collection.source_blake3(),"bytes":len,"format":if format==Format::Vcf{"vcf"}else{"bcf"},"original":matches!(args.format,VariantExportFormat::Original),"notices":exported.notices,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Dispose one exact private variant collection source-hash handle, releasing its retained session charge. Existing independent output artifacts/reference previews remain unchanged. A serialized collection summary cannot restore parser authority; reimport exact original bytes explicitly."
    )]
    async fn dispose_variants(
        &self,
        Parameters(args): Parameters<VariantQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let mut store = self.store.lock().unwrap();
        let collection = store.variants(&args.collection)?;
        if Arc::strong_count(&collection) > 2 {
            return Err(bad(
                "variant collection is busy; retry disposal after active requests finish",
            ));
        }
        let used = store
            .used
            .checked_sub(collection.retained_charge())
            .ok_or_else(|| bad("variant retained accounting mismatch"))?;
        store.variants.remove(&args.collection.id);
        store.used = used;
        output(json!({"disposed":true,"collection":args.collection}))
    }
    #[tool(
        description = "Import a complete explicit SAM or BAM/BGZF file as base64, up to4MiB decoded under finite source/record/retention limits. Strict EOF/framing/CRC completion publishes an immutable source-hash collection atomically; any parser/resource failure publishes nothing. All records, flags, exact header/dictionary, absent SEQ/QUAL, auxiliary evidence and original bytes remain retained. Parsing is not alignment, reference validation, coverage or variant evidence; no paths, URI/M5 lookup, network fetch or external aligner."
    )]
    async fn import_read_alignments(
        &self,
        Parameters(args): Parameters<ReadAlignmentImport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_source_bytes == 0
            || args.max_source_bytes > MAX_READ_ALIGNMENT_SOURCE
            || args.max_retained_bytes == 0
            || args.max_retained_bytes > MAX_STORE_BYTES
            || args.max_records == 0
            || args.max_records > 100_000
            || args.data_base64.len() > args.max_source_bytes.div_ceil(3) * 4
        {
            return Err(bad(
                "read alignment import exceeds explicit4MiB source/64MiB retention/100000 record ceilings",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            if bytes.len() > args.max_source_bytes {
                return Err(bad("read alignment source byte limit"));
            }
            let mut limits = collagene_core::read_alignment_store::ImportLimits::default();
            limits.max_source_bytes = args.max_source_bytes;
            limits.max_retained_bytes = args.max_retained_bytes;
            limits.max_records = args.max_records;
            limits.reader.max_input_bytes = args.max_source_bytes as u64;
            limits.reader.max_records = args.max_records as u64;
            let collection = collagene_core::read_alignment_store::Collection::import(
                &bytes,
                args.format.core(),
                &limits,
            )
            .map_err(bad)?;
            let id = collection.source_blake3().to_owned();
            let response = read_alignment_summary(&collection);
            let mut store = server.store.lock().unwrap();
            if let Some(existing) = store.read_alignments.get(&id) {
                if existing.original_bytes() != collection.original_bytes() {
                    return Err(bad("read alignment source identity conflict"));
                }
            } else {
                if store.read_alignments.len() >= MAX_ITEMS {
                    return Err(bad("at most128 read alignment collections"));
                }
                store.reserve(collection.retained_charge())?;
                store.read_alignments.insert(id, Arc::new(collection));
            }
            output(response)
        })
        .await
    }
    #[tool(
        description = "Inspect an immutable completed SAM/BAM collection with bounded record/dictionary pages, each at most100. Reports raw flags, reference positions, CIGAR previews, MAPQ/mate fields, absent sequence/quality and exact record hashes; original symbols/tags remain in source export. Sequence/quality cells and huge tag arrays are not returned. Dictionary/build/M5/UR and CIGAR/MD/NM are imported assertions, not actual reference validation or unique coverage. No automatic reference lookup."
    )]
    async fn inspect_read_alignments(
        &self,
        Parameters(args): Parameters<ReadAlignmentPage>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.limit) || !(1..=100).contains(&args.reference_limit) {
            return Err(bad("read alignment record/reference pages must be1..100"));
        }
        let server = self.clone();
        self.compute(move|| {
            let collection=server.store.lock().unwrap().read_alignments(&args.collection)?;
            let records=collection.page(args.start,args.limit).map_err(bad)?;
            let references=collection.header().references();
            if args.reference_start>references.len() {return Err(bad("reference dictionary page offset"));}
            let reference_end=args.reference_start.saturating_add(args.reference_limit).min(references.len());
            let page:Vec<_>=records.iter().enumerate().map(|(offset,record)| {
                let f=record.fields();
                json!({"index":args.start+offset,"record_blake3":record.original_blake3(),
                    "name":read_name(f.name.as_deref()),"name_truncated":f.name.as_ref().is_some_and(|s|s.chars().count()>256),
                    "flags":f.flags,"reverse":f.flags&16!=0,"unmapped":f.flags&4!=0,"secondary":f.flags&256!=0,"supplementary":f.flags&2048!=0,
                    "reference_name":read_name(f.reference.as_deref()),"reference_name_truncated":f.reference.as_ref().is_some_and(|s|s.chars().count()>256),
                    "reported_reference_start":f.start,"mapping_quality":f.mapping_quality,
                    "cigar_op_count":f.cigar.len(),"cigar_truncated":f.cigar.len()>64,
                    "cigar":f.cigar.iter().take(64).map(|op|json!({"length":op.len,"operation":op.kind.symbol().to_string()})).collect::<Vec<_>>(),
                    "sequence_present":f.sequence.is_some(),"sequence_length":f.sequence.as_ref().map(String::len),
                    "qualities_present":f.qualities.is_some(),"qualities_length":f.qualities.as_ref().map(Vec::len),
                    "mate_reference":read_name(f.mate_reference.as_deref()),"mate_start":f.mate_start,"template_length":f.template_length,
                    "aux_tag_count":f.aux.len(),"aux_tags":f.aux.iter().take(32).map(|a|String::from_utf8_lossy(&a.tag).into_owned()).collect::<Vec<_>>(),"aux_tags_truncated":f.aux.len()>32,
                    "virtual_offset":record.virtual_offset().map(|v|json!({"compressed":v.compressed.to_string(),"uncompressed":v.uncompressed})),
                })
            }).collect();
            let mut response=read_alignment_summary(&collection);
            response["start"]=json!(args.start);response["next_start"]=json!(args.start+records.len());response["done"]=json!(args.start+records.len()==collection.len());response["records"]=json!(page);
            response["reference_start"]=json!(args.reference_start);response["next_reference_start"]=json!(reference_end);response["references_done"]=json!(reference_end==references.len());
            response["references"]=json!(references[args.reference_start..reference_end].iter().map(|r|json!({"name":read_name(Some(&r.name)),"name_truncated":r.name.chars().count()>256,"length":r.length,"circular":r.circular,"authority":"imported dictionary assertion"})).collect::<Vec<_>>());
            output(response)
        }).await
    }
    #[tool(
        description = "Check one sealed imported SAM/BAM record against an explicitly selected immutable DNA DocumentRef and exact header reference name. Expected reference revision is taken from that DocumentRef; wrong/stale dictionary mappings are rejected. Returns a full source/record/reference-bound JSON artifact plus bounded column preview, distinguishing actual observed agreement/mismatch, reference-derived '=' symbols, insertions/deletions/skips/clips and missing sequence/quality. Does not validate MAPQ/M5/UR/MD/NM, infer unique placement, coverage, consensus, variants or experimental truth. No reference fetch."
    )]
    async fn check_read_alignment(
        &self,
        Parameters(args): Parameters<ReadAlignmentCheck>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.preview_columns > 100 || args.reference_name.is_empty() {
            return Err(bad(
                "explicit reference name and at most100 preview columns required",
            ));
        }
        let options = collagene_core::read_alignment_reference::Options {
            reference_name: args.reference_name.clone(),
            expected_reference_revision: args.reference.revision.clone(),
            max_columns: args.max_columns,
            max_work: args.max_work,
            max_report_bytes: args.max_report_bytes,
            ..Default::default()
        };
        options.validate().map_err(bad)?;
        let server = self.clone();
        self.compute(move|| {
            let (collection,reference)={
                let store=server.store.lock().unwrap();let collection=store.read_alignments(&args.collection)?;
                let reference=store.document(&args.reference)?;
                collagene_core::model::bounds::model_preflight(&[reference]).map_err(bad)?;
                let record=collection.record(args.record_index).ok_or_else(||bad("unknown read alignment record index"))?;
                let named=collection.header().reference(&args.reference_name).ok_or_else(||bad("explicit reference name is absent from imported dictionary"))?;
                if record.fields().reference.as_ref().is_some_and(|name|name!=&args.reference_name)
                    || named.length as usize!=reference.sequence.len()
                    || named.circular!=(reference.topology==collagene_core::Topology::Circular) {
                    return Err(bad("read alignment dictionary/reference mapping differs"));
                }
                (collection,reference.clone())
            };
            let checked=collection.check_reference(args.record_index,&reference,&options).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json(&checked)?,"application/json","read-alignment-reference-check.json")?;
            let report=&checked.report;
            output(json!({"collection":args.collection,"record_index":args.record_index,"reference":args.reference,
                "source_blake3":checked.source_blake3,"record_blake3":checked.record_blake3,"engine_version":report.engine_version,
                "revision":report.revision,"binding":report.binding,"reference_bound":report.reference_bound,"complete":report.complete,
                "incomplete":report.incomplete,"summary":report.summary,"warnings":report.warnings,
                "sequence_present":report.sequence_present,"qualities_present":report.qualities_present,"reverse":report.reverse,
                "run_count":report.runs.len(),"column_count":report.columns.len(),"columns":report.columns.iter().take(args.preview_columns).collect::<Vec<_>>(),
                "columns_truncated":report.columns.len()>args.preview_columns,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Export an immutable SAM/BAM collection as exact original bytes or through bounded SAM/BAM writers. Original retains lexical SAM and exact BAM/BGZF framing; conversions expose loss/framing notices and preserve missing SEQ/QUAL. Returns a standard binary artifact read in bounded base64 chunks. Conversion/capacity failure publishes no partial artifact; no arbitrary filesystem paths or external samtools."
    )]
    async fn export_read_alignments(
        &self,
        Parameters(args): Parameters<ReadAlignmentExport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_bytes == 0 || args.max_bytes > MAX_STORE_BYTES {
            return Err(bad("read alignment export cap must be1..64MiB"));
        }
        let server = self.clone();
        self.compute(move|| {
            use collagene_core::read_alignment_store::{Format,ImportLimits};
            let collection=server.store.lock().unwrap().read_alignments(&args.collection)?;
            let (bytes,format,notices)=match args.format {
                ReadAlignmentExportFormat::Original=>{
                    if collection.original_bytes().len()>args.max_bytes {return Err(bad("original read alignment export exceeds byte cap"));}
                    let format=if collection.header().original_bam_header().is_some(){Format::Bam}else{Format::Sam};
                    (collection.original_bytes().to_vec(),format,Vec::<String>::new())
                }
                mode=>{
                    let format=match mode{ReadAlignmentExportFormat::Sam=>Format::Sam,ReadAlignmentExportFormat::Bam=>Format::Bam,_=>unreachable!()};
                    let exported=collection.export(format,ImportLimits::default().reader,args.max_bytes).map_err(bad)?;
                    (exported.bytes,format,exported.notices.iter().map(|n|format!("{n:?}")).collect())
                }
            };
            let length=bytes.len();let format=match format{Format::Sam=>"sam",Format::Bam=>"bam"};
            let artifact=server.binary_artifact(bytes,"application/octet-stream",if format=="sam"{"read-alignments.sam"}else{"read-alignments.bam"})?;
            output(json!({"collection":args.collection,"source_blake3":collection.source_blake3(),"format":format,
                "original":matches!(args.format,ReadAlignmentExportFormat::Original),"bytes":length,"notices":notices,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Import aligned FASTA, CLUSTAL, Stockholm, strict/relaxed sequential/interleaved PHYLIP, MSF, NEXUS, PIR/NBRF, GDE flat or SELEX as exact base64 bytes with explicit DNA/RNA/protein alphabet. Strict bounded parsing retains original bytes, stable rows, gaps and supported metadata; no alignment engine is invoked. Creates immutable session evidence only."
    )]
    async fn import_alignment(
        &self,
        Parameters(args): Parameters<AlignmentImport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > alignment_formats::MAX_ALIGNMENT_BYTES.div_ceil(3) * 4 {
            return Err(bad("alignment exceeds8MiB decoded"));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            server.insert_alignment(
                alignment_formats::import_alignment(&bytes, args.format, args.alphabet, &args.name)
                    .map_err(bad)?,
            )
        })
        .await
    }
    #[tool(
        description = "Import a complete native alignment JSON attachment as base64, at most32MiB decoded. Source hashes, retained exact bytes, stable row identities and column metadata are validated. Does not edit browser/workspace documents."
    )]
    async fn import_alignment_document(
        &self,
        Parameters(args): Parameters<AlignmentDocumentImport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > (32 * 1024 * 1024usize).div_ceil(3) * 4 {
            return Err(bad("alignment attachment exceeds32MiB decoded"));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            server.insert_alignment(alignment_formats::load_alignment(&bytes).map_err(bad)?)
        })
        .await
    }
    #[tool(
        description = "List immutable session alignment evidence, at most128 collections. Summaries omit rows and raw source; read_alignment pages rows/columns and export_alignment retains full evidence."
    )]
    fn list_alignments(
        &self,
        Parameters(_): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.store.lock().unwrap();
        output(
            json!({"alignments":store.alignments.iter().map(|(id,a)|json!({"id":id,"revision":id,"name":a.name.chars().take(256).collect::<String>(),"alphabet":a.alphabet,"row_count":a.rows.len(),"column_count":a.column_count(),"source_blake3":a.source.blake3})).collect::<Vec<_>>()}),
        )
    }
    #[tool(
        description = "Read a bounded immutable alignment window: up to50 rows and1000 columns, zero-based offsets. Returns as-stored row residue indices or null gaps, never external reference coordinates. Imported consensus symbols are source evidence, not computed confidence. Exact alignment revision required."
    )]
    async fn read_alignment(
        &self,
        Parameters(args): Parameters<AlignmentRegionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let doc=server.store.lock().unwrap().alignment(&args.alignment)?;
            if args.max_rows==0||args.max_rows>50||args.max_columns==0||args.max_columns>1000||args.row_start>=doc.rows.len()||args.column_start>=doc.column_count(){return Err(bad("invalid alignment window (max50rows/1000columns)"));}
            let end=args.column_start.saturating_add(args.max_columns).min(doc.column_count());let row_end=args.row_start.saturating_add(args.max_rows).min(doc.rows.len());
            let rows=doc.rows[args.row_start..row_end].iter().map(|row|{
                let mut pos=row.aligned.as_bytes()[..args.column_start].iter().filter(|b|**b!=b'-').count();
                let positions=row.aligned.as_bytes()[args.column_start..end].iter().map(|b|if *b==b'-'{None}else{let p=pos;pos+=1;Some(p)}).collect::<Vec<_>>();
                json!({"id":row.id,"identifier":row.identifier.chars().take(256).collect::<String>(),"description":row.description.chars().take(256).collect::<String>(),"text_preview_truncated":row.identifier.chars().count()>256||row.description.chars().count()>256,"aligned":&row.aligned[args.column_start..end],"ungapped_positions":positions,"orientation":"AsStored"})
            }).collect::<Vec<_>>();
            output(json!({"alignment":args.alignment,"alphabet":doc.alphabet,"row_start":args.row_start,"row_end":row_end,"row_count":doc.rows.len(),"column_start":args.column_start,"column_end":end,"column_count":doc.column_count(),"rows":rows,"next_row_start":if row_end<doc.rows.len(){Some(row_end)}else{None},"next_column_start":if end<doc.column_count(){Some(end)}else{None},"clustal_consensus":doc.metadata.clustal_consensus.as_ref().map(|s|&s[args.column_start..end]),"source_blake3":doc.source.blake3,"coordinates":"zero-based alignment columns and ungapped row as stored; not an external reference map"}))
        }).await
    }
    #[tool(
        description = "Export a complete immutable alignment as aligned FASTA, CLUSTAL, Stockholm, strict/relaxed sequential/interleaved PHYLIP, MSF, NEXUS, PIR/NBRF, GDE flat or SELEX with explicit losses, Native validated JSON attachment, or exact Original source bytes. PHYLIP identifier truncation requires explicit options and reports renamings; collisions are rejected. Returns a UTF-8 artifact read with read_artifact chunks. No silent row renaming or external aligner. Native retains all supported evidence."
    )]
    async fn export_alignment(
        &self,
        Parameters(args): Parameters<AlignmentExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let doc=server.store.lock().unwrap().alignment(&args.alignment)?;
            let (bytes,losses,renamings,filename)=match args.format{
                AlignmentExportFormat::Native=>{
                    if args.options.is_some(){return Err(bad("Native alignment export has no exchange-layout options"));}
                    (alignment_formats::save_alignment(&doc).map_err(bad)?,vec![],std::collections::BTreeMap::new(),"alignment.cg2aln")
                },
                AlignmentExportFormat::Original=>{
                    if args.options.is_some(){return Err(bad("Original source export has no exchange-layout options"));}
                    (doc.original_bytes().map_err(bad)?,vec![],std::collections::BTreeMap::new(),"alignment-original.txt")
                },
                format=>{
                    use alignment_formats::AlignmentFormat as F;
                    let format=match format{
                        AlignmentExportFormat::AlignedFasta=>F::AlignedFasta,
                        AlignmentExportFormat::Clustal=>F::Clustal,
                        AlignmentExportFormat::Stockholm=>F::Stockholm,
                        AlignmentExportFormat::PhylipSequential=>F::PhylipSequential,
                        AlignmentExportFormat::PhylipInterleaved=>F::PhylipInterleaved,
                        AlignmentExportFormat::PhylipRelaxedSequential=>F::PhylipRelaxedSequential,
                        AlignmentExportFormat::PhylipRelaxedInterleaved=>F::PhylipRelaxedInterleaved,
                        AlignmentExportFormat::Msf=>F::Msf,
                        AlignmentExportFormat::Nexus=>F::Nexus,
                        AlignmentExportFormat::Pir=>F::Pir,
                        AlignmentExportFormat::Gde=>F::Gde,
                        AlignmentExportFormat::Selex=>F::Selex,
                        _=>unreachable!()
                    };
                    let report=if let Some(options)=&args.options{
                        collagene_core::alignment_exchange::export_alignment_with_options(&doc,format,options).map_err(bad)?
                    }else{alignment_formats::export_alignment(&doc,format).map_err(bad)?};
                    let filename=match format{F::AlignedFasta=>"alignment.fasta",F::Clustal=>"alignment.aln",F::Stockholm=>"alignment.sto",F::Msf=>"alignment.msf",F::Nexus=>"alignment.nex",F::Pir=>"alignment.pir",F::Gde=>"alignment.gde",F::Selex=>"alignment.slx",_=>"alignment.phy"};
                    (report.text.into_bytes(),report.losses,report.renamings,filename)
                }
            };
            let hash=blake3::hash(&bytes).to_hex().to_string();let size=bytes.len();let text=String::from_utf8(bytes).map_err(bad)?;
            let artifact=server.text_artifact(text,"text/plain",filename)?;
            output(json!({"alignment":args.alignment,"artifact_id":artifact,"bytes":size,"content_blake3":hash,"losses":losses,"renamings":renamings,"encoding":"utf-8","source_is_current":doc.source_is_current().map_err(bad)?}))
        }).await
    }
    #[tool(
        description = "Explicitly extract selected alignment rows into new immutable session molecules, removing gaps in as-stored orientation. Original annotations, consensus/confidence and external coordinates are not invented. All-gap, unknown or duplicate selections fail atomically. Does not modify browser/workspace documents."
    )]
    async fn extract_alignment_rows(
        &self,
        Parameters(args): Parameters<AlignmentRowsQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.row_ids.is_empty() || args.row_ids.len() > alignment_formats::MAX_ALIGNMENT_ROWS {
            return Err(bad("select1..1000 alignment rows"));
        }
        let server = self.clone();
        self.compute(move||{let doc=server.store.lock().unwrap().alignment(&args.alignment)?;let rows=doc.extract_rows(&args.row_ids).map_err(bad)?;server.insert_documents(rows,vec!["Extracted ungapped rows retain source/current alignment hashes; no external reference coordinates, alignment confidence or annotations were copied.".into()])}).await
    }
    #[tool(
        description = "Search an immutable document with explicit sequence, annotation, saved-primer, enzyme or all-source mode. DNA/RNA IUPAC and peptide search preserve document alphabet. Bounded hits/work/rows, strand, circular spans, exact source revision and completeness remain explicit. Readonly: selecting or modifying a document is a separate user action. Full report is a JSON artifact."
    )]
    async fn search_document(
        &self,
        Parameters(args): Parameters<SearchDocumentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let molecule=server.store.lock().unwrap().document(&args.document)?.clone();
            let report=collagene_core::search::search_document(&molecule,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&serde_json::to_value(&report).map_err(bad)?)?;
            output(json!({"document":args.document,"source_revision":report.source_revision,"engine_version":report.engine_version,"options":report.options,"complete":report.complete,"work":report.work,"rows":report.rows,"hit_count":report.hits.len(),"hits":report.hits.iter().take(100).map(|h|json!({"id":h.id,"kind":h.kind,"label":h.label.chars().take(256).collect::<String>(),"label_preview_truncated":h.label.chars().count()>256,"parts":h.parts.iter().take(64).collect::<Vec<_>>(),"parts_preview_truncated":h.parts.len()>64,"feature_index":h.feature_index,"enzyme":h.enzyme,"top_cut":h.top_cut,"bottom_cut":h.bottom_cut,"certain":h.certain})).collect::<Vec<_>>(),"preview_truncated":report.hits.len()>100,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Inspect bounded annotation, saved-primer, enzyme/cut and translation tracks for an exact immutable document and half-open viewport. Circular windows may cross the origin once. Uses the same deterministic track engine as the GUI; returns a bounded preview and full JSON artifact without changing selections or annotations."
    )]
    async fn view_tracks(
        &self,
        Parameters(args): Parameters<ViewTracksQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let molecule = server
                .store
                .lock()
                .unwrap()
                .document(&args.document)?
                .clone();
            let report =
                collagene_core::view_tracks::view_tracks(&molecule, args.window, &args.options)
                    .map_err(bad)?;
            let value = serde_json::to_value(&report).map_err(bad)?;
            let artifact = server.analysis_artifact(&value)?;
            let mut preview = value;
            if let Some(items) = preview.get_mut("items").and_then(Value::as_array_mut) {
                let count = items.len();
                items.truncate(100);
                let mut residue_budget = 300usize;
                for item in items {
                    if let Some(parts) = item.get_mut("parts").and_then(Value::as_array_mut) {
                        let truncated = parts.len() > 64;
                        parts.truncate(64);
                        item["parts_preview_truncated"] = json!(truncated);
                    }
                    let mut take = 0;
                    if let Some(residues) = item.get_mut("residues").and_then(Value::as_array_mut) {
                        take = residues.len().min(residue_budget);
                        residue_budget -= take;
                        let truncated = residues.len() > take;
                        residues.truncate(take);
                        item["residues_preview_truncated"] = json!(truncated);
                    }
                    if let Some(protein) = item.get("protein").and_then(Value::as_str) {
                        let text = protein.chars().take(take).collect::<String>();
                        item["protein"] = json!(text);
                    }
                }

                preview["item_count"] = json!(count);
                preview["preview_truncated"] = json!(count > 100);
            }
            preview["document"] = json!(args.document);
            preview["artifact_id"] = json!(artifact);
            output(preview)
        })
        .await
    }

    #[tool(
        description = "Create an immutable, versioned DNA feature library from explicitly selected features in your imported source documents. Biological extraction respects reverse/compound locations. Source revisions and review state remain evidence; no proprietary database or automatic target edits."
    )]
    async fn create_feature_library(
        &self,
        Parameters(args): Parameters<CreateFeatureLibrary>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.features.len() > 10_000 {
            return Err(bad("at most 10000 source feature selections"));
        }
        let server = self.clone();
        self.compute(move || {
            let sources = {
                let store = server.store.lock().unwrap();
                let mut sources = BTreeMap::new();
                let mut seen = std::collections::BTreeSet::new();
                for feature in &args.features {
                    let molecule = store.document(&feature.document)?;
                    if !seen.insert((feature.document.id.clone(), feature.feature_index)) {
                        return Err(bad("duplicate source feature selection"));
                    }
                    sources
                        .entry(feature.document.id.clone())
                        .or_insert_with(|| molecule.clone());
                }
                sources
            };
            let mut library =
                feature_library::FeatureLibrary::new(args.name, args.version).map_err(bad)?;
            for source in &args.features {
                library
                    .add_feature(&sources[&source.document.id], source.feature_index)
                    .map_err(bad)?;
            }
            server.insert_feature_library(library)
        })
        .await
    }
    #[tool(
        description = "Import a validated versioned feature-library JSON file encoded as base64 (at most 8 MiB decoded). Sequence, source identities, entry revisions and supported semantics are validated. Untrusted labels remain data; importing does not annotate documents."
    )]
    async fn import_feature_library(
        &self,
        Parameters(args): Parameters<ImportFeatureLibrary>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > (8 * 1024 * 1024usize).div_ceil(3) * 4 {
            return Err(bad("feature library import exceeds 8 MiB decoded"));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            if bytes.len() > 8 * 1024 * 1024 {
                return Err(bad("feature library import exceeds 8 MiB decoded"));
            }
            server.insert_feature_library(feature_library::load_library(&bytes).map_err(bad)?)
        })
        .await
    }
    #[tool(
        description = "List immutable session feature libraries with exact references and bounded source/version summaries. No filesystem or unsynchronized browser access."
    )]
    fn list_feature_libraries(
        &self,
        Parameters(_): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.store.lock().unwrap();
        output(
            json!({"libraries":store.feature_libraries.iter().map(|(id,l)|json!({"id":id,"revision":id,"name":l.name.chars().take(256).collect::<String>(),"version":l.version.chars().take(256).collect::<String>(),"entry_count":l.entries.len()})).collect::<Vec<_>>()}),
        )
    }
    #[tool(
        description = "Export a complete immutable feature library as validated JSON evidence. Exact library id/revision required; artifact chunks retain every source and entry rather than only previews."
    )]
    async fn export_feature_library(
        &self,
        Parameters(args): Parameters<FeatureLibraryQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let library = server
                .store
                .lock()
                .unwrap()
                .feature_library(&args.library)?
                .clone();
            let bytes = feature_library::save_library(&library).map_err(bad)?;
            let length = bytes.len();
            let text = String::from_utf8(bytes).map_err(bad)?;
            let artifact =
                server.text_artifact(text, "application/json", "feature-library.json")?;
            output(json!({"library":args.library,"bytes":length,"artifact_id":artifact}))
        })
        .await
    }
    #[tool(
        description = "Scan an exact DNA document against an immutable feature library on explicit strands with bounded mismatch/work/hit limits. Full-entry coverage, concrete versus ambiguity-compatible identities, repeats and completeness are explicit. Returns a revision-bound proposal set and JSON artifact; no annotations are automatically installed."
    )]
    async fn scan_feature_library(
        &self,
        Parameters(args): Parameters<FeatureLibraryScan>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (document,library)={let store=server.store.lock().unwrap();(store.document(&args.document)?.clone(),store.feature_library(&args.library)?.clone())};
            let report=feature_library::scan(&document,&library,&args.options).map_err(bad)?;
            let text=serde_json::to_string(&report).map_err(bad)?;let id=blake3::hash(text.as_bytes()).to_hex().to_string();
            let result=json!({"document":args.document,"library":args.library,"scan":{"id":id,"revision":id},"source_revision":report.source_revision,"library_revision":report.library_revision,"engine_version":report.engine_version,"complete":report.complete,"comparisons":report.comparisons,"options":report.options,"hit_count":report.hits.len(),"hits":report.hits.iter().take(100).map(|hit|json!({"id":hit.id,"entry_id":hit.entry_id,"entry_revision":hit.entry_revision,"start":hit.start,"length":hit.length,"reverse":hit.reverse,"concrete_matches":hit.concrete_matches,"compatible_ambiguous":hit.compatible_ambiguous,"mismatches":hit.mismatches,"identity_lower":hit.identity_lower,"identity_upper":hit.identity_upper,"coverage":hit.coverage,"already_annotated":hit.already_annotated,"name":hit.feature.label().chars().take(256).collect::<String>(),"kind":hit.feature.kind,"needs_review":hit.feature.needs_review})).collect::<Vec<_>>(),"hits_truncated":report.hits.len()>100,"warnings":report.warnings.iter().take(100).map(|w|w.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),"warnings_truncated":report.warnings.len()>100,"artifact_id":id,"applied":false});
            let mut store=server.store.lock().unwrap();let new_report=!store.feature_scans.contains_key(&id);let new_artifact=!store.artifacts.contains_key(&id);
            if new_report&&store.feature_scans.len()>=MAX_ITEMS{return Err(bad("at most 128 feature scans"));}if new_artifact&&store.artifacts.len()>=MAX_ITEMS{return Err(bad("at most 128 artifacts"));}
            let bytes=text.len().checked_mul(usize::from(new_report)+usize::from(new_artifact)).ok_or_else(||bad("session accounting overflow"))?;store.reserve(bytes)?;
            store.feature_scans.entry(id.clone()).or_insert(report);store.artifacts.entry(id).or_insert(Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"feature-library-scan.json"});output(result)
        }).await
    }
    #[tool(
        description = "Apply explicitly selected proposal IDs from an immutable feature-library scan to a NEW session document. Both exact library and document revisions are required; the core reconstructs the scan and rejects tampering, stale evidence and duplicates. Inferred annotations require review. Original documents and physical chemistry remain unchanged."
    )]
    async fn apply_feature_library(
        &self,
        Parameters(args): Parameters<FeatureLibraryApply>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{let (document,library,report)={let store=server.store.lock().unwrap();(store.document(&args.document)?.clone(),store.feature_library(&args.library)?.clone(),store.feature_scan(&args.scan)?.clone())};let annotated=feature_library::apply_selected(&document,&library,&report,&args.selected).map_err(bad)?;server.insert_documents(vec![annotated],vec!["Feature proposals installed only in a new immutable document; inferred annotations still require review.".into()])}).await
    }
    #[tool(
        description = "Parse GFF3, GTF or BED annotation evidence using an explicit source seqid and an exact imported reference revision. Coordinate/phase/strand/compound semantics are checked by the embedded parser; unrelated seqids remain explicit. Returns an immutable proposal dataset. No existing annotations are changed; <=8MiB source."
    )]
    async fn import_annotations(
        &self,
        Parameters(args): Parameters<AnnotationImport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.text.len() > 8 * 1024 * 1024 {
            return Err(bad("annotation text exceeds 8 MiB"));
        }
        let server = self.clone();
        self.compute(move||{
            let document=server.store.lock().unwrap().document(&args.document)?.clone();
            let report=annotation_formats::import_annotations(&document,&args.text,args.format,&args.seqid).map_err(bad)?;
            let bytes=serde_json::to_vec(&report).map_err(bad)?;let id=blake3::hash(&bytes).to_hex().to_string();
            let result=json!({"dataset":{"id":id,"revision":id},"document":args.document,"format":report.format,"seqid":report.seqid,"source_blake3":report.source_blake3,"reference_revision":report.reference_revision,"proposal_count":report.proposals.len(),"proposals":report.proposals.iter().take(100).enumerate().map(|(i,p)|json!({"index":i,"id":p.id,"name":p.feature.label().chars().take(256).collect::<String>(),"kind":p.feature.kind,"needs_review":p.feature.needs_review,"spans":p.feature.location.spans().iter().take(20).collect::<Vec<_>>(),"spans_truncated":p.feature.location.spans().len()>20,"source_rows":p.source_rows.iter().take(100).collect::<Vec<_>>(),"source_rows_truncated":p.source_rows.len()>100})).collect::<Vec<_>>(),"proposals_truncated":report.proposals.len()>100,"skipped_seqids":report.skipped_seqids.iter().take(100).map(|s|s.chars().take(256).collect::<String>()).collect::<Vec<_>>(),"skipped_seqids_truncated":report.skipped_seqids.len()>100,"warnings":report.warnings.iter().take(100).map(|w|w.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),"warnings_truncated":report.warnings.len()>100,"applied":false});
            let mut store=server.store.lock().unwrap();if !store.annotation_datasets.contains_key(&id){if store.annotation_datasets.len()>=MAX_ITEMS{return Err(bad("at most 128 annotation datasets"));}store.reserve(bytes.len())?;store.annotation_datasets.insert(id,report);}output(result)
        }).await
    }
    #[tool(
        description = "Apply explicitly selected indices from an immutable GFF3/GTF/BED dataset to a NEW document at its exact bound reference revision. Revalidates retained source and selected coordinates; no silent replacement, external write or inherited scientific verification."
    )]
    async fn apply_annotations(
        &self,
        Parameters(args): Parameters<AnnotationApply>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{let (document,report)={let store=server.store.lock().unwrap();(store.document(&args.document)?.clone(),store.annotation_dataset(&args.dataset)?.clone())};let annotated=annotation_formats::apply_annotations(&document,&report,&args.selected).map_err(bad)?;server.insert_documents(vec![annotated],vec!["Selected source annotations appended in a new immutable document; source evidence is retained and no cloning verification is implied.".into()])}).await
    }
    #[tool(
        description = "Export supported annotations from an exact document as GFF3, GTF or BED with an explicit output seqid. Phase, discontinuous/reverse locations and format losses follow the embedded strict profile; unsupported semantics fail rather than being invented. Full text artifact includes loss diagnostics separately."
    )]
    async fn export_annotations(
        &self,
        Parameters(args): Parameters<AnnotationExport>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{let document=server.store.lock().unwrap().document(&args.document)?.clone();let report=annotation_formats::export_annotations(&document,args.format,&args.seqid).map_err(bad)?;let length=report.text.len();let filename=match args.format { annotation_formats::AnnotationFormat::Gff3=>"annotations.gff3",annotation_formats::AnnotationFormat::Gtf=>"annotations.gtf",annotation_formats::AnnotationFormat::Bed=>"annotations.bed" };let artifact=server.text_artifact(report.text,"text/plain",filename)?;output(json!({"document":args.document,"format":args.format,"seqid":args.seqid,"bytes":length,"warnings":report.warnings,"artifact_id":artifact}))}).await
    }
    #[tool(
        description = "Export the complete immutable GFF3/GTF/BED proposal dataset and exact original source bytes as two artifacts. Includes explicit reference/seqid mapping and every row/proposal; this is evidence, not an automatic document edit or trusted verification."
    )]
    async fn export_annotation_dataset(
        &self,
        Parameters(args): Parameters<AnnotationDatasetQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let report=server.store.lock().unwrap().annotation_dataset(&args.dataset)?.clone();
            let text=serde_json::to_string(&report).map_err(bad)?;let original=report.original_bytes().map_err(bad)?;
            let text_id=blake3::hash(text.as_bytes()).to_hex().to_string();let original_id=blake3::hash(&original).to_hex().to_string();
            let mut store=server.store.lock().unwrap();let new_text=!store.artifacts.contains_key(&text_id);let new_original=!store.artifacts.contains_key(&original_id);
            if let Some(existing)=store.artifacts.get(&original_id){if !matches!(existing.data,ArtifactData::Binary(_)){return Err(bad("original source artifact already has another representation"));}}
            let added=usize::from(new_text)+usize::from(new_original);if store.artifacts.len()+added>MAX_ITEMS{return Err(bad("annotation artifacts exceed item limit; no outputs created"));}
            let bytes=(if new_text{text.len()}else{0usize}).checked_add(if new_original{original.len()}else{0}).ok_or_else(||bad("session accounting overflow"))?;store.reserve(bytes)?;
            store.artifacts.entry(text_id.clone()).or_insert(Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"annotation-dataset.json"});
            store.artifacts.entry(original_id.clone()).or_insert(Artifact{data:ArtifactData::Binary(original),mime_type:"application/octet-stream",encoding:"base64",filename:match report.format { annotation_formats::AnnotationFormat::Gff3=>"original-annotations.gff3",annotation_formats::AnnotationFormat::Gtf=>"original-annotations.gtf",annotation_formats::AnnotationFormat::Bed=>"original-annotations.bed" }});
            output(json!({"dataset":args.dataset,"dataset_artifact_id":text_id,"original_artifact_id":original_id}))
        }).await
    }
    #[tool(
        description = "Validate revision-bound original SCF/ABIF/explicit-encoding FASTQ reads against an independent DNA reference with embedded trimming, affine placement/uniqueness proof, coverage and variant evidence. Native trace scores are not invented Phred. Full JSON artifact; no physical/sample confirmation or automatic document mutation."
    )]
    async fn validate_reads(
        &self,
        Parameters(args): Parameters<ReadValidationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (reference,reads)={let store=server.store.lock().unwrap();(store.document(&args.reference)?.clone(),args.reads.iter().map(|r|store.document(r).cloned()).collect::<Result<Vec<_>,_>>()?)};
            let report=trace_validation::validate_reads(&reference,&reads,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"reference":args.reference,"reads":args.reads,"validation":report}))?;
            output(json!({"reference":args.reference,"read_revisions":report.read_revisions,"model":report.model,"engine_version":report.engine_version,"status":report.status,"eligible_coverage":report.eligible_coverage,"reference_length":report.consensus.len(),"reference_coverage_fraction":report.reference_coverage_fraction,"consensus_preview":report.consensus.chars().take(2000).collect::<String>(),"uncovered_ranges":report.uncovered_ranges.iter().take(100).collect::<Vec<_>>(),"uncovered_ranges_truncated":report.uncovered_ranges.len()>100,"differences":report.sites.iter().filter(|s|s.state!=trace_validation::SiteState::SupportedMatch).take(100).collect::<Vec<_>>(),"differences_truncated":report.sites.iter().filter(|s|s.state!=trace_validation::SiteState::SupportedMatch).count()>100,"reads":report.reads.iter().map(|r|json!({"name":r.name.chars().take(256).collect::<String>(),"read_format":r.read_format,"quality_source":r.trim.as_ref().map(|t|t.quality_source),"trim":r.trim.as_ref().map(|t|json!({"start":t.start,"end":t.end,"passing_calls":t.passing_calls,"low_or_unknown_count":t.low_or_unknown_calls.len()})),"placement_unique_proven":r.placement_unique_proven,"optimality_proven":r.optimality_proven,"search_complete":r.search_complete,"paired_reference_fraction":r.paired_reference_fraction,"error":r.error,"warnings":r.warnings})).collect::<Vec<_>>(),"insertion_observations":report.insertions.iter().take(100).map(|i|json!({"reference_boundary":i.reference_boundary,"sequence_preview":i.sequence.chars().take(100).collect::<String>(),"length":i.sequence.len(),"read_index":i.read_index,"eligible":i.eligible})).collect::<Vec<_>>(),"insertions_truncated":report.insertions.len()>100,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Import FASTQ caller evidence with required explicit Phred33 or Phred64 encoding. No ambiguous auto-detection, quality fabrication or automatic verification. Atomic immutable documents; text <=8 MiB."
    )]
    async fn import_fastq(
        &self,
        Parameters(args): Parameters<FastqImport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.text.len() > 8 * 1024 * 1024 {
            return Err(bad("FASTQ MCP import exceeds 8 MiB"));
        }
        let server = self.clone();
        self.compute(move || {
            let report =
                oligo_formats::import_fastq(args.text.as_bytes(), args.encoding).map_err(bad)?;
            server.insert_documents(report.molecules, report.warnings)
        })
        .await
    }
    #[tool(
        description = "Export exact called FASTQ with mapped caller qualities in a required explicit encoding. Stale qualities or unrepresentable scores fail instead of being invented/clamped. Original source attachment is preserved in native documents."
    )]
    async fn export_fastq(
        &self,
        Parameters(args): Parameters<FastqExport>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let molecule=server.store.lock().unwrap().document(&args.document)?.clone();
            let report=oligo_formats::export_fastq(&molecule,args.encoding).map_err(bad)?;
            let artifact=server.text_artifact(report.text.clone(),"text/plain","reads.fastq")?;
            output(json!({"document":args.document,"encoding":args.encoding,"bytes":report.text.len(),"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Import a CSV/TSV oligo collection using explicit column mapping. Retain the whole source table once, exact rows/modifications/units, and explicit selected source-row indices. Only supported selected primers become ssDNA reference documents; unsupported modifications stay evidence. Atomic; <=8 MiB source."
    )]
    async fn import_oligos(
        &self,
        Parameters(args): Parameters<OligoImport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.text.len() > oligo_formats::MAX_TABLE_BYTES {
            return Err(bad("oligo table exceeds 8 MiB"));
        }
        let server = self.clone();
        self.compute(move||{
            let table=oligo_formats::import_oligos(&args.text,&args.options).map_err(bad)?;
            let selected=table.selected_molecules(&args.selection).map_err(bad)?;
            let serialized=serde_json::to_vec(&table).map_err(bad)?;
            let id=blake3::hash(&serialized).to_hex().to_string();
            let mut pending=BTreeMap::new();for m in selected.molecules{pending.insert(m.revision(),m);}
            let selection=StoredSelection{collection_id:id.clone(),source_rows:args.selection.clone()};
            let selection_bytes=serde_json::to_vec(&selection).map_err(bad)?;
            let selection_id=blake3::hash(&selection_bytes).to_hex().to_string();
            let mut store=server.store.lock().unwrap();
            if !store.oligo_selections.contains_key(&selection_id)&&store.oligo_selections.len()>=MAX_ITEMS{return Err(bad("at most 128 oligo row selections"));}
            if !store.oligos.contains_key(&id)&&store.oligos.len()>=MAX_ITEMS{return Err(bad("at most 128 oligo collections"));}
            let (count,bytes)=pending.iter().filter(|(r,_)|!store.documents.contains_key(*r)).try_fold((0usize,(if store.oligos.contains_key(&id){0}else{serialized.len()})+(if store.oligo_selections.contains_key(&selection_id){0}else{selection_bytes.len()})),|(n,b),(_,m)|Ok::<_,ErrorData>((n+1,b.checked_add(serde_json::to_vec(m).map_err(bad)?.len()).ok_or_else(||bad("session accounting overflow"))?)))?;
            if store.documents.len()+count>MAX_DOCUMENT_REVISIONS{return Err(bad("at most 10000 session document revisions"));}
            let rows=table.records.iter().take(100).map(|r|json!({"row_index":r.row_index,"row_id":r.row_id,"name":r.name.chars().take(256).collect::<String>(),"dna_supported":r.primer.is_some(),"modification":r.modification.as_ref().map(|m|m.chars().take(512).collect::<String>()),"diagnostics":r.diagnostics.iter().take(10).map(|s|s.chars().take(512).collect::<String>()).collect::<Vec<_>>()})).collect::<Vec<_>>();
            let result=json!({"collection":{"id":id,"revision":id},"source_blake3":table.source_blake3,"row_count":table.records.len(),"rows":rows,"rows_truncated":table.records.len()>100,"selection":args.selection,"selection_id":selection_id,"documents":pending.iter().take(128).map(|(r,m)|summary(r,m)).collect::<Vec<_>>(),"documents_truncated":pending.len()>128,"warnings":selected.warnings,"table_warnings":table.warnings.iter().take(100).collect::<Vec<_>>(),"persistence":"session only"});
            store.reserve(bytes)?;store.documents.extend(pending);store.oligos.entry(id).or_insert(table);store.oligo_selections.entry(selection_id).or_insert(selection);
            output(result)
        }).await
    }
    #[tool(
        description = "Export an immutable complete oligo collection with exact original-table bytes, a normalized full CSV/TSV table and JSON evidence retaining explicit selected source rows. CSV contains ALL collection rows; selection is preserved separately in JSON, not used to fabricate a reduced source table. SpreadsheetSafe or MachineData must be explicit."
    )]
    async fn export_oligos(
        &self,
        Parameters(args): Parameters<OligoExport>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let table={let store=server.store.lock().unwrap();if args.collection.id!=args.collection.revision{return Err(bad("stale oligo collection revision"));}store.oligos.get(&args.collection.id).cloned().ok_or_else(||bad("unknown oligo collection"))?};
            if let Some(rows)=&args.selection{let mut unique=std::collections::BTreeSet::new();if rows.iter().any(|&r|r>=table.records.len()||!unique.insert(r)){return Err(bad("selected collection rows must be unique existing source rows"));}}
            let report=oligo_formats::export_oligos(&table,args.delimiter,args.mode).map_err(bad)?;
            // Preflight every output against aggregate limits before creating any artifact.
            let original=table.original_bytes().map_err(bad)?;
            let saved_selections={let store=server.store.lock().unwrap();store.oligo_selections.iter().filter(|(_,s)|s.collection_id==args.collection.id).map(|(id,s)|json!({"id":id,"source_rows":s.source_rows})).collect::<Vec<_>>()};
            let evidence=json!({"collection":args.collection,"selection":args.selection,"saved_selections":saved_selections,"table":table});
            let full=serde_json::to_string(&evidence).map_err(bad)?;
            let total=original.len().checked_add(full.len()).and_then(|v|v.checked_add(report.text.len())).ok_or_else(||bad("artifact accounting overflow"))?;
            {let store=server.store.lock().unwrap();if total>MAX_STORE_BYTES.saturating_sub(store.used)||store.artifacts.len()+3>MAX_ITEMS{return Err(bad("export artifacts exceed session size/item limit; no outputs created"));}}
            let original_id=server.binary_artifact(original,"application/octet-stream","original-oligos.table")?;
            let csv_id=server.text_artifact(report.text,"text/plain","oligos.table")?;
            let collection_id=server.text_artifact(full,"application/json","oligo-collection.json")?;
            output(json!({"collection":args.collection,"selected_source_rows":args.selection,"exports_all_rows":true,"original_artifact_id":original_id,"table_artifact_id":csv_id,"collection_artifact_id":collection_id,"warnings":report.warnings}))
        }).await
    }
    #[tool(
        description = "Import a full native .cg2workspace archive from base64 (<=8 MiB decoded convenience). Histories, candidates, source attachments and exact revisions are retained atomically. Candidates remain unverified data. For up to 64 MiB archives use bounded begin/append/finalize_workspace_upload."
    )]
    async fn import_workspace_archive(
        &self,
        Parameters(args): Parameters<ImportArchive>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > (8 * 1024 * 1024usize).div_ceil(3) * 4 {
            return Err(bad("use chunked upload for archives above 8 MiB decoded"));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            if bytes.len() > 8 * 1024 * 1024 {
                return Err(bad(
                    "small native archive import exceeds 8 MiB; use chunked upload",
                ));
            }
            server.import_archive_bytes(&bytes)
        })
        .await
    }
    #[tool(
        description = "Export a complete immutable native workspace with exact history, undo/redo, candidate data and attachments. Supply exact imported workspace id/revision. Returns a base64 artifact read in bounded read_artifact chunks; no path/write or trusted verification is granted."
    )]
    async fn export_workspace_archive(
        &self,
        Parameters(args): Parameters<WorkspaceRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let value={let store=server.store.lock().unwrap();if args.id!=args.revision{return Err(bad("stale workspace revision"));}store.workspaces.get(&args.id).ok_or_else(||bad("unknown workspace"))?.restore(&store)?};
            let bytes=workspace::save_workspace(&value).map_err(bad)?;let length=bytes.len();let digest=blake3::hash(&bytes).to_hex().to_string();
            let artifact=server.binary_artifact(bytes,"application/octet-stream","workspace.cg2workspace")?;
            output(json!({"workspace":args.id,"revision":args.revision,"bytes":length,"blake3":digest,"encoding":"base64","artifact_id":artifact,"filename":"workspace.cg2workspace","candidates_verified":false}))
        }).await
    }
    #[tool(
        description = "Finalize a saved cloning plan against explicit exact historical input/independent-target events in an immutable uploaded native workspace. Freshly reconstructs and retains ALL native step history; options control final construct name, prepared-output presentation, faithful saved oligo attachment and explicit incomplete acceptance. Atomically returns a NEW session workspace, native archive artifact and full receipt JSON; original workspace and sources are unchanged. Hiding preparations never deletes evidence. No browser/hosted writes, guessed DNA event mapping, experimental-material claim or imported verification authority."
    )]
    async fn finalize_cloning_outputs(
        &self,
        Parameters(args): Parameters<FinalizeCloningOutputs>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.input_bindings.len() > 128 {
            return Err(bad("at most128 explicit original input event bindings"));
        }
        args.options.validate().map_err(bad)?;
        let binding_bounded = |binding: &collagene_core::run_lineage::EventBinding| {
            binding.document_id.len() <= 128
                && binding.event_id.len() <= 128
                && binding.revision.len() == 64
        };
        if !binding_bounded(&args.target_binding)
            || args
                .input_bindings
                .iter()
                .any(|binding| binding.input_id.len() > 128 || !binding_bounded(&binding.event))
        {
            return Err(bad("oversized explicit event binding"));
        }
        let server = self.clone();
        self.compute(move || {
            let value = server.restore_workspace(&args.workspace)?;
            let plan = {
                let store = server.store.lock().unwrap();
                let saved = store
                    .plans
                    .get(&args.plan_revision)
                    .ok_or_else(|| bad("unknown saved plan revision"))?;
                if saved.target.revision() != args.target_binding.revision {
                    return Err(bad(
                        "chosen target event differs from the saved independent target",
                    ));
                }
                saved.plan.clone()
            };
            let expected = collagene_core::cloning_outputs::finalize_revision(
                &value,
                &plan,
                args.target_binding.clone(),
                args.input_bindings.clone(),
                &args.options,
            )
            .map_err(bad)?;
            let commit = collagene_core::cloning_outputs::finalize(
                &value,
                &plan,
                args.target_binding,
                args.input_bindings,
                &args.options,
            )
            .map_err(bad)?;
            commit.validate_binding(&expected).map_err(bad)?;
            let text = bounded_json_limit(&commit, 32 * 1024 * 1024)?;
            let receipt = serde_json::to_value(&commit.receipt).map_err(bad)?;
            let native = workspace::save_workspace(&commit.workspace).map_err(bad)?;
            server.commit_workspace_outputs(commit.workspace, Some((text, receipt)), Some(native))
        })
        .await
    }
    #[tool(
        description = "Freshly check sequencing reads against an exact construct reference event under explicit native-quality, trimming, placement and support rules, and retain the complete evidence in a NEW immutable native workspace. Reference and read event bindings are explicit; equal DNA never chooses events. Original sequences, heads and chromatograms are unchanged. Atomically returns workspace plus full JSON/native artifacts. Saved results are historical evidence, not restored current verification or experimental sample confirmation."
    )]
    async fn record_read_check(
        &self,
        Parameters(args): Parameters<RecordReadCheck>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.reads.is_empty() || args.reads.len() > 32 {
            return Err(bad("choose1..32 exact read event bindings"));
        }
        let server = self.clone();
        self.compute(move || {
            let frame=server.restore_workspace(&args.workspace)?;
            let expected=collagene_core::read_checks::capture_revision(&frame,&args.reference,&args.reads,&args.options).map_err(bad)?;
            let commit=collagene_core::read_checks::capture(&frame,args.reference,args.reads,&args.options).map_err(bad)?;
            commit.validate_binding(&expected).map_err(bad)?;
            let saved=collagene_core::read_checks::load(&commit.workspace,&commit.attachment_hash).map_err(bad)?;
            let text=bounded_json_limit(&saved,17*1024*1024)?;
            let receipt=json!({"attachment_hash":commit.attachment_hash,"context_revision":commit.context_revision,
                "reference":saved.reference,"read_count":saved.reads.len(),"recorded_status":saved.report.status,
                "recorded_engine":saved.report.engine_version,"state":"historical","verified":false});
            let native=workspace::save_workspace(&commit.workspace).map_err(bad)?;
            server.commit_workspace_outputs_named(commit.workspace,Some((text,receipt)),Some(native),"read-check.json",
                "Historical source-event-bound read evidence. Fresh replay reconstructs current evidence; no experimental sample identity or whole-construct confirmation is asserted.")
        }).await
    }
    #[tool(
        description = "List saved construct-owned read checks in an immutable native workspace. Every result retains its exact reference event and read count. Recorded engine/status are historical fields; this does not rerun reads or restore verification authority."
    )]
    async fn list_read_checks(
        &self,
        Parameters(args): Parameters<ListReadChecks>,
    ) -> Result<CallToolResult, ErrorData> {
        let frame = self.restore_workspace(&args.workspace)?;
        let checks = collagene_core::read_checks::list(&frame).map_err(bad)?;
        output(
            json!({"workspace":args.workspace,"checks":checks,"state":"historical","verified":false}),
        )
    }
    #[tool(
        description = "Replay one saved read check by resolving its exact retained reference/read events and freshly recomputing trimming, placement, coverage and native call evidence with the current engine. Does not change source heads, raw chromatograms, saved evidence or the native workspace. Current computational consistency remains distinct from empirical sample identity."
    )]
    async fn replay_read_check(
        &self,
        Parameters(args): Parameters<ReplayReadCheck>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.attachment_hash.len() != 64 {
            return Err(bad("an exact saved read-check attachment hash is required"));
        }
        let server = self.clone();
        self.compute(move || {
            let frame=server.restore_workspace(&args.workspace)?;
            let expected=collagene_core::read_checks::replay_revision(&frame,&args.attachment_hash).map_err(bad)?;
            let checked=collagene_core::read_checks::replay(&frame,&args.attachment_hash).map_err(bad)?;
            checked.validate_binding(&expected).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&checked,32*1024*1024)?,"application/json","read-check-replay.json")?;
            output(json!({"workspace":args.workspace,"attachment_hash":args.attachment_hash,"context_revision":checked.context_revision,
                "reference":checked.saved.reference,"reads":checked.saved.reads,"state":"fresh-computational-evidence",
                "engine_version":checked.report.engine_version,"status":checked.report.status,
                "eligible_coverage":checked.report.eligible_coverage,"reference_coverage_fraction":checked.report.reference_coverage_fraction,
                "artifact_id":artifact,"verified":false,"authority":"Fresh computational evidence under recorded rules; not experimental sample confirmation."}))
        }).await
    }
    #[tool(
        description = "Preview explicit decisions on a saved construct read check after fresh replay. Nominated bases, calls, deletion evidence and insertion read/boundary/sequence must match current eligible evidence exactly; index-only or unsupported edits are refused. Retain/unresolved decisions do not edit DNA. Original reads/reference remain unchanged; any revised reference is a proposal requiring a new read check, not measured sequence or experimental confirmation."
    )]
    async fn preview_read_review(
        &self,
        Parameters(args): Parameters<ReadReviewQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let frame=server.restore_workspace(&args.workspace)?;
            let preview=collagene_core::read_review::preview(&frame,&args.request).map_err(bad)?;
            let expected=collagene_core::read_review::preview_revision(&frame,&args.request).map_err(bad)?;
            preview.validate_binding(&expected).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&preview,32*1024*1024)?,"application/json","read-review-preview.json")?;
            let target=preview.target.as_ref().map(|m|json!({"name":m.name,"length":m.sequence.len(),"topology":m.topology,"revision":m.revision()}));
            output(json!({"state":"proposal","verified":false,"preview":{"context_revision":preview.context_revision,"request":preview.request,
                "reference":preview.reference,"reads":preview.reads,"target":target,"edits":preview.edits,"warnings":preview.warnings},"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Freshly replay a saved read check, validate explicit review decisions and atomically retain the decisions and optional derived reference in a NEW immutable native workspace. Exact original reference/read events become computational parents. Source heads and raw calls remain unchanged. Human choices never establish a measurement or current positive read verdict; a derived reference requires a new independent read check. Retain/unresolved-only requests save decisions without creating a reference."
    )]
    async fn commit_read_review(
        &self,
        Parameters(args): Parameters<ReadReviewQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let frame=server.restore_workspace(&args.workspace)?;
            let expected=collagene_core::read_review::commit_revision(&frame,&args.request).map_err(bad)?;
            let commit=collagene_core::read_review::commit(&frame,&args.request).map_err(bad)?;
            commit.validate_binding(&expected).map_err(bad)?;
            let saved=collagene_core::read_review::load(&commit.workspace,&commit.attachment_hash).map_err(bad)?;
            let text=bounded_json_limit(&saved,32*1024*1024)?;
            let receipt=json!({"attachment_hash":commit.attachment_hash,"context_revision":commit.context_revision,"target":commit.target,"state":"historical","verified":false});
            let native=workspace::save_workspace(&commit.workspace).map_err(bad)?;
            server.commit_workspace_outputs_named(commit.workspace,Some((text,receipt)),Some(native),"read-review.json",
                "Recorded caller decisions and optional proposed reference only. Unselected differences remain unreviewed; source calls are unchanged and the revised reference requires a new read check.")
        }).await
    }
    #[tool(
        description = "List immutable saved decisions on construct read checks and their optional derived-reference events. Historical review records do not restore current validation, experimental confirmation or an approved whole-construct sequence."
    )]
    async fn list_read_reviews(
        &self,
        Parameters(args): Parameters<ListReadChecks>,
    ) -> Result<CallToolResult, ErrorData> {
        let frame = self.restore_workspace(&args.workspace)?;
        let reviews = collagene_core::read_review::list(&frame).map_err(bad)?;
        output(
            json!({"workspace":args.workspace,"reviews":reviews,"state":"historical","verified":false}),
        )
    }
    #[tool(
        description = "Reconstruct an exact session plan against explicit historical target/input event bindings and retain all products in a NEW immutable session workspace. Does not alter the original, browser workspace or any file. Incomplete evidence requires explicit allow_incomplete; failed reconstructions never record products. Returns a full JSON artifact and new workspace reference."
    )]
    async fn record_workspace_run(
        &self,
        Parameters(args): Parameters<RecordWorkspaceRun>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let value = server.restore_workspace(&args.workspace)?;
            let plan = server
                .store
                .lock()
                .unwrap()
                .plans
                .get(&args.plan_revision)
                .ok_or_else(|| bad("unknown plan revision"))?
                .plan
                .clone();
            let commit = collagene_core::run_lineage::record_run(
                &value,
                &plan,
                args.target_binding,
                args.input_bindings,
                args.allow_incomplete,
            )
            .map_err(bad)?;
            let text = bounded_json(&commit)?;
            let receipt = serde_json::to_value(&commit.receipt).map_err(bad)?;
            server.commit_workspace_artifact(commit.workspace, Some((text, receipt)))
        })
        .await
    }
    #[tool(
        description = "List historical cloning runs and exact source events in an immutable workspace. Previews cap100 runs/events and20 products per run; full JSON artifact retains every binding. Stored runs never establish current verification authority."
    )]
    async fn list_workspace_runs(
        &self,
        Parameters(args): Parameters<WorkspaceRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let value = server.restore_workspace(&args)?;
            let runs = collagene_core::run_lineage::list_runs(&value).map_err(bad)?;
            let full = json!({"runs":runs,"events":value.history});
            let artifact = server.analysis_artifact(&full)?;
            let previews = runs.iter().take(100).map(|run|json!({"run_id":run.run_id,"attachment_hash":run.attachment_hash,"plan_name":run.plan_name.chars().take(256).collect::<String>(),"target":run.target,"product_count":run.products.len(),"products":run.products.iter().take(20).collect::<Vec<_>>(),"products_truncated":run.products.len()>20,"recorded_status":run.recorded_status,"recorded_engine":run.recorded_engine,"input_count":run.input_count})).collect::<Vec<_>>();
            let events = value.history.iter().take(100).map(|event|json!({"document_id":event.document_id,"event_id":event.id,"revision":event.revision,"operation":event.operation.chars().take(256).collect::<String>()})).collect::<Vec<_>>();
            output(json!({"workspace":args,"run_count":runs.len(),"runs":previews,"runs_truncated":runs.len()>100,"event_count":value.history.len(),"events":events,"events_truncated":value.history.len()>100,"artifact_id":artifact,"authority":"historical data; replay recomputes current evidence"}))
        }).await
    }
    #[tool(
        description = "Prepare one exact saved native run for current review using its immutable workspace/run ID. Resolves original input and independent-target events even when current heads changed. Atomically returns full documentary JSON plus a binary candidate artifact and unverified session plan usable with verify_plan. Does not replay/adopt historical status or material, replace current targets, alter workspace heads, or write a native database."
    )]
    async fn prepare_saved_run_review(
        &self,
        Parameters(args): Parameters<ReplayWorkspaceRun>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.run_id.is_empty() || args.run_id.len() > 128 {
            return Err(bad("explicit bounded native run ID required"));
        }
        let server = self.clone();
        self.compute(move || {
            let workspace = server.restore_workspace(&args.workspace)?;
            let expected = collagene_core::run_review::request_revision(&workspace, &args.run_id)
                .map_err(bad)?;
            let review =
                collagene_core::run_review::prepare(&workspace, &args.run_id).map_err(bad)?;
            review.validate_request_binding(&expected).map_err(bad)?;
            server.saved_review_artifacts(review, args.workspace)
        })
        .await
    }
    #[tool(
        description = "Export the canonical feature map of one explicit immutable session document revision as bounded native SVG. Returns exact source identity, renderer warnings and a complete image/svg+xml artifact. No arbitrary paths, external renderer or user-supplied SVG. This profile requests no restriction-enzyme overlay; map visibility is presentation only and never physical/material or verification authority."
    )]
    async fn export_map_svg(
        &self,
        Parameters(reference): Parameters<DocumentRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source={let store=server.store.lock().unwrap();let source=store.document(&reference)?;collagene_core::model::bounds::model_preflight(&[source]).map_err(bad)?;source.clone()};
            let rendered=collagene_core::svg_export::render(&collagene_core::svg_export::Request{source,source_revision:reference.revision.clone(),restriction_sites:None}).map_err(bad)?;
            rendered.validate_binding(&reference.revision).map_err(bad)?;
            let size=rendered.svg.len();let artifact=blake3::hash(rendered.svg.as_bytes()).to_hex().to_string();
            let result=output(json!({"document":reference,"source_revision":rendered.source_revision,"artifact_id":artifact,"filename":"map.svg","encoding":"utf-8","mime_type":"image/svg+xml","bytes":size,"warnings":rendered.warnings,"warning_count":rendered.warnings.len(),"warnings_truncated":false,"restriction_overlay":"none_requested","scope":"Canonical features of this exact native revision; no restriction-site overlay requested.","authority":"Documentary display only. Map export does not verify a cloning design, prove physical material, alter a source or write a native catalog."}))?;
            let mut store=server.store.lock().unwrap();
            if let Some(existing)=store.artifacts.get(&artifact){
                if !matches!(&existing.data,ArtifactData::Text(_))||existing.mime_type!="image/svg+xml"||existing.filename!="map.svg"{return Err(bad("map SVG artifact already has another representation"))}
            }else{
                if store.artifacts.len()>=MAX_ITEMS{return Err(bad("at most 128 artifacts; no map published"))}
                store.reserve(size)?;store.artifacts.insert(artifact,Artifact{data:ArtifactData::Text(rendered.svg),mime_type:"image/svg+xml",encoding:"utf-8",filename:"map.svg"});
            }Ok(result)
        }).await
    }
    #[tool(
        description = "Replay a retained run from exact historical sources using the current embedded verifier. Returns bounded product previews and a complete JSON evidence artifact. Engine/status/product comparisons are explicit. Does not edit a workspace, trust a stored report or invoke external programs."
    )]
    async fn replay_workspace_run(
        &self,
        Parameters(args): Parameters<ReplayWorkspaceRun>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let value = server.restore_workspace(&args.workspace)?;
            let report = collagene_core::run_lineage::replay_run(&value, &args.run_id).map_err(bad)?;
            report.validate_for_workspace(&value).map_err(bad)?;
            let full = serde_json::to_value(&report).map_err(bad)?;
            let artifact = server.analysis_artifact(&full)?;
            let products = report.current_report.products.iter().take(100).map(|p| json!({"revision":p.revision(),"name":p.name.chars().take(256).collect::<String>(),"length":p.sequence.len()})).collect::<Vec<_>>();
            output(json!({"workspace":args.workspace,"run_id":report.run_id,"source_workspace_fingerprint":report.source_workspace_fingerprint,"recorded_status":report.recorded_status,"recorded_engine":report.recorded_engine,"current_status":report.current_report.status,"current_engine":report.current_report.engine_version,"engine_matches":report.engine_matches,"products_match":report.products_match,"status_matches":report.status_matches,"products":products,"product_count":report.current_report.products.len(),"products_truncated":report.current_report.products.len()>100,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Read immutable historical reference-design and processed physical graph attachments for a retained run. Checks exact run/plan/engine/source/target/product-event bindings but does not replay chemistry or restore physical authority. Older runs may have no historical graph; fresh replay is a separate tool."
    )]
    async fn inspect_workspace_run_evidence(
        &self,
        Parameters(args): Parameters<ReplayWorkspaceRun>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let workspace=server.restore_workspace(&args.workspace)?;
            let record=collagene_core::run_lineage::load_run(&workspace,&args.run_id).map_err(bad)?;
            let evidence=collagene_core::run_evidence::download(&workspace,&args.run_id).map_err(bad)?;
            let artifact=evidence.bytes.as_ref().map(|bytes|server.text_artifact(String::from_utf8(bytes.clone()).map_err(bad)?,"application/json","historical-cloning-graphs.json")).transpose()?;
            output(json!({"state":"historical","verified":false,"run_id":record.run_id,"recorded_engine":record.assessment.engine_version,"scientific_step_count":evidence.step_count,"artifact_id":artifact,"authority":"Stored historical graph evidence only. Fresh replay of exact original inputs is required for current material state or verification."}))
        }).await
    }
    #[tool(
        description = "Reserve a session upload for exact native workspace bytes. Length 1..64 MiB, lower-case BLAKE3 required. At most four uploads; full expected length is charged to the shared 64 MiB session budget immediately. No filesystem or external programs."
    )]
    fn begin_workspace_upload(
        &self,
        Parameters(args): Parameters<BeginWorkspaceUpload>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.expected_bytes == 0
            || args.expected_bytes > workspace::MAX_ARCHIVE_BYTES
            || args.blake3.len() != 64
            || !args
                .blake3
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(bad("invalid archive length or lower-case BLAKE3"));
        }
        let mut store = self.store.lock().unwrap();
        if store.uploads.len() >= MAX_UPLOADS {
            return Err(bad("at most four outstanding workspace uploads"));
        }
        let serial = store
            .upload_serial
            .checked_add(1)
            .ok_or_else(|| bad("upload counter exhausted"))?;
        store.reserve(args.expected_bytes)?;
        let mut bytes = Vec::new();
        if let Err(error) = bytes.try_reserve_exact(args.expected_bytes) {
            store.used -= args.expected_bytes;
            return Err(bad(format!("upload allocation failed: {error}")));
        }
        store.upload_serial = serial;
        let id = format!("workspace-upload-{serial}");
        store.uploads.insert(
            id.clone(),
            WorkspaceUpload {
                expected: args.expected_bytes,
                hash: args.blake3,
                bytes,
            },
        );
        output(
            json!({"upload_id":id,"next_offset":0,"expected_bytes":args.expected_bytes,"max_decoded_chunk_bytes":MAX_UPLOAD_CHUNK}),
        )
    }
    #[tool(
        description = "Append one independently base64-encoded native workspace chunk, <=256 KiB decoded. Offset must equal the exact next decoded offset; gaps/overlaps/replays/overflow are rejected without changing the upload. No raw archive is parsed before finalize."
    )]
    fn append_workspace_upload(
        &self,
        Parameters(args): Parameters<AppendWorkspaceUpload>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > MAX_UPLOAD_CHUNK.div_ceil(3) * 4 {
            return Err(bad("upload chunk exceeds 256 KiB decoded"));
        }
        let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
        if bytes.is_empty() || bytes.len() > MAX_UPLOAD_CHUNK {
            return Err(bad("upload chunk must contain 1..256 KiB decoded bytes"));
        }
        let mut store = self.store.lock().unwrap();
        let upload = store
            .uploads
            .get_mut(&args.upload_id)
            .ok_or_else(|| bad("unknown upload"))?;
        if args.offset != upload.bytes.len()
            || bytes.len() > upload.expected.saturating_sub(upload.bytes.len())
        {
            return Err(bad(
                "incorrect next offset or chunk exceeds expected archive length",
            ));
        }
        upload.bytes.extend_from_slice(&bytes);
        output(
            json!({"upload_id":args.upload_id,"next_offset":upload.bytes.len(),"expected_bytes":upload.expected,"complete":upload.bytes.len()==upload.expected}),
        )
    }
    #[tool(
        description = "Finalize an exactly complete native workspace upload. Verifies length and BLAKE3 before parsing, then imports all revisions/history/candidates/attachments atomically under core64MiB archive/128MiB expanded and session64MiB limits. Incomplete uploads remain appendable; completed attempts consume/release upload storage even on hash/parse/storage failure. Imported candidates are never authoritative verification."
    )]
    async fn finalize_workspace_upload(
        &self,
        Parameters(args): Parameters<WorkspaceUploadRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let upload = {
                let mut store = server.store.lock().unwrap();
                let upload = store
                    .uploads
                    .get(&args.upload_id)
                    .ok_or_else(|| bad("unknown upload"))?;
                if upload.bytes.len() != upload.expected {
                    return Err(bad("upload incomplete; append remaining bytes or abort"));
                }
                store.uploads.remove(&args.upload_id).unwrap()
            };
            let decoded = if blake3::hash(&upload.bytes).to_hex().as_str() != upload.hash {
                Err(bad("archive BLAKE3 mismatch; upload consumed"))
            } else {
                workspace::load_workspace(&upload.bytes).map_err(bad)
            };
            // The reservation remains charged throughout hashing/decompression. Imports
            // are committed only after all core validation, with no partial documents.
            server.store.lock().unwrap().used -= upload.expected;
            drop(upload);
            server.commit_workspace(decoded?)
        })
        .await
    }
    #[tool(
        description = "Abort a pending workspace upload and release its full reserved session bytes. No documents or candidates are imported by abort."
    )]
    fn abort_workspace_upload(
        &self,
        Parameters(args): Parameters<WorkspaceUploadRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let mut store = self.store.lock().unwrap();
        let upload = store
            .uploads
            .remove(&args.upload_id)
            .ok_or_else(|| bad("unknown upload"))?;
        store.used -= upload.expected;
        output(json!({"upload_id":args.upload_id,"aborted":true,"released_bytes":upload.expected}))
    }
    #[tool(
        description = "List immutable imported workspace summaries with unverified candidate metadata. Exact workspace id/revision is required to export; molecule revisions are shared and accessible through paginated list_documents."
    )]
    fn list_workspaces(
        &self,
        Parameters(_args): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.store.lock().unwrap();
        output(
            json!({"workspaces":store.workspaces.iter().map(|(id,w)|w.summary(id)).collect::<Vec<_>>()}),
        )
    }
    #[tool(
        description = "Reconstruct explicit single-nucleotide tailing of BOTH DNA strand3-prime ends from a complete digest (optionally blunted) or a unique exact-PCR product. Requires actual blunt linear duplex and accessibleOH; retains five-prime chemistry and nicks. This idealendpoint transform does not predict enzyme/kit activity and cannot establish physical authority when imported as metadata. Full immutable source/recipe/strand report artifact included; use explicit Tail3Prime in a targetbound plan."
    )]
    async fn simulate_tailing(
        &self,
        Parameters(args): Parameters<TailingQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.store.lock().unwrap().document(&args.document)?.clone();
            let preview=verification::simulate_tailing(&source,&args.recipe,args.base).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"tailing":preview}))?;
            output(json!({"document":args.document,"source_revision":preview.source_revision,"recipe":preview.recipe,"base":preview.base,"engine_version":verification::ENGINE_VERSION,"top_length":preview.tailing.duplex.top.len(),"bottom_length":preview.tailing.duplex.bottom.len(),"bottom_start":preview.tailing.duplex.bottom_start,"additions":preview.tailing.additions,"warnings":preview.tailing.warnings,"artifact_id":artifact,"authority":"Explicit ideal tailing evidence only. Imported preview metadata cannot establish physical end state, template ancestry, enzyme activity or target verification."}))
        }).await
    }
    #[tool(
        description = "Check explicitly oriented ordered DNA overlap junctions against the declared Gibson/HiFi/In-Fusion-style profile and bounded repeat evidence. Source orientation, exact homology and length/GC findings are retained. Vendor recommendations are advisory: this neither simulates proprietary reagents nor establishes experimental compatibility, uniqueness or physical substrate chemistry. Full source/options-bound JSON artifact included; reconstruct ProfiledOverlap against an independent target to verify the candidate."
    )]
    async fn check_overlap_profile(
        &self,
        Parameters(args): Parameters<OverlapProfileQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let sources={let store=server.store.lock().unwrap();args.documents.iter().map(|r|store.document(r).cloned()).collect::<Result<Vec<_>,_>>()?};
            let report=collagene_core::overlap_profiles::check_oriented_fragments(&sources,&args.reverse,&args.overlaps,args.circular,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"documents":args.documents,"engine_version":verification::ENGINE_VERSION,"overlap_profile":report}))?;
            let mut preview=report.clone();preview.junctions.truncate(50);for junction in &mut preview.junctions{junction.occurrences.truncate(20);if junction.homology.len()>600{junction.homology.truncate(600);}}
            output(json!({"documents":args.documents,"reverse":args.reverse,"engine_version":verification::ENGINE_VERSION,"geometry_valid":report.geometry_valid,"recommendation_compliant":report.recommendation_compliant,"complete":report.complete,"exhaustive_alternative_search":false,"comparisons":report.comparisons,"junction_count":report.junctions.len(),"display_truncated":report.junctions.len()>50||report.junctions.iter().any(|j|j.occurrences.len()>20||j.homology.len()>600),"junctions":preview.junctions,"findings":report.findings.iter().take(100).collect::<Vec<_>>(),"warnings":report.warnings.iter().take(100).collect::<Vec<_>>(),"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Reconstruct a complete restriction digest and selected fragment, then apply explicit 3-prime fill-in/exonuclease end repair. Preserves 5-prime chemistry and internal nicks. Imported duplex metadata never grants physical authority. Returns bounded sequence/end summaries and full source-bound report artifact; use Digest then RepairEnds in a separately target-bound plan for verification."
    )]
    async fn simulate_end_repair(
        &self,
        Parameters(args): Parameters<EndRepairQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.store.lock().unwrap().document(&args.document)?.clone();
            let preview=verification::simulate_end_repair(&source,&args.enzymes,args.fragment_index,args.method).map_err(bad)?;
            let summary=preview.repair.duplex.summary().map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"repair":preview}))?;
            output(json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"length":preview.molecule.sequence.len(),"sequence_preview":preview.molecule.sequence.as_str().chars().take(200).collect::<String>(),"sequence_preview_truncated":preview.molecule.sequence.len()>200,"summary":summary,"top_added_3prime":preview.repair.top_added_3prime.chars().take(200).collect::<String>(),"top_added_count":preview.repair.top_added_3prime.len(),"top_removed_count":preview.repair.top_removed_3prime.len(),"bottom_added_count":preview.repair.bottom_added_3prime.len(),"bottom_removed_count":preview.repair.bottom_removed_3prime.len(),"annotation_count":preview.molecule.features.len(),"findings":preview.findings.iter().take(100).collect::<Vec<_>>(),"artifact_id":artifact,"authority":"Read-only preview. Physical intermediates must be reconstructed in the bound plan; this artifact is not an imported end-chemistry authority."}))
        }).await
    }
    #[tool(
        description = "Enumerate bounded exact antiparallel placements of two explicitly declared linear 5-prime-to-3-prime ssDNA oligos with supplied terminal chemistry. Does not assume automatic longest/unique pairing, phosphorylation, ligation or experimental success. Returns alternatives, completeness and a full immutable-source-bound artifact. Use an explicit AnnealOligos candidate in a separately bound plan and reverify."
    )]
    async fn simulate_oligo_annealing(
        &self,
        Parameters(args): Parameters<AnnealQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let(top,bottom)={let store=server.store.lock().unwrap();(store.document(&args.top)?.clone(),store.document(&args.bottom)?.clone())};
            let preview=verification::simulate_oligo_annealing(&top,&bottom,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"top":args.top,"bottom":args.bottom,"engine_version":verification::ENGINE_VERSION,"annealing":preview}))?;
            output(json!({"top":args.top,"bottom":args.bottom,"engine_version":verification::ENGINE_VERSION,"candidate_count":preview.annealing.candidates.len(),"complete":preview.annealing.complete,"unique":preview.annealing.unique,"comparisons":preview.annealing.comparisons,"options":preview.annealing.options,"warnings":preview.annealing.warnings,"display_truncated":preview.annealing.candidates.len()>50,"candidates":preview.annealing.candidates.iter().enumerate().take(50).map(|(index,c)|json!({"candidate_index":index,"id":c.id,"bottom_start":c.duplex.bottom_start,"top_length":c.duplex.top.len(),"bottom_length":c.duplex.bottom.len(),"paired_start":c.paired_start,"paired_end":c.paired_end,"paired_bases":c.paired_bases})).collect::<Vec<_>>(),"artifact_id":artifact,"authority":"Terminal chemistry is an explicit caller assumption. Base pairing creates no covalent bonds and imported metadata cannot establish a duplex state."}))
        }).await
    }
    #[tool(
        description = "Enumerate finite one-copy Golden Gate fragment assemblies using actual tracked-position cleavage checks. These bare originals have unknown modification state; unresolved cleavage withholds the requested pool. Use inspect_golden_gate_reaction for explicit original-generation or saved-material/PCR declarations. Distinguishes stable site-free candidates, recut paths, uncut carryover, source histories and budget completeness. This is not an exhaustive kinetic reaction-mixture search. Full JSON evidence artifact included; candidate index can be used in a separately target-bound GoldenGate plan step."
    )]
    async fn simulate_golden_gate(
        &self,
        Parameters(args): Parameters<GoldenGateQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let sources={let store=server.store.lock().unwrap();args.documents.iter().map(|r|store.document(r).cloned()).collect::<Result<Vec<_>,_>>()?};
            let report=collagene_core::golden_gate::simulate(&sources,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"documents":args.documents,"engine_version":verification::ENGINE_VERSION,"golden_gate":report}))?;
            output(json!({"documents":args.documents,"candidate_count":report.candidates.len(),"profile_revision":report.profile_revision,"products_determined":report.products_determined,"cleavage_comparisons":report.comparisons,"stable_physical_product_count":report.stable_physical_product_count,"enumeration_complete":report.enumeration_complete,"exhaustive_mixture_search":report.exhaustive_mixture_search,"states_examined":report.states_examined,"carryover_sources":report.carryover_sources,"truncation_reasons":report.truncation_reasons,"warnings":report.warnings,"display_truncated":report.candidates.len()>50,"candidates":report.candidates.iter().enumerate().take(50).map(|(i,c)|json!({"candidate_index":i,"length":c.duplex.top.len(),"topology":c.duplex.topology,"path":c.path,"restriction_site_free":c.restriction_site_free,"recut_site_count":c.recut_sites.len(),"recut_complete":c.cleavage.complete,"recut_products_determined":c.cleavage.products_determined,"directed_sequence_key":c.directed_sequence_key,"physical_sequence_key":c.physical_sequence_key})).collect::<Vec<_>>(),"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Reconstruct intentional-mismatch PCR from immutable template and actual tailed oligos under an explicit mismatch/3-prime/work policy. Returns introduced changes, competing primer roles, overlap/ambiguity diagnostics and search completeness. Overlapping-domain/uncertain products remain previews. Full report artifact is computed, not an AI assertion."
    )]
    async fn simulate_mismatch_pcr(
        &self,
        Parameters(args): Parameters<MismatchPcrQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.store.lock().unwrap().document(&args.document)?.clone();
            let report=collagene_core::mutagenesis::simulate_pcr(&source,&args.forward,&args.reverse,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"forward":args.forward,"reverse":args.reverse,"pcr":report}))?;
            output(json!({"document":args.document,"product_count":report.products.len(),"search_complete":report.search_complete,"enumeration_complete":report.enumeration_complete,"unsupported_geometries":report.unsupported_geometries,"conflicting_overlaps":report.conflicting_overlaps,"warnings":report.warnings,"display_truncated":report.products.len()>100,"products":report.products.iter().enumerate().take(100).map(|(i,p)|json!({"product_index":i,"length":p.sequence.len(),"forward_start":p.forward_start,"reverse_start":p.reverse_start,"forward_primer_index":p.forward_primer_index,"reverse_primer_index":p.reverse_primer_index,"introduced_changes":p.substitutions.iter().take(20).collect::<Vec<_>>(),"change_count":p.substitutions.len(),"certain":p.certain,"overlapping_domains":p.overlapping_domains})).collect::<Vec<_>>(),"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Design inverse-PCR replacement/insertion/deletion primers from an explicit circular source and independently supplied circular target. origin and origin_workspace are required nullable fields: null/null declares an original reference hypothesis; a saved Run/Derived selection requires its exact immutable native workspace and fresh source-event attestation. Reconstructs the origin when present, then actual PCR, phosphorylation and circular ligation. Nicked/unsupported templates refuse without repair or metadata authority. Atomically stores an unverified candidate and full source-bound artifact; inherited incomplete findings remain explicit. Verify the returned plan_revision before export/use; no persistent project writes or experimental success are implied."
    )]
    async fn design_inverse_pcr(
        &self,
        Parameters(args): Parameters<InversePcrQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            if args.origin.is_some()!=args.origin_workspace.is_some(){return Err(bad("Saved inverse-PCR origin and exact native workspace must be provided together."));}
            let frame=args.origin_workspace.as_ref().map(|w|server.restore_workspace(w)).transpose()?;
            let origin=match (&args.origin,&frame){
                (Some(selection),Some(frame))=>Some(server.resolve_material_origin(frame,selection.event(),selection)?),
                (None,None)=>None,
                _=>return Err(bad("Incomplete inverse-PCR origin binding.")),
            };
            let (source,target,expected)={
                let store=server.store.lock().unwrap();let source=store.document(&args.document)?;let target=store.document(&args.target)?;
                if let (Some(frame),Some(origin))=(&frame,&origin){if frame.revisions.get(&origin.event().revision)!=Some(source){return Err(bad("Inverse-PCR source differs from the exact retained material event."));}}
                let expected=collagene_core::mutagenesis::inverse_pcr_revision(source,target,origin.as_ref(),&args.options).map_err(bad)?;
                (source.clone(),target.clone(),expected)
            };
            let design=collagene_core::mutagenesis::design_inverse_pcr(&source,&target,origin.as_ref(),&args.options).map_err(bad)?;
            design.validate_binding(&expected).map_err(bad)?;
            transfer::Candidate::new(design.plan.clone()).validate().map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a>{document:&'a DocumentRef,target:&'a DocumentRef,origin:&'a Option<collagene_core::run_material::MaterialSelection>,origin_workspace:&'a Option<WorkspaceRef>,source_event:Option<&'a collagene_core::run_lineage::EventBinding>,engine_version:&'static str,inverse_pcr:&'a collagene_core::mutagenesis::InversePcrDesign}
            let text=bounded_json_limit(&Evidence{document:&args.document,target:&args.target,origin:&args.origin,origin_workspace:&args.origin_workspace,source_event:origin.as_ref().map(|o|o.event()),engine_version:verification::ENGINE_VERSION,inverse_pcr:&design},32*1024*1024)?;
            let artifact_id=blake3::hash(text.as_bytes()).to_hex().to_string();let plan_revision=verification::plan_revision(&design.plan);
            let plan_bytes=bounded_json_limit(&design.plan,32*1024*1024)?.len().checked_add(bounded_json_limit(&target,32*1024*1024)?.len()).ok_or_else(||bad("Inverse-PCR storage accounting overflow."))?;
            let result=output(json!({"document":args.document,"target":args.target,"origin":args.origin,"origin_workspace":args.origin_workspace,"source_event":origin.as_ref().map(|o|o.event()),"context_revision":expected,"plan_revision":plan_revision,"verified":false,"forward":design.forward,"reverse":design.reverse,"source_revision":design.source_revision,"target_revision":design.target_revision,"computed_status":design.verification.status,"findings":design.verification.findings,"warnings":design.warnings,"artifact_id":artifact_id}))?;
            let mut store=server.store.lock().unwrap();let new_plan=!store.plans.contains_key(&plan_revision);let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(existing)=store.artifacts.get(&artifact_id){if !matches!(existing.data,ArtifactData::Text(_))||existing.filename!="sequence-analysis.json"||existing.mime_type!="application/json"{return Err(bad("Inverse-PCR artifact has another representation; no candidate was published."));}}
            if (new_plan&&store.plans.len()>=MAX_ITEMS)||(new_artifact&&store.artifacts.len()>=MAX_ITEMS){return Err(bad("Inverse-PCR candidate/artifact capacity reached; nothing was published."));}
            let total=(if new_plan{plan_bytes}else{0}).checked_add(if new_artifact{text.len()}else{0}).ok_or_else(||bad("Inverse-PCR storage accounting overflow."))?;store.reserve(total)?;
            if new_plan{store.plans.insert(plan_revision,StoredPlan{plan:design.plan,target});}
            if new_artifact{store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"sequence-analysis.json"});}
            Ok(result)
        }).await
    }
    #[tool(
        description = "Predict a reference gel from revision-bound uncut, complete-digest or PCR lanes under an explicit ladder/calibration. CatalogDigest lanes require a native workspace reference proving the recipe's exact retained document/event/revision and a matching session document; historical sources need not be current heads. Native inspection is recomputed. Circular/unsupported profiles remain unpositioned. JSON retains declared lane recipes and references; SVG remains a prediction. No observed-gel/OpenGel integration."
    )]
    async fn predict_gel(
        &self,
        Parameters(args): Parameters<GelQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            if args.lanes.is_empty() || args.lanes.len() > 32 {
                return Err(bad("gel needs 1–32 source lanes"));
            }
            args.calibration.validate().map_err(bad)?;
            for lane in &args.lanes {
                if lane.label.len() > 256 || lane.label.chars().any(char::is_control)
                    || lane.mass_ng.is_some_and(|mass| !mass.is_finite() || mass <= 0. || mass > 1e9)
                {
                    return Err(bad("invalid gel lane label or loaded mass"));
                }
                if let collagene_core::gel::Preparation::CatalogDigest { recipe } = &lane.preparation {
                    recipe.validate().map_err(bad)?;
                    let reference = lane.workspace.as_ref().ok_or_else(|| bad(
                        "CatalogDigest requires a native workspace reference proving its exact retained source event"
                    ))?;
                    if reference.id.len() != 64 || reference.revision.len() != 64 {
                        return Err(bad("invalid pinned gel workspace reference"));
                    }
                }
            }
            // Resolve immutable session documents together, then release this
            // lock before restore_workspace acquires the same store mutex.
            let lanes = {
                let store = server.store.lock().unwrap();
                args.lanes.iter().map(|lane| Ok(collagene_core::gel::LaneInput {
                    label: lane.label.clone(),
                    molecule: store.document(&lane.document)?.clone(),
                    preparation: lane.preparation.clone(),
                    mass_ng: lane.mass_ng,
                })).collect::<Result<Vec<_>, ErrorData>>()?
            };
            let mut groups = std::collections::BTreeMap::<String, (WorkspaceRef, Vec<usize>)>::new();
            for (index, lane) in args.lanes.iter().enumerate() {
                if !matches!(&lane.preparation, collagene_core::gel::Preparation::CatalogDigest { .. }) { continue; }
                let reference = lane.workspace.as_ref().ok_or_else(|| bad("pinned gel workspace reference is unavailable"))?;
                if reference.id != reference.revision { return Err(bad("stale workspace revision")); }
                groups.entry(reference.id.clone()).or_insert_with(|| (reference.clone(), Vec::new())).1.push(index);
            }
            // Restore each distinct workspace once, then drop it before loading
            // another: shared large histories are never cloned once per lane.
            for (_, (reference, indices)) in groups {
                let frame = server.restore_workspace(&reference)?;
                for index in indices {
                    let lane = &lanes[index];
                    let collagene_core::gel::Preparation::CatalogDigest { recipe } = &lane.preparation else { continue; };
                    let source = &recipe.source;
                    if !frame.documents.iter().any(|document| document.id == source.document_id) {
                        return Err(bad("pinned gel source owner is not retained in this workspace"));
                    }
                    let event = frame.history.iter().find(|event| event.id == source.event_id)
                        .ok_or_else(|| bad("the exact pinned gel source event is not retained"))?;
                    if event.document_id != source.document_id || event.revision != source.revision {
                        return Err(bad("pinned gel source owner/event/revision mismatch"));
                    }
                    let original = frame.revisions.get(&source.revision)
                        .ok_or_else(|| bad("pinned gel historical molecule is not retained"))?;
                    if original != &lane.molecule || lane.molecule.revision() != source.revision {
                        return Err(bad("session lane document differs from the exact historical pinned gel source"));
                    }
                }
            }
            let report=collagene_core::gel::predict_gel(&lanes,&args.calibration).map_err(bad)?;
            let svg=collagene_core::gel::render_svg(&report).map_err(bad)?;
            let svg_artifact=server.text_artifact(svg,"image/svg+xml","predicted-gel.svg")?;
            let artifact=server.analysis_artifact(&json!({"documents":args.lanes.iter().map(|l|&l.document).collect::<Vec<_>>(),"lanes":args.lanes,"engine_version":verification::ENGINE_VERSION,"gel":report}))?;
            output(json!({"model":report.model,"lane_count":report.lanes.len(),"lanes":report.lanes.iter().map(|l|json!({"label":l.label,"source_revision":l.source_revision,"fragment_count":l.fragments.len(),"resolved_band_count":l.bands.len(),"unpositioned_fragment_count":l.unpositioned_fragments.len(),"complete_under_model":l.complete_under_model,"error":l.error,"warnings":l.warnings})).collect::<Vec<_>>(),"artifact_id":artifact,"svg_artifact_id":svg_artifact}))
        }).await
    }

    #[tool(
        description = "Rank bounded single/double restriction digests that distinguish immutable candidate plasmids under an explicit gel calibration/resolution. Reports distinguished, co-migrating and unknown pairs separately, with evaluated/work bounds and exhaustive-search flag. Unsupported profiles do not become negative specificity evidence. Full JSON evidence artifact included."
    )]
    async fn design_diagnostic_digest(
        &self,
        Parameters(args): Parameters<DiagnosticQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let sources={let store=server.store.lock().unwrap();args.documents.iter().map(|r|store.document(r).cloned()).collect::<Result<Vec<_>,_>>()?};
            let report=collagene_core::gel::diagnostic_digests(&sources,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"documents":args.documents,"engine_version":verification::ENGINE_VERSION,"diagnostic_digests":report}))?;
            output(json!({"documents":args.documents,"source_revisions":report.source_revisions,"evaluated_digests":report.evaluated_digests,"work_bases":report.work_bases,"exhaustive_search":report.exhaustive_search,"warnings":report.warnings,"display_truncated":report.candidates.len()>20,"candidates":report.candidates.iter().take(20).map(|c|json!({"enzymes":c.enzymes,"distinguished_pairs":c.distinguished_pairs,"indistinguishable_pairs":c.indistinguishable_pairs,"unknown_pairs":c.unknown_pairs,"minimum_diagnostic_separation":c.minimum_diagnostic_separation})).collect::<Vec<_>>(),"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Resolve a saved primer_bind annotation at its exact immutable document revision and zero-based feature index. Uses shared strict oligo/annealing/5-prime-tail consistency checks; missing fields may be inferred only from exact supported strand-oriented geometry, never fuzzy/unresolved/reviewed geometry. Returns the actual ordered Primer, sequence origin and warnings, preserving deliberate mismatches/IUPAC and historical source qualifiers as evidence. No binding, PCR, thermodynamic or experimental proof; pass the returned primer explicitly to simulate_pcr/analyze_primer."
    )]
    async fn resolve_primer(
        &self,
        Parameters(args): Parameters<SavedPrimerQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let resolved={let store=server.store.lock().unwrap();let source=store.document(&args.document)?;
                collagene_core::primer_reuse::resolve_at_revision(source,args.feature_index,&args.document.revision).map_err(bad)?};
            output(json!({"document":args.document,"resolved":resolved,
                "authority":"Saved oligo identity only. Current template binding and exact PCR products must be recomputed independently; source location and historical analysis are not amplification evidence."}))
        }).await
    }
    #[tool(
        description = "Analyze a tailed DNA primer under explicit concentration/sodium conditions, screen contiguous self/partner complementarity, and find bounded binding sites on an immutable template. Certain 3-prime matches, mismatches and ambiguous compatibility remain separate. This is not Primer3, empirical structure prediction or genome-wide specificity. Full JSON evidence artifact included."
    )]
    async fn analyze_primer(
        &self,
        Parameters(args): Parameters<PrimerQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.store.lock().unwrap().document(&args.document)?.clone();
            let analysis=primers::analyze_primer(&args.primer,&args.conditions).map_err(bad)?;
            let bindings=primers::find_binding_sites(&source,&args.primer,&args.binding_options).map_err(bad)?;
            let heterodimer=args.partner.as_ref().map(|p|primers::heterodimer(&args.primer,p)).transpose().map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"analysis":analysis,"bindings":bindings,"partner":args.partner,"heterodimer":heterodimer}))?;
            output(json!({"document":args.document,"analysis":analysis,"heterodimer":heterodimer,"bindings":bindings.hits.iter().take(100).collect::<Vec<_>>(),"binding_count":bindings.hits.len(),"search_truncated":bindings.truncated,"display_truncated":bindings.hits.len()>100,"binding_options":bindings.options,"warnings":bindings.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Design bounded PCR primer pairs anchored to an explicit reference interval (circular end may be unwrapped). Uses recorded nearest-neighbor conditions, GC/length/homopolymer/complementarity filters and separated synthetic tails. Candidates require binding-alternative/PCR and laboratory review; no automatic verification status or document edits. Full JSON artifact retains source/options and ranked pairs."
    )]
    async fn design_primers(
        &self,
        Parameters(args): Parameters<PrimerDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.store.lock().unwrap().document(&args.document)?.clone();
            let report=primers::design_primers(&source,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"design":report}))?;
            output(json!({"document":args.document,"options":report.options,"pairs":report.pairs.iter().take(20).collect::<Vec<_>>(),"count":report.pairs.len(),"display_truncated":report.pairs.len()>20,"evaluated_pairs":report.evaluated_pairs,"rejected_pairs":report.rejected_pairs,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "List the authenticated Addgene catalog's sequence choices for a plasmid URL/ID: deposited versus Addgene-verified, full versus partial, sequence ID and declared length. Configure COLLAGENE_ADDGENE_API_TOKEN in the server environment. Does not import or silently select a record; pass sequence_id to import_url after review. No raw sequence, secret or temporary download URL is returned."
    )]
    async fn list_source_sequences(
        &self,
        Parameters(args): Parameters<SourceQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let _permit = self
            .compute_gate
            .clone()
            .acquire_owned()
            .await
            .map_err(bad)?;
        let client = self
            .source_client
            .get_or_try_init(|| async {
                collagene_sources::SourceClient::new(
                    collagene_sources::Config::from_env().map_err(bad)?,
                )
                .map_err(bad)
            })
            .await?;
        let catalog = client.addgene_catalog(&args.url).await.map_err(bad)?;
        let artifact = self.analysis_artifact(&serde_json::to_value(&catalog).map_err(bad)?)?;
        output(
            json!({"plasmid_id":catalog.plasmid_id(),"name":catalog.name(),"response_blake3":catalog.response_blake3(),"count":catalog.choices().len(),"display_truncated":catalog.choices().len()>100,"choices":catalog.choices().iter().take(100).map(|choice|json!({"sequence_id":choice.sequence_id,"description":choice.description.chars().take(512).collect::<String>(),"description_truncated":choice.description.chars().count()>512,"origin":choice.origin,"completeness":choice.completeness,"length":choice.length,"catalog_sequence_blake3":choice.catalog_sequence_blake3})).collect::<Vec<_>>(),"artifact_id":artifact}),
        )
    }
    #[tool(
        description = "Translate a revision-bound DNA/RNA sequence with an explicit NCBI genetic code, frame, orientation and initiation policy. Returns bounded protein preview and a full JSON evidence artifact; ambiguous codons and possible stops remain explicit."
    )]
    async fn translate_sequence(
        &self,
        Parameters(args): Parameters<TranslateQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let molecule=server.store.lock().unwrap().document(&args.document)?.clone();
            let report=translation::translate(&molecule.sequence,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"translation":report}))?;
            output(json!({"document":args.document,"options":report.options,"code_name":report.code_name,"amino_acids":report.protein.len(),"protein_preview":report.protein.chars().take(2000).collect::<String>(),"preview_truncated":report.protein.len()>2000,"unresolved_codons":report.unresolved_codons,"possible_stop_codons":report.possible_stop_codons,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Translate a specific CDS annotation using its compound location, codon_start and transl_table. Unsupported exceptions/fuzzy locations fail explicitly. Stored-translation agreement is evidence, not experimental confirmation. Full report available as JSON artifact."
    )]
    async fn translate_cds(
        &self,
        Parameters(args): Parameters<CdsQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let molecule=server.store.lock().unwrap().document(&args.document)?.clone();
            let feature=molecule.features.get(args.feature_index).ok_or_else(||bad("CDS annotation index does not exist"))?;
            let report=translation::translate_cds(&molecule,feature).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"cds":report}))?;
            output(json!({"document":args.document,"feature_index":args.feature_index,"location":report.location,"code":report.translation.options.genetic_code,"protein_preview":report.protein.chars().take(2000).collect::<String>(),"amino_acids":report.protein.len(),"preview_truncated":report.protein.len()>2000,"complete_under_model":report.complete_under_model,"stored_translation_matches":report.stored_translation_matches,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Find bounded six-frame ORF candidates with explicit start/partial/code policies. Circular searches cover at most one source turn. Predictions carry exact/partial/possible-stop evidence; they are not validated genes. Returns first100 summaries plus full JSON artifact."
    )]
    async fn find_orfs(
        &self,
        Parameters(args): Parameters<OrfQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let molecule=server.store.lock().unwrap().document(&args.document)?.clone();
            let report=translation::find_orfs(&molecule,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"orfs":report}))?;
            output(json!({"document":args.document,"options":report.options,"count":report.orfs.len(),"truncated_search":report.truncated,"display_truncated":report.orfs.len()>100,"orfs":report.orfs.iter().take(100).map(|o|json!({"location":o.location,"reverse":o.reverse,"frame":o.frame,"amino_acids":o.protein.len(),"protein_preview":o.protein.chars().take(80).collect::<String>(),"crosses_origin":o.crosses_origin,"partial_5prime":o.partial_5prime,"partial_3prime":o.partial_3prime,"termination":o.termination})).collect::<Vec<_>>(),"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Run embedded affine-gap pairwise alignment on two immutable nucleotide revisions. Explicit global/local/query-in-reference mode, scoring and matrix budget; no circular-origin/orientation search in this tool. Ambiguous-compatible columns are distinct from certain matches. Full CIGAR/differences in JSON artifact."
    )]
    async fn align_sequences(
        &self,
        Parameters(args): Parameters<AlignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (reference,query)={let store=server.store.lock().unwrap();(store.document(&args.reference)?.clone(),store.document(&args.query)?.clone())};
            let report=alignment::align(&reference.sequence,&query.sequence,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"reference":args.reference,"query":args.query,"engine_version":verification::ENGINE_VERSION,"alignment":report}))?;
            output(json!({"reference":args.reference,"query":args.query,"algorithm":report.algorithm,"options":report.options,"score":report.score,"reference_range":[report.reference_start,report.reference_end],"query_range":[report.query_start,report.query_end],"stats":report.stats,"cigar_preview":report.cigar.chars().take(2000).collect::<String>(),"cigar_truncated":report.cigar.len()>2000,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Compare immutable molecule sequences with bounded circular-origin and reverse-query candidate search. Returns tested counts, completeness and optimality flags separately; seed/budget-limited search is never labeled exhaustive. Exact identity does not establish biological function or uniqueness."
    )]
    async fn compare_molecules(
        &self,
        Parameters(args): Parameters<CompareQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (reference,query)={let store=server.store.lock().unwrap();(store.document(&args.reference)?.clone(),store.document(&args.query)?.clone())};
            let report=alignment::compare_molecules(&reference,&query,&args.options).map_err(bad)?;
            let artifact=server.analysis_artifact(&json!({"reference":args.reference,"query":args.query,"engine_version":verification::ENGINE_VERSION,"comparison":report}))?;
            output(json!({"reference":args.reference,"query":args.query,"reference_origin":report.reference_origin,"query_reversed":report.query_reversed,"score":report.alignment.score,"stats":report.alignment.stats,"evaluated_candidates":report.evaluated_candidates,"exhaustive_reference_origins":report.exhaustive_reference_origins,"exhaustive_query_orientations":report.exhaustive_query_orientations,"optimality_proven":report.optimality_proven,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Import one bounded public nucleotide record using embedded Rust HTTP: NCBI accession/GenBank URL, explicit ena:ACCESSION.VERSION or official ENA URL, explicit ddbj:ACCESSION.VERSION or official DDBJ URL, or direct HTTPS sequence file. Bare accessions remain NCBI and numeric IDs remain Addgene. ENA/DDBJ record identity, declared lengths, original source bytes and provenance are checked by the shared source adapter; lists/search/ranges are not whole-record imports. Addgene requires an approved server token and explicit catalog sequence_id, which is Addgene-only. Public HTTPS/redirect/DNS limits remain enforced; no filesystem paths, cookie forwarding or private targets. Existing 32 MiB download and ordinary molecule limits remain; this tool does not stream large remote genomes. NCBI/Addgene credentials stay provider-scoped and out of metadata."
    )]
    async fn import_url(
        &self,
        Parameters(args): Parameters<ImportUrl>,
    ) -> Result<CallToolResult, ErrorData> {
        let _permit = self
            .compute_gate
            .clone()
            .acquire_owned()
            .await
            .map_err(bad)?;
        let client = self
            .source_client
            .get_or_try_init(|| async {
                let config = collagene_sources::Config::from_env().map_err(bad)?;
                collagene_sources::SourceClient::new(config).map_err(bad)
            })
            .await?;
        let report = client
            .import_url_selected(&args.url, args.sequence_id)
            .await
            .map_err(bad)?;
        self.insert_documents(report.molecules, report.warnings)
    }
    #[tool(
        description = "Export one exact session document revision as a native .cg2 artifact, preserving annotations and source evidence. Use this after native import_url when browser CORS prevents direct retrieval. Read base64 chunks, concatenate and decode before saving; opening the file in the browser is an explicit document import, not a project mutation."
    )]
    async fn export_document(
        &self,
        Parameters(reference): Parameters<DocumentRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move ||{
            let molecule=server.store.lock().unwrap().document(&reference)?.clone();
            let bytes=collagene_core::protocol::save_document(&molecule).map_err(bad)?;
            let id=server.binary_artifact(bytes.clone(),"application/octet-stream","sequence.cg2")?;
            output(json!({"artifact_id":id,"filename":"sequence.cg2","encoding":"base64","mime_type":"application/octet-stream","bytes":bytes.len(),"document":reference}))
        }).await
    }
    #[tool(
        description = "Review/export one immutable session revision as annotated GenBank, EMBL, GenPept or UniProt text using the native writers. GenBank/EMBL require explicit DNA/RNA documents; GenPept/UniProt require linear Protein. Protein profiles accept bounded explicit record identifiers and GenPept calendar export date; source accessions, review status, molecular mass and biological history are never invented. Returns exact writer warnings and no artifact until accept_omissions=true when warnings exist. Unsupported chemistry/locations/options remain errors even with acknowledgement. Output is bounded32MiB and session-only; retain native files/workspaces for exact source evidence and history."
    )]
    async fn export_annotated_sequence(
        &self,
        Parameters(args): Parameters<AnnotatedSequenceExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            use collagene_core::{Alphabet,protein_formats};
            let molecule=server.store.lock().unwrap().document(&args.document)?.clone();
            let protein=molecule.sequence.alphabet()==Alphabet::Protein;
            let (report,profile,filename)=match args.format {
                AnnotatedSequenceFormat::GenBank=> {
                    if protein{return Err(bad("GenBank requires DNA/RNA; use the explicit GenPept protein profile"));}
                    (formats::export_with_report(&molecule,formats::TextFormat::GenBank).map_err(bad)?,"genbank","sequence.gb")
                },
                AnnotatedSequenceFormat::Embl=> {
                    if protein{return Err(bad("EMBL requires DNA/RNA; use the explicit UniProt protein profile"));}
                    (formats::export_with_report(&molecule,formats::TextFormat::Embl).map_err(bad)?,"embl","sequence.embl")
                },
                profile=> {
                    if !protein{return Err(bad("GenPept/UniProt require an explicitly Protein document; no alphabet is inferred"));}
                    let (profile,options,label,filename)=match profile {
                        AnnotatedSequenceFormat::GenPept{locus_date,entry_name}=>{
                            if locus_date.as_ref().is_some_and(|v|v.len()!=11)||entry_name.as_ref().is_some_and(|v|v.len()>20){return Err(bad("GenPept date/entry identifier exceeds the supported profile bounds"));}
                            let date=locus_date.as_deref().map(protein_formats::InsdcDate::new).transpose().map_err(bad)?;
                            (protein_formats::Profile::GenPept,protein_formats::ExportOptions{locus_date:date,entry_name,uniprot_accession:None},"genpept","sequence.gp")
                        },
                        AnnotatedSequenceFormat::UniProt{entry_name,accession}=>{
                            if entry_name.as_ref().is_some_and(|v|v.len()>20)||accession.as_ref().is_some_and(|v|v.len()>32){return Err(bad("UniProt entry/accession exceeds the supported profile bounds"));}
                            (protein_formats::Profile::UniProt,protein_formats::ExportOptions{locus_date:None,entry_name,uniprot_accession:accession},"uniprot","sequence.dat")
                        },_=>unreachable!(),
                    };
                    (protein_formats::export_with_options(std::slice::from_ref(&molecule),profile,&options).map_err(bad)?,label,filename)
                },
            };
            if report.text.len()>32*1024*1024{return Err(bad("Annotated sequence export exceeds32MiB; use native export"));}
            if !args.accept_omissions&&!report.warnings.is_empty(){return output(json!({"document":args.document,"state":"review_required","profile":profile,"warnings":report.warnings,"artifact_id":null}));}
            let bytes=report.text.len();
            let artifact=server.text_artifact(report.text,"text/plain",filename)?;
            output(json!({"document":args.document,"source_revision":molecule.revision(),"state":"prepared","profile":profile,"warnings":report.warnings,"artifact_id":artifact,"bytes":bytes,"filename":filename,"mime_type":"text/plain","encoding":"utf8","persistence":"session only","authority":"Annotated exchange text is not a native history archive or new scientific verification."}))
        }).await
    }
    #[tool(
        description = "Import explicit FASTA, GCG/Wisconsin single-sequence or original Collagene XML profiles, at most 4 MiB of text. FASTA alphabet=null retains per-record declarations (unmarked records are DNA); an explicit alphabet handles unmarked RNA/protein and must agree with every declaration. GCG requires a caller-selected alphabet and verifies the declared type, length and checksum. XML lower-only migration requires AllowLowerOnly, retains exact original rows/bytes and never reconstructs physical cloning chemistry. Imports the entire selected project atomically as immutable session documents; no paths or external programs."
    )]
    async fn import_compat_sequence(
        &self,
        Parameters(args): Parameters<sequence_exchange_tools::ImportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            if args.text.len() > 4 * 1024 * 1024 {
                return Err(bad("Compatibility text imports are limited to 4 MiB"));
            }
            let report = match args.profile {
                sequence_exchange_tools::ImportProfile::Gcg(options) => {
                    collagene_core::gcg_sequence::import_gcg(args.text.as_bytes(), &options)
                }
                sequence_exchange_tools::ImportProfile::Fasta { alphabet } => {
                    let text = args.text.trim_start_matches('\u{feff}');
                    match alphabet {
                        Some(alphabet) => formats::import_fasta_with_alphabet(text, alphabet),
                        None => formats::import_fasta(text),
                    }
                }
                sequence_exchange_tools::ImportProfile::CollageneXml { policy } => {
                    collagene_core::collagene_xml_compat::import_with_policy(
                        args.text.as_bytes(),
                        policy,
                    )
                }
            }
            .map_err(bad)?;
            server.insert_documents(report.molecules, report.warnings)
        })
        .await
    }
    #[tool(
        description = "Review/export an exact immutable document revision under the independently implemented SnapGene revision-15 reference profile (.dna or .prot). Returns structured omissions and no artifact until accept_omissions=true when loss exists. Refuses unsupported physical/reference-only material under both policies. Emits fresh supported packets, never stale original opaque chunks or cloned experimental history; artifact is session-only and source remains unchanged."
    )]
    async fn export_snapgene(
        &self,
        Parameters(args): Parameters<sequence_exchange_tools::SnapGeneExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            use collagene_core::snapgene_export::{self,Policy};
            let molecule=server.store.lock().unwrap().document(&args.document)?.clone();
            let review=snapgene_export::review(&molecule).map_err(bad)?;
            if !review.losses.is_empty()&&!args.accept_omissions{return output(json!({"document":args.document,"state":"review_required","review":review,"artifact_id":null}));}
            let report=snapgene_export::export_snapgene(&molecule,if args.accept_omissions{Policy::AllowOmissions}else{Policy::Strict}).map_err(bad)?;
            let filename=match report.extension.as_str(){".dna"=>"sequence.dna",".prot"=>"sequence.prot",_=>return Err(bad("Unsupported SnapGene artifact extension"))};let size=report.bytes.len();
            let artifact=server.binary_artifact(report.bytes,"application/octet-stream",filename)?;
            output(json!({"document":args.document,"state":"prepared","profile":report.profile,"source_revision":report.source_revision,"losses":report.losses,"warnings":report.warnings,"artifact_id":artifact,"bytes":size,"filename":filename,"mime_type":report.mime,"encoding":"base64","persistence":"session only"}))
        }).await
    }
    #[tool(
        description = "Review/export an exact immutable DNA, RNA or protein revision as a GCG/Wisconsin single-sequence artifact with checksum. Requires explicit record ID/date for synthetic records; retained IDs/dates may be reused as documentary source data. Returns omission warnings without an artifact until accept_omissions=true when loss exists. GCG cannot carry annotations, circular topology, source attachments or cloning history. Unsupported physical/reference-only material is refused; no source mutation or host path."
    )]
    async fn export_gcg(
        &self,
        Parameters(args): Parameters<sequence_exchange_tools::GcgExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let molecule=server.store.lock().unwrap().document(&args.document)?.clone();
            let report=collagene_core::gcg_sequence::export_gcg(&molecule,&args.options).map_err(bad)?;
            if !report.warnings.is_empty()&&!args.accept_omissions{return output(json!({"document":args.document,"state":"review_required","warnings":report.warnings,"artifact_id":null}));}
            let bytes=report.text.into_bytes();let size=bytes.len();
            let artifact=server.binary_artifact(bytes,"text/plain","sequence.gcg")?;
            output(json!({"document":args.document,"state":"prepared","warnings":report.warnings,"artifact_id":artifact,"bytes":size,"filename":"sequence.gcg","mime_type":"text/plain","encoding":"base64","persistence":"session only"}))
        }).await
    }
    #[tool(
        description = "Import GenBank/EMBL/FASTA/Collagene XML or plain sequence text as immutable session documents. Atomic and deduplicated; returns IDs and exact revisions, without bulk sequence output."
    )]
    fn import_document(
        &self,
        Parameters(args): Parameters<ImportText>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.text.len() > 4 * 1024 * 1024 {
            return Err(bad("MCP text imports are limited to 4 MiB"));
        }
        let molecules = formats::import_text(&args.text).map_err(bad)?;
        let warnings = molecules.warnings;
        self.insert_documents(molecules.molecules, warnings)
    }
    fn insert_documents(
        &self,
        molecules: Vec<Molecule>,
        source_warnings: Vec<String>,
    ) -> Result<CallToolResult, ErrorData> {
        let mut pending = BTreeMap::new();
        let warnings_truncated = source_warnings.len() > 100;
        let warnings: Vec<_> = source_warnings
            .into_iter()
            .take(100)
            .map(|w| w.chars().take(1024).collect::<String>())
            .collect();
        for m in molecules {
            m.validate().map_err(bad)?;
            pending.insert(m.revision(), m);
        }
        let mut store = self.store.lock().unwrap();
        let mut bytes = 0;
        let mut count = 0;
        for (id, m) in &pending {
            if !store.documents.contains_key(id) {
                bytes += serde_json::to_vec(m).map_err(bad)?.len();
                count += 1;
            }
        }
        if store.documents.len() + count > MAX_DOCUMENT_REVISIONS {
            return Err(bad("at most 10000 session document revisions"));
        }
        store.reserve(bytes)?;
        let result: Vec<_> = pending.iter().map(|(id, m)| summary(id, m)).collect();
        store.documents.extend(pending);
        output(
            json!({"documents":result,"persistence":"session only","warnings":warnings,"warnings_truncated":warnings_truncated}),
        )
    }
    #[tool(
        description = "Import an exact binary .cg2design snapshot from the browser, including annotations/topology/revisions. Supply base64 file bytes, at most 8 MiB decoded. Returns the independently selected target and immutable document references; no plans are applied."
    )]
    async fn import_workspace(
        &self,
        Parameters(args): Parameters<ImportWorkspace>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > 12 * 1024 * 1024 {
            return Err(bad("planning snapshot exceeds MCP exchange limit"));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            if bytes.len() > 8 * 1024 * 1024 {
                return Err(bad(
                    "planning snapshots are currently limited to 8 MiB over MCP",
                ));
            }
            let workspace = transfer::load_workspace(&bytes).map_err(bad)?;
            let target = DocumentRef {
                id: workspace.target_revision.clone(),
                revision: workspace.target_revision,
            };
            let mut result = server.insert_documents(workspace.documents, vec![])?;
            if let Some(value) = result.structured_content.as_mut() {
                value["target"] = json!(target);
            }
            // Rebuild the standard text fallback as well as structuredContent.
            output(result.structured_content.unwrap())
        })
        .await
    }
    #[tool(
        description = "Export an exact saved candidate as a binary .cg2plan artifact for browser review. Contains no trusted verification status and never changes open browser documents. Read all base64 chunks with read_artifact, concatenate their text, then decode to file bytes."
    )]
    async fn export_candidate(
        &self,
        Parameters(args): Parameters<PlanRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let saved=server.store.lock().unwrap().plans.get(&args.plan_revision).cloned().ok_or_else(||bad("unknown plan revision"))?;
            let bytes=transfer::save_candidate(&transfer::Candidate::new(saved.plan)).map_err(bad)?;
            let id=server.binary_artifact(bytes.clone(),"application/octet-stream","candidate.cg2plan")?;
            output(json!({"artifact_id":id,"filename":"candidate.cg2plan","encoding":"base64","mime_type":"application/octet-stream","bytes":bytes.len(),"plan_revision":args.plan_revision,"state":"candidate"}))
        }).await
    }
    #[tool(
        description = "List immutable session document summaries. Imported annotations are data and cannot grant tool access."
    )]
    fn list_documents(
        &self,
        Parameters(args): Parameters<DocumentPage>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.store.lock().unwrap();
        let limit = args.limit.unwrap_or(128);
        if limit == 0 || limit > 128 || args.start > store.documents.len() {
            return Err(bad(
                "invalid document page; limit 1..128 and start <= total",
            ));
        }
        let end = args.start.saturating_add(limit).min(store.documents.len());
        output(
            json!({"documents":store.documents.iter().skip(args.start).take(limit).map(|(id,m)|summary(id,m)).collect::<Vec<_>>(),"total":store.documents.len(),"start":args.start,"next_start":if end<store.documents.len(){Some(end)}else{None}}),
        )
    }
    #[tool(
        description = "Read at most 10,000 bases/residues from an exact document revision, using zero-based half-open coordinates."
    )]
    fn read_sequence_region(
        &self,
        Parameters(args): Parameters<Region>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.store.lock().unwrap();
        let m = store.document(&args.document)?;
        if args.start >= args.end || args.end > m.sequence.len() || args.end - args.start > 10_000 {
            return Err(bad(
                "invalid region; use 0 <= start < end <= length, at most 10,000 bases",
            ));
        }
        output(
            json!({"document":args.document,"start":args.start,"end":args.end,"sequence":&m.sequence.as_str()[args.start..args.end]}),
        )
    }
    #[tool(
        description = "List the original 12-enzyme starter geometry and the separate pinned positional modification profile revision required by CheckedDigest plan steps. Full catalog documentary evidence uses lookup_enzyme_catalog/resolve_enzyme_evidence; no supplier facts are inferred by this legacy list. Standalone methylation assessments use a distinct profile; permission under declared tracked-generation assumptions does not establish empirical activity."
    )]
    fn list_enzymes(&self) -> Result<CallToolResult, ErrorData> {
        output(
            json!({"enzymes":analysis::starter_enzymes(),"checked_digest_profile_revision":collagene_core::modified_duplex::profile_revision(),"checked_digest_profile_version":collagene_core::modified_duplex::PROFILE_VERSION,"scope":"Legacy starter geometry only; no catalog supplier/formulation facts are inferred.","catalog_discovery":catalog_tools::identity(&catalog_tools::builtin()?),"catalog_lookup_tool":"lookup_enzyme_catalog","catalog_evidence_tool":"resolve_enzyme_evidence"}),
        )
    }
    #[tool(
        description = "Scan an immutable DNA source for restriction recognition spans and reference top/bottom cut boundaries without constructing duplexes, fragments or products. Select ordered exact catalog IDs or uniquely resolved names/aliases, optionally pinning edition/hash/full geometry. Certain and possible recognition are separate; unavailable cuts retain raw boundary diagnostics. Returns bounded site pages, validated source/geometry/options evidence and full JSON artifact. Per-enzyme observed counts are lower bounds when scan_complete=false; exact counts and single-cutter classification are null on incomplete scans, and possible recognition withholds single-cutter certainty. Scan completion describes sequence search only, not enzyme activity, methylation, buffer compatibility or empirical cleavage. Use inspect_restriction separately only when assumed-complete physical digest summaries are wanted."
    )]
    async fn scan_restriction(
        &self,
        Parameters(args): Parameters<RestrictionScanQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.max_sites) {
            return Err(bad("restriction site page size must be 1..100"));
        }
        if args.enzymes.is_empty() || args.enzymes.len() > 64 {
            return Err(bad("select 1..64 embedded enzymes"));
        }
        let catalog = catalog_tools::builtin()?;
        let ids = args
            .enzymes
            .iter()
            .map(|name| match catalog.lookup_enzyme(name) {
                collagene_enzyme_catalog::Lookup::Unique(enzyme) => Ok(enzyme.id.clone()),
                collagene_enzyme_catalog::Lookup::Missing => {
                    Err(bad("unknown enzyme identity; call lookup_enzyme_catalog"))
                }
                collagene_enzyme_catalog::Lookup::Ambiguous(_) => {
                    Err(bad("ambiguous enzyme alias; select exact catalog IDs"))
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        let selection = collagene_core::catalog_digest::select(&catalog, &ids).map_err(bad)?;
        if let Some(pin) = &args.catalog_pin {
            pin.validate_for(&catalog).map_err(bad)?;
            if pin != &selection {
                return Err(bad("catalog pin and ordered enzyme selection differ"));
            }
        }
        let enzymes = selection.resolved_enzymes(&catalog).map_err(bad)?;
        let server = self.clone();
        self.compute(move|| {
            let molecule={
                let store=server.store.lock().unwrap();
                let molecule=store.document(&args.document)?;
                collagene_core::restriction_inspection::source_preflight(&args.source_role,molecule,&enzymes,&args.options).map_err(bad)?;
                molecule.clone()
            };
            let report=collagene_core::restriction_inspection::scan(&args.source_role,&molecule,&enzymes,&args.options).map_err(bad)?;
            report.validate_bindings_for(&args.source_role,&molecule,&enzymes,&args.options).map_err(bad)?;
            if args.site_start>report.sites.len() {return Err(bad("restriction page offset exceeds retained evidence count"));}
            let mut counts=vec![(0usize,0usize);enzymes.len()];
            for site in &report.sites {
                let entry=&mut counts[site.enzyme_index];
                if site.certain {entry.0+=1;}else{entry.1+=1;}
            }
            let counts:Vec<_>=counts.into_iter().enumerate().map(|(i,(certain,possible))|json!({
                "enzyme_index":i,"catalog_id":selection.enzymes[i].id,"enzyme":enzymes[i].name,
                "counts_complete":report.scan_complete,"certain_observed":certain,"possible_observed":possible,
                "observed_site_count":certain+possible,
                "exact_site_count":(report.scan_complete && possible==0).then_some(certain),
                "single_cutter":if report.scan_complete && possible==0 {Some(certain==1)}else{None},
            })).collect();
            let authority="Recognition/cut-coordinate search evidence only; no duplex or fragment construction, no experimental cleavage or reaction-condition authority. Observed counts are lower bounds when scan is incomplete; possible IUPAC matches do not establish exact biological site counts.";
            let artifact_id=server.text_artifact(bounded_json(&json!({"catalog_pin":selection,"report":report,"enzyme_counts":counts,"authority":authority}))?,"application/json","restriction-scan.json")?;
            let end=args.site_start.saturating_add(args.max_sites).min(report.sites.len());
            output(json!({
                "document":args.document,"source_role":report.source_role,"source_revision":report.source_revision,
                "source_length":report.source_length,"topology":report.topology,"engine_version":report.engine_version,
                "batch_revision":report.batch_revision,"geometry_revision":report.geometry_revision,
                "options_revision":report.options_revision,"options":report.options,"enzymes":report.enzymes,
                "catalog_pin":selection,"scan_complete":report.scan_complete,"products_determined":false,
                "fragment_count":0,"ref_top_lengths":null,"comparisons":report.comparisons,
                "limits":report.limits,"diagnostics":report.diagnostics,"site_count":report.sites.len(),
                "enzyme_counts":counts,"site_start":args.site_start,"next_site_start":end,
                "sites_done":end==report.sites.len(),"sites":&report.sites[args.site_start..end],
                "artifact_id":artifact_id,"authority":authority,
            }))
        }).await
    }
    #[tool(
        description = "Inspect catalog-pinned cleavage of the actual outer-PCR product from an explicit four-step overlap-extension preparation. Freshly reconstructs both original-source PCR branches, selected physical strand extension and outer PCR before scanning/cutting its actual product; no intended target or claimed product is an input. Retains each stage's geometry counts, terminal-pair uncertainty, source/origin pins and full immutable evidence. Optional saved origin requires its exact native workspace and fresh replay; imported physical summaries never grant substrate authority. Returns a bounded fragment page and explicit selectable/incomplete status; it creates no stock, candidate, target or chosen cloning fragment."
    )]
    async fn inspect_overlap_extension_digest(
        &self,
        Parameters(args): Parameters<OverlapExtensionDigestInspectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.max_fragments) {
            return Err(bad(
                "overlap-extension digest fragment page size must be1..100",
            ));
        }
        if args.origin.is_some() && args.workspace.is_none() {
            return Err(bad(
                "saved overlap-extension source requires its exact native workspace reference",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let workspace = args.workspace.as_ref().map(|reference|server.restore_workspace(reference)).transpose()?;
            let source = {
                let store = server.store.lock().unwrap();
                let source = store.document(&args.document)?;
                collagene_core::overlap_extension_digest_inspection::source_preflight(source, args.origin.as_ref(), &args.request, workspace.as_ref()).map_err(bad)?;
                source.clone()
            };
            let expected = collagene_core::overlap_extension_digest_inspection::inspection_revision(&source, args.origin.as_ref(), &args.request, workspace.as_ref()).map_err(bad)?;
            let preview = collagene_core::overlap_extension_digest_inspection::inspect(&source, args.origin.as_ref(), &args.request, workspace.as_ref()).map_err(bad)?;
            preview.validate_binding(&expected).map_err(bad)?;
            let count = preview.inspection.as_ref().map_or(0, |inspection|inspection.fragments.len());
            if args.fragment_start > count { return Err(bad("overlap-extension digest fragment offset exceeds retained evidence")); }
            let end = args.fragment_start.saturating_add(args.max_fragments).min(count);
            let artifact = server.text_artifact(bounded_json_limit(&preview, 32 * 1024 * 1024)?, "application/json", "overlap-extension-digest-inspection.json")?;
            output(json!({
                "document":args.document,"mode":"overlap-extension-pcr","engine_version":preview.engine_version,
                "context_revision":preview.context_revision,"source_revision":preview.source_revision,
                "origin_revision":preview.origin_revision,"origin_incomplete":preview.origin_incomplete,
                "selectable":preview.selectable(),"incomplete":preview.incomplete,
                "amplicon":preview.amplicon.as_ref().map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})),
                "pcr_stages":preview.pcr_stages.iter().take(3).map(|stage|json!({
                    "stage":stage.stage,"template_revision":stage.template_revision,"geometry_count":stage.geometries.len(),
                    "product_revision":stage.product_revision,"product_length":stage.product_length
                })).collect::<Vec<_>>(),
                "extension":preview.extension.as_ref().map(|extension|json!({
                    "homology_length":extension.homology.len(),"terminal_pairings":extension.terminal_pairings.overlaps.iter().take(100).collect::<Vec<_>>(),
                    "pairings_truncated":extension.terminal_pairings.overlaps.len()>100,
                    "complete":extension.terminal_pairings.complete,"unique":extension.terminal_pairings.unique,
                    "comparisons":extension.terminal_pairings.comparisons
                })),
                "products_determined":preview.inspection.as_ref().is_some_and(|i|i.products_determined),
                "scan_complete":preview.inspection.as_ref().map(|i|i.scan_complete),
                "fragment_count":count,"fragment_start":args.fragment_start,"next_fragment_start":end,"fragments_done":end==count,
                "fragments":preview.inspection.as_ref().map(|i|&i.fragments[args.fragment_start..end]),
                "findings":preview.findings.iter().take(100).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "findings_truncated":preview.findings.len()>100||preview.findings.iter().take(100).any(|s|s.chars().count()>1024),
                "artifact_id":artifact,
                "authority":"Fresh source-bound four-step preparation and declared complete-cleavage inspection only. No target agreement, fragment isolation, stock, candidate or experimental success is inferred. Keep original-source bindings and incomplete evidence when explicitly selecting a fragment for a separate cloning recipe."
            }))
        }).await
    }
    #[tool(
        description = "Inspect actual restriction fragments of a freshly reconstructed unique exact PCR product from an immutable original source, explicit ordered primers and a full catalog pin. Primer 5′ extensions participate in the actual product and catalog scan. Ambiguous/no PCR geometry withholds physical fragment inspection; uncertain or partial cleavage retains diagnostics. Returns a paged actual end-chemistry summary and immutable full JSON evidence. Optional saved-run origin requires an explicit native workspace proving exact original run/source/target events before fresh reconstruction. Optional pcr_options enables the explicit mismatch-PCR engine with its own complete typed artifact and visible incomplete/conditional-singleton evidence; omitting it retains exact-PCR behavior. It never chooses a fragment, creates a plan, copies a result into an independent target, publishes a product document or infers amplification/cleavage success."
    )]
    async fn inspect_pcr_digest(
        &self,
        Parameters(args): Parameters<PcrDigestInspectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.max_fragments) {
            return Err(bad("PCR digest fragment page size must be1..100"));
        }
        if args.origin.is_some() && args.workspace.is_none() {
            return Err(bad(
                "saved-run PCR digest inspection requires its exact native workspace reference",
            ));
        }
        if let Some(pcr_options) = args.pcr_options {
            let request = collagene_core::mismatch_pcr_digest_inspection::Request {
                base: args.request,
                pcr_options,
            };
            let server = self.clone();
            return self.compute(move||{
                let workspace=args.workspace.as_ref().map(|reference|server.restore_workspace(reference)).transpose()?;
                let source={let store=server.store.lock().unwrap();let source=store.document(&args.document)?;
                    collagene_core::mismatch_pcr_digest_inspection::source_preflight(source,args.origin.as_ref(),&request,workspace.as_ref()).map_err(bad)?;source.clone()};
                let expected=collagene_core::mismatch_pcr_digest_inspection::inspection_revision(&source,args.origin.as_ref(),&request,workspace.as_ref()).map_err(bad)?;
                let preview=collagene_core::mismatch_pcr_digest_inspection::inspect(&source,args.origin.as_ref(),&request,workspace.as_ref()).map_err(bad)?;
                preview.validate_binding(&expected).map_err(bad)?;
                let count=preview.inspection.as_ref().map_or(0,|r|r.fragments.len());
                if args.fragment_start>count{return Err(bad("mismatch-PCR digest fragment offset exceeds retained evidence"));}
                let end=args.fragment_start.saturating_add(args.max_fragments).min(count);
                let result=json!({"document":args.document,"mode":"mismatch-pcr","engine_version":preview.engine_version,"context_revision":preview.context_revision,"source_revision":preview.source_revision,
                    "origin_revision":preview.origin_revision,"origin_incomplete":preview.origin_incomplete,"origin_findings_total":preview.origin_findings_total,"origin_findings_truncated":preview.origin_findings_truncated,"pcr_options":preview.request.pcr_options,"incomplete":preview.incomplete,
                    "pcr_candidate_count":preview.pcr.report.products.len(),"pcr_search_complete":preview.pcr.report.search_complete,"pcr_enumeration_complete":preview.pcr.report.enumeration_complete,
                    "pcr_candidates":preview.pcr.report.products.iter().take(100).enumerate().map(|(index,p)|json!({"index":index,"forward_start":p.forward_start,"reverse_start":p.reverse_start,"length":p.sequence.len(),"certain":p.certain,"overlapping_domains":p.overlapping_domains,"substitution_count":p.substitutions.len()})).collect::<Vec<_>>(),"pcr_candidates_truncated":preview.pcr.report.products.len()>100,
                    "amplicon":preview.amplicon.as_ref().map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})),
                    "products_determined":preview.inspection.as_ref().is_some_and(|i|i.products_determined),"scan_complete":preview.inspection.as_ref().map(|i|i.scan_complete),
                    "fragment_count":count,"fragment_start":args.fragment_start,"next_fragment_start":end,"fragments_done":end==count,"fragments":preview.inspection.as_ref().map(|i|&i.fragments[args.fragment_start..end]),
                    "findings":preview.findings.iter().take(100).collect::<Vec<_>>(),"findings_truncated":preview.findings.len()>100,
                    "authority":"Explicit mismatch-PCR policy and fresh catalog cleavage inspection only. Conditional singleton and incomplete-search evidence remain explicit. The report creates no independent target, selected cloning product, physical stock or experimental success."});
                let artifact=server.text_artifact(bounded_json_limit(&preview,32*1024*1024)?,"application/json","mismatch-pcr-digest-inspection.json")?;
                let mut result=result;result["artifact_id"]=json!(artifact);output(result)
            }).await;
        }
        let server = self.clone();
        self.compute(move||{
            let frame=args.workspace.as_ref().map(|reference|server.restore_workspace(reference)).transpose()?;
            let source={let store=server.store.lock().unwrap();let source=store.document(&args.document)?;
                if let Some(origin)=&args.origin{collagene_core::pcr_digest_inspection::preflight_from_material(source,origin,&args.request,frame.as_ref().ok_or_else(||bad("Saved PCR digest source workspace is unavailable"))?).map_err(bad)?;}
                else{collagene_core::pcr_digest_inspection::source_preflight(source,&args.request).map_err(bad)?;}source.clone()};
            let (ordinary,from_material)=if let Some(origin)=&args.origin {
                let frame=frame.as_ref().ok_or_else(||bad("Saved PCR digest source workspace is unavailable"))?;
                let expected=collagene_core::pcr_digest_inspection::inspection_revision_from_material(&source,origin,&args.request,frame).map_err(bad)?;
                let preview=collagene_core::pcr_digest_inspection::inspect_from_material(&source,origin,&args.request,frame).map_err(bad)?;
                preview.validate_binding(&expected).map_err(bad)?;(None,Some(preview))
            }else{
                let expected=collagene_core::pcr_digest_inspection::inspection_revision(&source,&args.request).map_err(bad)?;
                let preview=collagene_core::pcr_digest_inspection::inspect(&source,&args.request).map_err(bad)?;
                preview.validate_binding(&expected).map_err(bad)?;(Some(preview),None)
            };
            let preview=from_material.as_ref().map(|r|&r.preview).or(ordinary.as_ref()).ok_or_else(||bad("PCR digest inspection unavailable"))?;
            let context=from_material.as_ref().map(|r|r.context_revision.as_str()).unwrap_or(&preview.context_revision);
            let origin_revision=from_material.as_ref().map(|r|&r.origin_revision);
            let count=preview.inspection.as_ref().map_or(0,|i|i.fragments.len());
            if args.fragment_start>count{return Err(bad("PCR digest fragment offset exceeds retained evidence"));}
            let end=args.fragment_start.saturating_add(args.max_fragments).min(count);
            let result=json!({"document":args.document,"engine_version":preview.engine_version,"context_revision":context,"origin_revision":origin_revision,"source_revision":preview.source_revision,
                "origin_incomplete":from_material.as_ref().map(|r|r.origin_incomplete),"origin_findings_total":from_material.as_ref().map(|r|r.origin_findings_total),"origin_findings_truncated":from_material.as_ref().map(|r|r.origin_findings_truncated),
                "pcr_geometry_count":preview.pcr_geometries.len(),"pcr_geometries":preview.pcr_geometries.iter().take(100).collect::<Vec<_>>(),"pcr_geometries_truncated":preview.pcr_geometries.len()>100,
                "amplicon_length":preview.amplicon.as_ref().map(|m|m.sequence.len()),"products_determined":preview.inspection.as_ref().is_some_and(|i|i.products_determined),
                "scan_complete":preview.inspection.as_ref().map(|i|i.scan_complete),"limits":preview.inspection.as_ref().map(|i|&i.limits),
                "fragment_count":count,"fragment_start":args.fragment_start,"next_fragment_start":end,"fragments_done":end==count,
                "fragments":preview.inspection.as_ref().map(|i|&i.fragments[args.fragment_start..end]),
                "findings":preview.findings.iter().take(100).collect::<Vec<_>>(),"findings_truncated":preview.findings.len()>100,
                "authority":"Fresh explicitly bound original or saved-material PCR and assumed-complete catalog digestion evidence only. Inspect actual outer PCR and cut-derived termini; no selected product, target acceptance, reaction yield or imported physical authority is created."});
            let (text,filename)=if let Some(material)=&from_material{(bounded_json_limit(material,32*1024*1024)?,"material-pcr-digest-inspection.json")}else{(bounded_json_limit(preview,32*1024*1024)?,"pcr-digest-inspection.json")};
            let artifact=server.text_artifact(text,"application/json",filename)?;
            let mut result=result;result["artifact_id"]=json!(artifact);output(result)
        }).await
    }
    #[tool(
        description = "Preview or install target-free native exact PCR, intentional-mismatch PCR, catalog digestion or saved-material end repair from an exact current source event in an immutable workspace. Preview tasks return bounded metadata plus complete reference/geometry artifacts; they create no material documents. Install freshly reconstructs actual source strands, validates the literal selected PCR geometry or actual digest fragment and atomically returns a NEW session workspace, receipt and native archive, with optional newly performed PCR primers. No independent target, isolation, yield, experimental success, database write or filesystem authority is fabricated. Pass saved-origin MaterialSelection handles and explicit PCR preview index/revision only; closed capsules and whole amplicon geometry are resolved natively, never copied into requests. Saved Run/Derived handles must match their exact retained source events; marked reference-only documents cannot bypass replay as plain originals."
    )]
    async fn material_derivation(
        &self,
        Parameters(args): Parameters<MaterialDerivationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            use collagene_core::material_derivation as native;
            #[derive(Serialize)]
            struct MaterialPreviewEvidence<'a,T:Serialize,P:Serialize>{workspace:&'a WorkspaceRef,source:&'a collagene_core::run_lineage::EventBinding,parameters:P,preview:&'a T}
            let mut workspace=server.restore_workspace(&args.workspace)?;
            match args.task {
                MaterialDerivationTask::PreviewPcr{source,origin,forward,reverse,template_strand}=>{
                    let origin=server.resolve_material_derivation_origin(&workspace,&source,&origin)?;
                    let molecule=server.material_derivation_source(&workspace,&source,&origin)?;
                    let context=native::pcr_preview_revision(&molecule,&origin,&forward,&reverse,template_strand).map_err(bad)?;
                    let documents=native::preview_pcr(&molecule,&origin,&forward,&reverse,template_strand).map_err(bad)?;
                    #[derive(Serialize)]
                    struct PcrEvidence<'a>{workspace:&'a WorkspaceRef,source:&'a collagene_core::run_lineage::EventBinding,context_revision:&'a str,forward:&'a Primer,reverse:&'a Primer,template_strand:collagene_core::chemistry::Strand,documents:&'a[verification::PcrDocument]}
                    let text=bounded_json_limit(&PcrEvidence{workspace:&args.workspace,source:&source,context_revision:&context,forward:&forward,reverse:&reverse,template_strand,documents:&documents},32*1024*1024)?;
                    let artifact=server.text_artifact(text,"application/json","material-pcr-preview.json")?;
                    let rows=documents.iter().take(100).enumerate().map(|(index,d)|json!({"index":index,"preview_revision":d.molecule.revision(),"length":d.molecule.sequence.len(),"reusable":documents.len()==1 && d.molecule.sequence.as_str().bytes().all(|b|b"ACGT".contains(&b)) && forward.ordered_sequence().bytes().chain(reverse.ordered_sequence().bytes()).all(|b|b"ACGT".contains(&b)),"template_strand":d.geometry.template_strand,"forward_start":d.geometry.forward_start,"reverse_start":d.geometry.reverse_start,"crosses_origin":d.geometry.crosses_origin,"forward_primer":d.geometry.forward_primer,"reverse_primer":d.geometry.reverse_primer})).collect::<Vec<_>>();
                    output(json!({"workspace":args.workspace,"source":source,"context_revision":context,"count":documents.len(),"previews":rows,"previews_truncated":documents.len()>100,"artifact_id":artifact,"authority":"Reference PCR previews only. Use the explicit preview index and revision for Install; the server recomputes full geometry internally. Full geometry remains optional artifact evidence. Multiple products do not establish isolation or reusable material."}))
                },
                MaterialDerivationTask::PreviewMismatchPcr{source,origin,forward,reverse,options}=>{
                    let origin=server.resolve_material_derivation_origin(&workspace,&source,&origin)?;
                    let molecule=server.material_derivation_source(&workspace,&source,&origin)?;
                    let context=native::mismatch_pcr_preview_revision(&molecule,&origin,&forward,&reverse,&options).map_err(bad)?;
                    let preview=native::preview_mismatch_pcr(&molecule,&origin,&forward,&reverse,&options).map_err(bad)?;
                    preview.validate_binding(&context).map_err(bad)?;
                    let artifact=server.text_artifact(bounded_json_limit(&MaterialPreviewEvidence{workspace:&args.workspace,source:&source,parameters:json!({"forward":forward,"reverse":reverse,"options":options}),preview:&preview},32*1024*1024)?,"application/json","material-mismatch-pcr-preview.json")?;
                    let report=&preview.preview.report;
                    let reusable=preview.reusable();
                    let rows=preview.preview.documents.iter().zip(&report.products).take(100).enumerate().map(|(index,(document,geometry))|json!({"index":index,"preview_revision":document.revision(),"length":document.sequence.len(),"certain":geometry.certain,"overlapping_domains":geometry.overlapping_domains,"forward_start":geometry.forward_start,"reverse_start":geometry.reverse_start,"crosses_origin":geometry.crosses_origin,"forward_primer":geometry.forward_primer_index,"reverse_primer":geometry.reverse_primer_index,"substitution_count":geometry.substitutions.len(),"reusable":reusable&&index==0})).collect::<Vec<_>>();
                    output(json!({"workspace":args.workspace,"source":source,"context_revision":context,"count":preview.preview.documents.len(),"reusable":reusable,"previews":rows,"previews_truncated":preview.preview.documents.len()>100,"search_complete":report.search_complete,"enumeration_complete":report.enumeration_complete,"incomplete":preview.incomplete,"findings":preview.findings,"findings_total":preview.findings_total,"findings_truncated":preview.findings_truncated,"artifact_id":artifact,"authority":"Source-bound mismatch-PCR reference preview with full literal ordered primers and explicit mismatch/3-prime/search policy. Install freshly recomputes the selected geometry and independently refuses unsupported or competing physical products; this preview grants no material or isolation."}))
                },
                MaterialDerivationTask::PreviewRepair{source,origin,method}=>{
                    let origin=server.resolve_material_derivation_origin(&workspace,&source,&origin)?;
                    let molecule=server.material_derivation_source(&workspace,&source,&origin)?;
                    let context=native::repair_preview_revision(&molecule,&origin,method).map_err(bad)?;
                    let preview=native::preview_repair(&molecule,&origin,method).map_err(bad)?;
                    preview.validate_binding(&context).map_err(bad)?;
                    let artifact=server.text_artifact(bounded_json_limit(&MaterialPreviewEvidence{workspace:&args.workspace,source:&source,parameters:json!({"method":method}),preview:&preview},32*1024*1024)?,"application/json","material-end-repair-preview.json")?;
                    output(json!({"workspace":args.workspace,"source":source,"context_revision":context,"method":method,"length":preview.molecule.sequence.len(),"preview_revision":preview.molecule.revision(),"source_duplex":preview.source_duplex,"source_modifications":preview.source_modifications,"repair":preview.repair,"incomplete":preview.incomplete,"findings":preview.findings,"findings_total":preview.findings_total,"findings_truncated":preview.findings_truncated,"artifact_id":artifact,"authority":"Fresh saved-material strand repair under the explicit method; actual synthesized/removed bases and inherited findings remain evidence. No independent target, isolated stock, experimental success or saved product is asserted by this preview."}))
                },
                MaterialDerivationTask::PreviewDigest{source,origin,selection,options}=>{
                    let origin=server.resolve_material_derivation_origin(&workspace,&source,&origin)?;
                    let molecule=server.material_derivation_source(&workspace,&source,&origin)?;
                    let role=collagene_core::restriction_history::source_role(&source).map_err(bad)?;
                    let context=native::digest_inspection_revision(&role,&molecule,&origin,&selection,&options).map_err(bad)?;
                    let preview=native::inspect_digest(&role,&molecule,&origin,&selection,&options).map_err(bad)?;
                    preview.validate_binding(&context).map_err(bad)?;
                    let artifact=server.text_artifact(bounded_json_limit(&preview,32*1024*1024)?,"application/json","material-digest-preview.json")?;
                    output(json!({"workspace":args.workspace,"source":source,"source_role":role,"context_revision":context,"products_determined":preview.inspection.products_determined,"scan_complete":preview.inspection.scan_complete,"fragment_count":preview.inspection.fragments.len(),"fragments":preview.inspection.fragments.iter().take(100).collect::<Vec<_>>(),"fragments_truncated":preview.inspection.fragments.len()>100,"incomplete":preview.incomplete,"findings":preview.findings,"findings_total":preview.findings_total,"findings_truncated":preview.findings_truncated,"artifact_id":artifact,"authority":"Fresh actual source-strand catalog inspection under declared complete cleavage. This creates no material document, selected physical stock or independent target assessment."}))
                },
                MaterialDerivationTask::Install{request}=>{
                    if request.name.trim().is_empty()||request.name.trim()!=request.name||request.name.len()>256||request.name.chars().any(char::is_control){return Err(bad("Declare a bounded nonempty material product name."));}
                    let origin=server.resolve_material_derivation_origin(&workspace,&request.source,&request.origin)?;
                    let molecule=server.material_derivation_source(&workspace,&request.source,&origin)?;
                    let derivation=match request.derivation {
                        MaterialDerivationOperation::ExactPcr{forward,reverse,selection,template_strand}=>{
                            if selection.selected_index>=10000||selection.expected_preview_revision.len()!=64||!selection.expected_preview_revision.bytes().all(|b|b.is_ascii_digit()||(b'a'..=b'f').contains(&b)){return Err(bad("Invalid explicit PCR preview index/revision."));}
                            let previews=native::preview_pcr(&molecule,&origin,&forward,&reverse,template_strand).map_err(bad)?;
                            let chosen=previews.into_iter().nth(selection.selected_index).ok_or_else(||bad("The explicitly selected PCR preview no longer exists."))?;
                            if chosen.molecule.revision()!=selection.expected_preview_revision {return Err(bad("The explicitly selected PCR preview changed; preview again."));}
                            native::Derivation::ExactPcr{forward,reverse,template_strand,selection:native::PcrSelection{selected_index:selection.selected_index,expected_preview_revision:selection.expected_preview_revision,geometry:chosen.geometry}}
                        },
                        MaterialDerivationOperation::MismatchPcr{forward,reverse,options,selection}=>{
                            if selection.selected_index>=10000||selection.expected_preview_revision.len()!=64||!selection.expected_preview_revision.bytes().all(|b|b.is_ascii_digit()||(b'a'..=b'f').contains(&b)){return Err(bad("Invalid explicit mismatch-PCR preview index/revision."));}
                            let preview=native::preview_mismatch_pcr(&molecule,&origin,&forward,&reverse,&options).map_err(bad)?;
                            let chosen=preview.preview.documents.get(selection.selected_index).ok_or_else(||bad("The explicitly selected mismatch-PCR preview no longer exists."))?;
                            if chosen.revision()!=selection.expected_preview_revision{return Err(bad("The explicitly selected mismatch-PCR preview changed; preview again."));}
                            let geometry=preview.preview.report.products.get(selection.selected_index).cloned().ok_or_else(||bad("The selected mismatch-PCR geometry is unavailable."))?;
                            native::Derivation::MismatchPcr{forward,reverse,options,selection:native::MismatchPcrSelection{selected_index:selection.selected_index,expected_preview_revision:selection.expected_preview_revision,geometry}}
                        },
                        MaterialDerivationOperation::EndRepair{method}=>native::Derivation::EndRepair{method},
                        MaterialDerivationOperation::CatalogDigest{selection,options,fragment_index}=>native::Derivation::CatalogDigest{selection,options,fragment_index},
                    };
                    let request=native::Request{source:request.source,origin,derivation,name:request.name,save_primers:request.save_primers};
                    // Selecting the requested retained source in this local copy
                    // changes presentation only; the stored input is immutable.
                    workspace.active_document=request.source.document_id.clone();
                    let expected=native::revision(&workspace,&request).map_err(bad)?;
                    let installed=native::prepare_install(&workspace,&request).map_err(bad)?;
                    installed.receipt.validate_binding(&expected).map_err(bad)?;
                    installed.receipt.validate_archive(&installed.workspace).map_err(bad)?;
                    let receipt=serde_json::to_value(&installed.receipt).map_err(bad)?;
                    let text=bounded_json_limit(&installed.receipt,32*1024*1024)?;
                    let archive=workspace::save_workspace(&installed.workspace).map_err(bad)?;
                    let result=server.commit_workspace_outputs_named(installed.workspace,Some((text,receipt)),Some(archive),"material-derivation.json","Target-free material derivation in a NEW immutable session workspace only; every reuse reconstructs native source history again.")?;
                    let mut value=result.structured_content.ok_or_else(||bad("missing installed workspace receipt"))?;
                    value["source_workspace"]=json!(args.workspace);
                    value["source_selection"]=json!("Only the derived working copy selects the requested source tab before installation; all source biological heads and the input session workspace remain unchanged.");
                    value["authority"]=json!("Target-free declared-model material derivation saved to a NEW immutable session workspace/archive only. No independent target was invented or assessed, and no experimental isolation/yield, local catalog write, filesystem write or external publication occurred. Reuse freshly replays the complete native record.");
                    output(value)
                },
            }
        }).await
    }
    #[tool(
        description = "Discover or explicitly replay a saved run product or generic target-free material derivation origin for an exact event in an immutable native workspace. No origin is automatically selected. Returns bounded origin choices, findings and complete immutable inspection JSON including any selected closed documentary origin. Every later material preparation must freshly replay again; imported summaries and equal DNA are not material authority. Creates no documents, candidates, targets or stock."
    )]
    async fn inspect_material_origin(
        &self,
        Parameters(args): Parameters<MaterialOriginQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            if args.selection.as_ref().is_some_and(|selection|selection.event()!=&args.event) {
                return Err(bad("selected origin belongs to a different exact event"));
            }
            let workspace=server.restore_workspace(&args.workspace)?;
            let expected=collagene_core::run_material::material_inspection_revision(&workspace,&args.event,&args.selection).map_err(bad)?;
            let inspection=collagene_core::run_material::inspect_material(&workspace,&args.event,&args.selection).map_err(bad)?;
            inspection.validate_binding(&expected).map_err(bad)?;
            if let Some(selected)=&inspection.selected {
                selected.validate_for_workspace(&workspace,&args.event).map_err(bad)?;
                if args.selection.as_ref()!=Some(&selected.selection()) { return Err(bad("replayed origin differs from the explicit selection")); }
            }
            let text=bounded_json_limit(&inspection,32*1024*1024)?;
            let artifact=server.text_artifact(text,"application/json","material-origin-inspection.json")?;
            output(json!({"workspace":args.workspace,"event":inspection.event,"context_revision":inspection.context_revision,"engine_version":inspection.engine_version,"origins":inspection.origins,"selection":args.selection,"origin_available":inspection.selected.is_some(),"findings":inspection.findings,"artifact_id":artifact,"authority":"A closed documentary replay request is returned only after fresh native reconstruction. Its serialized bytes confer no private material authority; every later preparation independently replays against exact retained sources."}))
        }).await
    }
    #[tool(
        description = "Discover documentary originating run/step choices for an exact saved DNA event in an explicitly supplied immutable native workspace. Returns every bounded origin selection and immutable full inspection JSON. No origin is automatically selected, no chemistry is restored and no document, candidate or target is created. Equal DNA content is not an event identity."
    )]
    async fn list_run_material_origins(
        &self,
        Parameters(args): Parameters<RunMaterialOriginsQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let workspace=server.restore_workspace(&args.workspace)?;
            let expected=collagene_core::run_material::inspection_revision(&workspace,&args.event,&None).map_err(bad)?;
            let inspection=collagene_core::run_material::inspect(&workspace,&args.event,&None).map_err(bad)?;
            inspection.validate_binding(&expected).map_err(bad)?;
            let text=bounded_json_limit(&inspection,32*1024*1024)?;
            let artifact=server.text_artifact(text,"application/json","run-material-origins.json")?;
            output(json!({"workspace":args.workspace,"event":inspection.event,"context_revision":inspection.context_revision,"engine_version":inspection.engine_version,"origins":inspection.origins,"selected":false,"findings":inspection.findings,"artifact_id":artifact,"authority":"Documentary origin discovery only. No private material, successful replay or independent target acceptance is imported."}))
        }).await
    }
    #[tool(
        description = "Explicitly select an exact originating run/step/event in an immutable native workspace and freshly build its bounded durable reconstruction capsule. Rechecks original source/target/product event ownership and current deterministic replay. Returns documentary capsule JSON, not transferable private tracked material or scientific authority. Creates no document/candidate/target and performs no automatic origin choice. Existing source-bound recipes can consume the capsule through explicit FromRun preparation and freshly verify again."
    )]
    async fn prepare_run_material(
        &self,
        Parameters(args): Parameters<RunMaterialPreparationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let workspace=server.restore_workspace(&args.workspace)?;
            let selected=Some(args.selection.clone());
            let expected=collagene_core::run_material::inspection_revision(&workspace,&args.selection.event,&selected).map_err(bad)?;
            let inspection=collagene_core::run_material::inspect(&workspace,&args.selection.event,&selected).map_err(bad)?;
            inspection.validate_binding(&expected).map_err(bad)?;
            if inspection.selected.as_ref().is_none_or(|capsule|capsule.selection!=args.selection){return Err(bad("explicit selected origin produced no matching reconstruction capsule"));}
            let text=bounded_json_limit(&inspection,32*1024*1024)?;
            let artifact=server.text_artifact(text,"application/json","run-material-preparation.json")?;
            output(json!({"workspace":args.workspace,"event":inspection.event,"context_revision":inspection.context_revision,"engine_version":inspection.engine_version,"selection":args.selection,"capsule_available":true,"findings":inspection.findings,"artifact_id":artifact,"authority":"Freshly built documentary reconstruction request only; the JSON capsule does not carry private tracked material. Every later recipe independently reconstructs from the exact retained original run inputs and compares its separate target."}))
        }).await
    }
    #[tool(
        description = "Inspect an immutable DNA source with one or more exact catalog enzyme identities or uniquely resolved names. One artifact_id contains complete tagged Original or SavedRun evidence with full edition/hash/ordered-identity/geometry binding; inline pages remain bounded summaries. Preserves certain/possible recognition sites, exact origin-split spans, cuts outside substrate, bounded-scan limits and separate scan_complete/products_determined flags. Returns paged site and duplex-fragment summaries plus a complete immutable JSON artifact, bound to source role/revision, exact enzyme geometry and options. Products assume complete cleavage without blocking methylation or conditions; reports are evidence, not experimental activity or permission to release imported products."
    )]
    async fn inspect_restriction(
        &self,
        Parameters(args): Parameters<RestrictionInspectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.max_sites) || !(1..=100).contains(&args.max_fragments) {
            return Err(bad("restriction site/fragment page sizes must be 1..100"));
        }
        if args.enzymes.is_empty() || args.enzymes.len() > 64 {
            return Err(bad("select 1..64 embedded enzymes"));
        }
        if args.origin.is_some() && args.workspace.is_none() {
            return Err(bad(
                "saved-run catalog inspection requires its exact native workspace reference",
            ));
        }
        let catalog = catalog_tools::builtin()?;
        let ids = args
            .enzymes
            .iter()
            .map(|name| match catalog.lookup_enzyme(name) {
                collagene_enzyme_catalog::Lookup::Unique(enzyme) => Ok(enzyme.id.clone()),
                collagene_enzyme_catalog::Lookup::Missing => {
                    Err(bad("unknown enzyme identity; call lookup_enzyme_catalog"))
                }
                collagene_enzyme_catalog::Lookup::Ambiguous(_) => {
                    Err(bad("ambiguous enzyme alias; select exact catalog IDs"))
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        let selection = collagene_core::catalog_digest::select(&catalog, &ids).map_err(bad)?;
        if let Some(pin) = &args.catalog_pin {
            pin.validate_for(&catalog).map_err(bad)?;
            if pin != &selection {
                return Err(bad("catalog pin and ordered enzyme selection differ"));
            }
        }
        let server = self.clone();
        self.compute(move || {
            let molecule = {
                let store = server.store.lock().unwrap();
                let molecule = store.document(&args.document)?;
                if let Some(origin)=&args.origin {collagene_core::run_catalog_digest::source_preflight(&catalog,&selection,&args.source_role,molecule,origin,&args.options).map_err(bad)?;}
                else{collagene_core::catalog_digest::source_preflight(&catalog, &selection, &args.source_role, molecule, &args.options).map_err(bad)?;}
                molecule.clone()
            };
            let(ordinary,from_run)=if let Some(origin)=&args.origin {
                let workspace=args.workspace.as_ref().ok_or_else(||bad("saved-run catalog source workspace unavailable"))?;
                server.attest_run_source(workspace,origin,&molecule)?;
                let report=collagene_core::run_catalog_digest::inspect(&catalog,&selection,&args.source_role,&molecule,origin,&args.options).map_err(bad)?;
                report.validate_bindings_for(&catalog,&selection,&args.source_role,&molecule,origin,&args.options).map_err(bad)?;
                (None,Some(report))
            }else{
                let pinned=collagene_core::catalog_digest::inspect(&catalog,&selection,&args.source_role,&molecule,&args.options).map_err(bad)?;
                pinned.validate_bindings_for(&catalog,&selection,&args.source_role,&molecule,&args.options).map_err(bad)?;(Some(pinned),None)
            };
            let pinned=from_run.as_ref().map(|r|&r.inspection).or(ordinary.as_ref()).ok_or_else(||bad("catalog inspection unavailable"))?;
            let report=&pinned.inspection;
            if args.site_start > report.sites.len() || args.fragment_start > report.fragments.len() {
                return Err(bad("restriction page offset exceeds retained evidence count"));
            }
            let site_end = args.site_start.saturating_add(args.max_sites).min(report.sites.len());
            let fragment_end = args.fragment_start.saturating_add(args.max_fragments).min(report.fragments.len());
            let mut response=json!({
                "document":args.document,"engine_version":report.engine_version,"source_role":report.source_role,
                "source_revision":report.source_revision,"source_length":report.source_length,"topology":report.topology,
                "batch_revision":report.batch_revision,"geometry_revision":report.geometry_revision,
                "options_revision":report.options_revision,"enzymes":report.enzymes,"options":report.options,
                "scan_complete":report.scan_complete,"products_determined":report.products_determined,
                "comparisons":report.comparisons,"limits":report.limits,"diagnostics":report.diagnostics,
                "ref_top_lengths":report.ref_top_lengths,"site_count":report.sites.len(),
                "site_start":args.site_start,"next_site_start":site_end,"sites_done":site_end==report.sites.len(),
                "sites":&report.sites[args.site_start..site_end],"fragment_count":report.fragments.len(),
                "fragment_start":args.fragment_start,"next_fragment_start":fragment_end,
                "fragments_done":fragment_end==report.fragments.len(),"fragments":&report.fragments[args.fragment_start..fragment_end],
                "catalog_pin":pinned.selection,"catalog_batch_revision":pinned.batch_revision,
                "catalog_engine_version":pinned.engine_version,
                "authority":"Recognition evidence and assumed-complete cleavage summaries only; possible/partial/edge evidence cannot authorize a digest product. Reconstruct selected products from current source inputs."
            });
            if let Some(run)=&from_run {
                response["context_revision"]=json!(run.context_revision);
                response["origin_revision"]=json!(run.origin_revision);
                response["origin_status"]=json!(run.origin_status);
                response["origin_findings"]=json!(run.findings);
                response["authority"]=json!("Fresh saved-run tracked material and catalog cleavage inspection only. Serialized evidence never grants material, isolation, target acceptance or permission to install a reference-only digest product.");
            }
            let complete_evidence = match (ordinary, from_run) {
                (Some(inspection), None) => RestrictionInspectionEvidence::Original {
                    inspection,
                },
                (None, Some(inspection)) => RestrictionInspectionEvidence::SavedRun {
                    inspection,
                },
                _ => return Err(bad("catalog inspection source evidence is inconsistent")),
            };
            let artifact_id = server.text_artifact(
                bounded_json_limit(&complete_evidence, 32 * 1024 * 1024)?,
                "application/json",
                "catalog-restriction-inspection.json",
            )?;
            response["artifact_id"] = json!(artifact_id);
            output(response)
        }).await
    }
    #[tool(
        description = "Compute restriction sites and reference-strand fragment lengths for one embedded enzyme and exact document revision."
    )]
    async fn analyze_restriction(
        &self,
        Parameters(args): Parameters<Restriction>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.analyze_restriction_inner(args))
            .await
    }
    fn analyze_restriction_inner(&self, args: Restriction) -> Result<CallToolResult, ErrorData> {
        let m = self.store.lock().unwrap().document(&args.document)?.clone();
        let enzyme = analysis::starter_enzymes()
            .into_iter()
            .find(|e| e.name.eq_ignore_ascii_case(&args.enzyme))
            .ok_or_else(|| bad("unknown enzyme; call list_enzymes"))?;
        let sites = analysis::restriction_sites(&m, &enzyme).map_err(bad)?;
        let lengths = analysis::digest_lengths(&m, &sites).map_err(bad)?;
        // Sites are bounded separately to avoid flooding a model context.
        output(
            json!({"document":args.document,"enzyme":enzyme,"site_count":sites.len(),"sites":sites.iter().take(100).collect::<Vec<_>>(),"sites_truncated":sites.len()>100,"fragment_count":lengths.len(),"fragment_lengths":lengths.iter().take(100).collect::<Vec<_>>(),"fragments_truncated":lengths.len()>100,"model":"reference-strand lengths; full duplex end chemistry and methylation are not modeled"}),
        )
    }
    #[tool(
        description = "Compute complete duplex digestion using embedded enzyme names. Reports separate top/bottom lengths, 5′/3′ overhangs and terminal chemistry in ascending top-cut order. Imported unresolved chemical evidence is not trusted; methylation/conditions are not modeled. Returns at most 100 summaries with an explicit truncation flag."
    )]
    async fn digest_fragments(
        &self,
        Parameters(args): Parameters<DigestQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let molecule=server.store.lock().unwrap().document(&args.document)?.clone();
            let catalog=analysis::starter_enzymes();
            let enzymes=args.enzymes.iter().map(|name|catalog.iter().find(|e|e.name.eq_ignore_ascii_case(name)).cloned().ok_or_else(||bad("unknown enzyme; call list_enzymes"))).collect::<Result<Vec<_>,_>>()?;
            let fragments=chemistry::digest(&chemistry::Duplex::from_molecule(&molecule).map_err(bad)?,&enzymes).map_err(bad)?;
            let summaries=fragments.iter().take(100).enumerate().map(|(index,f)|f.duplex.summary().map(|summary|json!({"fragment_index":index,"top_start":f.top_start,"bottom_start":f.bottom_start,"duplex":summary}))).collect::<Result<Vec<_>,_>>().map_err(bad)?;
            output(json!({"document":args.document,"fragment_count":fragments.len(),"fragments":summaries,"truncated":fragments.len()>100,"assumption":"complete cleavage, no blocking methylation; unknown free-end chemistry remains unknown"}))
        }).await
    }
    #[tool(
        description = "Compute exact PCR products from an immutable template. Returns bounded summaries, considers competing oligo orientations and single-oligo products; mismatches/secondary structure are outside this model."
    )]
    async fn simulate_pcr(
        &self,
        Parameters(args): Parameters<Pcr>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.simulate_pcr_inner(args)).await
    }
    fn simulate_pcr_inner(&self, args: Pcr) -> Result<CallToolResult, ErrorData> {
        let m = self.store.lock().unwrap().document(&args.document)?.clone();
        let products = analysis::simulate_pcr(&m, &args.forward, &args.reverse).map_err(bad)?;
        let summaries: Vec<_> = products.iter().take(100).map(|p|json!({"length":p.sequence.len(),"sequence_hash":blake3::hash(p.sequence.as_str().as_bytes()).to_hex().to_string()})).collect();
        output(
            json!({"document":args.document,"product_count":products.len(),"products":summaries,"truncated":products.len()>100,"model":"exact non-overlapping binding sites within one template turn; no experimental guarantee"}),
        )
    }
    #[tool(
        description = "Build an intended reference from explicit regions of actual immutable source documents, their order, orientation and desired topology. No pre-existing target is required and no PCR, digest, ligation or assembly simulator supplies expected DNA. Atomically retains an unverified intended-reference target document and full source-placement JSON. Selection is design intent, not evidence a preparation can produce it or chemical/experimental material. Pass returned target DocumentRef with the same inputs to design_overlap_primers and propose_overlap_assembly for separate fresh preparation and verification."
    )]
    async fn build_assembly_intent(
        &self,
        Parameters(args): Parameters<AssemblyIntentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad("assembly intent requires 1..128 exact source bindings"));
        }
        let server = self.clone();
        self.compute(move || {
            let origin_workspace=args.origin_workspace.as_ref().map(|r|server.restore_workspace(r)).transpose()?;
            let native=args.request.native(origin_workspace.as_ref())?;

            let inputs={let store=server.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty()||bound.id.len()>128||!seen.insert(bound.id.as_str()) {return Err(bad("empty, oversized or duplicate intent source ID"));}
                    sources.push((bound.id.as_str(),store.document(&bound.document)?));
                }
                collagene_core::assembly_intent::preflight(&sources,&native,origin_workspace.as_ref()).map_err(bad)?;
                sources.into_iter().map(|(id,molecule)|PlanInput{id:id.into(),revision:molecule.revision(),molecule:molecule.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::assembly_intent::build(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            let expected=collagene_core::assembly_intent::intent_revision(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a>{origin_workspace:&'a Option<WorkspaceRef>,inputs:&'a [BoundInput],report:&'a collagene_core::assembly_intent::Report,interpretation:&'static str}
            let evidence=Evidence{origin_workspace:&args.origin_workspace,inputs:&args.inputs,report:&report,
                interpretation:"Declared source-region/order/orientation intent. Saved source strands are freshly reconstructed from exact native origins; the requested new preparation is never simulated to choose expected DNA. This target is not prepared material, an observed sequence or an independently reconstructed cloning product. Verify actual molecular steps separately against this frozen intent."};
            let (artifact,target_document)=server.candidate_document_artifact_bounded(&evidence,Some(report.target.clone()),"assembly-intent.json",32*1024*1024)?;
            let target_document=target_document.ok_or_else(||bad("intended-reference document was not retained"))?;
            output(json!({"state":"intended-reference","verified":false,"preparation_simulated":false,
                "inputs":args.inputs,"context_revision":report.context_revision,"request":args.request,
                "target":{"id":target_document["id"],"revision":target_document["revision"]},
                "target_summary":target_document,"intended_length":report.target.sequence.len(),"topology":report.target.topology,
                "fragment_count":report.fragments.len(),"fragments":report.fragments,
                "warnings":report.warnings.iter().take(100).map(|w|w.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "warnings_truncated":report.warnings.len()>100||report.warnings.iter().take(100).any(|w|w.chars().count()>1024),
                "artifact_id":artifact,
                "authority":"Unverified intended-reference document only. Original source snapshots are unchanged. Region selection does not imply PCR amplification, restriction cuts, compatible ends, sealing, chemistry or experimental success; actual preparation must be reconstructed separately."}))
        }).await
    }
    #[tool(
        description = "Build an independently declared assembly reference from exact immutable source snapshots, explicit region/order/orientation, junction trims/inserts and expected-base role-relative substitutions. Mutation intent comes from the declaration, never primers or simulated products. Checks original expected bases before copying/building and atomically retains an unverified intended-reference session document plus complete mapping JSON. Repeated source contributions retain distinct roles. Original/historical source snapshots remain unchanged. Material selectors are resolved from explicit native event/state bindings in origin_workspace; reference-only rows remain independent source-document declarations. No new preparation, material or successful verification is created."
    )]
    async fn build_substituted_assembly_intent(
        &self,
        Parameters(args): Parameters<SubstitutedAssemblyIntentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad(
                "substituted intent requires1..32 exact original source bindings",
            ));
        }
        let server = self.clone();
        self.compute(move||{
            let origin_workspace=args.origin_workspace.as_ref().map(|r|server.restore_workspace(r)).transpose()?;
            let native=args.request.native(origin_workspace.as_ref())?;

            let inputs={let store=server.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs{if bound.id.trim().is_empty()||bound.id.len()>128||!seen.insert(bound.id.as_str()){return Err(bad("invalid or duplicate original substitution source ID"));}sources.push((bound.id.as_str(),store.document(&bound.document)?));}
                collagene_core::assembly_substitution_intent::source_preflight(&sources,&native,origin_workspace.as_ref()).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::assembly_substitution_intent::build(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            let expected=collagene_core::assembly_substitution_intent::intent_revision(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            #[derive(Serialize)]struct Evidence<'a>{origin_workspace:&'a Option<WorkspaceRef>,inputs:&'a[BoundInput],report:&'a collagene_core::assembly_substitution_intent::Report,binding_scope:&'static str}
            let scope="Exact immutable source roles/revisions; material rows additionally attest explicit native events/state through the pinned origin_workspace. Reference-only rows declare correspondence, never material ancestry.";
            let evidence=Evidence{origin_workspace:&args.origin_workspace,inputs:&args.inputs,report:&report,binding_scope:scope};
            let(artifact,target)=server.candidate_document_artifact_bounded(&evidence,Some(report.target.clone()),"substituted-assembly-intent.json",32*1024*1024)?;
            let target=target.ok_or_else(||bad("substituted intended reference was not retained"))?;
            output(json!({"state":"intended-reference","verified":false,"preparation_simulated":false,"inputs":args.inputs,"context_revision":report.context_revision,"request":args.request,
                "target":{"id":target["id"],"revision":target["revision"]},"target_summary":target,"intended_length":report.target.sequence.len(),"topology":report.target.topology,
                "fragments":report.fragments,"applied":report.applied,"warnings":report.warnings,"artifact_id":artifact,"binding_scope":scope,
                "authority":"Independent source-bound design requirement only; original sources are unchanged. No primer-derived target, physical preparation, experimental success or native event ancestry is inferred."}))
        }).await
    }
    #[tool(
        description = "Build an independently declared assembly reference from exact immutable source snapshots, explicit region/order/orientation, junction trims/inserts, expected-base substitutions and literal contribution insertions/deletions/replacements. Mutation intent comes from the declaration, never primers or simulated products. Checks original expected literal bases before copying/building and atomically retains an unverified intended-reference session document plus complete mapping JSON. Repeated source contributions retain distinct roles. Original/historical source snapshots remain unchanged. Material selectors are resolved from explicit native event/state bindings in origin_workspace; reference-only rows remain independent source-document declarations. No new preparation, material or successful verification is created."
    )]
    async fn build_contribution_edit_intent(
        &self,
        Parameters(args): Parameters<ContributionEditIntentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad(
                "contribution edit intent requires1..32 exact original source bindings",
            ));
        }
        let server = self.clone();
        self.compute(move||{
            let origin_workspace=args.origin_workspace.as_ref().map(|r|server.restore_workspace(r)).transpose()?;
            let native=args.request.native(origin_workspace.as_ref())?;

            let inputs={let store=server.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs{if bound.id.trim().is_empty()||bound.id.len()>128||!seen.insert(bound.id.as_str()){return Err(bad("invalid or duplicate original contribution edit source ID"));}sources.push((bound.id.as_str(),store.document(&bound.document)?));}
                collagene_core::assembly_edit_intent::source_preflight(&sources,&native,origin_workspace.as_ref()).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::assembly_edit_intent::build(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            let expected=collagene_core::assembly_edit_intent::intent_revision(&inputs,&native,origin_workspace.as_ref()).map_err(bad)?;
            report.validate_sources(&inputs,&expected,origin_workspace.as_ref()).map_err(bad)?;
            #[derive(Serialize)]struct Evidence<'a>{origin_workspace:&'a Option<WorkspaceRef>,inputs:&'a[BoundInput],report:&'a collagene_core::assembly_edit_intent::Report,binding_scope:&'static str}
            let scope="Exact immutable source roles/revisions; material rows additionally attest explicit native events/state through the pinned origin_workspace. Reference-only rows declare correspondence, never material ancestry.";
            let evidence=Evidence{origin_workspace:&args.origin_workspace,inputs:&args.inputs,report:&report,binding_scope:scope};
            let(artifact,target)=server.candidate_document_artifact_bounded(&evidence,Some(report.target.clone()),"contribution-edit-intent.json",32*1024*1024)?;
            let target=target.ok_or_else(||bad("contribution-edited intended reference was not retained"))?;
            output(json!({"state":"intended-reference","verified":false,"preparation_simulated":false,"inputs":args.inputs,"context_revision":report.context_revision,"request":args.request,
                "target":{"id":target["id"],"revision":target["revision"]},"target_summary":target,"intended_length":report.target.sequence.len(),"topology":report.target.topology,
                "fragments":report.fragments,"applied":report.applied,"warnings":report.warnings,"artifact_id":artifact,"binding_scope":scope,
                "authority":"Independent source-bound design requirement only; original sources are unchanged. No primer-derived target, physical preparation, experimental success or native event ancestry is inferred."}))
        }).await
    }
    #[tool(
        description = "Build a declared edited assembly reference from exact original source regions and explicit junction deletions/insertions, independently of PCR or assembly simulation. Each edit trims bases from the upstream oriented core and downstream oriented core, then inserts literal DNA at that junction. Retains original order/orientation and an explicit circular origin. Atomically returns an unverified intended target document and complete JSON. Source documents stay unchanged; physical preparation must be reconstructed separately."
    )]
    async fn build_edited_assembly_intent(
        &self,
        Parameters(args): Parameters<EditedAssemblyIntentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 actual original source bindings"));
        }
        let server = self.clone();
        self.compute(move||{
            let origin_workspace=args.origin_workspace.as_ref().map(|r|server.restore_workspace(r)).transpose()?;
            let native=args.base.native(origin_workspace.as_ref())?;

            let inputs={let store=server.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs{if input.id.trim().is_empty()||input.id.len()>128||!seen.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::assembly_junction_edit::intent_preflight(&sources,&native,&args.edits,origin_workspace.as_ref()).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let intent=collagene_core::assembly_junction_edit::build_intent(&inputs,&native,&args.edits,origin_workspace.as_ref()).map_err(bad)?;
            let expected=collagene_core::assembly_junction_edit::intent_revision(&inputs,&native,&args.edits,origin_workspace.as_ref()).map_err(bad)?;intent.validate_binding(&expected).map_err(bad)?;
            let evidence=json!({"inputs":args.inputs,"origin_workspace":args.origin_workspace,"intent":intent});
            let (artifact,target)=server.candidate_document_artifact_bounded(&evidence,Some(intent.target.clone()),"edited-assembly-intent.json",32*1024*1024)?;
            let target=target.ok_or_else(||bad("edited intended reference was not retained"))?;
            output(json!({"state":"intended-reference","verified":false,"preparation_simulated":false,"inputs":args.inputs,
                "context_revision":intent.context_revision,"base":args.base,"edits":intent.edits,
                "target":{"id":target["id"],"revision":target["revision"]},"target_summary":target,"fragments":intent.fragments,
                "warnings":intent.warnings,"artifact_id":artifact,"authority":"Declared source-region edits only. Original sources are unchanged; primer production, assembly outcome and experimental success remain separate evidence."}))
        }).await
    }
    #[tool(
        description = "Propose original-source PCR primers for an explicitly junction-edited overlap assembly. Uses the declared trimmed source cores and literal inserted DNA to form the independent intended reference before designing primers. Added primer tails include the requested junction DNA plus next-core homology; full oligo metrics and existing filters are recomputed. Returns an advisory explicit-PCR recipe, no automatic application or verified cloning result. Reconstruct that recipe against the independently built edited target using propose_overlap_assembly."
    )]
    async fn design_edited_assembly_primers(
        &self,
        Parameters(args): Parameters<EditedAssemblyPrimersQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 actual original source bindings"));
        }
        let server = self.clone();
        self.compute(move||{
            let inputs={let store=server.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs{if input.id.trim().is_empty()||input.id.len()>128||!seen.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::assembly_junction_edit::preflight(&sources,&args.request).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::assembly_junction_edit::design(&inputs,&args.request).map_err(bad)?;
            let expected=collagene_core::assembly_junction_edit::design_revision(&inputs,&args.request).map_err(bad)?;report.validate_binding(&expected).map_err(bad)?;
            let evidence=json!({"inputs":args.inputs,"report":report});
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","edited-assembly-primers.json")?;
            let ordered=report.recipe.fragments.iter().map(|f|match &f.preparation{collagene_core::assembly_recipe::Preparation::ExactPcr{forward,reverse, .. }=>Ok(json!({"role":f.role,"input_id":f.input_id,"reverse_fragment":f.reverse,"forward":forward,"reverse":reverse,"forward_ordered":forward.ordered_sequence(),"reverse_ordered":reverse.ordered_sequence()})),_=>Err(bad("edited primer proposal lost its exact original-source PCR"))}).collect::<Result<Vec<_>,_>>()?;
            output(json!({"state":"advisory","verified":false,"context_revision":report.context_revision,"intent_context_revision":report.intent.context_revision,
                "target_revision":report.intent.target.revision(),"recipe":report.recipe,"ordered_oligos":ordered,"design_count":report.designs.len(),"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Propose original-source PCR primers for explicitly selected rows of a junction-edited overlap assembly, preserving all retained preparation DTOs exactly, including manual primers and pinned digests. Retained preparations may be incompatible: fresh full-recipe verification is required. Uses the declared trimmed source cores and literal inserted DNA to form the independent intended reference before designing primers. Added primer tails include the requested junction DNA plus next-core homology; full oligo metrics and existing filters are recomputed. Returns an advisory explicit-PCR recipe, no automatic application or verified cloning result. Reconstruct that recipe against the independently built edited target using propose_overlap_assembly."
    )]
    async fn design_selective_assembly_primers(
        &self,
        Parameters(args): Parameters<SelectiveAssemblyPrimersQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 actual original source bindings"));
        }
        let server = self.clone();
        self.compute(move||{
            let inputs={let store=server.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs{if input.id.trim().is_empty()||input.id.len()>128||!seen.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::selective_junction_edit::preflight(&sources,&args.request).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::selective_junction_edit::design(&inputs,&args.request).map_err(bad)?;
            let expected=collagene_core::selective_junction_edit::design_revision(&inputs,&args.request).map_err(bad)?;report.validate_binding(&expected).map_err(bad)?;
            let evidence=json!({"inputs":args.inputs,"report":report});
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","selective-assembly-primers.json")?;
            let ordered=report.recipe.fragments.iter().zip(&report.request.rows).filter_map(|(f,row)|match row {
                collagene_core::selective_junction_edit::Row::Retain(_) => None,
                collagene_core::selective_junction_edit::Row::DesignPcr{..} => Some(match &f.preparation {
                    collagene_core::assembly_recipe::Preparation::ExactPcr{forward,reverse, .. }=>Ok(json!({"role":f.role,"input_id":f.input_id,"reverse_fragment":f.reverse,"forward":forward,"reverse":reverse,"forward_ordered":forward.ordered_sequence(),"reverse_ordered":reverse.ordered_sequence()})),
                    _=>Err(bad("selected primer proposal lost its exact original-source PCR"))
                })
            }).collect::<Result<Vec<_>,_>>()?;
            output(json!({"state":"advisory","verified":false,"context_revision":report.context_revision,"intent_context_revision":report.intent.context_revision,
                "target_revision":report.intent.target.revision(),"recipe":report.recipe,"ordered_oligos":ordered,"design_count":ordered.len(),"retained_rows":report.retained_rows,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Propose explicitly selected original-source PCR and catalog-pinned Type IIS digest preparations for a junction-edited Golden Gate design. Builds independent source-core trim/insert intent before any PCR/digest. Actual tails distinguish fusion and inserted DNA, and fresh original-source PCR/cut/end geometry checks reject extra sites or nonunique binding. Retained rows, named oligos, catalog pins, treatments and repairs stay exact. This is an advisory ordered preparation/ligation recipe, not a cycling, domestication, fusion-optimization or experimental result. Build the edited intent independently and freshly check the returned complete recipe using propose_ligation_assembly."
    )]
    async fn design_golden_gate_junction_primers(
        &self,
        Parameters(args): Parameters<GoldenGateJunctionPrimersQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 actual original source bindings"));
        }
        let server = self.clone();
        self.compute(move||{
            let inputs={let store=server.store.lock().unwrap();let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs{if input.id.trim().is_empty()||input.id.len()>128||!seen.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::golden_gate_junction_edit::preflight(&sources,&args.request).map_err(bad)?;
                sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>()};
            let report=collagene_core::golden_gate_junction_edit::design(&inputs,&args.request).map_err(bad)?;
            let expected=collagene_core::golden_gate_junction_edit::design_revision(&inputs,&args.request).map_err(bad)?;report.validate_binding(&expected).map_err(bad)?;
            let evidence=json!({"inputs":args.inputs,"report":report});
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","golden-gate-junction-primers.json")?;
            let ordered=report.recipe.fragments.iter().zip(&report.request.rows).filter_map(|(fragment,row)|{let f=&fragment.source;match row {
                collagene_core::golden_gate_junction_edit::Row::Retain(_) => None,
                collagene_core::golden_gate_junction_edit::Row::DesignPcr{..} => Some(match &f.preparation {
                    collagene_core::assembly_recipe::Preparation::PcrCatalogDigest{forward,reverse,..}=>Ok(json!({"role":f.role,"input_id":f.input_id,"reverse_fragment":f.reverse,"forward":forward,"reverse":reverse,"forward_ordered":forward.ordered_sequence(),"reverse_ordered":reverse.ordered_sequence()})),
                    _=>Err(bad("selected Golden Gate proposal lost its original-source PCR/digest chain"))
                })
            }}).collect::<Result<Vec<_>,_>>()?;
            output(json!({"state":"advisory","verified":false,"context_revision":report.context_revision,"intent_context_revision":report.intent.context_revision,
                "target_revision":report.intent.target.revision(),"recipe":report.recipe,"ordered_oligos":ordered,"design_count":ordered.len(),"retained_rows":report.retained_rows,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Search bounded signed Golden Gate fusion-placement windows while retaining the original independently declared intended DNA and fixed manual preparations. Mathematical palindrome/end-distance constraints are not ligation fidelity or PCR feasibility. Returns complete candidate evidence as an artifact; explicitly select and materialize a candidate before fresh recipe verification."
    )]
    async fn search_golden_gate_fusions(
        &self,
        Parameters(args): Parameters<fusion_design_tools::SearchQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.fusion_search(args).await
    }
    #[tool(
        description = "Freshly materialize an explicitly selected source-bound Golden Gate fusion placement using original-template PCR and pinned Type IIS digests. Retains original intended DNA, records primer-added prefixes/suffixes, and preserves fixed preparations. Advisory only: verify the full returned ligation recipe against the original independent target before publication."
    )]
    async fn materialize_golden_gate_fusions(
        &self,
        Parameters(args): Parameters<fusion_design_tools::MaterializeQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.fusion_materialize(args).await
    }
    #[tool(
        description = "Estimate conditional ligation fidelity from an explicitly supplied, provenance-bound directed count matrix and selected fusion ends. Reports conditions and missing-data uncertainty; does not infer empirical data from exact complementarity or claim experimental yield. No external optimizer or runtime software is called."
    )]
    async fn score_ligation_fidelity(
        &self,
        Parameters(args): Parameters<fusion_design_tools::FidelityQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.fidelity_score(args).await
    }
    #[tool(
        description = "List exact versioned published ligation-frequency profiles, including formulation, conditions, source/conversion hashes, citation and separate CC-BY attribution. These data do not establish the actual preparation conditions or reaction efficiency."
    )]
    async fn list_ligation_profiles(
        &self,
        Parameters(args): Parameters<published_ligation_tools::ListQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.published_profiles_list(args).await
    }
    #[tool(
        description = "Score an explicit set of fusion ends using one exact published symmetric pair-frequency profile. Preserves literal published weights, zeros, palindromic/unknown cases and full conditions. No raw independent directional reads, measured yield or chemical verification is inferred."
    )]
    async fn score_published_ligation_profile(
        &self,
        Parameters(args): Parameters<published_ligation_tools::ScoreQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.published_profile_score(args).await
    }
    #[tool(
        description = "Freshly search source-bound Golden Gate fusion placements once and rank retained sets using one exact published profile prepared once. Requires an explicit conditional formulation assumption; never substitutes formulation from a generic enzyme name, selects a reaction winner, changes intended DNA or applies primers. Complete source/profile-bound evidence is returned as a bounded artifact."
    )]
    async fn rank_golden_gate_fusions(
        &self,
        Parameters(args): Parameters<published_ligation_tools::RankQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.published_fusion_rank(args).await
    }
    #[tool(
        description = "Freshly inspect TOPO insert PCR from one deliberately selected saved native material origin. Requires immutable source references, an explicit template role, a small explicit Run/Derived material selection and its exact retained native workspace. Replays the origin before actual new PCR; graph metadata or matching DNA cannot substitute for material authority. Select an explicit orientations, projection, independent-target Check, actions inspection or processed-actions task. Check resolves an exact immutable session target reference, never replaces it with a projected product. Returns bounded full typed evidence without saving any candidate or document; kit activation, transfers, processing and acceptance remain explicit assumptions."
    )]
    async fn inspect_saved_topo(
        &self,
        Parameters(args): Parameters<SavedTopoInspectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.len() > 32 {
            return Err(bad("bounded explicit TOPO input roles required"));
        }
        bounded_json_limit(&args.task, 2 * 1024 * 1024)?;
        if args
            .template
            .as_ref()
            .is_some_and(|t| t.role.trim().is_empty() || t.role.len() > 256)
        {
            return Err(bad("bounded saved PCR-template role required"));
        }
        if args.template.is_none() && args.task.preparation().selections()?.is_empty() {
            return Err(bad("choose an explicit saved graph or PCR-template origin"));
        }
        let server = self.clone();
        self.compute(move||{
            let selected=args.template.iter().map(|t|t.selection.clone()).collect::<Vec<_>>();
            let (preparation,mut materials)=server.hydrate_topo(args.task.preparation(),Some(&args.origin_workspace),&selected)?;
            let template=if let Some(t)=&args.template {
                if t.role.trim().is_empty()||t.role.len()>256||materials.len()!=1{return Err(bad("one exact saved PCR-template role required"));}
                Some(collagene_core::topo_run_preparation::BoundTemplate{role:t.role.clone(),origin:materials.remove(0)})
            }else{None};
            let (inputs,target)={let store=server.store.lock().unwrap();let mut ids=BTreeSet::new();let mut refs=vec![];let mut sources=vec![];
                for b in &args.inputs {if b.id.trim().is_empty()||b.id.len()>128||!ids.insert(b.id.as_str()){return Err(bad("invalid or duplicate saved TOPO role"));}refs.push(store.document(&b.document)?);sources.push(verification::TopoSource{role:b.id.clone(),source:verification::Source::Input(b.id.clone())});}
                let target=if let SavedTopoTask::Check{target,..}=&args.task {let m=store.document(target)?;collagene_core::topo_worker::target_preflight(m).map_err(bad)?;Some(m.clone())}else{None};
                let inputs=collagene_core::topo_run_preparation::source_inputs(&sources,&refs,&preparation,template.as_ref()).map_err(bad)?;
                (inputs,target)
            };
            let task=args.task.native(target,preparation)?;
            let context=collagene_core::topo_run_preparation::revision(&inputs,template.as_ref(),&task).map_err(bad)?;
            let report=collagene_core::topo_run_preparation::inspect(&inputs,template.as_ref(),&task).map_err(bad)?;
            report.validate_for(&inputs,template.as_ref(),&task,&context).map_err(bad)?;
            #[derive(Serialize)]struct Evidence<'a>{inputs:&'a[BoundInput],template:&'a Option<topo_graph_tools::BoundTemplate>,origin_workspace:&'a WorkspaceRef,task:&'a SavedTopoTask,inspection:&'a collagene_core::topo_run_preparation::Inspection}
            let evidence=Evidence{inputs:&args.inputs,template:&args.template,origin_workspace:&args.origin_workspace,task:&args.task,inspection:&report};
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","saved-topo-inspection.json")?;
            output(json!({"context_revision":context,"engine_version":collagene_core::topo_run_preparation::ENGINE_VERSION,"template":args.template,"origin_workspace":args.origin_workspace,"origin_incomplete":report.origin_incomplete,"origin_findings":report.origin_findings,"origin_findings_total":report.origin_findings_total,"origin_findings_truncated":report.origin_findings_truncated,"artifact_id":artifact,"verified":false,"authority":"Exact retained native graph and optional PCR-template origins freshly reconstructed; whole graph/prefix/byproducts retained. No duplicate insert synthesis, imported report authority, isolation or experimental success."}))
        }).await
    }
    #[tool(
        description = "List, freshly inspect, or explicitly capture a whole native TOPO graph in a pinned workspace. Source/event/record selections are small and session-attested; no supplied graph report or metadata becomes material. Capture returns a NEW immutable archive/receipt without changing heads or live application state. Retains partial-transfer history, activations, nicks, byproducts and inherited limitations. Retains the exact ordered transfer, flap-cleavage and nick-sealing transcript; no inventory, activation manufacture, isolation or experiment claim."
    )]
    async fn topo_graph(
        &self,
        Parameters(args): Parameters<topo_graph_tools::Query>,
    ) -> Result<CallToolResult, ErrorData> {
        self.topo_graph_dispatch(args).await
    }
    #[tool(
        description = "Compile explicit source-backed TA or directional CACC TOPO preparation declarations into the shared native recipe. Preserves actual ordered primers, vector strand geometry, terminal states and activation assumptions. Does not simulate activation or repair, add a tail silently, select transfers or invent a closed clone. Reuse inspect_topo_orientations and existing projection/target checks with the returned fields."
    )]
    async fn compile_topo_preparation(
        &self,
        Parameters(args): Parameters<published_ligation_tools::TopoQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.topo_compile_preparation(args).await
    }

    #[tool(
        description = "Freshly reconstruct source-bound chemical actions and prove a selected paired contour. Uses exact processing context/component/source seed and preserves the excised fragment, all parent/child graphs, created end chemistry and Unknown modifications. An annotated sequence view cannot restore private material; use the TopoActions plan operation with an independent target for full verification."
    )]
    async fn project_topo_actions(
        &self,
        Parameters(args): Parameters<topo_actions_tools::ProductQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.project_actions(args).await
    }
    #[tool(
        description = "Inspect source-bound restriction sites, explicit retained/protected ranges and coding-feature constraints for internal-site domestication. Uses exact catalog/feature/source pins and reports unsupported geometry or translation semantics rather than assuming a silent change. Reference inspection only; no physical stock is made."
    )]
    async fn inspect_domestication(
        &self,
        Parameters(args): Parameters<domestication_tools::DesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.domestication_inspect(args).await
    }
    #[tool(
        description = "Propose bounded synonymous source-coordinate substitutions that remove scoped unwanted restriction sites while preserving declared sites, protected bases and affected supported coding features. Positional minimum-edit constraints are inside the shared native search. Pruning remains explicit; a candidate is reference intent, not prepared DNA."
    )]
    async fn propose_domestication(
        &self,
        Parameters(args): Parameters<domestication_tools::DesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.domestication_propose(args).await
    }
    #[tool(
        description = "Independently check explicit source-bound domestication substitutions, preserving original annotation coordinates and qualifiers and rescanning affected proteins and the entire DNA for scoped unwanted/new sites. Retains the independently declared target and evidence atomically. This check confers no physical-preparation or experimental authority."
    )]
    async fn check_domestication(
        &self,
        Parameters(args): Parameters<domestication_tools::CheckQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.domestication_check(args).await
    }
    #[tool(
        description = "Materialize a literal checked domestication target using an explicit original-source route: inverse PCR on a circle, or whole-source overlap-extension PCR on linear/circular DNA. OE performs two actual original-template PCR branches, physical extension and outer PCR; circles then require explicit phosphorylation and ligation. All substitutions, coding/protected constraints and scoped sites are checked before planning. Atomically retains the independently frozen target, one canonical native plan and full route/oligo/search evidence for fresh verify_plan and HTML/PDF/order export. Unresolved PCR geometry retains a nonready diagnostic candidate; unsupported sources, settings and failed coverage are refused. No invented recoded stock or imported report authority."
    )]
    async fn prepare_domestication(
        &self,
        Parameters(args): Parameters<domestication_tools::PrepareQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        self.domestication_prepare(args).await
    }
    #[tool(
        description = "Build an unverified intended reference by explicitly replacing, inserting or deleting a single source interval. Uses only immutable source DNA and the requested edit, independently of any PCR or cloning simulator. Supports linear DNA and at most one circular turn with explicit origin shift. Atomically retains the intended target and full source-bound JSON. Does not mutate the source or claim physical material, successful cloning or observed DNA; design and verify inverse-PCR or another preparation separately against this frozen reference."
    )]
    async fn build_mutation_intent(
        &self,
        Parameters(args): Parameters<MutationIntentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let source = {
                let store = server.store.lock().unwrap();
                let source = store.document(&args.document)?;
                collagene_core::mutation_intent::source_preflight(source,&args.request).map_err(bad)?;
                source.clone()
            };
            let report = collagene_core::mutation_intent::build(&source,&args.request).map_err(bad)?;
            let expected = collagene_core::mutation_intent::intent_revision(&source,&args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a> { document:&'a DocumentRef, report:&'a collagene_core::mutation_intent::Report }
            let evidence = Evidence{document:&args.document,report:&report};
            let (artifact, target) = server.candidate_document_artifact_bounded(&evidence,Some(report.target.clone()),"mutation-intent.json",32*1024*1024)?;
            let target = target.ok_or_else(||bad("intended mutation reference was not retained"))?;
            output(json!({"state":"intended-reference","verified":false,"preparation_simulated":false,
                "document":args.document,"context_revision":report.context_revision,"request":report.request,
                "source_revision":report.source_revision,"target":{"id":target["id"],"revision":target["revision"]},
                "target_summary":target,"deleted_length":report.deleted.len(),"source_origin_shift":report.origin_shift,
                "warnings":report.warnings,"artifact_id":artifact,
                "authority":"Declared reference edit only. The source is unchanged; physical preparation, target agreement and experimental success require separate evidence."}))
        }).await
    }
    #[tool(
        description = "Plan bounded anchored PCR primers for explicitly ordered source intervals and orientations, using exact immutable DNA documents and an independently supplied target. Generates one-sided homology tails from the next actual source interval, never from target DNA, plus a practical overlap recipe and ordered oligos. Existing primer heuristic/conditions and Gibson/HiFi/InFusion advisory profile are explicit. Planner output is advisory and unverified, not a saved candidate or proof of amplification, specificity or target agreement. Pass its recipe and exact inputs/target to propose_overlap_assembly for fresh preparation and verification."
    )]
    async fn design_overlap_primers(
        &self,
        Parameters(args): Parameters<OverlapPrimerDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad(
                "overlap-primer design requires 1..128 actual source bindings",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target)={let store=server.store.lock().unwrap();let target=store.document(&args.target)?;
                let mut seen=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty()||bound.id.len()>128||!seen.insert(bound.id.as_str()) {return Err(bad("empty, oversized or duplicate template input ID"));}
                    sources.push((bound.id.as_str(),store.document(&bound.document)?));
                }
                collagene_core::overlap_primer_design::preflight(&sources,target,&args.request).map_err(bad)?;
                let inputs=sources.into_iter().map(|(id,molecule)|PlanInput{id:id.into(),revision:molecule.revision(),molecule:molecule.clone()}).collect::<Vec<_>>();
                (inputs,target.clone())};
            let report=collagene_core::overlap_primer_design::design(&inputs,&target,&args.request).map_err(bad)?;
            let expected=collagene_core::overlap_primer_design::design_revision(&inputs,&target,&args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            let ordered_oligos=report.recipe.fragments.iter().map(|fragment| {
                match &fragment.preparation {
                    collagene_core::assembly_recipe::Preparation::ExactPcr{forward,reverse, .. }=>Ok(json!({
                        "role":fragment.role,"input_id":fragment.input_id,"reverse_fragment":fragment.reverse,
                        "forward":{"primer":forward,"ordered_sequence":forward.ordered_sequence()},
                        "reverse":{"primer":reverse,"ordered_sequence":reverse.ordered_sequence()}})),
                    _=>Err(bad("planner did not retain explicit exact-PCR preparation"))
                }
            }).collect::<Result<Vec<_>,_>>()?;
            #[derive(Serialize)]
            struct Evidence<'a>{inputs:&'a [BoundInput],target:&'a DocumentRef,report:&'a collagene_core::overlap_primer_design::Report,interpretation:&'static str}
            let evidence=Evidence{inputs:&args.inputs,target:&args.target,report:&report,
                interpretation:"Advisory anchored primer heuristic and source-derived one-sided homology only. Independent target identity is retained but agreement, amplification and uniqueness require fresh compiler/verifier reconstruction; no source or candidate is promoted."};
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","overlap-primer-design.json")?;
            output(json!({"state":"advisory","verified":false,"inputs":args.inputs,"target":args.target,
                "context_revision":report.context_revision,"source_revisions":report.source_revisions,
                "recipe":report.recipe,"ordered_oligos":ordered_oligos,
                "design_count":report.designs.len(),"designs":report.designs.iter().take(32).map(|design|json!({
                    "source_revision":design.source_revision,"candidate_pair_count":design.pairs.len(),
                    "evaluated_pairs":design.evaluated_pairs,"rejected_pairs":design.rejected_pairs,
                    "warnings":design.warnings.iter().take(10).map(|w|w.chars().take(512).collect::<String>()).collect::<Vec<_>>(),
                    "warnings_truncated":design.warnings.len()>10||design.warnings.iter().take(10).any(|w|w.chars().count()>512)})).collect::<Vec<_>>(),
                "design_preview_truncated":report.designs.len()>32,
                "warnings":report.warnings.iter().take(100).map(|w|w.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "warnings_truncated":report.warnings.len()>100||report.warnings.iter().take(100).any(|w|w.chars().count()>1024),
                "artifact_id":artifact,"authority":"Advisory design only. Use the returned recipe with propose_overlap_assembly against the same exact inputs and independent target; method guidance and primer heuristics do not verify a cloning result."}))
        }).await
    }
    #[tool(
        description = "Propose primers for an entire explicitly edited overlap assembly from immutable original source documents and a typed independent contribution/junction/substitution edit intent. Unedited rows use exact PCR in their explicitly bound original or native strand frame; internally edited rows use actual two-branch overlap-extension PCR. Outgoing tails come from literal junction inserts plus the next edited source-core prefix, never a supplied product sequence. Freshly verifies the whole composed recipe against the separately built edit intent and retains failures/uncertainty. The single intent material-binding roster contains compact exact input/event Run/Derived selections, including actual Top/Bottom strands; all bindings require one exact origin_workspace. Closed recipes are hydrated jointly without replay, then the native proposal freshly reconstructs them. Returns a read-only advisory recipe, ordered oligos and bounded status summaries with full evidence; no source, target document, candidate or material status is installed. Retain the intended target separately with build_contribution_edit_intent using the same intent, then use propose_overlap_assembly and normal review/finalization."
    )]
    async fn design_edited_overlap_primers(
        &self,
        Parameters(args): Parameters<EditedOverlapPrimerDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 exact sources"));
        }
        bounded_json_limit(&args.request, 2 * 1024 * 1024)?;
        if !args.request.intent.base.base.material_bindings.is_empty()
            && args.origin_workspace.is_none()
        {
            return Err(bad(
                "material source bindings require one exact origin_workspace",
            ));
        }
        let server = self.clone();
        self.compute(move||{
            let frame=args.origin_workspace.as_ref().map(|wr|server.restore_workspace(wr)).transpose()?;
            let request=collagene_core::assembly_edited_primer_design::Request {
                intent:args.request.intent.native(frame.as_ref())?,overlap:args.request.overlap,
                options:args.request.options.clone(),profile_options:args.request.profile_options.clone(),finishing:args.request.finishing.clone(),
            };
            let inputs = {
                let store = server.store.lock().unwrap();
                let mut ids = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for input in &args.inputs {
                    if input.id.trim().is_empty() || input.id.len() > 128 || !ids.insert(input.id.as_str()) {
                        return Err(bad("invalid or duplicate original source ID"));
                    }
                    sources.push((input.id.as_str(), store.document(&input.document)?));
                }
                // Native source_preflight jointly attests the single input/event
                // roster against borrowed immutable session documents before cloning.
                collagene_core::assembly_edited_primer_design::source_preflight(&sources, &request, frame.as_ref()).map_err(bad)?;
                sources.into_iter().map(|(id,molecule)|PlanInput { id:id.into(),revision:molecule.revision(),molecule:molecule.clone() }).collect::<Vec<_>>()
            };
            let expected = collagene_core::assembly_edited_primer_design::design_revision(&inputs, &request, frame.as_ref()).map_err(bad)?;
            let report = collagene_core::assembly_edited_primer_design::design(&inputs, &request, frame.as_ref()).map_err(bad)?;
            report.validate_sources(&inputs, &expected, frame.as_ref()).map_err(bad)?;
            let mut ordered = Vec::new();
            for fragment in &report.recipe.fragments {
                use collagene_core::assembly_recipe::{Preparation,RunPreparation};
                let pairs:Vec<(&str,&Primer,&Primer,collagene_core::chemistry::Strand)> = match &fragment.preparation {
                    Preparation::ExactPcr {forward,reverse,template_strand} => vec![("exact_pcr",forward,reverse,*template_strand)],
                    Preparation::OverlapExtensionPcr {request} => vec![("left_pcr",&request.left.forward,&request.left.reverse,request.left.template_strand),("right_pcr",&request.right.forward,&request.right.reverse,request.right.template_strand),("outer_pcr",&request.outer.forward,&request.outer.reverse,request.outer.template_strand)],
                    Preparation::FromRun {preparation:RunPreparation::ExactPcr {forward,reverse,template_strand},..} | Preparation::FromDerived {preparation:RunPreparation::ExactPcr {forward,reverse,template_strand},..} => vec![("exact_pcr",forward,reverse,*template_strand)],
                    Preparation::FromRun {preparation:RunPreparation::OverlapExtensionPcr {request},..} | Preparation::FromDerived {preparation:RunPreparation::OverlapExtensionPcr {request},..} => vec![("left_pcr",&request.left.forward,&request.left.reverse,request.left.template_strand),("right_pcr",&request.right.forward,&request.right.reverse,request.right.template_strand),("outer_pcr",&request.outer.forward,&request.outer.reverse,request.outer.template_strand)],
                    _ => return Err(bad("edited assembly primer proposal lost its explicit exact/OE PCR preparation")),
                };
                for (stage,forward,reverse,template_strand) in pairs {
                    for (direction,primer) in [("forward",forward),("reverse",reverse)] {
                        ordered.push(json!({"role":fragment.role,"input_id":fragment.input_id,"reverse_fragment":fragment.reverse,"stage":stage,"template_strand":template_strand,"direction":direction,"primer":primer,"ordered_sequence":primer.ordered_sequence()}));
                    }
                }
            }
            #[derive(Serialize)]
            struct Evidence<'a> { inputs:&'a [BoundInput],origin_workspace:&'a Option<WorkspaceRef>,request:&'a EditedOverlapPrimerRequest,report:&'a collagene_core::assembly_edited_primer_design::Report,interpretation:&'static str }
            let evidence = Evidence { inputs:&args.inputs,origin_workspace:&args.origin_workspace,request:&args.request,report:&report,
                interpretation:"Read-only whole-construct primer proposal. Independent intended DNA is built only from exact original source contributions and explicit edit declarations. Primer tails and selected native source geometry never select or replace that intended target; whole-recipe reconstruction remains diagnostic until separate review/finalization." };
            let artifact = server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","edited-overlap-primer-design.json")?;
            output(json!({
                "state":"advisory","verified":false,"inputs":args.inputs,"origin_workspace":args.origin_workspace,
                "engine_version":report.engine_version,"context_revision":report.context_revision,"source_revisions":report.source_revisions,
                "intent_context_revision":report.intent.context_revision,
                "target_summary":{"revision":report.intent.target.revision(),"length":report.intent.target.sequence.len(),"topology":report.intent.target.topology},
                "recipe":report.recipe,"ordered_oligos":ordered,"tails":report.tails,
                "verification":{"status":report.preview.verification.status,"plan_revision":report.preview.verification.plan_revision,
                    "target_revision":report.preview.verification.target_revision,"product_count":report.preview.verification.products.len(),
                    "findings":report.preview.verification.findings.iter().take(50).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                    "findings_truncated":report.preview.verification.findings.len()>50||report.preview.verification.findings.iter().take(50).any(|s|s.chars().count()>1024)},
                "rows":report.rows.iter().take(32).map(|row|json!({"role":row.role,"input_id":row.input_id,"input_revision":row.input_revision,
                    "edit_count":row.edits.len(),"design_count":row.designs.len(),"template_strand":row.template_strand,"search":row.search})).collect::<Vec<_>>(),
                "warnings":report.warnings.iter().take(100).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "warnings_truncated":report.warnings.len()>100||report.warnings.iter().take(100).any(|s|s.chars().count()>1024),
                "artifact_id":artifact,
                "authority":"Read-only advisory recipe. Retain the independently declared target with build_contribution_edit_intent using the same inputs/request.intent, then freshly verify this recipe against that target through propose_overlap_assembly. No source, target document, candidate or physical material is installed."
            }))
        }).await
    }
    #[tool(
        description = "Propose six ordered oligos for a source-bound overlap-extension PCR preparation from an immutable original or exact saved Run/Derived source, an explicit stored-strand interval/orientation and independently declared literal edits. Saved sources require an exact material binding and owning origin_workspace; actual strand frames are reconstructed natively. Uses anchored primer heuristics on unchanged native source domains, then freshly reconstructs two PCR branches, selected reciprocal extension and outer-primer PCR against the separately constructed edit intent. Returns an advisory preparation, bounded oligo/status summaries and full source/request/intent/reconstruction evidence. Prepared DNA stays in stored-source orientation; apply request.fragment.reverse at subsequent assembly. No intended product is accepted as a source, no stock is changed, and no candidate, target document or material status is installed. Physical-state metadata and historical reports cannot bootstrap an original PCR template. Reuse the returned preparation in the ordinary assembly recipe and verify the full independent target."
    )]
    async fn design_overlap_extension_primers(
        &self,
        Parameters(args): Parameters<OverlapExtensionPrimerDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let workspace = args.origin_workspace.as_ref().map(|reference| server.restore_workspace(reference)).transpose()?;
            let source = {
                let store = server.store.lock().unwrap();
                let source = store.document(&args.document)?;
                collagene_core::overlap_extension_primer_design::preflight(source, &args.request, workspace.as_ref()).map_err(bad)?;
                source.clone()
            };
            let report = collagene_core::overlap_extension_primer_design::design(&source, &args.request, workspace.as_ref()).map_err(bad)?;
            let expected = collagene_core::overlap_extension_primer_design::design_revision(&source, &args.request, workspace.as_ref()).map_err(bad)?;
            report.validate_sources(&source, &expected, workspace.as_ref()).map_err(bad)?;
            let mut ordered_oligos = Vec::with_capacity(6);
            for (stage, pair) in [("left_pcr", &report.preparation.left), ("right_pcr", &report.preparation.right), ("outer_pcr", &report.preparation.outer)] {
                for (direction, primer) in [("forward", &pair.forward), ("reverse", &pair.reverse)] {
                    ordered_oligos.push(json!({"stage":stage,"template_strand":pair.template_strand,"direction":direction,"primer":primer,"ordered_sequence":primer.ordered_sequence()}));
                }
            }
            #[derive(Serialize)]
            struct Evidence<'a> {
                document: &'a DocumentRef,
                report: &'a collagene_core::overlap_extension_primer_design::Report,
                interpretation: &'static str,
            }
            let evidence = Evidence { document: &args.document, report: &report,
                interpretation: "Read-only native primer proposal from explicit immutable source coordinates and independent literal edits. Expected DNA comes from edit intent, not simulated products. The retained local reconstruction is diagnostic; no stock, candidate, target document or material authority is installed. Reuse the preparation in a full independently targeted assembly before review/finalization." };
            let artifact = server.text_artifact(bounded_json_limit(&evidence, 32 * 1024 * 1024)?, "application/json", "overlap-extension-primer-design.json")?;
            output(json!({
                "state":"advisory", "verified":false, "document":args.document,
                "context_revision":report.context_revision,"source_revision":report.source_revision,"engine_version":report.engine_version,
                "fragment":report.request.fragment,"reverse_at_assembly":report.request.fragment.reverse,
                "preparation":report.preparation,"ordered_oligos":ordered_oligos,
                "expected":{"revision":report.expected.revision(),"length":report.expected.sequence.len(),"topology":report.expected.topology,"orientation":"stored_source"},
                "verification":{"status":report.verification.status,"plan_revision":report.verification.plan_revision,"target_revision":report.verification.target_revision,
                    "product_count":report.verification.products.len(),"products":report.verification.products.iter().take(4).map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})).collect::<Vec<_>>(),
                    "findings":report.verification.findings.iter().take(30).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                    "findings_truncated":report.verification.findings.len()>30||report.verification.findings.iter().take(30).any(|s|s.chars().count()>1024)},
                "design_count":report.designs.len(),
                "search":report.search.iter().take(3).map(|search|json!({
                    "call":search.call,"candidates":search.candidates,"examined":search.examined,"rejected":search.rejected,
                    "selected_index":search.selected_index,"capped":search.capped,"unique_geometry":search.unique_geometry,
                    "notes":search.notes.iter().take(10).map(|s|s.chars().take(512).collect::<String>()).collect::<Vec<_>>(),
                    "notes_truncated":search.notes.len()>10||search.notes.iter().take(10).any(|s|s.chars().count()>512)
                })).collect::<Vec<_>>(),
                "warnings":report.warnings.iter().take(100).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "warnings_truncated":report.warnings.len()>100||report.warnings.iter().take(100).any(|s|s.chars().count()>1024),
                "artifact_id":artifact,
                "authority":"Advisory source-bound proposal only. The source and independent edits remain unchanged; use its explicit preparation in a full intended assembly for fresh native verification and review."
            }))
        }).await
    }
    #[tool(
        description = "Propose Type IIS site-bearing primers from explicit original source intervals, orientations, domain lengths and right fusion ends. Supports actual BsaI/BsmBI/BbsI/SapI geometry, including three-base SapI ends. Clamps and cleavage spacers are explicit DNA; target DNA does not select primers. Returns original-source PCR followed by pinned catalog digestion preparations and inspectable oligo segments. Advisory proposals do not install a candidate or establish experimental success or finite-mixture uniqueness. Pass fragments to propose_ligation_assembly for fresh reconstruction against the same independent target."
    )]
    async fn design_golden_gate_primers(
        &self,
        Parameters(args): Parameters<GoldenGatePrimerDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad(
                "Type IIS primer design requires 1..128 original source bindings",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs, target) = {
                let store = server.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty() || bound.id.len() > 128 || !seen.insert(bound.id.as_str()) { return Err(bad("empty, oversized or duplicate source ID")); }
                    sources.push((bound.id.as_str(), store.document(&bound.document)?));
                }
                collagene_core::golden_gate_primer_design::source_preflight(&sources, target, &args.request).map_err(bad)?;
                let inputs = sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>();
                (inputs,target.clone())
            };
            let report = collagene_core::golden_gate_primer_design::design(&inputs,&target,&args.request).map_err(bad)?;
            let expected = collagene_core::golden_gate_primer_design::preview_revision(&inputs,&target,&args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a> { inputs:&'a [BoundInput], target:&'a DocumentRef, report:&'a collagene_core::golden_gate_primer_design::Report }
            let evidence = Evidence{inputs:&args.inputs,target:&args.target,report:&report};
            let artifact = server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","golden-gate-primer-design.json")?;
            output(json!({"state":"advisory","verified":false,"inputs":args.inputs,"target":args.target,
                "context_revision":report.context_revision,"request":report.request,"fragments":report.fragments,
                "rows":report.rows.iter().map(|row|json!({
                    "role":row.role,"input_id":row.input_id,"input_revision":row.input_revision,
                    "forward":row.forward,"reverse":row.reverse,"forward_tail":row.forward_tail,"reverse_tail":row.reverse_tail,
                    "selection":row.selection,"fragment_index":row.fragment_index,
                    "oriented_length":row.oriented_top.len(),"oriented_top_preview":row.oriented_top.chars().take(160).collect::<String>(),
                    "oriented_top_truncated":row.oriented_top.len()>160,"left":row.left,"right":row.right,
                })).collect::<Vec<_>>(),"warnings":report.warnings,"artifact_id":artifact,
                "authority":"Explicit source-derived oligo proposals and preparation requests. Fresh ordered ligation verification against the independent target is separate; no claim of reaction yield or exhaustive mixture uniqueness."}))
        }).await
    }
    #[tool(
        description = "Preview an explicit physical homology-assembly recipe from original immutable sources and a separately bound intended target. Reconstructs actual preparations, oriented source strands, declared5-prime resection, complementary fill and sealing under an explicit finite policy; a vendor Gibson/HiFi/InFusion label supplies advisory sequence guidance only. Recipe overlaps/topology must exactly match the physical policy, never be silently overwritten. Retains products, rejected/ambiguous chemistry and bounded-search findings. Every FromRun role requires its exact native workspace and fresh closed-origin replay. Returns complete bounded plan/physical/junction evidence and concise summaries; it saves no candidate, source, target or product document. To keep a design, propose the retained typed operations against these same inputs/target and use normal verification/review."
    )]
    async fn preview_homology_recipe(
        &self,
        Parameters(args): Parameters<HomologyRecipePreviewQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty()
            || args.inputs.len() > collagene_core::assembly_recipe::MAX_ROWS
            || args.origin_workspaces.len() > collagene_core::assembly_recipe::MAX_ROWS
        {
            return Err(bad(
                "homology recipe sources/workspace bindings exceed their bounded row scope",
            ));
        }
        let origins = args
            .recipe
            .fragments
            .iter()
            .filter(|fragment| {
                matches!(
                    fragment.preparation,
                    collagene_core::assembly_recipe::Preparation::FromRun { .. }
                )
            })
            .count();
        if origins != args.origin_workspaces.len() {
            return Err(bad(
                "supply exactly one native workspace for each saved-source role",
            ));
        }
        let mut seen = BTreeSet::new();
        for binding in &args.origin_workspaces {
            if binding.role.is_empty()
                || binding.role.len() > 128
                || !seen.insert(binding.role.as_str())
                || !args.recipe.fragments.iter().any(|fragment| {
                    fragment.role == binding.role
                        && matches!(
                            fragment.preparation,
                            collagene_core::assembly_recipe::Preparation::FromRun { .. }
                        )
                })
            {
                return Err(bad(
                    "homology origin workspace role is invalid, duplicated or unused",
                ));
            }
        }
        let server = self.clone();
        self.compute(move||{
            let(inputs,target)={let store=server.store.lock().unwrap();let target=store.document(&args.target)?;let mut ids=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs {if input.id.trim().is_empty()||input.id.len()>128||!ids.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::assembly_recipe::homology_source_preflight(&sources,target,&args.recipe,&args.policy).map_err(bad)?;
                (sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>(),target.clone())};
            for fragment in &args.recipe.fragments {
                if let collagene_core::assembly_recipe::Preparation::FromRun{origin,..}=&fragment.preparation {
                    let workspace=args.origin_workspaces.iter().find(|binding|binding.role==fragment.role).ok_or_else(||bad("saved homology source has no exact workspace binding"))?;
                    let source=inputs.iter().find(|input|input.id==fragment.input_id).ok_or_else(||bad("saved homology source input unavailable"))?;
                    server.attest_run_source(&workspace.workspace,origin,&source.molecule)?;
                }
            }
            let preview=collagene_core::assembly_recipe::preview_homology(&inputs,&target,&args.recipe,&args.policy).map_err(bad)?;
            let context=collagene_core::assembly_recipe::homology_preview_revision(&inputs,&target,&args.recipe,&args.policy).map_err(bad)?;
            preview.validate_binding(&context).map_err(bad)?;
            let diagnostics=collagene_core::assembly_diagnostics::overlap(&preview).map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a>{inputs:&'a[BoundInput],target:&'a DocumentRef,recipe:&'a collagene_core::assembly_recipe::Recipe,policy:&'a collagene_core::homology_assembly::Options,origin_workspaces:&'a[MaterialOriginWorkspace],preview:&'a collagene_core::assembly_recipe::Preview,diagnostics:&'a collagene_core::assembly_diagnostics::Report}
            let evidence=Evidence{inputs:&args.inputs,target:&args.target,recipe:&args.recipe,policy:&args.policy,origin_workspaces:&args.origin_workspaces,preview:&preview,diagnostics:&diagnostics};
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","physical-homology-recipe.json")?;
            let product=(preview.verification.products.len()==preview.compiled.plan.steps.len()).then(||preview.verification.products.last()).flatten();
            let physical=preview.verification.physical_processing_reports.iter().find_map(|step|match &step.evidence{verification::PhysicalProcessingEvidence::HomologyAssembly(report) if step.step_index==preview.compiled.assembly_step=>Some(report),_=>None});
            output(json!({"state":"preview","verified":false,"candidate_saved":false,"inputs":args.inputs,"target":args.target,
                "context_revision":context,"plan_revision":preview.verification.plan_revision,"target_revision":preview.verification.target_revision,
                "engine_version":preview.verification.engine_version,"status":preview.verification.status,"preparation_complete":preview.preparation_complete,
                "product":product.map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})),
                "physical":physical.map(|report|json!({"engine_version":report.engine_version,"paired_product_available":report.product.is_some(),"reusable_material_available":report.product.as_ref().is_some_and(|product|product.nicks.is_empty()),"terminal_search_complete":report.terminal_search_complete,"incomplete":report.incomplete,"comparisons":report.comparisons,"source_geometry_revisions":report.source_geometry_revisions})),
                "findings":preview.verification.findings.iter().take(100).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "findings_truncated":preview.verification.findings.len()>100||preview.verification.findings.iter().take(100).any(|s|s.chars().count()>1024),
                "policy":args.policy,"diagnostics":diagnostics,"artifact_id":artifact,"plan_artifact_path":"preview.compiled.plan",
                "authority":"Fresh explicit physical-model preview and independent-target comparison only. Paired reference DNA can remain diagnostic when joining termini cannot seal; it does not grant reusable material. Eligible material requires fresh native reconstruction under the same policy and explicit saved-origin replay for reuse. This response saves no candidate or product and makes no experimental-success, isolation or mixture-uniqueness claim."
            }))
        }).await
    }
    #[tool(
        description = "Preview explicit 3-prime resection and homology annealing from original immutable sources and an independent target. Requires actual preparations, per-source left/right strand cuts, explicit junction lengths and bounded pairing search. Preserves remaining strands, source coordinates, unfilled gaps and unsealed nicks. The joined reference always remains incomplete and cannot confer reusable assembled material; no polymerase, ligase or host processing is inferred. FromRun roles require exact native workspace attestations and fresh closed-origin replay. Returns complete bounded typed plan, physical and junction evidence without saving any candidate or document."
    )]
    async fn preview_annealing_recipe(
        &self,
        Parameters(args): Parameters<AnnealingRecipePreviewQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty()
            || args.inputs.len() > collagene_core::assembly_recipe::MAX_ROWS
            || args.origin_workspaces.len() > collagene_core::assembly_recipe::MAX_ROWS
        {
            return Err(bad(
                "homology recipe sources/workspace bindings exceed their bounded row scope",
            ));
        }
        let origins = args
            .recipe
            .fragments
            .iter()
            .filter(|fragment| {
                matches!(
                    fragment.preparation,
                    collagene_core::assembly_recipe::Preparation::FromRun { .. }
                )
            })
            .count();
        if origins != args.origin_workspaces.len() {
            return Err(bad(
                "supply exactly one native workspace for each saved-source role",
            ));
        }
        let mut seen = BTreeSet::new();
        for binding in &args.origin_workspaces {
            if binding.role.is_empty()
                || binding.role.len() > 128
                || !seen.insert(binding.role.as_str())
                || !args.recipe.fragments.iter().any(|fragment| {
                    fragment.role == binding.role
                        && matches!(
                            fragment.preparation,
                            collagene_core::assembly_recipe::Preparation::FromRun { .. }
                        )
                })
            {
                return Err(bad(
                    "homology origin workspace role is invalid, duplicated or unused",
                ));
            }
        }
        let server = self.clone();
        self.compute(move||{
            let(inputs,target)={let store=server.store.lock().unwrap();let target=store.document(&args.target)?;let mut ids=BTreeSet::new();let mut sources=Vec::with_capacity(args.inputs.len());
                for input in &args.inputs {if input.id.trim().is_empty()||input.id.len()>128||!ids.insert(input.id.as_str()){return Err(bad("invalid or duplicate original source ID"));}sources.push((input.id.as_str(),store.document(&input.document)?));}
                collagene_core::assembly_recipe::annealing_source_preflight(&sources,target,&args.recipe,&args.policy).map_err(bad)?;
                (sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>(),target.clone())};
            for fragment in &args.recipe.fragments {
                if let collagene_core::assembly_recipe::Preparation::FromRun{origin,..}=&fragment.preparation {
                    let workspace=args.origin_workspaces.iter().find(|binding|binding.role==fragment.role).ok_or_else(||bad("saved homology source has no exact workspace binding"))?;
                    let source=inputs.iter().find(|input|input.id==fragment.input_id).ok_or_else(||bad("saved homology source input unavailable"))?;
                    server.attest_run_source(&workspace.workspace,origin,&source.molecule)?;
                }
            }
            let preview=collagene_core::assembly_recipe::preview_annealing(&inputs,&target,&args.recipe,&args.policy).map_err(bad)?;
            let context=collagene_core::assembly_recipe::annealing_preview_revision(&inputs,&target,&args.recipe,&args.policy).map_err(bad)?;
            preview.validate_binding(&context).map_err(bad)?;
            let diagnostics=collagene_core::assembly_diagnostics::overlap(&preview).map_err(bad)?;
            #[derive(Serialize)]
            struct Evidence<'a>{inputs:&'a[BoundInput],target:&'a DocumentRef,recipe:&'a collagene_core::assembly_recipe::Recipe,policy:&'a collagene_core::homology_annealing::Options,origin_workspaces:&'a[MaterialOriginWorkspace],preview:&'a collagene_core::assembly_recipe::Preview,diagnostics:&'a collagene_core::assembly_diagnostics::Report}
            let evidence=Evidence{inputs:&args.inputs,target:&args.target,recipe:&args.recipe,policy:&args.policy,origin_workspaces:&args.origin_workspaces,preview:&preview,diagnostics:&diagnostics};
            let artifact=server.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","three-prime-annealing-recipe.json")?;
            let product=(preview.verification.products.len()==preview.compiled.plan.steps.len()).then(||preview.verification.products.last()).flatten();
            let physical=preview.verification.physical_processing_reports.iter().find_map(|step|match &step.evidence{verification::PhysicalProcessingEvidence::HomologyAnnealing(report) if step.step_index==preview.compiled.assembly_step=>Some(report),_=>None});
            output(json!({"state":"preview","verified":false,"candidate_saved":false,"inputs":args.inputs,"target":args.target,
                "context_revision":context,"plan_revision":preview.verification.plan_revision,"target_revision":preview.verification.target_revision,
                "engine_version":preview.verification.engine_version,"status":preview.verification.status,"preparation_complete":preview.preparation_complete,
                "product":product.map(|m|json!({"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology})),
                "physical":physical.map(|report|json!({"engine_version":report.engine_version,"paired_diagnostic_available":report.paired_diagnostic.is_some(),"reusable_material_available":false,"terminal_search_complete":report.terminal_search_complete,"incomplete":report.incomplete,"comparisons":report.comparisons,"source_geometry_revisions":report.source_geometry_revisions})),
                "findings":preview.verification.findings.iter().take(100).map(|s|s.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "findings_truncated":preview.verification.findings.len()>100||preview.verification.findings.iter().take(100).any(|s|s.chars().count()>1024),
                "policy":args.policy,"diagnostics":diagnostics,"artifact_id":artifact,"plan_artifact_path":"preview.compiled.plan",
                "authority":"Fresh 3-prime resection and homology-annealing preview against an independent target. The joined reference remains incomplete and documentary; all gaps and nicks remain, no fill, sealing or host activity is inferred, and no reusable assembled physical material is created. This response saves no candidate or product and makes no experimental-success, isolation or mixture-uniqueness claim."
            }))
        }).await
    }
    #[tool(
        description = "Compile an explicit ordered overlap cloning recipe from immutable source document revisions and an independent target. Each fragment declares WholeLinear, an explicit DeclaredDuplex original-source terminal assumption, complete Digest or exact CatalogDigest selection with chosen fragment/optional repair, or ExactPcr with actual ordered oligos, or PcrCatalogDigest chaining actual original-source PCR into a pinned catalog digest; orientations, homology lengths, topology and published-method advisory profile are explicit. Freshly reconstructs preparations and verifies the resulting typed plan. Returns additive typed diagnostics locating failed/blocked preparation rows, junction evidence and final target mismatch, without parsing findings. Atomically saves an immutable candidate and full preparation/junction/verification JSON; no guessed tails, automatic winner, source mutation or laboratory guarantee. Reuse verify_plan and HTML/PDF exports with its plan_revision. Explicit MismatchPcr retains both ordered primers and its declared mismatch/3′-exact/search policy. OverlapExtensionPcr retains two actual PCR branches from the same source, selected physical overlap extension and separate outer-primer PCR; every ordered oligo and terminal-pair search policy is explicit. OverlapExtensionPcrCatalogDigest appends exact catalog-pinned cleavage of that actual outer-PCR product, explicit fragment index and optional repair; inspect_overlap_extension_digest exposes the same real product without requiring a target. FromRun preserves a deliberately selected saved origin before downstream preparation. Mutation intent and the target remain independently declared; no mutation is selected to fit a simulated product."
    )]
    async fn propose_overlap_assembly(
        &self,
        Parameters(args): Parameters<OverlapAssemblyProposal>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad("overlap recipe requires 1..128 explicit source inputs"));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target) = {
                let store = server.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty() || bound.id.len()>128 || !seen.insert(bound.id.as_str()) {
                        return Err(bad("empty, oversized or duplicate recipe input ID"));
                    }
                    sources.push((bound.id.as_str(),store.document(&bound.document)?));
                }
                collagene_core::assembly_recipe::source_preflight(&sources,target,&args.recipe).map_err(bad)?;
                let inputs = sources.into_iter().map(|(id,molecule)|PlanInput {
                    id:id.into(),revision:molecule.revision(),molecule:molecule.clone()
                }).collect::<Vec<_>>();
                (inputs,target.clone())
            };
            let preview = collagene_core::assembly_recipe::preview(&inputs,&target,&args.recipe).map_err(bad)?;
            let expected = collagene_core::assembly_recipe::preview_revision(&inputs,&target,&args.recipe).map_err(bad)?;
            preview.validate_binding(&expected).map_err(bad)?;
            // Typed locations come from the checked execution prefix, never finding strings.
            let diagnostics = collagene_core::assembly_diagnostics::overlap(&preview).map_err(bad)?;
            if diagnostics.context_revision != expected { return Err(bad("assembly diagnostics belong to another preview context")); }
            // Preserve all existing top-level preview fields, adding diagnostics only.
            #[derive(Serialize)]
            struct Evidence<'a> {
                #[serde(flatten)]
                preview: &'a collagene_core::assembly_recipe::Preview,
                diagnostics: &'a collagene_core::assembly_diagnostics::Report,
            }
            // Serialize all outputs and capacity-check before either is published.
            let text = bounded_json_limit(&Evidence { preview: &preview, diagnostics: &diagnostics },32*1024*1024)?;
            let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
            let revision = verification::plan_revision(&preview.compiled.plan);
            if revision.is_empty() {return Err(bad("compiled plan has no bounded canonical revision"));}
            let plan_bytes = bounded_json(&preview.compiled.plan)?.len().checked_add(bounded_json(&target)?.len()).ok_or_else(||bad("candidate storage accounting overflow"))?;
            let result = json!({"plan_revision":revision,"state":"candidate","verified":false,
                "target":args.target,"inputs":args.inputs,"context_revision":preview.context_revision,
                "preparation_complete":preview.preparation_complete,"diagnostics":diagnostics,"assembly_step":preview.compiled.assembly_step,
                "verification":{"status":preview.verification.status,"engine_version":preview.verification.engine_version,
                    "findings":preview.verification.findings.iter().take(100).collect::<Vec<_>>(),
                    "findings_truncated":preview.verification.findings.len()>100,
                    "assumptions":preview.verification.assumptions.iter().take(100).collect::<Vec<_>>(),
                    "assumptions_truncated":preview.verification.assumptions.len()>100,
                    "exhaustive_alternative_search":preview.verification.exhaustive_alternative_search},
                "artifact_id":artifact_id,
                "authority":"Fresh deterministic preparation and declared ordered-plan evidence; saving remains a candidate, method guidance is advisory, and sequence agreement is not laboratory success or exhaustive assembly uniqueness. Full rows, prepared products and junction evidence are in the artifact."});
            let response = output(result)?;
            let mut store=server.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&revision);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(existing)=store.artifacts.get(&artifact_id) {
                if !matches!(existing.data,ArtifactData::Text(_)) || existing.mime_type!="application/json" || existing.filename!="overlap-assembly-recipe.json" {
                    return Err(bad("assembly evidence already has another artifact representation"));
                }
            }
            if (new_plan && store.plans.len()>=MAX_ITEMS) || (new_artifact && store.artifacts.len()>=MAX_ITEMS) {
                return Err(bad("candidate/artifact item limit; no assembly outputs created"));
            }
            let bytes=(if new_plan {plan_bytes}else{0}).checked_add(if new_artifact{text.len()}else{0}).ok_or_else(||bad("session storage accounting overflow"))?;
            store.reserve(bytes)?;
            if new_plan {store.plans.insert(revision,StoredPlan{plan:preview.compiled.plan,target});}
            if new_artifact {store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"overlap-assembly-recipe.json"});}
            Ok(response)
        }).await
    }
    #[tool(
        description = "Compile an explicit ordered ligation cloning recipe from immutable DNA source revisions and an independent target. Each fragment declares WholeLinear, explicit DeclaredDuplex terminal assumptions on an original intact linear DNA input, complete Digest or catalog-pinned CatalogDigest with optional end repair, or actual ExactPcr preparation, explicit terminal treatment and orientation; circular and allow_nicked policies are explicit. Freshly reconstructs physical strands, oriented ends, junctions and retained nicks through the deterministic verifier. Returns additive typed row/step/junction diagnostics without treating available endpoints as proof of compatibility. Atomically stores a candidate and full preparation/physical/verification JSON, including failed drafts. No metadata-derived chemistry authority, automatic treatment/repair, experimental guarantee or source mutation. Reuse verify_plan and HTML/PDF exports with plan_revision. Explicit MismatchPcr retains both ordered primers and its declared mismatch/3′-exact/search policy. OverlapExtensionPcr retains two actual PCR branches from the same source, selected physical overlap extension and separate outer-primer PCR; every ordered oligo and terminal-pair search policy is explicit. OverlapExtensionPcrCatalogDigest appends exact catalog-pinned cleavage of that actual outer-PCR product, explicit fragment index and optional repair; inspect_overlap_extension_digest exposes the same real product without requiring a target. FromRun preserves a deliberately selected saved origin before downstream preparation. Mutation intent and the target remain independently declared; no mutation is selected to fit a simulated product."
    )]
    async fn propose_ligation_assembly(
        &self,
        Parameters(args): Parameters<LigationAssemblyProposal>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad(
                "ligation recipe requires 1..128 explicit source inputs",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target) = {
                let store = server.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty() || bound.id.len()>128 || !seen.insert(bound.id.as_str()) {
                        return Err(bad("empty, oversized or duplicate recipe input ID"));
                    }
                    sources.push((bound.id.as_str(),store.document(&bound.document)?));
                }
                collagene_core::ligation_recipe::source_preflight(&sources,target,&args.recipe).map_err(bad)?;
                let inputs = sources.into_iter().map(|(id,molecule)|PlanInput {
                    id:id.into(),revision:molecule.revision(),molecule:molecule.clone()
                }).collect::<Vec<_>>();
                (inputs,target.clone())
            };
            let preview = collagene_core::ligation_recipe::preview(&inputs,&target,&args.recipe, collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete).map_err(bad)?;
            let expected = collagene_core::ligation_recipe::preview_revision(&inputs,&target,&args.recipe, collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete).map_err(bad)?;
            preview.validate_binding(&expected).map_err(bad)?;
            // Typed locations come from the checked execution prefix, never finding strings.
            let diagnostics = collagene_core::assembly_diagnostics::ligation(&preview).map_err(bad)?;
            if diagnostics.context_revision != expected { return Err(bad("assembly diagnostics belong to another preview context")); }
            // Preserve all existing top-level preview fields, adding diagnostics only.
            #[derive(Serialize)]
            struct Evidence<'a> {
                #[serde(flatten)]
                preview: &'a collagene_core::ligation_recipe::Preview,
                diagnostics: &'a collagene_core::assembly_diagnostics::Report,
            }
            // Serialize all outputs and capacity-check before either is published.
            let text = bounded_json_limit(&Evidence { preview: &preview, diagnostics: &diagnostics },32*1024*1024)?;
            let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
            let revision = verification::plan_revision(&preview.compiled.plan);
            if revision.is_empty() {return Err(bad("compiled plan has no bounded canonical revision"));}
            let plan_bytes = bounded_json(&preview.compiled.plan)?.len().checked_add(bounded_json(&target)?.len()).ok_or_else(||bad("candidate storage accounting overflow"))?;
            let result = json!({"plan_revision":revision,"state":"candidate","verified":false,
                "target":args.target,"inputs":args.inputs,"context_revision":preview.context_revision,
                "preparation_complete":preview.preparation_complete,"diagnostics":diagnostics,"ligation_step":preview.compiled.ligation_step,
                "fragment_count":preview.compiled.rows.len(),"junction_count":preview.junctions.len(),
                "product_physical":preview.product_physical.as_ref().map(|physical|json!({"top_length":physical.top_length,"bottom_length":physical.bottom_length,"topology":physical.topology,"nick_count":physical.nicks.len(),"covalently_closed_in_model":physical.covalently_closed_in_model})),
                "verification":{"status":preview.verification.status,"engine_version":preview.verification.engine_version,
                    "findings":preview.verification.findings.iter().take(100).collect::<Vec<_>>(),
                    "findings_truncated":preview.verification.findings.len()>100,
                    "assumptions":preview.verification.assumptions.iter().take(100).collect::<Vec<_>>(),
                    "assumptions_truncated":preview.verification.assumptions.len()>100,
                    "exhaustive_alternative_search":preview.verification.exhaustive_alternative_search},
                "artifact_id":artifact_id,
                "authority":"Fresh deterministic strand/end/junction reconstruction under explicit preparation and terminal-treatment assumptions. Saving remains a candidate; sequence agreement is not laboratory success or exhaustive mixture uniqueness. Full physical components, junctions and nick evidence are retained in the artifact."});
            let response = output(result)?;
            let mut store=server.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&revision);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(existing)=store.artifacts.get(&artifact_id) {
                if !matches!(existing.data,ArtifactData::Text(_)) || existing.mime_type!="application/json" || existing.filename!="ligation-assembly-recipe.json" {
                    return Err(bad("assembly evidence already has another artifact representation"));
                }
            }
            if (new_plan && store.plans.len()>=MAX_ITEMS) || (new_artifact && store.artifacts.len()>=MAX_ITEMS) {
                return Err(bad("candidate/artifact item limit; no assembly outputs created"));
            }
            let bytes=(if new_plan {plan_bytes}else{0}).checked_add(if new_artifact{text.len()}else{0}).ok_or_else(||bad("session storage accounting overflow"))?;
            store.reserve(bytes)?;
            if new_plan {store.plans.insert(revision,StoredPlan{plan:preview.compiled.plan,target});}
            if new_artifact {store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"ligation-assembly-recipe.json"});}
            Ok(response)
        }).await
    }
    #[tool(
        description = "Compile explicit source-bound TA/GC cloning against an independently supplied target. Each row retains actual prepared ends or explicitly adds one A/T/G/C nucleotide at both 3′ ends, with separate optional end treatment and intended ligation order/orientation. The selected TA/GC profile is assessed separately from generic ligation verification: blunt ligation success does not establish a matching TA profile. Saved origins are freshly reconstructed, retaining original primers, modifications and incomplete evidence. No implicit polymerase, kit, phosphate, host repair, probabilities or experimental isolation is assumed. Retains an immutable candidate and full evidence; verify_plan and HTML/PDF export use its plan_revision."
    )]
    async fn propose_tailing_assembly(
        &self,
        Parameters(args): Parameters<TailingAssemblyProposal>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 128 {
            return Err(bad("tailing recipe requires 1..128 explicit source inputs"));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target) = {
                let store = server.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for bound in &args.inputs {
                    if bound.id.trim().is_empty() || bound.id.len()>128 || !seen.insert(bound.id.as_str()) {
                        return Err(bad("empty, oversized or duplicate recipe input ID"));
                    }
                    sources.push((bound.id.as_str(),store.document(&bound.document)?));
                }
                collagene_core::tailing_recipe::source_preflight(&sources,target,&args.recipe).map_err(bad)?;
                let inputs = sources.into_iter().map(|(id,molecule)|PlanInput {
                    id:id.into(),revision:molecule.revision(),molecule:molecule.clone()
                }).collect::<Vec<_>>();
                (inputs,target.clone())
            };
            let preview = collagene_core::tailing_recipe::preview(&inputs,&target,&args.recipe).map_err(bad)?;
            let expected = collagene_core::tailing_recipe::preview_revision(&inputs,&target,&args.recipe).map_err(bad)?;
            preview.validate_for(&inputs,&target,&args.recipe,&expected).map_err(bad)?;
            // Serialize bounded complete evidence before publishing either output.
            let text = bounded_json_limit(&preview,32*1024*1024)?;
            let artifact_id = blake3::hash(text.as_bytes()).to_hex().to_string();
            let revision = verification::plan_revision(&preview.compiled.plan);
            if revision.is_empty() {return Err(bad("compiled plan has no bounded canonical revision"));}
            let plan_bytes = bounded_json(&preview.compiled.plan)?.len().checked_add(bounded_json(&target)?.len()).ok_or_else(||bad("candidate storage accounting overflow"))?;
            let result = json!({"plan_revision":revision,"state":"candidate","verified":false,
                "target":args.target,"inputs":args.inputs,"context_revision":preview.context_revision,
                "profile_status":preview.profile_status,"profile_findings":preview.profile_findings,
                "preparation_complete":preview.preparation_complete,"ligation_step":preview.compiled.ligation_step,
                "fragment_count":preview.compiled.rows.len(),"junction_count":preview.junctions.len(),
                "product_physical":preview.product_physical.as_ref().map(|physical|json!({"top_length":physical.top_length,"bottom_length":physical.bottom_length,"topology":physical.topology,"nick_count":physical.nicks.len(),"covalently_closed_in_model":physical.covalently_closed_in_model})),
                "verification":{"status":preview.verification.status,"engine_version":preview.verification.engine_version,
                    "findings":preview.verification.findings.iter().take(100).collect::<Vec<_>>(),
                    "findings_truncated":preview.verification.findings.len()>100,
                    "assumptions":preview.verification.assumptions.iter().take(100).collect::<Vec<_>>(),
                    "assumptions_truncated":preview.verification.assumptions.len()>100,
                    "exhaustive_alternative_search":preview.verification.exhaustive_alternative_search},
                "artifact_id":artifact_id,
                "authority":"Fresh deterministic strand/end/junction reconstruction under explicit preparation and terminal-treatment assumptions. Saving remains a candidate; profile matching does not replace ligation verification, and sequence agreement is not laboratory success or exhaustive mixture uniqueness. Full physical components, junctions and nick evidence are retained in the artifact."});
            let response = output(result)?;
            let mut store=server.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&revision);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(existing)=store.artifacts.get(&artifact_id) {
                if !matches!(existing.data,ArtifactData::Text(_)) || existing.mime_type!="application/json" || existing.filename!="tailing-assembly-recipe.json" {
                    return Err(bad("assembly evidence already has another artifact representation"));
                }
            }
            if (new_plan && store.plans.len()>=MAX_ITEMS) || (new_artifact && store.artifacts.len()>=MAX_ITEMS) {
                return Err(bad("candidate/artifact item limit; no assembly outputs created"));
            }
            let bytes=(if new_plan {plan_bytes}else{0}).checked_add(if new_artifact{text.len()}else{0}).ok_or_else(||bad("session storage accounting overflow"))?;
            store.reserve(bytes)?;
            if new_plan {store.plans.insert(revision,StoredPlan{plan:preview.compiled.plan,target});}
            if new_artifact {store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"tailing-assembly-recipe.json"});}
            Ok(response)
        }).await
    }
    #[tool(
        description = "Inspect finite Golden Gate alternatives from freshly reconstructed explicit preparations and an independent target. Reuses actual private PCR/digest/treatment strand state; no serialized chemistry summary supplies a substrate. Redigests the prepared pool under the stated complete-digest, one-copy, orientation and search limits. Includes subsets/self-background and recutting candidates, preserves incomplete search and original ordered findings, and withholds enumeration when required physical preparation is missing. Returns all compact candidate summaries and an immutable full report, never an automatic winner or exhaustive real-mixture/yield claim."
    )]
    async fn inspect_prepared_golden_gate_pool(
        &self,
        Parameters(args): Parameters<PreparedGoldenGatePoolQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 original pool source bindings"));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target) = {
                let store = server.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for input in &args.inputs {
                    if input.id.trim().is_empty() || input.id.len()>128 || !seen.insert(input.id.as_str()) { return Err(bad("invalid or duplicate pool source ID")); }
                    sources.push((input.id.as_str(),store.document(&input.document)?));
                }
                collagene_core::prepared_golden_gate::preflight(&sources,target,&args.request).map_err(bad)?;
                (sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>(),target.clone())
            };
            server.attest_material_origins(&inputs,&args.origin_workspaces,&collagene_core::prepared_golden_gate::material_origins(&args.request))?;
            let report = collagene_core::prepared_golden_gate::inspect(&inputs,&target,&args.request).map_err(bad)?;
            let expected = collagene_core::prepared_golden_gate::inspect_revision(&inputs,&target,&args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            let candidates = report.assessments.iter().map(|a| {
                let c=report.pool.as_ref().and_then(|pool|pool.candidates.get(a.candidate_index)).ok_or_else(||bad("pool assessment has no bound candidate"))?;
                Ok(json!({"candidate_index":a.candidate_index,"directed_sequence_key":a.directed_sequence_key,"physical_sequence_key":a.physical_sequence_key,
                    "target_matches":a.target_matches,"source_roles":a.source_roles,"length":c.duplex.top.len(),"topology":c.duplex.topology,
                    "restriction_site_free":c.restriction_site_free,"recut_sites":c.recut_sites.len(),"recut_complete":c.cleavage.complete,"recut_products_determined":c.cleavage.products_determined}))
            }).collect::<Result<Vec<_>,ErrorData>>()?;
            #[derive(Serialize)]
            struct Evidence<'a> {kind:&'static str,inputs:&'a [BoundInput],target:&'a DocumentRef,origin_workspaces:&'a[MaterialOriginWorkspace],report:&'a collagene_core::prepared_golden_gate::Report}
            let text = bounded_json_limit(&Evidence {kind:"prepared-golden-gate-pool-v1",inputs:&args.inputs,target:&args.target,origin_workspaces:&args.origin_workspaces,report:&report},32*1024*1024)?;
            let artifact = server.text_artifact(text,"application/json","prepared-golden-gate-pool.json")?;
            output(json!({"state":"advisory","verified":false,"context_revision":report.context_revision,
                "plan_revision":report.preview.verification.plan_revision,"target_revision":report.preview.verification.target_revision,
                "inputs":args.inputs,"target":args.target,"preparation_complete":report.preview.preparation_complete,"missing_physical":report.missing_physical,"products_determined":report.pool.as_ref().map(|p|p.products_determined),"profile_revision":report.pool.as_ref().map(|p|p.profile_revision.clone()),"cleavage_comparisons":report.pool.as_ref().map(|p|p.comparisons),"cleavage_budget_scope":"Preparative checked cuts plus Golden Gate source/candidate scans in this native plan; saved-ancestry replay has separate bounded scopes.",
                "ordered_verification":{"status":report.preview.verification.status,"engine_version":report.preview.verification.engine_version,
                    "findings":report.preview.verification.findings.iter().take(100).collect::<Vec<_>>(),"findings_truncated":report.preview.verification.findings.len()>100,
                    "assumptions":report.preview.verification.assumptions.iter().take(100).collect::<Vec<_>>(),"assumptions_truncated":report.preview.verification.assumptions.len()>100},
                "enumeration_complete":report.pool.as_ref().is_some_and(|p|p.enumeration_complete),"exhaustive_mixture_search":false,
                "states_examined":report.pool.as_ref().map_or(0,|p|p.states_examined),"truncation_reasons":report.pool.as_ref().map(|p|&p.truncation_reasons),
                "stable_physical_product_count":report.pool.as_ref().map_or(0,|p|p.stable_physical_product_count),"candidate_count":candidates.len(),"display_truncated":false,
                "roles":report.preview.compiled.rows.iter().map(|r|json!({"index":r.preparation.row_index,"role":r.preparation.role,"input_id":r.preparation.input_id,"reverse_in_ordered_recipe":r.preparation.reverse})).collect::<Vec<_>>(),
                "candidates":candidates,"warnings":report.warnings,"artifact_id":artifact,
                "authority":"Fresh preparation and bounded one-copy pool inspection only; no automatic candidate, laboratory outcome or exhaustive mixture claim. Select an exact index and both keys, then freshly verify its plan."}))
        }).await
    }
    #[tool(
        description = "Select one exact site-free candidate from an immutable prepared Golden Gate pool artifact by context, index and both sequence keys. Preserves original inputs and every preparation/treatment step, replacing only final ordered ligation with the existing GoldenGate operation. Atomically saves an unverified plan and binary cg2plan artifact. Does not trust serialized physical state, install product DNA or assert target agreement. verify_plan must freshly reconstruct the selected plan against the same independent target."
    )]
    async fn select_prepared_golden_gate_candidate(
        &self,
        Parameters(args): Parameters<PreparedGoldenGateSelectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let (evidence,target) = {
                let store = server.store.lock().unwrap();
                let artifact=store.artifacts.get(&args.artifact_id).ok_or_else(||bad("unknown prepared pool artifact"))?;
                if artifact.filename!="prepared-golden-gate-pool.json" || artifact.mime_type!="application/json" { return Err(bad("artifact is not prepared pool evidence")); }
                let ArtifactData::Text(text)=&artifact.data else {return Err(bad("prepared pool evidence must be JSON"));};
                if text.len()>32*1024*1024 {return Err(bad("prepared pool evidence exceeds32MiB"));}
                let evidence: PreparedGoldenGateEvidence=serde_json::from_str(text).map_err(bad)?;
                if evidence.kind!="prepared-golden-gate-pool-v1" || evidence.inputs.len()!=evidence.report.preview.compiled.plan.inputs.len() { return Err(bad("prepared pool source context changed")); }
                evidence.report.validate_binding(&args.selection.context_revision).map_err(bad)?;
                let target=store.document(&evidence.target)?;
                let mut sources=Vec::with_capacity(evidence.inputs.len());
                for (bound,input) in evidence.inputs.iter().zip(&evidence.report.preview.compiled.plan.inputs) {
                    let source=store.document(&bound.document)?;
                    if input.id!=bound.id || input.revision!=bound.document.revision || source.revision()!=input.revision { return Err(bad("prepared pool original source differs from retained evidence")); }
                    sources.push((bound.id.as_str(),source));
                }
                collagene_core::prepared_golden_gate::preflight(&sources,target,&evidence.report.request).map_err(bad)?;
                if target.revision()!=evidence.report.preview.verification.target_revision { return Err(bad("prepared pool independent target changed")); }
                let target=target.clone();
                (evidence,target)
            };
            server.attest_material_origins(&evidence.report.preview.compiled.plan.inputs,&evidence.origin_workspaces,&collagene_core::prepared_golden_gate::material_origins(&evidence.report.request))?;
            let selected=collagene_core::prepared_golden_gate::select(&evidence.report,&args.selection).map_err(bad)?;
            let expected=collagene_core::prepared_golden_gate::selection_revision(&evidence.report,&args.selection).map_err(bad)?;
            selected.validate_binding(&expected).map_err(bad)?;
            let revision=verification::plan_revision(&selected.plan);
            let bytes=transfer::save_candidate(&transfer::Candidate::new(selected.plan.clone())).map_err(bad)?;
            let artifact_id=blake3::hash(&bytes).to_hex().to_string();
            let plan_bytes=bounded_json(&selected.plan)?.len().checked_add(bounded_json(&target)?.len()).ok_or_else(||bad("selected candidate storage overflow"))?;
            let response=output(json!({"state":"candidate","verified":false,"plan_revision":revision,"context_revision":selected.context_revision,"parent_context":selected.parent_context,
                "target":evidence.target,"inputs":evidence.inputs,"selection":selected.selection,"warnings":selected.warnings,
                "candidate_artifact_id":artifact_id,"artifact_id":artifact_id,"filename":"candidate.cg2plan","encoding":"base64","mime_type":"application/octet-stream","bytes":bytes.len()}))?;
            let mut store=server.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&revision);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(old)=store.artifacts.get(&artifact_id) {
                if !matches!(old.data,ArtifactData::Binary(_)) || old.filename!="candidate.cg2plan" || old.mime_type!="application/octet-stream" { return Err(bad("candidate bytes already have another artifact representation")); }
            }
            if new_plan&&store.plans.len()>=MAX_ITEMS || new_artifact&&store.artifacts.len()>=MAX_ITEMS {return Err(bad("candidate/artifact item limit; no pool selection created"));}
            let total=(if new_plan {plan_bytes}else{0}).checked_add(if new_artifact {bytes.len()}else{0}).ok_or_else(||bad("selected candidate storage overflow"))?;
            store.reserve(total)?;
            if new_plan {store.plans.insert(revision,StoredPlan{plan:selected.plan,target});}
            if new_artifact {store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Binary(bytes),mime_type:"application/octet-stream",encoding:"base64",filename:"candidate.cg2plan"});}
            Ok(response)
        }).await
    }
    #[tool(
        description = "Inspect an explicitly declared Golden Gate reaction from whole original circular/linear DNA, complete PCR amplicons or selected preparations. Checked positional cleavage includes all established fragments of every declared substrate and reports proven intact circular carryover separately. Unknown or impaired source cleavage withholds the entire requested pool; candidate recut evidence retains uncertainty. Actual preparation strands and enumeration come from one native pass; metadata cannot supply physical state. Retains recuttable products, backgrounds and search limits. Index0 probe diagnostics do not nominate a candidate or establish kinetics, repeated-copy stoichiometry, methylation clearance or experimental yield."
    )]
    async fn inspect_golden_gate_reaction(
        &self,
        Parameters(args): Parameters<GoldenGateReactionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.inputs.is_empty() || args.inputs.len() > 32 {
            return Err(bad("choose1..32 original pool source bindings"));
        }
        let server = self.clone();
        self.compute(move || {
            let (inputs,target) = {
                let store = server.store.lock().unwrap();
                let target = store.document(&args.target)?;
                let mut seen = BTreeSet::new();
                let mut sources = Vec::with_capacity(args.inputs.len());
                for input in &args.inputs {
                    if input.id.trim().is_empty() || input.id.len()>128 || !seen.insert(input.id.as_str()) { return Err(bad("invalid or duplicate pool source ID")); }
                    sources.push((input.id.as_str(),store.document(&input.document)?));
                }
                collagene_core::golden_gate_reaction::preflight(&sources,target,&args.request).map_err(bad)?;
                (sources.into_iter().map(|(id,m)|PlanInput{id:id.into(),revision:m.revision(),molecule:m.clone()}).collect::<Vec<_>>(),target.clone())
            };
            server.attest_material_origins(&inputs,&args.origin_workspaces,&collagene_core::golden_gate_reaction::material_origins(&args.request))?;
            let report = collagene_core::golden_gate_reaction::inspect(&inputs,&target,&args.request).map_err(bad)?;
            let expected = collagene_core::golden_gate_reaction::inspect_revision(&inputs,&target,&args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            let candidates = report.assessments.iter().map(|a| {
                let c=report.pool.as_ref().and_then(|pool|pool.candidates.get(a.candidate_index)).ok_or_else(||bad("pool assessment has no bound candidate"))?;
                Ok(json!({"candidate_index":a.candidate_index,"directed_sequence_key":a.directed_sequence_key,"physical_sequence_key":a.physical_sequence_key,
                    "target_matches":a.target_matches,"source_roles":a.source_roles,"length":c.duplex.top.len(),"topology":c.duplex.topology,
                    "restriction_site_free":c.restriction_site_free,"recut_sites":c.recut_sites.len(),"recut_complete":c.cleavage.complete,"recut_products_determined":c.cleavage.products_determined}))
            }).collect::<Result<Vec<_>,ErrorData>>()?;
            #[derive(Serialize)]
            struct Evidence<'a> {kind:&'static str,inputs:&'a [BoundInput],target:&'a DocumentRef,origin_workspaces:&'a[MaterialOriginWorkspace],report:&'a collagene_core::golden_gate_reaction::Report}
            let text = bounded_json_limit(&Evidence {kind:"golden-gate-reaction-v1",inputs:&args.inputs,target:&args.target,origin_workspaces:&args.origin_workspaces,report:&report},32*1024*1024)?;
            let artifact = server.text_artifact(text,"application/json","golden-gate-reaction.json")?;
            output(json!({"state":"advisory","verified":false,"context_revision":report.context_revision,
                "probe_plan_revision":report.verification.plan_revision,"target_revision":report.verification.target_revision,
                "inputs":args.inputs,"target":args.target,"preparation_complete":(report.verification.products.len()+1>=report.probe_plan.steps.len()),"pool_available":report.pool.is_some(),"products_determined":report.pool.as_ref().map(|p|p.products_determined),"profile_revision":report.pool.as_ref().map(|p|p.profile_revision.clone()),"cleavage_comparisons":report.pool.as_ref().map(|p|p.comparisons),"cleavage_budget_scope":"Preparative checked cuts plus Golden Gate source/candidate scans in this native plan; saved-ancestry replay has separate bounded scopes.",
                "probe_verification":{"status":report.verification.status,"engine_version":report.verification.engine_version,
                    "findings":report.verification.findings.iter().take(100).collect::<Vec<_>>(),"findings_truncated":report.verification.findings.len()>100,
                    "assumptions":report.verification.assumptions.iter().take(100).collect::<Vec<_>>(),"assumptions_truncated":report.verification.assumptions.len()>100},
                "enumeration_complete":report.pool.as_ref().is_some_and(|p|p.enumeration_complete),"exhaustive_mixture_search":false,
                "states_examined":report.pool.as_ref().map_or(0,|p|p.states_examined),"truncation_reasons":report.pool.as_ref().map(|p|&p.truncation_reasons),
                "stable_physical_product_count":report.pool.as_ref().map_or(0,|p|p.stable_physical_product_count),"candidate_count":candidates.len(),"display_truncated":false,
                "roles":report.materials.iter().map(|r|json!({"index":r.index,"role":r.role,"input_id":r.input_id,"source":r.source})).collect::<Vec<_>>(),
                "carryover_sources":report.pool.as_ref().map(|p|&p.carryover_sources),
                "candidates":candidates,"warnings":report.warnings,"artifact_id":artifact,
                "authority":"Fresh declared-reaction inspection only; the probe selects no saved winner and establishes no laboratory outcome or exhaustive mixture claim. Select an exact index and both keys, then freshly verify its plan."}))
        }).await
    }
    #[tool(
        description = "Select one exact site-free candidate from an immutable whole-reaction Golden Gate artifact by context, index and both keys. Preserves original inputs, all explicit preparations and the declared reaction substrates; changes only the final GoldenGate candidate index. Atomically saves an unverified plan and binary cg2plan artifact. Does not trust serialized physical state, install product DNA or assert target agreement. verify_plan must freshly reconstruct the selected plan against the same independent target."
    )]
    async fn select_golden_gate_reaction_candidate(
        &self,
        Parameters(args): Parameters<PreparedGoldenGateSelectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let (evidence,target) = {
                let store = server.store.lock().unwrap();
                let artifact=store.artifacts.get(&args.artifact_id).ok_or_else(||bad("unknown prepared pool artifact"))?;
                if artifact.filename!="golden-gate-reaction.json" || artifact.mime_type!="application/json" { return Err(bad("artifact is not prepared pool evidence")); }
                let ArtifactData::Text(text)=&artifact.data else {return Err(bad("prepared pool evidence must be JSON"));};
                if text.len()>32*1024*1024 {return Err(bad("prepared pool evidence exceeds32MiB"));}
                let evidence: GoldenGateReactionEvidence=serde_json::from_str(text).map_err(bad)?;
                if evidence.kind!="golden-gate-reaction-v1" || evidence.inputs.len()!=evidence.report.probe_plan.inputs.len() { return Err(bad("prepared pool source context changed")); }
                evidence.report.validate_binding(&args.selection.context_revision).map_err(bad)?;
                let target=store.document(&evidence.target)?;
                let mut sources=Vec::with_capacity(evidence.inputs.len());
                for (bound,input) in evidence.inputs.iter().zip(&evidence.report.probe_plan.inputs) {
                    let source=store.document(&bound.document)?;
                    if input.id!=bound.id || input.revision!=bound.document.revision || source.revision()!=input.revision { return Err(bad("prepared pool original source differs from retained evidence")); }
                    sources.push((bound.id.as_str(),source));
                }
                collagene_core::golden_gate_reaction::preflight(&sources,target,&evidence.report.request).map_err(bad)?;
                if target.revision()!=evidence.report.verification.target_revision { return Err(bad("prepared pool independent target changed")); }
                let target=target.clone();
                (evidence,target)
            };
            server.attest_material_origins(&evidence.report.probe_plan.inputs,&evidence.origin_workspaces,&collagene_core::golden_gate_reaction::material_origins(&evidence.report.request))?;
            let selected=collagene_core::golden_gate_reaction::select(&evidence.report,&args.selection).map_err(bad)?;
            let expected=collagene_core::golden_gate_reaction::selection_revision(&evidence.report,&args.selection).map_err(bad)?;
            selected.validate_binding(&expected).map_err(bad)?;
            let revision=verification::plan_revision(&selected.plan);
            let bytes=transfer::save_candidate(&transfer::Candidate::new(selected.plan.clone())).map_err(bad)?;
            let artifact_id=blake3::hash(&bytes).to_hex().to_string();
            let plan_bytes=bounded_json(&selected.plan)?.len().checked_add(bounded_json(&target)?.len()).ok_or_else(||bad("selected candidate storage overflow"))?;
            let response=output(json!({"state":"candidate","verified":false,"plan_revision":revision,"context_revision":selected.context_revision,"parent_context":selected.parent_context,
                "target":evidence.target,"inputs":evidence.inputs,"selection":selected.selection,"warnings":selected.warnings,
                "candidate_artifact_id":artifact_id,"artifact_id":artifact_id,"filename":"candidate.cg2plan","encoding":"base64","mime_type":"application/octet-stream","bytes":bytes.len()}))?;
            let mut store=server.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&revision);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if let Some(old)=store.artifacts.get(&artifact_id) {
                if !matches!(old.data,ArtifactData::Binary(_)) || old.filename!="candidate.cg2plan" || old.mime_type!="application/octet-stream" { return Err(bad("candidate bytes already have another artifact representation")); }
            }
            if new_plan&&store.plans.len()>=MAX_ITEMS || new_artifact&&store.artifacts.len()>=MAX_ITEMS {return Err(bad("candidate/artifact item limit; no pool selection created"));}
            let total=(if new_plan {plan_bytes}else{0}).checked_add(if new_artifact {bytes.len()}else{0}).ok_or_else(||bad("selected candidate storage overflow"))?;
            store.reserve(total)?;
            if new_plan {store.plans.insert(revision,StoredPlan{plan:selected.plan,target});}
            if new_artifact {store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Binary(bytes),mime_type:"application/octet-stream",encoding:"base64",filename:"candidate.cg2plan"});}
            Ok(response)
        }).await
    }
    #[tool(
        description = "Save an immutable candidate plan binding an independent target and input revisions. DeclareSourceModifications accepts only an explicit Input source with Dam/Dcm/CpG generation assumptions; CheckedDigest pins an exact CatalogSelection and explicit inspection limits, consumes freshly reconstructed tracked strands and reports the fixed current positional profile. DeclareSourceTermini requires an original intact linear DNA Input and explicit left/right5-prime/3-prime assumptions; it cannot promote Step, ssDNA, circular, imported physical or reference-only states, infer chemistry from metadata or repair a graph. OverlapExtension requires two freshly reconstructed intact linear duplexes, explicit orientations and terminal pairing/search options; it reciprocally extends A.top/B.bottom and can feed actual outer-primer PCR. A sequence-only Overlap/ProfiledOverlap cannot supply that physical template. Saving never verifies or applies a plan. Duplicate submissions are idempotent."
    )]
    fn propose_plan(
        &self,
        Parameters(args): Parameters<Proposal>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.name.trim().is_empty()
            || args.name.len() > 1024
            || args.steps.is_empty()
            || args.steps.len() > 100
            || args.inputs.len() > 128
            || serde_json::to_vec(&args).map_err(bad)?.len() > 2 * 1024 * 1024
        {
            return Err(bad(
                "plan requires a name, 1–100 steps, at most 128 inputs and at most 2 MiB",
            ));
        }
        let mut store = self.store.lock().unwrap();
        let mut estimated = serde_json::to_vec(&args).map_err(bad)?.len();
        estimated = estimated.saturating_add(
            serde_json::to_vec(store.document(&args.target)?)
                .map_err(bad)?
                .len(),
        );
        for input in &args.inputs {
            estimated = estimated.saturating_add(
                serde_json::to_vec(store.document(&input.document)?)
                    .map_err(bad)?
                    .len(),
            );
            if estimated > MAX_STORE_BYTES {
                return Err(bad(
                    "candidate input materialization exceeds session storage budget",
                ));
            }
        }
        let target = store.document(&args.target)?.clone();
        let mut seen = std::collections::BTreeSet::new();
        let mut inputs = Vec::new();
        for input in args.inputs {
            if input.id.trim().is_empty() || input.id.len() > 128 || !seen.insert(input.id.clone())
            {
                return Err(bad("empty or duplicate input ID"));
            }
            inputs.push(PlanInput {
                id: input.id,
                revision: input.document.revision.clone(),
                molecule: store.document(&input.document)?.clone(),
            });
        }
        let plan = DesignPlan {
            name: args.name,
            target_revision: args.target.revision.clone(),
            inputs,
            steps: args.steps,
        };
        let revision = verification::plan_revision(&plan);
        if !store.plans.contains_key(&revision) {
            if store.plans.len() >= MAX_ITEMS {
                return Err(bad("at most 128 candidate plans"));
            }
            let bytes = serde_json::to_vec(&plan).map_err(bad)?.len()
                + serde_json::to_vec(&target).map_err(bad)?.len();
            store.reserve(bytes)?;
            store
                .plans
                .insert(revision.clone(), StoredPlan { plan, target });
        }
        output(
            json!({"plan_revision":revision,"state":"candidate","verified":false,"target":args.target}),
        )
    }
    #[tool(
        description = "Freshly inspect DNA junctions of an exact saved overlap, ordered ligation or selected Golden Gate candidate and its stored independent target. Golden Gate boundaries and source mappings come from its actual selected digest path, captured in the same verifier pass; recipe order and target DNA never determine them. Explicit bounded flanks, genetic code, six-frame choices and CDS inspection options. Reports actual product boundary/homology/orientation/circular-origin windows, boundary-crossing codons, conditional source-CDS evidence and physical ends/nicks where supported. Full target verification and unsupported final-method diagnostics remain visible. Immutable JSON evidence only: no source edits, candidate promotion, inferred fusion function or experimental guarantee."
    )]
    async fn inspect_assembly_junctions(
        &self,
        Parameters(args): Parameters<AssemblyJunctionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            // Borrow actual saved source/target snapshots. No source or plan copy
            // precedes the core's bounded native-model preflight.
            let (plan,target,expected)={let store=server.store.lock().unwrap();
                let saved=store.plans.get(&args.plan_revision).ok_or_else(||bad("unknown plan revision"))?;
                collagene_core::assembly_junctions::source_preflight(&saved.plan,&saved.target).map_err(bad)?;
                let expected=collagene_core::assembly_junctions::inspection_revision(&saved.plan,&saved.target,&args.options).map_err(bad)?;
                (saved.plan.clone(),saved.target.clone(),expected)};
            let report=collagene_core::assembly_junctions::inspect_plan(&plan,&target,&args.options).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&report,32*1024*1024)?,"application/json","assembly-junctions.json")?;
            let junctions=report.junctions.iter().take(10).map(|junction|json!({
                "index":junction.index,"upstream":junction.upstream,"downstream":junction.downstream,
                "closing":junction.closing,"boundary":junction.boundary,"homology":junction.homology,
                "window":junction.window,
                "frames":junction.frames.iter().take(6).map(|frame|json!({"spec":frame.spec,
                    "crossing_codon":frame.crossing_codon,
                    "protein_preview":frame.translation.protein.chars().take(100).collect::<String>(),
                    "protein_truncated":frame.translation.protein.len()>100})).collect::<Vec<_>>(),
                "cds_count":junction.cds.len(),"cds_omitted":junction.cds_omitted,
                "cds":junction.cds.iter().take(10).map(|cds|json!({"feature_index":cds.feature_index,
                    "label":cds.label.chars().take(256).collect::<String>(),"label_truncated":cds.label.chars().count()>256,
                    "needs_review":cds.needs_review,"spans_boundary":cds.spans_boundary,
                    "warning":cds.warning.as_ref().map(|warning|warning.chars().take(1024).collect::<String>()),
                    "warning_truncated":cds.warning.as_ref().is_some_and(|warning|warning.chars().count()>1024)})).collect::<Vec<_>>(),
                "cds_preview_truncated":junction.cds.len()>10,
                "physical":junction.physical.as_ref().map(|physical|json!({"upstream_end":physical.upstream_end,"downstream_end":physical.downstream_end,"nick_count":physical.product_nicks.len()}))
            })).collect::<Vec<_>>();
            output(json!({"plan_revision":report.plan_revision,"target_revision":report.target_revision,
                "context_revision":report.context_revision,"engine_version":report.engine_version,"options":report.options,
                "verification_status":report.verification.status,"verification_findings":report.verification.findings.iter().take(100).map(|finding|finding.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "verification_findings_truncated":report.verification.findings.len()>100 || report.verification.findings.iter().take(100).any(|f|f.chars().count()>1024),
                "product_revision":report.product_revision,"product_length":report.product_length,"product_topology":report.product_topology,
                "method":report.method,"complete":report.complete,"omitted_junctions":report.omitted_junctions,
                "junction_count":report.junctions.len(),"junctions":junctions,"junction_preview_truncated":report.junctions.len()>10,
                "warnings":report.warnings.iter().take(100).map(|w|w.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),
                "warnings_truncated":report.warnings.len()>100 || report.warnings.iter().take(100).any(|w|w.chars().count()>1024),
                "artifact_id":artifact,
                "authority":"Fresh declared-model junction DNA and conditional annotation/frame evidence only. Local codons do not prove a functional fusion or alter full target verification; saved evidence never promotes this candidate or repairs physical strands."}))
        }).await
    }
    #[tool(
        description = "Reconstruct an exact saved plan revision against its separately bound target with the scientific core. Returns findings, assumptions and bounded product summaries. OverlapExtension freshly reconstructs the selected reciprocal extension; competing or unsearched terminal pairings remain explicit findings. Never trusts AI claimed products or status."
    )]
    async fn verify_plan(
        &self,
        Parameters(args): Parameters<PlanRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.verify_plan_inner(args)).await
    }
    fn verify_plan_inner(&self, args: PlanRef) -> Result<CallToolResult, ErrorData> {
        let saved = self
            .store
            .lock()
            .unwrap()
            .plans
            .get(&args.plan_revision)
            .cloned()
            .ok_or_else(|| bad("unknown plan revision"))?;
        let report = verification::verify(&saved.plan, &saved.target);
        let scientific_step_count = report.physical_processing_reports.len()
            + report.reference_step_reports.len()
            + report
                .step_reports
                .iter()
                .filter(|step| step.topo.is_some())
                .count();
        let reference_step_count = report.reference_step_reports.len();
        let artifact = if scientific_step_count > 0
            || reference_step_count > 0
            || !report.step_reports.is_empty()
        {
            Some(self.text_artifact(
                bounded_json_limit(&report, 32 * 1024 * 1024)?,
                "application/json",
                "verified-plan-scientific-evidence.json",
            )?)
        } else {
            None
        };
        let products:Vec<_>=report.products.iter().enumerate().map(|(i,m)|json!({"step":i,"revision":m.revision(),"length":m.sequence.len(),"topology":m.topology,"duplex":m.metadata.get("duplex_model").and_then(|s|serde_json::from_str::<Value>(s).ok()),"tracked_modifications":m.metadata.get("tracked_modifications").and_then(|s|serde_json::from_str::<Value>(s).ok())})).collect();
        output(
            json!({"plan_revision":report.plan_revision,"target_revision":report.target_revision,"engine_version":report.engine_version,"status":report.status,"findings":report.findings,"assumptions":report.assumptions,"products":products,"exhaustive_alternative_search":report.exhaustive_alternative_search,"scientific_step_count":scientific_step_count,"reference_step_count":reference_step_count,"artifact_id":artifact}),
        )
    }
    #[tool(
        description = "Generate an HTML protocol artifact from an exact saved plan revision. Always re-verifies at export, visibly includes invalid/incomplete findings and assumptions. Read the artifact in bounded chunks; no files are written."
    )]
    async fn export_protocol(
        &self,
        Parameters(args): Parameters<PlanRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.export_protocol_inner(args))
            .await
    }
    #[tool(
        description = "Export actual declared oligos from an exact saved plan and its stored independent target as spreadsheet-safe CSV or TSV plus complete use/domain/tail/mismatch evidence JSON. Exact ordered DNA/modifications may deduplicate; every original use and name remains separate. No PCR/target verification is run, no primer is designed, and downstream phosphorylation never becomes a supplier order modification. Missing/unrepresentable original order data is explicit complete:false/unavailable; do not order from a partial sheet. No quantities, stocks or supplier codes are invented."
    )]
    async fn export_plan_oligos(
        &self,
        Parameters(args): Parameters<PlanOligoExport>,
    ) -> Result<CallToolResult, ErrorData> {
        let delimiter = match args.format {
            PlanOligoFormat::Csv => oligo_formats::Delimiter::Comma,
            PlanOligoFormat::Tsv => oligo_formats::Delimiter::Tab,
        };
        let filename = match args.format {
            PlanOligoFormat::Csv => "plan-oligos.csv",
            PlanOligoFormat::Tsv => "plan-oligos.tsv",
        };
        let server = self.clone();
        self.compute(move || {
            let exported={let store=server.store.lock().unwrap();let saved=store.plans.get(&args.plan_revision).ok_or_else(||bad("unknown plan revision"))?;
                let expected=collagene_core::plan_oligos::export_revision(&saved.plan,&saved.target,&args.options,delimiter).map_err(bad)?;
                let exported=collagene_core::plan_oligos::export(&saved.plan,&saved.target,&args.options,delimiter).map_err(bad)?;
                exported.validate_binding(&expected).map_err(bad)?;exported};
            if exported.report.plan_revision!=args.plan_revision{return Err(bad("stored plan revision changed"));}
            let result=json!({"plan_revision":exported.report.plan_revision,"target_revision":exported.report.target_revision,
                "context_revision":exported.context_revision,"options":exported.report.options,"verified":false,
                "complete":exported.report.complete,"unavailable":exported.report.unavailable,
                "oligo_count":exported.report.oligos.records.len(),"use_count":exported.report.uses.len(),
                "uses":exported.report.uses.iter().take(100).collect::<Vec<_>>(),"uses_truncated":exported.report.uses.len()>100,
                "warnings":exported.report.warnings,"export_warnings":exported.export.warnings,
                "authority":"Declared order/use evidence only. PCR products, independent target agreement and experimental success require separate verification. Incomplete sheets are review evidence, not all-required orders; declared starting termini are not supplier modification requests."});
            let report_text=bounded_json_limit(&exported.report,32*1024*1024)?;
            let (sheet,report)=server.paired_text_artifacts([(exported.export.text,"text/plain",filename),(report_text,"application/json","plan-oligos.json")])?;
            let mut result=result;result["table_artifact_id"]=json!(sheet);result["report_artifact_id"]=json!(report);result["filename"]=json!(filename);output(result)
        }).await
    }
    #[tool(
        description = "Generate an embedded vector PDF protocol from an exact saved plan revision, re-verifying all steps and the independent target at export. Invalid/incomplete drafts remain conspicuous. Returns a bounded base64 artifact; read chunks, concatenate and decode to a PDF. No files, browser or external renderer are invoked."
    )]
    async fn export_protocol_pdf(
        &self,
        Parameters(args): Parameters<PlanRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let saved=server.store.lock().unwrap().plans.get(&args.plan_revision).cloned().ok_or_else(||bad("unknown plan revision"))?;
            let report=collagene_core::pdf::export_protocol_pdf(&saved.plan,&saved.target).map_err(bad)?;
            let id=server.binary_artifact(report.bytes.clone(),"application/pdf","cloning-protocol.pdf")?;
            output(json!({"artifact_id":id,"filename":"cloning-protocol.pdf","mime_type":"application/pdf","encoding":"base64","bytes":report.bytes.len(),"pages":report.pages,"status":report.status,"plan_revision":report.plan_revision,"target_revision":report.target_revision,"engine_version":report.engine_version,"warnings":report.warnings}))
        }).await
    }
    fn export_protocol_inner(&self, args: PlanRef) -> Result<CallToolResult, ErrorData> {
        let saved = self
            .store
            .lock()
            .unwrap()
            .plans
            .get(&args.plan_revision)
            .cloned()
            .ok_or_else(|| bad("unknown plan revision"))?;
        let html = verification::export_protocol_html(&saved.plan, &saved.target);
        let id = blake3::hash(html.as_bytes()).to_hex().to_string();
        let length = html.len();
        let mut store = self.store.lock().unwrap();
        if !store.artifacts.contains_key(&id) {
            if store.artifacts.len() >= MAX_ITEMS {
                return Err(bad("at most 128 artifacts"));
            }
            store.reserve(length)?;
            store.artifacts.insert(
                id.clone(),
                Artifact {
                    data: ArtifactData::Text(html),
                    mime_type: "text/html; charset=utf-8",
                    encoding: "utf-8",
                    filename: "cloning-protocol.html",
                },
            );
        }
        output(
            json!({"artifact_id":id,"mime_type":"text/html; charset=utf-8","bytes":length,"plan_revision":args.plan_revision}),
        )
    }
    #[tool(
        description = "Read an artifact in at most 16,384 text bytes. Follow next_start to done; UTF-8 text is literal, base64 text must be concatenated and decoded before saving binary files."
    )]
    fn read_artifact(
        &self,
        Parameters(args): Parameters<ArtifactRegion>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.store.lock().unwrap();
        let artifact = store
            .artifacts
            .get(&args.artifact_id)
            .ok_or_else(|| bad("unknown artifact"))?;
        if args.max_bytes == 0 || args.max_bytes > 16_384 {
            return Err(bad("artifact chunk size must be 1..16384"));
        }
        let (chunk, end, total) = match &artifact.data {
            ArtifactData::Text(text) => {
                if args.start > text.len() || !text.is_char_boundary(args.start) {
                    return Err(bad("invalid UTF-8 artifact offset"));
                }
                let mut end = args.start.saturating_add(args.max_bytes).min(text.len());
                while !text.is_char_boundary(end) {
                    end -= 1;
                }
                if end == args.start && end < text.len() {
                    return Err(bad("chunk too small for next UTF-8 character"));
                }
                (text[args.start..end].to_owned(), end, text.len())
            }
            ArtifactData::Binary(bytes) => {
                let total = bytes.len().div_ceil(3) * 4;
                if args.start > total {
                    return Err(bad("invalid base64 artifact offset"));
                }
                let end = args.start.saturating_add(args.max_bytes).min(total);
                let group_start = args.start / 4;
                let raw_start = (group_start * 3).min(bytes.len());
                let raw_end = (end.div_ceil(4) * 3).min(bytes.len());
                let encoded = STANDARD.encode(&bytes[raw_start..raw_end]);
                let skip = args.start - group_start * 4;
                (
                    encoded[skip..skip + (end - args.start)].to_owned(),
                    end,
                    total,
                )
            }
        };
        output(
            json!({"artifact_id":args.artifact_id,"start":args.start,"next_start":end,"done":end==total,"total_text_bytes":total,"text":chunk,"encoding":artifact.encoding,"mime_type":artifact.mime_type,"filename":artifact.filename}),
        )
    }
}

mod genome_store;

mod unit_tools;

mod assembly_library_history_tools;
mod catalog_tools;
mod combinatorial_tools;
mod domestication_tools;
mod fusion_design_tools;
mod multisite_recombination_tools;
mod published_ligation_tools;
mod recombination_tools;
mod topo_actions_tools;

#[cfg(test)]
mod finalized_workspace_atomic_capacity {
    use super::*;
    #[test]
    fn both_finalization_artifacts_are_reserved_before_any_workspace_publication() {
        let server = CollageneServer::default();
        {
            let mut store = server.store.lock().unwrap();
            for index in 0..MAX_ITEMS - 1 {
                store.artifacts.insert(
                    format!("reserved-{index}"),
                    Artifact {
                        data: ArtifactData::Text(String::new()),
                        mime_type: "text/plain",
                        encoding: "utf-8",
                        filename: "reserved.txt",
                    },
                );
            }
        }
        let source = Molecule::dna("source", "ACGT", collagene_core::Topology::Linear).unwrap();
        let workspace =
            workspace::Workspace::from_documents("Atomic fixture", vec![source], 0).unwrap();
        let native = workspace::save_workspace(&workspace).unwrap();
        let native_id = blake3::hash(&native).to_hex().to_string();
        let text = "{\"fixture\":\"receipt\"}".to_string();
        let report_id = blake3::hash(text.as_bytes()).to_hex().to_string();
        assert!(
            server
                .commit_workspace_outputs(workspace, Some((text, json!({}))), Some(native))
                .is_err()
        );
        let store = server.store.lock().unwrap();
        assert_eq!(store.artifacts.len(), MAX_ITEMS - 1);
        assert!(store.workspaces.is_empty());
        assert!(store.documents.is_empty());
        assert_eq!(store.used, 0);
        assert!(!store.artifacts.contains_key(&native_id));
        assert!(!store.artifacts.contains_key(&report_id));
    }
}

/// Serve the embedded server in this process. stdout belongs exclusively to MCP;
/// this helper neither opens a GUI nor starts an external executable.
pub async fn serve_stdio(
    server: CollageneServer,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    use rmcp::ServiceExt;
    server
        .serve(rmcp::transport::stdio())
        .await?
        .waiting()
        .await?;
    Ok(())
}

mod sequence_exchange_tools;

mod quantity_tools;
mod topo_graph_tools;
