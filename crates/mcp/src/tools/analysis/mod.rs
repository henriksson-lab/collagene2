mod base_editing;
mod codon;
mod design_checks;
mod ivt;

use crate::tools::default_inspection_page;

use crate::server::CollageneServer;
use crate::server::EmptyRequest;
use crate::server::response::bad;
use crate::server::response::bounded_json;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::session::ArtifactData;
use crate::session::MAX_ITEMS;
use crate::session::MAX_STORE_BYTES;
use crate::session::Store;
use crate::session::references::DocumentRef;
use crate::session::references::EvidenceRef;
use std::sync::Arc;

use base64::engine::general_purpose::STANDARD;
use collagene_core::alignment;
use collagene_core::formats::alignment as alignment_formats;
use collagene_model::Molecule;
use rmcp::ErrorData;
use rmcp::handler::server::wrapper::Parameters;
use rmcp::model::CallToolResult;
use schemars::JsonSchema;
use serde::Deserialize;
use serde::Serialize;
use serde_json::Value;

use base64::Engine as _;
use collagene_core::formats::oligos as oligo_formats;
use collagene_core::review::reads::trace_validation;
use collagene_core::sequence::translation;
use collagene_core::verification;
use rmcp::tool;
use rmcp::tool_router;
use serde_json::json;

pub(crate) const MAX_READ_ALIGNMENT_SOURCE: usize = 4 * 1024 * 1024;
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AlignmentImport {
    /// Exact file bytes as base64, at most8MiB decoded. No paths, URLs or aligner invocation.
    pub data_base64: String,
    pub format: alignment_formats::AlignmentFormat,
    pub alphabet: collagene_model::Alphabet,
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
    pub(crate) fn core(self) -> collagene_core::review::variants::reference::exchange::Format {
        match self {
            Self::Vcf => collagene_core::review::variants::reference::exchange::Format::Vcf,
            Self::Bcf => collagene_core::review::variants::reference::exchange::Format::Bcf,
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
    pub(crate) fn core(self) -> collagene_core::review::reads::store::Format {
        match self {
            Self::Sam => collagene_core::review::reads::store::Format::Sam,
            Self::Bam => collagene_core::review::reads::store::Format::Bam,
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
pub(crate) fn read_source_cap() -> usize {
    MAX_READ_ALIGNMENT_SOURCE
}
pub(crate) fn read_retained_cap() -> usize {
    32 * 1024 * 1024
}
pub(crate) fn read_record_cap() -> usize {
    10_000
}
pub(crate) fn read_preview() -> usize {
    20
}
pub(crate) fn read_export_cap() -> usize {
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
pub(crate) fn default_read_check_columns() -> usize {
    100_000
}
pub(crate) fn default_read_check_work() -> u64 {
    128_000_000
}
pub(crate) fn default_read_check_report() -> usize {
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
    pub options: Option<collagene_core::formats::alignment::profiles::ExportOptions>,
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
pub(crate) fn default_crispr_page() -> usize {
    20
}
fn crispr_hit_preview(hit: &collagene_core::sequence::crispr::OffTargetHit) -> serde_json::Value {
    json!({
        "id":hit.id,"reference_index":hit.reference_index,"reference_guide":hit.reference_guide,
        "pam_sequence":hit.pam_sequence,"guide_positions":hit.guide_positions,"pam_positions":hit.pam_positions,
        "reverse":hit.reverse,"min_mismatches":hit.min_mismatches,"max_mismatches":hit.max_mismatches,
        "pam_certain":hit.pam_certain,"design_pam":hit.design_pam,"definite":hit.definite,
        "is_intended":hit.is_intended,"nominal_top_cut":hit.nominal_top_cut,"nominal_bottom_cut":hit.nominal_bottom_cut,
        "alignment":{
            "complete":hit.alignment.complete,"limit_reached":hit.alignment.limit_reached,"work":hit.alignment.work,
            "alignments_retained":hit.alignment.alignments.len(),"preview_truncated":hit.alignment.alignments.len()>8,
            "alignments":hit.alignment.alignments.iter().take(8).map(|path|json!({
                "min_mismatches":path.min_mismatches,"max_mismatches":path.max_mismatches,
                "dna_bulge_bases":path.dna_bulge_bases,"rna_bulge_bases":path.rna_bulge_bases,"definite":path.definite,
            })).collect::<Vec<_>>(),
            "column_evidence":"Complete retained columns are in the immutable full report artifact."
        }
    })
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CrisprGuidesQuery {
    pub document: DocumentRef,
    #[serde(default)]
    pub options: collagene_core::sequence::crispr::GuideOptions,
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
    pub query: collagene_core::sequence::crispr::GuideQuery,
    pub references: Vec<CrisprReferenceRef>,
    #[serde(default)]
    pub options: collagene_core::sequence::crispr::OffTargetOptions,
    #[serde(default)]
    pub hit_start: usize,
    #[serde(default = "default_crispr_page")]
    pub max_hits: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct MethylationQuery {
    pub document: DocumentRef,
    pub options: collagene_core::chemistry::methylation::MethylationOptions,
    #[serde(default)]
    pub site_start: usize,
    #[serde(default = "default_inspection_page")]
    pub max_sites: usize,
}
pub(crate) fn default_contig_page() -> usize {
    5
}
pub(crate) fn default_overlap_page() -> usize {
    10
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AssembleContigsQuery {
    pub reads: Vec<ContigReadRef>,
    pub options: collagene_core::review::contigs::ContigOptions,
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
    pub options: collagene_core::protein::properties::ProteinOptions,
}
pub(crate) fn default_dotplot_page() -> usize {
    50
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DotPlotQuery {
    pub reference: DocumentRef,
    pub query: DocumentRef,
    pub options: collagene_core::alignment::dotplot::DotPlotOptions,
    #[serde(default)]
    pub hit_start: usize,
    #[serde(default = "default_dotplot_page")]
    pub max_hits: usize,
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
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SimilarityQuery {
    pub query: DocumentRef,
    pub subjects: Vec<DocumentRef>,
    #[serde(default)]
    pub options: collagene_core::sequence::similarity::Options,
}

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AlignmentAntibodyQuery {
    pub alignment: EvidenceRef,
    pub row_id: String,
    pub options: collagene_core::antibody_numbering::domain_candidates::Options,
    pub limits: collagene_core::antibody_numbering::domain_candidates::Limits,
    pub selection: collagene_core::antibody_numbering::domain_tracks::Selection,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AntibodyQuery {
    pub document: DocumentRef,
    pub options: collagene_core::antibody_numbering::domain_candidates::Options,
    pub limits: collagene_core::antibody_numbering::domain_candidates::Limits,
    pub selection: collagene_core::antibody_numbering::domain_tracks::Selection,
}

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProteinLiabilityQuery {
    pub document: DocumentRef,
    pub options: collagene_core::protein_liabilities::Options,
}

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RnaSourceQuery {
    pub document: DocumentRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RnaChemistryQuery {
    pub document: DocumentRef,
    pub declaration: collagene_core::rna_chemistry::Declaration,
}

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct StructureImport {
    /// Complete ASCII PDB/mmCIF file, at most 4 MiB. Original line endings remain evidence.
    pub text: String,
    pub format: collagene_core::protein::structure::coordinates::Format,
    /// Null chooses the lowest model ID, all chains/alternates and deposited H/D.
    #[serde(default)]
    pub selection: Option<collagene_core::protein::structure::coordinates::scene::Selection>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct StructureInspect {
    /// Exact source handle returned by import_structure; not a molecule or report handle.
    pub source: EvidenceRef,
    #[serde(default)]
    pub selection: Option<collagene_core::protein::structure::coordinates::scene::Selection>,
}
pub(crate) fn capture_variant(
    store: &Store,
    args: &VariantCheck,
) -> Result<
    (
        Arc<collagene_core::review::variants::reference::exchange::Collection>,
        Molecule,
        collagene_core::review::variants::reference::Options,
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
    collagene_model::bounds::model_preflight(&[reference]).map_err(bad)?;
    if collection.record(args.record_index).is_none() {
        return Err(bad("variant record index unavailable"));
    }
    let options = collagene_core::review::variants::reference::Options {
        source_blake3: collection.source_blake3().into(),
        record_index: args.record_index,
        reference_name: args.reference_name.clone(),
        expected_reference_revision: args.reference.revision.clone(),
        alternate_index: args.alternate_index,
    };
    Ok((collection, reference.clone(), options))
}
pub(crate) fn variant_summary(
    collection: &collagene_core::review::variants::reference::exchange::Collection,
) -> Value {
    json!({"collection":{"id":collection.source_blake3(),"revision":collection.source_blake3()},
        "source_blake3":collection.source_blake3(),"source_bytes":collection.original_bytes().len(),
        "format":match collection.format(){collagene_core::review::variants::reference::exchange::Format::Vcf=>"vcf",collagene_core::review::variants::reference::exchange::Format::Bcf=>"bcf"},
        "record_count":collection.len(),"sample_column_count":collection.header().sample_names().len(),
        "sample_columns":collection.header().sample_names().iter().take(32).map(|s|read_name(Some(s))).collect::<Vec<_>>(),
        "sample_columns_truncated":collection.header().sample_names().len()>32,
        "authority":"Declared file evidence only; header/build/sample/genotype/INFO assertions are unverified. No implicit DNA reference mapping, variant calling or clinical inference."})
}
pub(crate) fn read_alignment_summary(
    collection: &collagene_core::review::reads::store::Collection,
) -> Value {
    let summary = collection.summary();
    json!({"collection":{"id":summary.source_blake3,"revision":summary.source_blake3},
        "source_blake3":summary.source_blake3,"source_bytes":summary.source_bytes,"format":summary.format,
        "record_count":summary.records,"reference_count":summary.references.len(),
        "flag_unmapped":summary.flag_unmapped,"flag_secondary":summary.flag_secondary,"flag_supplementary":summary.flag_supplementary,
        "retained_charge":collection.retained_charge(),"file_complete":true,
        "warnings":summary.warnings.iter().take(20).map(|s|s.chars().take(512).collect::<String>()).collect::<Vec<_>>(),
        "warnings_omitted":summary.warnings.len().saturating_sub(20),"interpretation":summary.interpretation})
}
pub(crate) fn read_name(value: Option<&str>) -> Option<String> {
    value.map(|s| s.chars().take(256).collect())
}
#[tool_router(router = analysis_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Design source-bound CRISPR guide/PAM candidates under an explicit canonical SpCas9, Cas12a, SaCas9, SpCas9-VRER or custom nuclease geometry and bounded GC/homopolymer/search filters. Default options use SpCas9; returns available builtin profiles, paged exact positions/strands, nominal cuts, ambiguity/filter evidence and complete immutable JSON. Model geometry and nominal cuts do not establish activity, editing outcome or genome-wide specificity. No sequence edits or persistent application writes."
    )]
    pub(crate) async fn design_crispr_guides(
        &self,
        Parameters(args): Parameters<CrisprGuidesQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_candidates == 0 || args.max_candidates > 100 {
            return Err(bad("guide pages require max_candidates 1..100"));
        }
        let server = self.clone();
        self.compute(move||{
            let molecule={let store=server.session.store.lock().unwrap();let molecule=store.document(&args.document)?;collagene_core::sequence::crispr::guide_preflight(molecule,&args.options).map_err(bad)?;molecule.clone()};
            let report=collagene_core::sequence::crispr::design_guides(&molecule,&args.options).map_err(bad)?;
            report.validate_bindings_for(&molecule,&args.options).map_err(bad)?;
            if args.candidate_start>report.candidates.len(){return Err(bad("guide page start exceeds retained candidate count"));}
            let end=args.candidate_start.saturating_add(args.max_candidates).min(report.candidates.len());
            let artifact_id=server.session.text_artifact(bounded_json(&report)?,"application/json","crispr-guides.json")?;
            output(json!({"document":args.document,"engine_version":report.engine_version,"source_revision":report.source_revision,
                "batch_revision":report.batch_revision,"options_revision":report.options_revision,"profile_revision":report.profile_revision,
                "options":report.options,"source_length":report.source_length,"circular":report.circular,
                "complete":report.complete,"comparisons":report.comparisons,"limits":report.limits,"warnings":report.warnings,
                "candidate_count":report.candidates.len(),"candidate_start":args.candidate_start,"next_candidate_start":end,"done":end==report.candidates.len(),
                "candidates":&report.candidates[args.candidate_start..end],"artifact_id":artifact_id,
                "available_profiles":collagene_core::sequence::crispr::catalogue(),
                "authority":"Source-bound PAM/guide geometry under the declared model; no empirical activity, editing outcome or whole-genome specificity."}))
        }).await
    }
    #[tool(
        description = "Search an explicit immutable selected reference set for bounded full-input CRISPR guide/target alignment evidence with exact PAM/strand/coordinate binding. Borrowed reference preflight precedes owned copies. Build/accession/scope labels are user declarations, not independently verified genome coverage. Every retained hit remains visible; is_intended requires exact role, revision, guide/PAM positions and orientation. Returns bounded hit pages and full JSON evidence. Explicit DNA/RNA bulge caps include terminal/mixed gaps; zero caps use the same ungapped policy. A separate declared screening PAM retains its rationale/sources and design-PAM evidence. Per-locus pages contain bounded path summaries; full retained columns remain in the artifact. No cleavage, genome-wide safety, sequence editing or persistent writes."
    )]
    pub(crate) async fn search_crispr_off_targets(
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
            let references={let store=server.session.store.lock().unwrap();
                let borrowed=args.references.iter().map(|r|Ok(collagene_core::sequence::crispr::ReferenceBorrow{id:&r.id,molecule:store.document(&r.document)?,build:r.build.as_deref(),accession:r.accession.as_deref()})).collect::<Result<Vec<_>,ErrorData>>()?;
                collagene_core::sequence::crispr::reference_batch_preflight(&args.query,&borrowed,&args.options).map_err(bad)?;
                borrowed.into_iter().map(|r|collagene_core::sequence::crispr::ReferenceInput{id:r.id.to_owned(),molecule:r.molecule.clone(),build:r.build.map(str::to_owned),accession:r.accession.map(str::to_owned)}).collect::<Vec<_>>()};
            let report=collagene_core::sequence::crispr::search_off_targets(&args.query,&references,&args.options).map_err(bad)?;
            report.validate_bindings_for(&args.query,&references,&args.options).map_err(bad)?;
            if args.hit_start>report.hits.len(){return Err(bad("off-target page start exceeds retained hit count"));}
            let end=args.hit_start.saturating_add(args.max_hits).min(report.hits.len());
            let artifact_id=server.session.text_artifact(bounded_json(&report)?,"application/json","crispr-off-targets.json")?;
            output(json!({"reference_documents":args.references,"engine_version":report.engine_version,"alignment_model":report.alignment_model,"batch_revision":report.batch_revision,
                "options_revision":report.options_revision,"profile_revision":report.profile_revision,"query":report.query,"options":report.options,
                "references":report.references,"complete":report.complete,"exhaustive_for_selected_records":report.exhaustive_for_selected_records,
                "comparisons":report.comparisons,"index_entries":report.index_entries,"index_used":report.index_used,"limits":report.limits,"warnings":report.warnings,
                "hit_count":report.hits.len(),"hit_start":args.hit_start,"next_hit_start":end,"done":end==report.hits.len(),"hits":report.hits[args.hit_start..end].iter().map(crispr_hit_preview).collect::<Vec<_>>(),"artifact_id":artifact_id,
                "authority":"Declared global alignment and screening-PAM evidence for the explicit immutable selected reference set only; declared build labels do not prove complete genome coverage, activity or empirical safety."}))
        }).await
    }
    #[tool(
        description = "Assess starter-enzyme recognition sites under explicitly declared Dam, Dcm and CpG modification assumptions using versioned supplier profiles and bounded sequence-context matching. Returns source/options/profile-bound site evidence, explicit ambiguity/unsupported/partial-search warnings, paged site summaries and a complete immutable JSON artifact. Hemi/mixed and unknown states cannot silently become unmodified. This assessment does not modify DNA, predict digest products, prove empirical activity or grant cloning verification authority."
    )]
    pub(crate) async fn assess_methylation(
        &self,
        Parameters(args): Parameters<MethylationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_sites == 0 || args.max_sites > 100 {
            return Err(bad("methylation site pages require max_sites 1..100"));
        }
        let server = self.clone();
        self.compute(move || {
            let molecule = {
                let store = server.session.store.lock().unwrap();
                let molecule = store.document(&args.document)?;
                collagene_core::chemistry::methylation::source_preflight(molecule, &args.options).map_err(bad)?;
                molecule.clone()
            };
            let report = collagene_core::chemistry::methylation::assess(&molecule, &args.options).map_err(bad)?;
            report.validate_bindings_for(&molecule, &args.options).map_err(bad)?;
            if args.site_start > report.sites.len() {
                return Err(bad("methylation page start exceeds retained site count"));
            }
            let end = args.site_start.saturating_add(args.max_sites).min(report.sites.len());
            let text = bounded_json(&report)?;
            let artifact_id = server.session.text_artifact(text, "application/json", "methylation-assessment.json")?;
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
    pub(crate) async fn assemble_contigs(
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
                let store = server.session.store.lock().unwrap();
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
                args.reads.iter().map(|r| Ok(collagene_core::review::contigs::ReadInput {
                    id: r.id.clone(), molecule: store.document(&r.document)?.clone(),
                })).collect::<Result<Vec<_>, ErrorData>>()?
            };
            let batch_revision = collagene_core::review::contigs::batch_revision(&reads, &args.options).map_err(bad)?;
            let cached = {
                let store = server.session.store.lock().unwrap();
                store.contig_reports.get(&batch_revision)
                    .and_then(|id| store.artifacts.get(id))
                    .and_then(|a| match &a.data { ArtifactData::Text(text) if a.mime_type=="application/json" && a.filename=="contig-assembly.json" => Some(text.clone()), _=>None })
            };
            let reused_result = cached.is_some();
            let report = match cached {
                Some(text) => serde_json::from_str::<collagene_core::review::contigs::ContigReport>(&text).map_err(bad)?,
                None => collagene_core::review::contigs::assemble(&reads, &args.options).map_err(bad)?,
            };
            report.validate_bindings_for(&reads, &args.options).map_err(bad)?;
            if args.contig_start > report.contigs.len() || args.overlap_start > report.overlaps.len() {
                return Err(bad("assembly page start exceeds the retained result count"));
            }
            let end = args.contig_start.saturating_add(args.max_contigs).min(report.contigs.len());
            let overlap_end = args.overlap_start.saturating_add(args.max_overlaps).min(report.overlaps.len());
            let molecules = report.contigs.iter().filter_map(|c| c.molecule.clone()).collect();
            let (artifact_id, documents) = server.session.contig_artifacts(&report, molecules, &batch_revision)?;
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
    pub(crate) async fn analyze_dotplot(
        &self,
        Parameters(args): Parameters<DotPlotQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_hits == 0 || args.max_hits > 200 {
            return Err(bad("dot plot hit pages require max_hits1..200"));
        }
        let server = self.clone();
        self.compute(move || {
            let (reference, query) = { let store=server.session.store.lock().unwrap(); (store.document(&args.reference)?.clone(), store.document(&args.query)?.clone()) };
            let report = collagene_core::alignment::dotplot::analyze(&reference,&query,&args.options).map_err(bad)?;
            report.validate_bindings_for(&reference,&query,&args.options).map_err(bad)?;
            if args.hit_start > report.hits.len() { return Err(bad("dot plot hit_start exceeds retained hit count")); }
            let end = args.hit_start.saturating_add(args.max_hits).min(report.hits.len());
            let svg=collagene_core::alignment::dotplot::export_svg(&report).map_err(bad)?;
            let (artifact_id,svg_artifact_id)=server.session.dotplot_artifacts(&report,svg)?;
            output(json!({"reference":args.reference,"query":args.query,"engine_version":report.engine_version,"reference_revision":report.reference_revision,"query_revision":report.query_revision,"options_revision":report.options_revision,"reference_length":report.reference_length,"query_length":report.query_length,"alphabet":report.alphabet,"reference_circular":report.reference_circular,"query_circular":report.query_circular,"strategy":report.strategy,"work_used":report.work_used,"examined_windows":report.examined_windows,"total_windows":report.total_windows,"grid_complete":report.grid_complete,"complete":report.complete,"sampled":report.sampled,"hit_count":report.hits.len(),"hit_start":args.hit_start,"next_hit_start":end,"done":end==report.hits.len(),"hits":&report.hits[args.hit_start..end],"limits":report.limits,"warnings":report.warnings,"artifact_id":artifact_id,"svg_artifact_id":svg_artifact_id,"authority":"Window-matching evidence only; not an alignment, homology proof or physical ancestry."}))
        }).await
    }
    #[tool(
        description = "Analyze an immutable protein document with explicit terminal-stop, ambiguity, disulfide and modification assumptions. Embedded composition, average/monoisotopic mass, charge/pI and extinction models expose unavailable or bounded values rather than inventing exact ambiguity/modification chemistry. Full revision-bound report is a JSON artifact; no sequence edits or experimental claims."
    )]
    pub(crate) async fn analyze_protein(
        &self,
        Parameters(args): Parameters<ProteinPropertiesQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let molecule = server
                .session
                .store
                .lock()
                .unwrap()
                .document(&args.document)?
                .clone();
            let report =
                collagene_core::protein::properties::analyze_protein(&molecule, &args.options)
                    .map_err(bad)?;
            let value = serde_json::to_value(&report).map_err(bad)?;
            let artifact = server.session.analysis_artifact(&value)?;
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
    pub(crate) async fn import_variants(
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
            let collection =
                collagene_core::review::variants::reference::exchange::Collection::import(
                    &bytes,
                    args.format.core(),
                    &collagene_core::review::variants::session::mcp_limits(),
                )
                .map_err(bad)?;
            let id = collection.source_blake3().to_owned();
            let response = variant_summary(&collection);
            let mut store = server.session.store.lock().unwrap();
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
    pub(crate) async fn inspect_variants(
        &self,
        Parameters(args): Parameters<VariantPage>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.limit) {
            return Err(bad("variant pages must be 1..100"));
        }
        let server = self.clone();
        self.compute(move||{
            let collection=server.session.store.lock().unwrap().variants(&args.collection)?;
            let rows=collection.page(args.start,args.limit).map_err(bad)?;
            let full=rows.iter().map(|r|json!({"index":r.index,"reference_name":r.reference_name,"position_1based":r.position_1based,"vcf":r.vcf})).collect::<Vec<_>>();
            let evidence=json!({"collection":args.collection,"source_blake3":collection.source_blake3(),"start":args.start,"records":full,"interpretation":"Normalized declared variant rows; exact original lexical bytes remain separately available."});
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","variant-record-page.json")?;
            let mut response=variant_summary(&collection);response["start"]=json!(args.start);response["next_start"]=json!(args.start+rows.len());response["done"]=json!(args.start+rows.len()==collection.len());
            response["records"]=json!(rows.iter().map(|r|json!({"index":r.index,"reference_name":read_name(Some(&r.reference_name)),"reference_name_truncated":r.reference_name.chars().count()>256,"position_1based":r.position_1based,"vcf_preview":r.vcf.chars().take(1024).collect::<String>(),"vcf_preview_truncated":r.vcf.chars().count()>1024})).collect::<Vec<_>>());response["artifact_id"]=json!(artifact);output(response)
        }).await
    }
    #[tool(
        description = "Freshly compare one explicitly selected 1-based ALT and declared REF against actual immutable DNA DocumentRef with exact file contig mapping. Requires private collection source hash and current document revision; wrong contig/stale identity rejected. Returns bounded source/record/reference-bound JSON evidence. Match only establishes REF correspondence, never genotype, build identity, sample status or variant validity; no name-based automatic reference lookup."
    )]
    pub(crate) async fn check_variant_reference(
        &self,
        Parameters(args): Parameters<VariantCheck>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (collection,reference,options)=capture_variant(&server.session.store.lock().unwrap(),&args)?;
            let report=collagene_core::review::variants::reference::check_record(&collection,&reference,&options).map_err(bad)?;
            let evidence=json!({"collection":args.collection,"reference":args.reference,"report":report});
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","variant-reference-check.json")?;
            output(json!({"collection":args.collection,"reference":args.reference,"record_index":args.record_index,"alternate_index":args.alternate_index,"revision":report.revision,"outcome":report.outcome,"interpretation":report.interpretation,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Freshly check a literal selected variant REF/ALT against an explicitly mapped immutable DNA reference, then create an edited reference-only preview. Symbolic/SV/breakend/ambiguous or mismatching alleles cannot produce a preview. Atomically installs reference-only document plus bounded source/reference-bound JSON artifact; never physical material, confirmed genotype or target-guided variant selection. Original sources remain unchanged."
    )]
    pub(crate) async fn preview_variant_allele(
        &self,
        Parameters(args): Parameters<VariantCheck>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (collection,reference,options)=capture_variant(&server.session.store.lock().unwrap(),&args)?;
            let preview=collagene_core::review::variants::reference::preview(&collection,&reference,&options).map_err(bad)?;
            let evidence=json!({"collection":args.collection,"reference":args.reference,"report":preview.report,"interpretation":"Reference-only explicit selected-allele design; no measured genotype or physical-material authority."});
            let (artifact,document)=server.session.candidate_document_artifact_bounded(&evidence,Some(preview.molecule),"variant-allele-preview.json",32*1024*1024)?;
            output(json!({"collection":args.collection,"reference":args.reference,"record_index":args.record_index,"alternate_index":args.alternate_index,"revision":preview.report.revision,"outcome":preview.report.outcome,"reference_only":true,"document":document,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Replay a retained variant reference-only preview from its bounded unverified recipe against the original private collection source_id and actual immutable reference/design DocumentRefs. Freshly reconstructs and compares literal DNA plus topology only. Returns full bounded source/reference/design-bound JSON evidence and claimed/fresh revisions; installs nothing and restores no genotype, physical chemistry, biological function or historical verification authority. Designs without retained recipes require explicit original record/ALT selection instead."
    )]
    pub(crate) async fn replay_variant_preview(
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
                let store=server.session.store.lock().unwrap();
                let collection=store.variants(&EvidenceRef{id:args.source_id.clone(),revision:args.source_id.clone()})?;
                let reference=store.document(&args.reference)?;let design=store.document(&args.design)?;
                collagene_model::bounds::model_preflight(&[reference,design]).map_err(bad)?;
                (collection,reference.clone(),design.clone())
            };
            let replay=collagene_core::review::variants::reference::replay_preview(&collection,&reference,&design).map_err(bad)?;
            let fresh_revision=replay.fresh.molecule.revision();
            let evidence=json!({"source_id":args.source_id,"reference":args.reference,"design":args.design,"fresh_revision":fresh_revision,"replay":replay});
            let artifact=server.session.text_artifact(bounded_json_limit(&evidence,32*1024*1024)?,"application/json","variant-preview-replay.json")?;
            output(json!({"source_id":args.source_id,"reference":args.reference,"design":args.design,
                "literal_dna_and_topology_match":replay.literal_dna_and_topology_match,"claimed_revision":replay.claimed_revision,"fresh_revision":fresh_revision,"interpretation":replay.interpretation,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Export a private completed variant collection as exact original bytes or normalized VCF/BCF, with explicit conversion notices/errors for unsupported representations. Returns bounded binary artifact read via read_artifact; any codec/size/storage failure publishes no partial artifact. No paths or external bcftools. Original is the byte-exact evidence mode."
    )]
    pub(crate) async fn export_variants(
        &self,
        Parameters(args): Parameters<VariantExport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_bytes == 0 || args.max_bytes > MAX_STORE_BYTES {
            return Err(bad("variant export byte cap must be 1..64 MiB"));
        }
        let server = self.clone();
        self.compute(move||{
            use collagene_core::review::variants::reference::exchange::{Format, ExportFormat};
            let collection=server.session.store.lock().unwrap().variants(&args.collection)?;
            let mode=match args.format{VariantExportFormat::Original=>ExportFormat::Original,VariantExportFormat::Vcf=>ExportFormat::Vcf,VariantExportFormat::Bcf=>ExportFormat::Bcf};
            let format=match args.format{VariantExportFormat::Original=>collection.format(),VariantExportFormat::Vcf=>Format::Vcf,VariantExportFormat::Bcf=>Format::Bcf};
            let exported=collection.export(mode,args.max_bytes).map_err(bad)?;let len=exported.bytes.len();
            let artifact=server.session.binary_artifact(exported.bytes,"application/octet-stream",if format==Format::Vcf{"variants.vcf"}else{"variants.bcf"})?;
            output(json!({"collection":args.collection,"source_blake3":collection.source_blake3(),"bytes":len,"format":if format==Format::Vcf{"vcf"}else{"bcf"},"original":matches!(args.format,VariantExportFormat::Original),"notices":exported.notices,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Dispose one exact private variant collection source-hash handle, releasing its retained session charge. Existing independent output artifacts/reference previews remain unchanged. A serialized collection summary cannot restore parser authority; reimport exact original bytes explicitly."
    )]
    pub(crate) async fn dispose_variants(
        &self,
        Parameters(args): Parameters<VariantQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let mut store = self.session.store.lock().unwrap();
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
    pub(crate) async fn import_read_alignments(
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
            let mut limits = collagene_core::review::reads::store::ImportLimits::default();
            limits.max_source_bytes = args.max_source_bytes;
            limits.max_retained_bytes = args.max_retained_bytes;
            limits.max_records = args.max_records;
            limits.reader.max_input_bytes = args.max_source_bytes as u64;
            limits.reader.max_records = args.max_records as u64;
            let collection = collagene_core::review::reads::store::Collection::import(
                &bytes,
                args.format.core(),
                &limits,
            )
            .map_err(bad)?;
            let id = collection.source_blake3().to_owned();
            let response = read_alignment_summary(&collection);
            let mut store = server.session.store.lock().unwrap();
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
    pub(crate) async fn inspect_read_alignments(
        &self,
        Parameters(args): Parameters<ReadAlignmentPage>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.limit) || !(1..=100).contains(&args.reference_limit) {
            return Err(bad("read alignment record/reference pages must be1..100"));
        }
        let server = self.clone();
        self.compute(move|| {
            let collection=server.session.store.lock().unwrap().read_alignments(&args.collection)?;
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
    pub(crate) async fn check_read_alignment(
        &self,
        Parameters(args): Parameters<ReadAlignmentCheck>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.preview_columns > 100 || args.reference_name.is_empty() {
            return Err(bad(
                "explicit reference name and at most100 preview columns required",
            ));
        }
        let options = collagene_core::review::reads::reference::Options {
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
                let store=server.session.store.lock().unwrap();let collection=store.read_alignments(&args.collection)?;
                let reference=store.document(&args.reference)?;
                collagene_model::bounds::model_preflight(&[reference]).map_err(bad)?;
                let record=collection.record(args.record_index).ok_or_else(||bad("unknown read alignment record index"))?;
                let named=collection.header().reference(&args.reference_name).ok_or_else(||bad("explicit reference name is absent from imported dictionary"))?;
                if record.fields().reference.as_ref().is_some_and(|name|name!=&args.reference_name)
                    || named.length as usize!=reference.sequence.len()
                    || named.circular!=(reference.topology==collagene_model::Topology::Circular) {
                    return Err(bad("read alignment dictionary/reference mapping differs"));
                }
                (collection,reference.clone())
            };
            let checked=collection.check_reference(args.record_index,&reference,&options).map_err(bad)?;
            let artifact=server.session.text_artifact(bounded_json(&checked)?,"application/json","read-alignment-reference-check.json")?;
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
    pub(crate) async fn export_read_alignments(
        &self,
        Parameters(args): Parameters<ReadAlignmentExport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.max_bytes == 0 || args.max_bytes > MAX_STORE_BYTES {
            return Err(bad("read alignment export cap must be1..64MiB"));
        }
        let server = self.clone();
        self.compute(move|| {
            use collagene_core::review::reads::store::{Format, ImportLimits};
            let collection=server.session.store.lock().unwrap().read_alignments(&args.collection)?;
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
            let artifact=server.session.binary_artifact(bytes,"application/octet-stream",if format=="sam"{"read-alignments.sam"}else{"read-alignments.bam"})?;
            output(json!({"collection":args.collection,"source_blake3":collection.source_blake3(),"format":format,
                "original":matches!(args.format,ReadAlignmentExportFormat::Original),"bytes":length,"notices":notices,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Import aligned FASTA, CLUSTAL, Stockholm, strict/relaxed sequential/interleaved PHYLIP, MSF, NEXUS, PIR/NBRF, GDE flat or SELEX as exact base64 bytes with explicit DNA/RNA/protein alphabet. Strict bounded parsing retains original bytes, stable rows, gaps and supported metadata; no alignment engine is invoked. Creates immutable session evidence only."
    )]
    pub(crate) async fn import_alignment(
        &self,
        Parameters(args): Parameters<AlignmentImport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > alignment_formats::MAX_ALIGNMENT_BYTES.div_ceil(3) * 4 {
            return Err(bad("alignment exceeds8MiB decoded"));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            server.session.insert_alignment(
                alignment_formats::import_alignment(&bytes, args.format, args.alphabet, &args.name)
                    .map_err(bad)?,
            )
        })
        .await
    }
    #[tool(
        description = "Import a complete native alignment JSON attachment as base64, at most32MiB decoded. Source hashes, retained exact bytes, stable row identities and column metadata are validated. Does not edit browser/workspace documents."
    )]
    pub(crate) async fn import_alignment_document(
        &self,
        Parameters(args): Parameters<AlignmentDocumentImport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > (32 * 1024 * 1024usize).div_ceil(3) * 4 {
            return Err(bad("alignment attachment exceeds32MiB decoded"));
        }
        let server = self.clone();
        self.compute(move || {
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            server
                .session
                .insert_alignment(alignment_formats::load_alignment(&bytes).map_err(bad)?)
        })
        .await
    }
    #[tool(
        description = "List immutable session alignment evidence, at most128 collections. Summaries omit rows and raw source; read_alignment pages rows/columns and export_alignment retains full evidence."
    )]
    pub(crate) fn list_alignments(
        &self,
        Parameters(_): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.session.store.lock().unwrap();
        output(
            json!({"alignments":store.alignments.iter().map(|(id,a)|json!({"id":id,"revision":id,"name":a.name.chars().take(256).collect::<String>(),"alphabet":a.alphabet,"row_count":a.rows.len(),"column_count":a.column_count(),"source_blake3":a.source.blake3})).collect::<Vec<_>>()}),
        )
    }
    #[tool(
        description = "Read a bounded immutable alignment window: up to50 rows and1000 columns, zero-based offsets. Returns as-stored row residue indices or null gaps, never external reference coordinates. Imported consensus symbols are source evidence, not computed confidence. Exact alignment revision required."
    )]
    pub(crate) async fn read_alignment(
        &self,
        Parameters(args): Parameters<AlignmentRegionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let doc=server.session.store.lock().unwrap().alignment(&args.alignment)?;
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
    pub(crate) async fn export_alignment(
        &self,
        Parameters(args): Parameters<AlignmentExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let doc=server.session.store.lock().unwrap().alignment(&args.alignment)?;
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
                    use collagene_core::formats::alignment::AlignmentFormat as F;
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
                        collagene_core::formats::alignment::profiles::export_alignment_with_options(&doc,format,options).map_err(bad)?
                    }else{alignment_formats::export_alignment(&doc,format).map_err(bad)?};
                    let filename=match format{F::AlignedFasta=>"alignment.fasta",F::Clustal=>"alignment.aln",F::Stockholm=>"alignment.sto",F::Msf=>"alignment.msf",F::Nexus=>"alignment.nex",F::Pir=>"alignment.pir",F::Gde=>"alignment.gde",F::Selex=>"alignment.slx",_=>"alignment.phy"};
                    (report.text.into_bytes(),report.losses,report.renamings,filename)
                }
            };
            let hash=blake3::hash(&bytes).to_hex().to_string();let size=bytes.len();let text=String::from_utf8(bytes).map_err(bad)?;
            let artifact=server.session.text_artifact(text,"text/plain",filename)?;
            output(json!({"alignment":args.alignment,"artifact_id":artifact,"bytes":size,"content_blake3":hash,"losses":losses,"renamings":renamings,"encoding":"utf-8","source_is_current":doc.source_is_current().map_err(bad)?}))
        }).await
    }
    #[tool(
        description = "Explicitly extract selected alignment rows into new immutable session molecules, removing gaps in as-stored orientation. Original annotations, consensus/confidence and external coordinates are not invented. All-gap, unknown or duplicate selections fail atomically. Does not modify browser/workspace documents."
    )]
    pub(crate) async fn extract_alignment_rows(
        &self,
        Parameters(args): Parameters<AlignmentRowsQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.row_ids.is_empty() || args.row_ids.len() > alignment_formats::MAX_ALIGNMENT_ROWS {
            return Err(bad("select1..1000 alignment rows"));
        }
        let server = self.clone();
        self.compute(move||{let doc=server.session.store.lock().unwrap().alignment(&args.alignment)?;let rows=doc.extract_rows(&args.row_ids).map_err(bad)?;server.session.insert_documents(rows,vec!["Extracted ungapped rows retain source/current alignment hashes; no external reference coordinates, alignment confidence or annotations were copied.".into()])}).await
    }
    #[tool(
        description = "Validate revision-bound original SCF/ABIF/explicit-encoding FASTQ reads against an independent DNA reference with embedded trimming, affine placement/uniqueness proof, coverage and variant evidence. Native trace scores are not invented Phred. Full JSON artifact; no physical/sample confirmation or automatic document mutation."
    )]
    pub(crate) async fn validate_reads(
        &self,
        Parameters(args): Parameters<ReadValidationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (reference,reads)={let store=server.session.store.lock().unwrap();(store.document(&args.reference)?.clone(),args.reads.iter().map(|r|store.document(r).cloned()).collect::<Result<Vec<_>,_>>()?)};
            let report=trace_validation::validate_reads(&reference,&reads,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"reference":args.reference,"reads":args.reads,"validation":report}))?;
            output(json!({"reference":args.reference,"read_revisions":report.read_revisions,"model":report.model,"engine_version":report.engine_version,"status":report.status,"eligible_coverage":report.eligible_coverage,"reference_length":report.consensus.len(),"reference_coverage_fraction":report.reference_coverage_fraction,"consensus_preview":report.consensus.chars().take(2000).collect::<String>(),"uncovered_ranges":report.uncovered_ranges.iter().take(100).collect::<Vec<_>>(),"uncovered_ranges_truncated":report.uncovered_ranges.len()>100,"differences":report.sites.iter().filter(|s|s.state!=trace_validation::SiteState::SupportedMatch).take(100).collect::<Vec<_>>(),"differences_truncated":report.sites.iter().filter(|s|s.state!=trace_validation::SiteState::SupportedMatch).count()>100,"reads":report.reads.iter().map(|r|json!({"name":r.name.chars().take(256).collect::<String>(),"read_format":r.read_format,"quality_source":r.trim.as_ref().map(|t|t.quality_source),"trim":r.trim.as_ref().map(|t|json!({"start":t.start,"end":t.end,"passing_calls":t.passing_calls,"low_or_unknown_count":t.low_or_unknown_calls.len()})),"placement_unique_proven":r.placement_unique_proven,"optimality_proven":r.optimality_proven,"search_complete":r.search_complete,"paired_reference_fraction":r.paired_reference_fraction,"error":r.error,"warnings":r.warnings})).collect::<Vec<_>>(),"insertion_observations":report.insertions.iter().take(100).map(|i|json!({"reference_boundary":i.reference_boundary,"sequence_preview":i.sequence.chars().take(100).collect::<String>(),"length":i.sequence.len(),"read_index":i.read_index,"eligible":i.eligible})).collect::<Vec<_>>(),"insertions_truncated":report.insertions.len()>100,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Import FASTQ caller evidence with required explicit Phred33 or Phred64 encoding. No ambiguous auto-detection, quality fabrication or automatic verification. Atomic immutable documents; text <=8 MiB."
    )]
    pub(crate) async fn import_fastq(
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
            server
                .session
                .insert_documents(report.molecules, report.warnings)
        })
        .await
    }
    #[tool(
        description = "Export exact called FASTQ with mapped caller qualities in a required explicit encoding. Stale qualities or unrepresentable scores fail instead of being invented/clamped. Original source attachment is preserved in native documents."
    )]
    pub(crate) async fn export_fastq(
        &self,
        Parameters(args): Parameters<FastqExport>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let report=oligo_formats::export_fastq(&molecule,args.encoding).map_err(bad)?;
            let artifact=server.session.text_artifact(report.text.clone(),"text/plain","reads.fastq")?;
            output(json!({"document":args.document,"encoding":args.encoding,"bytes":report.text.len(),"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Translate a revision-bound DNA/RNA sequence with an explicit NCBI genetic code, frame, orientation and initiation policy. Returns bounded protein preview and a full JSON evidence artifact; ambiguous codons and possible stops remain explicit."
    )]
    pub(crate) async fn translate_sequence(
        &self,
        Parameters(args): Parameters<TranslateQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let report=translation::translate(&molecule.sequence,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"translation":report}))?;
            output(json!({"document":args.document,"options":report.options,"code_name":report.code_name,"amino_acids":report.protein.len(),"protein_preview":report.protein.chars().take(2000).collect::<String>(),"preview_truncated":report.protein.len()>2000,"unresolved_codons":report.unresolved_codons,"possible_stop_codons":report.possible_stop_codons,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Translate a specific CDS annotation using its compound location, codon_start and transl_table. Unsupported exceptions/fuzzy locations fail explicitly. Stored-translation agreement is evidence, not experimental confirmation. Full report available as JSON artifact."
    )]
    pub(crate) async fn translate_cds(
        &self,
        Parameters(args): Parameters<CdsQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let feature=molecule.features.get(args.feature_index).ok_or_else(||bad("CDS annotation index does not exist"))?;
            let report=translation::translate_cds(&molecule,feature).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"cds":report}))?;
            output(json!({"document":args.document,"feature_index":args.feature_index,"location":report.location,"code":report.translation.options.genetic_code,"protein_preview":report.protein.chars().take(2000).collect::<String>(),"amino_acids":report.protein.len(),"preview_truncated":report.protein.len()>2000,"complete_under_model":report.complete_under_model,"stored_translation_matches":report.stored_translation_matches,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Find bounded six-frame ORF candidates with explicit start/partial/code policies. Circular searches cover at most one source turn. Predictions carry exact/partial/possible-stop evidence; they are not validated genes. Returns first100 summaries plus full JSON artifact."
    )]
    pub(crate) async fn find_orfs(
        &self,
        Parameters(args): Parameters<OrfQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let report=translation::find_orfs(&molecule,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"document":args.document,"engine_version":verification::ENGINE_VERSION,"orfs":report}))?;
            output(json!({"document":args.document,"options":report.options,"count":report.orfs.len(),"truncated_search":report.truncated,"display_truncated":report.orfs.len()>100,"orfs":report.orfs.iter().take(100).map(|o|json!({"location":o.location,"reverse":o.reverse,"frame":o.frame,"amino_acids":o.protein.len(),"protein_preview":o.protein.chars().take(80).collect::<String>(),"crosses_origin":o.crosses_origin,"partial_5prime":o.partial_5prime,"partial_3prime":o.partial_3prime,"termination":o.termination})).collect::<Vec<_>>(),"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Analyze one exact protein alignment row against the explicit current species/locus V/J reference model. Requires immutable alignment reference, stable row ID, explicit native options and work/retention limits. Preserves gaps via a bounded ungapped-index to alignment-column map; gaps never acquire numbering. Returns conditional track summaries with explicit numbering/region choices and full report artifact without alignment, document, history or material mutation. Not calibrated domain classification or independent numbering certification."
    )]
    pub(crate) async fn analyze_alignment_antibody(
        &self,
        Parameters(args): Parameters<AlignmentAntibodyQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            use collagene_core::protein::antibody_alignment as a;
            let doc={server.session.store.lock().unwrap().alignment(&args.alignment)?};
            let source=a::source(&doc,&args.row_id).map_err(bad)?;
            let r=a::analyze(&doc,&source,&args.options,&args.limits,&args.selection,&mut || false).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&serde_json::to_value(&r).map_err(bad)?)?;
            output(json!({"alignment":args.alignment,"source":r.source,"references":args.options.references,"models":r.analysis.evidence.v_evidence.mappings.iter().map(|m|json!({"locus":m.locus,"model_revision":m.model_revision,"training_allele_count":m.training_alleles.len(),"omitted_alleles":m.omitted_alleles})).collect::<Vec<_>>(),"engine":r.engine,"context_revision":r.context_revision,"selection":r.analysis.selection,"columns":r.column_count,"ungapped_residues":r.residue_columns.len(),"qualified_pairs_in_retained_v":r.analysis.evidence.qualified_pairs_in_retained_v,"retained_candidates":r.analysis.tracks.len(),"v_retention_truncated":r.analysis.evidence.v_retention_truncated,"pair_retention_truncated":r.analysis.evidence.pair_retention_truncated,"preview_truncated":r.analysis.tracks.len()>20,"candidates":r.analysis.tracks.iter().take(20).map(|t|json!({"candidate_id":t.candidate_id,"locus":t.locus,"unassigned_residues":t.unassigned_residues,"unassigned_regions":t.unassigned_regions,"issues":t.issues,"first_supported_row_residue":t.first_supported_source_index,"row_residue_end":t.source_end})).collect::<Vec<_>>(),"artifact_id":artifact,"authority":"Conditional explicitly selected reference/rule evidence. Alignment gaps have no residue label; no calibrated domain/species classification or independent biological numbering certification."}))
        }).await
    }
    #[tool(
        description = "Describe the current explicit Human or Mouse C57BL/6 IGH reference scopes, IMGT numbering/region rules, dataset provenance and explicit resource limits. Human scope includes supported IMGT/Kabat/Chothia/North rules; mouse C57BL/6 is heavy-only IMGT/IMGT. Other species, mouse light chains and unsupported mouse schemes are refused. Not calibrated species/domain classification."
    )]
    pub(crate) async fn list_antibody_profiles(&self) -> Result<CallToolResult, ErrorData> {
        use collagene_core::antibody_numbering as n;
        let selections = [
            n::v_mapping::ReferenceSelection::human_all(),
            n::v_mapping::ReferenceSelection {
                species: n::references::Species::MouseC57Bl6,
                loci: vec![n::v_mapping::Locus::Heavy],
            },
        ];
        let mut profiles = Vec::new();
        for references in selections {
            let corpus = n::references::load(&references.species).map_err(bad)?;
            let human = references.species == n::references::Species::Human;
            profiles.push(json!({
                "references":references,"corpus_revision":corpus.revision().map_err(bad)?,
                "source_policy":corpus.source_policy(),"j_label_rule":corpus.j_label_rule(),"sources":corpus.sources(),
                "v_references":corpus.v_entries().len(),"j_references":corpus.j_entries().len(),
                "numbering_schemes":if human {vec!["Imgt","Kabat1991","Chothia1997"]}else{vec!["Imgt"]},
                "region_definitions":if human {vec!["Imgt","Kabat1991","ChothiaCanonical","North2011"]}else{vec!["Imgt"]},
            }));
        }
        output(json!({
            "engine":n::domain_tracks::TRACK_ENGINE,"rule_revision":n::domain_tracks::catalog_revision(),
            "profile_count":profiles.len(),"profiles":profiles,
            "model_identity_scope":"Actual model revisions depend on declared scoring and are returned with each analysis; these are explicit reference choices, not species predictions.",
            "max_source_residues":n::MAX_SOURCE_RESIDUES,"max_candidates":n::domain_candidates::MAX_CANDIDATES,
            "max_endpoints_per_locus":n::profile::MAX_ENDPOINTS,"max_cells":n::profile::MAX_CELLS,
            "max_trace_states":n::profile::MAX_TRACE_STATES,"max_scan_windows":n::domain_candidates::MAX_SCAN_WINDOWS,
            "authority":"Conditional selected-reference/rule evidence only. Alternatives and unresolved positions remain explicit. No calibrated classification, biological numbering certification or cloning/material authority."
        }))
    }
    #[tool(
        description = "Analyze an exact immutable protein document with explicit V scoring, CDR3/FR4 scope and work/retention limits. Returns conditional source-coordinate tracks under explicit numbering/region choices, all-optimal invariant framework evidence and explicit unresolved residues/CDRs. Per-endpoint optimal V paths only; truncation is not absence. Full source/options-bound report artifact plus at most 20 candidate summaries; no source/history/material edits. Not independent biological numbering certification or calibrated domain/species classification."
    )]
    pub(crate) async fn analyze_antibody(
        &self,
        Parameters(args): Parameters<AntibodyQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            use collagene_core::antibody_numbering::domain_tracks as n;
            let molecule={let store=server.session.store.lock().unwrap();let m=store.document(&args.document)?;n::request_revision(m,&args.options,&args.limits,&args.selection).map_err(bad)?;m.clone()};
            let r=n::analyze(&molecule,&args.options,&args.limits,&args.selection,&mut || false).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&serde_json::to_value(&r).map_err(bad)?)?;
            let previews=r.tracks.iter().take(20).map(|t|json!({"candidate_id":t.candidate_id,"locus":t.locus,"numbering_scheme":t.numbering_scheme,"region_definition":t.region_definition,"invariant_position_one":t.invariant_position_one,"first_supported_source_index":t.first_supported_source_index,"source_end":t.source_end,"unassigned_residues":t.unassigned_residues,"unassigned_regions":t.unassigned_regions,"issues":t.issues})).collect::<Vec<_>>();
            output(json!({"document":args.document,"references":args.options.references,"models":r.evidence.v_evidence.mappings.iter().map(|m|json!({"locus":m.locus,"model_revision":m.model_revision,"training_allele_count":m.training_alleles.len(),"omitted_alleles":m.omitted_alleles})).collect::<Vec<_>>(),"engine":r.engine,"context_revision":r.context_revision,"selection":r.selection,"source_revision":r.source_revision,"corpus_revision":r.evidence.v_evidence.corpus_revision,"qualified_pairs_in_retained_v":r.evidence.qualified_pairs_in_retained_v,"retained_candidates":r.tracks.len(),"v_retention_truncated":r.evidence.v_retention_truncated,"pair_retention_truncated":r.evidence.pair_retention_truncated,"preview_truncated":r.tracks.len()>20,"candidates":previews,"artifact_id":artifact,"authority":"Conditional explicitly selected reference/rule evidence. Unassigned residues and unsupported loops remain unresolved. Not independent biological numbering certification or calibrated domain/species classification."}))
        }).await
    }
    #[tool(
        description = "List the current embedded protein liability-screen rules, literal masks, citations, interpretation and resource bounds. Patterns are candidate sequence reviews, not measured PTMs, antibody identity, risk probabilities or automatic domain/numbering/CDR predictions. Optional supplied CDR intervals must bind an exact source revision."
    )]
    pub(crate) async fn list_protein_liability_profiles(
        &self,
    ) -> Result<CallToolResult, ErrorData> {
        use collagene_core::protein_liabilities as native;
        output(
            json!({"engine":native::ENGINE,"catalog_revision":native::catalog_revision(),"rules":native::rules(),"max_residues":native::MAX_RESIDUES,"max_hits":native::MAX_HITS,"max_provided_cdrs":native::MAX_CDRS,"indexing":"zero-based half-open stored-source coordinates; CDR membership uses the motif anchor","authority":"Literal certain/possible pattern matches only; no measured chemical state, antibody numbering, probabilities, activity or cloning/material authority."}),
        )
    }
    #[tool(
        description = "Freshly screen an exact immutable protein document using the embedded versioned liability profile. Requires explicit retained-hit limit (0 is counts-only); exact counts continue across the full source even after the retention cap. B/Z/J/X uncertainty produces possible/certain sequence-pattern matches, not PTM probabilities. Optional supplied CDRs have exact source, label/provenance and nonoverlapping interval validation; they are not validated numbering and do not assert complete CDR coverage. Returns at most 100 hit summaries plus a full source/options-bound JSON artifact. No source/history edits, saved declarations, candidate or material installs."
    )]
    pub(crate) async fn screen_protein_liabilities(
        &self,
        Parameters(args): Parameters<ProteinLiabilityQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            use collagene_core::protein_liabilities as native;
            let molecule={let store=server.session.store.lock().unwrap();let m=store.document(&args.document)?;native::request_revision(m,&args.options).map_err(bad)?;m.clone()};
            let report=native::screen(&molecule,&args.options,&mut||false).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&serde_json::to_value(&report).map_err(bad)?)?;
            output(json!({"document":args.document,"engine":report.engine,"catalog_revision":report.catalog_revision,"context_revision":report.context_revision,"source_length":report.source_length,"counts":report.counts,"retained_hit_count":report.hits.len(),"hits_truncated":report.hits_truncated,"display_truncated":report.hits.len()>100,"hits":report.hits.iter().take(100).collect::<Vec<_>>(),"ambiguous_residues":report.ambiguous_residues,"stop_symbols":report.stop_symbols,"cdr_rule_unassessed":report.cdr_rule_unassessed,"artifact_id":artifact,"authority":"Sequence-pattern candidates only. Supplied CDR context is documentary, not validated antibody numbering or complete domain coverage. No measured PTM or biochemical risk prediction."}))
        }).await
    }
    #[tool(
        description = "List the current embedded RNA covalent-chemistry model, supported base and imported annotation profiles, declaration version, catalog binding and bounds. Coverage and termini must be explicitly supplied; no default unmodified/closed assumption. This models declared neutral-acid composition, not measured sample state or enzyme behavior."
    )]
    pub(crate) async fn list_rna_chemistry_profiles(&self) -> Result<CallToolResult, ErrorData> {
        use collagene_core::rna_chemistry as rna;
        output(
            json!({"engine":rna::ENGINE,"declaration_version":rna::DECLARATION_VERSION,"catalog_revision":rna::catalog_revision(),"nucleosides":rna::profiles(),"annotation_profiles":rna::annotation_profiles(),"mass_model":rna::MASS_MODEL,"max_bases":rna::MAX_BASES,"max_sites":rna::MAX_SITES,"max_links":rna::MAX_LINKS,"max_declaration_bytes":rna::MAX_DECLARATION_BYTES,"indexing":"zero-based stored strand; link.after names upstream residue","authority":"Declared covalent hypothesis only; neutral fully protonated acid without counterions/adducts. Partial coverage and unresolved chemistry withhold aggregate quantities. No material, measured-state or cloning authority."}),
        )
    }
    #[tool(
        description = "Inspect retained RNA modification features on an exact immutable session document. Returns source-feature fingerprints and typed candidate residue declarations for a bounded INSDC modified_base profile. Suggestions are not accepted reviews. Unsupported codes, uncertain/compound/complement locations, edited features and conflicting parents remain unresolved. Preview bounded to 100 entries/80 code characters; complete evidence is a JSON artifact. No document/history/material mutation."
    )]
    pub(crate) async fn inspect_rna_annotations(
        &self,
        Parameters(args): Parameters<RnaSourceQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            use collagene_core::rna_chemistry as rna;
            let molecule={let store=server.session.store.lock().unwrap();let m=store.document(&args.document)?;rna::source_revision(m).map_err(bad)?;m.clone()};
            let candidates=rna::annotation_candidates(&molecule).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"document":args.document,"engine":rna::ENGINE,"catalog_revision":rna::catalog_revision(),"candidates":candidates}))?;
            let preview=candidates.iter().take(100).map(|c| json!({"feature_index":c.feature_index,"feature_revision":c.feature_revision,"position":c.position,"raw_code_preview":c.raw_code.as_ref().map(|s|s.chars().take(80).collect::<String>()),"raw_code_truncated":c.raw_code.as_ref().is_some_and(|s|s.chars().count()>80),"suggested_site":c.suggested_site,"unsupported_reason":c.unsupported_reason})).collect::<Vec<_>>();
            output(json!({"document":args.document,"engine":rna::ENGINE,"count":candidates.len(),"preview_truncated":candidates.len()>100,"candidates":preview,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Freshly analyze an explicit versioned RNA chemistry declaration against an exact immutable RNA source using embedded Rust. Requires explicit coverage and termini. Checks exact feature review fingerprints and base/sugar constraints, source parents, link positions and chemical contradictions before deriving CHNOPS/F composition and neutral-acid mass. Partial, ambiguous, unknown or unreviewed evidence withholds quantities; PS stereochemistry remains separate. Full source-bound report artifact and bounded diagnostics; no saved declaration, source edit, candidate or material install. Not experimental proof, activity prediction or cloning verification."
    )]
    pub(crate) async fn analyze_rna_chemistry(
        &self,
        Parameters(args): Parameters<RnaChemistryQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            use collagene_core::rna_chemistry as rna;
            let molecule={let store=server.session.store.lock().unwrap();let m=store.document(&args.document)?;rna::declaration_revision(m,&args.declaration).map_err(bad)?;m.clone()};
            let report=rna::analyze(&molecule,&args.declaration).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&serde_json::to_value(&report).map_err(bad)?)?;
            output(json!({"document":args.document,"engine":report.engine,"catalog_revision":report.catalog_revision,"declaration_revision":report.declaration_revision,"context_revision":report.context_revision,"formula":report.formula,"mass":report.mass,"formula_complete":report.formula_complete,"identity_complete":report.identity_complete,"issues":report.issues,"unspecified_stereocenter_count":report.unspecified_stereocenters.len(),"unspecified_stereocenter_preview":report.unspecified_stereocenters.iter().take(100).collect::<Vec<_>>(),"stereocenter_preview_truncated":report.unspecified_stereocenters.len()>100,"retained_evidence_count":report.retained_evidence.len(),"unreviewed_evidence_count":report.unreviewed_evidence.len(),"artifact_id":artifact,"authority":"Declared covalent hypothesis; not measured sample, material or cloning authority."}))
        }).await
    }
    #[tool(
        description = "Import a complete explicit ASCII PDB/mmCIF coordinate file up to 4 MiB using embedded Rust. Keeps exact original text and a full source/options-bound atom/scene report as immutable session artifacts, reserved atomically. Preserves separate model, label/author, alternate/insertion and deposited connection evidence. Optional typed display selection; null chooses the lowest model ID. Returns a source handle and bounded original-coordinate preview. Does not create molecule documents, infer a sequence link, predict a fold or verify cloning. No host paths, URLs or external software."
    )]
    pub(crate) async fn import_structure(
        &self,
        Parameters(args): Parameters<StructureImport>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            if args.text.len() > 4 * 1024 * 1024 {
                return Err(bad("MCP structure source exceeds 4 MiB"));
            }
            let report = collagene_core::protein::structure::prepare(
                args.text.as_bytes(),
                args.format,
                args.selection,
            )
            .map_err(bad)?;
            server.session.structure_output(args.text, &report)
        })
        .await
    }
    #[tool(
        description = "Freshly inspect an exact retained structure source handle with a typed model/chain/alternate/H-D selection. Reimports original PDB/mmCIF bytes and recomputes bounded scene geometry in Rust; a cached display report is never replayed as fresh evidence. Returns a bounded atom preview plus immutable original/full-report artifacts. Refuses stale, unknown, foreign-role source handles, absent model/selection and unsupported syntax. Does not mutate sequence/history, install a candidate/material or infer structure-to-sequence correspondence."
    )]
    pub(crate) async fn inspect_structure(
        &self,
        Parameters(args): Parameters<StructureInspect>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            if args.source.id != args.source.revision {
                return Err(bad("stale structure source revision"));
            }
            let (text, format) = {
                use collagene_core::protein::structure::coordinates::Format;
                let store = server.session.store.lock().unwrap();
                let a = store
                    .artifacts
                    .get(&args.source.id)
                    .ok_or_else(|| bad("unknown structure source"))?;
                let format = match (a.mime_type, a.filename) {
                    ("chemical/x-pdb", "structure.pdb") => Format::Pdb,
                    ("chemical/x-mmcif", "structure.cif") => Format::Mmcif,
                    _ => return Err(bad("artifact is not a retained structure source")),
                };
                let ArtifactData::Text(text) = &a.data else {
                    return Err(bad("artifact is not an original text structure"));
                };
                if text.len() > 4 * 1024 * 1024
                    || blake3::hash(text.as_bytes()).to_hex().as_str() != args.source.revision
                {
                    return Err(bad("structure source integrity/bound mismatch"));
                }
                (text.clone(), format)
            };
            let report = collagene_core::protein::structure::prepare(
                text.as_bytes(),
                format,
                args.selection,
            )
            .map_err(bad)?;
            server.session.structure_output(text, &report)
        })
        .await
    }
    #[tool(
        description = "Search an explicit bounded collection of immutable DNA/RNA document references using native affine local similarity. Returns one best observed score per subject, original reverse/circular coordinates, ambiguity and explicit score completeness. Exact canonical matches prove a score upper bound; capped origins/work remain incomplete. This is not BLAST statistics, sample identity, cloning verification or a material/candidate install. Full source-bound evidence is retained as a JSON artifact."
    )]
    pub(crate) async fn search_sequence_library(
        &self,
        Parameters(args): Parameters<SimilarityQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            use collagene_core::sequence::similarity as native;
            if args.subjects.is_empty()||args.subjects.len()>native::MAX_SUBJECTS{return Err(bad("Choose 1–128 immutable subject references."));}
            let (query,subjects)={
                let store=server.session.store.lock().unwrap();let query=store.document(&args.query)?;
                let sources=args.subjects.iter().map(|r|store.document(r).map(|m|(r.id.as_str(),m))).collect::<Result<Vec<_>,_>>()?;
                native::preflight_borrowed((&args.query.id,query),&sources,&args.options).map_err(bad)?;
                (native::Source{key:args.query.id.clone(),molecule:query.clone()},sources.into_iter().map(|(id,m)|native::Source{key:id.into(),molecule:m.clone()}).collect::<Vec<_>>())
            };
            let report=native::search(&query,&subjects,&args.options,&mut||false).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"query":args.query,"subjects":args.subjects,"report":report}))?;
            output(json!({"engine_version":report.engine_version,"context_revision":report.context_revision,"query":args.query,"complete":report.complete,"hits_truncated":report.hits_truncated,"work_used":report.work_used,"hits":report.hits,"subjects":report.subjects.iter().map(|s|json!({"key":s.key,"revision":s.revision,"best_score_proven":s.best_score_proven,"alignments_evaluated":s.alignments_evaluated,"limits":s.limits})).collect::<Vec<_>>(),"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Run embedded affine-gap pairwise alignment on two immutable nucleotide revisions. Explicit global/local/query-in-reference mode, scoring and matrix budget; no circular-origin/orientation search in this tool. Ambiguous-compatible columns are distinct from certain matches. Full CIGAR/differences in JSON artifact."
    )]
    pub(crate) async fn align_sequences(
        &self,
        Parameters(args): Parameters<AlignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (reference,query)={let store=server.session.store.lock().unwrap();(store.document(&args.reference)?.clone(),store.document(&args.query)?.clone())};
            let report=alignment::align(&reference.sequence,&query.sequence,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"reference":args.reference,"query":args.query,"engine_version":verification::ENGINE_VERSION,"alignment":report}))?;
            output(json!({"reference":args.reference,"query":args.query,"algorithm":report.algorithm,"options":report.options,"score":report.score,"reference_range":[report.reference_start,report.reference_end],"query_range":[report.query_start,report.query_end],"stats":report.stats,"cigar_preview":report.cigar.chars().take(2000).collect::<String>(),"cigar_truncated":report.cigar.len()>2000,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Compare immutable molecule sequences with bounded circular-origin and reverse-query candidate search. Returns tested counts, completeness and optimality flags separately; seed/budget-limited search is never labeled exhaustive. Exact identity does not establish biological function or uniqueness."
    )]
    pub(crate) async fn compare_molecules(
        &self,
        Parameters(args): Parameters<CompareQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (reference,query)={let store=server.session.store.lock().unwrap();(store.document(&args.reference)?.clone(),store.document(&args.query)?.clone())};
            let report=alignment::compare_molecules(&reference,&query,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&json!({"reference":args.reference,"query":args.query,"engine_version":verification::ENGINE_VERSION,"comparison":report}))?;
            output(json!({"reference":args.reference,"query":args.query,"reference_origin":report.reference_origin,"query_reversed":report.query_reversed,"score":report.alignment.score,"stats":report.alignment.stats,"evaluated_candidates":report.evaluated_candidates,"exhaustive_reference_origins":report.exhaustive_reference_origins,"exhaustive_query_orientations":report.exhaustive_query_orientations,"optimality_proven":report.optimality_proven,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
}
