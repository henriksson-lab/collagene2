use crate::session::presentation::summary;
pub mod annotation_merge;
pub mod reference_edit;
pub mod sbol_tools;
pub mod sequence_exchange_tools;

use crate::server::CollageneServer;
use crate::server::EmptyRequest;
use crate::server::response::bad;
use crate::server::response::bounded_json;
use crate::server::response::output;
use crate::session::Artifact;
use crate::session::ArtifactData;
use crate::session::MAX_DOCUMENT_REVISIONS;
use crate::session::MAX_ITEMS;
use crate::session::MAX_STORE_BYTES;
use crate::session::StoredSelection;
use crate::session::references::CollectionRef;
use crate::session::references::DocumentRef;
use crate::session::references::EvidenceRef;
use crate::tools::catalog::catalog_tools;

use std::collections::BTreeMap;

use base64::Engine as _;
use base64::engine::general_purpose::STANDARD;
use collagene_core::formats::annotations as annotation_formats;
use collagene_core::formats::oligos as oligo_formats;
use collagene_core::formats::sequence as formats;
use collagene_core::sequence::analysis;
use collagene_core::sequence::annotations::library as feature_library;
use rmcp::ErrorData;
use rmcp::handler::server::wrapper::Parameters;
use rmcp::model::CallToolResult;
use rmcp::tool;
use rmcp::tool_router;
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::Value;
use serde_json::json;

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
    pub options: collagene_core::sequence::search::SearchOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ViewTracksQuery {
    pub document: DocumentRef,
    pub window: collagene_core::tracks::view::TrackWindow,
    pub options: collagene_core::tracks::view::TrackOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditTemplateQuery {
    pub source: DocumentRef,
    pub target: DocumentRef,
    pub options: collagene_core::sequence::template::EditTemplateOptions,
    #[serde(default)]
    pub guide: Option<collagene_core::sequence::template::GuideBinding>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EditTemplateCheckQuery {
    pub source: DocumentRef,
    pub target: DocumentRef,
    pub donor: DocumentRef,
    pub options: collagene_core::sequence::template::EditTemplateOptions,
    #[serde(default)]
    pub guide: Option<collagene_core::sequence::template::GuideBinding>,
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
pub struct AnnotationImport {
    pub document: DocumentRef,
    pub text: String,
    /// Explicit annotation syntax: Gff3, Gtf or Bed. No name-based reference mapping.
    pub format: annotation_formats::AnnotationFormat,
    pub seqid: String,
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
pub struct DocumentPage {
    #[serde(default)]
    pub start: usize,
    /// Default 128; at most 128 summaries per response.
    #[serde(default)]
    pub limit: Option<usize>,
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
pub struct ArtifactRegion {
    pub artifact_id: String,
    /// UTF-8 byte offset, must be on a character boundary.
    pub start: usize,
    /// Maximum 16,384 UTF-8 bytes; next_start is returned for pagination.
    pub max_bytes: usize,
}
#[tool_router(router = sequence_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Edit one immutable reference document: reverse-complement, rotate a circular origin, extract or replace a zero-based half-open interval, set/delete an annotation, flip annotation strands, rename or change reference topology. Uses validated native coordinate/qualifier projection and annotation-edit guards. Publishes a new session-only reference revision; the source stays unchanged. This is reference intent only, never physical cleavage/ligation, recorded workspace ancestry, or verified cloning material. Unsupported coordinate semantics fail rather than silently moving source coordinates."
    )]
    pub(crate) async fn edit_reference_document(
        &self,
        Parameters(args): Parameters<reference_edit::Request>,
    ) -> Result<CallToolResult, ErrorData> {
        if matches!(
            &args.operation,
            collagene_core::documents::reference_edit::Operation::MergeAnnotations { .. }
        ) {
            return Err(bad(
                "Use apply_annotation_merge with immutable session source references; raw merge payloads do not resolve source ownership.",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let source = {
                let store = server.session.store.lock().unwrap();
                let source = store.document(&args.document)?;
                collagene_model::bounds::model_preflight(&[source]).map_err(bad)?;
                source.validate().map_err(bad)?;
                args.operation.preflight().map_err(bad)?;
                source.clone()
            };
            let operation = args.operation.name();
            let product = args.operation.apply(&source).map_err(bad)?;
            collagene_model::bounds::model_preflight(&[&product]).map_err(bad)?;
            let review_count = product.features.iter().filter(|feature| feature.needs_review).count();
            let mut warnings = vec!["Reference edit only: no physical product, reusable material origin, native history or current scientific verification is established.".into()];
            if review_count > 0 {
                warnings.push(format!("Result contains {review_count} annotations marked for review; inspect their locations and source-evidence notes before using their biological claims."));
            }
            let published = server.session.insert_documents(vec![product], warnings)?;
            let Some(mut receipt) = published.structured_content else {
                return Err(ErrorData::internal_error("Missing session document receipt", None));
            };
            receipt["source"] = json!(args.document);
            receipt["operation"] = json!(operation);
            receipt["reference_only"] = json!(true);
            receipt["review_required"] = json!(review_count > 0);
            receipt["annotations_requiring_review"] = json!(review_count);
            output(receipt)
        }).await
    }
    #[tool(
        description = "Design a literal source-derived homology editing donor against an independently supplied immutable target. Explicit edit/deletion/replacement/arm/orientation options and an optional exact source-bound guide are independently checked base by base. The target is a noncausal requirement, never donor source DNA. Returns bounded donor/edited-reference previews and full immutable JSON; installs an immutable donor document only when the independent check is consistent, atomically with the artifact. Sequence consistency is not cellular HDR, nuclease activity or editing success. No physical operation or persistent application writes."
    )]
    pub(crate) async fn design_edit_template(
        &self,
        Parameters(args): Parameters<EditTemplateQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let(source,target)={let store=server.session.store.lock().unwrap();let source=store.document(&args.source)?;let target=store.document(&args.target)?;
                collagene_core::sequence::template::source_preflight(source,target,&args.options,args.guide.as_ref()).map_err(bad)?;(source.clone(),target.clone())};
            let report=collagene_core::sequence::template::design(&source,&target,&args.options,args.guide.as_ref()).map_err(bad)?;
            report.validate_bindings_for(&source,&target,&args.options,args.guide.as_ref()).map_err(bad)?;
            let donor=report.check.consistent.then(||report.donor.clone());
            let(artifact_id,donor_document)=server.session.edit_template_artifact(&report,donor)?;
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
    pub(crate) async fn check_edit_template(
        &self,
        Parameters(args): Parameters<EditTemplateCheckQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let(source,target,donor)={let store=server.session.store.lock().unwrap();let source=store.document(&args.source)?;let target=store.document(&args.target)?;
                collagene_core::sequence::template::source_preflight(source,target,&args.options,args.guide.as_ref()).map_err(bad)?;
                let donor=store.document(&args.donor)?;
                if donor.sequence.alphabet()!=collagene_model::Alphabet::Dna || donor.topology!=collagene_model::Topology::Linear || donor.sequence.len()>collagene_core::sequence::template::MAX_DONOR{return Err(bad("donor candidate must be bounded linear DNA"));}
                (source.clone(),target.clone(),donor.sequence.as_str().to_owned())};
            let report=collagene_core::sequence::template::check(&source,&target,&args.options,args.guide.as_ref(),&donor).map_err(bad)?;
            report.validate_bindings_for(&source,&target,&args.options,args.guide.as_ref(),&donor).map_err(bad)?;
            let artifact_id=server.session.text_artifact(bounded_json(&report)?,"application/json","edit-template-check.json")?;
            output(json!({"source":args.source,"target":args.target,"donor_document":args.donor,"check":report,"artifact_id":artifact_id,
                "authority":"Independent literal source/edit/candidate/target sequence check; no cellular repair, empirical activity, physical product or inherited verification authority."}))
        }).await
    }
    #[tool(
        description = "Search an immutable document with explicit sequence, annotation, saved-primer, enzyme or all-source mode. DNA/RNA IUPAC and peptide search preserve document alphabet. Bounded hits/work/rows, strand, circular spans, exact source revision and completeness remain explicit. Readonly: selecting or modifying a document is a separate user action. Full report is a JSON artifact."
    )]
    pub(crate) async fn search_document(
        &self,
        Parameters(args): Parameters<SearchDocumentQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let report=collagene_core::sequence::search::search_document(&molecule,&args.options).map_err(bad)?;
            let artifact=server.session.analysis_artifact(&serde_json::to_value(&report).map_err(bad)?)?;
            output(json!({"document":args.document,"source_revision":report.source_revision,"engine_version":report.engine_version,"options":report.options,"complete":report.complete,"work":report.work,"rows":report.rows,"hit_count":report.hits.len(),"hits":report.hits.iter().take(100).map(|h|json!({"id":h.id,"kind":h.kind,"label":h.label.chars().take(256).collect::<String>(),"label_preview_truncated":h.label.chars().count()>256,"parts":h.parts.iter().take(64).collect::<Vec<_>>(),"parts_preview_truncated":h.parts.len()>64,"feature_index":h.feature_index,"enzyme":h.enzyme,"top_cut":h.top_cut,"bottom_cut":h.bottom_cut,"certain":h.certain})).collect::<Vec<_>>(),"preview_truncated":report.hits.len()>100,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    #[tool(
        description = "Inspect bounded annotation, saved-primer, enzyme/cut and translation tracks for an exact immutable document and half-open viewport. Circular windows may cross the origin once. Uses the same deterministic track engine as the GUI; returns a bounded preview and full JSON artifact without changing selections or annotations."
    )]
    pub(crate) async fn view_tracks(
        &self,
        Parameters(args): Parameters<ViewTracksQuery>,
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
                collagene_core::tracks::view::view_tracks(&molecule, args.window, &args.options)
                    .map_err(bad)?;
            let value = serde_json::to_value(&report).map_err(bad)?;
            let artifact = server.session.analysis_artifact(&value)?;
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
    pub(crate) async fn create_feature_library(
        &self,
        Parameters(args): Parameters<CreateFeatureLibrary>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.features.len() > 10_000 {
            return Err(bad("at most 10000 source feature selections"));
        }
        let server = self.clone();
        self.compute(move || {
            let sources = {
                let store = server.session.store.lock().unwrap();
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
            server.session.insert_feature_library(library)
        })
        .await
    }
    #[tool(
        description = "Import a validated versioned feature-library JSON file encoded as base64 (at most 8 MiB decoded). Sequence, source identities, entry revisions and supported semantics are validated. Untrusted labels remain data; importing does not annotate documents."
    )]
    pub(crate) async fn import_feature_library(
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
            server
                .session
                .insert_feature_library(feature_library::load_library(&bytes).map_err(bad)?)
        })
        .await
    }
    #[tool(
        description = "List immutable session feature libraries with exact references and bounded source/version summaries. No filesystem or unsynchronized browser access."
    )]
    pub(crate) fn list_feature_libraries(
        &self,
        Parameters(_): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.session.store.lock().unwrap();
        output(
            json!({"libraries":store.feature_libraries.iter().map(|(id,l)|json!({"id":id,"revision":id,"name":l.name.chars().take(256).collect::<String>(),"version":l.version.chars().take(256).collect::<String>(),"entry_count":l.entries.len()})).collect::<Vec<_>>()}),
        )
    }
    #[tool(
        description = "Export a complete immutable feature library as validated JSON evidence. Exact library id/revision required; artifact chunks retain every source and entry rather than only previews."
    )]
    pub(crate) async fn export_feature_library(
        &self,
        Parameters(args): Parameters<FeatureLibraryQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let library = server
                .session
                .store
                .lock()
                .unwrap()
                .feature_library(&args.library)?
                .clone();
            let bytes = feature_library::save_library(&library).map_err(bad)?;
            let length = bytes.len();
            let text = String::from_utf8(bytes).map_err(bad)?;
            let artifact =
                server
                    .session
                    .text_artifact(text, "application/json", "feature-library.json")?;
            output(json!({"library":args.library,"bytes":length,"artifact_id":artifact}))
        })
        .await
    }
    #[tool(
        description = "Scan an exact DNA document against an immutable feature library on explicit strands with bounded mismatch/work/hit limits. Full-entry coverage, concrete versus ambiguity-compatible identities, repeats and completeness are explicit. Returns a revision-bound proposal set and JSON artifact; no annotations are automatically installed."
    )]
    pub(crate) async fn scan_feature_library(
        &self,
        Parameters(args): Parameters<FeatureLibraryScan>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let (document,library)={let store=server.session.store.lock().unwrap();(store.document(&args.document)?.clone(),store.feature_library(&args.library)?.clone())};
            let report=feature_library::scan(&document,&library,&args.options).map_err(bad)?;
            let text=serde_json::to_string(&report).map_err(bad)?;let id=blake3::hash(text.as_bytes()).to_hex().to_string();
            let result=json!({"document":args.document,"library":args.library,"scan":{"id":id,"revision":id},"source_revision":report.source_revision,"library_revision":report.library_revision,"engine_version":report.engine_version,"complete":report.complete,"comparisons":report.comparisons,"options":report.options,"hit_count":report.hits.len(),"hits":report.hits.iter().take(100).map(|hit|json!({"id":hit.id,"entry_id":hit.entry_id,"entry_revision":hit.entry_revision,"start":hit.start,"length":hit.length,"reverse":hit.reverse,"concrete_matches":hit.concrete_matches,"compatible_ambiguous":hit.compatible_ambiguous,"mismatches":hit.mismatches,"identity_lower":hit.identity_lower,"identity_upper":hit.identity_upper,"coverage":hit.coverage,"same_label_geometry_hint":hit.same_label_geometry_hint,"name":hit.feature.label().chars().take(256).collect::<String>(),"kind":hit.feature.kind,"needs_review":hit.feature.needs_review})).collect::<Vec<_>>(),"hits_truncated":report.hits.len()>100,"warnings":report.warnings.iter().take(100).map(|w|w.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),"warnings_truncated":report.warnings.len()>100,"artifact_id":id,"applied":false});
            let mut store=server.session.store.lock().unwrap();let new_report=!store.feature_scans.contains_key(&id);let new_artifact=!store.artifacts.contains_key(&id);
            if new_report&&store.feature_scans.len()>=MAX_ITEMS{return Err(bad("at most 128 feature scans"));}if new_artifact&&store.artifacts.len()>=MAX_ITEMS{return Err(bad("at most 128 artifacts"));}
            let bytes=text.len().checked_mul(usize::from(new_report)+usize::from(new_artifact)).ok_or_else(||bad("session accounting overflow"))?;store.reserve(bytes)?;
            store.feature_scans.entry(id.clone()).or_insert(report);store.artifacts.entry(id).or_insert(Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"feature-library-scan.json"});output(result)
        }).await
    }
    #[tool(
        description = "Parse GFF3, GTF or BED annotation evidence using an explicit source seqid and an exact imported reference revision. Coordinate/phase/strand/compound semantics are checked by the embedded parser; unrelated seqids remain explicit. Returns an immutable proposal dataset. No existing annotations are changed; <=8MiB source."
    )]
    pub(crate) async fn import_annotations(
        &self,
        Parameters(args): Parameters<AnnotationImport>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.text.len() > 8 * 1024 * 1024 {
            return Err(bad("annotation text exceeds 8 MiB"));
        }
        let server = self.clone();
        self.compute(move||{
            let document=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let report=annotation_formats::import_annotations(&document,&args.text,args.format,&args.seqid).map_err(bad)?;
            let bytes=serde_json::to_vec(&report).map_err(bad)?;let id=blake3::hash(&bytes).to_hex().to_string();
            let result=json!({"dataset":{"id":id,"revision":id},"document":args.document,"format":report.format,"seqid":report.seqid,"source_blake3":report.source_blake3,"reference_revision":report.reference_revision,"proposal_count":report.proposals.len(),"proposals":report.proposals.iter().take(100).enumerate().map(|(i,p)|json!({"index":i,"id":p.id,"name":p.feature.label().chars().take(256).collect::<String>(),"kind":p.feature.kind,"needs_review":p.feature.needs_review,"spans":p.feature.location.spans().iter().take(20).collect::<Vec<_>>(),"spans_truncated":p.feature.location.spans().len()>20,"source_rows":p.source_rows.iter().take(100).collect::<Vec<_>>(),"source_rows_truncated":p.source_rows.len()>100})).collect::<Vec<_>>(),"proposals_truncated":report.proposals.len()>100,"skipped_seqids":report.skipped_seqids.iter().take(100).map(|s|s.chars().take(256).collect::<String>()).collect::<Vec<_>>(),"skipped_seqids_truncated":report.skipped_seqids.len()>100,"warnings":report.warnings.iter().take(100).map(|w|w.chars().take(1024).collect::<String>()).collect::<Vec<_>>(),"warnings_truncated":report.warnings.len()>100,"applied":false});
            let mut store=server.session.store.lock().unwrap();if !store.annotation_datasets.contains_key(&id){if store.annotation_datasets.len()>=MAX_ITEMS{return Err(bad("at most 128 annotation datasets"));}store.reserve(bytes.len())?;store.annotation_datasets.insert(id,report);}output(result)
        }).await
    }
    #[tool(
        description = "Export supported annotations from an exact document as GFF3, GTF or BED with an explicit output seqid. Phase, discontinuous/reverse locations and format losses follow the embedded strict profile; unsupported semantics fail rather than being invented. Full text artifact includes loss diagnostics separately."
    )]
    pub(crate) async fn export_annotations(
        &self,
        Parameters(args): Parameters<AnnotationExport>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{let document=server.session.store.lock().unwrap().document(&args.document)?.clone();let report=annotation_formats::export_annotations(&document,args.format,&args.seqid).map_err(bad)?;let length=report.text.len();let filename=match args.format { annotation_formats::AnnotationFormat::Gff3=>"annotations.gff3",annotation_formats::AnnotationFormat::Gtf=>"annotations.gtf",annotation_formats::AnnotationFormat::Bed=>"annotations.bed" };let artifact=server.session.text_artifact(report.text,"text/plain",filename)?;output(json!({"document":args.document,"format":args.format,"seqid":args.seqid,"bytes":length,"warnings":report.warnings,"artifact_id":artifact}))}).await
    }
    #[tool(
        description = "Export the complete immutable GFF3/GTF/BED proposal dataset and exact original source bytes as two artifacts. Includes explicit reference/seqid mapping and every row/proposal; this is evidence, not an automatic document edit or trusted verification."
    )]
    pub(crate) async fn export_annotation_dataset(
        &self,
        Parameters(args): Parameters<AnnotationDatasetQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let report=server.session.store.lock().unwrap().annotation_dataset(&args.dataset)?.clone();
            let text=serde_json::to_string(&report).map_err(bad)?;let original=report.original_bytes().map_err(bad)?;
            let text_id=blake3::hash(text.as_bytes()).to_hex().to_string();let original_id=blake3::hash(&original).to_hex().to_string();
            let mut store=server.session.store.lock().unwrap();let new_text=!store.artifacts.contains_key(&text_id);let new_original=!store.artifacts.contains_key(&original_id);
            if let Some(existing)=store.artifacts.get(&original_id){if !matches!(existing.data,ArtifactData::Binary(_)){return Err(bad("original source artifact already has another representation"));}}
            let added=usize::from(new_text)+usize::from(new_original);if store.artifacts.len()+added>MAX_ITEMS{return Err(bad("annotation artifacts exceed item limit; no outputs created"));}
            let bytes=(if new_text{text.len()}else{0usize}).checked_add(if new_original{original.len()}else{0}).ok_or_else(||bad("session accounting overflow"))?;store.reserve(bytes)?;
            store.artifacts.entry(text_id.clone()).or_insert(Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"annotation-dataset.json"});
            store.artifacts.entry(original_id.clone()).or_insert(Artifact{data:ArtifactData::Binary(original),mime_type:"application/octet-stream",encoding:"base64",filename:match report.format { annotation_formats::AnnotationFormat::Gff3=>"original-annotations.gff3",annotation_formats::AnnotationFormat::Gtf=>"original-annotations.gtf",annotation_formats::AnnotationFormat::Bed=>"original-annotations.bed" }});
            output(json!({"dataset":args.dataset,"dataset_artifact_id":text_id,"original_artifact_id":original_id}))
        }).await
    }
    #[tool(
        description = "Import a CSV/TSV oligo collection using explicit column mapping. Retain the whole source table once, exact rows/modifications/units, and explicit selected source-row indices. Only supported selected primers become ssDNA reference documents; unsupported modifications stay evidence. Atomic; <=8 MiB source."
    )]
    pub(crate) async fn import_oligos(
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
            let mut store=server.session.store.lock().unwrap();
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
    pub(crate) async fn export_oligos(
        &self,
        Parameters(args): Parameters<OligoExport>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let table={let store=server.session.store.lock().unwrap();if args.collection.id!=args.collection.revision{return Err(bad("stale oligo collection revision"));}store.oligos.get(&args.collection.id).cloned().ok_or_else(||bad("unknown oligo collection"))?};
            if let Some(rows)=&args.selection{let mut unique=std::collections::BTreeSet::new();if rows.iter().any(|&r|r>=table.records.len()||!unique.insert(r)){return Err(bad("selected collection rows must be unique existing source rows"));}}
            let report=oligo_formats::export_oligos(&table,args.delimiter,args.mode).map_err(bad)?;
            // Preflight every output against aggregate limits before creating any artifact.
            let original=table.original_bytes().map_err(bad)?;
            let saved_selections={let store=server.session.store.lock().unwrap();store.oligo_selections.iter().filter(|(_,s)|s.collection_id==args.collection.id).map(|(id,s)|json!({"id":id,"source_rows":s.source_rows})).collect::<Vec<_>>()};
            let evidence=json!({"collection":args.collection,"selection":args.selection,"saved_selections":saved_selections,"table":table});
            let full=serde_json::to_string(&evidence).map_err(bad)?;
            let total=original.len().checked_add(full.len()).and_then(|v|v.checked_add(report.text.len())).ok_or_else(||bad("artifact accounting overflow"))?;
            {let store=server.session.store.lock().unwrap();if total>MAX_STORE_BYTES.saturating_sub(store.used)||store.artifacts.len()+3>MAX_ITEMS{return Err(bad("export artifacts exceed session size/item limit; no outputs created"));}}
            let original_id=server.session.binary_artifact(original,"application/octet-stream","original-oligos.table")?;
            let csv_id=server.session.text_artifact(report.text,"text/plain","oligos.table")?;
            let collection_id=server.session.text_artifact(full,"application/json","oligo-collection.json")?;
            output(json!({"collection":args.collection,"selected_source_rows":args.selection,"exports_all_rows":true,"original_artifact_id":original_id,"table_artifact_id":csv_id,"collection_artifact_id":collection_id,"warnings":report.warnings}))
        }).await
    }
    #[tool(
        description = "List the authenticated Addgene catalog's sequence choices for a plasmid URL/ID: deposited versus Addgene-verified, full versus partial, sequence ID and declared length. Configure COLLAGENE_ADDGENE_API_TOKEN in the server environment. Does not import or silently select a record; pass sequence_id to import_url after review. No raw sequence, secret or temporary download URL is returned."
    )]
    pub(crate) async fn list_source_sequences(
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
        let artifact = self
            .session
            .analysis_artifact(&serde_json::to_value(&catalog).map_err(bad)?)?;
        output(
            json!({"plasmid_id":catalog.plasmid_id(),"name":catalog.name(),"response_blake3":catalog.response_blake3(),"count":catalog.choices().len(),"display_truncated":catalog.choices().len()>100,"choices":catalog.choices().iter().take(100).map(|choice|json!({"sequence_id":choice.sequence_id,"description":choice.description.chars().take(512).collect::<String>(),"description_truncated":choice.description.chars().count()>512,"origin":choice.origin,"completeness":choice.completeness,"length":choice.length,"catalog_sequence_blake3":choice.catalog_sequence_blake3})).collect::<Vec<_>>(),"artifact_id":artifact}),
        )
    }
    #[tool(
        description = "Import one bounded public nucleotide record using embedded Rust HTTP: NCBI accession/GenBank URL, explicit ena:ACCESSION.VERSION or official ENA URL, explicit ddbj:ACCESSION.VERSION or official DDBJ URL, or direct HTTPS sequence file. Bare accessions remain NCBI and numeric IDs remain Addgene. ENA/DDBJ record identity, declared lengths, original source bytes and provenance are checked by the shared source adapter; lists/search/ranges are not whole-record imports. Addgene requires an approved server token and explicit catalog sequence_id, which is Addgene-only. Public HTTPS/redirect/DNS limits remain enforced; no filesystem paths, cookie forwarding or private targets. Existing 32 MiB download and ordinary molecule limits remain; this tool does not stream large remote genomes. NCBI/Addgene credentials stay provider-scoped and out of metadata."
    )]
    pub(crate) async fn import_url(
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
        self.session
            .insert_documents(report.molecules, report.warnings)
    }
    #[tool(
        description = "Export one exact session document revision as a native .cg2 artifact, preserving annotations and source evidence. Use this after native import_url when browser CORS prevents direct retrieval. Read base64 chunks, concatenate and decode before saving; opening the file in the browser is an explicit document import, not a project mutation."
    )]
    pub(crate) async fn export_document(
        &self,
        Parameters(reference): Parameters<DocumentRef>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move ||{
            let molecule=server.session.store.lock().unwrap().document(&reference)?.clone();
            let bytes=collagene_core::protocol::save_document(&molecule).map_err(bad)?;
            let id=server.session.binary_artifact(bytes.clone(),"application/octet-stream","sequence.cg2")?;
            output(json!({"artifact_id":id,"filename":"sequence.cg2","encoding":"base64","mime_type":"application/octet-stream","bytes":bytes.len(),"document":reference}))
        }).await
    }
    #[tool(
        description = "Review/export one immutable session revision as annotated GenBank, EMBL, GenPept or UniProt text using the native writers. GenBank/EMBL require explicit DNA/RNA documents; GenPept/UniProt require linear Protein. Protein profiles accept bounded explicit record identifiers and GenPept calendar export date; source accessions, review status, molecular mass and biological history are never invented. Returns exact writer warnings and no artifact until accept_omissions=true when warnings exist. Unsupported chemistry/locations/options remain errors even with acknowledgement. Output is bounded32MiB and session-only; retain native files/workspaces for exact source evidence and history."
    )]
    pub(crate) async fn export_annotated_sequence(
        &self,
        Parameters(args): Parameters<AnnotatedSequenceExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            use collagene_core::formats::protein as protein_formats;
use collagene_model::Alphabet;
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
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
            let artifact=server.session.text_artifact(report.text,"text/plain",filename)?;
            output(json!({"document":args.document,"source_revision":molecule.revision(),"state":"prepared","profile":profile,"warnings":report.warnings,"artifact_id":artifact,"bytes":bytes,"filename":filename,"mime_type":"text/plain","encoding":"utf8","persistence":"session only","authority":"Annotated exchange text is not a native history archive or new scientific verification."}))
        }).await
    }
    #[tool(
        description = "Import explicit FASTA, GCG/Wisconsin single-sequence or original Collagene XML profiles, at most 4 MiB of text. FASTA alphabet=null retains per-record declarations (unmarked records are DNA); an explicit alphabet handles unmarked RNA/protein and must agree with every declaration. GCG requires a caller-selected alphabet and verifies the declared type, length and checksum. XML lower-only migration requires AllowLowerOnly, retains exact original rows/bytes and never reconstructs physical cloning chemistry. Imports the entire selected project atomically as immutable session documents; no paths or external programs."
    )]
    pub(crate) async fn import_compat_sequence(
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
                    collagene_core::formats::sequence::gcg::import_gcg(
                        args.text.as_bytes(),
                        &options,
                    )
                }
                sequence_exchange_tools::ImportProfile::Fasta { alphabet } => {
                    let text = args.text.trim_start_matches('\u{feff}');
                    match alphabet {
                        Some(alphabet) => formats::import_fasta_with_alphabet(text, alphabet),
                        None => formats::import_fasta(text),
                    }
                }
                sequence_exchange_tools::ImportProfile::CollageneXml { policy } => {
                    collagene_core::formats::original::collagene_lower_strand::import_with_policy(
                        args.text.as_bytes(),
                        policy,
                    )
                }
            }
            .map_err(bad)?;
            server
                .session
                .insert_documents(report.molecules, report.warnings)
        })
        .await
    }
    #[tool(
        description = "Import up to 4 MiB of SBOL 3.1 structural JSON-LD as immutable documentary graph evidence. No remote contexts/URLs are fetched. Returns an exact original-source handle, bounded Component/linked-Sequence previews and a full typed graph/issue artifact. Does not choose a sequence, flatten hierarchy, create documents or restore cloning/material/verification authority. Use inspect_sbol for further pages and project_sbol for explicitly reviewed reference creation."
    )]
    pub(crate) async fn import_sbol(
        &self,
        Parameters(args): Parameters<sbol_tools::ImportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.sbol_import(args)).await
    }
    #[tool(
        description = "Freshly inspect a retained SBOL graph from an exact imported-source handle or exact native workspace attachment. Pages Components with their actual linked Sequence choices; full typed objects, foreign properties and issues remain in an artifact. No network, implicit choices, source edits or scientific authority. Page bounds are explicit."
    )]
    pub(crate) async fn inspect_sbol(
        &self,
        Parameters(args): Parameters<sbol_tools::InspectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.sbol_inspect(args)).await
    }
    #[tool(
        description = "Project one explicitly selected SBOL Component and one of its linked Sequences into a native reference. Requires exact graph source and native selection/topology policy; never assembles child sequences. Unsupported/unsafe physical types still refuse. Representation losses produce review_required with a full documentary report until accept_omissions=true. Creation atomically retains the reference, original graph attachment and full projection receipt in a NEW immutable native workspace. No target fitting, source edits or inherited verification/material authority. Preserve the workspace for source-graph fidelity."
    )]
    pub(crate) async fn project_sbol(
        &self,
        Parameters(args): Parameters<sbol_tools::ProjectQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.sbol_project(args)).await
    }
    #[tool(
        description = "Export explicitly requested Original or Expanded bytes from an exact retained SBOL graph source or exact native workspace attachment. Original preserves supplied bytes; Expanded preserves the complete supported retained graph including foreign properties, without RDF canonicalization. This is full local graph exchange, not a selected lab-library publication grant. No network, sequence flattening, source edits or scientific verification."
    )]
    pub(crate) async fn export_sbol_graph(
        &self,
        Parameters(args): Parameters<sbol_tools::GraphExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.sbol_export_graph(args)).await
    }
    #[tool(
        description = "Export an exact immutable native reference as SBOL 3.1 structural JSON-LD using explicit namespace/Component/Sequence identities. Reuses the existing core adapter and the shared reference-export restriction guard used by the browser. Losses require accept_omissions=true; physical, tracked, modified or reference-only state cannot be erased through acknowledgment. Emits reviewed design data only, never imported experimental history or material/verification authority. Source remains unchanged."
    )]
    pub(crate) async fn export_sbol(
        &self,
        Parameters(args): Parameters<sbol_tools::NativeExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || server.sbol_export_native(args)).await
    }
    #[tool(
        description = "Import exact local sequence-file bytes as base64, at most 8 MiB decoded, with an explicit SnapGene DNA/protein, XDNA v0, ABIF, SCF or current native .cg2 profile. Uses the same bounded core parsers as the browser; format mismatch refuses without text sniffing or fallback. Preserves supported annotations, original file/trace evidence and parser warnings in immutable session documents. XDNA requires a bounded fallback name. No paths, URLs, external programs, browser edits or restored physical material authority. Use import_compat_sequence for explicit text profiles."
    )]
    pub(crate) async fn import_sequence_file(
        &self,
        Parameters(args): Parameters<sequence_exchange_tools::FileImportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.data_base64.len() > sequence_exchange_tools::MAX_FILE_BYTES.div_ceil(3) * 4 {
            return Err(bad("Sequence file exceeds 8 MiB decoded"));
        }
        let server = self.clone();
        self.compute(move || {
            let report = sequence_exchange_tools::import_file(args).map_err(bad)?;
            server
                .session
                .insert_documents(report.molecules, report.warnings)
        })
        .await
    }
    #[tool(
        description = "Review/export an exact immutable DNA or selected linear-protein document under the current XDNA v0 reference profile. Returns structured losses without an artifact until accept_omissions=true when fields cannot be represented. Refuses unsupported modified/physical/reference-only evidence even with acknowledgment. Reuses the bounded core writer; no stale opaque chunks, source mutation or cloned experimental history. Artifact is session-only."
    )]
    pub(crate) async fn export_xdna(
        &self,
        Parameters(args): Parameters<sequence_exchange_tools::ReferenceExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            use collagene_core::formats::original::xdna_export::{self, Policy};
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let review=xdna_export::review(&molecule).map_err(bad)?;
            if !args.accept_omissions&&!review.losses.is_empty() {
                return output(json!({"document":args.document,"state":"review_required","review":review,"artifact_id":null}));
            }
            let report=xdna_export::export_xdna(&molecule,if args.accept_omissions{Policy::AllowOmissions}else{Policy::Strict}).map_err(bad)?;
            let size=report.bytes.len();
            let artifact=server.session.binary_artifact(report.bytes,"application/octet-stream","sequence.xdna")?;
            output(json!({"document":args.document,"state":"prepared","profile":report.profile,"source_revision":report.source_revision,"losses":report.losses,"warnings":report.warnings,"artifact_id":artifact,"bytes":size,"filename":"sequence.xdna","mime_type":report.mime,"encoding":"base64","persistence":"session only"}))
        }).await
    }
    #[tool(
        description = "Review/export regenerated SCF 3.00 from an exact immutable linear DNA document's retained SCF attachment and explicit precision. Freshly checks source hash, actual signal/calls/peaks/four-channel accuracies and declared orientation. Sequence edits, missing trace evidence, ABIF scalar scores, fabricated signal or precision loss refuse even with accept_omissions=true. Returns documentary losses without an artifact until acknowledged; no basecalling, rescaling or invented qualities. Source and original bytes remain unchanged."
    )]
    pub(crate) async fn export_trace_scf(
        &self,
        Parameters(args): Parameters<sequence_exchange_tools::ScfExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            use collagene_core::{formats::scf::{self as scf_export, Options, LossPolicy}, review::reads::traces};
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let trace=traces::trace_from_molecule(&molecule).map_err(bad)?;
            let report=scf_export::export(&trace,&Options{precision:args.precision,loss_policy:LossPolicy::AllowDocumentaryOmissions}).map_err(bad)?;
            if !args.accept_omissions&&!report.losses.is_empty() {
                return output(json!({"document":args.document,"state":"review_required","source_revision":molecule.revision(),"source_hash":report.source_hash,"orientation":report.orientation,"context_revision":report.context_revision,"sample_size":report.sample_size,"losses":report.losses,"warnings":report.warnings,"artifact_id":null}));
            }
            let size=report.bytes.len();
            let artifact=server.session.binary_artifact(report.bytes,"application/octet-stream","trace.scf")?;
            output(json!({"document":args.document,"state":"prepared","source_revision":molecule.revision(),"source_hash":report.source_hash,"orientation":report.orientation,"context_revision":report.context_revision,"sample_size":report.sample_size,"losses":report.losses,"warnings":report.warnings,"artifact_id":artifact,"bytes":size,"filename":"trace.scf","mime_type":"application/octet-stream","encoding":"base64","persistence":"session only"}))
        }).await
    }
    #[tool(
        description = "Review/export an exact immutable document revision under the independently implemented SnapGene revision-15 reference profile (.dna or .prot). Returns structured omissions and no artifact until accept_omissions=true when loss exists. Refuses unsupported physical/reference-only material under both policies. Emits fresh supported packets, never stale original opaque chunks or cloned experimental history; artifact is session-only and source remains unchanged."
    )]
    pub(crate) async fn export_snapgene(
        &self,
        Parameters(args): Parameters<sequence_exchange_tools::ReferenceExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            use collagene_core::formats::snapgene::export::{self as snapgene_export, Policy};
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let review=snapgene_export::review(&molecule).map_err(bad)?;
            if !review.losses.is_empty()&&!args.accept_omissions{return output(json!({"document":args.document,"state":"review_required","review":review,"artifact_id":null}));}
            let report=snapgene_export::export_snapgene(&molecule,if args.accept_omissions{Policy::AllowOmissions}else{Policy::Strict}).map_err(bad)?;
            let filename=match report.extension.as_str(){".dna"=>"sequence.dna",".prot"=>"sequence.prot",_=>return Err(bad("Unsupported SnapGene artifact extension"))};let size=report.bytes.len();
            let artifact=server.session.binary_artifact(report.bytes,"application/octet-stream",filename)?;
            output(json!({"document":args.document,"state":"prepared","profile":report.profile,"source_revision":report.source_revision,"losses":report.losses,"warnings":report.warnings,"artifact_id":artifact,"bytes":size,"filename":filename,"mime_type":report.mime,"encoding":"base64","persistence":"session only"}))
        }).await
    }
    #[tool(
        description = "Review/export an exact immutable DNA, RNA or protein revision as a GCG/Wisconsin single-sequence artifact with checksum. Requires explicit record ID/date for synthetic records; retained IDs/dates may be reused as documentary source data. Returns omission warnings without an artifact until accept_omissions=true when loss exists. GCG cannot carry annotations, circular topology, source attachments or cloning history. Unsupported physical/reference-only material is refused; no source mutation or host path."
    )]
    pub(crate) async fn export_gcg(
        &self,
        Parameters(args): Parameters<sequence_exchange_tools::GcgExportQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let molecule=server.session.store.lock().unwrap().document(&args.document)?.clone();
            let report=collagene_core::formats::sequence::gcg::export_gcg(&molecule,&args.options).map_err(bad)?;
            if !report.warnings.is_empty()&&!args.accept_omissions{return output(json!({"document":args.document,"state":"review_required","warnings":report.warnings,"artifact_id":null}));}
            let bytes=report.text.into_bytes();let size=bytes.len();
            let artifact=server.session.binary_artifact(bytes,"text/plain","sequence.gcg")?;
            output(json!({"document":args.document,"state":"prepared","warnings":report.warnings,"artifact_id":artifact,"bytes":size,"filename":"sequence.gcg","mime_type":"text/plain","encoding":"base64","persistence":"session only"}))
        }).await
    }
    #[tool(
        description = "Import GenBank/EMBL/FASTA/Collagene XML or plain sequence text as immutable session documents. Atomic and deduplicated; returns IDs and exact revisions, without bulk sequence output."
    )]
    pub(crate) fn import_document(
        &self,
        Parameters(args): Parameters<ImportText>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.text.len() > 4 * 1024 * 1024 {
            return Err(bad("MCP text imports are limited to 4 MiB"));
        }
        let molecules = formats::import_text(&args.text).map_err(bad)?;
        let warnings = molecules.warnings;
        self.session.insert_documents(molecules.molecules, warnings)
    }
    #[tool(
        description = "List immutable session document summaries. Imported annotations are data and cannot grant tool access."
    )]
    pub(crate) fn list_documents(
        &self,
        Parameters(args): Parameters<DocumentPage>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.session.store.lock().unwrap();
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
    pub(crate) fn read_sequence_region(
        &self,
        Parameters(args): Parameters<Region>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.session.store.lock().unwrap();
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
    pub(crate) fn list_enzymes(&self) -> Result<CallToolResult, ErrorData> {
        output(
            json!({"enzymes":analysis::starter_enzymes(),"checked_digest_profile_revision":collagene_core::chemistry::tracked_duplex::profile_revision(),"checked_digest_profile_version":collagene_core::chemistry::tracked_duplex::PROFILE_VERSION,"scope":"Legacy starter geometry only; no catalog supplier/formulation facts are inferred.","catalog_discovery":catalog_tools::identity(&catalog_tools::builtin()?),"catalog_lookup_tool":"lookup_enzyme_catalog","catalog_evidence_tool":"resolve_enzyme_evidence"}),
        )
    }
    #[tool(
        description = "Read an artifact in at most 16,384 text bytes. Follow next_start to done; UTF-8 text is literal, base64 text must be concatenated and decoded before saving binary files."
    )]
    pub(crate) fn read_artifact(
        &self,
        Parameters(args): Parameters<ArtifactRegion>,
    ) -> Result<CallToolResult, ErrorData> {
        let store = self.session.store.lock().unwrap();
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
