//! Native MCP graph exchange. Source handles and workspace attachments remain
//! documentary; only explicit projection creates a new native reference.
use crate::{
    ArtifactData, CollageneServer, DocumentRef, EvidenceRef, WorkspaceRef, bad, bounded_json_limit,
    output, summary,
};
use collagene_core::workspace;
use collagene_design_exchange::{self as sbol, Kind, ProjectionPolicy};
use collagene_model::Molecule;
use rmcp::{ErrorData, model::CallToolResult};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::collections::BTreeMap;
const MAX_SOURCE_BYTES: usize = 4 * 1024 * 1024;
const MIME: &str = "application/ld+json";
const FILENAME: &str = "design.sbol.jsonld";
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum Source {
    Imported(EvidenceRef),
    Workspace {
        workspace: WorkspaceRef,
        source_hash: String,
    },
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportQuery {
    pub text: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InspectQuery {
    pub source: Source,
    pub component_start: usize,
    pub max_components: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProjectQuery {
    pub source: Source,
    pub selection: sbol::Selection,
    #[serde(default)]
    pub accept_omissions: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub enum Representation {
    Original,
    Expanded,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GraphExportQuery {
    pub source: Source,
    pub representation: Representation,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NativeExportQuery {
    pub document: DocumentRef,
    pub ids: sbol::ExportIds,
    #[serde(default)]
    pub accept_omissions: bool,
}
fn source_text(bytes: &[u8]) -> Result<String, ErrorData> {
    if bytes.len() > MAX_SOURCE_BYTES {
        return Err(bad("MCP SBOL source exceeds 4 MiB"));
    }
    String::from_utf8(bytes.to_vec()).map_err(bad)
}
impl CollageneServer {
    fn sbol_source(&self, source: &Source) -> Result<String, ErrorData> {
        let store = self.store.lock().unwrap();
        let (bytes, expected): (&[u8], &str) = match source {
            Source::Imported(reference) => {
                if reference.id != reference.revision || !store.sbol_sources.contains(&reference.id)
                {
                    return Err(bad("unknown, stale or foreign-role SBOL source"));
                }
                let artifact = store
                    .artifacts
                    .get(&reference.id)
                    .ok_or_else(|| bad("missing SBOL source"))?;
                let ArtifactData::Text(text) = &artifact.data else {
                    return Err(bad("SBOL source is not text"));
                };
                if artifact.mime_type != MIME || artifact.filename != FILENAME {
                    return Err(bad("SBOL source artifact role changed"));
                }
                (text.as_bytes(), &reference.revision)
            }
            Source::Workspace {
                workspace,
                source_hash,
            } => {
                if workspace.id != workspace.revision {
                    return Err(bad("stale SBOL workspace"));
                }
                let value = store
                    .workspaces
                    .get(&workspace.id)
                    .ok_or_else(|| bad("unknown SBOL workspace"))?;
                let attachment = value
                    .attachments
                    .get(source_hash)
                    .ok_or_else(|| bad("SBOL source is not attached to this exact workspace"))?;
                if attachment.kind != collagene_library_sharing::SBOL_SOURCE_ATTACHMENT_KIND
                    || attachment.media_type != MIME
                {
                    return Err(bad("workspace attachment is not an original SBOL graph"));
                }
                (&attachment.bytes, source_hash)
            }
        };
        if blake3::hash(bytes).to_hex().as_str() != expected {
            return Err(bad("SBOL source hash mismatch"));
        }
        source_text(bytes)
    }
    fn sbol_inspection_output(
        &self,
        source: Source,
        doc: &sbol::Document,
        start: usize,
        maximum: usize,
        original: Option<String>,
    ) -> Result<CallToolResult, ErrorData> {
        let components = doc
            .objects()
            .values()
            .filter_map(|o| {
                if let Kind::Component(c) = &o.kind {
                    Some((o, c))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        if maximum == 0 || maximum > 50 || start > components.len() {
            return Err(bad(
                "SBOL component page requires max_components 1..50 and a valid start",
            ));
        }
        #[derive(Serialize)]
        struct Inspection<'a> {
            profile: &'static str,
            source: &'a Source,
            source_hash: String,
            objects: &'a BTreeMap<String, sbol::Object>,
            issues: &'a [sbol::Issue],
        }
        let text = bounded_json_limit(
            &Inspection {
                profile: sbol::PROFILE,
                source: &source,
                source_hash: doc.revision(),
                objects: doc.objects(),
                issues: doc.issues(),
            },
            sbol::MAX_BYTES,
        )?;
        let artifact = if let Some(original) = original {
            let (source_id, report_id) = self.paired_text_artifacts([
                (original, MIME, FILENAME),
                (text, "application/json", "sbol-inspection.json"),
            ])?;
            self.store.lock().unwrap().sbol_sources.insert(source_id);
            report_id
        } else {
            self.text_artifact(text, "application/json", "sbol-inspection.json")?
        };
        let rows=components.iter().skip(start).take(maximum).map(|(o,c)|json!({
            "component":o.identity.iri,"name":o.identity.name.as_ref().map(|n|n.chars().take(256).collect::<String>()),
            "types":c.types.iter().take(4).collect::<Vec<_>>(),"types_truncated":c.types.len()>4,
            "sequence_count":c.sequences.len(),"sequence_preview_truncated":c.sequences.len()>4,
            "sequences":c.sequences.iter().take(4).map(|id|{
                let record=doc.objects().get(id).and_then(|o|if let Kind::Sequence(s)=&o.kind{Some(s)}else{None});
                json!({"iri":id,"length":record.and_then(|s|s.elements.as_ref().map(String::len)),"encoding":record.and_then(|s|s.encoding.as_ref())})
            }).collect::<Vec<_>>(),"features":c.features.len(),"constraints":c.constraints.len(),
        })).collect::<Vec<_>>();
        output(
            json!({"source":source,"profile":sbol::PROFILE,"source_hash":doc.revision(),"object_count":doc.objects().len(),"component_count":components.len(),"component_start":start,"components":rows,"components_truncated":start.saturating_add(maximum)<components.len(),"issue_count":doc.issues().len(),"issue_preview":doc.issues().iter().take(20).collect::<Vec<_>>(),"issues_truncated":doc.issues().len()>20,"artifact_id":artifact,"authority":"Retained design data only; no selected reference, cloning material or verification is created."}),
        )
    }
    pub(crate) fn sbol_import(&self, args: ImportQuery) -> Result<CallToolResult, ErrorData> {
        if args.text.len() > MAX_SOURCE_BYTES {
            return Err(bad("MCP SBOL source exceeds 4 MiB"));
        }
        let doc = sbol::import(args.text.as_bytes()).map_err(bad)?;
        let id = doc.revision();
        self.sbol_inspection_output(
            Source::Imported(EvidenceRef {
                id: id.clone(),
                revision: id,
            }),
            &doc,
            0,
            20,
            Some(args.text),
        )
    }
    pub(crate) fn sbol_inspect(&self, args: InspectQuery) -> Result<CallToolResult, ErrorData> {
        let text = self.sbol_source(&args.source)?;
        let doc = sbol::import(text.as_bytes()).map_err(bad)?;
        self.sbol_inspection_output(
            args.source,
            &doc,
            args.component_start,
            args.max_components,
            None,
        )
    }
    pub(crate) fn sbol_project(&self, args: ProjectQuery) -> Result<CallToolResult, ErrorData> {
        let text = self.sbol_source(&args.source)?;
        let doc = sbol::import(text.as_bytes()).map_err(bad)?;
        let projection = doc
            .project(&args.selection, ProjectionPolicy::AllowOmissions)
            .map_err(bad)?;
        let document = summary(&projection.molecule.revision(), &projection.molecule);
        let receipt = json!({"state":"created","profile":sbol::PROFILE,"source":args.source,"source_hash":doc.revision(),"selection":args.selection,"document":document,"issue_count":projection.issues.len(),"issue_preview":projection.issues.iter().take(32).collect::<Vec<_>>(),"issues_truncated":projection.issues.len()>32,"warnings":projection.warnings});
        #[derive(Serialize)]
        struct Evidence<'a> {
            profile: &'static str,
            source: &'a Source,
            source_hash: String,
            selection: &'a sbol::Selection,
            molecule: &'a Molecule,
            issues: &'a [sbol::Issue],
            warnings: &'a [String],
        }
        let report = bounded_json_limit(
            &Evidence {
                profile: sbol::PROFILE,
                source: &args.source,
                source_hash: doc.revision(),
                selection: &args.selection,
                molecule: &projection.molecule,
                issues: &projection.issues,
                warnings: &projection.warnings,
            },
            sbol::MAX_BYTES,
        )?;
        if !args.accept_omissions && !projection.issues.is_empty() {
            let artifact =
                self.text_artifact(report, "application/json", "sbol-projection.json")?;
            return output(
                json!({"state":"review_required","source":args.source,"selection":args.selection,"document":null,"issue_count":projection.issues.len(),"issue_preview":projection.issues.iter().take(32).collect::<Vec<_>>(),"issues_truncated":projection.issues.len()>32,"artifact_id":artifact,"authority":"Read-only proposed reference and full omission review; no document or workspace has been created."}),
            );
        }
        let mut frame =
            workspace::Workspace::from_documents("SBOL reference", vec![projection.molecule], 0)
                .map_err(bad)?;
        frame.attachments.insert(
            doc.revision(),
            workspace::WorkspaceAttachment {
                kind: collagene_library_sharing::SBOL_SOURCE_ATTACHMENT_KIND.into(),
                name: FILENAME.into(),
                media_type: MIME.into(),
                bytes: text.into_bytes(),
            },
        );
        frame.validate().map_err(bad)?;
        self.commit_workspace_outputs_named(frame,Some((report,receipt)),None,"sbol-projection.json","Documentary SBOL reference with its original graph only; no imported material, historical run or verification authority.")
    }
    pub(crate) fn sbol_export_graph(
        &self,
        args: GraphExportQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let original = self.sbol_source(&args.source)?;
        let doc = sbol::import(original.as_bytes()).map_err(bad)?;
        let text = match args.representation {
            Representation::Original => original,
            Representation::Expanded => {
                String::from_utf8(doc.expanded_bytes().map_err(bad)?).map_err(bad)?
            }
        };
        let length = text.len();
        let hash = blake3::hash(text.as_bytes()).to_hex().to_string();
        let artifact = self.text_artifact(text, MIME, FILENAME)?;
        output(
            json!({"source":args.source,"source_hash":doc.revision(),"representation":args.representation,"artifact_id":artifact,"bytes":length,"blake3":hash,"filename":FILENAME,"mime_type":MIME,"encoding":"utf-8","authority":"Full local documentary graph exchange only, not a lab publication or material/verification grant."}),
        )
    }
    pub(crate) fn sbol_export_native(
        &self,
        args: NativeExportQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let molecule = self.store.lock().unwrap().document(&args.document)?.clone();
        collagene_core::format_fidelity::reference_export_guard(&molecule).map_err(bad)?;
        let exported =
            sbol::export_molecule(&molecule, &args.ids, ProjectionPolicy::AllowOmissions)
                .map_err(bad)?;
        if !args.accept_omissions && !exported.issues.is_empty() {
            return output(
                json!({"state":"review_required","document":args.document,"ids":args.ids,"issues":exported.issues,"artifact_id":null}),
            );
        }
        let length = exported.bytes.len();
        let hash = blake3::hash(&exported.bytes).to_hex().to_string();
        let artifact = self.text_artifact(
            String::from_utf8(exported.bytes).map_err(bad)?,
            MIME,
            FILENAME,
        )?;
        output(
            json!({"state":"prepared","profile":sbol::PROFILE,"document":args.document,"ids":args.ids,"issues":exported.issues,"artifact_id":artifact,"bytes":length,"blake3":hash,"filename":FILENAME,"mime_type":MIME,"encoding":"utf-8","authority":"Reviewed sequence-reference design exchange; no physical or verification authority. Source is unchanged."}),
        )
    }
}
