//! Documentary image tools. Registered handlers use the native immutable owner;
//! neither uploaded pixels nor an archived preview can confer scientific authority.
use crate::{
    server::{
        CollageneServer,
        response::{bad, bounded_json_limit, output},
    },
    session::references::WorkspaceRef,
};
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_core::{
    documents::image_evidence::{self, Binding, ImageEvidence, Prepared},
    material::lineage::workspace_fingerprint,
};
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

const MAX_INLINE_BYTES: usize = 8 * 1024 * 1024;
const AUTHORITY: &str = "Documentary image attachment only. No sequence, feature, material, history event or verification authority is created by pixels or an image preview.";

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AttachImageEvidence {
    pub workspace: WorkspaceRef,
    /// Exact current document identity and head revision in this workspace.
    pub document: Binding,
    /// Fresh fingerprint returned by list_image_evidence for this exact workspace.
    pub expected_source_fingerprint: String,
    /// Display label only; it is never a filesystem path or artifact filename.
    pub name: String,
    /// Exact original JPEG/TIFF/PNG bytes; at most 8 MiB decoded.
    pub data_base64: String,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ListImageEvidence {
    pub workspace: WorkspaceRef,
    /// Exact retained selection. All images for this document are discoverable,
    /// including records associated with its earlier retained revisions.
    pub document: Binding,
    #[serde(default)]
    pub start: usize,
    #[serde(default = "default_limit")]
    pub limit: usize,
}
fn default_limit() -> usize {
    20
}
#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
pub enum Representation {
    Original,
    Preview,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ExportImageEvidence {
    pub workspace: WorkspaceRef,
    /// Must equal this image record's original associated document/revision.
    pub document: Binding,
    pub manifest_hash: String,
    pub representation: Representation,
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn bounded_selection(workspace: &WorkspaceRef, document: &Binding) -> Result<(), ErrorData> {
    if !digest(&workspace.id)
        || !digest(&workspace.revision)
        || document.document_id.is_empty()
        || document.document_id.len() > 128
        || !digest(&document.revision)
    {
        return Err(bad(
            "bounded exact workspace and document references are required",
        ));
    }
    Ok(())
}
fn original_filename(format: image_evidence::Format) -> &'static str {
    match format {
        image_evidence::Format::Jpeg => "documentary-image.jpg",
        image_evidence::Format::Tiff => "documentary-image.tif",
        image_evidence::Format::Png => "documentary-image.png",
    }
}

#[tool_router(router = image_evidence_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Attach exact JPEG/TIFF/PNG original bytes (<=8 MiB decoded base64) to an explicitly selected current document/revision in an exact immutable session workspace. Requires its fresh source fingerprint. The private native compute owner decodes bounded pixels and prepares a documentary PNG; original bytes remain exact. Atomically returns a NEW workspace and image-specific JSON receipt. Original workspace, sequences, features, history and material stay unchanged; no laboratory, gel-comparison or verification claim, database access or filesystem path."
    )]
    pub(crate) async fn attach_image_evidence(
        &self,
        Parameters(args): Parameters<AttachImageEvidence>,
    ) -> Result<CallToolResult, ErrorData> {
        bounded_selection(&args.workspace, &args.document)?;
        if !digest(&args.expected_source_fingerprint)
            || args.name.is_empty()
            || args.name.len() > 256
            || args.data_base64.len() > MAX_INLINE_BYTES.div_ceil(3) * 4
        {
            return Err(bad(
                "image fingerprint/name/base64 input exceeds its bounded profile",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let frame = server.session.restore_workspace(&args.workspace)?;
            args.document.validate_selected(&frame).map_err(bad)?;
            let source_fingerprint = workspace_fingerprint(&frame).map_err(bad)?;
            if source_fingerprint != args.expected_source_fingerprint {
                return Err(bad(
                    "source workspace fingerprint changed; no image was attached",
                ));
            }
            let bytes = STANDARD.decode(args.data_base64).map_err(bad)?;
            if bytes.is_empty() || bytes.len() > MAX_INLINE_BYTES {
                return Err(bad("image input exceeds 8 MiB decoded or is empty"));
            }
            let prepared = Prepared::new(&args.name, args.document, bytes).map_err(bad)?;
            let captured = prepared.apply_to(&frame).map_err(bad)?;
            if captured.source_workspace_fingerprint != source_fingerprint {
                return Err(bad("image preparation source fingerprint differs"));
            }
            let receipt = json!({
                "source_workspace":args.workspace,"source_workspace_fingerprint":source_fingerprint,
                "image":captured.summary,"state":"documentary","verified":false,
                "original_workspace_unchanged":true,"authority":AUTHORITY
            });
            let text = bounded_json_limit(&receipt, 64 * 1024)?;
            server.session.commit_workspace_outputs_named(
                captured.workspace,
                Some((text, receipt)),
                None,
                "image-evidence.json",
                AUTHORITY,
            )
        })
        .await
    }

    #[tool(
        description = "List documentary images for one explicitly selected retained document/revision in an exact immutable session workspace. Earlier image association revisions remain explicit and discoverable after sequence edits. Returns a source fingerprint for attach CAS and pages of at most25 metadata records; never returns bulk source/preview pixels or restores verification authority."
    )]
    pub(crate) async fn list_image_evidence(
        &self,
        Parameters(args): Parameters<ListImageEvidence>,
    ) -> Result<CallToolResult, ErrorData> {
        bounded_selection(&args.workspace, &args.document)?;
        if args.limit == 0 || args.limit > 25 || args.start > image_evidence::MAX_IMAGES {
            return Err(bad(
                "image metadata page limit is1..25; start must fit retained records",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let frame = server.session.restore_workspace(&args.workspace)?;
            args.document.validate_retained(&frame).map_err(bad)?;
            let images: Vec<_> = ImageEvidence::new(&frame)
                .map_err(bad)?
                .list()
                .map_err(bad)?
                .into_iter()
                .filter(|summary| {
                    summary.record.association.document_id == args.document.document_id
                })
                .collect();
            if args.start > images.len() {
                return Err(bad(
                    "image metadata start exceeds selected document's retained records",
                ));
            }
            let end = args.start.saturating_add(args.limit).min(images.len());
            output(json!({
                "workspace":args.workspace,"selected_document":args.document,
                "source_workspace_fingerprint":workspace_fingerprint(&frame).map_err(bad)?,
                "images":&images[args.start..end],"total":images.len(),"start":args.start,
                "next_start":if end < images.len() { Some(end) } else { None },
                "state":"documentary","verified":false,"authority":AUTHORITY
            }))
        })
        .await
    }

    #[tool(
        description = "Export one documentary image from an exact immutable session workspace and its exact original document/revision association. Original returns byte-exact retained JPEG/TIFF/PNG without claiming a fresh decode; Preview freshly reopens original bytes with the current bounded decoder and compares the stored descriptor/preview before exporting PNG. Returns a paginated base64 artifact reference, never a filesystem write or scientific/gel-comparison/verification authority."
    )]
    pub(crate) async fn export_image_evidence(
        &self,
        Parameters(args): Parameters<ExportImageEvidence>,
    ) -> Result<CallToolResult, ErrorData> {
        bounded_selection(&args.workspace, &args.document)?;
        if !digest(&args.manifest_hash) {
            return Err(bad("exact image manifest hash is required"));
        }
        let server = self.clone();
        self.compute(move || {
            let frame = server.session.restore_workspace(&args.workspace)?;
            args.document.validate_retained(&frame).map_err(bad)?;
            let owner = ImageEvidence::new(&frame).map_err(bad)?;
            let summary = owner
                .list()
                .map_err(bad)?
                .into_iter()
                .find(|s| s.manifest_hash == args.manifest_hash)
                .ok_or_else(|| bad("unknown image manifest in this exact workspace"))?;
            if summary.record.association != args.document {
                return Err(bad(
                    "image association differs from explicitly selected document/revision",
                ));
            }
            let (bytes, media_type, filename, representation, redecoded) = match args.representation
            {
                Representation::Original => {
                    let original = owner.original(&args.manifest_hash).map_err(bad)?;
                    (
                        original.bytes.to_vec(),
                        original.media_type,
                        original_filename(summary.record.format),
                        "Original",
                        false,
                    )
                }
                Representation::Preview => {
                    let preview = owner
                        .prepared(&args.manifest_hash)
                        .map_err(bad)?
                        .reopen()
                        .map_err(bad)?
                        .into_preview();
                    preview.validate_binding(&summary).map_err(bad)?;
                    (
                        preview.png,
                        "image/png",
                        "documentary-image.png",
                        "Preview",
                        true,
                    )
                }
            };
            let length = bytes.len();
            let content_hash = blake3::hash(&bytes).to_hex().to_string();
            let artifact_id = server
                .session
                .binary_artifact(bytes, media_type, filename)?;
            output(json!({
                "workspace":args.workspace,"image":summary,"representation":representation,
                "artifact_id":artifact_id,"bytes":length,"content_hash":content_hash,
                "mime_type":media_type,"filename":filename,"encoding":"base64",
                "freshly_redecoded":redecoded,"calibrated":false,
                "state":"documentary","verified":false,"authority":AUTHORITY
            }))
        })
        .await
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::session::{Artifact, ArtifactData, MAX_ITEMS, MAX_STORE_BYTES};
    use collagene_core::documents::workspace::Workspace;
    use collagene_model::{Molecule, Topology};
    const PNG: &[u8] = include_bytes!("../../../tests/data/image-evidence/known_2x3.png");

    #[tokio::test]
    async fn failed_image_capacity_cannot_publish_workspace_or_receipt() {
        for exhaust_items in [true, false] {
            let server = CollageneServer::default();
            let frame = Workspace::from_documents(
                "Image capacity",
                vec![Molecule::dna("Reference", "ACGT", Topology::Linear).unwrap()],
                0,
            )
            .unwrap();
            let summary = server
                .session
                .commit_workspace(frame.clone())
                .unwrap()
                .structured_content
                .unwrap();
            let id = summary["workspace"]["id"].as_str().unwrap().to_owned();
            let reference = WorkspaceRef {
                id: id.clone(),
                revision: id,
            };
            let before = {
                let mut store = server.session.store.lock().unwrap();
                if exhaust_items {
                    for index in 0..MAX_ITEMS {
                        store.artifacts.insert(
                            format!("occupied-{index}"),
                            Artifact {
                                data: ArtifactData::Text(String::new()),
                                mime_type: "text/plain",
                                encoding: "utf-8",
                                filename: "occupied.txt",
                            },
                        );
                    }
                } else {
                    store.used = MAX_STORE_BYTES - 1;
                }
                (
                    store.workspaces.len(),
                    store.documents.len(),
                    store.artifacts.len(),
                    store.used,
                )
            };
            assert!(
                server
                    .attach_image_evidence(Parameters(AttachImageEvidence {
                        workspace: reference.clone(),
                        document: Binding {
                            document_id: frame.documents[0].id.clone(),
                            revision: frame.documents[0].head.clone()
                        },
                        expected_source_fingerprint: workspace_fingerprint(&frame).unwrap(),
                        name: "Documentary image".into(),
                        data_base64: STANDARD.encode(PNG),
                    }))
                    .await
                    .is_err()
            );
            let after = {
                let store = server.session.store.lock().unwrap();
                (
                    store.workspaces.len(),
                    store.documents.len(),
                    store.artifacts.len(),
                    store.used,
                )
            };
            assert_eq!(after, before);
            let restored = server.session.restore_workspace(&reference).unwrap();
            assert!(
                ImageEvidence::new(&restored)
                    .unwrap()
                    .list()
                    .unwrap()
                    .is_empty()
            );
            assert_eq!(
                serde_json::to_value(restored).unwrap(),
                serde_json::to_value(frame).unwrap()
            );
        }
    }

    #[tokio::test]
    async fn oversized_inline_image_is_rejected_before_compute_or_decode() {
        let server = CollageneServer::default();
        server.close_compute();
        let result = server
            .attach_image_evidence(Parameters(AttachImageEvidence {
                workspace: WorkspaceRef {
                    id: "0".repeat(64),
                    revision: "0".repeat(64),
                },
                document: Binding {
                    document_id: "document-1".into(),
                    revision: "0".repeat(64),
                },
                expected_source_fingerprint: "0".repeat(64),
                name: "Documentary".into(),
                data_base64: "A".repeat(MAX_INLINE_BYTES.div_ceil(3) * 4 + 4),
            }))
            .await;
        assert!(result.unwrap_err().message.contains("bounded profile"));
        let store = server.session.store.lock().unwrap();
        assert!(store.workspaces.is_empty());
        assert!(store.artifacts.is_empty());
        assert_eq!(store.used, 0);
    }
}
