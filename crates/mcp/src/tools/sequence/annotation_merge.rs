//! Explicit annotation decisions over immutable session sources; no native ancestry is invented.
use crate::{
    server::{
        CollageneServer,
        response::{bad, bounded_json_limit, bounded_json_size, output},
    },
    session::references::{DocumentRef, EvidenceRef, WorkspaceRef},
};
use collagene_core::{material::lineage::EventBinding, sequence::annotations::merge as native};
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::json;

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum SourceQuery {
    Library {
        library: EvidenceRef,
        scan: EvidenceRef,
    },
    Dataset {
        dataset: EvidenceRef,
    },
    Donor {
        workspace: WorkspaceRef,
        event: EventBinding,
        transform: native::Transform,
        indices: Vec<usize>,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Query {
    pub document: DocumentRef,
    pub source: SourceQuery,
    pub selections: Vec<native::Selection>,
    #[serde(default)]
    pub limits: native::Limits,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProvenanceQuery {
    pub document: DocumentRef,
    pub feature_index: usize,
}
impl Query {
    fn resolve(
        &self,
        server: &CollageneServer,
    ) -> Result<(collagene_model::Molecule, native::Request), ErrorData> {
        bounded_json_size(self, native::MAX_REQUEST_BYTES)?;
        if self.selections.len() > native::MAX_SELECTIONS {
            return Err(bad("Too many annotation merge decisions"));
        }
        let source = match &self.source {
            SourceQuery::Library { library, scan } => {
                let store = server.session.store.lock().unwrap();
                let library = store.feature_library(library)?;
                let report = store.feature_scan(scan)?;
                native::Source::from_library(library, report).map_err(bad)?
            }
            SourceQuery::Dataset { dataset } => {
                let store = server.session.store.lock().unwrap();
                let report = store.annotation_dataset(dataset)?;
                native::Source::from_dataset(report).map_err(bad)?
            }
            SourceQuery::Donor {
                workspace,
                event,
                transform,
                indices,
            } => {
                let frame = server.session.restore_workspace(workspace)?;
                native::Source::from_donor(&frame, event, transform.clone(), indices)
                    .map_err(bad)?
            }
        };
        let target = {
            let store = server.session.store.lock().unwrap();
            let molecule = store.document(&self.document)?;
            collagene_model::bounds::native_frame_preflight(
                std::iter::once(molecule),
                native::MAX_REQUEST_BYTES,
            )
            .map_err(bad)?;
            molecule.clone()
        };
        let request = native::Request {
            schema_version: 1,
            target_revision: target.revision(),
            source,
            selections: self.selections.clone(),
            limits: self.limits.clone(),
        };
        request.validate(&target).map_err(bad)?;
        Ok((target, request))
    }
}

#[tool_router(router = annotation_merge_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Preview explicit Append, Skip or AddProvenance annotation decisions from a retained library scan, imported annotation dataset or exact donor event in an immutable workspace. Reconstructs source evidence once; distinguishes overlap, conflicts, orientation and equivalent annotations. Attribution is documentary, not biological verification or native ancestry. Returns a bounded preview and full JSON artifact; changes no document."
    )]
    pub(crate) async fn preview_annotation_merge(
        &self,
        Parameters(args): Parameters<Query>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let (target, request) = args.resolve(&server)?;
            let context = native::revision(&target, &request).map_err(bad)?;
            let preview = native::preview(&target, &request).map_err(bad)?;
            preview.validate_binding(&context).map_err(bad)?;
            preview.validate_request(&target, &request).map_err(bad)?;
            let candidates = preview.candidates.iter().take(20).map(|candidate| json!({
                "proposal_id":candidate.proposal_id,
                "decision":candidate.decision,
                "kind":candidate.feature.kind,
                "label":candidate.feature.label().chars().take(256).collect::<String>(),
                "relations":candidate.relations.iter().take(10).collect::<Vec<_>>(),
                "relations_preview_truncated":candidate.relations.len()>10,
                "relations_complete":candidate.relations_complete,
            })).collect::<Vec<_>>();
            let mut summary = json!({"document":args.document,"context_revision":context,"report_revision":preview.report_revision,"candidate_count":preview.candidates.len(),"candidates":candidates,"candidates_preview_truncated":preview.candidates.len()>20,"complete":preview.complete,"source_complete":preview.source_complete,"comparisons":preview.comparisons,"warnings":preview.warnings,"applied":false,"authority":"Documentary attribution and selected coordinate correspondence; no function, material or experiment is verified."});
            bounded_json_size(&summary, 1024*1024-1024)?;
            let artifact = server.session.text_artifact(bounded_json_limit(&preview,native::MAX_REPORT_BYTES)?,"application/json","annotation-merge-preview.json")?;
            summary["artifact_id"]=json!(artifact);
            output(summary)
        }).await
    }
    #[tool(
        description = "Atomically execute explicit annotation merge decisions in a NEW immutable session document. Revalidates the actual target and source evidence, preserves DNA and existing feature identity, and never overwrites biological qualifiers or clears review. Append creates review-required features; Skip is explicit; AddProvenance targets an exact annotation revision with compatible geometry and meaning. Returns actual decision receipt and full bounded JSON. Donor source ownership is resolved from its selected workspace; no native history or physical cloning authority is created."
    )]
    pub(crate) async fn apply_annotation_merge(
        &self,
        Parameters(args): Parameters<Query>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let (target, request)=args.resolve(&server)?;
            let context=native::revision(&target,&request).map_err(bad)?;
            let applied=native::apply(&target,&request).map_err(bad)?;
            applied.validate_binding(&context).map_err(bad)?;
            applied.validate_request(&target, &request).map_err(bad)?;
            let mut summary=json!({"source":args.document,"context_revision":context,"receipt":applied.receipt,"persistence":"session only","warnings":applied.review.warnings,"native_history_created":false,"physical_cloning_verified":false});
            bounded_json_size(&summary,1024*1024-4096)?;
            let (artifact,document)=server.session.candidate_document_artifact_bounded(&applied,Some(applied.molecule.clone()),"annotation-merge-applied.json",native::MAX_REPORT_BYTES)?;
            summary["artifact_id"]=json!(artifact);
            summary["documents"]=json!([document.ok_or_else(||bad("Missing annotation merge document receipt"))?]);
            output(summary)
        }).await
    }
    #[tool(
        description = "Inspect the single typed accumulated annotation-source contribution ledger on one exact immutable document feature. Returns documentary library, source-file and donor-event claims; validates the current ledger and rejects unsupported versions. Imported claims do not establish authenticated source ownership, annotation inheritance, biological function or experimental validation."
    )]
    pub(crate) async fn inspect_annotation_provenance(
        &self,
        Parameters(args): Parameters<ProvenanceQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let store=server.session.store.lock().unwrap();
            let molecule=store.document(&args.document)?;
            let feature=molecule.features.get(args.feature_index).ok_or_else(||bad("Annotation index is out of range"))?;
            let ledger=native::Ledger::from_feature(feature).map_err(bad)?;
            let result=json!({"document":args.document,"feature_index":args.feature_index,"ledger":ledger,"authority":"Documentary contributions only; no authenticated biological or native ancestry claim."});
            bounded_json_size(&result,1024*1024)?;
            output(result)
        }).await
    }
}
