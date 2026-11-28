//! Immutable-source IVT declaration checks; no RNA production or material authority.
use crate::{
    server::{
        CollageneServer, EmptyRequest,
        response::{bad, bounded_json_limit, output},
    },
    session::references::DocumentRef,
};
use collagene_core::rna_chemistry::ivt;
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct IvtTemplateQuery {
    pub document: DocumentRef,
    pub request: ivt::Request,
}

#[tool_router(router = ivt_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Discover exact installed in-vitro-transcription DNA-template check profiles, profile identities, primary source citations, model scope and resource bounds. Current scope is the declared T7 consensus promoter and original CleanCap AG N-7113 initiating-base compatibility. Other polymerases and cap analogues are unsupported. No transcription yield, capping efficiency, end homogeneity, stability or experimental/material authority. Read-only session discovery."
    )]
    pub(crate) async fn list_ivt_profiles(
        &self,
        Parameters(_args): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        self.compute(move || {
            let profiles = ivt::profiles().map_err(bad)?;
            let revision = ivt::catalogue_revision().map_err(bad)?;
            bounded_json_limit(&profiles, 64 * 1024)?;
            output(json!({
                "engine_version": ivt::ENGINE,
                "catalogue_revision": revision,
                "profile_count": profiles.len(),
                "profiles": profiles,
                "bounds": {
                    "template_bases": ivt::MAX_TEMPLATE_BASES,
                    "request_bytes": ivt::MAX_REQUEST_BYTES,
                    "report_bytes": ivt::MAX_REPORT_BYTES,
                    "retained_poly_a_runs": ivt::MAX_RETAINED_RUNS
                },
                "coordinates": "Original zero-based first coding-strand base; reverse coding proceeds towards decreasing indices. Cuts are supplied bond offsets along the coding direction from +1, not evidence cleavage occurred.",
                "authority": "Scoped sequence and supplied-geometry checks only; no measured transcription/capping outcome or material authority."
            }))
        }).await
    }

    #[tool(
        description = "Check an explicit source-revision/profile-bound IVT request against one exact immutable DNA document using embedded Rust. Declare transcript start/length/coding strand, cap choice, optional additional initiator, nullable supplied cut geometry and nullable encoded poly(A) region with provenance. Reports typed promoter/initiator/end/runoff/poly(A) checks, original-coordinate spans, complete counts and bounded segmented-run preview plus complete JSON artifact. Unknown/ambiguous facts remain explicit; no overall biological pass, RNA synthesis, physical cleavage, capping efficiency, yield, stability or edits/history/material writes."
    )]
    pub(crate) async fn assess_ivt_template(
        &self,
        Parameters(args): Parameters<IvtTemplateQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let molecule = {
                let store = server.session.store.lock().unwrap();
                let molecule = store.document(&args.document)?;
                args.request.validate_for(molecule).map_err(bad)?;
                molecule.clone()
            };
            let context = ivt::revision(&molecule, &args.request).map_err(bad)?;
            let report = ivt::assess(&molecule, &args.request).map_err(bad)?;
            report.validate_binding(&context).map_err(bad)?;
            let artifact_id = server.session.text_artifact(
                bounded_json_limit(&report, ivt::MAX_REPORT_BYTES)?,
                "application/json",
                "ivt-template.json",
            )?;
            output(json!({
                "document": args.document,
                "engine_version": report.engine_version,
                "context_revision": report.context_revision,
                "source_revision": report.source_revision,
                "source_length": report.source_length,
                "source_topology": report.source_topology,
                "profile": report.profile.reference,
                "request": report.request,
                "report_revision": report.report_revision,
                "transcript_spans": report.transcript_spans,
                "promoter": report.promoter,
                "initiator_context": report.initiator_context,
                "linearization": report.linearization,
                "poly_a": report.poly_a.as_ref().map(|encoding| json!({
                    "region": encoding.region,
                    "adenosines": encoding.adenosines,
                    "other_concrete_bases": encoding.other_concrete_bases,
                    "ambiguous_bases": encoding.ambiguous_bases,
                    "longest_exact_a_run": encoding.longest_exact_a_run,
                    "terminal": encoding.terminal,
                    "runs_total": encoding.runs_total,
                    "retained_run_count": encoding.runs.len(),
                    "segmentation_complete": encoding.segmentation_complete
                })),
                "findings": report.findings,
                "sources": report.profile.sources,
                "assumptions": report.assumptions,
                "artifact_id": artifact_id,
                "persistence": "session only",
                "authority": "Individual sequence/supplied-geometry checks only; no transcription, capping, isolation or experimental-success verdict."
            }))
        }).await
    }
}
