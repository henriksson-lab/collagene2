//! Immutable-source alternative IVT restriction hypotheses, with no automatic winner.
use crate::{
    server::{
        CollageneServer,
        response::{bad, bounded_json_limit, output},
    },
    session::references::DocumentRef,
};
use collagene_core::cloning::restriction::ivt_linearization as native;
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct IvtLinearizationQuery {
    pub document: DocumentRef,
    pub request: native::Request,
}

#[tool_router(router = ivt_linearization_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Inspect alternative single-enzyme IVT template linearization hypotheses on one immutable DNA document. Reuse list_ivt_profiles and restriction catalogue discovery for exact installed identities. Supply the complete typed source-bound IVT declaration (linearization must be null), exact enzyme selection, original-generation assumptions OR an existing saved MaterialOrigin, nullable exact reaction conditions and explicit eligibility/unknown-history/resource policies. Native Rust replays saved material and reuses positional checked cleavage; it retains each enzyme's completeness, unknown/blocked sites, multiple cuts, all fragment geometries and exact forward/reverse/circular coding/template cut offsets. Unique-cut/runoff/promoter-retention eligibility is conditional on declarations, never an automatic winner, completed cleavage, isolation, transcription/capping yield or experimental success. Returns compact alternatives and a complete bounded JSON artifact; no edits or material/history writes."
    )]
    pub(crate) async fn assess_ivt_linearization(
        &self,
        Parameters(args): Parameters<IvtLinearizationQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            // Validate source/request/capsule bounds before cloning the immutable document.
            // Fresh native replay and the actual checked digest run after releasing the store lock.
            let molecule = {
                let store = server.session.store.lock().unwrap();
                let source = store.document(&args.document)?;
                args.request.validate(source).map_err(bad)?;
                source.clone()
            };
            let expected = native::revision(&molecule, &args.request).map_err(bad)?;
            let report = native::assess(&molecule, &args.request).map_err(bad)?;
            report.validate_binding(&expected).map_err(bad)?;
            let artifact_id = server.session.text_artifact(
                bounded_json_limit(&report, native::MAX_REPORT_BYTES)?,
                "application/json",
                "ivt-linearization.json",
            )?;
            let alternatives = report
                .candidates
                .iter()
                .map(|candidate| {
                    json!({
                        "enzyme_id": candidate.enzyme_id,
                        "state": candidate.state,
                        "issues": candidate.issues,
                        "scan_complete": candidate.cleavage.complete,
                        "products_determined": candidate.cleavage.products_determined,
                        "comparisons": candidate.cleavage.comparisons,
                        "site_count": candidate.cleavage.sites.len(),
                        "sites": candidate.cleavage.sites.iter().map(|site| json!({
                            "recognition_start": site.recognition_start,
                            "reverse": site.reverse,
                            "top_cut": site.top_cut,
                            "bottom_cut": site.bottom_cut,
                            "outcome": site.outcome
                        })).collect::<Vec<_>>(),
                        "limits": candidate.cleavage.limits,
                        "fragment_count": candidate.fragments.len(),
                        "retained_fragment": candidate.retained_fragment,
                        "retained_fragment_geometry": candidate.retained_fragment
                            .and_then(|index| candidate.fragments.get(index)),
                        "mapping": candidate.mapping,
                        "conditions": candidate.conditions.as_ref().map(|conditions| json!({
                            "context_revision": conditions.context_revision,
                            "source_revision": conditions.source_revision,
                            "material_revision": conditions.material_revision,
                            "status": conditions.report.documentary_status,
                            "star_activity": "unassessed",
                            "incubation_duration_compatibility": "unassessed"
                        })),
                        "ivt": candidate.ivt.as_ref().map(|ivt| json!({
                            "context_revision": ivt.context_revision,
                            "report_revision": ivt.report_revision,
                            "linearization": ivt.linearization,
                            "findings": ivt.findings
                        }))
                    })
                })
                .collect::<Vec<_>>();
            let summary = json!({
                "document": args.document,
                "engine_version": report.engine_version,
                "context_revision": report.context_revision,
                "source_revision": report.source_revision,
                "source_topology": report.source_topology,
                "material_revision": report.source_modifications.state_revision,
                "source_history_incomplete": report.source_history_incomplete,
                "source_route": match &report.request.source {
                    native::SourceInput::DeclaredOriginal { .. } => "declared_original",
                    native::SourceInput::Saved { .. } => "fresh_saved_origin_replay"
                },
                "ivt_profile": report.request.ivt.profile,
                "transcript": report.request.ivt.transcript,
                "selection": report.request.selection,
                "reaction_choices": report.request.conditions,
                "eligibility_policy": report.request.policy,
                "limits": report.request.limits,
                "comparisons": report.comparisons,
                "alternative_count": alternatives.len(),
                "alternatives": alternatives,
                "report_revision": report.report_revision,
                "artifact_id": artifact_id,
                "persistence": "session only",
                "authority": report.statement
            });
            bounded_json_limit(&summary, 2 * 1024 * 1024)?;
            output(summary)
        })
        .await
    }
}
