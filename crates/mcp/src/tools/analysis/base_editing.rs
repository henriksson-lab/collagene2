//! Conditional finite base-editing states on immutable document revisions.
use crate::{
    server::{
        CollageneServer, EmptyRequest,
        response::{bad, bounded_json_limit, output},
    },
    session::references::DocumentRef,
};
use collagene_core::sequence::crispr::base_editing;
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

fn default_states() -> usize {
    25
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BaseEditingQuery {
    pub document: DocumentRef,
    pub request: base_editing::Request,
    /// Page in retained conditional states, never in all biological outcomes.
    #[serde(default)]
    pub state_start: usize,
    #[serde(default = "default_states")]
    pub max_states: usize,
}

#[tool_router(router = base_editing_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Discover exact installed base-editor window profiles, editor/guide geometry revisions, conversion direction and primary citations. Current builtin is the published typical BE1 C4-8 window in PAM-distal guide orientation. This is a declared finite reference-substitution model, not empirical efficiency, repair chemistry or an absolute biological activity boundary. Custom declarations retain distinct identities and do not inherit builtin authority. Read-only discovery."
    )]
    pub(crate) async fn list_base_editor_profiles(
        &self,
        Parameters(_args): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        self.compute(move || {
            let profiles = base_editing::catalogue().map_err(bad)?;
            let catalogue_revision = base_editing::catalogue_revision().map_err(bad)?;
            bounded_json_limit(&profiles, 128 * 1024)?;
            output(json!({
                "engine_version": base_editing::ENGINE_VERSION,
                "catalogue_revision": catalogue_revision,
                "profile_count": profiles.len(),
                "profiles": profiles,
                "bounds": {
                    "request_bytes": base_editing::MAX_REQUEST_BYTES,
                    "report_bytes": base_editing::MAX_REPORT_BYTES,
                    "window_width": 12,
                    "conditional_states": 4096,
                    "display_page": 100
                },
                "coordinates": "One-based biological guide positions, counted from the PAM-distal end for the supported three-prime-PAM geometry; editable sites retain original zero-based source positions, both oriented and source-forward substitutions.",
                "statement": base_editing::STATEMENT
            }))
        }).await
    }

    #[tool(
        description = "Assess declared base-editor substitutions on one exact immutable DNA document using the native source/guide-bound kernel. Supply current guide options/candidate, exact editor profile, explicit intended guide positions, finite state cap and required nullable CDS assumption. Returns actual intended/bystander sites and paged retained conditional reference states, including the unedited state, with explicit finite-enumeration completeness and a complete bounded JSON artifact. Optional ordinary codon consequences use the independently bound CDS/genetic-code assumption. Finite completeness is not biological outcome completeness; no calibrated efficiency, probability, physical editing/repair, overall biological pass, candidate/history/material writes or sequence edits."
    )]
    pub(crate) async fn assess_base_editing(
        &self,
        Parameters(args): Parameters<BaseEditingQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.max_states) {
            return Err(bad("state page size must be 1..100"));
        }
        let server = self.clone();
        self.compute(move || {
            let molecule = {
                let store = server.session.store.lock().unwrap();
                let molecule = store.document(&args.document)?;
                // Native validation admits the exact source/guide and resource
                // bounds before cloning an immutable document snapshot.
                args.request.validate(molecule).map_err(bad)?;
                molecule.clone()
            };
            let context = base_editing::revision(&molecule, &args.request).map_err(bad)?;
            let report = base_editing::assess(&molecule, &args.request).map_err(bad)?;
            report.validate_binding(&context).map_err(bad)?;
            if args.state_start > report.states.len() {
                return Err(bad("state page offset exceeds retained conditional states"));
            }
            let end = args
                .state_start
                .saturating_add(args.max_states)
                .min(report.states.len());
            let artifact_id = server.session.text_artifact(
                bounded_json_limit(&report, base_editing::MAX_REPORT_BYTES)?,
                "application/json",
                "base-editing-reference-states.json",
            )?;
            output(json!({
                "document": args.document,
                "engine_version": report.engine_version,
                "context_revision": report.context_revision,
                "source_revision": report.source_revision,
                "source_length": report.source_length,
                "source_topology": report.source_topology,
                "request": report.request,
                "editor_profile_revision": report.request.profile.revision,
                "guide_profile_revision": report.request.profile.guide_profile_revision,
                "report_revision": report.report_revision,
                "sites": report.sites,
                "intent_declared": !report.request.intended_positions.is_empty(),
                "finite_state_count": report.finite_state_count,
                "retained_state_count": report.states.len(),
                "finite_enumeration_complete": report.finite_enumeration_complete,
                "biological_outcomes_exhaustive": report.biological_outcomes_exhaustive,
                "state_start": args.state_start,
                "next_start": end,
                "done": end == report.states.len(),
                "display_truncated": args.state_start != 0 || end != report.states.len(),
                "states": &report.states[args.state_start..end],
                "cds": report.cds,
                "findings": report.findings,
                "statement": report.statement,
                "artifact_id": artifact_id,
                "persistence": "session only"
            }))
        })
        .await
    }
}
