//! Declared research-vector graph checks over exact immutable document handles.
use crate::{
    server::{
        CollageneServer, EmptyRequest,
        response::{bad, bounded_json_limit, output},
    },
    session::references::DocumentRef,
};
use collagene_core::construct_qc::vector_class;
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

fn default_findings() -> usize {
    25
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VectorClassQuery {
    pub document: DocumentRef,
    pub request: vector_class::Request,
    #[serde(default)]
    pub finding_start: usize,
    #[serde(default = "default_findings")]
    pub max_findings: usize,
}

#[tool_router(router = vector_class_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Discover installed declared research-vector graph profiles with exact identities, multiplicity/order/orientation data, boundary-role conventions, descriptive study context, primary sources and limitations. Current scopes are conventional rAAV transfer-cassette and SIN-labelled lentiviral transfer-plasmid annotation graphs. No biological reference DNA, inferred vector class, packaging/helper-system design, replication-competence, tropism, yield or safety verdict. Read-only discovery."
    )]
    pub(crate) async fn list_vector_class_profiles(
        &self,
        Parameters(_args): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        self.compute(move || {
            let profiles = vector_class::catalogue().map_err(bad)?;
            bounded_json_limit(&profiles, 256 * 1024)?;
            output(json!({
                "engine_version": vector_class::ENGINE_VERSION,
                "catalogue_revision": vector_class::catalogue_revision().map_err(bad)?,
                "profile_count": profiles.len(),
                "profiles": profiles,
                "bounds": {
                    "request_bytes": vector_class::MAX_REQUEST_BYTES,
                    "report_bytes": vector_class::MAX_REPORT_BYTES,
                    "bindings": 64,
                    "spacing_constraints": 32,
                    "display_findings": 100
                },
                "coordinates": "Original zero-based half-open source intervals. The user declares exact start/end binding IDs and graph axis; boundary size includes the two selected elements, while interior size excludes them. Neither is replaced by total plasmid length.",
                "statement": vector_class::STATEMENT
            }))
        }).await
    }

    #[tool(
        description = "Assess an explicit source/profile-bound research-vector graph against one immutable DNA document. Supply current annotation assumptions, pinned user reference sequences with provenance, role bindings, required nullable boundary selection, explicit spacing and optional bounded literal-copy/forbidden-literal requirements. Reuses native reference identity and literal search. Returns paged typed findings, actual identity verdicts separately from dependent role premises, original-coordinate boundary measurements and descriptive study context plus complete bounded JSON artifact. Missing/unverified roles stay indeterminate; a literal mismatch stays failed without implying biological absence. No overall biological pass, packaging/titre/propagation prediction, source edits, history/material writes or helper-system design."
    )]
    pub(crate) async fn assess_vector_class(
        &self,
        Parameters(args): Parameters<VectorClassQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=100).contains(&args.max_findings) {
            return Err(bad("finding page size must be 1..100"));
        }
        let server = self.clone();
        self.compute(move || {
            let molecule = {
                let store = server.session.store.lock().unwrap();
                let molecule = store.document(&args.document)?;
                args.request.validate(molecule).map_err(bad)?;
                molecule.clone()
            };
            let context = vector_class::revision(&molecule, &args.request).map_err(bad)?;
            let report = vector_class::assess(&molecule, &args.request).map_err(bad)?;
            report.validate_binding(&context).map_err(bad)?;
            if args.finding_start > report.findings.len() {
                return Err(bad("finding offset exceeds retained vector-class findings"));
            }
            let end = args
                .finding_start
                .saturating_add(args.max_findings)
                .min(report.findings.len());
            let mut summary = json!({
                "document": args.document,
                "engine_version": report.engine_version,
                "context_revision": report.context_revision,
                "source_revision": report.source_revision,
                "source_length": report.source_length,
                "source_feature_count": report.source_feature_count,
                "source_topology": report.source_topology,
                "annotations_revision": report.annotations_revision,
                "assumptions_revision": report.assumptions_revision,
                "references_revision": report.references_revision,
                "report_revision": report.report_revision,
                "profile": report.profile.reference,
                "scope": report.profile.definition.scope,
                "sources": report.profile.definition.sources,
                "limitations": report.profile.definition.limitations,
                "bindings": report.bindings,
                "boundaries": report.boundaries,
                "numeric_context": report.profile.definition.numeric_context,
                "finding_count": report.findings.len(),
                "finding_start": args.finding_start,
                "next_start": end,
                "done": end == report.findings.len(),
                "display_truncated": args.finding_start != 0 || end != report.findings.len(),
                "findings": &report.findings[args.finding_start..end],
                "literal_parts": report.literal_parts,
                "literal_search_work": report.literal_search_work,
                "literal_searches": report.literal_searches.iter().map(|search| json!({
                    "reference_id": search.reference_id,
                    "reference_revision": search.reference_revision,
                    "complete": search.complete,
                    "work": search.work,
                    "retained_matches": search.matches.len(),
                    "reason": search.reason
                })).collect::<Vec<_>>(),
                "literal_findings": report.literal_findings,
                "statement": report.statement,
                "persistence": "session only"
            });
            // Reserve the compact response budget before committing its immutable
            // artifact. Full user reference sequences remain in that artifact.
            bounded_json_limit(&summary, 1024 * 1024 - 1024)?;
            let artifact_id = server.session.text_artifact(
                bounded_json_limit(&report, vector_class::MAX_REPORT_BYTES)?,
                "application/json",
                "vector-class-qc.json",
            )?;
            summary["artifact_id"] = json!(artifact_id);
            output(summary)
        })
        .await
    }
}
