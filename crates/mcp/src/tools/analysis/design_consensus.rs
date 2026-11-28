//! Explicit design/supplied-consensus roles on retained native workspace events.
use crate::{
    server::{
        CollageneServer,
        response::{bad, bounded_json_limit, bounded_json_size, output},
    },
    session::references::WorkspaceRef,
};
use collagene_core::alignment::design_consensus as native;
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

fn page_size() -> usize {
    4
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DesignConsensusQuery {
    pub workspace: WorkspaceRef,
    pub request: native::Request,
    #[serde(default)]
    pub cds_start: usize,
    #[serde(default = "page_size")]
    pub max_cds: usize,
}

#[tool_router(router = design_consensus_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Compare explicitly named design and supplied-consensus retained events in one exact immutable workspace.
        Supply whole-construct alignment options, explicit consensus-axis assumptions and optional design-CDS annotation/genetic-code declarations.
        Computes one native pairwise report and conditional projected-CDS consequences on that selected alignment path.
        Returns exact event/settings/search-evidence pins, paged CDS completeness/findings and a full bounded JSON artifact.
        Supplied consensus is a reference sequence, not reads, quality, a measured sample or biological verification.
        Alternative alignments, ambiguous bases, incomplete search/mapping and unsupported CDS remain explicit.
        No edits, history writes, function/expression prediction or material authority."
    )]
    pub(crate) async fn assess_design_consensus(
        &self,
        Parameters(args): Parameters<DesignConsensusQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=4).contains(&args.max_cds) {
            return Err(bad("CDS page size must be 1..4"));
        }
        bounded_json_size(&args.request, native::MAX_REQUEST_BYTES)?;
        let server = self.clone();
        self.compute(move || {
            let workspace = server.session.restore_workspace(&args.workspace)?;
            args.request.validate(&workspace).map_err(bad)?;
            let context = native::revision(&workspace, &args.request).map_err(bad)?;
            let report = native::assess(&workspace, &args.request).map_err(bad)?;
            report.validate_binding(&context).map_err(bad)?;
            report
                .validate_request(&workspace, &args.request)
                .map_err(bad)?;
            if args.cds_start > report.evidence.cds.len() {
                return Err(bad("CDS offset exceeds the retained result count"));
            }
            let evidence = &report.evidence;
            let end = args
                .cds_start
                .saturating_add(args.max_cds)
                .min(evidence.cds.len());
            let cds = report.evidence.cds[args.cds_start..end]
                .iter()
                .map(|row| {
                    json!({
                        "id": row.id,
                        "design_cds": row.design_cds,
                        "code": row.code,
                        "mapping_complete": row.mapping_complete,
                        "register_complete": row.register_complete,
                        "consequences_complete": row.consequences_complete,
                        "register_segment_count": row.register_segments.len(),
                        "consequence_count": row.consequences.len(),
                        "finding_count": row.findings.len(),
                        "findings": row.findings.iter().take(8).collect::<Vec<_>>(),
                        "findings_preview_truncated": row.findings.len() > 8,
                        "design_translation_available": row.design_translation.is_some(),
                        "projected_translation_available": row.projected_translation.is_some()
                    })
                })
                .collect::<Vec<_>>();
            let comparison = &report.pairwise.comparison;
            let mut summary = json!({
                "workspace": args.workspace,
                "engine_version": evidence.engine_version,
                "context_revision": evidence.context_revision,
                "workspace_fingerprint": evidence.workspace_fingerprint,
                "request_revision": evidence.request_revision,
                "evidence_revision": evidence.evidence_revision,
                "design": evidence.request.design,
                "consensus": evidence.request.consensus,
                "consensus_axis": evidence.request.consensus_axis,
                "comparison_options": evidence.request.comparison,
                "comparison_revision": report.pairwise.context_revision,
                "source_provenance": evidence.provenance,
                "reference_origin": comparison.reference_origin,
                "query_reversed": comparison.query_reversed,
                "evaluated_candidates": comparison.evaluated_candidates,
                "exhaustive_reference_origins": comparison.exhaustive_reference_origins,
                "exhaustive_query_orientations": comparison.exhaustive_query_orientations,
                "optimality_proven": comparison.optimality_proven,
                "paired_input_coverage": {
                    "design_fraction": comparison.alignment.stats.reference_coverage,
                    "consensus_fraction": comparison.alignment.stats.query_coverage,
                    "scope": "Paired original-input fractions, not read coverage or confidence"
                },
                "evaluated_cells": comparison.evaluated_cells,
                "candidate_evidence": comparison.candidate_evidence,
                "alignment_search": comparison.alignment.search,
                "difference_count": comparison.alignment.differences.len(),
                "cds_count": report.evidence.cds.len(),
                "cds_start": args.cds_start,
                "next_start": end,
                "done": end == report.evidence.cds.len(),
                "page_truncated": args.cds_start != 0 || end != report.evidence.cds.len(),
                "cds": cds,
                "coding_complete": evidence.coding_complete,
                "selected_path_conditional": evidence.selected_path_conditional,
                "work": evidence.work,
                "statement": evidence.statement,
                "persistence": "session evidence only"
            });
            // Reject a large compact response before publishing an artifact.
            bounded_json_size(&summary, 1024 * 1024 - 1024)?;
            let artifact = server.session.text_artifact(
                bounded_json_limit(&report, native::MAX_REPORT_BYTES)?,
                "application/json",
                "design-consensus.json",
            )?;
            summary["artifact_id"] = json!(artifact);
            output(summary)
        })
        .await
    }
}
