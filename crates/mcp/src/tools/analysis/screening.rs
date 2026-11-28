//! Named screening hypotheses over immutable retained workspace events.
use crate::{
    server::{
        CollageneServer,
        response::{bad, bounded_json_limit, output},
    },
    session::references::{MaterialDerivationOrigin, WorkspaceRef},
};
use collagene_core::{
    chemistry::tracked_duplex::GenerationAssumptions,
    cloning::screening,
    material::{capsule, derivation::SourceOrigin, lineage::EventBinding},
};
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::{Deserialize, Deserializer, Serialize};
use serde_json::json;

fn required_generation<'de, D: Deserializer<'de>>(
    d: D,
) -> Result<Option<GenerationAssumptions>, D::Error> {
    Option::deserialize(d)
}
fn default_page_size() -> usize {
    8
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ScreeningComponent {
    pub id: String,
    pub label: String,
    pub event: EventBinding,
    pub origin: MaterialDerivationOrigin,
    #[serde(deserialize_with = "required_generation")]
    #[schemars(required)]
    pub original_generation: Option<GenerationAssumptions>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ScreeningQuery {
    pub workspace: WorkspaceRef,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, ScreeningComponent, 16>"
    )]
    pub components: Vec<ScreeningComponent>,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, screening::Subject, 16>"
    )]
    pub subjects: Vec<screening::Subject>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, screening::Test, 64>")]
    pub tests: Vec<screening::Test>,
    pub readout: screening::IdealSizeReadout,
    pub limits: screening::Limits,
    #[serde(default)]
    pub test_start: usize,
    #[serde(default = "default_page_size")]
    pub max_tests: usize,
}
impl ScreeningQuery {
    fn native(
        &self,
        workspace: &collagene_core::documents::workspace::Workspace,
    ) -> Result<screening::Request, ErrorData> {
        self.limits.validate().map_err(bad)?;
        if self.components.is_empty()
            || self.components.len() > screening::MAX_COMPONENTS
            || !(2..=screening::MAX_SUBJECTS).contains(&self.subjects.len())
            || self.tests.is_empty()
            || self.tests.len() > screening::MAX_TESTS
        {
            return Err(bad(
                "bounded explicit components, subjects and screening tests required",
            ));
        }
        let mut bases = 0usize;
        let mut borrowed = Vec::with_capacity(self.components.len());
        let mut selections = Vec::new();
        for component in &self.components {
            if component.id.is_empty()
                || component.id.len() > 96
                || component.label.is_empty()
                || component.label.len() > 256
                || component
                    .id
                    .chars()
                    .chain(component.label.chars())
                    .any(char::is_control)
            {
                return Err(bad("bounded component identifiers and labels required"));
            }
            let event = workspace
                .history
                .iter()
                .find(|event| event.id == component.event.event_id)
                .ok_or_else(|| bad("screening component event is unavailable"))?;
            if event.document_id != component.event.document_id
                || event.revision != component.event.revision
            {
                return Err(bad(
                    "screening component document/event/revision binding changed",
                ));
            }
            let source = workspace
                .revisions
                .get(&component.event.revision)
                .ok_or_else(|| bad("screening component revision is unavailable"))?;
            collagene_model::bounds::model_preflight(&[source]).map_err(bad)?;
            source.validate().map_err(bad)?;
            bases = bases
                .checked_add(source.sequence.len())
                .ok_or_else(|| bad("component base bound"))?;
            if bases > 1_000_000 {
                return Err(bad(
                    "screening components exceed one million aggregate bases",
                ));
            }
            match &component.origin {
                MaterialDerivationOrigin::OriginalReference => {}
                MaterialDerivationOrigin::Saved(selection) => {
                    selection.validate().map_err(bad)?;
                    if selection.event() != &component.event
                        || component.original_generation.is_some()
                    {
                        return Err(bad(
                            "saved selection must name the exact component event and cannot declare original generation",
                        ));
                    }
                    selections.push(selection.clone());
                }
            }
            if let Some(generation) = &component.original_generation {
                generation.validate().map_err(bad)?;
            }
            borrowed.push(source);
        }
        collagene_model::bounds::model_preflight(&borrowed).map_err(bad)?;
        // Hydrate closed declarations once, without performing scientific replay
        // here. The native screening owner freshly replays each actual material.
        let mut origins = capsule::documentary_origins(workspace, &selections, &[])
            .map_err(bad)?
            .materials
            .into_iter();
        let mut components = Vec::with_capacity(self.components.len());
        for (component, source) in self.components.iter().zip(borrowed) {
            let origin = match &component.origin {
                MaterialDerivationOrigin::OriginalReference => SourceOrigin::OriginalReference,
                MaterialDerivationOrigin::Saved(_) => SourceOrigin::Saved(
                    origins
                        .next()
                        .ok_or_else(|| bad("native material hydration count differs"))?,
                ),
            };
            origin.validate_shallow_for(&component.event).map_err(bad)?;
            components.push(screening::Component {
                id: component.id.clone(),
                label: component.label.clone(),
                revision: component.event.revision.clone(),
                molecule: source.clone(),
                origin,
                original_generation: component.original_generation.clone(),
            });
        }
        if origins.next().is_some() {
            return Err(bad("unexpected native material hydration"));
        }
        let request = screening::Request {
            components,
            subjects: self.subjects.clone(),
            tests: self.tests.clone(),
            readout: self.readout.clone(),
            limits: self.limits.clone(),
        };
        request.validate().map_err(bad)?;
        Ok(request)
    }
}
#[tool_router(router = screening_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Assess explicitly named intended and decoy products or mixtures using ideal native exact-PCR and checked positional digest size signatures. Every component names an immutable workspace event and an explicit original-reference or retained native-material selection; caller-supplied capsules or synthetic decoy generation are not accepted. Required nullable original_generation declares chemistry only for original references. Returns paged actual cells, possible PCR geometries, signatures, distinction witnesses, ranking and bounded panel search plus a complete immutable JSON artifact. A named empty-vector decoy is an explicit supplied molecule, whereas an empty complete PCR signature means no product in that exact model. Limited/unsupported evidence never proves absence; page omissions never imply complete display. No band intensity, gel resolution calibration, success probability, biological winner, source edit, material installation or implicit database write."
    )]
    pub(crate) async fn assess_screening_strategy(
        &self,
        Parameters(args): Parameters<ScreeningQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if !(1..=8).contains(&args.max_tests) {
            return Err(bad("screening test page size must be 1..8"));
        }
        bounded_json_limit(&args, 512 * 1024)?;
        let server = self.clone();
        self.compute(move || {
            let workspace = server.session.restore_workspace(&args.workspace)?;
            let request = args.native(&workspace)?;
            let context = screening::revision(&request).map_err(bad)?;
            let report = screening::assess(&request).map_err(bad)?;
            report.validate_binding(&context).map_err(bad)?;
            if args.test_start > report.tests.len() { return Err(bad("test offset exceeds retained screening tests")); }
            let end = args.test_start.saturating_add(args.max_tests).min(report.tests.len());
            let mut summary = json!({
                "workspace": args.workspace,
                "component_events": args.components.iter().map(|component| json!({"id": component.id,"event": component.event,"origin": component.origin})).collect::<Vec<_>>(),
                "engine_version": report.engine_version,"context_revision": report.context_revision,"report_revision": report.report_revision,
                "components": report.components,"subjects": report.subjects,"readout": report.readout,"limits": report.limits,
                "test_count": report.tests.len(),"test_start": args.test_start,"next_start": end,"done": end == report.tests.len(),
                "page_truncated": args.test_start != 0 || end != report.tests.len(),
                "tests": &report.tests[args.test_start..end],
                "cells": report.cells.iter().filter(|cell| args.test_start <= cell.test_index && cell.test_index < end).collect::<Vec<_>>(),
                "ranked_test_indices": report.ranked_test_indices,"panel": report.panel,
                "cleavage_work": report.cleavage_work,"pcr_work": report.pcr_work,"statement": report.statement,
                "persistence": "immutable session evidence only"
            });
            bounded_json_limit(&summary, 1024 * 1024 - 1024)?;
            let artifact = server.session.text_artifact(bounded_json_limit(&report, screening::MAX_REPORT_BYTES)?, "application/json", "screening-strategy.json")?;
            summary["artifact_id"] = json!(artifact);
            output(summary)
        }).await
    }
}
