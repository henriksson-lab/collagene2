//! Immutable-source construct QC, assembly-standard and oligo-collection tools.
use crate::{
    server::{
        CollageneServer, EmptyRequest,
        response::{bad, bounded_json_limit, output},
    },
    session::references::DocumentRef,
};
use collagene_core::cloning::{
    assembly::standards,
    pcr::{collection_binding, primers},
};
use collagene_core::construct_qc;
use collagene_model::Molecule;
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::json;

fn default_page() -> usize {
    25
}
fn page(total: usize, start: usize, count: usize) -> Result<usize, ErrorData> {
    if !(1..=100).contains(&count) || start > total {
        return Err(bad(
            "page size must be 1..100 and offset cannot exceed the retained result count",
        ));
    }
    Ok(start.saturating_add(count).min(total))
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ConstructQcQuery {
    pub document: DocumentRef,
    pub options: construct_qc::Options,
    #[serde(default)]
    pub finding_start: usize,
    #[serde(default = "default_page")]
    pub max_findings: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct StandardPart {
    pub id: String,
    pub document: DocumentRef,
    pub start: usize,
    pub length: usize,
    pub reverse: bool,
    pub slot: String,
    pub body_start: usize,
    pub body_end: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct StandardQuery {
    pub profile_id: String,
    pub profile_identity: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, StandardPart, 128>")]
    pub parts: Vec<StandardPart>,
    pub limits: collagene_assembly_standards::assessment::Limits,
    #[serde(default)]
    pub finding_start: usize,
    #[serde(default = "default_page")]
    pub max_findings: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OligoCollectionQuery {
    pub document: DocumentRef,
    pub collection: collagene_oligo_binding::Collection,
    pub options: collagene_oligo_binding::Options,
    pub budget: collagene_oligo_binding::BudgetLimits,
    pub conditions: primers::ThermoConditions,
    #[serde(default)]
    pub hit_start: usize,
    #[serde(default = "default_page")]
    pub max_hits: usize,
}
/// Snapshot admission before owned copies. Scientific interpretation remains in
/// the core owners; repeated document references share this aggregate charge.
fn source_preflight(sources: &[&Molecule]) -> Result<(), ErrorData> {
    collagene_model::bounds::model_preflight(sources).map_err(bad)?;
    let mut bases = 0_usize;
    for source in sources {
        bases = bases
            .checked_add(source.sequence.len())
            .ok_or_else(|| bad("source snapshot size overflow"))?;
        if bases > crate::session::MAX_STORE_BYTES {
            return Err(bad(
                "source snapshots exceed the ordinary session byte bound",
            ));
        }
    }
    Ok(())
}

#[tool_router(router = design_checks_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "Check declared construct structure against an exact immutable nucleotide document, annotations, genetic code and explicit source/revision-bound annotation/reference/fusion/junction assumptions, declared initiation-context patterns and optional bounded ORF scans. Returns paged coordinates/findings and a complete bounded JSON artifact. Pass/fail/indeterminate are individual structural checks, never expression, packaging, titre, experimental QC or annotation truth. No sequence edits or persistent writes."
    )]
    pub(crate) async fn check_construct_qc(
        &self,
        Parameters(args): Parameters<ConstructQcQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        page(usize::MAX, args.finding_start, args.max_findings)?;
        let server = self.clone();
        self.compute(move || {
            let molecule = {
                let store = server.session.store.lock().unwrap();
                let molecule = store.document(&args.document)?;
                args.options.validate(molecule).map_err(bad)?;
                source_preflight(&[molecule])?;
                molecule.clone()
            };
            let report = construct_qc::assess(&molecule, &args.options).map_err(bad)?;
            let end = page(report.findings.len(), args.finding_start, args.max_findings)?;
            let artifact_id = server.session.text_artifact(
                bounded_json_limit(&report, construct_qc::MAX_REPORT_BYTES)?,
                "application/json", "construct-qc.json",
            )?;
            output(json!({"document":args.document,"engine_version":report.engine_version,
                "context_revision":report.context_revision,"source_revision":report.source_revision,
                "annotations_revision":report.annotations_revision,"options_revision":report.options_revision,
                "references_revision":report.references_revision,"report_revision":report.report_revision,
                "finding_count":report.findings.len(),"finding_start":args.finding_start,
                "next_start":end,"done":end==report.findings.len(),
                "findings":&report.findings[args.finding_start..end],
                "counts":{"pass":report.findings.iter().filter(|f|f.verdict==construct_qc::Verdict::Pass).count(),
                    "fail":report.findings.iter().filter(|f|f.verdict==construct_qc::Verdict::Fail).count(),
                    "indeterminate":report.findings.iter().filter(|f|f.verdict==construct_qc::Verdict::Indeterminate).count(),
                    "informational":report.findings.iter().filter(|f|f.verdict==construct_qc::Verdict::Informational).count()},
                "statement":report.statement,"artifact_id":artifact_id,"persistence":"session only"}))
        }).await
    }

    #[tool(
        description = "List the exact installed assembly-standard catalogue, profile/catalog identities, actual supported part-sequence scope and source evidence. The sole built-in is the sourced RFC10 part-sequence (2007) profile; this is not full vector/host/registry compliance or a guessed MoClo/SEVA/RFC1000 profile. No document changes or assembly simulation."
    )]
    pub(crate) async fn list_assembly_standards(
        &self,
        Parameters(_args): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        self.compute(move || {
            let profiles = collagene_core::cloning::assembly::standards::catalogue().map_err(bad)?;
            if profiles.len() > 32 { return Err(bad("assembly-standard catalogue exceeds the discovery bound")); }
            bounded_json_limit(&profiles, 512 * 1024)?;
            output(json!({"profile_count":profiles.len(),"profiles":profiles,
                "authority":"Installed documentary part-sequence profiles only; no experimental success or complete-standard compliance"}))
        }).await
    }

    #[tool(
        description = "Check an ordered set of immutable document fragments against an exact discovered assembly-standard profile identity and declared slot/body boundaries. Fragment start/length use original forward source coordinates; reverse and a single circular wrap are explicit; body coordinates index the oriented extracted part. Returns bounded findings, mapped original-source site intervals, local grammar successors and complete JSON evidence. Sole installed profile is RFC10 part-sequence (2007), not full registry/vector/host compliance, functional ordering or physical assembly verification. No sequence edits or persistence."
    )]
    pub(crate) async fn check_assembly_standard(
        &self,
        Parameters(args): Parameters<StandardQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        page(usize::MAX, args.finding_start, args.max_findings)?;
        if args.parts.len() > 128 {
            return Err(bad("at most128 declared standard parts"));
        }
        let server = self.clone();
        self.compute(move || {
            let parts = {
                let store = server.session.store.lock().unwrap();
                let sources = args.parts.iter().map(|p| store.document(&p.document)).collect::<Result<Vec<_>,_>>()?;
                source_preflight(&sources)?;
                args.parts.iter().zip(sources).map(|(p, molecule)| standards::SourcePart {
                    id:p.id.clone(),molecule:molecule.clone(),start:p.start,length:p.length,reverse:p.reverse,
                    slot:p.slot.clone(),body_start:p.body_start,body_end:p.body_end,
                }).collect()
            };
            let request = standards::Request { profile_id:args.profile_id,profile_identity:args.profile_identity,parts,limits:args.limits };
            let report = standards::assess(&request).map_err(bad)?;
            let end = page(report.assessment.findings.len(), args.finding_start, args.max_findings)?;
            let mut findings = Vec::with_capacity(end-args.finding_start);
            for finding in &report.assessment.findings[args.finding_start..end] {
                let source_spans = match (&finding.kind, finding.part_index) {
                    (collagene_assembly_standards::assessment::FindingKind::AvoidedSite {start,end,..},Some(index)) => {
                        let binding = report.sources.get(index).ok_or_else(||bad("standard finding source index is out of bounds"))?;
                        Some(binding.map_interval(*start,*end).map_err(bad)?)
                    },
                    _ => None,
                };
                findings.push(json!({"finding":finding,"source_spans":source_spans}));
            }
            let artifact_id = server.session.text_artifact(bounded_json_limit(&report,standards::MAX_REPORT_BYTES)?,"application/json","assembly-standard.json")?;
            output(json!({"engine_version":report.engine_version,"context_revision":report.context_revision,
                "profile":report.profile,"sources":report.sources,"parts":args.parts,
                "report_revision":report.report_revision,"assessment_revision":report.assessment.request_revision,
                "outcome":report.assessment.outcome,"complete":report.assessment.complete,
                "violations_observed":report.assessment.violations_observed,
                "indeterminate_observed":report.assessment.indeterminate_observed,
                "advisories_observed":report.assessment.advisories_observed,
                "budgets_hit":report.assessment.budgets_hit,
                "finding_count":report.assessment.findings.len(),"finding_start":args.finding_start,
                "next_start":end,"done":end==report.assessment.findings.len(),"findings":findings,
                "successors":report.successors,"artifact_id":artifact_id,"persistence":"session only",
                "authority":report.assessment.authority}))
        }).await
    }

    #[tool(
        description = "Search a declared bounded oligo collection against one exact immutable source document using the embedded ungapped 3-prime-anchored annealing-suffix model. Returns paged original-coordinate footprints and score/mismatch/certainty evidence, native work/completeness limits, explicit solution conditions, separate annealing/whole-oligo thermodynamics and full JSON artifact. Annealing Tm uses a perfect duplex or the admitted single internal DNA-mismatch model with pH7, sodium-only correction and two distinct equal-concentration strands. Declared tails stay separate; unsupported mismatch/ambiguous thermodynamics are unavailable, never silent perfect matches. Whole-oligo Tm is a hypothetical perfect duplex, not actual binding; no PCR, kinetics, calibrated probability, sequence edits or persistence."
    )]
    pub(crate) async fn search_oligo_collection(
        &self,
        Parameters(args): Parameters<OligoCollectionQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        page(usize::MAX, args.hit_start, args.max_hits)?;
        args.collection.validate().map_err(bad)?;
        args.options.validate().map_err(bad)?;
        let server = self.clone();
        self.compute(move || {
            let source = {
                let store = server.session.store.lock().unwrap();
                let source = store.document(&args.document)?;
                source_preflight(&[source])?;
                source.clone()
            };
            let request = collection_binding::Request {source,collection:args.collection,options:args.options,budget:args.budget,conditions:args.conditions};
            let report = collection_binding::assess(&request).map_err(bad)?;
            let end = page(report.search.hits.len(),args.hit_start,args.max_hits)?;
            let hits = report.search.hits[args.hit_start..end].iter().enumerate().map(|(page_index,hit)| {
                let index = args.hit_start + page_index;
                let thermodynamics = report.thermodynamics.iter().find(|t|t.hit_index==index);
                json!({"hit_index":index,"oligo_index":hit.oligo_index,"oligo_id":hit.oligo_id,
                    "orientation":hit.orientation,"three_prime_source_position":hit.three_prime_source_position,
                    "source_start":hit.source_start,"length":hit.length,"crosses_origin":hit.crosses_origin,
                    "declared_tail_length":hit.declared_tail_length,"extra_unpaired_length":hit.extra_unpaired_length,
                    "annealing_start":hit.annealing_start,"minimum_mismatches":hit.minimum_mismatches,
                    "maximum_mismatches":hit.maximum_mismatches,"minimum_score":hit.minimum_score,
                    "maximum_score":hit.maximum_score,"three_prime_certain":hit.three_prime_certain,
                    "three_prime_possible":hit.three_prime_possible,"certainly_admissible":hit.certainly_admissible,
                    "score_scope":hit.score_scope,"thermodynamics":thermodynamics})
            }).collect::<Vec<_>>();
            let artifact_id = server.session.text_artifact(bounded_json_limit(&report,collection_binding::MAX_REPORT_BYTES)?,"application/json","oligo-collection-binding.json")?;
            output(json!({"document":args.document,"engine_version":report.engine_version,
                "context_revision":report.context_revision,"source_revision":report.source_revision,
                "report_revision":report.report_revision,"collection_hash":report.search.collection_hash,
                "options_hash":report.search.options_hash,"source_content_hash":report.search.source_content_hash,
                "model":report.search.model,"complete":report.search.complete,"conditions":report.conditions,
                "search_engine_version":report.search.engine_version,
                "budget_limits":report.search.budget_limits,"budget_before":report.search.budget_before,
                "budget_after":report.search.budget_after,"hit_count":report.search.hits.len(),
                "oligo_count":report.search.oligos.len(),
                "oligo_progress":report.search.oligos.iter().take(100).collect::<Vec<_>>(),
                "oligo_progress_truncated":report.search.oligos.len()>100,
                "hit_start":args.hit_start,"next_start":end,"done":end==report.search.hits.len(),
                "hits":hits,"warnings":report.warnings,"artifact_id":artifact_id,"persistence":"session only",
                "authority":"Sequence search and explicitly scoped perfect-duplex or isolated internal mismatch equilibrium thermodynamics; no PCR extension, experimental binding, calibrated probability or cloning success"}))
        }).await
    }
}
