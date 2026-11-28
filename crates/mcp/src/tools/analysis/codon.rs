//! Source-bound codon design and the pinned observed-host preference catalogue.
use crate::{
    server::{
        CollageneServer, EmptyRequest,
        response::{bad, bounded_json_limit, output},
    },
    session::references::DocumentRef,
};
use collagene_core::sequence::codon_design::{CodonPreferences, hosts};
use rmcp::{
    ErrorData, handler::server::wrapper::Parameters, model::CallToolResult, tool, tool_router,
};
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::json;

const HOST_AUTHORITY: &str = "Observed internal reference-CDS codon preferences under the pinned admission policy. These are not CAI, calibrated gene expression, host fitness, experimental success, physical synthesis or cloning verification.";

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CodonDesignQuery {
    pub document: DocumentRef,
    pub options: collagene_core::sequence::codon_design::DesignOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CodonCheckQuery {
    pub source: DocumentRef,
    pub candidate: DocumentRef,
    pub options: collagene_core::sequence::codon_design::DesignOptions,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CodonHostQuery {
    pub reference: hosts::Reference,
}

#[tool_router(router = codon_router, vis = "pub(crate)")]
impl CollageneServer {
    #[tool(
        description = "List compact descriptors for the embedded curated reference-host codon catalogue. Each selection has an exact profile id and identity, organism/assembly/code, admitted/excluded loci, internal-codon count and admission policy. This reads small pinned data, never parses a runtime genome, contacts a provider or predicts expression/experimental success."
    )]
    pub(crate) async fn list_codon_hosts(
        &self,
        Parameters(_args): Parameters<EmptyRequest>,
    ) -> Result<CallToolResult, ErrorData> {
        self.compute(move || {
            let catalogue = hosts::catalogue().map_err(bad)?;
            output(json!({"host_count":catalogue.len(),"hosts":catalogue,
                "scope":"Embedded curated reference-host catalogue only", "expression_prediction":false,
                "calibrated":false,"authority":HOST_AUTHORITY}))
        }).await
    }

    #[tool(
        description = "Open one exact embedded reference-host codon profile by pinned id and identity. Returns a compact descriptor, a ReferenceHost preferences value accepted by design_codons/check_codon_design, and a full JSON artifact of observed internal/start/stop counts, explicit codon ordering, resolved integer preferences, source hashes/URLs and derivation policy. No user-supplied organism label/profile is trusted; no runtime genome parsing, provider request, expression/CAI/experimental prediction or molecular authority."
    )]
    pub(crate) async fn codon_host_profile(
        &self,
        Parameters(args): Parameters<CodonHostQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        if args.reference.id.is_empty()
            || args.reference.id.len() > 128
            || args.reference.identity.len() != 64
            || !args
                .reference
                .identity
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(bad(
                "bounded pinned codon-host id and identity are required",
            ));
        }
        let server = self.clone();
        self.compute(move || {
            let profile = hosts::resolve(&args.reference).map_err(bad)?;
            let descriptor = hosts::catalogue()
                .map_err(bad)?
                .into_iter()
                .find(|host| {
                    host.reference.id == args.reference.id
                        && host.reference.identity == args.reference.identity
                })
                .ok_or_else(|| bad("host reference is absent from the current pinned catalogue"))?;
            let weights = profile.preferences().map_err(bad)?;
            let preferences = CodonPreferences::ReferenceHost {
                reference: args.reference.clone(),
            };
            let full = json!({"reference":args.reference,"profile":profile,
                "codon_order":weights.iter().map(|weight|weight.codon.as_str()).collect::<Vec<_>>(),
                "resolved_preferences":weights,"preferences":preferences,
                "runtime_source":"Embedded small pinned profile; no runtime genome parsing",
                "expression_prediction":false,"calibrated":false,"authority":HOST_AUTHORITY});
            let text = bounded_json_limit(&full, 64 * 1024)?;
            let artifact_id = server.session.text_artifact(
                text,
                "application/json",
                "codon-host-profile.json",
            )?;
            output(
                json!({"host":descriptor,"preferences":preferences,"artifact_id":artifact_id,
                "encoding":"utf-8","expression_prediction":false,"calibrated":false,
                "authority":HOST_AUTHORITY,"persistence":"session only"}),
            )
        })
        .await
    }

    #[tool(
        description = "Design synonymous DNA for an immutable concrete protein or exact DNA coding source using explicit genetic code, Uniform/Explicit/pinned ReferenceHost codon preferences and synthesis constraints. Deterministic bounded search independently checks every returned candidate against the source protein; incomplete search never claims unsatisfiability or a global optimum. Contexts constrain junctions but are excluded from the candidate CDS. Returns a bounded preview, full JSON artifact and immutable candidate document reference when successful. No organism expression, physical synthesis or cloning verification claim; no persistent application writes."
    )]
    pub(crate) async fn design_codons(
        &self,
        Parameters(args): Parameters<CodonDesignQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let source = server.session.store.lock().unwrap().document(&args.document)?.clone();
            let report = collagene_core::sequence::codon_design::design(&source, &args.options).map_err(bad)?;
            report.validate_for_source(&source, &args.options).map_err(bad)?;
            let (artifact_id, candidate_document) = server.session.codon_artifact(&report, report.candidate.clone())?;
            output(json!({"document":args.document,"engine_version":report.engine_version,"source_revision":report.source_revision,"options_revision":report.options_revision,"preferences":args.options.preferences,"target_protein_preview":report.target_protein.chars().take(200).collect::<String>(),"target_protein_length":report.target_protein.len(),"target_preview_truncated":report.target_protein.len()>200,"status":report.status,"search_complete":report.search_complete,"optimality_proven":report.optimality_proven,"work_used":report.work_used,"pruned_states":report.pruned_states,"score":report.score,"candidate":candidate_document,"candidate_sequence_preview":report.candidate.as_ref().map(|m|m.sequence.as_str().chars().take(200).collect::<String>()),"candidate_preview_truncated":report.candidate.as_ref().is_some_and(|m|m.sequence.len()>200),"candidate_check":report.candidate_check.as_ref().map(|check|json!({"complete":check.complete,"satisfies":check.satisfies,"candidate_revision":check.candidate_revision,"gc_count":check.gc_count,"coding_length":check.coding_length,"score":check.score})),"limits":report.limits.iter().take(100).collect::<Vec<_>>(),"limits_truncated":report.limits.len()>100,"warnings":report.warnings.iter().take(100).collect::<Vec<_>>(),"warnings_truncated":report.warnings.len()>100,"artifact_id":artifact_id,"persistence":"session only","authority":"Candidate sequence constraints are checked against the immutable required protein. No organism expression, physical synthesis or cloning verification is established."}))
        }).await
    }
    #[tool(
        description = "Independently check a supplied immutable DNA candidate against a separately bound concrete protein or exact DNA coding source and explicit synthesis constraints. Re-derives the required protein, translation, codon preference score, GC, forbidden motifs on both strands and supplied context junctions. Search is not used. Incomplete bounded checks never grant a positive pass. Returns bounded evidence and a full JSON artifact; no edits or persistent writes."
    )]
    pub(crate) async fn check_codon_design(
        &self,
        Parameters(args): Parameters<CodonCheckQuery>,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move || {
            let (source, candidate) = { let store = server.session.store.lock().unwrap();
                (store.document(&args.source)?.clone(), store.document(&args.candidate)?.clone()) };
            let report = collagene_core::sequence::codon_design::check(&source, &candidate, &args.options).map_err(bad)?;
            report.validate_for(&source, &candidate, &args.options).map_err(bad)?;
            let (artifact_id, _) = server.session.codon_artifact(&report, None)?;
            output(json!({"source":args.source,"candidate":args.candidate,"engine_version":report.engine_version,"source_revision":report.source_revision,"candidate_revision":report.candidate_revision,"options_revision":report.options_revision,"preferences":args.options.preferences,"target_protein_preview":report.target_protein.chars().take(200).collect::<String>(),"translated_protein_preview":report.translated_protein.chars().take(200).collect::<String>(),"target_protein_length":report.target_protein.len(),"translated_protein_length":report.translated_protein.len(),"protein_previews_truncated":report.target_protein.len()>200||report.translated_protein.len()>200,"complete":report.complete,"satisfies":report.satisfies,"gc_count":report.gc_count,"coding_length":report.coding_length,"score":report.score,"violations":report.violations.iter().take(100).collect::<Vec<_>>(),"violations_truncated":report.violations.len()>100,"work_used":report.work_used,"warnings":report.warnings.iter().take(100).collect::<Vec<_>>(),"warnings_truncated":report.warnings.len()>100,"artifact_id":artifact_id}))
        }).await
    }
}
