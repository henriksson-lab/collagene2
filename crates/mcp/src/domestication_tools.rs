//! Source-bound reference constraints remain distinct from physical preparation.
use super::*;
use collagene_core::{domestication as design, domestication_preparation as preparation};
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DesignQuery {
    pub source: DocumentRef,
    pub request: design::Request,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckQuery {
    pub source: DocumentRef,
    pub request: design::Request,
    pub selection: design::Selection,
}
#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PrepareQuery {
    pub source: DocumentRef,
    pub request: preparation::Request,
}
impl CollageneServer {
    fn domestication_source(
        &self,
        source: &DocumentRef,
        request: &design::Request,
    ) -> Result<Molecule, ErrorData> {
        let store = self.store.lock().unwrap();
        let actual = store.document(source)?;
        design::preflight(actual, request).map_err(bad)?;
        Ok(actual.clone())
    }
    pub(super) async fn domestication_inspect(
        &self,
        args: DesignQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.domestication_source(&args.source,&args.request)?;
            let report=design::inspect(&source,&args.request).map_err(bad)?;
            report.validate_binding(&design::inspect_revision(&source,&args.request).map_err(bad)?).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&json!({"source":args.source,"inspection":report}),32*1024*1024)?,"application/json","domestication-inspection.json")?;
            output(json!({"state":"reference-constraint-inspection","verified":false,"context_revision":report.context_revision,
                "site_count":report.sites.len(),"sites":report.sites.iter().take(40).collect::<Vec<_>>(),
                "feature_count":report.features.len(),"features":report.features.iter().take(20).collect::<Vec<_>>(),
                "blocking":report.blocking,"warnings":report.warnings,"artifact_id":artifact}))
        }).await
    }
    pub(super) async fn domestication_propose(
        &self,
        args: DesignQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.domestication_source(&args.source,&args.request)?;
            let report=design::propose(&source,&args.request).map_err(bad)?;
            report.validate_binding(&design::propose_revision(&source,&args.request).map_err(bad)?).map_err(bad)?;
            let artifact=server.text_artifact(bounded_json_limit(&json!({"source":args.source,"proposal":report}),32*1024*1024)?,"application/json","domestication-proposal.json")?;
            output(json!({"state":"advisory-synonymous-proposal","verified":false,"context_revision":report.context_revision,
                "status":report.status,"search_complete":report.search_complete,"candidate":report.candidate,
                "searches":report.searches,"blocking":report.inspection.blocking,"warnings":report.warnings,"artifact_id":artifact,
                "authority":"Reference design only. Search pruning cannot prove a minimum or impossibility. Independently check the literal patches and explicitly materialize a real preparation; no edited stock is invented."}))
        }).await
    }
    pub(super) async fn domestication_check(
        &self,
        args: CheckQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.domestication_source(&args.source,&args.request)?;
            let checked=design::check(&source,&args.request,&args.selection).map_err(bad)?;
            checked.report().validate_binding(&design::check_revision(&source,&args.request,&args.selection).map_err(bad)?).map_err(bad)?;
            let evidence=json!({"source":args.source,"report":checked.report(),"target":checked.target()});
            let (artifact,target)=server.candidate_document_artifact_bounded(&evidence,Some(checked.target().clone()),"domestication-check.json",32*1024*1024)?;
            output(json!({"state":"checked-reference-intent","verified":false,"context_revision":checked.report().context_revision,
                "source":args.source,"selection":args.selection,"target":target,"artifact_id":artifact,
                "authority":"Literal substitutions, all affected supported coding features, protected bases and scoped sites were freshly checked. This independently declared reference is not physically prepared DNA or experimental verification."}))
        }).await
    }
    pub(super) async fn domestication_prepare(
        &self,
        args: PrepareQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let source=server.domestication_source(&args.source,&args.request.design)?;
            let input=PlanInput{id:"template".into(),revision:source.revision(),molecule:source};
            let prepared=preparation::materialize(&input,&args.request).map_err(bad)?;
            prepared.validate_binding(&preparation::revision(&input,&args.request).map_err(bad)?).map_err(bad)?;
            let target=prepared.target.clone();
            let plan=prepared.plan.clone();
            let plan_revision=verification::plan_revision(&plan);
            if plan_revision.is_empty()||plan.target_revision!=target.revision(){return Err(bad("prepared plan lost its independent target"));}
            #[derive(Serialize)]
            struct Evidence<'a>{source:&'a DocumentRef,prepared:&'a preparation::Prepared}
            let text=bounded_json_limit(&Evidence{source:&args.source,prepared:&prepared},32*1024*1024)?;
            let method=match &prepared.method{preparation::MethodEvidence::InversePcr{..}=>"InversePcr",preparation::MethodEvidence::OverlapExtension{..}=>"OverlapExtension",preparation::MethodEvidence::InversePcrRounds{..}=>"InversePcrRounds"};
            let artifact_id=blake3::hash(text.as_bytes()).to_hex().to_string();
            let target_id=target.revision();
            let target_bytes=bounded_json(&target)?.len();
            let plan_bytes=bounded_json(&plan)?.len().checked_add(target_bytes).ok_or_else(||bad("prepared candidate size overflow"))?;
            let result=output(json!({"state":if prepared.proposal_ready{"candidate"}else{"diagnostic-candidate"},"verified":false,"source":args.source,
                "proposal_ready":prepared.proposal_ready,"proposal_findings":prepared.proposal_findings.iter().take(40).collect::<Vec<_>>(),"proposal_findings_total":prepared.proposal_findings.len(),"proposal_findings_truncated":prepared.proposal_findings.len()>40,"proposal_scope":"Requested primer-construction conditions only; separate from molecular verification, isolation or experimental success.",
                "context_revision":prepared.context_revision,"plan_revision":plan_revision,"method":method,"operation_count":plan.steps.len(),
                "target":summary(&target_id,&target),"verification":{"status":prepared.verification.status,
                    "engine_version":prepared.verification.engine_version,"findings":prepared.verification.findings.iter().take(40).collect::<Vec<_>>(),"findings_total":prepared.verification.findings.len(),"findings_truncated":prepared.verification.findings.len()>40},
                "warnings":prepared.warnings,"artifact_id":artifact_id,
                "authority":"Native original-source preparation proposal against the separately frozen reference. A diagnostic fallback is not ready under the requested primer-construction conditions. Actual verification remains separate. Use verify_plan, export_protocol, export_protocol_pdf and export_plan_oligos with plan_revision. No efficiency, template clearance or experimental success is established."}))?;
            let mut store=server.store.lock().unwrap();
            let new_plan=!store.plans.contains_key(&plan_revision);
            let new_doc=!store.documents.contains_key(&target_id);
            let new_artifact=!store.artifacts.contains_key(&artifact_id);
            if new_plan&&store.plans.len()>=MAX_ITEMS||new_doc&&store.documents.len()>=MAX_DOCUMENT_REVISIONS||new_artifact&&store.artifacts.len()>=MAX_ITEMS{return Err(bad("domestication session capacity; no outputs published"));}
            if let Some(old)=store.artifacts.get(&artifact_id){if !matches!(old.data,ArtifactData::Text(_))||old.mime_type!="application/json"||old.filename!="domestication-preparation.json"{return Err(bad("prepared evidence has another artifact representation"));}}
            let bytes=(if new_plan{plan_bytes}else{0}).checked_add(if new_doc{target_bytes}else{0}).and_then(|n|n.checked_add(if new_artifact{text.len()}else{0})).ok_or_else(||bad("prepared output size overflow"))?;
            store.reserve(bytes)?;
            if new_plan{store.plans.insert(plan_revision,StoredPlan{plan,target:target.clone()});}
            if new_doc{store.documents.insert(target_id,target);}
            if new_artifact{store.artifacts.insert(artifact_id,Artifact{data:ArtifactData::Text(text),mime_type:"application/json",encoding:"utf-8",filename:"domestication-preparation.json"});}
            Ok(result)
        }).await
    }
}
