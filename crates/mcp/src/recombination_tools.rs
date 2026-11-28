//! Explicit reference-only BP/LR requests; source documents are bounded before copying.
use super::*;
use collagene_core::recombination::{self, reference};

pub const PROFILE_URI: &str = "collagene://recombination/profile/manufacturer-context-v1";
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SourceRef {
    pub role: String,
    pub document: DocumentRef,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InspectQuery {
    /// Exactly two ordered, explicitly named source roles. Equal DNA does not merge roles.
    pub sources: [SourceRef; 2],
    /// Explicit full-site profile snapshot, reaction and resource/annotation limits.
    pub request: recombination::Request,
    /// Optional typed plan source slots. When supplied, exact core plan role derivation is used; clients never compute private role hashes.
    #[serde(default)]
    pub plan_sources: Option<[verification::Source; 2]>,
    #[serde(default)]
    pub offset: usize,
    #[serde(default = "page_size")]
    pub limit: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckQuery {
    pub sources: [SourceRef; 2],
    pub request: recombination::Request,
    /// Optional typed plan source slots. When supplied, exact core plan role derivation is used; clients never compute private role hashes.
    #[serde(default)]
    pub plan_sources: Option<[verification::Source; 2]>,
    /// Independent requirement; never copied into a product or used as causal DNA.
    pub target: DocumentRef,
    pub policy: reference::TargetPolicy,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SelectQuery {
    pub sources: [SourceRef; 2],
    pub request: recombination::Request,
    /// Optional typed plan source slots. When supplied, exact core plan role derivation is used; clients never compute private role hashes.
    #[serde(default)]
    pub plan_sources: Option<[verification::Source; 2]>,
    pub selection: recombination::Selection,
    pub target: DocumentRef,
    pub policy: reference::TargetPolicy,
}
fn page_size() -> usize {
    10
}
pub fn profile() -> Result<Value, ErrorData> {
    let snapshot = recombination::ProfileSnapshot::manufacturer_context().map_err(bad)?;
    Ok(
        json!({"snapshot":snapshot,"authority":"Explicit opt-in deposited full-site contexts only; not a universal site definition, chemistry, kinetics or manufacturer activity guarantee."}),
    )
}
impl CollageneServer {
    fn recombination_sources(
        &self,
        sources: &[SourceRef; 2],
        request: &recombination::Request,
        plan_sources: Option<&[verification::Source; 2]>,
        target: Option<&DocumentRef>,
    ) -> Result<(Vec<reference::Input>, Option<Molecule>), ErrorData> {
        let store = self.store.lock().unwrap();
        let a = store.document(&sources[0].document)?;
        let b = store.document(&sources[1].document)?;
        let target = target.map(|r| store.document(r)).transpose()?;
        recombination::borrowed_preflight(
            [(sources[0].role.as_str(), a), (sources[1].role.as_str(), b)],
            request,
            target,
        )
        .map_err(bad)?;
        let inputs = match plan_sources {
            Some(slots) => {
                verification::recombination_inputs(slots, [a, b], request).map_err(bad)?
            }
            None => vec![
                reference::Input {
                    role: sources[0].role.clone(),
                    molecule: a.clone(),
                },
                reference::Input {
                    role: sources[1].role.clone(),
                    molecule: b.clone(),
                },
            ],
        };
        Ok((inputs, target.cloned()))
    }
    pub(super) async fn recombination_inspect(
        &self,
        args: InspectQuery,
    ) -> Result<CallToolResult, ErrorData> {
        if args.limit == 0 || args.limit > 25 {
            return Err(bad("candidate pages require 1..25 rows"));
        }
        let server = self.clone();
        self.compute(move||{
            let(inputs,_)=server.recombination_sources(&args.sources,&args.request,args.plan_sources.as_ref(),None)?;
            let report=recombination::simulate(&inputs,&args.request).map_err(bad)?;
            report.validate_bindings_for(&inputs,&args.request).map_err(bad)?;
            if args.offset>report.reference.candidates.len(){return Err(bad("candidate offset beyond report"));}
            let rows=report.reference.candidates.iter().skip(args.offset).take(args.limit).map(|c|
                json!({"candidate_id":c.id,"sites":c.sites,"components":c.products.iter().enumerate().map(|(i,p)|json!({"component_id":i,"length":p.molecule.sequence.len(),"topology":p.molecule.topology,"reference_revision":p.molecule.revision()})).collect::<Vec<_>>()})).collect::<Vec<_>>();
            let next_offset=args.offset+rows.len();
            let artifact_id=server.text_artifact(bounded_json(&report)?,"application/json","recombination-inspection.json")?;
            output(json!({"sources":args.sources,"plan_sources":args.plan_sources,"request":args.request,"batch_revision":report.batch_revision,"complete":report.reference.complete,"candidate_count":report.reference.candidates.len(),"offset":args.offset,"next_offset":next_offset,"done":next_offset==report.reference.candidates.len(),"candidates":rows,"diagnostics":report.reference.diagnostics,"warnings":report.warnings,"artifact_id":artifact_id,"authority":"Reference-sequence alternatives under the declared full-site profile only. All reciprocal components remain in the full artifact; no preferred biological winner or chemical activity is inferred."}))
        }).await
    }
    pub(super) async fn recombination_check(
        &self,
        args: CheckQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let(inputs,target)=server.recombination_sources(&args.sources,&args.request,args.plan_sources.as_ref(),Some(&args.target))?;
            let target=target.ok_or_else(||bad("independent target required"))?;
            let assessment=recombination::check_target(&inputs,&args.request,&target,&args.policy).map_err(bad)?;
            assessment.validate_bindings_for(&inputs,&args.request,&target,&args.policy).map_err(bad)?;
            let artifact_id=server.text_artifact(bounded_json(&assessment)?,"application/json","recombination-target-check.json")?;
            output(json!({"sources":args.sources,"plan_sources":args.plan_sources,"target":args.target,"request":args.request,"assessment":assessment,"artifact_id":artifact_id,"authority":"Fresh independent target comparison under the declared profile, retaining ambiguity and incomplete-search flags; no experimental or chemical verification."}))
        }).await
    }
    pub(super) async fn recombination_select(
        &self,
        args: SelectQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let(inputs,target)=server.recombination_sources(&args.sources,&args.request,args.plan_sources.as_ref(),Some(&args.target))?;
            let target=target.ok_or_else(||bad("independent target required"))?;
            let checked=recombination::checked_product(&inputs,&args.request,&args.selection,&target,&args.policy).map_err(bad)?;
            checked.validate_bindings_for(&inputs,&args.request,&args.selection,&target,&args.policy).map_err(bad)?;
            let(artifact_id,document)=server.candidate_document_artifact(&checked,Some(checked.product.molecule.clone()),"recombination-component.json")?;
            output(json!({"sources":args.sources,"plan_sources":args.plan_sources,"target":args.target,"request":args.request,"selection":args.selection,"batch_revision":checked.batch_revision,"check":checked.check,"reciprocal_components":checked.product.reciprocal_components,"annotation_evidence":checked.product.annotation_evidence,"warnings":checked.product.warnings,"document":document,"artifact_id":artifact_id,"authority":"Freshly reconstructed selected reference component matching the independent target. This is not physical reaction validation; product chemistry/state is unmapped. Target is noncausal. Native causal recording requires the typed plan/run workflow."}))
        }).await
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn report_and_reference_document_capacity_failure_is_atomic() {
        let server = CollageneServer::default();
        let m = Molecule::dna(
            "reference-only capacity fixture",
            "ACGT",
            collagene_model::Topology::Linear,
        )
        .unwrap();
        {
            let mut store = server.store.lock().unwrap();
            store.used = MAX_STORE_BYTES - 1;
        }
        assert!(
            server
                .candidate_document_artifact(
                    &json!({"reference":"bounded evidence"}),
                    Some(m),
                    "recombination-component.json"
                )
                .is_err()
        );
        let store = server.store.lock().unwrap();
        assert_eq!(store.used, MAX_STORE_BYTES - 1);
        assert!(store.documents.is_empty());
        assert!(store.artifacts.is_empty());
    }
}
