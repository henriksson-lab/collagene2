use crate::server::CollageneServer;
use crate::server::response::bad;
use crate::server::response::bounded_json_limit;
use crate::server::response::output;
use crate::session::references::DocumentRef;
use crate::session::references::WorkspaceRef;

use collagene_model::Molecule;
use rmcp::ErrorData;
use rmcp::model::CallToolResult;
use schemars::JsonSchema;
use serde::Deserialize;
use serde::Serialize;
use serde_json::Value;
use serde_json::json;

use collagene_core::cloning::recombination::{multisite as design, reference};

pub const THREE_ENTRY_PROFILE_URI: &str =
    "collagene://recombination/profile/tol2kit-228-383-302-vectorii-lr-reference/2026-11-25.1";
pub fn three_entry_profile() -> Result<Value, ErrorData> {
    let profile = reference::multisite::three_entry_source_contexts().map_err(bad)?;
    let revision = profile.revision().map_err(bad)?;
    Ok(
        json!({"profile":profile,"profile_revision":revision,"recommended_reaction":"LR","specificities":[4,1,2,3],
        "authority":"Finite source-edition context facts, eight LR contexts plus full P3 competitor. Review exact source editions and limitations before use. No universal kit recognition, physical chemistry or experimental success is established."}),
    )
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct NativePin {
    pub workspace: WorkspaceRef,
    pub event: collagene_core::material::lineage::EventBinding,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SourceRef {
    pub role: String,
    pub document: DocumentRef,
    /// Optional exact retained native event. Ownership and literal molecule equality are checked.
    #[serde(default)]
    pub native: Option<NativePin>,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InspectQuery {
    /// Two through five ordered immutable source roles; identical DNA does not merge roles.
    pub sources: Vec<SourceRef>,
    /// Complete caller-supplied finite profile, declared specificity set and work limits.
    pub request: design::Request,
    #[serde(default)]
    pub offset: usize,
    #[serde(default = "page_size")]
    pub limit: usize,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProductQuery {
    pub sources: Vec<SourceRef>,
    pub request: design::Request,
    pub selection: collagene_core::cloning::recombination::Selection,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckQuery {
    pub sources: Vec<SourceRef>,
    pub request: design::Request,
    pub selection: collagene_core::cloning::recombination::Selection,
    /// Independent noncausal requirement; never determines site or product selection.
    pub target: DocumentRef,
}
fn page_size() -> usize {
    10
}
impl CollageneServer {
    fn multisite_sources(
        &self,
        sources: &[SourceRef],
        request: &design::Request,
        target: Option<&DocumentRef>,
    ) -> Result<(Vec<reference::Input>, Option<Molecule>), ErrorData> {
        design::request_preflight(request).map_err(bad)?;
        if !(2..=5).contains(&sources.len()) {
            return Err(bad("multisite sources require 2..5 ordered roles"));
        }
        let mut roles = std::collections::BTreeSet::new();
        for s in sources {
            if s.role.trim().is_empty() || s.role.len() > 256 || !roles.insert(&s.role) {
                return Err(bad("source roles must be nonempty, bounded and unique"));
            }
            if let Some(pin) = &s.native {
                if pin.workspace.id.len() != 64
                    || pin.workspace.revision.len() != 64
                    || pin.event.document_id.is_empty()
                    || pin.event.document_id.len() > 128
                    || pin.event.event_id.is_empty()
                    || pin.event.event_id.len() > 128
                    || pin.event.revision.len() != 64
                {
                    return Err(bad("invalid bounded native source pin"));
                }
            }
        }
        let (inputs, target) = {
            let store = self.session.store.lock().unwrap();
            let molecules = sources
                .iter()
                .map(|s| store.document(&s.document))
                .collect::<Result<Vec<_>, _>>()?;
            design::source_preflight(&molecules, request).map_err(bad)?;
            let target = target.map(|r| store.document(r)).transpose()?;
            if let Some(t) = target {
                reference::model_preflight(&[t]).map_err(bad)?;
                t.validate().map_err(bad)?;
                if t.sequence.alphabet() != collagene_model::Alphabet::Dna {
                    return Err(bad("independent target must be DNA"));
                }
            }
            (
                sources
                    .iter()
                    .zip(molecules)
                    .map(|(s, m)| reference::Input {
                        role: s.role.clone(),
                        molecule: m.clone(),
                    })
                    .collect::<Vec<_>>(),
                target.cloned(),
            )
        };
        // Group exact native pins to avoid restoring one potentially large frame per role.
        let mut groups = std::collections::BTreeMap::<String, (&WorkspaceRef, Vec<usize>)>::new();
        for (i, s) in sources.iter().enumerate() {
            if let Some(pin) = &s.native {
                if pin.workspace.id != pin.workspace.revision {
                    return Err(bad("stale native workspace revision"));
                }
                groups
                    .entry(pin.workspace.id.clone())
                    .or_insert_with(|| (&pin.workspace, Vec::new()))
                    .1
                    .push(i);
            }
        }
        for (_, (workspace, indices)) in groups {
            let frame = self.session.restore_workspace(workspace)?;
            for i in indices {
                let pin = &sources[i]
                    .native
                    .as_ref()
                    .ok_or_else(|| bad("missing native source pin"))?
                    .event;
                if !frame.documents.iter().any(|d| d.id == pin.document_id) {
                    return Err(bad("native source owner is not retained"));
                }
                let event = frame
                    .history
                    .iter()
                    .find(|e| e.id == pin.event_id)
                    .ok_or_else(|| bad("exact native source event is not retained"))?;
                if event.document_id != pin.document_id || event.revision != pin.revision {
                    return Err(bad("native source owner/event/revision mismatch"));
                }
                if pin.revision != inputs[i].molecule.revision()
                    || frame.revisions.get(&pin.revision) != Some(&inputs[i].molecule)
                {
                    return Err(bad(
                        "immutable source differs from pinned historical molecule",
                    ));
                }
            }
        }
        Ok((inputs, target))
    }
    pub(crate) async fn multisite_inspect(
        &self,
        args: InspectQuery,
    ) -> Result<CallToolResult, ErrorData> {
        if args.limit == 0 || args.limit > 25 || args.offset > 256 {
            return Err(bad("candidate pages require limit1..25 and bounded offset"));
        }
        let server = self.clone();
        self.compute(move||{
            let(inputs,_)=server.multisite_sources(&args.sources,&args.request,None)?;
            let report=design::simulate(&inputs,&args.request).map_err(bad)?;
            report.validate_binding(&design::revision(&inputs,&args.request).map_err(bad)?).map_err(bad)?;
            if args.offset>report.reference.candidates.len(){return Err(bad("candidate offset beyond report"));}
            let rows=report.reference.candidates.iter().skip(args.offset).take(args.limit).map(|c|json!({"candidate_id":c.id,"sites":c.sites,"pairs":c.pairs,"components":c.products.iter().enumerate().map(|(i,p)|json!({"component_id":i,"length":p.molecule.sequence.len(),"topology":p.molecule.topology,"reference_revision":p.molecule.revision()})).collect::<Vec<_>>()})).collect::<Vec<_>>();
            let next=args.offset+rows.len();
            let artifact=server.session.text_artifact(bounded_json_limit(&json!({"sources":args.sources,"batch":report}),design::MAX_REPORT_BYTES)?,"application/json","multisite-recombination-batch.json")?;
            output(json!({"sources":args.sources,"batch_revision":report.batch_revision,"engine_version":report.engine_version,"complete":report.reference.complete,"candidate_count":report.reference.candidates.len(),"offset":args.offset,"next_offset":next,"done":next==report.reference.candidates.len(),"candidates":rows,"unsupported":report.reference.unsupported,"limits":report.reference.limits,"warnings":report.warnings,"artifact_id":artifact,"authority":"Finite caller-profile reference alternatives only. Full typed artifact retains every reciprocal component and annotation correspondence. No preferred winner, chemical activity or inherited physical state is established."}))
        }).await
    }
    pub(crate) async fn multisite_select(
        &self,
        args: ProductQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let(inputs,_)=server.multisite_sources(&args.sources,&args.request,None)?;
            let product=design::product(&inputs,&args.request,&args.selection).map_err(bad)?;
            product.validate_binding(&args.selection.batch_revision).map_err(bad)?;
            let evidence=json!({"sources":args.sources,"request":args.request,"product":product});
            let(artifact,document)=server.session.candidate_document_artifact_bounded(&evidence,Some(product.molecule.clone()),"multisite-recombination-component.json",design::MAX_REPORT_BYTES)?;
            output(json!({"sources":args.sources,"selection":args.selection,"document":document,"reciprocal_components":product.reciprocal_components,"annotation_evidence":product.annotation_evidence,"warnings":product.warnings,"artifact_id":artifact,"verified":false,"authority":"Fresh selected reference-only projection; explicit choice is not independent target agreement or physical reaction validation. Exact source pins remain in the artifact; native causal lineage requires a typed plan/run."}))
        }).await
    }
    pub(crate) async fn multisite_check(
        &self,
        args: CheckQuery,
    ) -> Result<CallToolResult, ErrorData> {
        let server = self.clone();
        self.compute(move||{
            let(inputs,target)=server.multisite_sources(&args.sources,&args.request,Some(&args.target))?;
            let target=target.ok_or_else(||bad("independent target required"))?;
            let report=design::check_target(&inputs,&args.request,&args.selection,&target).map_err(bad)?;
            report.validate_binding(&args.selection.batch_revision).map_err(bad)?;
            let artifact=server.session.text_artifact(bounded_json_limit(&json!({"sources":args.sources,"target":args.target,"request":args.request,"checked":report}),design::MAX_REPORT_BYTES)?,"application/json","multisite-recombination-target-check.json")?;
            output(json!({"sources":args.sources,"target":args.target,"selection":args.selection,"context_revision":report.context_revision().map_err(bad)?,"assessment":report.assessment,"reciprocal_components":report.product.reciprocal_components,"artifact_id":artifact,"authority":"Fresh comparison of the explicitly selected reference component against independent oriented DNA/topology. Target is noncausal; mismatch and incomplete search remain visible. No document is installed by this check and no chemical success is asserted."}))
        }).await
    }
}
