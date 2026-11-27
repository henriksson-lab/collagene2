//! Explicit edits of candidate declarations. Scientific evidence is recomputed
//! after Apply; graph edits never promote a displayed or imported product.
use crate::digest::SourceInfo;
use collagene_core::{
    analysis::Primer,
    chemistry::{AnnealOptions, EndRepairMethod, EndTreatment, FivePrime, ThreePrime},
    plan_editing,
    run_lineage::EventBinding,
    tailing::TailingBase,
    transfer::Candidate,
    verification::{self, DesignPlan, Operation, OrientedSource, PlanInput, Source},
};
use leptos::prelude::*;
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

fn title(op: &Operation) -> &'static str {
    match op {
        Operation::ReconstructRunProduct { .. } => "Reuse saved cloning material",
        Operation::ReconstructDigestFragment { .. } => "Reuse saved digest fragment",
        Operation::ReconstructMaterial { .. } => "Reuse target-free derived material",
        Operation::Pcr { .. } => "PCR",
        Operation::OverlapExtension { .. } => "Reciprocal overlap extension",
        Operation::PcrWithMismatches { .. } => "PCR with mismatches",
        Operation::Digest { .. } => "Restriction digest",
        Operation::CatalogDigest { .. } => "Catalog-pinned digest",
        Operation::CheckedDigest { .. } => "Modification-aware digest",
        Operation::TreatEnds { .. } => "End treatment",
        Operation::RepairEnds { .. } => "End repair",
        Operation::Tail3Prime { .. } => "3′ tailing",
        Operation::Ligate { .. } => "Ligation with bounded alternatives",
        Operation::AnnealOligos { .. } => "Oligo annealing",
        Operation::HomologyAssembly { .. } => "Physical homology assembly",
        Operation::ResectThreePrime { .. } => "Explicit 3′ source resection",
        Operation::HomologyAnnealing { .. } => {
            "3′ resection and unsealed 5′-tail annealing · incomplete reference"
        }
        Operation::Overlap { .. } | Operation::ProfiledOverlap { .. } => "Overlap assembly",
        Operation::GoldenGate { .. } | Operation::SelectedGoldenGate { .. } => "Golden Gate",
        Operation::DeclareSourceModifications { .. } => "Source modification assumptions",
        Operation::DeclareSourceTermini { .. } => "Source terminal assumptions",
        Operation::MultiSiteRecombineReference { .. } => "Multi-parent reference exchange",
        Operation::RecombineReference { .. } => "Reference recombination",
        Operation::Topo { .. } | Operation::CheckedTopo { .. } => "TOPO reference inspection",
        Operation::ProcessedTopo { .. } => "Processed TOPO",
        Operation::TopoActions { .. } => "TOPO processing actions",
    }
}
fn source_key(source: &Source) -> String {
    match source {
        Source::Input(id) => format!("input:{id}"),
        Source::Step(i) => format!("step:{i}"),
    }
}
fn parse_source(key: &str) -> Result<Source, String> {
    if let Some(id) = key.strip_prefix("input:") {
        if !id.is_empty() {
            return Ok(Source::Input(id.into()));
        }
    }
    if let Some(index) = key
        .strip_prefix("step:")
        .and_then(|s| s.parse::<usize>().ok())
    {
        return Ok(Source::Step(index));
    }
    Err("Choose an input or an earlier step product.".into())
}
fn binding_key(binding: &EventBinding) -> String {
    format!(
        "{}|{}|{}",
        binding.document_id, binding.event_id, binding.revision
    )
}
fn short(value: &str) -> String {
    value.chars().take(12).collect()
}
fn error(e: impl std::fmt::Display) -> String {
    e.to_string()
}
fn number(value: &str, label: &str) -> Result<usize, String> {
    value
        .trim()
        .parse()
        .map_err(|_| format!("{label} must be a nonnegative integer."))
}
fn first_source(ids: &[String], index: usize) -> Option<Source> {
    ids.first()
        .map(|id| Source::Input(id.clone()))
        .or_else(|| index.checked_sub(1).map(Source::Step))
}
fn new_operation(kind: &str, source: Source) -> Operation {
    match kind {
        "resect-three" => Operation::ResectThreePrime {
            source,
            cuts: collagene_core::three_prime_resection::Cuts { left: 0, right: 0 },
        },
        "digest" => Operation::Digest {
            template: source,
            enzymes: vec![],
            fragment_index: 0,
        },
        "treatment" => Operation::TreatEnds {
            template: source,
            treatment: EndTreatment::Phosphorylate,
        },
        "repair" => Operation::RepairEnds {
            template: source,
            method: EndRepairMethod::Blunt,
        },
        "tail" => Operation::Tail3Prime {
            template: source,
            base: TailingBase::A,
        },
        "ligate" => Operation::Ligate {
            options: Default::default(),
            fragments: vec![OrientedSource {
                source,
                reverse: false,
            }],
            circular: false,
            allow_nicked: false,
        },
        "anneal" => Operation::AnnealOligos {
            top: source.clone(),
            bottom: source,
            options: AnnealOptions {
                top_chemistry: collagene_core::chemistry::Terminus::UNKNOWN,
                bottom_chemistry: collagene_core::chemistry::Terminus::UNKNOWN,
                ..Default::default()
            },
            candidate_index: 0,
        },
        _ => Operation::Pcr {
            template: source,
            forward: Primer::new("Forward", "N", "").expect("literal primer"),
            reverse: Primer::new("Reverse", "N", "").expect("literal primer"),
        },
    }
}
#[derive(Clone)]
struct EditFrame {
    generation: u64,
    index: usize,
    insert: bool,
    operation: Operation,
}

#[component]
pub fn PlanEditingView(
    candidate: Candidate,
    epoch: RwSignal<u64>,
    input_sources: Memo<Vec<SourceInfo>>,
    capture_input: Callback<(String, EventBinding), Result<PlanInput, String>>,
    on_apply: Callback<DesignPlan, Result<(), String>>,
    on_close: Callback<()>,
    #[prop(optional)] bound_input_ids: Vec<String>,
) -> impl IntoView {
    let opened_epoch = epoch.get_untracked();
    let mut base = candidate.plan;
    let inputs = RwSignal::new(
        std::mem::take(&mut base.inputs)
            .into_iter()
            .map(Arc::new)
            .collect::<Vec<_>>(),
    );
    let steps = RwSignal::new(std::mem::take(&mut base.steps));
    let name = RwSignal::new(std::mem::take(&mut base.name));
    let target = base.target_revision.clone();
    let base = StoredValue::new(base);
    let ids = Memo::new(move |_| {
        inputs.with(|rows| rows.iter().map(|row| row.id.clone()).collect::<Vec<_>>())
    });
    let bound = RwSignal::new(bound_input_ids.into_iter().collect::<BTreeSet<_>>());
    let source_choices = RwSignal::new(BTreeMap::<String, String>::new());
    let new_id = RwSignal::new(String::new());
    let new_source = RwSignal::new(String::new());
    let selected = RwSignal::new(0usize);
    let page = RwSignal::new(0usize);
    let input_page = RwSignal::new(0usize);
    let kind = RwSignal::new("pcr".to_string());
    let status = RwSignal::new(String::new());
    let frame = RwSignal::new(None::<Arc<EditFrame>>);
    let generation = RwSignal::new(0u64);
    let stale = Memo::new(move |_| epoch.get() != opened_epoch);
    let input_labels = Memo::new(move |_| {
        inputs.with(|rows| {
            rows.iter()
                .map(|row| {
                    (
                        source_key(&Source::Input(row.id.clone())),
                        format!(
                            "{} · {} · {}…",
                            row.id,
                            row.molecule.name.chars().take(80).collect::<String>(),
                            short(&row.revision)
                        ),
                    )
                })
                .collect::<Vec<_>>()
        })
    });
    let start_edit = Callback::new(move |(index, insert): (usize, bool)| {
        if stale.get_untracked() || frame.with_untracked(Option::is_some) {
            return;
        }
        let operation = if insert {
            let Some(source) = first_source(&ids.get_untracked(), index) else {
                status.set("Add an input before inserting the first step.".into());
                return;
            };
            new_operation(&kind.get_untracked(), source)
        } else {
            let Some(op) = steps.with_untracked(|rows| rows.get(index).cloned()) else {
                return;
            };
            op
        };
        generation.update(|g| *g = g.wrapping_add(1));
        frame.set(Some(Arc::new(EditFrame {
            generation: generation.get_untracked(),
            index,
            insert,
            operation,
        })));
        status.set(String::new());
    });
    let move_step = Callback::new(move |to: usize| {
        if stale.get_untracked() || frame.with_untracked(Option::is_some) {
            return;
        }
        let from = selected.get_untracked();
        let result = steps
            .with_untracked(|rows| plan_editing::moved(&ids.get_untracked(), rows, from, to))
            .map_err(error);
        match result {
            Ok(rows) => {
                steps.set(rows);
                selected.set(to);
                page.set(to / 10);
                status.set("Step order changed; source dependencies were preserved. Verify after applying.".into());
            }
            Err(e) => status.set(e),
        }
    });
    let remove = move |_| {
        if stale.get_untracked() || frame.with_untracked(Option::is_some) {
            return;
        }
        let index = selected.get_untracked();
        match steps
            .with_untracked(|rows| plan_editing::removed(&ids.get_untracked(), rows, index))
            .map_err(error)
        {
            Ok(rows) => {
                let next = index.min(rows.len().saturating_sub(1));
                steps.set(rows);
                selected.set(next);
                page.set(next / 10);
                status.set(
                    "Step removed. Remaining dependencies were preserved; verify after applying."
                        .into(),
                );
            }
            Err(e) => status.set(e),
        }
    };
    let add_input = move |_| {
        if stale.get_untracked() {
            return;
        }
        let id = new_id.get_untracked().trim().to_string();
        if id.is_empty()
            || id.len() > 128
            || ids.with_untracked(|ids| ids.contains(&id))
            || inputs.with_untracked(|rows| rows.len() >= 1000)
        {
            status.set("Choose a new nonempty input ID of at most 128 bytes, within the 1,000-input limit.".into());
            return;
        }
        let key = new_source.get_untracked();
        let Some(source) = input_sources.with_untracked(|sources| {
            sources
                .iter()
                .find(|s| binding_key(&s.binding) == key)
                .cloned()
        }) else {
            status.set("Choose a current source version again.".into());
            return;
        };
        match capture_input.run((id.clone(), source.binding)) {
            Ok(input) if input.id == id => {
                inputs.update(|rows| rows.push(Arc::new(input)));
                bound.update(|ids| {
                    ids.insert(id);
                });
                new_id.set(String::new());
                status.set(
                    "Input version pinned. It is now available in step source selectors.".into(),
                );
            }
            Ok(_) => status.set("Captured input identity differs from the requested ID.".into()),
            Err(e) => status.set(e),
        }
    };
    let bind_input = Callback::new(move |id: String| {
        if stale.get_untracked() {
            return;
        }
        let expected = inputs.with_untracked(|rows| {
            rows.iter()
                .find(|row| row.id == id)
                .map(|row| row.revision.clone())
        });
        let Some(expected) = expected else {
            return;
        };
        let key = source_choices
            .with_untracked(|choices| choices.get(&id).cloned())
            .unwrap_or_default();
        let source = input_sources.with_untracked(|sources| {
            sources
                .iter()
                .find(|s| s.binding.revision == expected && binding_key(&s.binding) == key)
                .cloned()
        });
        let Some(source) = source else {
            status.set("Choose the exact matching input version again.".into());
            return;
        };
        match capture_input.run((id.clone(), source.binding)) {
            Ok(input) if input.id == id && input.revision == expected => {
                inputs.update(|rows| {
                    if let Some(row) = rows.iter_mut().find(|row| row.id == id) {
                        *row = Arc::new(input);
                    }
                });
                bound.update(|ids| {
                    ids.insert(id);
                });
                status.set("Exact input source version pinned.".into());
            }
            Ok(_) => status.set("Captured input differs from the original candidate input.".into()),
            Err(e) => status.set(e),
        }
    });
    let all_bound = Memo::new(move |_| {
        inputs.with(|rows| bound.with(|bound| rows.iter().all(|row| bound.contains(&row.id))))
    });
    let apply = move |_| {
        if stale.get_untracked() || frame.with_untracked(Option::is_some) {
            return;
        }
        if steps.with_untracked(|rows| rows.is_empty()) {
            status.set("Add at least one step before applying the candidate.".into());
            return;
        }
        if !all_bound.get_untracked() {
            status.set(
                "Bind every original input to an exact source version before applying.".into(),
            );
            return;
        }
        let name = name.get_untracked();
        if name.trim().is_empty() || name.len() > 1024 {
            status.set("Enter a plan name of 1–1,024 bytes.".into());
            return;
        }
        if let Err(e) =
            steps.with_untracked(|rows| plan_editing::validate_steps(&ids.get_untracked(), rows))
        {
            status.set(e.to_string());
            return;
        }
        let mut plan = base.get_value();
        plan.name = name;
        plan.steps = steps.get_untracked();
        plan.inputs =
            inputs.with_untracked(|rows| rows.iter().map(|row| row.as_ref().clone()).collect());
        match on_apply.run(plan) {
            Ok(()) => {}
            Err(e) => status.set(e),
        }
    };
    let editing_rows = move || frame.get().into_iter().collect::<Vec<_>>();
    let editing_key = |frame: &Arc<EditFrame>| frame.generation;
    let editing_child = move |edit: Arc<EditFrame>| {
        let limit = edit.index;
        let at = edit.index;
        let inserting = edit.insert;
        let choices = Memo::new(move |_| {
            let mut choices = input_labels.get();
            steps.with(|rows| {
                choices.extend(rows.iter().take(limit).enumerate().map(|(i, op)| {
                    (
                        source_key(&Source::Step(i)),
                        format!("Step {} product · {}", i + 1, title(op)),
                    )
                }))
            });
            choices
        });
        view! {<OperationEditor edit choices stale on_cancel=Callback::new(move|()|frame.set(None)) on_commit=Callback::new(move|operation:Operation|{
            if stale.get_untracked(){return Err("Workspace changed; reopen the plan editor.".into());}
            let result=steps.with_untracked(|rows|if inserting{plan_editing::inserted(&ids.get_untracked(),rows,at,operation)}else{plan_editing::replaced(&ids.get_untracked(),rows,at,operation)}).map_err(error)?;
            steps.set(result);selected.set(at);page.set(at/10);frame.set(None);status.set("Step declaration updated. Apply the plan and verify again.".into());Ok(())
        })/>}
    };
    view! {<section class="plan-editing-view">
        <style>".plan-editing-view .plan-edit-actions{display:flex;gap:6px;flex-wrap:wrap;margin:8px 0}.plan-editing-view .plan-edit-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,240px),1fr));gap:10px}.plan-editing-view label{display:block;margin:8px 0}.plan-editing-view input,.plan-editing-view select{max-width:100%;box-sizing:border-box}.plan-editing-view pre{white-space:pre-wrap;overflow-wrap:anywhere;max-height:260px;overflow:auto}.plan-editing-view .plan-edit-row{border-bottom:1px solid #dce3e8;padding:8px 0}.plan-editing-view .plan-edit-row.active{border-left:3px solid #2b8775;padding-left:8px}.plan-editing-view fieldset{min-width:0;border:1px solid #dce3e8;padding:10px}.plan-editing-view .plan-edit-source select{width:100%}"</style>
        <h2>"Edit cloning plan"</h2><p class="helper">"Edit declarations, then apply and verify the complete candidate. Existing scientific results are not carried over to the edited plan."</p>
        <p class="helper">{format!("Independent target remains fixed · {}…",short(&target))}</p>
        <Show when=move||stale.get()><p role="alert">"The workspace changed. Close and reopen this editor before applying changes."</p></Show>
        <fieldset disabled=move||stale.get()>
        <label>"Plan name"<input aria-label="Edited plan name" maxlength="1024" prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label>
        <details><summary>"Pinned input versions"</summary>
            <For each=move||inputs.with(|rows|rows.iter().skip(input_page.get()*10).take(10).map(|row|row.id.clone()).collect::<Vec<_>>()) key=|id|id.clone() children=move|id|{
                let row_id=StoredValue::new(id);view!{<div class="plan-edit-row">
                    <p>{move||inputs.with(|rows|rows.iter().find(|row|row.id==row_id.get_value()).map(|row|format!("{} · {} · {}…",row.id,row.molecule.name,short(&row.revision))).unwrap_or_default())}</p>
                    <Show when=move||bound.with(|ids|ids.contains(&row_id.get_value())) fallback=move||view!{
                        <label>"Choose matching source version"<select aria-label="Bind candidate input source" prop:value=move||source_choices.with(|choices|choices.get(&row_id.get_value()).cloned().unwrap_or_default()) on:change=move|e|source_choices.update(|choices|{choices.insert(row_id.get_value(),event_target_value(&e));})><option value="">"Choose exact source…"</option>{move||{
                            let revision=inputs.with(|rows|rows.iter().find(|row|row.id==row_id.get_value()).map(|row|row.revision.clone()));
                            input_sources.with(|sources|sources.iter().filter(|source|Some(&source.binding.revision)==revision.as_ref()).map(|source|view!{<option value=binding_key(&source.binding)>{format!("{} · {} · {}",source.name,source.binding.document_id,short(&source.binding.event_id))}</option>}).collect_view())
                        }}</select></label><button class="button secondary" on:click=move|_|bind_input.run(row_id.get_value())>"Bind source version"</button>
                    }><p class="helper">"Exact source version is pinned"</p></Show>
                </div>}
            }/>
            <div class="plan-edit-actions"><button class="text-tool" disabled=move||input_page.get()==0 on:click=move|_|input_page.update(|p|*p=p.saturating_sub(1))>"Previous inputs"</button><button class="text-tool" disabled=move||{(input_page.get()+1)*10>=inputs.with(|rows|rows.len())} on:click=move|_|input_page.update(|p|*p+=1)>"Next inputs"</button></div>
        </details>
        <details><summary>"Add source input"</summary><p class="helper">"Choose the original material and a new input ID. The independent target is not source material."</p><div class="plan-edit-grid"><label>"New input ID"<input aria-label="New plan input ID" maxlength="128" prop:value=move||new_id.get() on:input=move|e|new_id.set(event_target_value(&e))/></label><label>"Source document version"<select aria-label="New plan input source" prop:value=move||new_source.get() on:change=move|e|new_source.set(event_target_value(&e))><option value="">"Choose source…"</option>{move||input_sources.with(|sources|sources.iter().map(|source|view!{<option value=binding_key(&source.binding)>{format!("{} · {}…",source.name,short(&source.binding.revision))}</option>}).collect_view())}</select></label></div><button class="button secondary" on:click=add_input>"Add input version"</button></details>
        <h3>"Ordered steps"</h3>
        <For each=move||steps.with(|rows|(page.get()*10..rows.len().min((page.get()+1)*10)).collect::<Vec<_>>()) key=|i|*i children=move|i|view!{<div class="plan-edit-row" class:active=move||selected.get()==i><button class="text-tool" disabled=move||frame.with(Option::is_some) on:click=move|_|selected.set(i)>{move||steps.with(|rows|rows.get(i).map(|op|format!("{}. {}",i+1,title(op))).unwrap_or_default())}</button></div>}/>
        <div class="plan-edit-actions"><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous steps"</button><button class="text-tool" disabled=move||{(page.get()+1)*10>=steps.with(|rows|rows.len())} on:click=move|_|page.update(|p|*p+=1)>"Next steps"</button></div>
        <div class="plan-edit-actions"><button class="button secondary" disabled=move||frame.with(Option::is_some)||steps.with(|rows|rows.is_empty()) on:click=move|_|start_edit.run((selected.get_untracked(),false))>"Edit selected step"</button><button class="text-tool" disabled=move||frame.with(Option::is_some)||selected.get()==0 on:click=move|_|move_step.run(selected.get_untracked().saturating_sub(1))>"Move earlier"</button><button class="text-tool" disabled=move||frame.with(Option::is_some)||{selected.get()+1>=steps.with(|rows|rows.len())} on:click=move|_|move_step.run(selected.get_untracked()+1)>"Move later"</button><button class="text-tool" disabled=move||frame.with(Option::is_some)||steps.with(|rows|rows.is_empty()) on:click=remove>"Remove selected step"</button></div>
        <label>"New operation"<select aria-label="New plan operation" prop:value=move||kind.get() on:change=move|e|kind.set(event_target_value(&e))>{[("pcr","PCR"),("digest","Restriction digest"),("resect-three","Explicit 3′ source resection"),("anneal","Oligo annealing"),("treatment","End treatment"),("repair","End repair"),("tail","3′ tailing"),("ligate","Ordered ligation")].into_iter().map(|(value,label)|view!{<option value=value>{label}</option>}).collect_view()}</select></label>
        <div class="plan-edit-actions"><button class="button secondary" disabled=move||frame.with(Option::is_some)||steps.with(|rows|rows.len()>=1000) on:click=move|_|start_edit.run((selected.get_untracked().min(steps.with_untracked(|rows|rows.len())),true))>"Insert before selected step"</button><button class="button secondary" disabled=move||frame.with(Option::is_some)||steps.with(|rows|rows.len()>=1000) on:click=move|_|start_edit.run((steps.with_untracked(|rows|rows.len()),true))>"Append step"</button></div>
        <For each=editing_rows key=editing_key children=editing_child/>
        <p class="helper">"Moves preserve dependencies and are refused if a source would follow its consumer. A referenced step must be disconnected before removal."</p>
        <div class="plan-edit-actions"><button class="button primary" disabled=move||frame.with(Option::is_some)||!all_bound.get()||steps.with(|rows|rows.is_empty()) on:click=apply>"Apply edited candidate"</button></div>
        </fieldset><button class="button secondary" on:click=move|_|on_close.run(())>"Cancel plan editing"</button><p role="status">{move||status.get()}</p>
    </section>}
}

#[derive(Clone, Default)]
struct Fields {
    forward_name: String,
    forward_domain: String,
    forward_tail: String,
    reverse_name: String,
    reverse_domain: String,
    reverse_tail: String,
    enzymes: String,
    index: String,
    mode: String,
    circular: bool,
    nicked: bool,
    min_pair: String,
    top_five: String,
    top_three: String,
    bottom_five: String,
    bottom_three: String,
    resection: crate::plan_processing_model::ResectionDraft,
}
impl Fields {
    fn capture(op: &Operation, new: bool) -> Self {
        let mut f = Self::default();
        match op {
            Operation::Pcr {
                forward, reverse, ..
            }
            | Operation::PcrWithMismatches {
                forward, reverse, ..
            } => {
                f.forward_name = forward.name.clone();
                f.forward_domain = if new {
                    String::new()
                } else {
                    forward.annealing.as_str().into()
                };
                f.forward_tail = forward.tail.clone();
                f.reverse_name = reverse.name.clone();
                f.reverse_domain = if new {
                    String::new()
                } else {
                    reverse.annealing.as_str().into()
                };
                f.reverse_tail = reverse.tail.clone();
            }
            Operation::Digest {
                enzymes,
                fragment_index,
                ..
            } => {
                f.enzymes = enzymes.join(", ");
                f.index = fragment_index.to_string();
            }
            Operation::TreatEnds { treatment, .. } => {
                f.mode = if *treatment == EndTreatment::Phosphorylate {
                    "phosphorylate"
                } else {
                    "dephosphorylate"
                }
                .into()
            }
            Operation::RepairEnds { method, .. } => {
                f.mode = match method {
                    EndRepairMethod::FillFivePrimeOverhangs => "fill",
                    EndRepairMethod::RemoveThreePrimeOverhangs => "remove",
                    EndRepairMethod::Blunt => "blunt",
                }
                .into()
            }
            Operation::ResectThreePrime { cuts, .. } => {
                // Insertion starts unfinished; neither end is silently chosen.
                if !new {
                    f.resection = crate::plan_processing_model::ResectionDraft::capture(cuts);
                }
            }
            Operation::Tail3Prime { base, .. } => f.mode = base.as_str().into(),
            Operation::Ligate {
                circular,
                allow_nicked,
                ..
            } => {
                f.circular = *circular;
                f.nicked = *allow_nicked;
            }
            Operation::AnnealOligos {
                options,
                candidate_index,
                ..
            } => {
                f.min_pair = options.min_paired_bases.to_string();
                f.index = candidate_index.to_string();
                f.top_five = five_key(options.top_chemistry.five_prime).into();
                f.top_three = three_key(options.top_chemistry.three_prime).into();
                f.bottom_five = five_key(options.bottom_chemistry.five_prime).into();
                f.bottom_three = three_key(options.bottom_chemistry.three_prime).into();
            }
            _ => {}
        }
        f
    }
}
fn edited_primer(name: &str, domain: &str, tail: &str) -> Result<Primer, String> {
    let p = crate::primer_pair_model::OligoDraft {
        name: name.into(),
        annealing: domain.into(),
        tail: tail.into(),
        anchor: None,
    }
    .primer()?;
    if p.annealing.is_empty() {
        return Err(
            "Enter a nonempty annealing domain for each primer before inserting the PCR step."
                .into(),
        );
    }
    Ok(p)
}
fn five_key(value: FivePrime) -> &'static str {
    match value {
        FivePrime::Phosphate => "P",
        FivePrime::Hydroxyl => "OH",
        FivePrime::Unknown => "unknown",
    }
}
fn three_key(value: ThreePrime) -> &'static str {
    match value {
        ThreePrime::Hydroxyl => "OH",
        ThreePrime::Blocked => "blocked",
        ThreePrime::Unknown => "unknown",
    }
}
fn five(value: &str) -> Result<FivePrime, String> {
    match value {
        "P" => Ok(FivePrime::Phosphate),
        "OH" => Ok(FivePrime::Hydroxyl),
        "unknown" => Ok(FivePrime::Unknown),
        _ => Err("Choose explicit 5′ chemistry.".into()),
    }
}
fn three(value: &str) -> Result<ThreePrime, String> {
    match value {
        "OH" => Ok(ThreePrime::Hydroxyl),
        "blocked" => Ok(ThreePrime::Blocked),
        "unknown" => Ok(ThreePrime::Unknown),
        _ => Err("Choose explicit 3′ chemistry.".into()),
    }
}
fn text_field(
    fields: RwSignal<Fields>,
    label: &'static str,
    get: fn(&Fields) -> &str,
    set: fn(&mut Fields, String),
) -> AnyView {
    view!{<label>{label}<input aria-label=label maxlength="1600" prop:value=move||fields.with(|f|get(f).to_owned()) on:input=move|e|fields.update(|f|set(f,event_target_value(&e)))/></label>}.into_any()
}
fn choice_field(
    fields: RwSignal<Fields>,
    label: &'static str,
    get: fn(&Fields) -> &str,
    set: fn(&mut Fields, String),
    options: Vec<(&'static str, &'static str)>,
) -> AnyView {
    view!{<label>{label}<select aria-label=label prop:value=move||fields.with(|f|get(f).to_owned()) on:change=move|e|fields.update(|f|set(f,event_target_value(&e)))>{options.into_iter().map(|(value,label)|view!{<option value=value>{label}</option>}).collect_view()}</select></label>}.into_any()
}
#[component]
fn OperationEditor(
    edit: Arc<EditFrame>,
    choices: Memo<Vec<(String, String)>>,
    stale: Memo<bool>,
    on_commit: Callback<Operation, Result<(), String>>,
    on_cancel: Callback<()>,
) -> impl IntoView {
    let original = StoredValue::new(edit.clone());
    let fields = RwSignal::new(Fields::capture(&edit.operation, edit.insert));
    let sources = RwSignal::new(
        verification::operation_sources(&edit.operation)
            .into_iter()
            .cloned()
            .collect::<Vec<_>>(),
    );
    let reverse = RwSignal::new(match &edit.operation {
        Operation::Ligate { fragments, .. } => {
            fragments.iter().map(|f| f.reverse).collect::<Vec<_>>()
        }
        _ => vec![],
    });
    let source_page = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let is_ligation = matches!(edit.operation, Operation::Ligate { .. });
    let params=match &edit.operation {
        Operation::Pcr{..}|Operation::PcrWithMismatches{..}=>view!{<div class="plan-edit-grid">
            {text_field(fields,"Forward primer name",|f|&f.forward_name,|f,v|f.forward_name=v)}{text_field(fields,"Forward annealing domain (5′→3′)",|f|&f.forward_domain,|f,v|f.forward_domain=v)}{text_field(fields,"Forward 5′ tail",|f|&f.forward_tail,|f,v|f.forward_tail=v)}
            {text_field(fields,"Reverse primer name",|f|&f.reverse_name,|f,v|f.reverse_name=v)}{text_field(fields,"Reverse annealing domain (5′→3′)",|f|&f.reverse_domain,|f,v|f.reverse_domain=v)}{text_field(fields,"Reverse 5′ tail",|f|&f.reverse_tail,|f,v|f.reverse_tail=v)}
        </div>}.into_any(),
        Operation::ReconstructMaterial{origin,..}=>view!{<div class="helper"><p>{format!("Exact reusable material event {}",origin.product.event_id)}</p><p>"Its retained target-free PCR/digest derivation is replayed. Imported reference summaries confer no material authority."</p></div>}.into_any(),
        Operation::ReconstructDigestFragment{origin,..}=>view!{<div class="helper"><p>{format!("Exact saved digest event {} · originating run {}",origin.record.product.event_id,origin.record.request.origin.selection.run_id)}</p><p>"The original run and selected digest are replayed. This editor changes only the outer dependency; documentary end summaries do not confer material authority."</p></div>}.into_any(),
        Operation::ReconstructRunProduct{origin,..}=>view!{<div class="helper"><p>{format!("Immutable originating run {} · step {} · selected source event {}",origin.selection.run_id,origin.selection.step_index+1,origin.selection.event.event_id)}</p><p>{format!("Original independent target: {} · event {}",origin.target.name,origin.target_event.event_id)}</p><p>"The retained origin plan is reconstructed against its own target. Origin parameters and historical selections are preserved; this editor changes only the outer dependency. A different event cannot acquire material authority through matching DNA or imported summaries. Historical primers are not added to a new order."</p></div>}.into_any(),
        Operation::Digest{..}=>view!{<div>{text_field(fields,"Restriction enzymes (comma separated)",|f|&f.enzymes,|f,v|f.enzymes=v)}{text_field(fields,"Digest fragment index (zero based)",|f|&f.index,|f,v|f.index=v)}</div>}.into_any(),
        Operation::ResectThreePrime{..}=>view!{<div>
            {text_field(fields,"Left bottom-strand 3′ cut (bases)",|f|&f.resection.left,|f,v|f.resection.left=v)}
            {text_field(fields,"Right top-strand 3′ cut (bases)",|f|&f.resection.right,|f,v|f.resection.right=v)}
            <p class="helper">"Enter both cuts explicitly, including zero for an end left unchanged. Values shorten stored 3′ suffixes of the selected physical duplex; original 5′ groups remain. Choose an earlier eligible physical preparation (such as PCR or digest), or an explicitly declared source-termini step. Raw sequence letters alone do not establish accessible 3′ hydroxyls. Native verification rejects circular, nicked, blocked, unknown or exhausted substrates as applicable; no enzyme, unblocking, fill, sealing or host repair is inferred."</p>
        </div>}.into_any(),
        Operation::TreatEnds{..}=>choice_field(fields,"End treatment",|f|&f.mode,|f,v|f.mode=v,vec![("phosphorylate","Phosphorylate 5′ ends"),("dephosphorylate","Dephosphorylate 5′ ends")]),
        Operation::RepairEnds{..}=>choice_field(fields,"End repair",|f|&f.mode,|f,v|f.mode=v,vec![("fill","Fill 5′ overhangs"),("remove","Remove 3′ overhangs"),("blunt","Blunt ends")]),
        Operation::Tail3Prime{..}=>choice_field(fields,"Added 3′ base",|f|&f.mode,|f,v|f.mode=v,vec![("A","A"),("C","C"),("G","G"),("T","T")]),
        Operation::Ligate{..}=>view!{<div><label class="checkbox-label"><input type="checkbox" prop:checked=move||fields.with(|f|f.circular) on:change=move|e|fields.update(|f|f.circular=event_target_checked(&e))/>"Circular ligation"</label><label class="checkbox-label"><input type="checkbox" prop:checked=move||fields.with(|f|f.nicked) on:change=move|e|fields.update(|f|f.nicked=event_target_checked(&e))/>"Allow retained nicks"</label></div>}.into_any(),
        Operation::AnnealOligos{..}=>view!{<div class="plan-edit-grid">
            {text_field(fields,"Minimum paired bases",|f|&f.min_pair,|f,v|f.min_pair=v)}{text_field(fields,"Annealing candidate index (zero based)",|f|&f.index,|f,v|f.index=v)}
            {choice_field(fields,"Top oligo 5′ chemistry",|f|&f.top_five,|f,v|f.top_five=v,vec![("unknown","Unknown"),("P","Phosphate"),("OH","Hydroxyl")])}
            {choice_field(fields,"Top oligo 3′ chemistry",|f|&f.top_three,|f,v|f.top_three=v,vec![("unknown","Unknown"),("OH","Hydroxyl"),("blocked","Blocked")])}
            {choice_field(fields,"Bottom oligo 5′ chemistry",|f|&f.bottom_five,|f,v|f.bottom_five=v,vec![("unknown","Unknown"),("P","Phosphate"),("OH","Hydroxyl")])}
            {choice_field(fields,"Bottom oligo 3′ chemistry",|f|&f.bottom_three,|f,v|f.bottom_three=v,vec![("unknown","Unknown"),("OH","Hydroxyl"),("blocked","Blocked")])}
        </div>}.into_any(),
        _=>view!{<p class="helper">"Existing preparation parameters are preserved. Use the method's preparation panel for detailed parameter changes, then inspect the method again after rebinding a source."</p>}.into_any(),
    };
    let commit = move |_| {
        if stale.get_untracked() {
            return;
        }
        let result = (|| -> Result<Operation, String> {
            let f = fields.get_untracked();
            let refs = sources.get_untracked();
            let mut op = original.with_value(|e| e.operation.clone());
            if is_ligation {
                let flags = reverse.get_untracked();
                if refs.is_empty() || refs.len() > 1000 {
                    return Err("Choose 1–1,000 ordered ligation fragments.".into());
                }
                let fragments = refs
                    .into_iter()
                    .enumerate()
                    .map(|(i, source)| OrientedSource {
                        source,
                        reverse: flags.get(i).copied().unwrap_or(false),
                    })
                    .collect();
                let Operation::Ligate { options, .. } = &op else {
                    return Err("Expected a ligation operation.".into());
                };
                return Ok(Operation::Ligate {
                    fragments,
                    circular: f.circular,
                    allow_nicked: f.nicked,
                    options: options.clone(),
                });
            }
            let mut slots = verification::operation_sources_mut(&mut op);
            if slots.len() != refs.len() {
                return Err("Operation source count changed.".into());
            }
            for (slot, source) in slots.iter_mut().zip(refs) {
                **slot = source;
            }
            drop(slots);
            match &mut op {
                Operation::Pcr {
                    forward, reverse, ..
                }
                | Operation::PcrWithMismatches {
                    forward, reverse, ..
                } => {
                    *forward = edited_primer(&f.forward_name, &f.forward_domain, &f.forward_tail)?;
                    *reverse = edited_primer(&f.reverse_name, &f.reverse_domain, &f.reverse_tail)?;
                }
                Operation::Digest {
                    enzymes,
                    fragment_index,
                    ..
                } => {
                    *enzymes = f
                        .enzymes
                        .split(',')
                        .map(str::trim)
                        .filter(|s| !s.is_empty())
                        .map(str::to_owned)
                        .collect();
                    if enzymes.is_empty()
                        || enzymes.len() > 32
                        || enzymes.iter().any(|s| s.len() > 256)
                    {
                        return Err("Enter 1–32 enzyme names, separated by commas.".into());
                    }
                    *fragment_index = number(&f.index, "Digest fragment index")?;
                }
                Operation::TreatEnds { treatment, .. } => {
                    *treatment = match f.mode.as_str() {
                        "phosphorylate" => EndTreatment::Phosphorylate,
                        "dephosphorylate" => EndTreatment::Dephosphorylate,
                        _ => return Err("Choose an end treatment.".into()),
                    }
                }
                Operation::RepairEnds { method, .. } => {
                    *method = match f.mode.as_str() {
                        "fill" => EndRepairMethod::FillFivePrimeOverhangs,
                        "remove" => EndRepairMethod::RemoveThreePrimeOverhangs,
                        "blunt" => EndRepairMethod::Blunt,
                        _ => return Err("Choose an end-repair method.".into()),
                    }
                }
                Operation::Tail3Prime { base, .. } => {
                    *base = match f.mode.as_str() {
                        "A" => TailingBase::A,
                        "C" => TailingBase::C,
                        "G" => TailingBase::G,
                        "T" => TailingBase::T,
                        _ => return Err("Choose a single tailing base.".into()),
                    }
                }
                Operation::AnnealOligos {
                    options,
                    candidate_index,
                    ..
                } => {
                    options.min_paired_bases = number(&f.min_pair, "Minimum pairing")?;
                    if options.min_paired_bases == 0
                        || options.min_paired_bases
                            > collagene_core::chemistry::MAX_ANNEAL_OLIGO_BASES
                    {
                        return Err("Minimum pairing must be between 1 and 10,000 bases.".into());
                    }
                    *candidate_index = number(&f.index, "Annealing candidate index")?;
                    options.top_chemistry.five_prime = five(&f.top_five)?;
                    options.top_chemistry.three_prime = three(&f.top_three)?;
                    options.bottom_chemistry.five_prime = five(&f.bottom_five)?;
                    options.bottom_chemistry.three_prime = three(&f.bottom_three)?;
                }
                _ => {}
            }
            if matches!(&op, Operation::ResectThreePrime { .. }) {
                return f.resection.apply(&op);
            }
            Ok(op)
        })()
        .and_then(|op| on_commit.run(op));
        if let Err(e) = result {
            status.set(e);
        }
    };
    let source_labels = StoredValue::new(match &edit.operation {
        Operation::AnnealOligos { .. } => vec!["Top oligo".to_string(), "Bottom oligo".to_string()],
        Operation::Topo { sources, .. }
        | Operation::CheckedTopo { sources, .. }
        | Operation::ProcessedTopo { sources, .. }
        | Operation::TopoActions { sources, .. } => {
            sources.iter().map(|source| source.role.clone()).collect()
        }
        Operation::ReconstructMaterial { .. }
        | Operation::ReconstructDigestFragment { .. }
        | Operation::ReconstructRunProduct { .. }
        | Operation::Pcr { .. }
        | Operation::PcrWithMismatches { .. }
        | Operation::Digest { .. }
        | Operation::CatalogDigest { .. }
        | Operation::CheckedDigest { .. }
        | Operation::TreatEnds { .. }
        | Operation::RepairEnds { .. }
        | Operation::Tail3Prime { .. }
        | Operation::ResectThreePrime { .. }
        | Operation::DeclareSourceModifications { .. }
        | Operation::DeclareSourceTermini { .. } => vec!["Template".to_string()],
        _ => vec![],
    });
    let detail = crate::app::describe_step(&edit.operation);
    let truncated = detail.chars().count() > 50_000;
    let detail = detail.chars().take(50_000).collect::<String>();
    view! {<fieldset class="plan-edit-operation" disabled=move||stale.get()><legend>{format!("{} step {} · {}",if edit.insert{"Insert"}else{"Edit"},edit.index+1,title(&edit.operation))}</legend>
        <h4>"Sources"</h4><p class="helper">"Select exact bound inputs or products of earlier steps. Source selection does not establish compatible chemistry."</p>
        <For each=move||sources.with(|rows|(source_page.get()*6..rows.len().min((source_page.get()+1)*6)).collect::<Vec<_>>()) key=|i|*i children=move|i|view!{<div class="plan-edit-source">
            <label>{source_labels.with_value(|labels|labels.get(i).cloned().unwrap_or_else(||format!("Fragment {}",i+1)))}<select aria-label=format!("Step source {}",i+1) prop:value=move||sources.with(|rows|rows.get(i).map(source_key).unwrap_or_default()) on:change=move|e|match parse_source(&event_target_value(&e)){Ok(source)=>sources.update(|rows|if let Some(row)=rows.get_mut(i){*row=source}),Err(e)=>status.set(e)}>{move||choices.with(|rows|rows.iter().map(|(value,label)|view!{<option value=value.clone()>{label.clone()}</option>}).collect_view())}</select></label>
            <Show when=move||is_ligation><div class="plan-edit-actions"><label class="checkbox-label"><input type="checkbox" prop:checked=move||reverse.with(|rows|rows.get(i).copied().unwrap_or(false)) on:change=move|e|reverse.update(|rows|if let Some(row)=rows.get_mut(i){*row=event_target_checked(&e)})/>"Reverse fragment"</label><button class="text-tool" disabled=i==0 on:click=move|_|{if i>0{sources.update(|rows|rows.swap(i,i-1));reverse.update(|rows|rows.swap(i,i-1));}}>"Earlier fragment"</button><button class="text-tool" disabled=move||{i+1>=sources.with(|rows|rows.len())} on:click=move|_|{if i+1<sources.with_untracked(|rows|rows.len()){sources.update(|rows|rows.swap(i,i+1));reverse.update(|rows|rows.swap(i,i+1));}}>"Later fragment"</button><button class="text-tool" disabled=move||sources.with(|rows|rows.len()<=1) on:click=move|_|{sources.update(|rows|{rows.remove(i);});reverse.update(|rows|{rows.remove(i);});source_page.update(|page|*page=(*page).min(sources.with_untracked(|rows|rows.len().saturating_sub(1)/6)));}>"Remove fragment"</button></div></Show>
        </div>}/>
        <div class="plan-edit-actions"><button class="text-tool" disabled=move||source_page.get()==0 on:click=move|_|source_page.update(|p|*p=p.saturating_sub(1))>"Previous source references"</button><button class="text-tool" disabled=move||{(source_page.get()+1)*6>=sources.with(|rows|rows.len())} on:click=move|_|source_page.update(|p|*p+=1)>"Next source references"</button><Show when=move||is_ligation><button class="button secondary" disabled=move||sources.with(|rows|rows.len()>=1000) on:click=move|_|{if let Some((key,_))=choices.with_untracked(|rows|rows.first().cloned()){if let Ok(source)=parse_source(&key){sources.update(|rows|rows.push(source));reverse.update(|rows|rows.push(false));source_page.set(sources.with_untracked(|rows|rows.len().saturating_sub(1)/6));}}}>"Add ligation fragment"</button></Show></div>
        {params}<details><summary>"Existing operation details"</summary><pre>{detail}</pre>{truncated.then(||view!{<p class="helper">"Long preparation details are abbreviated here. Use the existing method panel for its full preparation controls."</p>})}</details>
        <div class="plan-edit-actions"><button class="button primary" on:click=commit>{if edit.insert{"Insert step declaration"}else{"Replace step declaration"}}</button><button class="button secondary" on:click=move|_|on_cancel.run(())>"Discard step edits"</button></div><p role="status">{move||status.get()}</p>
    </fieldset>}
}
