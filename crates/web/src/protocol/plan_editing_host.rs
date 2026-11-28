//! Explicit candidate editing preserves document-event choices, never chemistry.
use crate::protocol::plan_editing_sources;
use crate::protocol::plan_editing_view::PlanEditingView;
use collagene_core::{
    documents::{history::DocumentHistory, transfer::Candidate},
    material::lineage::{EventBinding, InputEventBinding},
    planning::plan::{DesignPlan, PlanInput},
};
use collagene_model::{Alphabet, bounds::model_preflight};
use leptos::prelude::*;
use std::{collections::BTreeMap, sync::Arc};

pub type Preferences = (String, u64, EventBinding, Vec<InputEventBinding>);
struct Frame {
    original: Candidate,
    expected_candidate: Option<String>,
    target: EventBinding,
    epoch: u64,
}
fn event(doc: &DocumentHistory) -> EventBinding {
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event().into(),
        revision: doc.head_revision().into(),
    }
}
fn retained(docs: &[DocumentHistory], binding: &EventBinding) -> bool {
    plan_editing_sources::resolve(docs, binding).is_ok()
}
#[component]
pub fn ManualPlanEditor(
    candidate: RwSignal<Option<Candidate>>,
    #[prop(optional)] creating: bool,
    preferences: RwSignal<Option<Preferences>>,
    documents: RwSignal<Vec<DocumentHistory>>,
    active: RwSignal<usize>,
    epoch: RwSignal<u64>,
    locked: Memo<bool>,
    on_apply: Callback<(Candidate, EventBinding, Vec<InputEventBinding>), Result<(), String>>,
    on_close: Callback<()>,
) -> impl IntoView {
    let prepared = (|| -> Result<(Frame, BTreeMap<String, EventBinding>), String> {
        if locked.get_untracked() {
            return Err("Wait for the workspace update to finish.".into());
        }
        let initial = candidate.get_untracked();
        let expected_candidate = initial.as_ref().map(|c| c.plan_revision.clone());
        let pinned = if creating {
            None
        } else {
            preferences
                .get_untracked()
                .filter(|p| initial.as_ref().is_some_and(|c| p.0 == c.plan_revision))
        };
        if pinned
            .as_ref()
            .is_some_and(|p| p.1 != epoch.get_untracked())
        {
            return Err("The candidate edit scope changed; prepare it again.".into());
        }
        let target = match &pinned {
            Some(p) => p.2.clone(),
            None => documents
                .with_untracked(|docs| docs.get(active.get_untracked()).map(event))
                .ok_or("Select the candidate's independent target.")?,
        };
        let original = if creating {
            Candidate::new(DesignPlan {
                standard: None,
                name: "Manual cloning plan".into(),
                target_revision: target.revision.clone(),
                inputs: vec![],
                steps: vec![],
            })
        } else {
            let value = initial.ok_or("Open a candidate before editing its steps.")?;
            value.validate().map_err(|e| e.to_string())?;
            value
        };
        documents.with_untracked(|docs| {
            plan_editing_sources::resolve_target(docs, &target, &original.plan.target_revision)
                .map(|_| ())
        })?;
        let mut bindings = BTreeMap::new();
        if let Some((revision, version, preferred_target, preferred)) = preferences
            .get_untracked()
            .filter(|p| p.0 == original.plan_revision)
        {
            if revision != original.plan_revision
                || version != epoch.get_untracked()
                || preferred_target != target
            {
                return Err("The candidate's pinned target changed. Reopen it against the original target before editing.".into());
            }
            documents.with_untracked(|docs| {
                for binding in preferred {
                    if original
                        .plan
                        .inputs
                        .iter()
                        .any(|i| i.id == binding.input_id && i.revision == binding.event.revision)
                        && retained(docs, &binding.event)
                    {
                        bindings.insert(binding.input_id, binding.event);
                    }
                }
            });
        }
        Ok((
            Frame {
                original,
                expected_candidate,
                target,
                epoch: epoch.get_untracked(),
            },
            bindings,
        ))
    })();
    let (frame,initial)=match prepared {
        Ok(value)=>value,
        Err(error)=>return view!{<p role="alert">{error}</p><button class="button secondary" on:click=move|_|on_close.run(())>"Close editor"</button>}.into_any(),
    };
    let bound_input_ids = initial.keys().cloned().collect::<Vec<_>>();
    let frame = StoredValue::new(Arc::new(frame));
    let bindings = RwSignal::new(initial);
    let guard = Callback::new(move |()| -> Result<(), String> {
        if locked.get_untracked() || epoch.get_untracked() != frame.with_value(|f| f.epoch) {
            return Err("The workspace changed; reopen this candidate editor.".into());
        }
        if candidate.with_untracked(|c| c.as_ref().map(|c| c.plan_revision.clone()))
            != frame.with_value(|f| f.expected_candidate.clone())
        {
            return Err("The active candidate changed; reopen the editor.".into());
        }
        documents.with_untracked(|docs| {
            frame.with_value(|f| {
                plan_editing_sources::resolve_target(
                    docs,
                    &f.target,
                    &f.original.plan.target_revision,
                )
                .map(|_| ())
            })
        })?;
        Ok(())
    });
    let capture_input = Callback::new(
        move |(id, binding): (String, EventBinding)| -> Result<PlanInput, String> {
            guard.run(())?;
            if id.trim().is_empty() || id.len() > 128 || id.chars().any(char::is_control) {
                return Err("Choose an input ID of 1–128 characters.".into());
            }
            if binding == frame.with_value(|f| f.target.clone()) {
                return Err(
                    "Choose a source document separately from the independent target.".into(),
                );
            }
            if let Some(original) = frame.with_value(|f| {
                f.original
                    .plan
                    .inputs
                    .iter()
                    .find(|input| input.id == id)
                    .map(|i| i.revision.clone())
            }) {
                if binding.revision != original {
                    return Err("Existing inputs keep their original revision. Add a new input ID to use another source.".into());
                }
            }
            if bindings.with_untracked(|b| b.get(&id).is_some_and(|previous| *previous != binding))
            {
                return Err(
                    "This input ID already has an exact source event. Add a different input ID."
                        .into(),
                );
            }
            let input = documents.with_untracked(|docs| -> Result<PlanInput, String> {
                let molecule = plan_editing_sources::resolve(docs, &binding)?;
                if molecule.sequence.alphabet() != Alphabet::Dna {
                    return Err("Cloning inputs must be DNA.".into());
                }
                model_preflight(&[molecule]).map_err(|e| e.to_string())?;
                if molecule.sequence.len() > 32_000_000 {
                    return Err("The source exceeds the cloning-plan size limit.".into());
                }
                Ok(PlanInput {
                    id: id.clone(),
                    revision: binding.revision.clone(),
                    molecule: molecule.clone(),
                })
            })?;
            bindings.update(|b| {
                b.insert(id, binding);
            });
            Ok(input)
        },
    );
    let apply = Callback::new(move |plan: DesignPlan| -> Result<(), String> {
        guard.run(())?;
        if plan.target_revision != frame.with_value(|f| f.original.plan.target_revision.clone()) {
            return Err("The independent target cannot be changed by editing steps.".into());
        }
        let ids = plan.inputs.iter().map(|i| i.id.clone()).collect::<Vec<_>>();
        collagene_core::planning::editing::validate_steps(&ids, &plan.steps)
            .map_err(|e| e.to_string())?;
        if plan.steps.is_empty() {
            return Err("Add at least one cloning step before applying.".into());
        }
        model_preflight(&plan.inputs.iter().map(|i| &i.molecule).collect::<Vec<_>>())
            .map_err(|e| e.to_string())?;
        plan.inputs.iter().try_fold(0usize, |sum, input| {
            sum.checked_add(input.molecule.sequence.len())
                .filter(|n| *n <= 32_000_000)
                .ok_or("The combined inputs exceed the cloning-plan size limit.")
        })?;
        let preserved = frame.with_value(|f| {
            f.original.plan.inputs.iter().all(|old| {
                plan.inputs
                    .iter()
                    .any(|new| new.id == old.id && new.revision == old.revision)
            })
        });
        if !preserved {
            return Err(
                "Original input identities must remain unchanged; add another source explicitly."
                    .into(),
            );
        }
        let exact = bindings.with_untracked(|bindings| {
            plan.inputs
                .iter()
                .map(|input| {
                    let event = bindings
                        .get(&input.id)
                        .filter(|b| b.revision == input.revision)
                        .cloned()
                        .ok_or_else(|| {
                            format!("Choose the exact source event for input {}.", input.id)
                        })?;
                    Ok(InputEventBinding {
                        input_id: input.id.clone(),
                        event,
                    })
                })
                .collect::<Result<Vec<_>, String>>()
        })?;
        if !documents.with_untracked(|docs| exact.iter().all(|b| retained(docs, &b.event))) {
            return Err(
                "A selected source event is no longer available; reopen the editor.".into(),
            );
        }
        let value = Candidate::new(plan);
        value.validate().map_err(|e| e.to_string())?;
        guard.run(())?;
        on_apply.run((value, frame.with_value(|f| f.target.clone()), exact))
    });
    let required = StoredValue::new(frame.with_value(|f| {
        f.original
            .plan
            .inputs
            .iter()
            .map(|i| i.revision.clone())
            .collect::<std::collections::BTreeSet<_>>()
    }));
    let catalog = Memo::new(move |_| {
        documents.with(|docs| {
            required.with_value(|required| plan_editing_sources::catalog(docs, required, 4096))
        })
    });
    let input_sources = Memo::new(move |_| {
        catalog.with(|result| {
            result
                .as_ref()
                .map(|catalog| {
                    catalog
                        .choices
                        .iter()
                        .filter(|source| source.binding != frame.with_value(|f| f.target.clone()))
                        .cloned()
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default()
        })
    });
    let original = frame.with_value(|f| f.original.clone());
    let target_label=frame.with_value(|f|format!("Independent target: document {} · event {} · revision {}. Viewer navigation does not replace this target.",f.target.document_id,f.target.event_id,f.target.revision));
    view! {<p class="helper">{target_label}</p><p class="helper">"Choose an exact event for each unbound original input. Equal DNA in another document does not choose its ancestry."</p>
        <p class="helper" role="status">{move||catalog.with(|result|match result {Err(error)=>error.clone(),Ok(catalog) if catalog.omitted>0=>format!("{} additional historical source choices omitted. Activate the wanted event in its document history and reopen this editor to include it as a current source.",catalog.omitted),_=>String::new()})}</p>
        <PlanEditingView candidate=original epoch input_sources capture_input bound_input_ids on_apply=apply on_close/>}.into_any()
}
