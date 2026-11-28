//! Explicit source-bound site removal: reference constraints and physical plans stay separate.
use crate::cloning::methods::domestication_form::Form;
use crate::cloning::methods::domestication_form::RouteForm;
use crate::cloning::methods::domestication_review::Open;
use crate::history::run_origin_model::PinnedOrigin;
use crate::platform::worker::CachedWorker;
use crate::sources::identity::SourceInfo;
use collagene_core::{
    cloning::golden_gate::domestication::{self as design, preparation as physical},
    documents::workspace::Workspace,
    planning::plan::PlanInput,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::Molecule;
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

struct Scope {
    source: SourceInfo,
    molecule: Arc<Molecule>,
    epoch: u64,
    form: Form,
    request: design::Request,
}
struct Inspection {
    scope: Arc<Scope>,
    report: design::Inspection,
}
struct Proposed {
    scope: Arc<Scope>,
    report: design::Proposal,
}
struct Checked {
    scope: Arc<Scope>,
    selection: design::Selection,
    report: design::CheckReport,
    target: Molecule,
}
struct MaterialSource {
    origin: Option<PinnedOrigin>,
    workspace: Option<Arc<Workspace>>,
}
struct Prepared {
    material: Arc<MaterialSource>,
    checked: Arc<Checked>,
    route: RouteForm,
    report: physical::Prepared,
}
#[derive(Clone)]
enum Work {
    Inspect,
    Propose,
    Check(Arc<Proposed>),
    Prepare(Arc<Checked>, RouteForm),
}
const LIST_PAGE: usize = 20;
const SITE_PAGE: usize = 40;
const PEPTIDE_PAGE: usize = 600;
fn pages(total: usize, limit: usize) -> usize {
    total.max(1).div_ceil(limit)
}
fn feature_name(feature: &collagene_model::Feature, index: usize) -> String {
    feature
        .qualifiers
        .iter()
        .find(|(key, _)| matches!(key.as_str(), "label" | "gene" | "product"))
        .and_then(|(_, value)| value.as_deref())
        .filter(|value| !value.is_empty())
        .map(str::to_string)
        .unwrap_or_else(|| format!("CDS {}", index + 1))
}
fn source_window(source: &Molecule, position: usize) -> String {
    let dna = source.sequence.as_str();
    if position >= dna.len() {
        return "Source position unavailable".into();
    }
    let start = position.saturating_sub(12);
    let end = position.saturating_add(13).min(dna.len());
    format!(
        "Original source {}–{} (linearized display): {}[{}]{}",
        start + 1,
        end,
        &dna[start..position],
        &dna[position..position + 1],
        &dna[position + 1..end]
    )
}
#[component]
fn Pager(page: RwSignal<usize>, total: usize, limit: usize, label: &'static str) -> impl IntoView {
    let count = pages(total, limit);
    view! {<p>{move||format!("{} {} · page {} of {}",total,label,page.get().min(count-1)+1,count)}</p>
    <button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|value|*value=value.saturating_sub(1))>{format!("Previous {label}")}</button>
    <button class="text-tool" disabled={move||page.get()+1>=count} on:click=move|_|page.update(|value|*value=(*value+1).min(count-1))>{format!("Next {label}")}</button>}
}

#[component]
pub fn DomesticationView(
    source: Memo<SourceInfo>,
    molecule: Memo<Molecule>,
    selected_feature: RwSignal<Option<usize>>,
    epoch: RwSignal<u64>,
    get_workspace: Callback<(), Result<Workspace, String>>,
    on_open: Callback<Open, Result<(), String>>,
    on_navigate: Callback<(usize, usize), Result<(), String>>,
) -> impl IntoView {
    let form =
        RwSignal::new(molecule.with_untracked(|m| Form::new(m, selected_feature.get_untracked())));
    let route = RwSignal::new(RouteForm::default());
    let origin = RwSignal::new(None::<PinnedOrigin>);
    let inspection = RwSignal::new(None::<Arc<Inspection>>);
    let proposed = RwSignal::new(None::<Arc<Proposed>>);
    let checked = RwSignal::new(None::<Arc<Checked>>);
    let prepared = RwSignal::new(None::<Arc<Prepared>>);
    let pending = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let feature_page = RwSignal::new(0usize);
    let site_page = RwSignal::new(0usize);
    let diagnostic_page = RwSignal::new(0usize);
    let retained_page = RwSignal::new(0usize);
    let search_page = RwSignal::new(0usize);
    let inspection_note_page = RwSignal::new(0usize);
    let preparation_note_page = RwSignal::new(0usize);
    let patch_page = RwSignal::new(0usize);
    let check_page = RwSignal::new(0usize);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
    });
    let cancel = Callback::new(move |clear_sites: bool| {
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        pending.set(false);
        proposed.set(None);
        checked.set(None);
        prepared.set(None);
        patch_page.set(0);
        check_page.set(0);
        if clear_sites {
            inspection.set(None);
            site_page.set(0);
            diagnostic_page.set(0);
        }
    });
    Effect::new(move |_| {
        let _ = (source.get(), epoch.get());
        cancel.run(true);
        form.set(molecule.with_untracked(|m| Form::new(m, selected_feature.get_untracked())));
        route.set(RouteForm::default());
        origin.set(None);
        feature_page.set(0);
        retained_page.set(0);
        search_page.set(0);
        inspection_note_page.set(0);
        preparation_note_page.set(0);
    });
    Effect::new(move |_| {
        let current = form.get();
        cancel.run(false);
        if inspection.with_untracked(|saved| {
            saved
                .as_ref()
                .is_some_and(|saved| saved.scope.form.enzymes != current.enzymes)
        }) {
            inspection.set(None);
            site_page.set(0);
            diagnostic_page.set(0);
        }
    });
    Effect::new(move |_| {
        let _ = route.get();
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        pending.set(false);
        prepared.set(None);
    });
    Effect::new(move |_| {
        let _ = origin.get();
        generation.update(|value| *value = value.wrapping_add(1));
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        pending.set(false);
        prepared.set(None);
    });
    let attest_material=Callback::new(move |(saved,selected):(Arc<Scope>,Option<PinnedOrigin>)|->Result<Arc<MaterialSource>,String>{
        let workspace=if selected.is_some(){Some(Arc::new(get_workspace.run(())?))}else{None};
        if let Some(workspace)=workspace.as_deref(){crate::cloning::methods::domestication_review::validate_origin(workspace,&saved.source.binding,&saved.molecule,selected.as_ref())?;}
        Ok(Arc::new(MaterialSource{origin:selected,workspace}))
    });
    let guard = Callback::new(move |saved: Arc<Scope>| -> Result<(), String> {
        if source.get_untracked() != saved.source
            || epoch.get_untracked() != saved.epoch
            || form.get_untracked() != saved.form
            || !molecule.with_untracked(|current| current == saved.molecule.as_ref())
        {
            return Err(
                "The exact source event or declarations changed; inspect, propose and check again."
                    .into(),
            );
        }
        Ok(())
    });
    let capture = Callback::new(move |()| -> Result<Arc<Scope>, String> {
        let source = source.get_untracked();
        let settings = form.get_untracked();
        let (original, request) = molecule.with_untracked(|original| -> Result<_, String> {
            collagene_model::bounds::model_preflight(&[original])
                .map_err(|error| error.to_string())?;
            let request = settings.request(original)?;
            if request.source_revision != source.binding.revision {
                return Err("The source event and original DNA do not agree.".into());
            }
            Ok((original.clone(), request))
        })?;
        Ok(Arc::new(Scope {
            source,
            molecule: Arc::new(original),
            epoch: epoch.get_untracked(),
            form: settings,
            request,
        }))
    });
    let run_alive = alive.clone();
    let run = Callback::new(move |work: Work| {
        if pending.get_untracked() {
            return;
        }
        let result = (|| -> Result<_, String> {
            let mut material = None;
            let (scope, revision, job) = match &work {
                Work::Inspect => {
                    let scope = capture.run(())?;
                    let revision = design::inspect_revision(&scope.molecule, &scope.request)
                        .map_err(|e| e.to_string())?;
                    let job = AnalysisJob::InspectDomestication {
                        source: scope.molecule.as_ref().clone(),
                        request: scope.request.clone(),
                    };
                    (scope, revision, job)
                }
                Work::Propose => {
                    let scope = capture.run(())?;
                    let revision = design::propose_revision(&scope.molecule, &scope.request)
                        .map_err(|e| e.to_string())?;
                    let job = AnalysisJob::ProposeDomestication {
                        source: scope.molecule.as_ref().clone(),
                        request: scope.request.clone(),
                    };
                    (scope, revision, job)
                }
                Work::Check(saved) => {
                    guard.run(saved.scope.clone())?;
                    if !proposed.with_untracked(|current| {
                        current
                            .as_ref()
                            .is_some_and(|current| Arc::ptr_eq(current, saved))
                    }) {
                        return Err("The proposed candidate is stale.".into());
                    }
                    let selection = saved
                        .report
                        .candidate
                        .as_ref()
                        .ok_or("No concrete candidate was proposed.")?;
                    let revision = design::check_revision(
                        &saved.scope.molecule,
                        &saved.scope.request,
                        selection,
                    )
                    .map_err(|e| e.to_string())?;
                    let job = AnalysisJob::CheckDomestication {
                        source: saved.scope.molecule.as_ref().clone(),
                        request: saved.scope.request.clone(),
                        selection: selection.clone(),
                    };
                    (saved.scope.clone(), revision, job)
                }
                Work::Prepare(saved, settings) => {
                    guard.run(saved.scope.clone())?;
                    if !checked.with_untracked(|current| {
                        current
                            .as_ref()
                            .is_some_and(|current| Arc::ptr_eq(current, saved))
                    }) || route.get_untracked() != *settings
                    {
                        return Err("The checked substitutions or physical route changed.".into());
                    }
                    let captured =
                        attest_material.run((saved.scope.clone(), origin.get_untracked()))?;
                    let request = settings.request(
                        &saved.scope.molecule,
                        saved.scope.request.clone(),
                        saved.selection.clone(),
                        captured.origin.as_ref().map(|pin| pin.as_ref()),
                        captured.workspace.as_deref(),
                    )?;
                    let input = PlanInput {
                        id: physical::INPUT_ID.into(),
                        revision: saved.scope.source.binding.revision.clone(),
                        molecule: saved.scope.molecule.as_ref().clone(),
                    };
                    let revision =
                        physical::revision(&input, &request, captured.workspace.as_deref())
                            .map_err(|e| e.to_string())?;
                    let workspace = captured.workspace.as_deref().cloned();
                    material = Some(captured);
                    (
                        saved.scope.clone(),
                        revision,
                        AnalysisJob::PrepareDomestication {
                            input,
                            request,
                            workspace,
                        },
                    )
                }
            };
            Ok((
                scope,
                revision,
                job,
                Rc::new(CachedWorker::new()?),
                material,
            ))
        })();
        let (scope, revision, job, w, material) = match result {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        worker.set_value(Some(w.clone()));
        pending.set(true);
        match &work {
            Work::Inspect => inspection.set(None),
            Work::Propose => {
                proposed.set(None);
                checked.set(None);
                prepared.set(None);
            }
            Work::Check(_) => {
                checked.set(None);
                prepared.set(None);
            }
            Work::Prepare(_, _) => prepared.set(None),
        }
        status.set(
            match &work {
                Work::Inspect => "Inspecting actual source sites and coding constraints…",
                Work::Propose => "Searching source-bound synonymous substitutions in Rust…",
                Work::Check(_) => "Independently rescanning sites and affected coding sequences…",
                Work::Prepare(_, _) => {
                    "Reconstructing the selected source and declared PCR preparation…"
                }
            }
            .into(),
        );
        let alive = run_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let outcome = (|| -> Result<&'static str, String> {
                guard.run(scope.clone())?;
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Domestication response differs from its exact request.".into());
                }
                match (work, response.result) {
                    (_, AnalysisResult::Error(error)) => Err(error),
                    (Work::Inspect, AnalysisResult::DomesticationInspected(report)) => {
                        report
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        inspection.set(Some(Arc::new(Inspection { scope, report })));
                        Ok(
                            "Source inspection ready. Keeping a site is an explicit declaration; changing scope or coding selections requires a new proposal.",
                        )
                    }
                    (Work::Propose, AnalysisResult::DomesticationProposed(report)) => {
                        report
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        inspection.set(Some(Arc::new(Inspection {
                            scope: scope.clone(),
                            report: report.inspection.clone(),
                        })));
                        proposed.set(Some(Arc::new(Proposed { scope, report })));
                        patch_page.set(0);
                        Ok(
                            "Proposal ready for review. Explicitly Check the candidate before opening a reference or preparing any material.",
                        )
                    }
                    (
                        Work::Check(previous),
                        AnalysisResult::DomesticationChecked { report, target },
                    ) => {
                        report
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        let selection = previous
                            .report
                            .candidate
                            .as_ref()
                            .ok_or("The original candidate is unavailable.")?
                            .clone();
                        if report.selection.candidate_id != selection.candidate_id
                            || report.target_revision != target.revision()
                        {
                            return Err("Checked result belongs to another candidate or independent target.".into());
                        }
                        checked.set(Some(Arc::new(Checked {
                            scope,
                            selection,
                            report,
                            target,
                        })));
                        check_page.set(0);
                        Ok(
                            "Reference constraints checked. This target is reference-only; a physical preparation remains a separate explicit action.",
                        )
                    }
                    (
                        Work::Prepare(previous, settings),
                        AnalysisResult::DomesticationPrepared(report),
                    ) => {
                        if settings != route.get_untracked() {
                            return Err("The physical preparation route changed.".into());
                        }
                        let material =
                            material.ok_or("Missing captured physical source context.")?;
                        if origin.get_untracked() != material.origin {
                            return Err(
                                "The selected saved material changed; prepare again.".into()
                            );
                        }
                        let _ = attest_material
                            .run((previous.scope.clone(), material.origin.clone()))?;
                        report
                            .validate_binding(&revision)
                            .map_err(|e| e.to_string())?;
                        if report.target != previous.target
                            || report.check.context_revision != previous.report.context_revision
                        {
                            return Err("Physical preparation differs from the independently checked original substitutions.".into());
                        }
                        prepared.set(Some(Arc::new(Prepared {
                            material,
                            checked: previous,
                            route: settings,
                            report,
                        })));
                        Ok(
                            "Actual primer and reaction reconstruction ready for review. Open the plan for its normal verification and output workflow.",
                        )
                    }
                    _ => Err("Unexpected domestication result type.".into()),
                }
            })();
            match outcome {
                Ok(message) => status.set(message.into()),
                Err(error) => status.set(error),
            }
        });
    });
    let open_alive = alive.clone();
    let open_reference = Callback::new(move |saved: Arc<Checked>| {
        let result = (|| -> Result<(), String> {
            guard.run(saved.scope.clone())?;
            if !checked.with_untracked(|current| {
                current
                    .as_ref()
                    .is_some_and(|current| Arc::ptr_eq(current, &saved))
            }) {
                return Err("The checked reference is stale.".into());
            }
            let selected = origin.get_untracked();
            let _ = attest_material.run((saved.scope.clone(), selected.clone()))?;
            on_open.run(Open {
                origin: selected,
                source: saved.scope.source.binding.clone(),
                epoch: saved.scope.epoch,
                request: saved.scope.request.clone(),
                selection: saved.selection.clone(),
                check: saved.report.clone(),
                target: saved.target.clone(),
                prepared: None,
            })
        })();
        if let Err(error) = result {
            if open_alive.load(Ordering::Relaxed) {
                status.set(error);
            }
        }
    });
    let plan_alive = alive.clone();
    let open_plan = Callback::new(move |saved: Arc<Prepared>| {
        let result = (|| -> Result<(), String> {
            guard.run(saved.checked.scope.clone())?;
            if route.get_untracked() != saved.route
                || !prepared.with_untracked(|current| {
                    current
                        .as_ref()
                        .is_some_and(|current| Arc::ptr_eq(current, &saved))
                })
            {
                return Err("The physical preparation is stale.".into());
            }
            if origin.get_untracked() != saved.material.origin {
                return Err("The selected material origin changed; prepare again.".into());
            }
            let _ = attest_material
                .run((saved.checked.scope.clone(), saved.material.origin.clone()))?;
            on_open.run(Open {
                origin: saved.material.origin.clone(),
                source: saved.checked.scope.source.binding.clone(),
                epoch: saved.checked.scope.epoch,
                request: saved.checked.scope.request.clone(),
                selection: saved.checked.selection.clone(),
                check: saved.checked.report.clone(),
                target: saved.checked.target.clone(),
                prepared: Some(saved.report.clone()),
            })
        })();
        if let Err(error) = result {
            if plan_alive.load(Ordering::Relaxed) {
                status.set(error);
            }
        }
    });
    let navigate_alive = alive.clone();
    let navigate = Callback::new(move |range: (usize, usize)| {
        if let Err(error) = on_navigate.run(range) {
            if navigate_alive.load(Ordering::Relaxed) {
                status.set(error);
            }
        }
    });
    let check_current = Memo::new(move |_| {
        checked.with(|value| {
            value.as_ref().is_some_and(|saved| {
                source.get() == saved.scope.source
                    && epoch.get() == saved.scope.epoch
                    && form.get() == saved.scope.form
            })
        })
    });
    let capture_preparation = Callback::new(move |()| -> Result<_, String> {
        let saved = prepared
            .get_untracked()
            .ok_or("Reconstruct a physical preparation first")?;
        guard.run(saved.checked.scope.clone())?;
        if route.get_untracked() != saved.route {
            return Err("Preparation route or settings changed.".into());
        }
        let input = PlanInput {
            id: physical::INPUT_ID.into(),
            revision: saved.checked.scope.molecule.revision(),
            molecule: saved.checked.scope.molecule.as_ref().clone(),
        };
        if origin.get_untracked() != saved.material.origin {
            return Err("The selected material origin changed; prepare again.".into());
        }
        let current =
            attest_material.run((saved.checked.scope.clone(), saved.material.origin.clone()))?;
        let context =
            physical::revision(&input, &saved.report.request, current.workspace.as_deref())
                .map_err(|e| e.to_string())?;
        saved
            .report
            .validate_binding(&context)
            .map_err(|e| e.to_string())?;
        Ok((saved.report.plan.clone(), saved.report.target.clone()))
    });
    let preparation_identity = Memo::new(move |_| {
        prepared
            .get()
            .map(|saved| {
                (
                    collagene_core::planning::plan::DesignPlan::revision(&saved.report.plan),
                    saved.report.target.revision(),
                )
            })
            .unwrap_or_default()
    });
    let preparation_current = Memo::new(move |_| {
        prepared.with(|value| {
            value.as_ref().is_some_and(|saved| {
                check_current.get()
                    && route.get() == saved.route
                    && origin.get() == saved.material.origin
            })
        })
    });
    view!{<section class="domestication-design"><h2>"Remove unwanted restriction sites"</h2>
        <p>{move||source.with(|source|format!("Source: {}",source.name))}</p>
        {move||source.with(|source|view!{<details><summary>"Exact source event"</summary><p>{source.binding.document_id.clone()}</p><code>{source.binding.event_id.clone()}</code><p>"Revision"</p><code>{source.binding.revision.clone()}</code></details>})}
        <crate::history::run_origin_view::SourceOriginPicker source=Memo::new(move|_|Some(source.get())) epoch get_workspace selected=Memo::new(move|_|origin.get()) on_select=Callback::new(move|selected:Option<PinnedOrigin>|{
            if let Some(pin)=selected.as_ref(){
                let frame=get_workspace.run(())?;let current=source.get_untracked();
                molecule.with_untracked(|m|crate::cloning::methods::domestication_review::validate_origin(&frame,&current.binding,m,Some(pin)))?;
            }
            origin.set(selected);Ok(())
        })/>
        <p class="helper">"Saved material changes the preparation source, not the independently checked substitutions. Unsupported source geometry stays diagnostic."</p>
        <p class="helper">"Select the sites to remove and the coding sequences whose protein must be preserved. Every affected CDS is checked, including unselected overlaps. A checked reference is design intent; only a separately reconstructed physical route can propose material."</p>
        <fieldset disabled=move||pending.get()><legend>"Source declarations"</legend>
            <label>"Enzyme catalog IDs (comma-separated)"<input maxlength="2048" prop:value=move||form.with(|value|value.enzymes.clone()) on:input=move|event|form.update(|value|value.enzymes=event_target_value(&event))/></label>
            <label>"Scope start (1-based)"<input type="number" min="1" prop:value=move||form.with(|value|value.scope_start.clone()) on:input=move|event|form.update(|value|value.scope_start=event_target_value(&event))/></label>
            <label>"Scope end (inclusive)"<input type="number" min="1" prop:value=move||form.with(|value|value.scope_end.clone()) on:input=move|event|form.update(|value|value.scope_end=event_target_value(&event))/></label>
            <p class="helper">"On circular DNA an end before the start declares an origin-spanning interval."</p>
            <label>"Protected ranges (start-end, comma-separated)"<input maxlength="8192" prop:value=move||form.with(|value|value.protected.clone()) on:input=move|event|form.update(|value|value.protected=event_target_value(&event))/></label>
            <details><summary>"Search limits"</summary>
                <label>"Beam width"<input type="number" min="1" prop:value=move||form.with(|value|value.beam.clone()) on:input=move|event|form.update(|value|value.beam=event_target_value(&event))/></label>
                <label>"Work budget"<input type="number" min="1" prop:value=move||form.with(|value|value.work.clone()) on:input=move|event|form.update(|value|value.work=event_target_value(&event))/></label>
                <label>"Maximum substitutions"<input type="number" min="1" prop:value=move||form.with(|value|value.max_patches.clone()) on:input=move|event|form.update(|value|value.max_patches=event_target_value(&event))/></label>
                <p class="helper">"A pruned or incomplete search does not prove infeasibility or a global optimum."</p>
            </details>
            <h3>"Coding sequences"</h3>
            {move||molecule.with(|molecule|{
                let cds=molecule.features.iter().enumerate().filter(|(_,feature)|feature.kind.eq_ignore_ascii_case("CDS")).collect::<Vec<_>>();let count=cds.len();
                view!{<section>{cds.into_iter().skip(feature_page.get().min(pages(count,LIST_PAGE)-1)*LIST_PAGE).take(LIST_PAGE).map(|(index,feature)|{
                    let first=feature.location.spans().first().map(|part|(part.0,part.1));let name=feature_name(feature,index);
                    view!{<article><label><input type="checkbox" prop:checked=move||form.with(|value|value.features.contains(&index)) on:change=move|event|form.update(|value|{if event_target_checked(&event){if !value.features.contains(&index){value.features.push(index);value.features.sort_unstable();}}else{value.features.retain(|selected|*selected!=index);}})/>{format!("{} · feature {}{}",name,index+1,if feature.needs_review{" · needs review"}else{""})}</label>
                        {first.map(|range|view!{<button class="text-tool" on:click=move|_|navigate.run(range)>"View first source part"</button>})}
                    </article>}
                }).collect_view()}<Pager page=feature_page total=count limit=LIST_PAGE label="coding features"/></section>}
            })}
        </fieldset>
        <details><summary>"Deliberately retained recognition occurrences"</summary>
            <p class="helper">"These are explicit unchanged-site declarations. If you change the catalog selection, remove obsolete declarations explicitly before searching again."</p>
            {move||{let retained=form.with(|value|value.retained.clone());let count=retained.len();view!{<section>
                {retained.into_iter().skip(retained_page.get().min(pages(count,SITE_PAGE)-1)*SITE_PAGE).take(SITE_PAGE).map(|pin|{
                    let label=format!("{} · source position {} · {} strand",pin.enzyme_id,pin.start+1,if pin.reverse{"reverse"}else{"forward"});
                    view!{<p>{label}<button class="text-tool" disabled=move||pending.get() on:click=move|_|form.update(|value|value.retained.retain(|kept|kept!=&pin))>"Remove keep declaration"</button></p>}
                }).collect_view()}<Pager page=retained_page total=count limit=SITE_PAGE label="retained occurrences"/>
            </section>}}}
        </details>
        <button class="button secondary" disabled=move||pending.get() on:click=move|_|run.run(Work::Inspect)>"Inspect source sites and coding constraints"</button>
        <button class="button primary" disabled=move||pending.get() on:click=move|_|run.run(Work::Propose)>"Propose substitutions"</button>
        <button class="text-tool" disabled=move||!pending.get() on:click=move|_|{cancel.run(false);status.set("Domestication computation cancelled. Original source unchanged.".into());}>"Cancel computation"</button>
        {move||inspection.get().map(|saved|{
            let count=saved.report.sites.len();let diagnostic_count=saved.report.features.len();let note_count=saved.report.blocking.len()+saved.report.warnings.len();
            view!{<section><h3>"Last source inspection"</h3><p class="helper">"Site positions remain from this exact source and catalog. Scope and coding-selection labels below describe the last inspection; inspect again after changing declarations."</p>
                {saved.report.blocking.iter().chain(&saved.report.warnings).skip(inspection_note_page.get().min(pages(note_count,LIST_PAGE)-1)*LIST_PAGE).take(LIST_PAGE).map(|message|view!{<p>{message.clone()}</p>}).collect_view()}
                <Pager page=inspection_note_page total=note_count limit=LIST_PAGE label="inspection findings"/>
                {saved.report.sites.iter().skip(site_page.get().min(pages(count,SITE_PAGE)-1)*SITE_PAGE).take(SITE_PAGE).map(|site|{
                    let pin=site.pin.clone();let recognition=site.recognition.clone();let position=pin.start;let length=recognition.len();
                    let guard_source=saved.scope.clone();
                    view!{<article><p>{format!("{} · position {} · {} strand · {} · in last scope {}",pin.enzyme_id,position+1,if pin.reverse{"reverse"}else{"forward"},recognition,site.in_scope)}</p>
                        <label><input type="checkbox" disabled=move||pending.get() prop:checked={let pin=pin.clone();move||form.with(|value|value.retained.contains(&pin))} on:change=move|event|{
                            if source.get_untracked()!=guard_source.source||epoch.get_untracked()!=guard_source.epoch||form.with_untracked(|value|value.enzymes!=guard_source.form.enzymes){return;}
                            let keep=event_target_checked(&event);form.update(|value|{if keep{if !value.retained.contains(&pin){value.retained.push(pin.clone());}}else{value.retained.retain(|value|value!=&pin);}});
                        }/>"Keep this exact recognition occurrence unchanged"</label>
                        <button class="text-tool" on:click=move|_|molecule.with_untracked(|molecule|{let end=position.saturating_add(length).min(molecule.sequence.len());navigate.run((position,end));})>"View source recognition"</button>
                    </article>}
                }).collect_view()}<Pager page=site_page total=count limit=SITE_PAGE label="recognition sites"/>
                <details><summary>"Coding support and translation diagnostics"</summary>
                    {saved.report.features.iter().enumerate().skip(diagnostic_page.get().min(pages(diagnostic_count,LIST_PAGE)-1)*LIST_PAGE).take(LIST_PAGE).map(|(index,_)|view!{<InspectionCds saved=saved.clone() index/>}).collect_view()}
                    <Pager page=diagnostic_page total=diagnostic_count limit=LIST_PAGE label="coding diagnostics"/>
                </details>
            </section>}
        })}
        {move||proposed.get().map(|saved|{
            let candidate=saved.report.candidate.as_ref();let count=candidate.map_or(0,|candidate|candidate.patches.len());let check=saved.clone();
            view!{<section><h3>"Proposed reference substitutions"</h3><p>{format!("Search result {:?} · complete within declared model {}",saved.report.status,saved.report.search_complete)}</p>
                {saved.report.searches.iter().skip(search_page.get().min(pages(saved.report.searches.len(),LIST_PAGE)-1)*LIST_PAGE).take(LIST_PAGE).map(|search|view!{<p>{format!("Feature {} · {:?} · work {} · search complete {} · minimum changed-base count proven {}",search.feature_index+1,search.status,search.work_used,search.search_complete,search.minimum_changes_proven)}</p>}).collect_view()}
                <Pager page=search_page total=saved.report.searches.len() limit=LIST_PAGE label="coding searches"/>
                {saved.report.warnings.iter().map(|message|view!{<p>{message.clone()}</p>}).collect_view()}
                {candidate.map(|candidate|candidate.patches.iter().skip(patch_page.get().min(pages(count,LIST_PAGE)-1)*LIST_PAGE).take(LIST_PAGE).map(|patch|{
                    let position=patch.position;let scope=saved.scope.clone();
                    view!{<article><p>{format!("Source position {}: {} → {}",patch.position+1,patch.expected,patch.replacement)}</p>
                        <pre>{source_window(&saved.scope.molecule,position)}</pre>
                        <button class="text-tool" on:click=move|_|if guard.run(scope.clone()).is_ok(){navigate.run((position,position+1))}>"View this source substitution"</button>
                    </article>}
                }).collect_view())}
                <Pager page=patch_page total=count limit=LIST_PAGE label="substitutions"/>
                <button class="button primary" disabled=move||pending.get()||count==0 on:click=move|_|run.run(Work::Check(check.clone()))>"Check this explicit candidate"</button>
                <details><summary>"Exact proposal metadata"</summary><code>{saved.report.context_revision.clone()}</code>{candidate.map(|candidate|view!{<p>"Candidate identity"</p><code>{candidate.candidate_id.clone()}</code>})}</details>
            </section>}
        })}
        {move||checked.get().map(|saved|{
            let open=saved.clone();let suggest=StoredValue::new(saved.clone());let count=saved.report.features.len();
            view!{<section><h3>"Checked intended reference — no material created"</h3><p>{format!("{} · {} bases · {} substitutions",saved.target.name,saved.target.sequence.len(),saved.selection.patches.len())}</p>
                {saved.report.warnings.iter().map(|message|view!{<p>{message.clone()}</p>}).collect_view()}
                <details><summary>"Preserved coding sequences"</summary>{saved.report.features.iter().enumerate().skip(check_page.get().min(pages(count,LIST_PAGE)-1)*LIST_PAGE).take(LIST_PAGE).map(|(index,_)|view!{<CheckedCds saved=saved.clone() index/>}).collect_view()}<Pager page=check_page total=count limit=LIST_PAGE label="checked coding sequences"/></details>
                <crate::assembly::review::assembly_product_view::ReferenceCanvas product=saved.target.clone() identity=saved.report.target_revision.clone() current=check_current on_message=Callback::new(move|message|status.set(message))/>
                <button class="button secondary" disabled=move||pending.get() on:click=move|_|open_reference.run(open.clone())>"Open reference-only design"</button>
                <details><summary>"Choose a physical preparation"</summary>
                    <p class="helper">"Choose the preparation independently after checking source, feature and restriction constraints. The intended target never becomes an amplification template."</p>
                    <label>"Physical preparation route"<select prop:value=move||match route.get().kind{crate::cloning::methods::domestication_form::PreparationKind::InversePcr=>"inverse",crate::cloning::methods::domestication_form::PreparationKind::InversePcrRounds=>"rounds",crate::cloning::methods::domestication_form::PreparationKind::OverlapExtension=>"extension"} on:change=move|e|route.update(|r|r.kind=match event_target_value(&e).as_str(){"extension"=>crate::cloning::methods::domestication_form::PreparationKind::OverlapExtension,"rounds"=>crate::cloning::methods::domestication_form::PreparationKind::InversePcrRounds,_=>crate::cloning::methods::domestication_form::PreparationKind::InversePcr})><option value="inverse">"Inverse PCR — circular source"</option><option value="rounds">"Sequential inverse PCR — scattered changes"</option><option value="extension">"Overlap-extension PCR — whole linear or circular source"</option></select></label>
                    <Show when=move||route.get().kind==crate::cloning::methods::domestication_form::PreparationKind::InversePcr><p class="helper">"Every changed base must fit the declared replacement interval and actual primer tails. Complete oligos are limited to200 bases. Scattered edits can require a different route or reviewed rounds; no recoded template is invented."</p>
                    <button class="text-tool" disabled=move||pending.get() on:click=move|_|{let suggest=suggest.get_value();if guard.run(suggest.scope.clone()).is_ok(){if let (Some(first),Some(last))=(suggest.selection.patches.first(),suggest.selection.patches.last()){route.update(|value|{value.start=(first.position+1).to_string();value.end=(last.position+1).to_string();});}}}>"Suggest straight interval covering these substitutions"</button>
                    <p class="helper">"This suggestion preserves original numbering. Enter an origin-spanning interval manually when appropriate; suggestion is not preparation approval."</p>
                    <fieldset disabled=move||pending.get()>
                        <label>"Replacement start (1-based)"<input type="number" min="1" prop:value=move||route.with(|value|value.start.clone()) on:input=move|event|route.update(|value|value.start=event_target_value(&event))/></label>
                        <label>"Replacement end (inclusive)"<input type="number" min="1" prop:value=move||route.with(|value|value.end.clone()) on:input=move|event|route.update(|value|value.end=event_target_value(&event))/></label>
                        <label>"Forward annealing length"<input type="number" min="8" max="80" prop:value=move||route.with(|value|value.forward.clone()) on:input=move|event|route.update(|value|value.forward=event_target_value(&event))/></label>
                        <label>"Reverse annealing length"<input type="number" min="8" max="80" prop:value=move||route.with(|value|value.reverse.clone()) on:input=move|event|route.update(|value|value.reverse=event_target_value(&event))/></label>
                        <label>"Replacement bases assigned to reverse tail"<input type="number" min="0" prop:value=move||route.with(|value|value.reverse_tail.clone()) on:input=move|event|route.update(|value|value.reverse_tail=event_target_value(&event))/></label>
                    </fieldset>
                    </Show>
                    <Show when=move||route.get().kind==crate::cloning::methods::domestication_form::PreparationKind::InversePcrRounds>
                        <label>"Ordered rounds in original source coordinates"<textarea prop:value=move||route.get().rounds maxlength="4096" placeholder="21-21:18:18:0
75-75:18:18:0" on:input=move|e|route.update(|r|r.rounds=event_target_value(&e))/></label>
                        <p class="helper">"One start-end:forward_length:reverse_length:reverse_tail_bases declaration per line. Intervals are 1-based inclusive original source coordinates, including circular origin crossings. Each selected change must belong to exactly one round. Later PCRs consume the preceding closed product; every primer pair is reviewed before opening the plan."</p>
                    </Show>
                    <Show when=move||route.get().kind==crate::cloning::methods::domestication_form::PreparationKind::OverlapExtension><OverlapControls route pending/><p class="helper">"Amplify the original source in two PCRs, extend its explicitly edited overlap, then amplify with outer primers. Circular output additionally declares end phosphorylation and ligation. Linear output remains linear. No added outer tails, automatic source repair or independent-target replacement."</p></Show>
                    <button class="button primary" disabled=move||pending.get() on:click=move|_|run.run(Work::Prepare(saved.clone(),route.get_untracked()))>"Reconstruct chosen physical preparation"</button>
                </details>
            </section>}
        })}
        {move||prepared.get().map(|saved|{
            let open=saved.clone();let physical=&saved.report;let actual=&physical.plan;
            let complete=physical.verification.products.len()==actual.steps.len();
            let note_count=physical.verification.findings.len()+physical.verification.assumptions.len()+physical.warnings.len();
            view!{<section><h3>"Actual physical preparation review"</h3><p>{format!("Native plan assessment: {:?}. Completed {} of {} operations.",physical.verification.status,physical.verification.products.len(),actual.steps.len())}</p>
                <h4>"Primer proposal"</h4><p>{if physical.proposal_ready{"Ready under the declared construction and primer constraints. This is separate from sequence verification and experimental success."}else{"Not ready under the declared primer constraints. The actual sequence verification remains unchanged; open this as a diagnostic candidate for review or manual editing."}}</p><ul>{physical.proposal_findings.iter().map(|finding|view!{<li>{finding.clone()}</li>}).collect_view()}</ul>
                <PreparationEvidence method=Arc::new(physical.method.clone())/>
                <crate::assembly::review::assembly_outputs::AssemblyOutputs capture=capture_preparation identity=preparation_identity epoch/>

                {physical.verification.findings.iter().chain(&physical.verification.assumptions).chain(&physical.warnings).skip(preparation_note_page.get().min(pages(note_count,LIST_PAGE)-1)*LIST_PAGE).take(LIST_PAGE).map(|message|view!{<p>{message.clone()}</p>}).collect_view()}
                <Pager page=preparation_note_page total=note_count limit=LIST_PAGE label="preparation findings"/>
                {if complete {physical.verification.products.last().map(|product|view!{<h4>"Reconstructed final product"</h4><crate::assembly::review::assembly_product_view::ReferenceCanvas product=product.clone() identity=product.revision() current=preparation_current on_message=Callback::new(move|message|status.set(message))/>}).into_any()}else{view!{<p>"A full physical product was not reconstructed. A completed preparation prefix is not the final product."</p>}.into_any()}}
                <button class="button primary" disabled=move||pending.get() on:click=move|_|open_plan.run(open.clone())>{if physical.proposal_ready{"Open plan for verification and outputs"}else{"Open diagnostic plan for review or editing"}}</button>
            </section>}
        })}
        <p role="status">{move||status.get()}</p>
    </section>}.into_any()
}

#[component]
fn InspectionCds(saved: Arc<Inspection>, index: usize) -> impl IntoView {
    let saved = StoredValue::new(saved);
    let page = RwSignal::new(0usize);
    let (label,detail,total)=saved.with_value(|saved|{
        let feature=&saved.report.features[index];
        (feature.label.clone(),format!("Feature {} · selected at inspection {} · supported {} · code {:?} · codon start {:?} · terminal stop {:?}",feature.pin.index+1,feature.selected,feature.supported,feature.genetic_code,feature.codon_start,feature.has_terminal_stop),feature.protein.as_ref().map_or(0,String::len))
    });
    view!{<article><h4>{label}</h4><p>{detail}</p>
        {saved.with_value(|saved|saved.report.features[index].reason.as_ref().map(|reason|view!{<p>{reason.clone()}</p>}))}
        {move||saved.with_value(|saved|saved.report.features[index].protein.as_ref().map(|protein|{let start=page.get().min(pages(total,PEPTIDE_PAGE)-1)*PEPTIDE_PAGE;view!{<pre>{protein.chars().skip(start).take(PEPTIDE_PAGE).collect::<String>()}</pre>}}))}
        <Pager page total limit=PEPTIDE_PAGE label="protein residues"/>
    </article>}.into_any()
}
#[component]
fn CheckedCds(saved: Arc<Checked>, index: usize) -> impl IntoView {
    let saved = StoredValue::new(saved);
    let page = RwSignal::new(0usize);
    let (label, detail, total) = saved.with_value(|saved| {
        let feature = &saved.report.features[index];
        (
            feature.label.clone(),
            format!(
                "Feature {} · code {} · codon start {} · terminal stop {}",
                feature.pin.index + 1,
                feature.genetic_code,
                feature.codon_start,
                feature.has_terminal_stop
            ),
            feature.before_protein.len(),
        )
    });
    view!{<article><h4>{label}</h4><p>{detail}</p>
        {move||saved.with_value(|saved|{let feature=&saved.report.features[index];let start=page.get().min(pages(total,PEPTIDE_PAGE)-1)*PEPTIDE_PAGE;let end=(start+PEPTIDE_PAGE).min(total);
            view!{<p>"Original protein"</p><pre>{feature.before_protein.chars().skip(start).take(end-start).collect::<String>()}</pre><p>"Checked target protein"</p><pre>{feature.after_protein.chars().skip(start).take(end-start).collect::<String>()}</pre>}
        })}<Pager page total limit=PEPTIDE_PAGE label="protein residues"/>
    </article>}.into_any()
}

#[component]
fn OverlapControls(route: RwSignal<RouteForm>, pending: RwSignal<bool>) -> impl IntoView {
    let labels = [
        "Minimum annealing length",
        "Maximum annealing length",
        "Target annealing Tm (°C)",
        "Maximum pair ΔTm (°C)",
        "Minimum GC fraction",
        "Maximum GC fraction",
        "Maximum homopolymer run",
        "Maximum 3′ complementarity run",
        "Ranked pairs per primer call",
        "Homology flank per edit (bp)",
        "Maximum ordered oligo length",
        "Minimum paired overlap bases",
        "Extension comparison limit",
        "Extension candidate limit",
        "Sodium concentration (M)",
        "Total strand concentration (M)",
    ];
    view! {<details><summary>"Overlap-extension primer and search settings"</summary><fieldset disabled=move||pending.get()>{labels.into_iter().enumerate().map(|(index,label)|view!{<label>{label}<input maxlength="64" prop:value=move||route.with(|r|r.overlap_extension.values[index].clone()) on:input=move|e|route.update(|r|r.overlap_extension.values[index]=event_target_value(&e))/></label>}).collect_view()}</fieldset><p class="helper">"Tm and complementarity are bounded sequence models, not empirical amplification predictions. Search limits and incomplete evidence remain visible. No declared outer tails are discarded."</p></details>}
}
fn primer_text(primer: &collagene_core::sequence::analysis::Primer) -> String {
    format!(
        "{}\nOrdered 5′→3′: {}{}\nAnnealing domain: {}\n5′ tail: {}",
        primer.name,
        primer.tail,
        primer.annealing.as_str(),
        primer.annealing.as_str(),
        primer.tail
    )
}
#[component]
fn PreparationEvidence(method: Arc<physical::MethodEvidence>) -> impl IntoView {
    match method.as_ref(){
  physical::MethodEvidence::InversePcr{forward,reverse,..}=>view!{<details open><summary>"Inverse-PCR primer pair"</summary><pre>{format!("{}\n\n{}",primer_text(forward),primer_text(reverse))}</pre></details>}.into_any(),
  physical::MethodEvidence::InversePcrRounds{rounds}=>view!{<section><h4>"Sequential inverse-PCR primer pairs"</h4>{rounds.iter().enumerate().map(|(index,round)|view!{<details><summary>{format!("Round {} · original {}–{}",index+1,round.round.edit_start+1,round.round.edit_end)}</summary><p>{format!("Template {:?}; intended source origin {} (zero-based original coordinate).",round.template,round.source_origin)}</p><pre>{format!("{}\n\n{}",primer_text(&round.forward),primer_text(&round.reverse))}</pre></details>}).collect_view()}<p class="helper">"The canonical graph retains every round, actual PCR alternatives, phosphorylation and circular ligation. Intermediate intended frames are primer declarations, never extra template stocks."</p></section>}.into_any(),
  physical::MethodEvidence::OverlapExtension{source_frame,preparation,designs,search,..}=>view!{<section><p>{format!("Template {:?} · {:?} · {} actual bases",source_frame.template_strand,source_frame.topology,source_frame.sequence.len())}</p><h4>"Three actual PCR primer pairs"</h4>{[("Left PCR",&preparation.left),("Right PCR",&preparation.right),("Outer PCR",&preparation.outer)].into_iter().map(|(label,pair)|view!{<details><summary>{label}</summary><pre>{format!("{}\n\n{}",primer_text(&pair.forward),primer_text(&pair.reverse))}</pre></details>}).collect_view()}{designs.iter().enumerate().map(|(i,report)|view!{<details><summary>{format!("Primer call {} diagnostics",i+1)}</summary><p>{format!("{} ranked pairs retained · {} evaluated · {} rejected",report.pairs.len(),report.evaluated_pairs,report.rejected_pairs)}</p>{report.warnings.iter().map(|w|view!{<p>{w.clone()}</p>}).collect_view()}{search.get(i).map(|s|view!{<p>{format!("{} · selected rank {} · {} examined · {} rejected · capped {} · unique geometry {}",s.call,s.selected_index+1,s.examined,s.rejected,s.capped,s.unique_geometry)}</p>{s.notes.iter().map(|n|view!{<p>{n.clone()}</p>}).collect_view()}{report.pairs.get(s.selected_index).map(|pair|view!{<p>{format!("Annealing Tm forward {:?} / reverse {:?} °C · ΔTm {:.2} °C · heterodimer max3′run {}",pair.forward.analysis.tm.as_ref().map(|t|t.tm_c),pair.reverse.analysis.tm.as_ref().map(|t|t.tm_c),pair.tm_difference_c,pair.heterodimer.max_three_prime_run)}</p>})}})}</details>}).collect_view()}<p class="helper">"Order sheets and HTML/PDF use the canonical reconstructed plan above, including all six ordered oligos and declared circular closure steps."</p></section>}.into_any(),
 }
}
