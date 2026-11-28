//! Single explicitly formatted source-part inspection against a pinned grammar.
use crate::{platform::worker::CachedWorker, sources::identity::SourceInfo};
use collagene_assembly_standards::assessment::{FindingKind, Limits, Side};
use collagene_core::{
    cloning::assembly::standards as native,
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
#[derive(Clone, PartialEq, Eq)]
struct Draft {
    profile: String,
    slot: String,
    start: String,
    length: String,
    reverse: bool,
    body_start: String,
    body_end: String,
    comparisons: String,
}
impl Default for Draft {
    fn default() -> Self {
        Self {
            profile: String::new(),
            slot: String::new(),
            start: "1".into(),
            length: String::new(),
            reverse: false,
            body_start: String::new(),
            body_end: String::new(),
            comparisons: "8000000".into(),
        }
    }
}
#[derive(Clone)]
struct Capture {
    source: SourceInfo,
    epoch: u64,
    generation: u32,
    draft: Draft,
    revision: String,
}
fn number(value: &str, label: &str) -> Result<usize, String> {
    value
        .trim()
        .parse()
        .map_err(|_| format!("Enter a whole {label}."))
}
#[component]
pub fn AssemblyStandardsPanel(
    source: Memo<SourceInfo>,
    get_molecule: Callback<(), Molecule>,
    epoch: RwSignal<u64>,
    on_select: Callback<Vec<(usize, usize)>>,
) -> impl IntoView {
    let catalogue = StoredValue::new(native::catalogue().map_err(|e| e.to_string()));
    let draft = RwSignal::new(Draft::default());
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let report = RwSignal::new(None::<(Capture, Arc<native::Report>)>);
    let page = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        let s = source.get();
        epoch.get();
        draft.set(Draft {
            length: s.length.to_string(),
            ..Default::default()
        });
    });
    Effect::new(move |_| {
        let _ = (source.get(), epoch.get(), draft.get());
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
        page.set(0);
        status.set("Choose a profile, slot, source interval and explicit body boundaries.".into());
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(w) = worker.get_value() {
                w.close();
            }
        });
    }
    let run = {
        let alive = alive.clone();
        Callback::new(move |(): ()| {
            if busy.get_untracked() {
                return;
            }
            let prepared = (|| {
                let d = draft.get_untracked();
                let c = catalogue.get_value()?;
                let p = c
                    .iter()
                    .find(|p| p.profile_id == d.profile)
                    .ok_or("Choose an installed profile.")?;
                let pin = source.get_untracked();
                let molecule = get_molecule.run(());
                if molecule.revision() != pin.binding.revision {
                    return Err("The current source changed. Select the document again.".into());
                }
                let start = number(&d.start, "source start")?
                    .checked_sub(1)
                    .ok_or("Source start is one-based (at least 1).")?;
                let body_start = number(&d.body_start, "body start")?
                    .checked_sub(1)
                    .ok_or("Body start is one-based (at least 1).")?;
                let request = native::Request {
                    profile_id: p.profile_id.clone(),
                    profile_identity: p.profile_identity.clone(),
                    parts: vec![native::SourcePart {
                        id: "selected-part".into(),
                        molecule,
                        start,
                        length: number(&d.length, "part length")?,
                        reverse: d.reverse,
                        slot: d.slot.clone(),
                        body_start,
                        body_end: number(&d.body_end, "body end")?,
                    }],
                    limits: Limits {
                        max_comparisons: number(&d.comparisons, "comparison budget")?,
                        ..Default::default()
                    },
                };
                request.validate().map_err(|e| e.to_string())?;
                let revision = native::revision(&request).map_err(|e| e.to_string())?;
                Ok::<_, String>((pin, d, request, revision))
            })();
            let (pin, d, request, revision) = match prepared {
                Ok(v) => v,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            let w = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            if let Some(old) = worker.get_value() {
                old.close();
            }
            worker.set_value(Some(w.clone()));
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            let capture = Capture {
                source: pin,
                epoch: epoch.get_untracked(),
                generation: token,
                draft: d,
                revision,
            };
            busy.set(true);
            report.set(None);
            status.set("Inspecting the declared formatted part…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: capture.revision.clone(),
                        job: AnalysisJob::InspectAssemblyStandard { request },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                busy.set(false);
                worker.set_value(None);
                if source.get_untracked() != capture.source
                    || epoch.get_untracked() != capture.epoch
                    || draft.get_untracked() != capture.draft
                    || get_molecule.run(()).revision() != capture.source.binding.revision
                {
                    status.set("Source or part declarations changed; discard this result.".into());
                    return;
                }
                let result = result.and_then(|r| match r.result {
                    AnalysisResult::AssemblyStandard(r) => {
                        r.validate_binding(&capture.revision)
                            .map_err(|e| e.to_string())?;
                        if r.sources
                            .first()
                            .is_none_or(|s| s.revision != capture.source.binding.revision)
                        {
                            return Err("Standard report belongs to another source.".into());
                        }
                        Ok(r)
                    }
                    AnalysisResult::Error(e) => Err(e),
                    _ => Err("Unexpected assembly-standard result.".into()),
                });
                match result {
                    Ok(r) => {
                        report.set(Some((capture, Arc::new(r))));
                        status.set("Documentary format findings ready. No assembly or experimental success is inferred.".into());
                    }
                    Err(e) => status.set(e),
                }
            });
        })
    };
    let use_flanks = Callback::new(move |(): ()| {
        let result = (|| {
            let d = draft.get_untracked();
            let profiles = catalogue.get_value()?;
            let p = profiles
                .iter()
                .find(|p| p.profile_id == d.profile)
                .ok_or("Choose a profile.")?;
            let s = p
                .definition
                .grammar
                .slots
                .iter()
                .find(|s| s.id == d.slot)
                .ok_or("Choose a slot.")?;
            let length = number(&d.length, "part length")?;
            let end = length
                .checked_sub(s.right_addition.len())
                .ok_or("The interval is shorter than its declared suffix.")?;
            if s.left_addition.len() >= end {
                return Err("Declared flanks leave no body bases.".into());
            }
            draft.update(|d| {
                d.body_start = (s.left_addition.len() + 1).to_string();
                d.body_end = end.to_string();
            });
            Ok::<_, String>(())
        })();
        if let Err(e) = result {
            status.set(e);
        }
    });
    view! {
        <section class="inspector-content assembly-standards" aria-label="Assembly standard format" style="overflow-wrap:anywhere">
            <h3>
                "Assembly standard format"
            </h3>
            <p class="helper">
                "Inspect one explicitly declared formatted part. This does not recognize a vector class, simulate an assembly or establish biological function. Whole-construct grammar requirements may remain unsatisfied for a single part."
            </p>
            <p class="helper" role="alert">
                {
                    catalogue.with_value(|r|r.as_ref().err().cloned().unwrap_or_default())
                }
            </p>
            <label>
                "Installed profile"
                <select
                    aria-label="Assembly standard profile"
                    prop:value=move ||draft.with(|d|d.profile.clone())
                    on:change=move|e| {
                        let id = event_target_value(&e);
                        draft.update(|d| {
                            d.profile = id;
                            d.slot.clear();
                            d.body_start.clear();
                            d.body_end.clear();
                        });
                    }
                    >
                    <option value="">
                        "Choose profile"
                    </option>
                    {
                        catalogue.get_value().unwrap_or_default().into_iter().map(|p| view! {
                            <option value=p.profile_id>
                                {p.name}
                            </option>
                        }).collect_view()
                    }
                </select>
            </label>
            {
                move || catalogue.with_value(|c|c.as_ref().ok().and_then(|c|c.iter().find(|p|p.profile_id == draft
                    .with(|d|d.profile.clone()))).map(|p| view! {
                    <p class="helper">
                        {p.scope.clone()}
                    </p>
                    {
                        p.definition.limitations.iter().cloned().map(|s| view! {
                            <p class="helper">
                                {s}
                            </p>
                        }).collect_view()
                    }
                    <details>
                        <summary>
                            "Profile sources"
                        </summary>
                        {
                            p.definition.sources.iter().map(|s| view! {
                                <p class="helper">
                                    {format!("{} · {} · {} · {}", s.author, s.url, s.revision, s.locator)}
                                </p>
                            }).collect_view()
                        }
                    </details>
                }))
            }
            <label>
                "Declared slot"
                <select
                    aria-label="Assembly standard slot"
                    prop:value=move ||draft.with(|d|d.slot.clone())
                    on:change=move|e| {
                        let slot = event_target_value(&e);
                        draft.update(|d| {
                            d.slot = slot;
                            d.body_start.clear();
                            d.body_end.clear();
                        });
                    }
                    >
                    <option value="">
                        "Choose slot"
                    </option>
                    {
                        move || catalogue.with_value(|c|c.as_ref().ok().and_then(|c|c.iter().find(|p|p
                            .profile_id == draft.with(|d|d.profile.clone()))).map(|p|p.definition.grammar
                            .slots.iter().map(|s| view! {
                            <option value=s.id.clone()>
                                {s.name.clone()}
                            </option>
                        }).collect_view()))
                    }
                </select>
            </label>
            {
                move || catalogue.with_value(|c|c.as_ref().ok().and_then(|c|c.iter().find(|p|p.profile_id == draft
                    .with(|d|d.profile.clone()))).and_then(|p|p.definition.grammar.slots.iter().find(|s|s
                    .id == draft.with(|d|d.slot.clone()))).map(|s| view! {
                    <p class="helper">
                        {
                            format!("Expected prefix: {} · expected suffix: {} · overhangs: {} / {}", s.left_addition, s.right_addition, s.left_overhang,
                            s.right_overhang)
                        }
                    </p>
                    {
                        s.body_rules.iter().map(|r| view! {
                            <p class="helper">
                                {format!("Body rule: {r:?}")}
                            </p>
                        }).collect_view()
                    }
                }))
            }
            <label>
                "Source start (one-based)"
                <input
                    inputmode="numeric"
                    prop:value=move ||draft.with(|d|d.start.clone())
                    on:input=move|e| {
                        let v = event_target_value(&e);
                        draft.update(|d|d.start = v);
                    }
                    />
            </label>
            <label>
                "Formatted interval length (bases)"
                <input
                    inputmode="numeric"
                    prop:value=move ||draft.with(|d|d.length.clone())
                    on:input=move|e| {
                        let v = event_target_value(&e);
                        draft.update(|d|d.length = v);
                    }
                    />
            </label>
            <p class="helper">
                "The source interval is forward; a circular source may wrap once. Reverse orientation is applied after extraction. Body coordinates refer to this oriented formatted interval."
            </p>
            <label class="check-row">
                <input
                    type="checkbox"
                    prop:checked=move ||draft.with(|d|d.reverse)
                    on:change=move|e| {
                        let v = event_target_checked(&e);
                        draft.update(|d|d.reverse = v);
                    }
                    />
                "Reverse-complement interval"
            </label>
            <label>
                "Body start (one-based)"
                <input
                    inputmode="numeric"
                    prop:value=move ||draft.with(|d|d.body_start.clone())
                    on:input=move|e| {
                        let v = event_target_value(&e);
                        draft.update(|d|d.body_start = v);
                    }
                    />
            </label>
            <label>
                "Body end (inclusive)"
                <input
                    inputmode="numeric"
                    prop:value=move ||draft.with(|d|d.body_end.clone())
                    on:input=move|e| {
                        let v = event_target_value(&e);
                        draft.update(|d|d.body_end = v);
                    }
                    />
            </label>
            <button class="text-tool" on:click=move |_|use_flanks.run(())>
                "Set body from declared flank lengths"
            </button>
            <p class="helper">
                "This button declares a boundary; sequence conformity is checked only when you run the inspection."
            </p>
            <details>
                <summary>
                    "Analysis limit"
                </summary>
                <label>
                    "Maximum comparisons"
                    <input
                        inputmode="numeric"
                        prop:value=move ||draft.with(|d|d.comparisons.clone())
                        on:input=move|e| {
                            let v = event_target_value(&e);
                            draft.update(|d|d.comparisons = v);
                        }
                        />
                </label>
            </details>
            <button class="button secondary full" disabled=move ||busy.get() on:click=move |_|run.run(())>
                {
                    move || if busy.get() {"Inspecting…"}
                    else {"Inspect declared part"}
                }
            </button>
            <p class="helper" role="status">
                {move || status.get()}
            </p>
            {
                move || report.with(|r|r.as_ref().map(|(capture, r)| {
                    let c = capture.clone();
                    let report = r.clone();
                    let click_alive = alive.clone();
                    let select = Callback::new(move |index: usize| {
                        if !click_alive.load(Ordering::Relaxed) {
                            return;
                        }
                        if generation.get_untracked() != c.generation
                            || source.get_untracked() != c.source
                            || epoch .get_untracked() != c.epoch
                            || draft.get_untracked() != c.draft
                            || get_molecule .run(()).revision() != c.source.binding.revision {
                            status.set("This finding is stale; inspect the current declarations again.".into());
                            return;
                        }
                        let selected = (|| {
                            let f = report.assessment.findings.get(index).ok_or("Finding no longer exists.")?;
                            let i = f.part_index.ok_or("This grammar finding has no single source interval.")?;
                            let part = report.assessment.parts.get(i).ok_or("No retained part binding.")?;
                            let (a, b) = match &f.kind {
                                FindingKind::AvoidedSite { start, end, .. } => (*start, *end),
                                FindingKind::Addition { side: Side::Left, .. } => (0, part.body_start),
                                FindingKind::Addition { side: Side::Right, .. } => (
                                    part.body_end,
                                    part.sequence_length
                                ),
                                FindingKind::BodyRule { .. } => (
                                    part.body_start,
                                    part.body_end
                                ),
                                FindingKind::Grammar(_) => return Err(
                                    "This grammar finding is a declaration-level constraint.".to_string()
                                )
                            };
                            if a == b {
                                return Err("This finding concerns a missing addition; no existing bases occupy its interval.".into());
                            }
                            report.sources.get(i).ok_or("Source binding is unavailable.")?.map_interval(a, b).map_err(|e|e.to_string())
                        })();
                        match selected {Ok(ranges) => on_select.run(ranges), Err(e) => status.set(e)}
                    });
                    view! {
                        <p class="helper">
                            {
                                format!("{:?} · complete={} · {} comparisons", r.assessment.outcome, r.assessment.complete,
                                r.assessment.comparisons)
                            }
                        </p>
                        <p class="helper">
                            {
                                format!("Observed violations: {} · indeterminate: {} · advisories: {}", r.assessment.violations_observed, r.assessment.indeterminate_observed,
                                r.assessment.advisories_observed)
                            }
                        </p>
                        {
                            r.assessment.limitations.iter().cloned().map(|s| view! {
                                <p class="helper">
                                    {s}
                                </p>
                            }).collect_view()
                        }
                        {
                            r.assessment.findings.iter().enumerate().skip(page.get() * 40).take(40).map(|(i,
                            f)| view! {
                                <details>
                                    <summary>
                                        {format!("{:?}", f.severity)}
                                    </summary>
                                    <p>
                                        {format!("{:?}", f.kind)}
                                    </p>
                                    <button
                                        class="text-tool"
                                        disabled=(f.part_index.is_none() || matches!(&f.kind,FindingKind::Grammar(_)))
                                        on:click=move |_|select.run(i)>
                                        "Select source coordinates"
                                    </button>
                                </details>
                            }).collect_view()
                        }
                        <p class="helper">
                            {
                                format!("Allowed next slots: {} · declared list may finish: {}", r.successors.allowed_slots.join(", "),
                                r.successors.may_finish)
                            }
                        </p>
                    }
                }))
            }
            <Show when=move ||report.with(Option::is_some)>
                <p class="helper">
                    {
                        move || report.with(|r|r.as_ref().map(|(_, r)|format!("{} retained findings · page {} · budgets reached: {:?}", r.assessment
                            .findings.len(), page.get() + 1, r.assessment.budgets_hit)).unwrap_or_default())
                    }
                </p>
                <button
                    class="text-tool"
                    disabled=move ||page.get()==0
                    on:click=move |_|page.update(|p|*p=p.saturating_sub(1))>
                    "Previous findings"
                </button>
                <button
                    class="text-tool"
                    disabled=move ||report.with(|r|r.as_ref().is_none_or(|(_,r)|(page.get()+1)*40>=r.assessment.findings.len()))
                    on:click=move |_|page.update(|p|*p+=1)>
                    "Next findings"
                </button>
            </Show>
        </section>
    }
}
