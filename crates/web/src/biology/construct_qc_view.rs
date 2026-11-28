//! Current-event construct checks, with explicit declaration controls and worker ownership.
use super::construct_qc_model::Draft;
use crate::{platform::worker::CachedWorker, sources::identity::SourceInfo};
use collagene_core::{
    construct_qc as native,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::translation,
};
use collagene_model::{Molecule, Sequence};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
#[derive(Clone)]
struct Capture {
    source: SourceInfo,
    epoch: u64,
    generation: u32,
    draft: Draft,
    revision: String,
}
pub(super) fn feature_choice(
    rows: Memo<Vec<(usize, String, bool)>>,
    value: RwSignal<String>,
    label: &'static str,
) -> impl IntoView {
    view! {
        <label>
            {label}
            <select
                aria-label=label
                prop:value=move ||value.get()
                on:change=move |e|value.set(event_target_value(&e))>
                <option value="">
                    "Choose annotation"
                </option>
                {
                    move || rows.get().into_iter().map(|(i, s, _)| view! {
                        <option value=i.to_string()>
                            {s}
                        </option>
                    }).collect_view()
                }
            </select>
        </label>
    }
}
fn index(value: RwSignal<String>) -> Result<usize, String> {
    value
        .get_untracked()
        .parse()
        .map_err(|_| "Choose an annotation.".into())
}
#[component]
pub fn ConstructQcPanel(
    source: Memo<SourceInfo>,
    get_molecule: Callback<(), Molecule>,
    epoch: RwSignal<u64>,
    on_select: Callback<Vec<(usize, usize)>>,
) -> impl IntoView {
    let draft = RwSignal::new(Draft::default());
    let rows = Memo::new(move |_| {
        source.get();
        get_molecule
            .run(())
            .features
            .iter()
            .take(256)
            .enumerate()
            .map(|(i, f)| {
                (
                    i,
                    format!(
                        "{} · {} · {}",
                        i + 1,
                        f.kind,
                        f.label().chars().take(90).collect::<String>()
                    ),
                    f.needs_review,
                )
            })
            .collect::<Vec<_>>()
    });
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let report = RwSignal::new(None::<(Capture, Arc<native::Report>)>);
    let page = RwSignal::new(0usize);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        source.get();
        epoch.get();
        draft.set(Draft::default());
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
        status.set("Run checks for the current source and committed declarations.".into());
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
                let pin = source.get_untracked();
                let m = get_molecule.run(());
                if m.revision() != pin.binding.revision {
                    return Err("Current source changed; choose the document again.".to_owned());
                }
                let d = draft.get_untracked();
                let options = d.options_for(&m)?;
                let revision = native::revision(&m, &options).map_err(|e| e.to_string())?;
                Ok((pin, m, d, options, revision))
            })();
            let (pin, m, d, options, revision) = match prepared {
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
            status.set("Checking declared structural facts…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let result = w
                    .ask(Request {
                        id: token,
                        revision: capture.revision.clone(),
                        job: AnalysisJob::AnalyzeConstructQc {
                            molecule: m,
                            options,
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                if source.get_untracked() != capture.source
                    || epoch.get_untracked() != capture.epoch
                    || draft.get_untracked() != capture.draft
                    || get_molecule.run(()).revision() != capture.source.binding.revision
                {
                    status.set("Source or declarations changed; discard this result.".into());
                    return;
                }
                let outcome = result.and_then(|r| match r.result {
                    AnalysisResult::ConstructQc(r) => {
                        r.validate_binding(&capture.revision)
                            .map_err(|e| e.to_string())?;
                        if r.source_revision != capture.source.binding.revision {
                            return Err("QC response belongs to another source.".into());
                        }
                        Ok(r)
                    }
                    AnalysisResult::Error(e) => Err(e),
                    _ => Err("Unexpected construct QC response.".into()),
                });
                match outcome {
                    Ok(r) => {
                        report.set(Some((capture, Arc::new(r))));
                        status.set("Structural findings ready; no expression or experimental outcome is inferred.".into());
                    }
                    Err(e) => status.set(e),
                }
            });
        })
    };
    let ref_feature = RwSignal::new(String::new());
    let ref_id = RwSignal::new(String::new());
    let ref_sequence = RwSignal::new(String::new());
    let ref_url = RwSignal::new(String::new());
    let ref_accession = RwSignal::new(String::new());
    let ref_date = RwSignal::new(String::new());
    let ref_license = RwSignal::new(String::new());
    let add_reference = Callback::new(move |(): ()| {
        let result = (|| {
            let feature = if ref_feature.get_untracked().is_empty() {
                None
            } else {
                Some(index(ref_feature)?)
            };
            let m = get_molecule.run(());
            let literal = ref_sequence.get_untracked();
            if literal.len() > 1_000_000 {
                return Err("Reference input is limited to one million characters.".into());
            }
            let sequence =
                Sequence::new(m.sequence.alphabet(), &literal).map_err(|e| e.to_string())?;
            let id = ref_id.get_untracked();
            let reference = native::ReferenceElement {
                id: id.clone(),
                sequence_revision: native::reference_revision(&sequence)
                    .map_err(|e| e.to_string())?,
                sequence,
                provenance: native::Provenance {
                    source_url: ref_url.get_untracked(),
                    accession: ref_accession.get_untracked(),
                    retrieved_on: ref_date.get_untracked(),
                    license: ref_license.get_untracked(),
                },
            };
            let mut d = draft.get_untracked();
            d.options.references.push(reference);
            if let Some(i) = feature {
                d.options.reference_claims.push(native::ReferenceClaim {
                    feature_index: i,
                    reference_id: id,
                });
            }
            d.options_for(&m)?;
            draft.set(d);
            Ok::<_, String>(())
        })();
        if let Err(e) = result {
            status.set(e);
        }
    });
    let fusion_cds = RwSignal::new(String::new());
    let fusion_part = RwSignal::new(String::new());
    let fusion_phase = RwSignal::new("0".to_string());
    let add_fusion = Callback::new(move |(): ()| {
        let result = (|| {
            let mut d = draft.get_untracked();
            d.options.fusions.push(native::FusionClaim {
                cds_feature: index(fusion_cds)?,
                part_feature: index(fusion_part)?,
                expected_offset_mod3: fusion_phase
                    .get_untracked()
                    .parse()
                    .map_err(|_| "Choose phase 0, 1 or 2.")?,
                rationale:
                    "Explicit component boundaries are declared to share the selected coding phase."
                        .into(),
            });
            d.options_for(&get_molecule.run(()))?;
            draft.set(d);
            Ok::<_, String>(())
        })();
        if let Err(e) = result {
            status.set(e);
        }
    });
    let junction_cds = RwSignal::new(String::new());
    let junction_boundary = RwSignal::new(String::new());
    let junction_phase = RwSignal::new("0".to_string());
    let add_junction = Callback::new(move |(): ()| {
        let result = (|| {
            let mut d = draft.get_untracked();
            d.options.junctions.push(native::JunctionClaim {
                cds_feature: index(junction_cds)?,
                boundary: junction_boundary.get_untracked().parse().map_err(|_|"Enter the boundary after a whole base count (0..source length).")?,
                expected_offset_mod3: junction_phase.get_untracked().parse().map_err(|_|"Choose phase 0, 1 or 2.")?,
                rationale: "User-declared boundary in the current construct; no causal edit history inferred.".into(),
                history_document_id: None,
                history_event_id: None
            });
            d.options_for(&get_molecule.run(()))?;
            draft.set(d);
            Ok::<_, String>(())
        })();
        if let Err(e) = result {
            status.set(e);
        }
    });
    let phase = |value: RwSignal<String>, label: &'static str| {
        view! {
            <label>
                {label}
                <select
                    aria-label=label
                    prop:value=move ||value.get()
                    on:change=move |e|value.set(event_target_value(&e))>
                    <option value="0">
                        "0 · codon boundary"
                    </option>
                    <option value="1">
                        "1"
                    </option>
                    <option value="2">
                        "2"
                    </option>
                </select>
            </label>
        }
    };
    Effect::new(move |_| {
        source.get();
        epoch.get();
        ref_feature.set(String::new());
        ref_id.set(String::new());
        ref_sequence.set(String::new());
        ref_url.set(String::new());
        ref_accession.set(String::new());
        ref_date.set(String::new());
        ref_license.set(String::new());
        fusion_cds.set(String::new());
        fusion_part.set(String::new());
        junction_cds.set(String::new());
        junction_boundary.set(String::new());
    });
    view! {
        <section class="inspector-content construct-qc" aria-label="Construct QC" style="overflow-wrap:anywhere">
            <h3>
                "Construct QC"
            </h3>
            <p class="helper">
                {native::STATEMENT}
            </p>
            <p class="helper">
                "Current checks: CDS integrity, stored translation, literal references, declared fusion/junction register, structural element chains and literal copy requirements. Explicit initiation-pattern and bounded ORF checks are optional; vector-class and history-derived junction checks remain unevaluated."
            </p>
            <label>
                "Genetic code when CDS has no /transl_table"
                <select
                    aria-label="Construct QC genetic code"
                    prop:value=move ||draft.with(|d|d.options.genetic_code.to_string())
                    on:change=move|e|if let Ok(v) = event_target_value(&e).parse() {
                        draft.update(|d|d.options.genetic_code = v);
                    }
                    >
                    {
                        translation::genetic_codes().into_iter().map(|c| view! {
                            <option value=c.id.to_string()>
                                {format!("{} · {}", c.id, c.name)}
                            </option>
                        }).collect_view()
                    }
                </select>
            </label>
            <details open>
                <summary>
                    "Declared annotation roles"
                </summary>
                <p class="helper">
                    "Check only annotations whose roles and coordinates you explicitly assume for these checks. This does not verify biological identity. Features marked for review remain indeterminate."
                </p>
                {
                    move || rows.get().into_iter().map(|(i, label, review)| view! {
                        <label class="check-row">
                            <input
                                type="checkbox"
                                disabled=review
                                prop:checked=move ||draft.with(|d|d.assumed.contains(&i))
                                on:change=move|e| {
                                    let checked = event_target_checked(&e);
                                    draft.update(|d| {
                                        d.assumed.retain(|v| * v != i);
                                        if checked {
                                            d.assumed.push(i);
                                        }
                                    });
                                }
                                />
                            {label}
                            {
                                if review {" · review required"}
                                else {""}
                            }
                        </label>
                    }).collect_view()
                }
                <Show when=move ||rows.with(Vec::is_empty)>
                    <p class="helper">
                        "No annotations are present. Missing CDS annotation yields indeterminate findings."
                    </p>
                </Show>
            </details>
            <details>
                <summary>
                    "Pinned reference elements"
                </summary>
                {feature_choice(rows, ref_feature, "Reference feature")}
                <p class="helper">"Leave the annotation unset to capture a reference for literal copy-count checks. Selecting an annotation also adds a feature identity comparison."</p>
                <label>
                    "Reference ID"
                    <input
                        maxlength="128"
                        prop:value=move ||ref_id.get()
                        on:input=move |e|ref_id.set(event_target_value(&e))/>
                </label>
                <label>
                    "Literal oriented reference sequence"
                    <textarea
                        maxlength="1000000"
                        prop:value=move ||ref_sequence.get()
                        on:input=move |e|ref_sequence.set(event_target_value(&e))/>
                </label>
                <label>
                    "Source URL (HTTPS)"
                    <input
                        maxlength="2048"
                        prop:value=move ||ref_url.get()
                        on:input=move |e|ref_url.set(event_target_value(&e))/>
                </label>
                <label>
                    "Accession / captured source ID"
                    <input
                        maxlength="256"
                        prop:value=move ||ref_accession.get()
                        on:input=move |e|ref_accession.set(event_target_value(&e))/>
                </label>
                <label>
                    "Retrieval date"
                    <input
                        type="date"
                        prop:value=move ||ref_date.get()
                        on:input=move |e|ref_date.set(event_target_value(&e))/>
                </label>
                <label>
                    "Reference license"
                    <input
                        maxlength="256"
                        prop:value=move ||ref_license.get()
                        on:input=move |e|ref_license.set(event_target_value(&e))/>
                </label>
                <button class="text-tool" on:click=move |_|add_reference.run(())>
                    "Add pinned reference"
                </button>
                {move || draft.with(|d| d.options.references.iter().map(|reference| {
                    let id = reference.id.clone();
                    let label = format!("{} · {} bases · {}", id, reference.sequence.len(), reference.sequence_revision);
                    view! {<p>{label}<button class="text-tool" on:click=move |_| {
                        let outcome = draft.try_update(|d| d.remove_reference(&id));
                        if let Some(Err(error)) = outcome { status.set(error); }
                    }>"Remove pinned reference"</button></p>}
                }).collect_view())}
                {
                    move || draft.with(|d|d.options.reference_claims.iter().map(|r| {
                        let label = format!("Feature {} ↔ {}", r.feature_index + 1, r.reference_id);
                        let needle = r.clone();
                        let remove = Callback::new(move |(): ()|draft.update(|d| {
                            d.options.reference_claims.retain(|r|r != &needle);
                        }));
                        view! {
                            <p>
                                {label}
                                <button class="text-tool" on:click=move |_|remove.run(())>
                                    "Remove"
                                </button>
                            </p>
                        }
                    }).collect_view())
                }
            </details>
            <details>
                <summary>
                    "Declared fusion register"
                </summary>
                {feature_choice(rows, fusion_cds, "Fusion CDS")}
                {feature_choice(rows, fusion_part, "Component annotation")}
                {phase(fusion_phase, "Expected phase at both component boundaries")}
                <button class="text-tool" on:click=move |_|add_fusion.run(())>
                    "Add fusion check"
                </button>
                {
                    move || draft.with(|d|d.options.fusions.iter().map(|f| {
                        let label = format!("CDS {} · component {} · phase {}", f.cds_feature + 1, f.part_feature + 1,
                        f.expected_offset_mod3);
                        let needle = f.clone();
                        let remove = Callback::new(move |(): ()|draft.update(|d|d.options.fusions.retain(|f|f != &needle)));
                        view! {
                            <p>
                                {label}
                                <button class="text-tool" on:click=move |_|remove.run(())>
                                    "Remove"
                                </button>
                            </p>
                        }
                    }).collect_view())
                }
            </details>
            <details>
                <summary>
                    "Declared current junction"
                </summary>
                <p class="helper">
                    "This checks current register only; it does not infer that an assembly or edit introduced a change."
                </p>
                {feature_choice(rows, junction_cds, "Junction CDS")}
                <label>
                    "Boundary after base (0..source length)"
                    <input
                        inputmode="numeric"
                        prop:value=move ||junction_boundary.get()
                        on:input=move |e|junction_boundary.set(event_target_value(&e))/>
                </label>
                {phase(junction_phase, "Expected junction phase")}
                <button class="text-tool" on:click=move |_|add_junction.run(())>
                    "Add junction check"
                </button>
                {
                    move || draft.with(|d|d.options.junctions.iter().map(|j| {
                        let label = format!("CDS {} · after base {} · phase {}", j.cds_feature + 1, j.boundary, j.expected_offset_mod3);
                        let needle = j.clone();
                        let remove = Callback::new(move |(): ()|draft.update(|d|d.options.junctions.retain(|j|j != &needle)));
                        view! {
                            <p>
                                {label}
                                <button class="text-tool" on:click=move |_|remove.run(())>
                                    "Remove"
                                </button>
                            </p>
                        }
                    }).collect_view())
                }
            </details>
            <super::construct_qc_elements::ElementDeclarations
                draft rows get_molecule source epoch status
            />
            <super::construct_qc_scans::ScanDeclarations draft rows get_molecule status/>
            <button class="button secondary full" disabled=move ||busy.get() on:click=move |_|run.run(())>
                {
                    move || if busy.get() {"Checking…"}
                    else {"Run structural checks"}
                }
            </button>
            <p class="helper" role="status">
                {move || status.get()}
            </p>
            {
                move || report.with(|r|r.as_ref().map(|(capture, r)| {
                    let capture = capture.clone();
                    r.findings.iter().skip(page.get() * 40).take(40).map(|f| {
                        let f = f.clone();
                        let ranges = f.coordinates.iter().flat_map(|l|l.spans().into_iter().map(|(a,
                        b, _)|(a, b))).collect::<Vec<_>>();
                        let selectable = !ranges.is_empty();
                        let c = capture.clone();
                        let click_alive = alive.clone();
                        let click = Callback::new(move |(): ()| {
                            if !click_alive.load(Ordering::Relaxed) {
                                return;
                            }
                            if generation.get_untracked() != c.generation
                                || source.get_untracked() != c .source
                                || epoch.get_untracked() != c.epoch
                                || draft.get_untracked() != c .draft
                                || get_molecule.run(()).revision() != c.source.binding.revision {
                                status.set("This finding is stale. Run the current checks again.".into());
                                return;
                            }
                            on_select.run(ranges.clone());
                        });
                        view! {
                            <details class="qc-finding">
                                <summary>
                                    {format!("{:?} · {:?}", f.verdict, f.rule)}
                                </summary>
                                <p>
                                    {f.reason}
                                </p>
                                <p class="helper">
                                    {
                                        if f.under_declared_annotation_assumptions {"Under your declared annotation assumptions; biological roles are not verified."}
                                        else {"No verified annotation-role premise is established."}
                                    }
                                </p>
                                <button class="text-tool" disabled=(!selectable) on:click=move |_|click.run(())>
                                    "Select coordinates"
                                </button>
                                {
                                    f.citations.into_iter().map(|c| view! {
                                        <p class="helper">
                                            {format!("{} · {}", c.source, c.locator)}
                                        </p>
                                    }).collect_view()
                                }
                            </details>
                        }
                    }).collect_view()
                }))
            }
            <Show when=move ||report.with(Option::is_some)>
                <p class="helper">
                    {
                        move || report.with(|r|r.as_ref().map(|(_, r)|format!("{} findings · page {}", r.findings.len(),
                        page.get() + 1)).unwrap_or_default())
                    }
                </p>
                <button
                    class="text-tool"
                    disabled=move ||page.get()==0
                    on:click=move |_|page.update(|p|*p=p.saturating_sub(1))>
                    "Previous"
                </button>
                <button
                    class="text-tool"
                    disabled=move ||report.with(|r|r.as_ref().is_none_or(|(_,r)|(page.get()+1)*40>=r.findings.len()))
                    on:click=move |_|page.update(|p|*p+=1)>
                    "Next"
                </button>
            </Show>
        </section>
    }
}
