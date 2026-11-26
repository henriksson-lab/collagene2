//! Inspect every returned placement and explicitly copy one into an editable pair.
//! Captured geometry is comparison context, not restored binding authority.
use crate::primer_pair_model::{self, OligoDraft, PairDraft};
use collagene_core::{
    Feature, Molecule, analysis::Primer, oligo_formats::OligoTable, protocol::PrimerEvidence,
    run_lineage::EventBinding,
};
use leptos::prelude::*;
use std::sync::Arc;

const PAGE_SIZE: usize = 50;

#[component]
pub fn PrimerEvidenceView(
    source: EventBinding,
    report: Arc<PrimerEvidence>,
    guard: Callback<(), Result<(), String>>,
    resolve_source: Callback<EventBinding, Result<Molecule, String>>,
    pair_draft: RwSignal<PairDraft>,
    on_select_binding: Callback<(EventBinding, Vec<(usize, usize)>), Result<(), String>>,
    on_annotate: Callback<(EventBinding, Feature), Result<(), String>>,
    on_collection: Callback<OligoTable, Result<(), String>>,
    on_message: Callback<String>,
    on_edit_domain: Callback<()>,
    on_edit_options: Callback<()>,
    on_rerun: Callback<()>,
) -> impl IntoView {
    let primer = match Primer::new(
        &report.analysis.name,
        &report.analysis.annealing,
        &report.analysis.tail,
    ) {
        Ok(primer)
            if report.bindings.source_revision == source.revision
                && report.bindings.primer_name == primer.name
                && report.bindings.annealing == primer.annealing.as_str()
                && report.analysis.ordered_sequence == primer.ordered_sequence() =>
        {
            primer
        }
        _ => {
            return view! {<p role="alert">"Primer evidence does not match its source or oligo. Run analysis again."</p>}.into_any();
        }
    };
    let count = report.bindings.hits.len();
    let page_count = count.div_ceil(PAGE_SIZE).max(1);
    let source = StoredValue::new(source);
    let primer = StoredValue::new(primer);
    let report = StoredValue::new(report);
    let page = RwSignal::new(0usize);
    let page_input = RwSignal::new("1".to_string());
    let move_page = Callback::new(move |next: usize| {
        page.set(next.min(page_count - 1));
        page_input.set((page.get_untracked() + 1).to_string());
    });
    let check = Callback::new(move |()| -> Result<Molecule, String> {
        guard.run(())?;
        let molecule = resolve_source.run(source.get_value())?;
        if molecule.revision() != source.with_value(|s| s.revision.clone()) {
            return Err("The retained source no longer matches this binding report.".into());
        }
        Ok(molecule)
    });
    let copy_unanchored = Callback::new(move |reverse: bool| {
        let result = guard.run(()).map(|()| {
            let value = primer.with_value(|p| OligoDraft {
                name: p.name.clone(),
                annealing: p.annealing.as_str().into(),
                tail: p.tail.clone(),
                anchor: None,
            });
            pair_draft.update(|pair| {
                if reverse {
                    pair.reverse = value;
                } else {
                    pair.forward = value;
                }
            });
        });
        on_message.run(match result {
            Ok(()) => "Copied the oligo without a selected placement. Apply the pair separately, then run PCR analysis.".into(),
            Err(error) => error,
        });
    });
    let use_hit = Callback::new(move |(index, reverse_slot): (usize, bool)| {
        let result = (|| -> Result<(), String> {
            let molecule = check.run(())?;
            let value = report.with_value(|report| {
                let hit = report
                    .bindings
                    .hits
                    .get(index)
                    .ok_or("Binding hit is unavailable.")?;
                primer.with_value(|primer| {
                    primer_pair_model::capture_binding_hit(
                        &molecule,
                        &source.get_value(),
                        primer,
                        hit,
                    )
                })
            })?;
            // Resolvers are host callbacks: recheck before changing either slot.
            guard.run(())?;
            pair_draft.update(|pair| {
                if reverse_slot {
                    pair.reverse = value;
                } else {
                    pair.forward = value;
                }
            });
            Ok(())
        })();
        on_message.run(match result {
            Ok(()) => format!("Copied binding {} into the {} draft slot with its exact source comparison. Oligo sequence and the other slot were preserved. Apply the pair separately and rerun PCR; other placements and incomplete enumeration remain relevant.",index+1,if reverse_slot{"reverse"}else{"forward"}),
            Err(error) => error,
        });
    });
    let navigate = Callback::new(move |index: usize| {
        let result = (|| -> Result<(), String> {
            let molecule = check.run(())?;
            let spans = report.with_value(|report| {
                let hit = report
                    .bindings
                    .hits
                    .get(index)
                    .ok_or("Binding hit is unavailable.")?;
                primer.with_value(|primer| {
                    primer_pair_model::capture_binding_hit(
                        &molecule,
                        &source.get_value(),
                        primer,
                        hit,
                    )
                })?;
                Ok::<_, String>(
                    hit.location
                        .spans()
                        .into_iter()
                        .map(|(a, b, _)| (a, b))
                        .collect::<Vec<_>>(),
                )
            })?;
            guard.run(())?;
            on_select_binding.run((source.get_value(), spans))
        })();
        on_message.run(match result {
            Ok(()) => format!("Selected every source span of binding {}. An origin-spanning hit retains both parts; the captured pair anchor retains the complete domain.",index+1),
            Err(error) => error,
        });
    });
    let save_hit = Callback::new(move |index: usize| {
        let result = (|| -> Result<(), String> {
            let molecule = check.run(())?;
            let feature = report.with_value(|report| {
                let hit = report
                    .bindings
                    .hits
                    .get(index)
                    .ok_or("Binding hit is unavailable.")?;
                primer.with_value(|primer| {
                    primer_pair_model::capture_binding_hit(
                        &molecule,
                        &source.get_value(),
                        primer,
                        hit,
                    )?;
                    Ok::<_, String>(crate::primers::annotation(
                        primer,
                        hit.location.clone(),
                        &report.analysis,
                        !hit.certain,
                        &source.with_value(|s| s.revision.clone()),
                    ))
                })
            })?;
            guard.run(())?;
            on_annotate.run((source.get_value(), feature))
        })();
        on_message.run(match result {
            Ok(()) => "Saved the primer annotation with its exact source revision. The source changed; rerun analysis before saving another placement.".into(),
            Err(error) => error,
        });
    });
    let export_alive = RwSignal::new(());
    let download = move |_| {
        let result = guard
            .run(())
            .and_then(|()| {
                report.with_value(|report| {
                    serde_json::to_vec_pretty(report.as_ref()).map_err(|e| e.to_string())
                })
            })
            .and_then(|bytes| {
                crate::app::prepare_download("primer-evidence.json", "application/json", &bytes)
            });
        match result {
            Ok(artifact) => {
                let check = move || {
                    if export_alive.try_get_untracked().is_none() {
                        return Err("Primer review closed.".into());
                    }
                    guard.run(())
                };
                artifact.complete_guarded(check, move |result| {
                    if export_alive.try_get_untracked().is_none() || guard.run(()).is_err() {
                        return;
                    }
                    if let Err(error) = result {
                        on_message.run(error);
                    }
                });
            }
            Err(error) => on_message.run(error),
        }
    };
    let metadata = report.with_value(|r| {
        let a = &r.analysis;
        let tm =
            a.tm.as_ref()
                .map(|t| format!("Tm {:.2} °C · {}", t.tm_c, t.model))
                .unwrap_or_else(|| {
                    a.tm_unavailable
                        .clone()
                        .unwrap_or_else(|| "Tm unavailable".into())
                });
        (
            a.ordered_sequence.clone(),
            format!(
                "Annealing {} nt · tail {} nt · GC {:.1}–{:.1}%",
                a.annealing.len(),
                a.tail.len(),
                a.annealing_gc.minimum_fraction * 100.,
                a.annealing_gc.maximum_fraction * 100.
            ),
            tm,
            format!(
                "Contiguous possible runs: hairpin {} · homodimer {} · homodimer at a 3′ end {}",
                a.hairpin.max_run, a.homodimer.max_run, a.homodimer.max_three_prime_run
            ),
            a.warnings
                .iter()
                .chain(&r.bindings.warnings)
                .cloned()
                .collect::<Vec<_>>(),
            r.bindings.truncated,
        )
    });
    view! {<section class="primer-evidence"><h3>"Oligo evidence"</h3>
        <pre>{metadata.0}</pre>
        <details><summary>"Exact binding source"</summary><p>{source.with_value(|s|format!("Document {} · event {} · revision {}",s.document_id,s.event_id,s.revision))}</p></details>
        <crate::primer_orders::PrimerOrdersPanel primers=vec![primer.get_value()] guard on_collection/>
        <div class="primer-actions"><button class="button secondary" on:click=move|_|copy_unanchored.run(false)>"Keep as forward primer without placement"</button><button class="button secondary" on:click=move|_|copy_unanchored.run(true)>"Keep as reverse primer without placement"</button></div>
        <p>{metadata.1}</p><p class="primer-tm">{metadata.2}</p><p>{metadata.3}</p>
        <ul>{metadata.4.into_iter().map(|warning|view!{<li>{warning}</li>}).collect_view()}</ul>
        <h4>{format!("{count} returned binding hits{}",if metadata.5{" · incomplete enumeration"}else{""})}</h4>
        <p class="helper">{if count==0{"No placements were returned under the declared mismatch/3′ policy. Inspect the domain and binding policy, then run again; the software will not relax these settings automatically."}else if count>1{"Multiple placements were returned. Inspect an explicit placement before copying it to a pair slot; choosing one does not establish uniqueness."}else{"One placement was returned. Inspect its strand and diagnostic calls; a single returned row does not establish specificity or amplification."}}</p>
        {metadata.5.then(||view!{<p class="primer-anchor-warning">"The search was incomplete. Missing placements remain unknown; choosing a returned hit does not complete enumeration."</p>})}
        <p class="helper">{report.with_value(|r|format!("Declared policy: at most {} disjoint mismatches, at least {} exact 3′ bases; {} placements examined. Synthetic tails are excluded from this template search.",r.bindings.options.max_mismatches,r.bindings.options.min_three_prime_exact,r.bindings.placements_examined))}</p>
        <div class="primer-actions"><button class="button secondary" on:click=move|_|on_edit_domain.run(())>"Edit domain or tail"</button><button class="button secondary" on:click=move|_|on_edit_options.run(())>"Review binding policy"</button><button class="button secondary" on:click=move|_|on_rerun.run(())>"Run analysis again"</button></div>
        <Show when=move||{count>PAGE_SIZE}><div class="primer-actions" aria-label="Binding hit pages">
            <button class="button secondary" disabled=move||page.get()==0 on:click=move|_|move_page.run(page.get_untracked().saturating_sub(1))>"Previous hits"</button>
            <span>{move||format!("Page {} of {}",page.get()+1,page_count)}</span>
            <button class="button secondary" disabled=move||{page.get()+1>=page_count} on:click=move|_|move_page.run(page.get_untracked()+1)>"Next hits"</button>
            <label>"Page"<input aria-label="Binding hit page" type="number" min="1" max=page_count prop:value=move||page_input.get() on:input=move|event|page_input.set(event_target_value(&event))/></label>
            <button class="button secondary" on:click=move|_|{match page_input.get_untracked().parse::<usize>(){Ok(value) if value>0&&value<=page_count=>move_page.run(value-1),_=>on_message.run(format!("Enter a binding page from 1 to {page_count}."))}}>"Go to page"</button>
        </div></Show>
        {move||{
            let first=page.get()*PAGE_SIZE;
            report.with_value(|report|report.bindings.hits.iter().enumerate().skip(first).take(PAGE_SIZE).map(|(index,hit)|{
                let spans=hit.location.spans();
                let offsets=|values:&[usize]|if values.is_empty(){"none".to_string()}else{values.iter().map(|v|(v+1).to_string()).collect::<Vec<_>>().join(", ")};
                view!{<article class="primer-hit"><strong>{format!("Binding {} · {}–{}{}{}",index+1,hit.start+1,hit.start+hit.length,if hit.reverse{" · reverse"}else{" · forward"},if hit.crosses_origin{" · crosses origin"}else{""})}</strong>
                    <p>{format!("{} mismatches · {} ambiguous positions · {} exact 3′ bases{}",hit.mismatches.len(),hit.ambiguous_positions.len(),hit.three_prime_exact,if hit.certain{" · certain sequence comparison"}else{" · conditional/ambiguous sequence comparison"})}</p>
                    <p class="helper">{format!("Own 5′→3′ annealing-call offsets (1-based): disjoint mismatches {}; compatible ambiguity {}.",offsets(&hit.mismatches),offsets(&hit.ambiguous_positions))}</p>
                    <p class="helper">{format!("Source spans: {}",spans.iter().map(|(a,b,_)|format!("{}–{}",a+1,b)).collect::<Vec<_>>().join(", "))}</p>
                    <div class="primer-actions"><button class="button secondary" on:click=move|_|navigate.run(index)>"Select binding"</button>
                        <button class="button secondary" on:click=move|_|use_hit.run((index,false))>"Use binding as forward primer"</button>
                        <button class="button secondary" on:click=move|_|use_hit.run((index,true))>"Use binding as reverse primer"</button>
                        <button class="button secondary" on:click=move|_|save_hit.run(index)>"Save primer annotation"</button>
                    </div>
                </article>}
            }).collect_view())
        }}
        <p class="helper">"Every returned hit is reachable through the pages. Copying a hit updates one editable pair slot and retains a source comparison, not a verified PCR result. The slot name does not change the selected hit's strand. Apply the pair separately; all alternative placements and search limits still apply."</p>
        <button class="button secondary full" on:click=download>"Download complete primer evidence"</button>
    </section>}.into_any()
}
