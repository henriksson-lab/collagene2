//! Explicit previews modify draft geometry only; pair application remains separate.
use crate::primers::primer_domain_edit;
use crate::primers::primer_domain_edit::Edit;
use crate::primers::primer_domain_edit::Preview;
use crate::primers::primer_pair_model::OligoDraft;
use crate::primers::primer_pair_model::PairDraft;
use collagene_core::material::lineage::EventBinding;
use collagene_model::Molecule;
use leptos::prelude::*;
use std::sync::Arc;

#[derive(Clone)]
struct Staged {
    pair: PairDraft,
    source: EventBinding,
    epoch: u64,
    source_length: usize,
    fields: (String, String, String),
    preview: Preview,
}
fn oligo(pair: &PairDraft, reverse_slot: bool) -> &OligoDraft {
    if reverse_slot {
        &pair.reverse
    } else {
        &pair.forward
    }
}
fn positions(values: &[usize]) -> String {
    if values.is_empty() {
        "none".into()
    } else {
        values
            .iter()
            .map(|n| (n + 1).to_string())
            .collect::<Vec<_>>()
            .join(", ")
    }
}
fn range_label(oligo: &OligoDraft, source_length: usize) -> String {
    let Some(anchor) = &oligo.anchor else {
        return "No captured source anchor".into();
    };
    let spans = if anchor.end > source_length {
        format!(
            "{}–{}, 1–{}",
            anchor.start + 1,
            source_length,
            anchor.end - source_length
        )
    } else {
        format!("{}–{}", anchor.start + 1, anchor.end)
    };
    let three_prime = if anchor.reverse {
        anchor.start + 1
    } else {
        (anchor.end - 1) % source_length + 1
    };
    format!(
        "Source-forward interval {}–{} (1-based inclusive, unrolled) · displayed spans {} · {} · 3′ template base {}",
        anchor.start + 1,
        anchor.end,
        spans,
        if anchor.reverse {
            "reverse-complement oligo"
        } else {
            "forward oligo"
        },
        three_prime
    )
}

#[component]
pub fn PrimerDomainEditView(
    draft: RwSignal<PairDraft>,
    reverse_slot: bool,
    resolve_source: Callback<EventBinding, Result<Molecule, String>>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let start = RwSignal::new(String::new());
    let end = RwSignal::new(String::new());
    let tail_length = RwSignal::new(String::new());
    let staged = RwSignal::new(None::<Arc<Staged>>);
    let status = RwSignal::new(String::new());
    let current_oligo = Memo::new(move |_| draft.with(|pair| oligo(pair, reverse_slot).clone()));
    Effect::new(move |_| {
        let current = current_oligo.get();
        start.set(
            current
                .anchor
                .as_ref()
                .map(|a| a.start.saturating_add(1).to_string())
                .unwrap_or_default(),
        );
        end.set(
            current
                .anchor
                .as_ref()
                .map(|a| a.end.to_string())
                .unwrap_or_default(),
        );
        tail_length.set(
            current
                .primer()
                .map(|p| p.tail.len().to_string())
                .unwrap_or_default(),
        );
    });
    Effect::new(move |_| {
        let pair = draft.get();
        let current_epoch = epoch.get();
        if staged.with_untracked(|saved| {
            saved
                .as_ref()
                .is_some_and(|saved| saved.pair != pair || saved.epoch != current_epoch)
        }) {
            staged.set(None);
            status.set("The pair or workspace changed; create a fresh preview.".into());
        }
    });
    let fields = Callback::new(move |()| {
        (
            start.get_untracked(),
            end.get_untracked(),
            tail_length.get_untracked(),
        )
    });
    let create_preview = Callback::new(move |boundary: bool| {
        let result = (|| -> Result<Arc<Staged>, String> {
            let pair = draft.get_untracked();
            let current = oligo(&pair, reverse_slot);
            let anchor = current.anchor.as_ref().ok_or(
                "Capture a source selection before editing its endpoints or tail boundary.",
            )?;
            let captured_epoch = epoch.get_untracked();
            let source = anchor.source.clone();
            let molecule = resolve_source.run(source.clone())?;
            let form = fields.run(());
            let edit = if boundary {
                Edit::TailBoundary {
                    tail_length: form
                        .2
                        .parse::<usize>()
                        .map_err(|_| "Enter a nonnegative tail length in ordered 5′→3′ bases.")?,
                }
            } else {
                let first = form
                    .0
                    .parse::<usize>()
                    .map_err(|_| "Enter a 1-based source start.")?;
                let start = first
                    .checked_sub(1)
                    .ok_or("Source start must be at least 1.")?;
                let end = form
                    .1
                    .parse::<usize>()
                    .map_err(|_| "Enter a 1-based inclusive, possibly unrolled end.")?;
                Edit::Endpoints { start, end }
            };
            let preview =
                primer_domain_edit::preview(&molecule, &source, captured_epoch, current, edit)?;
            if draft.get_untracked() != pair || epoch.get_untracked() != captured_epoch {
                return Err("The pair or workspace changed during preview.".into());
            }
            Ok(Arc::new(Staged {
                pair,
                source,
                epoch: captured_epoch,
                source_length: molecule.sequence.len(),
                fields: form,
                preview,
            }))
        })();
        match result {
            Ok(preview) => {
                staged.set(Some(preview));
                status.set("Preview only. Accept updates this draft slot; apply the pair separately before rerunning analysis.".into());
            }
            Err(error) => {
                staged.set(None);
                status.set(error);
            }
        }
    });
    let preview_current = Memo::new(move |_| {
        staged.with(|saved| {
            saved.as_ref().is_some_and(|saved| {
                saved.pair == draft.get()
                    && saved.epoch == epoch.get()
                    && saved.fields == (start.get(), end.get(), tail_length.get())
            })
        })
    });
    let accept_preview = Callback::new(move |()| {
        let result = (|| -> Result<(), String> {
            let saved = staged
                .get_untracked()
                .ok_or("Create an edit preview first.")?;
            if draft.get_untracked() != saved.pair
                || epoch.get_untracked() != saved.epoch
                || fields.run(()) != saved.fields
            {
                return Err(
                    "The pair, endpoints, boundary or workspace changed; preview again.".into(),
                );
            }
            let molecule = resolve_source.run(saved.source.clone())?;
            let current = oligo(&saved.pair, reverse_slot);
            let edited = primer_domain_edit::accept(
                &saved.preview,
                &molecule,
                &saved.source,
                saved.epoch,
                current,
            )?;
            if draft.get_untracked() != saved.pair || epoch.get_untracked() != saved.epoch {
                return Err("The pair or workspace changed; this preview was not accepted.".into());
            }
            draft.update(|pair| {
                if reverse_slot {
                    pair.reverse = edited
                } else {
                    pair.forward = edited
                }
            });
            staged.set(None);
            Ok(())
        })();
        status.set(match result{Ok(())=>"Accepted this slot's draft edit. Apply the pair explicitly and rerun the workflow before using any analysis result.".into(),Err(error)=>error});
    });
    let label = if reverse_slot {
        "Reverse slot"
    } else {
        "Forward slot"
    };
    view! {
        <details class="primer-domain-edit"><summary>"Preview endpoints and tail boundary"</summary>
            <p class="helper">"Edits use this exact retained source version, which may be historical. Applying the pair reruns the selected workflow and does not rebind the anchor to another source."</p>
            <p class="helper">{move||current_oligo.with(|current|current.anchor.as_ref().map(|a|format!("Captured strand: {}. Endpoint edits keep this direction; use a new source selection to change it.",if a.reverse{"reverse complement"}else{"forward"})).unwrap_or_else(||"Capture an exact source selection first. Freeform text without an anchor remains editable above.".into()))}</p>
            <label>"Source start (1-based inclusive)"<input type="number" min="1" aria-label=format!("{label} edit source start") prop:value=move||start.get() on:input=move|event|{start.set(event_target_value(&event));}/></label>
            <label>"Source end (1-based inclusive, unrolled)"<input type="number" min="1" aria-label=format!("{label} edit source end") prop:value=move||end.get() on:input=move|event|{end.set(event_target_value(&event));}/></label>
            <p class="helper">"Coordinates follow the forward reference, including for reverse primers. On circular DNA, start is within the source and end may continue beyond its length across the origin; at most one turn is allowed. Overlapping literal domain calls are preserved; extensions use actual source calls."</p>
            <button class="button secondary" on:click=move|_|{create_preview.run(false);}>"Preview endpoint change"</button>
            <label>"5′ tail length in ordered oligo bases"<input type="number" min="0" max="199" aria-label=format!("{label} edit tail boundary") prop:value=move||tail_length.get() on:input=move|event|{tail_length.set(event_target_value(&event));}/></label>
            <p class="helper">"Moving the tail boundary preserves the full ordered 5′→3′ oligo. Bases promoted into annealing remain literal calls and are compared against their newly captured template positions. At least one annealing base must remain."</p>
            <button class="button secondary" on:click=move|_|{create_preview.run(true);}>"Preview tail-boundary change"</button>
            {move||staged.get().map(|saved|{
                let plan=&saved.preview;
                let before=plan.before.primer().map(|p|p.ordered_sequence()).unwrap_or_else(|e|e);
                let after=plan.after.primer().map(|p|p.ordered_sequence()).unwrap_or_else(|e|e);
                let comparison=plan.comparison.clone();
                view!{<section class="primer-domain-preview"><h4>"Proposed draft edit · not yet accepted"</h4>
                    <details><summary>"Exact source event"</summary><p>{format!("Document {} · event {} · revision {} · workspace epoch {}",saved.source.document_id,saved.source.event_id,saved.source.revision,saved.epoch)}</p></details>
                    <p>{range_label(&plan.after,saved.source_length)}</p>
                    <pre>{format!("Before ordered 5′→3′: {before}\nAfter ordered 5′→3′:  {after}\n\nBefore tail: {}\nAfter tail:  {}\nBefore annealing: {}\nAfter annealing:  {}",if plan.before.tail.is_empty(){"(none)"}else{&plan.before.tail},if plan.after.tail.is_empty(){"(none)"}else{&plan.after.tail},plan.before.annealing,plan.after.annealing)}</pre>
                    <pre>{format!("New captured template, own 5′→3′: {}\nLiteral annealing calls, own 5′→3′: {}",comparison.reference,comparison.domain)}</pre>
                    <div class="primer-domain-calls" aria-label=format!("{label} proposed annealing comparisons")>{comparison.domain.chars().enumerate().map(|(index,call)|{
                        let mismatch=comparison.mismatches.contains(&index);let ambiguous=comparison.ambiguous.contains(&index);
                        view!{<span class:primer-mismatch=mismatch class:primer-ambiguous=ambiguous title=format!("Annealing call {}{}",index+1,if mismatch{" · disjoint mismatch"}else if ambiguous{" · compatible ambiguity"}else{" · exact call"})>{call.to_string()}</span>}
                    }).collect_view()}</div>
                    <p class="helper">{format!("Disjoint mismatches at domain calls: {}. Compatible ambiguous comparisons: {}. Domain indices are one-based and exclude the tail.",positions(&comparison.mismatches),positions(&comparison.ambiguous))}</p>
                    <ul>{plan.removed_edits.iter().map(|edit|view!{<li>{format!("Edited call leaving annealing: template base {} · old domain call {} · reference {} → literal {}{}{}",edit.source_position+1,edit.old_domain_index+1,edit.reference,edit.literal,if edit.mismatch{" · mismatch"}else{""},if edit.ambiguous{" · ambiguity"}else{""})}</li>}).collect_view()}</ul>
                    {plan.warnings.iter().map(|warning|view!{<p class="primer-anchor-warning">{warning.clone()}</p>}).collect_view()}
                    <p class="helper">"A new captured comparison describes source geometry and literal sequence only. It does not establish binding specificity, amplification or thermodynamic suitability. Calls leaving annealing may remain in the tail after a boundary move; inspect the ordered sequences."</p>
                    <Show when=move||{!preview_current.get()}><p class="primer-anchor-warning">"Inputs changed after this preview. Create a new preview before accepting."</p></Show>
                    <button class="button primary" disabled=move||{!preview_current.get()} on:click=move|_|{accept_preview.run(());}>"Accept preview into this draft slot"</button>
                    <button class="text-tool" on:click=move|_|{staged.set(None);}>"Discard preview"</button>
                </section>}
            })}
            <p class="helper" role="status">{move||status.get()}</p>
        </details>
    }.into_any()
}
