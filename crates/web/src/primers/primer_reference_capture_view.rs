//! Compare an unchanged existing oligo to an explicitly selected reference.
use crate::primers::primer_pair_model::OligoDraft;
use crate::primers::primer_pair_model::PairDraft;
use crate::primers::primer_pair_model::SelectionSeed;
use collagene_core::material::lineage::EventBinding;
use collagene_model::Molecule;
use leptos::prelude::*;
use std::sync::Arc;

#[derive(Clone)]
struct Preview {
    pair: PairDraft,
    seed: SelectionSeed,
    epoch: u64,
    reverse: bool,
    after: OligoDraft,
}
fn slot(pair: &PairDraft, reverse_slot: bool) -> &OligoDraft {
    if reverse_slot {
        &pair.reverse
    } else {
        &pair.forward
    }
}
fn resolve_selection(
    seed: &SelectionSeed,
    resolve_source: Callback<EventBinding, Result<Molecule, String>>,
) -> Result<(), String> {
    let molecule = resolve_source.run(seed.source.clone())?;
    let actual = SelectionSeed::capture(&molecule, &seed.source, seed.start, seed.end)?;
    if actual != *seed {
        return Err("The selection does not match its exact retained source event.".into());
    }
    Ok(())
}

#[component]
pub fn PrimerReferenceCaptureView(
    draft: RwSignal<PairDraft>,
    reverse_slot: bool,
    selection: Callback<(), Result<SelectionSeed, String>>,
    resolve_source: Callback<EventBinding, Result<Molecule, String>>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let staged = RwSignal::new(None::<Arc<Preview>>);
    let status = RwSignal::new(String::new());
    Effect::new(move |_| {
        let pair = draft.get();
        let current_epoch = epoch.get();
        if staged.with_untracked(|value| {
            value
                .as_ref()
                .is_some_and(|value| value.pair != pair || value.epoch != current_epoch)
        }) {
            staged.set(None);
            status.set("The pair or workspace changed; preview the comparison again.".into());
        }
    });
    let preview = Callback::new(move |reverse: bool| {
        let result = (|| -> Result<Preview, String> {
            let pair = draft.get_untracked();
            let captured_epoch = epoch.get_untracked();
            let seed = selection.run(())?;
            resolve_selection(&seed, resolve_source)?;
            let after = seed.compare_existing(slot(&pair, reverse_slot), reverse)?;
            if draft.get_untracked() != pair || epoch.get_untracked() != captured_epoch {
                return Err("The pair or workspace changed; preview again.".into());
            }
            Ok(Preview {
                pair,
                seed,
                epoch: captured_epoch,
                reverse,
                after,
            })
        })();
        match result {
            Ok(value) => {
                staged.set(Some(Arc::new(value)));
                status.set("Preview only. Your oligo text is unchanged; accept to retain this reference comparison.".into());
            }
            Err(error) => {
                staged.set(None);
                status.set(error);
            }
        }
    });
    let accept = move |_| {
        let result = (|| -> Result<(), String> {
            let saved = staged
                .get_untracked()
                .ok_or("Preview a comparison first.")?;
            if draft.get_untracked() != saved.pair || epoch.get_untracked() != saved.epoch {
                return Err("The pair or workspace changed; preview again.".into());
            }
            if selection.run(())? != saved.seed {
                return Err(
                    "The selected source or region changed; preview the intended comparison again."
                        .into(),
                );
            }
            resolve_selection(&saved.seed, resolve_source)?;
            let after = saved
                .seed
                .compare_existing(slot(&saved.pair, reverse_slot), saved.reverse)?;
            if after != saved.after
                || draft.get_untracked() != saved.pair
                || epoch.get_untracked() != saved.epoch
            {
                return Err("The comparison changed; preview again.".into());
            }
            draft.update(|pair| {
                if reverse_slot {
                    pair.reverse = after;
                } else {
                    pair.forward = after;
                }
            });
            staged.set(None);
            Ok(())
        })();
        status.set(match result {
            Ok(()) => "Reference comparison captured. Name, annealing text and tail are unchanged. Apply the pair and run analysis separately.".into(),
            Err(error) => error,
        });
    };
    view! {
        <details class="primer-reference-capture">
            <summary>"Compare existing oligo to selection"</summary>
            <p class="helper">"Select template DNA of the same length as this oligo's annealing domain, then choose its direction. This keeps every literal oligo call, including intentional mismatches. It replaces only the captured reference comparison; it does not search binding sites."</p>
            <div class="primer-selection-actions">
                <button class="button secondary" on:click=move|_|preview.run(false)>"Preview forward comparison"</button>
                <button class="button secondary" on:click=move|_|preview.run(true)>"Preview reverse-complement comparison"</button>
            </div>
            {move||staged.get().map(|saved|{
                let comparison=saved.after.comparison().ok().flatten();
                let source=&saved.seed.source;
                view!{<section class="primer-reference-preview">
                    <h4>"Proposed reference comparison · oligo unchanged"</h4>
                    <p>{format!("Source-forward interval {}–{} (1-based inclusive, unrolled) · {}",saved.seed.start+1,saved.seed.end,if saved.reverse{"reverse-complement oligo"}else{"forward oligo"})}</p>
                    <details><summary>"Exact source event"</summary><p>{format!("Document {} · event {} · revision {}",source.document_id,source.event_id,source.revision)}</p></details>
                    <pre>{format!("Unchanged name: {}\nUnchanged annealing text: {}\nUnchanged 5′ tail: {}",saved.after.name,saved.after.annealing,saved.after.tail)}</pre>
                    {comparison.map(|comparison|view!{<div>
                        <pre>{format!("Selected reference, own 5′→3′: {}\nLiteral domain, own 5′→3′:    {}",comparison.reference,comparison.domain)}</pre>
                        <div class="primer-domain-calls">{comparison.domain.chars().enumerate().map(|(index,call)|{
                            let mismatch=comparison.mismatches.contains(&index);
                            let ambiguous=comparison.ambiguous.contains(&index);
                            view!{<span class:primer-mismatch=mismatch class:primer-ambiguous=ambiguous title=format!("Annealing call {}{}",index+1,if mismatch{" · disjoint mismatch"}else if ambiguous{" · compatible ambiguity"}else{" · exact call"})>{call.to_string()}</span>}
                        }).collect_view()}</div>
                        <p>{format!("{} disjoint mismatches · {} compatible ambiguous comparisons. No binding or amplification result is inferred.",comparison.mismatches.len(),comparison.ambiguous.len())}</p>
                    </div>})}
                    <button class="button primary" on:click=accept>"Accept comparison; keep oligo unchanged"</button>
                    <button class="text-tool" on:click=move|_|staged.set(None)>"Discard comparison preview"</button>
                </section>}
            })}
            <p role="status">{move||status.get()}</p>
        </details>
    }.into_any()
}
