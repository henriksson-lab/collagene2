//! Shared explicit primer inputs. Local comparisons describe only the captured
//! selection; amplification, specificity and thermodynamics require fresh analysis.
use crate::primers::primer_pair_model::OligoDraft;
use crate::primers::primer_pair_model::PairDraft;
use crate::primers::primer_pair_model::SelectionSeed;
use crate::primers::primer_picker::PrimerPicker;
use crate::primers::primer_picker::ResolveRequest;
use crate::sources::identity::SourceInfo;
use collagene_core::{
    cloning::pcr::reuse as primer_reuse, material::lineage::EventBinding,
    sequence::analysis::Primer,
};
use collagene_model::Molecule;
use leptos::prelude::*;

fn slot(pair: &PairDraft, reverse_slot: bool) -> &OligoDraft {
    if reverse_slot {
        &pair.reverse
    } else {
        &pair.forward
    }
}
fn slot_mut(pair: &mut PairDraft, reverse_slot: bool) -> &mut OligoDraft {
    if reverse_slot {
        &mut pair.reverse
    } else {
        &mut pair.forward
    }
}
fn positions(values: &[usize]) -> String {
    if values.is_empty() {
        "none".into()
    } else {
        values
            .iter()
            .take(200)
            .map(|n| (n + 1).to_string())
            .collect::<Vec<_>>()
            .join(", ")
    }
}

#[component]
pub fn PrimerPairEditor(
    draft: RwSignal<PairDraft>,
    selection: Callback<(), Result<SelectionSeed, String>>,
    anchor_available: Callback<EventBinding, bool>,
    resolve_source: Callback<EventBinding, Result<Molecule, String>>,
    primer_sources: Memo<Vec<SourceInfo>>,
    list_primers: Callback<EventBinding, Result<primer_reuse::List, String>>,
    resolve_primer: Callback<ResolveRequest, Result<primer_reuse::Resolved, String>>,
    epoch: RwSignal<u64>,
    on_apply: Callback<(Primer, Primer), Result<(), String>>,
    apply_label: String,
) -> impl IntoView {
    let status = RwSignal::new(String::new());
    let apply = move |_| {
        let result = draft
            .with_untracked(PairDraft::primers)
            .and_then(|pair| on_apply.run(pair));
        match result {
            Ok(()) => status.set(
                "Primer pair applied. Run the relevant analysis again before using its result."
                    .into(),
            ),
            Err(error) => status.set(error),
        }
    };
    view! {<section class="primer-pair-editor" aria-label="Primer pair editor">
        <style>".primer-pair-editor .primer-pair-slots{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,280px),1fr));gap:12px}.primer-pair-editor fieldset{min-width:0;border:1px solid #dce3e8;padding:10px;border-radius:5px}.primer-pair-editor label{display:block;margin:8px 0}.primer-pair-editor input,.primer-pair-editor textarea{box-sizing:border-box;width:100%}.primer-pair-editor pre{white-space:pre-wrap;overflow-wrap:anywhere;font:12px/1.6 monospace;max-height:240px;overflow:auto}.primer-pair-editor .primer-selection-actions{display:flex;flex-wrap:wrap;gap:6px}.primer-pair-editor .primer-mismatch{background:#ffe0df;color:#8b1712}.primer-pair-editor .primer-ambiguous{background:#fff0bf;color:#6a4700}.primer-pair-editor .primer-domain-calls{white-space:pre-wrap;overflow-wrap:anywhere;font:12px/1.7 monospace}.primer-pair-editor details{overflow-wrap:anywhere}.primer-pair-editor .primer-anchor-warning{color:#8b4d24}"</style>
        <p class="helper">"This editable pair is saved with your workspace. Choose a selection or a saved primer explicitly, or enter the sequences. Every oligo is written 5′→3′; the 5′ tail is separate from the annealing domain. Selection direction is independent of which primer slot you are editing."</p>
        <div class="primer-pair-slots">
            <PrimerSlotEditor draft reverse_slot=false selection anchor_available resolve_source primer_sources list_primers resolve_primer epoch/>
            <PrimerSlotEditor draft reverse_slot=true selection anchor_available resolve_source primer_sources list_primers resolve_primer epoch/>
        </div>
        <p class="helper">"A selection comparison is not a binding-site search, an annealing-temperature calculation or proof of PCR. Manually edited mismatches remain explicit design inputs; changed domain lengths leave the captured comparison unresolved."</p>
        <button class="button primary" on:click=apply>{apply_label}</button>
        <p role="status">{move||status.get()}</p>
    </section>}.into_any()
}

#[component]
fn PrimerSlotEditor(
    draft: RwSignal<PairDraft>,
    reverse_slot: bool,
    selection: Callback<(), Result<SelectionSeed, String>>,
    anchor_available: Callback<EventBinding, bool>,
    resolve_source: Callback<EventBinding, Result<Molecule, String>>,
    primer_sources: Memo<Vec<SourceInfo>>,
    list_primers: Callback<EventBinding, Result<primer_reuse::List, String>>,
    resolve_primer: Callback<ResolveRequest, Result<primer_reuse::Resolved, String>>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let label = if reverse_slot { "Reverse" } else { "Forward" };
    let status = RwSignal::new(String::new());
    let oligo = Memo::new(move |_| draft.with(|pair| slot(pair, reverse_slot).clone()));
    let comparison = Memo::new(move |_| oligo.with(OligoDraft::comparison));
    let available = Memo::new(move |_| {
        let _ = epoch.get();
        oligo.with(|oligo| {
            oligo
                .anchor
                .as_ref()
                .is_none_or(|anchor| anchor_available.run(anchor.source.clone()))
        })
    });
    let use_selection = Callback::new(move |reverse: bool| {
        let captured_epoch = epoch.get_untracked();
        let result = (|| -> Result<OligoDraft, String> {
            let seed = selection.run(())?;
            if epoch.get_untracked() != captured_epoch || !anchor_available.run(seed.source.clone())
            {
                return Err("The selected source changed or is unavailable. Select its current region again.".into());
            }
            draft.with_untracked(|pair| {
                let current = slot(pair, reverse_slot);
                seed.oligo(&current.name, reverse, &current.tail)
            })
        })();
        match result {
            Ok(value) => {
                draft.update(|pair| *slot_mut(pair, reverse_slot) = value);
                status.set("Copied the selected DNA in the requested direction. Name and 5′ tail were preserved; review before applying.".into());
            }
            Err(error) => status.set(error),
        }
    });
    let on_copy = Callback::new(move |primer: Primer| {
        let value = OligoDraft {
            name: primer.name,
            annealing: primer.annealing.as_str().to_owned(),
            tail: primer.tail,
            anchor: None,
        };
        draft.update(|pair| *slot_mut(pair, reverse_slot) = value);
        status.set("Copied the saved primer's name, annealing domain and tail. It has no current sequence-selection anchor.".into());
    });
    view! {<fieldset class="primer-pair-slot"><legend>{format!("{label} primer")}</legend>
        <label>"Name"<input aria-label=format!("{label} primer name") maxlength="256" prop:value=move||oligo.with(|o|o.name.clone()) on:input=move|e|draft.update(|pair|slot_mut(pair,reverse_slot).name=event_target_value(&e))/></label>
        <label>"Annealing domain (5′→3′)"<textarea aria-label=format!("{label} primer annealing domain") rows="3" maxlength="1600" spellcheck="false" prop:value=move||oligo.with(|o|o.annealing.clone()) on:input=move|e|draft.update(|pair|slot_mut(pair,reverse_slot).annealing=event_target_value(&e))></textarea></label>
        <label>"5′ tail"<textarea aria-label=format!("{label} primer 5 prime tail") rows="2" maxlength="1600" spellcheck="false" prop:value=move||oligo.with(|o|o.tail.clone()) on:input=move|e|draft.update(|pair|slot_mut(pair,reverse_slot).tail=event_target_value(&e))></textarea></label>
        <div class="primer-selection-actions"><button class="button secondary" on:click=move|_|use_selection.run(false)>"Use selection forward"</button><button class="button secondary" on:click=move|_|use_selection.run(true)>"Use selection reverse complement"</button></div>
        <crate::primers::primer_reference_capture_view::PrimerReferenceCaptureView draft reverse_slot selection resolve_source epoch/>
        <PrimerPicker label=format!("{label} pair editor") sources=primer_sources list_primers resolve_primer epoch on_copy/>
        <pre>{move||oligo.with(|o|format!("Ordered oligo, 5′→3′: {}{}\nAnnealing domain: {}\n5′ tail: {}",o.tail,o.annealing,o.annealing,if o.tail.is_empty(){"(none)"}else{&o.tail}))}</pre>
        {move||oligo.with(|o|o.anchor.clone()).map(|anchor|view!{<section class="primer-anchor-comparison">
            <p>{format!("Captured selection {}–{} (1-based inclusive; end may continue past circular origin) · {}",anchor.start.saturating_add(1),anchor.end,if anchor.reverse{"reverse complement"}else{"forward"})}</p>
            <Show when=move||!available.get()><p class="primer-anchor-warning">"The captured source event is unavailable. Your edits and captured reference are retained; this is not evidence against another source."</p></Show>
            <details><summary>"Source details"</summary><p>{format!("Document {} · event {} · revision {}",anchor.source.document_id,anchor.source.event_id,anchor.source.revision)}</p></details>
            {move||{match comparison.get(){
                Ok(Some(value)) => view!{<div><pre>{format!("Captured reference, own 5′→3′: {}\nEdited annealing domain, own 5′→3′: {}",value.reference,value.domain)}</pre>
                    <div class="primer-domain-calls" aria-label=format!("{label} primer compared bases")>{value.domain.chars().enumerate().map(|(index,call)|{let mismatch=value.mismatches.contains(&index);let ambiguous=value.ambiguous.contains(&index);view!{<span class:primer-mismatch=mismatch class:primer-ambiguous=ambiguous title=format!("Annealing call {}{}",index+1,if mismatch{" · mismatch"}else if ambiguous{" · ambiguous comparison"}else{" · exact call"})>{call.to_string()}</span>}}).collect_view()}</div>
                    <p>{format!("Mismatches at annealing calls: {}. Ambiguous comparisons at calls: {}. Positions are 1-based within this oligo's annealing domain and exclude its tail.",positions(&value.mismatches),positions(&value.ambiguous))}</p>
                </div>}.into_any(),
                Ok(None) => view!{<p class="helper">"No captured selection comparison."</p>}.into_any(),
                Err(error) => view!{<p class="primer-anchor-warning">{format!("Selection comparison unresolved: {error}")}</p>}.into_any(),
            }}}
            <button class="text-tool" on:click=move|_|draft.update(|pair|slot_mut(pair,reverse_slot).anchor=None)>"Detach selection comparison"</button>
        </section>})}
        <crate::primers::primer_domain_edit_view::PrimerDomainEditView draft reverse_slot resolve_source epoch/>
        <p role="status">{move||status.get()}</p>
    </fieldset>}.into_any()
}
