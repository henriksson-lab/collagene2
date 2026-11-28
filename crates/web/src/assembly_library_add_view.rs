//! Explicit source-bin additions; immutable proposal review belongs to the host.
use crate::{
    assembly_draft_model::{CapturedSource, Draft, Slot},
    assembly_library_topology::{self as topology, ChoiceSeed},
    cloning_methods_model::OverlapInput,
    digest::SourceInfo,
};
use collagene_core::run_lineage::EventBinding;
use collagene_model::{Alphabet, Molecule};
use leptos::prelude::*;
use std::sync::Arc;
const PAGE: usize = 10;
struct SeedBudget(usize);
impl std::io::Write for SeedBudget {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|n| *n <= 32 << 20)
            .ok_or_else(|| std::io::Error::other("new source choices exceed32MiB"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn seed_name(seed: &ChoiceSeed) -> String {
    let source = seed
        .slot
        .captured
        .as_ref()
        .map(|v| v.source.name.as_str())
        .unwrap_or("Unassigned");
    format!(
        "{} · {} · {}",
        source,
        crate::assembly_draft_model::selector_label(&seed.slot.selector),
        seed.preparation.summary()
    )
}
#[component]
pub fn AddLibraryBinView(
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    busy: RwSignal<bool>,
    sources: Memo<Vec<SourceInfo>>,
    get_source: Callback<EventBinding, Result<Molecule, String>>,
    on_propose: Callback<(Arc<Draft>, u64, topology::Change)>,
) -> impl IntoView {
    let base = RwSignal::new(None::<(Arc<Draft>, u64)>);
    let label = RwSignal::new(String::new());
    let position = RwSignal::new(String::new());
    let seeds = RwSignal::new(Vec::<ChoiceSeed>::new());
    let assigned = RwSignal::new(Vec::<Option<usize>>::new());
    let substrates = RwSignal::new(Vec::<Option<u8>>::new());
    let row_page = RwSignal::new(0usize);
    let choice_page = RwSignal::new(0usize);
    let source_index = RwSignal::new(String::new());
    let start = RwSignal::new("1".to_owned());
    let end = RwSignal::new(String::new());
    let reverse = RwSignal::new(false);
    let all_seed = RwSignal::new(String::new());
    let all_substrate = RwSignal::new(String::new());
    let status = RwSignal::new(String::new());
    Effect::new(move |_| {
        let current = Arc::new(draft.get());
        let e = epoch.get();
        let n = current.alternatives.as_ref().map_or(0, |l| l.rows.len());
        position.set(current.slots.len().to_string());
        base.set(Some((current, e)));
        label.set(String::new());
        seeds.set(vec![]);
        assigned.set(vec![None; n]);
        substrates.set(vec![None; n]);
        row_page.set(0);
        choice_page.set(0);
        source_index.set(String::new());
        start.set("1".into());
        end.set(String::new());
        reverse.set(false);
        all_seed.set(String::new());
        all_substrate.set(String::new());
        status.set(String::new());
    });
    let guard = Callback::new(move |()| {
        if busy.get_untracked() {
            status.set("Wait for the workspace operation to finish.".into());
            return false;
        }
        let valid = base.with_untracked(|saved| {
            saved.as_ref().is_some_and(|(saved, e)| {
                *e == epoch.get_untracked()
                    && draft.with_untracked(|current| current == saved.as_ref())
            })
        });
        if !valid {
            status.set("The source library changed; configure this addition again.".into());
        }
        valid
    });
    let add_seed = Callback::new(move |seed: ChoiceSeed| {
        if !guard.run(()) {
            return;
        }
        let existing = base.with_untracked(|b| {
            b.as_ref()
                .and_then(|(d, _)| d.alternatives.as_ref())
                .map_or(0, |l| l.bins.iter().map(|b| b.choices.len()).sum::<usize>())
        });
        if seeds.with_untracked(Vec::len) >= 64
            || existing
                .saturating_add(seeds.with_untracked(Vec::len))
                .saturating_add(1)
                > 256
        {
            status.set("At most64 new choices and256 total library choices are allowed.".into());
            return;
        }
        let bounded = seeds.with_untracked(|current| {
            serde_json::to_writer(&mut SeedBudget(0), &(current, &seed)).is_ok()
        });
        if !bounded {
            status.set("New source choices exceed their32MiB editable input budget.".into());
            return;
        }
        seeds.update(|s| s.push(seed));
        status.set(
            "Choice added. Assign a source and preparation explicitly for every construct below."
                .into(),
        );
    });
    let capture = move |_| {
        if !guard.run(()) {
            return;
        }
        let result = (|| -> Result<ChoiceSeed, String> {
            let selected_event = source_index.get_untracked();
            if selected_event.is_empty() {
                return Err("Choose an open DNA source.".into());
            }
            let source = sources
                .with_untracked(|s| {
                    s.iter()
                        .find(|source| source.binding.event_id == selected_event)
                        .cloned()
                })
                .ok_or("The selected exact source event is no longer open.")?;
            if source.alphabet != Alphabet::Dna {
                return Err("Choose a DNA source.".into());
            }
            let molecule = get_source.run(source.binding.clone())?;
            if molecule.revision() != source.binding.revision
                || molecule.name != source.name
                || molecule.sequence.len() != source.length
                || molecule.topology != source.topology
                || molecule.sequence.alphabet() != source.alphabet
            {
                return Err("The exact source event changed; choose it again.".into());
            }
            let first = start
                .get_untracked()
                .trim()
                .parse::<usize>()
                .ok()
                .and_then(|v| v.checked_sub(1))
                .ok_or("Start is a1-based position.")?;
            let last = end
                .get_untracked()
                .trim()
                .parse::<usize>()
                .map_err(|_| "Enter an explicit inclusive end position.")?;
            let flipped = reverse.get_untracked();
            let slot = Slot {
                material_origin: None,
                selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                    start: first,
                    end: last,
                    reverse: flipped,
                },
                label: "New source choice".into(),
                captured: Some(Arc::new(CapturedSource { source, molecule })),
                digest_geometry: None,
                anneal_bottom: None,
            };
            let preparation = OverlapInput {
                interval_start: first.saturating_add(1).to_string(),
                interval_end: last.to_string(),
                reverse: flipped,
                ..Default::default()
            };
            crate::assembly_draft_model::validate_slot(&slot, &preparation, 0)?;
            Ok(ChoiceSeed { slot, preparation })
        })();
        match result {
            Ok(seed) => add_seed.run(seed),
            Err(error) => status.set(error),
        }
    };
    let submit = move |_| {
        if !guard.run(()) {
            return;
        }
        let result = (|| -> Result<(), String> {
            let (saved, e) = base
                .get_untracked()
                .ok_or("Library context is unavailable.")?;
            let choices = seeds.get_untracked();
            if choices.is_empty() {
                return Err("Add at least one explicit source choice.".into());
            }
            let selected = assigned.get_untracked();
            let pool = substrates.get_untracked();
            let rows = topology::addition_rows(&saved, &choices, &selected, &pool)?;
            let change = topology::Change::Add {
                position: position
                    .get_untracked()
                    .parse::<usize>()
                    .map_err(|_| "Choose the insertion position.")?,
                label: label.get_untracked(),
                choices,
                rows,
            };
            topology::propose(&saved, change.clone())?;
            on_propose.run((saved, e, change));
            Ok(())
        })();
        if let Err(error) = result {
            status.set(error);
        }
    };
    view! {
        <details class="assembly-library-add"><summary>"Add source bin"</summary>
            <p class="helper">"Choose source alternatives, then assign one source and its copied preparation to every construct. Review all effects before applying."</p>
            <label>"New bin label"<input maxlength="80" prop:value=move||label.get() on:input=move|ev|label.set(event_target_value(&ev))/></label>
            <label>"Insertion position"<select prop:value=move||position.get() on:change=move|ev|position.set(event_target_value(&ev))>
                {move||base.with(|b|b.as_ref().map(|(d,_)|(0..=d.slots.len()).map(|i|view!{<option value=i.to_string()>{format!("{} · {}",i+1,if i==d.slots.len(){"after last bin"}else{"before this bin"})}</option>}).collect_view()))}
            </select></label>
            <details><summary>"Copy an existing captured source and its preparation"</summary>
                {move|| {
                    let saved=base.get();let Some((saved,_))=saved else{return ().into_any();};
                    let Some(library)=saved.alternatives.as_ref() else{return ().into_any();};
                    let first=choice_page.get().saturating_mul(PAGE);
                    library.bins.iter().flat_map(|bin|bin.choices.iter()).filter(|choice|choice.slot.captured.is_some()).skip(first).take(PAGE).map(|choice| {
                        let seed=ChoiceSeed{slot:choice.slot.clone(),preparation:choice.row.clone()};
                        let text=seed_name(&seed);
                        view!{<p>{text}<button class="button secondary" disabled=move||busy.get() on:click=move |_|add_seed.run(seed.clone())>"Add this choice unchanged"</button></p>}
                    }).collect_view().into_any()
                }}
                <button class="button secondary" disabled=move||choice_page.get()==0 on:click=move |_|choice_page.update(|p|*p=p.saturating_sub(1))>"Previous alternatives"</button>
                <button class="button secondary" disabled=move||base.with(|b|b.as_ref().and_then(|(d,_)|d.alternatives.as_ref()).is_none_or(|l|(choice_page.get()+1).saturating_mul(PAGE)>=l.bins.iter().flat_map(|b|b.choices.iter()).filter(|c|c.slot.captured.is_some()).count())) on:click=move |_|choice_page.update(|p|*p+=1)>"Next alternatives"</button>
            </details>
            <details><summary>"Capture an open DNA source"</summary>
                <label>"Exact source event"<select prop:value=move||source_index.get() on:change=move|ev|source_index.set(event_target_value(&ev))><option value="">"Choose a DNA source…"</option>
                    {move||sources.get().into_iter().filter(|s|s.alphabet==Alphabet::Dna).map(|s|view!{<option value=s.binding.event_id.clone()>{format!("{} · {} bp",s.name,s.length)}</option>}).collect_view()}
                </select></label>
                <label>"Original contribution start (1-based)"<input type="number" min="1" prop:value=move||start.get() on:input=move|ev|start.set(event_target_value(&ev))/></label>
                <label>"Original contribution end (inclusive; at most one turn)"<input type="number" min="1" prop:value=move||end.get() on:input=move|ev|end.set(event_target_value(&ev))/></label>
                <label><input type="checkbox" prop:checked=move||reverse.get() on:change=move|ev|reverse.set(event_target_checked(&ev))/>"Reverse complement contribution"</label>
                <p class="helper">"New captures start with whole linear source preparation; circular DNA and PCR or digest preparation must be configured later. Contribution coordinates do not crop or amplify the physical source."</p>
                <button class="button secondary" disabled=move||busy.get() on:click=capture>"Capture as new choice"</button>
            </details>
            <p>{move||format!("{} new choices",seeds.with(Vec::len))}</p>
            {move||seeds.get().into_iter().enumerate().map(|(i,seed)|view!{<p>{format!("{}: {}",i+1,seed_name(&seed))}<button class="text-tool" disabled=move||busy.get() on:click=move |_|{if guard.run(()){seeds.update(|s|{s.remove(i);});assigned.update(|a|a.fill(None));all_seed.set(String::new());status.set("Choice removed; explicitly assign each construct again.".into());}}>"Remove choice"</button></p>}).collect_view()}
            <label>"Explicit source and copied preparation for all constructs"<select prop:value=move||all_seed.get() on:change=move|ev|all_seed.set(event_target_value(&ev))><option value="">"Choose…"</option>{move||seeds.get().into_iter().enumerate().map(|(i,seed)|view!{<option value=i.to_string()>{format!("{}: {}",i+1,seed_name(&seed))}</option>}).collect_view()}</select></label>
            <button class="button secondary" disabled=move||busy.get() on:click=move |_|{if guard.run(()){if let Ok(i)=all_seed.get_untracked().parse::<usize>(){if i<seeds.with_untracked(Vec::len){assigned.update(|a|a.fill(Some(i)));}}}}>"Assign this source and preparation to all constructs"</button>
            <label>"Reaction substrate for all pool constructs"<select prop:value=move||all_substrate.get() on:change=move|ev|all_substrate.set(event_target_value(&ev))><option value="">"Choose…"</option><option value="0">"Explicit prepared material"</option><option value="1">"Whole original source"</option><option value="2">"Whole PCR amplicon"</option></select></label>
            <button class="button secondary" disabled=move||busy.get() on:click=move |_|{if guard.run(()){if let Ok(i)=all_substrate.get_untracked().parse::<u8>(){if i<=2{substrates.update(|s|s.fill(Some(i)));}}}}>"Assign substrate to all pool constructs"</button>
            {move||{
                let Some((saved,_))=base.get() else{return ().into_any();};let Some(library)=saved.alternatives.as_ref() else{return ().into_any();};
                library.rows.iter().enumerate().skip(row_page.get().saturating_mul(PAGE)).take(PAGE).map(|(i,row)|{
                    let effective=topology::effective_edits(&saved,row);
                    let pool=effective.methods.golden_gate_pool.is_some();
                    let name=effective.name.to_string();
                    view!{<div class="assembly-library-row"><strong>{name}</strong><label>"Source and copied preparation"<select prop:value=move||assigned.with(|a|a.get(i).copied().flatten().map(|v|v.to_string()).unwrap_or_default()) on:change=move|ev|{let v=event_target_value(&ev).parse::<usize>().ok();if guard.run(()){assigned.update(|a|{if let Some(s)=a.get_mut(i){*s=v;}});}}><option value="">"Explicit choice required…"</option>{move||seeds.get().into_iter().enumerate().map(|(j,seed)|view!{<option value=j.to_string()>{format!("{}: {}",j+1,seed_name(&seed))}</option>}).collect_view()}</select></label>
                        <Show when=move||pool><label>"Reaction substrate"<select prop:value=move||substrates.with(|s|s.get(i).copied().flatten().map(|v|v.to_string()).unwrap_or_default()) on:change=move|ev|{let v=event_target_value(&ev).parse::<u8>().ok();if guard.run(()){substrates.update(|s|{if let Some(t)=s.get_mut(i){*t=v;}});}}><option value="">"Explicit substrate required…"</option><option value="0">"Explicit prepared material"</option><option value="1">"Whole original source"</option><option value="2">"Whole PCR amplicon"</option></select></label></Show>
                    </div>}
                }).collect_view().into_any()
            }}
            <button class="button secondary" disabled=move||row_page.get()==0 on:click=move |_|row_page.update(|p|*p=p.saturating_sub(1))>"Previous constructs"</button>
            <button class="button secondary" disabled=move||base.with(|b|b.as_ref().and_then(|(d,_)|d.alternatives.as_ref()).is_none_or(|l|(row_page.get()+1).saturating_mul(PAGE)>=l.rows.len())) on:click=move |_|row_page.update(|p|*p+=1)>"Next constructs"</button>
            <button class="button" disabled=move||busy.get() on:click=submit>"Review adding this bin"</button>
            <p role="status">{move||status.get()}</p>
        </details>
    }
}
