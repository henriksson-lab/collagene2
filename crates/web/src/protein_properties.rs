//! Explicit sequence-based protein estimates, isolated from biological verification authority.
use crate::worker::CachedWorker;
use collagene_core::{
    protein_properties::*,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Alphabet, Molecule};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

fn interval(value: &NumericInterval, places: usize) -> String {
    if value.minimum == value.maximum {
        format!("{:.*}", places, value.minimum)
    } else {
        format!("{:.*}–{:.*}", places, value.minimum, places, value.maximum)
    }
}

#[component]
pub fn ProteinPanel(molecule: Memo<Molecule>) -> impl IntoView {
    let stop = RwSignal::new(StopPolicy::AllowSingleTerminal);
    let ambiguity = RwSignal::new(AmbiguityPolicy::Unavailable);
    let disulfide_mode = RwSignal::new("reduced".to_string());
    let bond_count = RwSignal::new("0".to_string());
    let modifications = RwSignal::new(ModificationState::AssumeUnmodified);
    let ph = RwSignal::new("7.0".to_string());
    let limit = RwSignal::new("1000000".to_string());
    let options = Memo::new(move |_| {
        let charge_ph = ph
            .get()
            .trim()
            .parse::<f64>()
            .map_err(|_| "Enter a numeric protein charge pH.".to_string())?;
        let max_residues = limit
            .get()
            .trim()
            .parse::<usize>()
            .map_err(|_| "Enter a whole-number protein residue limit.".to_string())?;
        let disulfides =
            match disulfide_mode.get().as_str() {
                "maximum" => DisulfideState::MaximumIntrachain,
                "count" => {
                    DisulfideState::Count(bond_count.get().trim().parse::<usize>().map_err(
                        |_| "Enter a nonnegative whole intrachain bond count.".to_string(),
                    )?)
                }
                _ => DisulfideState::Reduced,
            };
        let o = ProteinOptions {
            stop_policy: stop.get(),
            ambiguity: ambiguity.get(),
            disulfides,
            modification_state: modifications.get(),
            charge_ph,
            max_residues,
        };
        o.validate().map_err(|e| e.to_string())?;
        Ok::<ProteinOptions, String>(o)
    });
    let source = Memo::new(move |_| molecule.with(Molecule::revision));
    let report = RwSignal::new(None::<ProteinReport>);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        let _ = (source.get(), options.get());
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(worker) = active.get_value() {
            worker.close();
        }
        active.set_value(None);
        report.set(None);
        busy.set(false);
        status.set(String::new());
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            if let Some(worker) = active.get_value() {
                worker.close();
            }
            active.set_value(None);
        });
    }
    let run = {
        let alive = alive.clone();
        move |_: web_sys::MouseEvent| {
            let opts = match options.get_untracked() {
                Ok(o) => o,
                Err(e) => {
                    status.set(e);
                    return;
                }
            };
            if molecule.with_untracked(|m| m.sequence.alphabet() != Alphabet::Protein) {
                status.set("Choose an explicit protein document; nucleotide translation is a separate step.".into());
                return;
            }
            if let Some(worker) = active.get_value() {
                worker.close();
            }
            let worker = match CachedWorker::new() {
                Ok(worker) => Rc::new(worker),
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            active.set_value(Some(worker.clone()));
            let revision = source.get_untracked();
            let token = generation.get_untracked().wrapping_add(1);
            generation.set(token);
            report.set(None);
            busy.set(true);
            status.set("Calculating protein properties…".into());
            let alive = alive.clone();
            let molecule = molecule.get_untracked();
            spawn_local(async move {
                let response = worker
                    .ask(Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::ProteinProperties {
                            molecule,
                            options: opts.clone(),
                        },
                    })
                    .await;
                worker.close();
                if !alive.load(Ordering::Relaxed) {
                    return;
                }
                if generation.get_untracked() != token
                    || source.get_untracked() != revision
                    || options.get_untracked() != Ok(opts.clone())
                {
                    return;
                }
                busy.set(false);
                active.set_value(None);
                match response {
                    Ok(response) => match response.result {
                        AnalysisResult::ProteinProperties(value)
                            if value.source_revision == revision && value.options == opts =>
                        {
                            status.set(if value.complete{"All requested estimates available under the declared model.".into()}else{"Some identities or quantities remain unresolved; review the reasons and assumptions.".into()});
                            report.set(Some(value));
                        }
                        AnalysisResult::Error(error) => status.set(error),
                        _ => status.set("Unexpected protein-properties response".into()),
                    },
                    Err(error) => status.set(error),
                }
            });
        }
    };
    view! {<section class="inspector-content protein-properties" style="overflow-wrap:anywhere" aria-label="Protein properties">
        <h2>"Protein properties"</h2>
        <Show when=move||molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein) fallback=move||view!{<p class="helper">"Open an explicit protein document to calculate mass, composition, theoretical pI and A280. Translate a nucleotide source separately; this panel never guesses a protein alphabet."</p>}>
            <p class="helper">"Sequence/model estimates for one linear chain. Free termini and the declared cysteine state are assumptions; folding, processing and experimental function are not inferred."</p>
            <label>"Protein stop policy"<select aria-label="Protein stop policy" prop:value=move||if stop.get()==StopPolicy::Reject{"reject"}else{"terminal"} on:change=move|e|stop.set(if event_target_value(&e)=="reject"{StopPolicy::Reject}else{StopPolicy::AllowSingleTerminal})><option value="terminal">"Allow one terminal stop"</option><option value="reject">"Reject all stops"</option></select></label>
            <label>"Protein ambiguity"<select aria-label="Protein ambiguity" prop:value=move||if ambiguity.get()==AmbiguityPolicy::Unavailable{"unavailable"}else{"bounds"} on:change=move|e|ambiguity.set(if event_target_value(&e)=="bounds"{AmbiguityPolicy::CanonicalBounds}else{AmbiguityPolicy::Unavailable})><option value="unavailable">"Leave unsupported quantities unavailable"</option><option value="bounds">"Explicit canonical-residue bounds"</option></select></label>
            <label>"Protein cysteine state"<select aria-label="Protein cysteine state" prop:value=move||disulfide_mode.get() on:change=move|e|disulfide_mode.set(event_target_value(&e))><option value="reduced">"All cysteines reduced"</option><option value="maximum">"Maximum intrachain pairs"</option><option value="count">"Declare intrachain bond count"</option></select></label>
            <Show when=move||disulfide_mode.get()=="count"><label>"Intrachain disulfide bonds"<input type="number" min="0" step="1" aria-label="Intrachain disulfide bonds" prop:value=move||bond_count.get() on:input=move|e|bond_count.set(event_target_value(&e))/></label></Show>
            <label>"Protein modification state"<select aria-label="Protein modification state" prop:value=move||if modifications.get()==ModificationState::AssumeUnmodified{"unmodified"}else{"unknown"} on:change=move|e|modifications.set(if event_target_value(&e)=="unknown"{ModificationState::KnownOrUnknownModifications}else{ModificationState::AssumeUnmodified})><option value="unmodified">"Assume no other modifications"</option><option value="unknown">"Known or unknown modifications"</option></select></label>
            <label>"Charge at pH"<input type="number" min="0" max="14" step="0.1" aria-label="Protein charge pH" prop:value=move||ph.get() on:input=move|e|ph.set(event_target_value(&e))/></label>
            <details><summary>"Protein calculation limit"</summary><label>"Maximum analyzed residues"<input type="number" min="1" max="10000000" step="1" aria-label="Protein residue limit" prop:value=move||limit.get() on:input=move|e|limit.set(event_target_value(&e))/></label><p class="helper">"A limit failure rejects the calculation; it never reports a partial molecular mass."</p></details>
            {move||options.get().err().map(|error|view!{<p role="alert" class="helper">{error}</p>})}
            <button class="button primary full" disabled=move||busy.get()||options.get().is_err() on:click=run.clone()>"Calculate protein properties"</button>
            <Show when=move||busy.get()><button class="text-tool" on:click=move|_|{generation.update(|n|*n=n.wrapping_add(1));if let Some(worker)=active.get_value(){worker.close();}active.set_value(None);busy.set(false);status.set("Protein calculation cancelled.".into());}>"Cancel protein calculation"</button></Show>
            <p class="helper protein-status" role="status">{move||status.get()}</p>
            {move||report.get().map(|r|{
                let mass=r.mass.as_ref().map(|m|format!("{} Da",interval(&m.average_da,3))).unwrap_or_else(||"Unavailable".into());
                let mono=r.mass.as_ref().map(|m|format!("{} Da",interval(&m.monoisotopic_da,6))).unwrap_or_else(||"Unavailable".into());
                let pi=r.isoelectric_point.map(|v|format!("{v:.4} (model)")).unwrap_or_else(||"Unavailable".into());
                let charge=r.net_charge_at_ph.map(|v|format!("{v:+.4} at pH {}",r.options.charge_ph)).unwrap_or_else(||"Unavailable".into());
                let epsilon=r.extinction_280.as_ref().map(|e|format!("{} L mol⁻¹ cm⁻¹",interval(&e.epsilon_molar_inverse_cm,0))).unwrap_or_else(||"Unavailable".into());
                let absorbance=r.extinction_280.as_ref().and_then(|e|e.absorbance_1mg_ml_1cm.as_ref()).map(|v|format!("{} (1 mg/mL, 1 cm)",interval(v,4))).unwrap_or_else(||"Unavailable".into());
                let complete=r.complete;let identity=r.mass.as_ref().is_some_and(|m|m.exact_identity);let source_label=r.source_revision.chars().take(16).collect::<String>();
                view!{<div class="protein-result" data-complete=complete.to_string()>
                    <h3>"Protein model report"</h3><p class="helper">{format!("{} residues · {} · source {source_label}",r.residue_count,if identity{"resolved identities"}else{"identities or chemistry unresolved"})}</p>
                    <dl style="font-size:13px;overflow-wrap:anywhere;line-height:1.6"><dt>"Average neutral mass"</dt><dd style="margin:0 0 10px" class="protein-average-mass">{mass}</dd><dt>"Monoisotopic neutral mass"</dt><dd style="margin:0 0 10px" class="protein-mono-mass">{mono}</dd><dt>"Theoretical pI"</dt><dd style="margin:0 0 10px" class="protein-pi">{pi}</dd><dt>"Model net charge"</dt><dd style="margin:0 0 10px" class="protein-charge">{charge}</dd><dt>"A280 molar extinction"</dt><dd style="margin:0 0 10px" class="protein-epsilon">{epsilon}</dd><dt>"Predicted absorbance"</dt><dd style="margin:0 0 10px" class="protein-absorbance">{absorbance}</dd><dt>"Declared intrachain disulfides"</dt><dd style="margin:0 0 10px">{interval(&r.disulfide_bonds,0)}</dd></dl>
                    <details open><summary>"Unavailable quantities and identity limits"</summary>{r.unavailable.into_iter().map(|reason|view!{<p class="helper">{reason}</p>}).collect_view()}<p class="helper">{if complete{"Every requested quantity is available within the declared model; this is not empirical validation."}else{"A point-valued number does not establish resolved identity or complete biological semantics."}}</p></details>
                    <details><summary>"Residue composition"</summary><table style="font-size:12px;width:100%"><thead><tr><th scope="col">"Symbol"</th><th scope="col">"Count"</th><th scope="col">"Fraction"</th></tr></thead><tbody>{r.composition.into_iter().map(|c|view!{<tr><th scope="row">{c.symbol}</th><td>{c.count}</td><td>{format!("{:.6}",c.fraction)}</td></tr>}).collect_view()}</tbody></table></details>
                    <details><summary>"Models and assumptions"</summary><p class="helper">{format!("{} · {} · {} · {}",r.engine_version,MASS_MODEL,r.charge_model,EXTINCTION_MODEL)}</p>{r.warnings.into_iter().map(|warning|view!{<p class="helper">{warning}</p>}).collect_view()}</details>
                    <button class="button secondary full" on:click=move|_|if let Some(value)=report.get_untracked(){match serde_json::to_vec_pretty(&value){Ok(bytes)=>crate::app::download_with_completion("protein-properties.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_none(){return;}if let Err(error)=result {status.set(error);}}),Err(error)=>status.set(error.to_string())}}>"Download protein report"</button>
                </div>}
            })}
        </Show>
    </section>}
}
