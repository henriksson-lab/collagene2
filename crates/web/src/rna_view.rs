//! Typed RNA declarations with cancellable worker analysis and exact-event saves.
use crate::{
    rna_model::{self, Capture},
    source_identity::SourceInfo,
    worker::CachedWorker,
};
use collagene_core::{
    Alphabet, Molecule, Topology,
    protocol::{AnalysisJob, AnalysisResult, Request},
    rna_chemistry as rna, rna_history,
    run_lineage::EventBinding,
};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
fn bases() -> Vec<(&'static str, rna::Nucleoside)> {
    use rna::Nucleoside::*;
    vec![
        ("Canonical", Canonical),
        ("N6-methyladenosine", N6Methyladenosine),
        ("5-methylcytidine", FiveMethylcytidine),
        ("5-methyluridine", FiveMethyluridine),
        ("Pseudouridine", Pseudouridine),
        ("Inosine", Inosine),
        (
            "Unknown",
            Unknown {
                label: "Unspecified nucleoside".into(),
            },
        ),
    ]
}
fn sugars() -> Vec<(&'static str, rna::Sugar)> {
    use rna::Sugar::*;
    vec![
        ("Ribose", Ribose),
        ("2′-O-methyl", TwoPrimeOMethyl),
        ("2′-deoxy", TwoPrimeDeoxy),
        ("2′-fluoro", TwoPrimeFluoro),
        (
            "Unknown",
            Unknown {
                label: "Unspecified sugar".into(),
            },
        ),
    ]
}
fn links() -> Vec<(&'static str, rna::Linkage)> {
    use rna::{Linkage::*, Stereo};
    vec![
        ("Phosphodiester", Phosphodiester),
        (
            "Phosphorothioate Rp",
            Phosphorothioate { stereo: Stereo::Rp },
        ),
        (
            "Phosphorothioate Sp",
            Phosphorothioate { stereo: Stereo::Sp },
        ),
        (
            "Phosphorothioate, stereo unresolved",
            Phosphorothioate {
                stereo: Stereo::Unspecified,
            },
        ),
        (
            "Unknown",
            Unknown {
                label: "Unspecified linkage".into(),
            },
        ),
    ]
}
fn five() -> Vec<(&'static str, rna::FivePrime)> {
    use rna::FivePrime::*;
    vec![
        (
            "Unknown",
            Unknown {
                label: "Unspecified 5′ end".into(),
            },
        ),
        ("Hydroxyl", Hydroxyl),
        ("Monophosphate", Monophosphate),
        ("Triphosphate", Triphosphate),
    ]
}
fn three() -> Vec<(&'static str, rna::ThreePrime)> {
    use rna::ThreePrime::*;
    vec![
        (
            "Unknown",
            Unknown {
                label: "Unspecified 3′ end".into(),
            },
        ),
        ("Hydroxyl", Hydroxyl),
        ("Monophosphate", Monophosphate),
        ("2′,3′ cyclic phosphate", TwoThreeCyclicPhosphate),
    ]
}
fn choice<T: Clone>(rows: Vec<(&str, T)>, value: String) -> Option<T> {
    rows.get(value.parse::<usize>().ok()?)
        .map(|(_, v)| v.clone())
}
fn index<T: PartialEq>(rows: Vec<(&str, T)>, value: &T) -> String {
    rows.iter()
        .position(|(_, v)| v == value)
        .or_else(|| {
            rows.iter()
                .position(|(_, v)| std::mem::discriminant(v) == std::mem::discriminant(value))
        })
        .expect("RNA form covers the current native enum")
        .to_string()
}
#[component]
pub fn RnaPanel(
    active_source: Memo<SourceInfo>,
    epoch: RwSignal<u64>,
    workspace_busy: RwSignal<bool>,
    resolve: Callback<EventBinding, Result<Molecule, String>>,
    saved: Memo<Result<Vec<rna_history::Summary>, String>>,
    on_save: Callback<(u64, rna_history::Record), Result<(), String>>,
    on_load: Callback<String, Result<rna_history::Record, String>>,
    on_delete: Callback<String, Result<(), String>>,
) -> impl IntoView {
    let pin = RwSignal::new(active_source.get_untracked());
    let draft = RwSignal::new(rna_model::blank(
        pin.get_untracked().binding.revision,
        pin.get_untracked().topology,
    ));
    let candidates = RwSignal::new(Vec::<rna::AnnotationCandidate>::new());
    let site_page = RwSignal::new(0usize);
    let link_page = RwSignal::new(0usize);
    let annotation_page = RwSignal::new(0usize);
    let inspection_done = RwSignal::new(false);
    let name = RwSignal::new(String::new());
    let status = RwSignal::new(String::new());
    let report = RwSignal::new(None::<(Capture, rna::Report)>);
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        draft.with(|_| ());
        let _ = (pin.get(), epoch.get(), workspace_busy.get());
        generation.update(|n| *n = n.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        report.set(None);
    });
    Effect::new(move |_| {
        let _ = (pin.get(), epoch.get());
        candidates.set(vec![]);
        inspection_done.set(false);
        site_page.set(0);
        link_page.set(0);
        annotation_page.set(0);
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
    let use_active = move |_| {
        let s = active_source.get_untracked();
        if s.alphabet != Alphabet::Rna {
            status.set("Select an RNA document first.".into());
            return;
        }
        draft.set(rna_model::blank(s.binding.revision.clone(), s.topology));
        pin.set(s);
        name.set(String::new());
        status
            .set("Pinned the active RNA history event. Enter an explicit model statement.".into());
    };
    let inspect = move |_| match resolve
        .run(pin.get_untracked().binding)
        .and_then(|m| rna::annotation_candidates(&m).map_err(|e| e.to_string()))
    {
        Ok(c) => {
            candidates.set(c);
            inspection_done.set(true);
            status.set(
                "Review each imported annotation explicitly; unsupported records stay unresolved."
                    .into(),
            );
        }
        Err(e) => status.set(e),
    };
    let run = {
        let alive = alive.clone();
        move |_| {
            if busy.get_untracked() || workspace_busy.get_untracked() {
                return;
            }
            let result = (|| {
                let s = pin.get_untracked();
                let d = draft.get_untracked();
                let m = resolve.run(s.binding.clone())?;
                let c = Capture::prepare(epoch.get_untracked(), s.binding, &m, d)?;
                Ok::<_, String>((m, c))
            })();
            let (m, c) = match result {
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
            busy.set(true);
            report.set(None);
            status.set("Calculating declared RNA composition…".into());
            let alive = alive.clone();
            wasm_bindgen_futures::spawn_local(async move {
                let reply = w
                    .ask(Request {
                        id: token,
                        revision: c.context_revision.clone(),
                        job: AnalysisJob::AnalyzeRnaChemistry {
                            molecule: m,
                            declaration: c.declaration.clone(),
                        },
                    })
                    .await;
                w.close();
                if !alive.load(Ordering::Relaxed) || generation.get_untracked() != token {
                    return;
                }
                busy.set(false);
                worker.set_value(None);
                if let Err(e) = c.guard(
                    epoch.get_untracked(),
                    &pin.get_untracked().binding,
                    &draft.get_untracked(),
                ) {
                    status.set(e);
                    return;
                }
                match reply.and_then(|r| match r.result {
                    AnalysisResult::RnaChemistry(value) => {
                        c.validate_report(&value)?;
                        Ok(value)
                    }
                    AnalysisResult::Error(e) => Err(e),
                    _ => Err("Unexpected RNA analysis response.".into()),
                }) {
                    Ok(r) => {
                        status.set(if r.formula_complete {
                            "Composition available under the explicit model.".into()
                        } else {
                            "Quantities unavailable; review unresolved chemistry below.".into()
                        });
                        report.set(Some((c, r)));
                    }
                    Err(e) => status.set(e),
                }
            });
        }
    };
    let site_position = RwSignal::new("1".to_owned());
    let base = RwSignal::new(rna::Nucleoside::Canonical);
    let sugar = RwSignal::new(rna::Sugar::Ribose);
    let site_note = RwSignal::new(String::new());
    let site_ref = RwSignal::new(String::new());
    let add_site = move |_| {
        let result = (|| {
            let i =
                rna_model::position(&site_position.get_untracked(), pin.get_untracked().length)?;
            let note = site_note.get_untracked();
            if note.trim().is_empty() {
                return Err("Describe the residue declaration's provenance.".into());
            }
            let reference = site_ref.get_untracked();
            let s = rna::Site {
                index: i,
                nucleoside: base.get_untracked(),
                sugar: sugar.get_untracked(),
                provenance: rna::Provenance {
                    statement: note,
                    reference: (!reference.trim().is_empty()).then_some(reference),
                },
            };
            draft
                .try_update(|d| rna_model::upsert_site(d, s))
                .ok_or("RNA panel is closed.")??;
            Ok::<_, String>(())
        })();
        match result {
            Ok(()) => {
                status.set("Residue entry updated; calculate to check the complete model.".into())
            }
            Err(e) => status.set(e),
        }
    };
    let bond_position = RwSignal::new("1".to_owned());
    let bond = RwSignal::new(rna::Linkage::Phosphodiester);
    let bond_note = RwSignal::new(String::new());
    let add_link = move |_| {
        let result = (|| {
            let p = pin.get_untracked();
            let count = if p.topology == Topology::Circular {
                p.length
            } else {
                p.length.saturating_sub(1)
            };
            let after = rna_model::position(&bond_position.get_untracked(), count)?;
            let note = bond_note.get_untracked();
            if note.trim().is_empty() {
                return Err("Describe the linkage declaration's provenance.".into());
            }
            let l = rna::Link {
                after,
                linkage: bond.get_untracked(),
                provenance: rna::Provenance {
                    statement: note,
                    reference: None,
                },
            };
            draft.update(|d| {
                if let Some(old) = d.links.iter_mut().find(|l| l.after == after) {
                    *old = l;
                } else {
                    d.links.push(l);
                }
            });
            Ok::<_, String>(())
        })();
        match result {
            Ok(()) => status.set("Linkage entry updated.".into()),
            Err(e) => status.set(e),
        }
    };
    view! {<section class="inspector-content rna-properties" aria-label="RNA chemistry" style="overflow-wrap:anywhere">
        <h2>"RNA chemistry"</h2><p class="helper">"Declare a covalent model for an exact RNA version. Calculations describe a neutral acid without counterions or adducts; they do not establish experimental chemistry or enzyme behavior."</p>
        <button class="button secondary full" disabled=move||workspace_busy.get()||active_source.with(|s|s.alphabet!=Alphabet::Rna) on:click=use_active>"Use active RNA version"</button>
        <p class="helper">{move||pin.with(|s|format!("Pinned: {} · {} residues · event {}",s.name,s.length,s.binding.event_id.chars().take(18).collect::<String>()))}</p>
        <Show when=move||pin.with(|s|s.alphabet==Alphabet::Rna) fallback=move||view!{<p class="helper">"Select an RNA document, then use its active version."</p>}>
        <fieldset disabled=move||workspace_busy.get() style="border:0;padding:0;min-width:0">
        <label>"Model statement"<textarea aria-label="RNA model statement" maxlength="1024" prop:value=move||draft.with(|d|d.provenance.statement.clone()) on:input=move|e|draft.update(|d|d.provenance.statement=event_target_value(&e)) /></label>
        <label>"Model reference (optional)"<input aria-label="RNA model reference" maxlength="4096" prop:value=move||draft.with(|d|d.provenance.reference.clone().unwrap_or_default()) on:input=move|e|{let value=event_target_value(&e);draft.update(|d|d.provenance.reference=(!value.trim().is_empty()).then_some(value));}/></label>
        <label>"Coverage"<select aria-label="RNA coverage" prop:value=move||if draft.with(|d|d.coverage==rna::Coverage::ClosedModel){"closed"}else{"partial"} on:change=move|e|draft.update(|d|d.coverage=if event_target_value(&e)=="closed"{rna::Coverage::ClosedModel}else{rna::Coverage::Partial})><option value="partial">"Partial — unspecified chemistry remains unresolved"</option><option value="closed">"Closed model — all unlisted residues/bonds canonical"</option></select></label>
        <Show when=move||pin.with(|s|s.topology==Topology::Linear) fallback=move||view!{<p class="helper">"Circular model: includes the closing 3′–5′ bond, with no linear end groups."</p>}>
            <label>"5′ terminus"<select aria-label="RNA five prime terminus" prop:value=move||draft.with(|d|if let rna::Ends::Linear{five:v,..}=&d.ends{index(five(),v)}else{"0".into()}) on:change=move|e|{if let Some(v)=choice(five(),event_target_value(&e)){draft.update(|d|if let rna::Ends::Linear{five,..}=&mut d.ends{*five=v;});}}>{five().into_iter().enumerate().map(|(i,(name,_))|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label>
            <Show when=move||draft.with(|d|matches!(d.ends,rna::Ends::Linear{five:rna::FivePrime::Unknown{..},..}))><label>"Unknown 5′ end description"<input aria-label="Unknown RNA five prime end" maxlength="256" prop:value=move||draft.with(|d|if let rna::Ends::Linear{five:rna::FivePrime::Unknown{label},..}=&d.ends{label.clone()}else{String::new()}) on:input=move|e|draft.update(|d|if let rna::Ends::Linear{five,..}=&mut d.ends{*five=rna::FivePrime::Unknown{label:event_target_value(&e)};})/></label></Show>
            <label>"3′ terminus"<select aria-label="RNA three prime terminus" prop:value=move||draft.with(|d|if let rna::Ends::Linear{three:v,..}=&d.ends{index(three(),v)}else{"0".into()}) on:change=move|e|{if let Some(v)=choice(three(),event_target_value(&e)){draft.update(|d|if let rna::Ends::Linear{three,..}=&mut d.ends{*three=v;});}}>{three().into_iter().enumerate().map(|(i,(name,_))|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label>
            <Show when=move||draft.with(|d|matches!(d.ends,rna::Ends::Linear{three:rna::ThreePrime::Unknown{..},..}))><label>"Unknown 3′ end description"<input aria-label="Unknown RNA three prime end" maxlength="256" prop:value=move||draft.with(|d|if let rna::Ends::Linear{three:rna::ThreePrime::Unknown{label},..}=&d.ends{label.clone()}else{String::new()}) on:input=move|e|draft.update(|d|if let rna::Ends::Linear{three,..}=&mut d.ends{*three=rna::ThreePrime::Unknown{label:event_target_value(&e)};})/></label></Show>
        </Show>
        <details><summary>"Residue declarations"</summary><p class="helper">"Positions are one-based on the stored strand. Base and sugar chemistry share one entry per residue."</p>
            <label>"Residue position"<input aria-label="RNA residue position" inputmode="numeric" prop:value=move||site_position.get() on:input=move|e|site_position.set(event_target_value(&e))/></label>
            <label>"Nucleoside"<select aria-label="RNA nucleoside" prop:value=move||index(bases(),&base.get()) on:change=move|e|if let Some(v)=choice(bases(),event_target_value(&e)){base.set(v);}>{bases().into_iter().enumerate().map(|(i,(name,_))|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label>
            <Show when=move||matches!(base.get(),rna::Nucleoside::Unknown{..})><label>"Unknown nucleoside description"<input aria-label="Unknown RNA nucleoside" maxlength="256" prop:value=move||if let rna::Nucleoside::Unknown{label}=base.get(){label}else{String::new()} on:input=move|e|base.set(rna::Nucleoside::Unknown{label:event_target_value(&e)})/></label></Show>
            <label>"Sugar"<select aria-label="RNA sugar" prop:value=move||index(sugars(),&sugar.get()) on:change=move|e|if let Some(v)=choice(sugars(),event_target_value(&e)){sugar.set(v);}>{sugars().into_iter().enumerate().map(|(i,(name,_))|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label>
            <Show when=move||matches!(sugar.get(),rna::Sugar::Unknown{..})><label>"Unknown sugar description"<input aria-label="Unknown RNA sugar" maxlength="256" prop:value=move||if let rna::Sugar::Unknown{label}=sugar.get(){label}else{String::new()} on:input=move|e|sugar.set(rna::Sugar::Unknown{label:event_target_value(&e)})/></label></Show>
            <label>"Residue provenance"<input aria-label="RNA residue provenance" maxlength="1024" prop:value=move||site_note.get() on:input=move|e|site_note.set(event_target_value(&e))/></label>
            <label>"Reference (optional)"<input aria-label="RNA residue reference" maxlength="4096" prop:value=move||site_ref.get() on:input=move|e|site_ref.set(event_target_value(&e))/></label>
            <button class="button secondary full" on:click=add_site>"Add or replace residue"</button>
            {move||draft.with(|d|d.sites.iter().enumerate().skip(site_page.get().saturating_mul(50)).take(50).map(|(i,s)|(i,s.clone())).collect::<Vec<_>>()).into_iter().map(|(i,s)|{let edit=s.clone();view!{<div class="rna-residue"><p class="helper">{format!("{} · {:?} · {:?}",s.index+1,s.nucleoside,s.sugar)}</p><button class="text-tool" on:click=move|_|{site_position.set((edit.index+1).to_string());base.set(edit.nucleoside.clone());sugar.set(edit.sugar.clone());site_note.set(edit.provenance.statement.clone());site_ref.set(edit.provenance.reference.clone().unwrap_or_default());}>"Edit residue"</button><button class="text-tool" on:click=move|_|{draft.update(|d|rna_model::remove_site(d,i));site_page.set(0);}>"Remove residue and its reviews"</button></div>}}).collect_view()}
            <p class="helper">{move||format!("{} residue entries · page {}",draft.with(|d|d.sites.len()),site_page.get()+1)}</p><button class="text-tool" disabled=move||site_page.get()==0 on:click=move|_|site_page.update(|p|*p=p.saturating_sub(1))>"Previous residues"</button><button class="text-tool" disabled=move||site_page.get().saturating_add(1).saturating_mul(50)>=draft.with(|d|d.sites.len()) on:click=move|_|site_page.update(|p|*p+=1)>"Next residues"</button>
        </details>
        <details><summary>"Backbone declarations"</summary><p class="helper">"Bond after a one-based residue; the final position is the circular seam only. Rp/Sp share a formula but remain distinct identities."</p>
            <label>"Bond after residue"<input aria-label="RNA linkage position" inputmode="numeric" prop:value=move||bond_position.get() on:input=move|e|bond_position.set(event_target_value(&e))/></label>
            <label>"Linkage"<select aria-label="RNA linkage" prop:value=move||index(links(),&bond.get()) on:change=move|e|if let Some(v)=choice(links(),event_target_value(&e)){bond.set(v);}>{links().into_iter().enumerate().map(|(i,(name,_))|view!{<option value=i.to_string()>{name}</option>}).collect_view()}</select></label>
            <Show when=move||matches!(bond.get(),rna::Linkage::Unknown{..})><label>"Unknown linkage description"<input aria-label="Unknown RNA linkage" maxlength="256" prop:value=move||if let rna::Linkage::Unknown{label}=bond.get(){label}else{String::new()} on:input=move|e|bond.set(rna::Linkage::Unknown{label:event_target_value(&e)})/></label></Show>
            <label>"Linkage provenance"<input aria-label="RNA linkage provenance" maxlength="1024" prop:value=move||bond_note.get() on:input=move|e|bond_note.set(event_target_value(&e))/></label><button class="button secondary full" disabled=move||pin.with(|s|s.length<2&&s.topology==Topology::Linear) on:click=add_link>"Add or replace linkage"</button>
            {move||draft.with(|d|d.links.iter().enumerate().skip(link_page.get().saturating_mul(50)).take(50).map(|(i,s)|(i,s.clone())).collect::<Vec<_>>()).into_iter().map(|(i,l)|{let edit=l.clone();view!{<div><p class="helper">{format!("After {} · {:?}",l.after+1,l.linkage)}</p><button class="text-tool" on:click=move|_|{bond_position.set((edit.after+1).to_string());bond.set(edit.linkage.clone());bond_note.set(edit.provenance.statement.clone());}>"Edit linkage"</button><button class="text-tool" on:click=move|_|{draft.update(|d|{if i<d.links.len(){d.links.remove(i);}});link_page.set(0);}>"Remove linkage"</button></div>}}).collect_view()}
            <p class="helper">{move||format!("{} linkage entries · page {}",draft.with(|d|d.links.len()),link_page.get()+1)}</p><button class="text-tool" disabled=move||link_page.get()==0 on:click=move|_|link_page.update(|p|*p=p.saturating_sub(1))>"Previous linkages"</button><button class="text-tool" disabled=move||link_page.get().saturating_add(1).saturating_mul(50)>=draft.with(|d|d.links.len()) on:click=move|_|link_page.update(|p|*p+=1)>"Next linkages"</button>
        </details>
        <details><summary>"Imported modification annotations"</summary><p class="helper">"Only the components named by an imported code constrain the combined residue entry. Accepting a review never supplies closed coverage."</p><button class="button secondary full" on:click=inspect>"Inspect imported annotations"</button><Show when=move||inspection_done.get()&&candidates.with(Vec::is_empty)><p class="helper">"No supported modification feature records were found. Opaque metadata, if present, still blocks quantities."</p></Show>
            {move||candidates.with(|rows|rows.iter().skip(annotation_page.get().saturating_mul(50)).take(50).cloned().collect::<Vec<_>>()).into_iter().map(|c|{let id=c.feature_index;let label=format!("Feature {} · {} · {}",id+1,c.raw_code.as_ref().map(|s|s.chars().take(80).collect::<String>()).unwrap_or_else(||"Missing or repeated code".into()),c.position.map(|p|format!("residue {}",p+1)).unwrap_or_else(||"uncertain location".into()));let supported=c.suggested_site.is_some();view!{<div><p class="helper">{label}</p><p class="helper">{c.suggested_site.map(|s|format!("Suggested residue: {:?} / {:?}",s.nucleoside,s.sugar)).or(c.unsupported_reason)}</p><button class="text-tool" disabled=move||!supported||draft.with(|d|d.annotation_reviews.iter().any(|r|r.feature_index==id)) on:click=move|_|{let result=resolve.run(pin.get_untracked().binding).and_then(|m|{let mut d=draft.get_untracked();rna_model::accept_annotation(&mut d,&m,id)?;draft.set(d);Ok(())});status.set(result.map(|_|"Feature reviewed against its exact source and residue entry.".into()).unwrap_or_else(|e|e));}>"Review and apply constraint"</button><button class="text-tool" disabled=move||!draft.with(|d|d.annotation_reviews.iter().any(|r|r.feature_index==id)) on:click=move|_|draft.update(|d|d.annotation_reviews.retain(|r|r.feature_index!=id))>"Withdraw review"</button></div>}}).collect_view()}
            <p class="helper">{move||format!("{} annotation records · page {}",candidates.with(Vec::len),annotation_page.get()+1)}</p><button class="text-tool" disabled=move||annotation_page.get()==0 on:click=move|_|annotation_page.update(|p|*p=p.saturating_sub(1))>"Previous annotations"</button><button class="text-tool" disabled=move||annotation_page.get().saturating_add(1).saturating_mul(50)>=candidates.with(Vec::len) on:click=move|_|annotation_page.update(|p|*p+=1)>"Next annotations"</button>
        </details>
        <button class="button primary full" disabled=move||busy.get()||draft.with(|d|d.provenance.statement.trim().is_empty()) on:click=run.clone()>"Calculate RNA composition"</button>
        <Show when=move||busy.get()><button class="text-tool" on:click=move|_|{generation.update(|n|*n=n.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);busy.set(false);status.set("RNA calculation cancelled.".into());}>"Cancel RNA calculation"</button></Show>
        </fieldset>
        {move||report.with(|item|item.as_ref().map(|(_,r)|view!{<div class="rna-result"><h3>"Declared composition"</h3><dl><dt>"Formula"</dt><dd>{r.formula.map(|f|f.hill()).unwrap_or_else(||"Unavailable".into())}</dd><dt>"Monoisotopic mass"</dt><dd>{r.mass.as_ref().map(|m|format!("{:.6} Da",m.monoisotopic_da)).unwrap_or_else(||"Unavailable".into())}</dd><dt>"Average mass"</dt><dd>{r.mass.as_ref().map(|m|format!("{:.3} Da",m.average_da)).unwrap_or_else(||"Unavailable".into())}</dd></dl><p class="helper">{format!("Formula complete: {} · identity complete: {} · unresolved PS stereocenters: {} · unreviewed source records: {}",r.formula_complete,r.identity_complete,r.unspecified_stereocenters.len(),r.unreviewed_evidence.len())}</p>{r.issues.clone().into_iter().map(|i|view!{<p class="helper">{format!("{:?}: {} · residue examples {}{}",i.kind,i.count,i.positions.iter().map(|p|(p+1).to_string()).collect::<Vec<_>>().join(", "),if i.positions_truncated{" …"}else{""})}</p>}).collect_view()}<button class="button secondary full" on:click=move|_|{if let Some((c,r))=report.get_untracked(){if let Err(e)=c.guard(epoch.get_untracked(),&pin.get_untracked().binding,&draft.get_untracked()){status.set(e);return;}match serde_json::to_vec_pretty(&r){Ok(bytes)=>crate::app::download_with_completion("rna-chemistry.json","application/json",&bytes,move|result|{if status.try_get_untracked().is_some(){if let Err(e)=result{status.set(e);}}}),Err(e)=>status.set(e.to_string())}}}>"Download RNA report"</button></div>}))}
        <label>"Saved declaration name"<input aria-label="RNA saved declaration name" maxlength="256" prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label><button class="button secondary full" disabled=move||workspace_busy.get()||busy.get()||report.with(Option::is_none)||name.with(|s|s.trim().is_empty()) on:click=move|_|{if let Some((c,_))=report.get_untracked(){let result=c.guard(epoch.get_untracked(),&pin.get_untracked().binding,&draft.get_untracked()).and_then(|()|on_save.run((c.epoch,rna_history::Record{version:rna_history::VERSION,name:name.get_untracked().trim().into(),source:c.source,declaration:c.declaration})));status.set(result.map(|_|"Saved the declaration with its exact source event; quantities are recomputed on reopening.".into()).unwrap_or_else(|e|e));}}>"Save RNA declaration"</button>
        </Show>
        <details><summary>"Saved RNA declarations"</summary>{move||match saved.get(){Err(e)=>view!{<p role="alert">{e}</p>}.into_any(),Ok(rows)=>view!{<div>{rows.into_iter().map(|s|{let hash=s.hash.clone();let delete=s.hash;view!{<div><p class="helper">{format!("{} · {} · event {}",s.name,s.source_name,s.source.event_id.chars().take(18).collect::<String>())}</p><button class="text-tool" disabled=move||workspace_busy.get() on:click=move|_|{match on_load.run(hash.clone()).and_then(|r|resolve.run(r.source.clone()).map(|m|(r,m))){Ok((r,m))=>{batch(||{pin.set(SourceInfo{binding:r.source,name:m.name.clone(),length:m.sequence.len(),alphabet:m.sequence.alphabet(),topology:m.topology});draft.set(r.declaration);name.set(r.name);report.set(None);});status.set("Opened the saved declaration and exact historical source. Calculate to obtain fresh quantities.".into());},Err(e)=>status.set(e)}}>"Open declaration"</button><button class="text-tool" disabled=move||workspace_busy.get() on:click=move|_|status.set(on_delete.run(delete.clone()).map(|_|"Saved declaration deleted.".into()).unwrap_or_else(|e|e))>"Delete declaration"</button></div>}}).collect_view()}</div>}.into_any()}}</details>
        <p role="status" class="helper">{move||status.get()}</p>
    </section>}
}
