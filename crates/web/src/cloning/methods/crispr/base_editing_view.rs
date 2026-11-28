//! Conditional base-editor reference states for the selected exact guide.
use super::*;
use collagene_core::{construct_qc, sequence::crispr::base_editing as editing};

#[derive(Clone)]
struct Capture {
    source: DocumentInfo,
    epoch: u64,
    request: editing::Request,
    revision: String,
}

fn intended_positions(text: &str) -> Result<Vec<usize>, String> {
    if text.trim().is_empty() {
        return Ok(vec![]);
    }
    if text.len() > 128 {
        return Err("Enter at most 12 intended guide positions.".into());
    }
    let positions = text
        .split(',')
        .map(|part| {
            let part = part.trim();
            if part.is_empty() || !part.bytes().all(|b| b.is_ascii_digit()) {
                return Err("Separate intended guide positions with commas.".into());
            }
            whole(part, "intended guide position")
        })
        .collect::<Result<Vec<usize>, String>>()?;
    if positions.len() > 12 {
        return Err("Enter at most 12 intended guide positions.".into());
    }
    Ok(positions)
}

#[component]
pub(super) fn BaseEditingPanel(
    source: Memo<DocumentInfo>,
    options: Memo<Result<GuideOptions, String>>,
    guide: Memo<Option<GuideCandidate>>,
    get_document: Callback<String, Option<Molecule>>,
    epoch: RwSignal<u64>,
    on_navigate: Callback<ReadNavigation>,
) -> impl IntoView {
    let intended = RwSignal::new(String::new());
    let limit = RwSignal::new("32".to_string());
    let cds = RwSignal::new(String::new());
    let genetic_code = RwSignal::new("1".to_string());
    let status = RwSignal::new(String::new());
    let report = RwSignal::new(None::<(Capture, Arc<editing::Report>)>);
    let page = RwSignal::new(0usize);
    let job = JobSlot::new();
    let alive = Arc::new(AtomicBool::new(true));
    Effect::new(move |_| {
        let _ = (
            source.get(),
            epoch.get(),
            guide.get(),
            options.get(),
            intended.get(),
            limit.get(),
            cds.get(),
            genetic_code.get(),
        );
        job.stop();
        report.set(None);
        page.set(0);
        status.set("Choose a guide and declare any intended positions before enumerating reference states.".into());
    });
    {
        let alive = alive.clone();
        on_cleanup(move || {
            alive.store(false, Ordering::Relaxed);
            job.stop();
        });
    }
    let prepare = Callback::new(move |(): ()| -> Result<(Molecule, Capture), String> {
        let pin = source.get_untracked();
        let molecule = get_document
            .run(pin.binding.document_id.clone())
            .ok_or("The selected source is no longer available.")?;
        if molecule.revision() != pin.binding.revision {
            return Err("The guide source changed; design guides again.".into());
        }
        let candidate = guide
            .get_untracked()
            .ok_or("Select a current guide first.")?;
        let selected_cds = if cds.get_untracked().is_empty() {
            None
        } else {
            let index = whole::<usize>(&cds.get_untracked(), "CDS annotation")?;
            let feature = molecule
                .features
                .get(index)
                .ok_or("The selected CDS annotation changed.")?;
            if !feature.kind.eq_ignore_ascii_case("CDS") || feature.needs_review {
                return Err("Choose a CDS annotation whose coordinates have been reviewed.".into());
            }
            Some(editing::CdsRequest {
                annotation: construct_qc::AnnotationAssumption {
                    feature_index: index,
                    source_revision: pin.binding.revision.clone(),
                    annotation_revision: construct_qc::annotation_revision(feature).map_err(|e|e.to_string())?,
                    rationale: "User explicitly assumed the selected CDS annotation role and coordinates for ordinary-codon consequences; biological annotation truth is not established.".into(),
                },
                genetic_code: whole(&genetic_code.get_untracked(), "genetic code")?,
            })
        };
        let request = editing::Request {
            guide_options: options.get_untracked()?,
            guide: candidate,
            profile: editing::be1_2016().map_err(|e| e.to_string())?,
            intended_positions: intended_positions(&intended.get_untracked())?,
            max_outcomes: whole(&limit.get_untracked(), "reference state limit")?,
            cds: selected_cds,
        };
        let revision = editing::revision(&molecule, &request).map_err(|e| e.to_string())?;
        Ok((
            molecule,
            Capture {
                source: pin,
                epoch: epoch.get_untracked(),
                request,
                revision,
            },
        ))
    });
    let guard = Callback::new(move |capture: Capture| -> Result<(), String> {
        if source.try_get_untracked().as_ref() != Some(&capture.source)
            || epoch.try_get_untracked() != Some(capture.epoch)
        {
            return Err("The source event or workspace changed; enumerate again.".into());
        }
        let (_, current) = prepare.run(())?;
        if current.request != capture.request || current.revision != capture.revision {
            return Err("The guide or base-editing settings changed; enumerate again.".into());
        }
        Ok(())
    });
    let run = {
        let alive = alive.clone();
        move |_| {
            if job.busy.get_untracked() {
                return;
            }
            let (molecule, capture) = match prepare.run(()) {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            let (worker, token) = match job.start() {
                Ok(value) => value,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
            status.set("Enumerating conditional base-editing reference states…".into());
            let alive = alive.clone();
            spawn_local(async move {
                let response = worker
                    .ask(Request {
                        id: token,
                        revision: capture.revision.clone(),
                        job: AnalysisJob::AnalyzeBaseEditing {
                            molecule,
                            request: capture.request.clone(),
                        },
                    })
                    .await;
                worker.close();
                if !alive.load(Ordering::Relaxed) || !job.current(token) {
                    return;
                }
                job.busy.set(false);
                job.worker.set_value(None);
                let checked = (|| {
                    guard.run(capture.clone())?;
                    match response?.result {
                        AnalysisResult::BaseEditing(value) => {
                            value
                                .validate_binding(&capture.revision)
                                .map_err(|e| e.to_string())?;
                            if value.request != capture.request
                                || value.source_revision != capture.source.binding.revision
                            {
                                return Err("The base-editing response binding changed.".into());
                            }
                            Ok(value)
                        }
                        AnalysisResult::Error(error) => Err(error),
                        _ => Err("Unexpected base-editing response.".into()),
                    }
                })();
                match checked {
                    Ok(value) => {
                        report.set(Some((capture, Arc::new(value))));
                        status.set("Reference states ready. Review intended edits, bystanders and unmodeled outcomes.".into());
                    }
                    Err(error) => status.set(error),
                }
            });
        }
    };
    let export_alive = alive.clone();
    let export = move |_| {
        let Some((capture, value)) = report.get_untracked() else {
            return;
        };
        let prepared = (|| {
            guard.run(capture.clone())?;
            let bytes = serde_json::to_vec_pretty(value.as_ref()).map_err(|e| e.to_string())?;
            ExportArtifact::prepare(
                "base-editing-reference-states.json",
                "application/json",
                &bytes,
            )
        })();
        match prepared {
            Ok(artifact) => {
                let alive = export_alive.clone();
                artifact.complete_guarded(
                    move || {
                        if !alive.load(Ordering::Relaxed) {
                            return Err("The base-editing panel was closed.".into());
                        }
                        guard.run(capture.clone())
                    },
                    move |result| {
                        if let Err(error) = result {
                            let _ = status.try_set(error);
                        }
                    },
                );
            }
            Err(error) => status.set(error),
        }
    };
    view! {
        <details class="crispr-base-editing"><summary>"Base editing · conditional reference states"</summary>
            <p class="helper">"Original BE1 cytosine-transition window, guide positions 4–8 from the PAM-distal end. Requires its supported SpCas9 geometry. No activity or efficiency prediction."</p>
            <label>"Intended guide positions (comma separated, optional)"<input maxlength="128"
                prop:value=move ||intended.get() on:input=move |e|intended.set(event_target_value(&e))/></label>
            <label>"Maximum retained reference states"<input inputmode="numeric" maxlength="4"
                prop:value=move ||limit.get() on:input=move |e|limit.set(event_target_value(&e))/></label>
            <details><summary>"Optional assumed CDS and codon consequences"</summary>
                <label>"Assume this CDS annotation for codon checking"
                    <select prop:value=move ||cds.get() on:change=move |e|cds.set(event_target_value(&e))>
                        <option value="">"No CDS assumption"</option>
                        {move ||get_document.run(source.get().binding.document_id).map(|m|m.features.iter().enumerate()
                            .filter(|(_,f)|f.kind.eq_ignore_ascii_case("CDS")&&!f.needs_review)
                            .map(|(i,_)|view! {<option value=i.to_string()>{format!("CDS annotation {}",i+1)}</option>}).collect_view())}
                    </select>
                </label>
                <label>"NCBI genetic code"<input inputmode="numeric" maxlength="3"
                    prop:value=move ||genetic_code.get() on:input=move |e|genetic_code.set(event_target_value(&e))/></label>
                <p class="helper">"Ordinary codon changes use the explicitly assumed annotation and code. Unsupported transcript paths and initiation effects remain unresolved."</p>
            </details>
            <button class="button secondary" disabled=move ||job.busy.get()||guide.get().is_none() on:click=run>"Enumerate reference states"</button>
            <p role="status" class="helper">{move ||status.get()}</p>
            {move ||report.get().map(|(capture,r)| {
                let total=r.states.len();
                let first=page.get().saturating_mul(25);
                let summary=format!("{} of {} finite reference states retained · finite enumeration complete: {} · biological outcomes exhaustive: {}",
                    total,r.finite_state_count,r.finite_enumeration_complete,r.biological_outcomes_exhaustive);
                let states=r.states.iter().enumerate().skip(first).take(25).map(|(index,state)| {
                    let consequence=state.codon_consequences.iter().map(|c|view! {<p class="helper">{
                        format!("CDS {} · codon {} · {} → {} · {} → {} · {:?}",c.feature_index+1,c.codon_number,c.before,c.after,c.before_amino_acid,c.after_amino_acid,c.kind)
                    }</p>}).collect_view();
                    view! {<article class="orf-entry"><strong>{format!("Reference state {}",index+1)}</strong>
                        <p>{format!("Edited guide positions {:?} · bystanders {:?} · intended changes satisfied: {}",state.edited_positions,state.bystander_positions,
                            match state.satisfies_intent {Some(true)=>"yes",Some(false)=>"no",None=>"not declared"})}</p>
                        <code>{state.reference_guide.clone()}</code>{consequence}
                    </article>}
                }).collect_view();
                let sites=r.sites.iter().map(|site| {
                    let pin=capture.source.binding.clone();let at=site.reference_position;let saved_epoch=capture.epoch;
                    view! {<p>{format!("Guide position {} · source base {} · {} → {}{}",site.guide_position,at+1,site.reference_from,site.reference_to,if site.intended{" · intended"}else if r.request.intended_positions.is_empty(){" · editable; no intent declared"}else{" · possible bystander"})}
                        <button class="text-tool" on:click=move |_|on_navigate.run(ReadNavigation{binding:pin.clone(),call_index:at,epoch:saved_epoch})>"Show base"</button></p>}
                }).collect_view();
                view! {<div><p>{summary}</p><p class="helper">{r.statement.clone()}</p>{sites}
                    {r.findings.iter().map(|finding|view! {<p>{format!("{:?} · {:?}: {}",finding.rule,finding.verdict,finding.message)}</p>}).collect_view()}
                    {r.cds.as_ref().map(|cds|view! {<p class="helper">{cds.reason.clone()}</p>})}
                    {states}
                    <div class="crispr-actions"><button class="text-tool" disabled=move ||page.get()==0 on:click=move |_|page.update(|p|*p=p.saturating_sub(1))>"Previous states"</button>
                        <button class="text-tool" disabled=move ||page.get().saturating_add(1).saturating_mul(25)>=total on:click=move |_|page.update(|p|*p+=1)>"Next states"</button>
                        <button class="text-tool" on:click=export.clone()>"Export complete evidence (JSON)"</button>
                    </div>
                </div>}
            })}
        </details>
    }
}
