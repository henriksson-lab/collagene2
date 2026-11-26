//! Source-bound PCR-product digest inspection. Only the requested fragment index
//! is applied; the complete recipe reconstructs material again before verification.
use crate::{
    cloning_methods_model::MethodsState, pcr_digest_model, source_identity::SourceInfo,
    worker::CachedWorker,
};
use collagene_core::{
    Molecule, catalog_digest,
    pcr_digest_inspection::{self, Preview},
    protocol::{AnalysisJob, AnalysisResult, Request},
};
use leptos::prelude::*;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
#[derive(Clone, PartialEq, Eq)]
struct Context {
    source: SourceInfo,
    document: usize,
    request_key: String,
    epoch: u64,
    origin: Option<crate::run_origin_model::PinnedOrigin>,
}
struct Saved {
    context: Context,
    revision: String,
    preview: Preview,
}
fn request_key(request: &pcr_digest_inspection::Request) -> Result<String, String> {
    serde_json::to_string(request).map_err(|e| e.to_string())
}
#[component]
pub fn PcrDigestPicker(
    config: RwSignal<MethodsState>,
    index: usize,
    get_query: Callback<usize, Option<Molecule>>,
    source_info: Callback<usize, Option<SourceInfo>>,
    epoch: RwSignal<u64>,
) -> impl IntoView {
    let row = move || config.with(|c| c.rows.get(index).cloned().unwrap_or_default());
    let catalog = StoredValue::new(catalog_digest::builtin_catalog().ok());
    let search = RwSignal::new(String::new());
    let status = RwSignal::new(String::new());
    let proof = RwSignal::new(None::<Arc<Saved>>);
    let busy = RwSignal::new(false);
    let page = RwSignal::new(0usize);
    let geometry_page = RwSignal::new(0usize);
    let generation = RwSignal::new(0u32);
    let active = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = active.get_value() {
            worker.close();
        }
    });
    let context = Memo::new(move |_| {
        let state = config.get();
        let row = state.rows.get(index)?;
        let request = pcr_digest_model::build_request(&state, index).ok()?;
        Some(Context {
            source: source_info.run(row.document)?,
            document: row.document,
            request_key: request_key(&request).ok()?,
            epoch: epoch.get(),
            origin: row.run_origin.clone(),
        })
    });
    Effect::new(move |_| {
        let _ = context.get();
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(worker) = active.get_value() {
            worker.close();
        }
        active.set_value(None);
        proof.set(None);
        busy.set(false);
        page.set(0);
        geometry_page.set(0);
    });
    let toggle = Callback::new(move |(id, checked): (String, bool)| {
        let selected = row().catalog_selection;
        let mut ids = selected
            .as_ref()
            .map(|s| s.enzymes.iter().map(|e| e.id.clone()).collect::<Vec<_>>())
            .unwrap_or_default();
        if checked && !ids.contains(&id) {
            ids.push(id);
        } else if !checked {
            ids.retain(|selected| selected != &id);
        }
        if ids.is_empty() {
            config.update(|c| {
                if let Some(row) = c.rows.get_mut(index) {
                    row.catalog_selection = None;
                }
            });
            return;
        }
        let selection=catalog.with_value(|catalog|catalog.as_ref().ok_or("The built-in enzyme catalog is unavailable.").and_then(|catalog|catalog_digest::select(catalog,&ids).map_err(|_|"The selected enzyme set exceeds the catalog bounds or contains an unsupported entry.")));
        match selection {
            Ok(selection) => config.update(|c| {
                if let Some(row) = c.rows.get_mut(index) {
                    row.catalog_selection = Some(selection);
                }
            }),
            Err(error) => status.set(error.into()),
        }
    });
    let compute = Callback::new(move |()| {
        let request =
            match config.with_untracked(|state| pcr_digest_model::build_request(state, index)) {
                Ok(request) => request,
                Err(error) => {
                    status.set(error);
                    return;
                }
            };
        let Some(bound) = context.get_untracked() else {
            status.set("Choose an available original source and exact PCR primers.".into());
            return;
        };
        if bound.source.length > 8 * 1024 * 1024 {
            status.set("The original source exceeds the8MiB preparation profile.".into());
            return;
        }
        let Some(source) = get_query.run(bound.document) else {
            status.set("The pinned original template is unavailable.".into());
            return;
        };
        if source.revision() != bound.source.binding.revision
            || request_key(&request).ok() != Some(bound.request_key.clone())
            || context.get_untracked() != Some(bound.clone())
        {
            status.set("The original source or PCR request changed; inspect again.".into());
            return;
        }
        if bound
            .origin
            .as_ref()
            .is_some_and(|origin| origin.selection.event != bound.source.binding)
        {
            status.set("The saved-material origin belongs to a different source event.".into());
            return;
        }
        let revision = match bound.origin.as_ref().map_or_else(
            || pcr_digest_inspection::inspection_revision(&source, &request),
            |origin| {
                pcr_digest_inspection::inspection_revision_from_run(
                    &source,
                    origin.as_ref(),
                    &request,
                )
            },
        ) {
            Ok(revision) => revision,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let worker = match CachedWorker::new() {
            Ok(worker) => Rc::new(worker),
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Some(old) = active.get_value() {
            old.close();
        }
        active.set_value(Some(worker.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        proof.set(None);
        status.set("Reconstructing exact PCR and the selected catalog digest…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let job = if let Some(origin) = &bound.origin {
                AnalysisJob::InspectRunPcrDigest {
                    molecule: source,
                    origin: origin.as_ref().clone(),
                    request,
                }
            } else {
                AnalysisJob::InspectPcrDigest {
                    molecule: source,
                    request,
                }
            };
            let response = worker
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
            worker.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || context.get_untracked() != Some(bound.clone())
            {
                return;
            }
            active.set_value(None);
            busy.set(false);
            let result = match response {
                Ok(response) if response.id == token && response.revision == revision => {
                    match response.result {
                        AnalysisResult::PcrDigestInspection(preview) if bound.origin.is_none() => {
                            Ok(preview)
                        }
                        AnalysisResult::RunPcrDigestInspection(report)
                            if bound.origin.is_some() =>
                        {
                            report
                                .validate_binding(&revision)
                                .map(|()| report.preview)
                                .map_err(|e| e.to_string())
                        }
                        AnalysisResult::Error(error) => Err(error),
                        _ => Err("Unexpected PCR digest response.".into()),
                    }
                }
                Ok(_) => Err("Stale PCR digest response withheld.".into()),
                Err(error) => Err(error),
            };
            match result {
                Ok(preview) => {
                    let preview_revision = if bound.origin.is_some() {
                        preview.context_revision.clone()
                    } else {
                        revision.clone()
                    };
                    if let Err(error) = preview.validate_binding(&preview_revision) {
                        status.set(error.to_string());
                        return;
                    }
                    if preview.source_revision != bound.source.binding.revision
                        || request_key(&preview.request).ok() != Some(bound.request_key.clone())
                    {
                        status.set("The inspection differs from the pinned original or primer/catalog request.".into());
                        return;
                    }
                    let selectable = preview.pcr_geometries.len() == 1
                        && preview
                            .inspection
                            .as_ref()
                            .is_some_and(|i| i.scan_complete && i.products_determined);
                    status.set(if selectable{"Choose a reconstructed digest fragment below. The complete assembly will independently reconstruct it again."}else{"Fragment selection is withheld. Review PCR alternatives and the reported limits or unavailable chemistry."}.into());
                    proof.set(Some(Arc::new(Saved {
                        context: bound,
                        revision: preview_revision,
                        preview,
                    })));
                }
                Err(error) => status.set(error),
            }
        });
    });
    let choose = Callback::new(move |fragment: usize| {
        let Some(saved) = proof.get_untracked() else {
            return;
        };
        if context.get_untracked() != Some(saved.context.clone()) {
            status.set(
                "Source, primer or catalog settings changed. Inspect the PCR digest again.".into(),
            );
            return;
        }
        let updated = config.with_untracked(|state| {
            pcr_digest_model::apply(
                state,
                index,
                &saved.preview.request,
                &saved.preview,
                &saved.revision,
                fragment,
            )
        });
        match updated {
            Ok(updated) => {
                config.set(updated);
                status.set(format!("Fragment {} requested. Reconstruct and verify the complete prepared assembly against the independent target.",fragment+1));
            }
            Err(error) => status.set(error),
        }
    });
    view!{<section class="pcr-digest-picker"><details><summary>{move||format!("PCR digest enzymes · {} selected",row().catalog_selection.as_ref().map(|s|s.enzymes.len()).unwrap_or(0))}</summary>
        <label>"Find enzymes"<input type="search" aria-label=format!("PCR digest enzyme search {}",index+1) prop:value=move||search.get() on:input=move|e|search.set(event_target_value(&e))/></label>
        {move||{
            let needle=search.get().to_ascii_lowercase();let selected=row().catalog_selection.map(|s|s.enzymes.iter().map(|e|e.id.clone()).collect::<Vec<_>>()).unwrap_or_default();
            catalog.with_value(|catalog|catalog.as_ref().map(|catalog|{
                let all=catalog.edition().enzymes.iter().filter(|e|needle.is_empty()&&selected.contains(&e.id)||!needle.is_empty()&&e.name.to_ascii_lowercase().contains(&needle)).take(61).collect::<Vec<_>>();
                let capped=all.len()>60;
                view!{<div class="pcr-digest-enzyme-choices">{all.into_iter().take(60).map(|entry|{let id=entry.id.clone();let checked=selected.contains(&id);view!{<label class="checkbox-label"><input type="checkbox" prop:checked=checked on:change=move|e|toggle.run((id.clone(),event_target_checked(&e)))/>{entry.name.clone()}</label>}}).collect_view()}<Show when=move||capped><p class="helper">"First60 matching enzymes shown. Refine the name search."</p></Show></div>}
            }))
        }}
        <button class="text-tool" disabled=move||row().catalog_selection.is_none() on:click=move|_|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.catalog_selection=None;})>"Clear enzyme selection"</button>
        <p class="helper">"Search by name and select each enzyme needed for the digest. The exact catalog identities and cleavage geometries are pinned together."</p></details>
        <p>{move||row().catalog_selection.map(|selection|selection.enzymes.iter().map(|e|e.enzyme.name.clone()).collect::<Vec<_>>().join(" + ")).unwrap_or_else(||"Choose the PCR-product digest enzymes.".into())}</p>
        <button class="button secondary" disabled=move||busy.get() on:click=move|_|compute.run(())>"Inspect PCR digest fragments"</button>
        <p class="helper" role="status">{move||status.get()}</p>
        {move||proof.get().map(|saved|{
            let preview=&saved.preview;let count=preview.pcr_geometries.len();let start=geometry_page.get().min(count.saturating_sub(1)/20)*20;
            let geometries=preview.pcr_geometries.iter().skip(start).take(20).map(|g|view!{<tr><td>{g.index+1}</td><td>{g.product_length}</td><td>{format!("{} / {}",g.forward_start+1,g.reverse_start+1)}</td></tr>}).collect_view();
            let selectable=count==1&&preview.amplicon.is_some()&&preview.inspection.as_ref().is_some_and(|i|i.scan_complete&&i.products_determined);
            view!{<section class="pcr-digest-result"><p>{format!("{} exact PCR geometries",count)}</p>
                {preview.findings.iter().map(|finding|view!{<p class="helper">{finding.clone()}</p>}).collect_view()}
                {(count>0).then(||view!{<details open={count>1}><summary>"PCR geometry review"</summary><table><thead><tr><th>"Geometry"</th><th>"Product nt"</th><th>"Template binding starts · 1-based"</th></tr></thead><tbody>{geometries}</tbody></table>
                    <div class="method-row-actions"><button class="text-tool" disabled=move||geometry_page.get()==0 on:click=move|_|geometry_page.update(|p|*p=p.saturating_sub(1))>"Previous PCR geometries"</button><button class="text-tool" disabled=move||(geometry_page.get()+1)*20>=count on:click=move|_|geometry_page.update(|p|*p=p.saturating_add(1))>"Next PCR geometries"</button></div>
                </details>})}
                {preview.amplicon.as_ref().map(|amplicon|view!{<p>{format!("Reconstructed amplicon · {} nt",amplicon.sequence.len())}</p>})}
                {preview.inspection.as_ref().map(|inspection|{
                    let count=inspection.fragments.len();let start=page.get().min(count.saturating_sub(1)/20)*20;
                    view!{<section><p>{format!("Recognition scan complete {} · products determined {} · {} sites · {} fragments",inspection.scan_complete,inspection.products_determined,inspection.sites.len(),count)}</p>
                        <p class="helper">{format!("Recognition comparisons: {} · limits: {}",inspection.comparisons,if inspection.limits.is_empty(){"none reached".into()}else{inspection.limits.iter().map(|limit|match limit{collagene_core::restriction_inspection::InspectionLimit::Comparisons=>"comparison budget",collagene_core::restriction_inspection::InspectionLimit::Sites=>"site count",collagene_core::restriction_inspection::InspectionLimit::Fragments=>"fragment count",collagene_core::restriction_inspection::InspectionLimit::ReportBytes=>"report bytes"}).collect::<Vec<_>>().join(", ")})}</p>
                        {inspection.diagnostics.iter().map(|message|view!{<p class="helper">{message.clone()}</p>}).collect_view()}
                        <p class="helper">"Ends and coordinates use the PCR product’s top-strand reference, before assembly orientation, end repair and terminal treatment. Full reconstruction checks those later steps."</p>
                        {inspection.fragments.iter().skip(start).take(20).map(|fragment|{let i=fragment.index;view!{<details><summary>{format!("Fragment {} · top {} nt / bottom {} nt · starts at amplicon top boundary {} (zero-based)",i+1,fragment.duplex.top_length,fragment.duplex.bottom_length,fragment.top_start)}</summary>{crate::end_view::summary(&fragment.duplex)}<button class="button secondary" disabled=!selectable on:click=move|_|choose.run(i)>{move||if config.with(|c|c.rows.get(index).and_then(|row|row.fragment.parse::<usize>().ok())==Some(i+1)){"Current fragment request"}else{"Use this fragment"}}</button></details>}}).collect_view()}
                        <div class="method-row-actions"><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous fragments"</button><button class="text-tool" disabled=move||(page.get()+1)*20>=count on:click=move|_|page.update(|p|*p=p.saturating_add(1))>"Next fragments"</button></div>
                    </section>}
                })}
            </section>}
        })}
        <details><summary>"Advanced declared fragment index"</summary><label>"Fragment number · 1-based"<input type="number" min="1" aria-label=format!("PCR catalog fragment {}",index+1) prop:value=move||row().fragment on:input=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.fragment=event_target_value(&e)})/></label><p class="helper">"This is a request for fresh reconstruction, including fragments pinned by applied Type IIS primer proposals."</p></details>
    </section>}.into_any()
}
