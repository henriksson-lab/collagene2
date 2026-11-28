//! Source-bound PCR-product digest inspection. Only the requested fragment index
//! is applied; the complete recipe reconstructs material again before verification.
use crate::cloning::preparation::cloning_methods_model::MethodsState;
use crate::cloning::preparation::pcr_digest_model;
use crate::cloning::preparation::pcr_digest_model::InspectionPreview;
use crate::cloning::preparation::pcr_digest_model::InspectionRequest;
use crate::platform::worker::CachedWorker;
use crate::sources::identity::SourceInfo;
use collagene_core::{
    cloning::restriction::catalog as catalog_digest,
    protocol::{request::Request, schema::AnalysisResult},
};
use collagene_model::Molecule;
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
    origin: Option<crate::history::run_origin_model::PinnedOrigin>,
}
struct Saved {
    context: Context,
    revision: String,
    preview: InspectionPreview,
}
fn request_key(request: &InspectionRequest) -> Result<String, String> {
    request.key()
}
#[component]
pub fn PcrDigestPicker(
    config: RwSignal<MethodsState>,
    index: usize,
    get_query: Callback<usize, Option<Molecule>>,
    get_workspace: Option<
        Callback<(), Result<collagene_core::documents::workspace::Workspace, String>>,
    >,
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
        let request = pcr_digest_model::build_inspection_request(&state, index).ok()?;
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
        let request = match config
            .with_untracked(|state| pcr_digest_model::build_inspection_request(state, index))
        {
            Ok(request) => request,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        let Some(bound) = context.get_untracked() else {
            status.set("Choose an available original source and explicit PCR primers.".into());
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
            .is_some_and(|origin| *origin.event() != bound.source.binding)
        {
            status.set("The saved-material origin belongs to a different source event.".into());
            return;
        }
        let origin = bound.origin.as_ref().map(|o| o.as_ref().clone());
        let workspace = if origin.is_some() {
            let Some(get_workspace) = get_workspace else {
                status.set("The exact saved-material workspace is unavailable.".into());
                return;
            };
            match get_workspace.run(()) {
                Ok(workspace) => Some(workspace),
                Err(error) => {
                    status.set(error);
                    return;
                }
            }
        } else {
            None
        };
        let revision = match request.revision(&source, origin.as_ref(), workspace.as_ref()) {
            Ok(revision) => revision,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let job = match request.clone().into_job(source, origin, workspace) {
            Ok(job) => job,
            Err(error) => {
                status.set(error);
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
        status.set("Reconstructing the declared preparation and pinned catalog digest…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let mismatch = matches!(&request, InspectionRequest::Mismatch(_));
            let extension = matches!(&request, InspectionRequest::Extension(_));
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
                        AnalysisResult::PcrDigestInspection(preview)
                            if !mismatch && !extension && bound.origin.is_none() =>
                        {
                            Ok(InspectionPreview::Exact(preview))
                        }
                        AnalysisResult::MaterialPcrDigestInspection(report)
                            if !mismatch && !extension && bound.origin.is_some() =>
                        {
                            report
                                .validate_binding(&revision)
                                .map(|()| InspectionPreview::Material(report))
                                .map_err(|e| e.to_string())
                        }
                        AnalysisResult::MismatchPcrDigestInspection(report) if mismatch => report
                            .validate_binding(&revision)
                            .map(|()| InspectionPreview::Mismatch(report))
                            .map_err(|e| e.to_string()),
                        AnalysisResult::OverlapExtensionDigestInspection(report)
                            if extension
                                && report.origin_revision.is_some() == bound.origin.is_some() =>
                        {
                            report
                                .validate_binding(&revision)
                                .map(|()| InspectionPreview::Extension(report))
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
                    let preview_revision = preview.context_revision().to_owned();
                    if let Err(error) = preview.validate_binding(&preview_revision) {
                        status.set(error.to_string());
                        return;
                    }
                    if preview.source_revision() != bound.source.binding.revision
                        || request_key(&preview.request()).ok() != Some(bound.request_key.clone())
                    {
                        status.set("The inspection differs from the pinned original or primer/catalog request.".into());
                        return;
                    }
                    let selectable = preview.selectable();
                    status.set(if !selectable {"Fragment selection is withheld. Review PCR alternatives and the reported limits or unavailable chemistry."} else if preview.incomplete() {"Conditional geometry only: the PCR or saved-origin search is incomplete. Choosing a fragment keeps these limits; final verification retains the unresolved alternatives."} else {"Choose a reconstructed digest fragment below. The complete assembly will independently reconstruct it again."}.into());
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
            pcr_digest_model::apply_inspection(
                state,
                index,
                &saved.preview.request(),
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
            let preview=&saved.preview;let count=preview.geometry_count();let start=geometry_page.get().min(count.saturating_sub(1)/20)*20;
            let geometries=preview.geometry_page(start).iter().map(|g|view!{<tr><td>{g.index+1}</td><td>{g.product_length}</td><td>{format!("{:?}: {} / {}",g.template_strand,g.forward_start+1,g.reverse_start+1)}</td></tr>}).collect_view();
            let selectable=preview.selectable();
            let extension=matches!(preview,InspectionPreview::Extension(_));
            view!{<section class="pcr-digest-result"><p>{if extension {format!("Overlap-extension preparation · {} reconstructed outer-PCR product",count)} else {format!("{} {} PCR geometries",count,if matches!(preview,InspectionPreview::Mismatch(_)){"intentional-mismatch"}else{"exact"})}}</p>
                {preview.incomplete().then(||view!{<p class="helper" role="status">"Conditional geometry: PCR alternatives or saved-origin evidence remain incomplete. Selecting a fragment does not resolve those limitations."</p>})}
                {preview.findings().iter().take(40).map(|finding|view!{<p class="helper">{finding.clone()}</p>}).collect_view()}
                {if let InspectionPreview::Mismatch(report)=preview {
                    let pcr=&report.pcr.report;
                    view!{<details><summary>"Mismatch policy and introduced bases"</summary>
                        <p>{format!("At most {} mismatches / {} exact 3′ bases per primer · search complete {} · enumeration complete {} · conditional {}",report.request.pcr_options.binding.max_mismatches,report.request.pcr_options.binding.min_three_prime_exact,pcr.search_complete,pcr.enumeration_complete,report.incomplete)}</p>
                        {pcr.products.iter().enumerate().skip(start).take(20).map(|(i,product)|view!{<details><summary>{format!("Geometry {} · {} introduced substitutions · certain {} · overlapping domains {}",i+1,product.substitutions.len(),product.certain,product.overlapping_domains)}</summary>
                            <p>{format!("Forward/reverse mismatch counts: {}/{} · actual exact 3′ matches: {}/{}",product.forward_binding.mismatches.len(),product.reverse_binding.mismatches.len(),product.forward_binding.three_prime_exact,product.reverse_binding.three_prime_exact)}</p>
                            {product.substitutions.iter().take(40).map(|change|view!{<p>{format!("Template {}: {} → {} at amplicon {}{}",change.source_position+1,change.template_base,change.product_base,change.product_position+1,if change.template_ambiguous{" (ambiguous template)"}else{""})}</p>}).collect_view()}
                            {(product.substitutions.len()>40).then(||view!{<p>"First40 introduced bases shown for this geometry."</p>})}
                        </details>}).collect_view()}
                    </details>}.into_any()
                }else{().into_any()}}
                {if let InspectionPreview::Extension(report)=preview {
                    view!{<details><summary>"Actual PCR stages and terminal pairing"</summary>
                        {report.extension.as_ref().map(|geometry|view!{<p>{format!("Selected terminal pairing {} bp · complete alternative search {} · unique within terminal model {}",geometry.options.overlap,geometry.terminal_pairings.complete,geometry.terminal_pairings.unique)}</p>})}
                        {report.pcr_stages.iter().map(|stage|view!{<details><summary>{format!("PCR stage {} · {} geometries · {}",stage.stage,stage.geometries.len(),stage.product_length.map(|n|format!("{n} nt product")).unwrap_or_else(||"product unavailable".into()))}</summary>
                            <p class="helper">"Binding starts address this PCR stage’s actual template, before outer digestion and assembly orientation."</p>
                            <p class="helper">{format!("Template revision {}",stage.template_revision)}</p>
                            <table><thead><tr><th>"Geometry"</th><th>"Product nt"</th><th>"Binding starts · 1-based"</th></tr></thead><tbody>{stage.geometries.iter().take(20).map(|geometry|view!{<tr><td>{geometry.index+1}</td><td>{geometry.product_length}</td><td>{format!("{:?}: {} / {}",geometry.template_strand,geometry.forward_start+1,geometry.reverse_start+1)}</td></tr>}).collect_view()}</tbody></table>
                            {(stage.geometries.len()>20).then(||view!{<p class="helper">"First20 geometries shown for this stage; the full retained stage evidence remains bound to this inspection."</p>})}
                        </details>}).collect_view()}
                    </details>}.into_any()
                } else {().into_any()}}
                {(!extension&&count>0).then(||view!{<details open={count>1}><summary>"PCR geometry review"</summary><table><thead><tr><th>"Geometry"</th><th>"Product nt"</th><th>"Template binding starts · 1-based"</th></tr></thead><tbody>{geometries}</tbody></table>
                    <div class="method-row-actions"><button class="text-tool" disabled=move||geometry_page.get()==0 on:click=move|_|geometry_page.update(|p|*p=p.saturating_sub(1))>"Previous PCR geometries"</button><button class="text-tool" disabled={move||(geometry_page.get()+1)*20>=count} on:click=move|_|geometry_page.update(|p|*p=p.saturating_add(1))>"Next PCR geometries"</button></div>
                </details>})}
                {preview.amplicon().map(|amplicon|view!{<p>{format!("Reconstructed amplicon · {} nt",amplicon.sequence.len())}</p>})}
                {preview.inspection().map(|inspection|{
                    let count=inspection.fragments.len();let start=page.get().min(count.saturating_sub(1)/20)*20;
                    view!{<section><p>{format!("Recognition scan complete {} · products determined {} · {} sites · {} fragments",inspection.scan_complete,inspection.products_determined,inspection.sites.len(),count)}</p>
                        <p class="helper">{format!("Recognition comparisons: {} · limits: {}",inspection.comparisons,if inspection.limits.is_empty(){"none reached".into()}else{inspection.limits.iter().map(|limit|match limit{collagene_core::cloning::restriction::inspection::InspectionLimit::Comparisons=>"comparison budget",collagene_core::cloning::restriction::inspection::InspectionLimit::Sites=>"site count",collagene_core::cloning::restriction::inspection::InspectionLimit::Fragments=>"fragment count",collagene_core::cloning::restriction::inspection::InspectionLimit::ReportBytes=>"report bytes"}).collect::<Vec<_>>().join(", ")})}</p>
                        {inspection.diagnostics.iter().map(|message|view!{<p class="helper">{message.clone()}</p>}).collect_view()}
                        <p class="helper">"Ends and coordinates use the PCR product’s top-strand reference, before assembly orientation, end repair and terminal treatment. Full reconstruction checks those later steps."</p>
                        {inspection.fragments.iter().skip(start).take(20).map(|fragment|{let i=fragment.index;let conditional=preview.incomplete();view!{<details><summary>{format!("Fragment {} · top {} nt / bottom {} nt · starts at amplicon top boundary {} (zero-based)",i+1,fragment.duplex.top_length,fragment.duplex.bottom_length,fragment.top_start)}</summary>{crate::cloning::preparation::end_view::summary(&fragment.duplex)}<button class="button secondary" disabled=!selectable on:click=move|_|choose.run(i)>{move||if config.with(|c|c.rows.get(index).and_then(|row|row.fragment.parse::<usize>().ok())==Some(i+1)){if conditional{"Current conditional fragment request"}else{"Current fragment request"}}else{if conditional{"Use conditional fragment"}else{"Use this fragment"}}}</button></details>}}).collect_view()}
                        <div class="method-row-actions"><button class="text-tool" disabled=move||page.get()==0 on:click=move|_|page.update(|p|*p=p.saturating_sub(1))>"Previous fragments"</button><button class="text-tool" disabled={move||(page.get()+1)*20>=count} on:click=move|_|page.update(|p|*p=p.saturating_add(1))>"Next fragments"</button></div>
                    </section>}
                })}
            </section>}
        })}
        <details><summary>"Advanced declared fragment index"</summary><label>"Fragment number · 1-based"<input type="number" min="1" aria-label=format!("PCR catalog fragment {}",index+1) prop:value=move||row().fragment on:input=move|e|config.update(|c|if let Some(row)=c.rows.get_mut(index){row.fragment=event_target_value(&e)})/></label><p class="helper">"This is a request for fresh reconstruction, including fragments pinned by applied Type IIS primer proposals."</p></details>
    </section>}.into_any()
}
