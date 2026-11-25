//! One source-bound recognition scan shared by sequence, map and enzyme picker.
use crate::worker::CachedWorker;
use collagene_core::{
    Alphabet, Molecule,
    analysis::Enzyme,
    catalog_digest, enzyme_sets,
    protocol::{AnalysisJob, AnalysisResult, Request},
    restriction_inspection::{self, InspectionOptions, InspectionSite},
    restriction_selection::{self, Basis, Region, Relation},
};
use leptos::prelude::*;
use std::{
    collections::BTreeMap,
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen::{JsCast, closure::Closure};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone, Default, PartialEq, Eq, serde::Serialize)]
pub struct SiteCount {
    pub certain: usize,
    pub possible: usize,
    pub complete: bool,
}
impl SiteCount {
    fn label(&self) -> String {
        if !self.complete {
            return if self.certain + self.possible == 0 {
                "Unknown".into()
            } else {
                format!("≥{} · incomplete", self.certain + self.possible)
            };
        }
        if self.possible > 0 {
            format!("{} + {} possible", self.certain, self.possible)
        } else {
            self.certain.to_string()
        }
    }
    fn accepts(&self, filter: &str) -> bool {
        if filter == "all" {
            return true;
        }
        if filter == "uncertain" {
            return !self.complete || self.possible > 0;
        }
        self.complete
            && self.possible == 0
            && match filter {
                "0" => self.certain == 0,
                "1" => self.certain == 1,
                "2" => self.certain == 2,
                "many" => self.certain >= 3,
                _ => false,
            }
    }
}
#[derive(Clone)]
pub struct Scan {
    pub context: String,
    pub revision: String,
    pub sites: Arc<Vec<InspectionSite>>,
    pub counts: BTreeMap<String, SiteCount>,
    pub reports: Arc<Vec<restriction_inspection::RestrictionInspection>>,
}
struct Timer {
    id: i32,
    _callback: Closure<dyn FnMut()>,
}
impl Drop for Timer {
    fn drop(&mut self) {
        if let Some(w) = web_sys::window() {
            w.clear_timeout_with_handle(self.id);
        }
    }
}
#[derive(Clone)]
struct PinnedRegion {
    context: String,
    revision: String,
    region: Region,
}

#[derive(Clone, Copy)]
pub struct Restrictions {
    pub visible: RwSignal<bool>,
    pub names: RwSignal<Vec<String>>,
    pub count_filter: RwSignal<&'static str>,
    pub region_filter: RwSignal<&'static str>,
    pub basis: RwSignal<Basis>,
    pub region_status: Memo<String>,
    scoped_counts: Memo<BTreeMap<String, SiteCount>>,
    pinned: RwSignal<Option<PinnedRegion>>,
    source: Memo<Molecule>,
    context: Memo<String>,
    selection: Memo<Option<Vec<(usize, usize)>>>,
    pub selected: RwSignal<Option<usize>>,
    pub scan: RwSignal<Option<Arc<Scan>>>,
    pub displayed: Memo<Option<Arc<Vec<InspectionSite>>>>,
    pub status: RwSignal<String>,
    pub busy: RwSignal<bool>,
    pub catalog: StoredValue<Vec<Enzyme>>,
    pub refresh: RwSignal<u64>,
}
impl Restrictions {
    pub fn count(&self, name: &str) -> SiteCount {
        self.scoped_counts
            .with(|counts| counts.get(name).cloned().unwrap_or_default())
    }
    pub fn capture_selection(&self) {
        let m = self.source.get_untracked();
        let region = self
            .selection
            .get_untracked()
            .and_then(|spans| Region::with_topology(m.sequence.len(), m.topology, spans).ok());
        self.pinned.set(region.map(|region| PinnedRegion {
            context: self.context.get_untracked(),
            revision: m.revision(),
            region,
        }));
    }

    pub fn accepts(&self, name: &str) -> bool {
        self.count(name).accepts(self.count_filter.get())
    }
}
pub fn use_restrictions(
    molecule: Memo<Molecule>,
    context: Memo<String>,
    locked: Memo<bool>,
    selection: Memo<Option<Vec<(usize, usize)>>>,
) -> Restrictions {
    let catalog_result = catalog_digest::builtin_catalog()
        .map(|c| {
            let mut enzymes = c
                .edition()
                .enzymes
                .iter()
                .map(|e| Enzyme {
                    name: e.name.clone(),
                    recognition: e.geometry.recognition.clone(),
                    top_cut: e.geometry.top_cut,
                    bottom_cut: e.geometry.bottom_cut,
                })
                .collect::<Vec<_>>();
            enzymes.sort_by(|a, b| a.name.cmp(&b.name));
            enzymes
        })
        .map_err(|e| e.to_string());
    let catalog = StoredValue::new(catalog_result.clone().unwrap_or_default());
    let visible = RwSignal::new(true);
    let names = RwSignal::new(vec![
        "EcoRI".into(),
        "BamHI".into(),
        "HindIII".into(),
        "BsaI".into(),
    ]);
    let count_filter = RwSignal::new("all");
    let region_filter = RwSignal::new("all");
    let basis = RwSignal::new(Basis::RecognitionBases);
    let pinned = RwSignal::new(None::<PinnedRegion>);
    let selected = RwSignal::new(None);
    let scan = RwSignal::new(None::<Arc<Scan>>);
    let status = RwSignal::new("Preparing restriction sites…".into());
    let busy = RwSignal::new(false);
    let refresh = RwSignal::new(0u64);
    let generation = RwSignal::new(0u64);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let timer = StoredValue::new_local(None::<Timer>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        timer.set_value(None);
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
    });
    let revision = Memo::new(move |_| molecule.with(Molecule::revision));
    Effect::new(move |_| {
        let current_context = context.get();
        let current_revision = revision.get();
        if pinned.with_untracked(|p| {
            p.as_ref()
                .is_some_and(|p| p.context != current_context || p.revision != current_revision)
        }) {
            pinned.set(None);
        }
    });
    let region = Memo::new(move |_| {
        pinned
            .get()
            .filter(|p| p.context == context.get() && p.revision == revision.get())
            .map(|p| p.region)
    });
    let region_status = Memo::new(move |_| {
        if region_filter.get() == "all" {
            return "Counts cover the whole source.".into();
        }
        match region.get() {
            Some(r) => format!("Pinned bases {}. Counts cover the chosen relation; site navigation keeps this region.", r.spans().iter().map(|(a,b)|format!("{}–{}",a+1,b)).collect::<Vec<_>>().join(", ")),
            None => "Select a nonempty exact range, then use current selection. Scoped counts are unknown.".into(),
        }
    });
    let scoped_sites = Memo::new(move |_| {
        let snapshot = scan.get()?;
        if snapshot.context != context.get() || snapshot.revision != revision.get() {
            return None;
        }
        let relation = region_filter.get();
        if relation == "all" {
            return Some(snapshot.sites.clone());
        }
        let Some(region) = region.get() else {
            return Some(Arc::new(Vec::new()));
        };
        let basis = basis.get();
        Some(Arc::new(
            snapshot
                .sites
                .iter()
                .filter(|site| {
                    let classified = restriction_selection::classify(site, &region, &basis)
                        .unwrap_or(Relation::Undetermined);
                    matches!(
                        (relation, classified),
                        ("inside", Relation::Inside)
                            | ("outside", Relation::Outside)
                            | ("crossing", Relation::Crossing)
                            | ("undetermined", Relation::Undetermined)
                    )
                })
                .cloned()
                .collect(),
        ))
    });
    let scoped_counts = Memo::new(move |_| {
        let Some(snapshot) = scan.get() else {
            return BTreeMap::new();
        };
        if snapshot.context != context.get() || snapshot.revision != revision.get() {
            return BTreeMap::new();
        }
        if region_filter.get() == "all" {
            return snapshot.counts.clone();
        }
        let valid = region.get().is_some();
        let mut counts = snapshot
            .counts
            .iter()
            .map(|(name, count)| {
                (
                    name.clone(),
                    SiteCount {
                        complete: valid && count.complete,
                        ..SiteCount::default()
                    },
                )
            })
            .collect::<BTreeMap<_, _>>();
        if valid && region_filter.get() != "undetermined" {
            if let Some(region) = region.get() {
                for site in snapshot.sites.iter() {
                    if restriction_selection::classify(site, &region, &basis.get())
                        .unwrap_or(Relation::Undetermined)
                        == Relation::Undetermined
                    {
                        if let Some(count) = counts.get_mut(&site.enzyme) {
                            count.complete = false;
                        }
                    }
                }
            }
        }
        if let Some(sites) = scoped_sites.get() {
            for site in sites.iter() {
                if let Some(count) = counts.get_mut(&site.enzyme) {
                    if site.certain {
                        count.certain += 1;
                    } else {
                        count.possible += 1;
                    }
                }
            }
        }
        counts
    });
    let displayed = Memo::new(move |_| {
        if !visible.get() {
            return None;
        }
        let sites = scoped_sites.get()?;
        let selected_names = names.get();
        let filter = count_filter.get();
        Some(Arc::new(
            sites
                .iter()
                .filter(|site| {
                    selected_names.contains(&site.enzyme)
                        && scoped_counts.with(|counts| {
                            counts.get(&site.enzyme).is_some_and(|n| n.accepts(filter))
                        })
                })
                .cloned()
                .collect(),
        ))
    });
    Effect::new(move |_| {
        let _ = displayed.get();
        selected.set(None);
    });
    let effect_alive = alive.clone();
    Effect::new(move |_| {
        let expected_context = context.get();
        let expected_revision = revision.get();
        let _ = refresh.get();
        let is_locked = locked.get();
        let next = generation.get_untracked().wrapping_add(1);
        generation.set(next);
        timer.set_value(None);
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        scan.set(None);
        selected.set(None);
        busy.set(false);
        if is_locked {
            return;
        }
        if let Err(error) = &catalog_result {
            status.set(error.clone());
            return;
        }
        if molecule.with_untracked(|m| m.sequence.alphabet() != Alphabet::Dna) {
            status.set("Restriction enzymes apply to DNA documents.".into());
            return;
        }
        status.set("Scanning restriction recognition and strand-cut positions…".into());
        busy.set(true);
        let active = effect_alive.clone();
        let callback = Closure::<dyn FnMut()>::new(move || {
            if !active.load(Ordering::Relaxed) || generation.try_get_untracked() != Some(next) {
                return;
            }
            let m = molecule.get_untracked();
            let enzymes = catalog.get_value();
            let requested_context = expected_context.clone();
            let requested_revision = expected_revision.clone();
            let w = match CachedWorker::new() {
                Ok(w) => Rc::new(w),
                Err(e) => {
                    status.set(e);
                    busy.set(false);
                    return;
                }
            };
            worker.set_value(Some(w.clone()));
            let active = active.clone();
            spawn_local(async move {
                let mut sites = Vec::new();
                let mut reports = Vec::new();
                let mut counts = enzymes
                    .iter()
                    .map(|e| (e.name.clone(), SiteCount::default()))
                    .collect::<BTreeMap<_, _>>();
                let mut complete = true;
                let mut failure = None;
                for (batch, chunk) in enzymes.chunks(64).enumerate() {
                    if !active.load(Ordering::Relaxed)
                        || generation.try_get_untracked() != Some(next)
                    {
                        return;
                    }
                    if sites.len() >= 10_000 {
                        complete = false;
                        break;
                    }
                    let options = InspectionOptions {
                        max_comparisons: 50_000_000,
                        max_sites: 10_000 - sites.len(),
                        max_fragments: 1,
                        max_report_bytes: 4 * 1024 * 1024,
                    };
                    let role = format!("editor-restrictions-{batch}");
                    let binding =
                        match restriction_inspection::batch_revision(&role, &m, chunk, &options) {
                            Ok(v) => v,
                            Err(e) => {
                                failure = Some(e.to_string());
                                break;
                            }
                        };
                    let response = w
                        .ask(Request {
                            id: batch as u32 + 1,
                            revision: binding.clone(),
                            job: AnalysisJob::ScanRestriction {
                                source_role: role.clone(),
                                molecule: m.clone(),
                                enzymes: chunk.to_vec(),
                                options: options.clone(),
                            },
                        })
                        .await;
                    if !active.load(Ordering::Relaxed)
                        || generation.try_get_untracked() != Some(next)
                        || context.try_get_untracked().as_ref() != Some(&requested_context)
                        || revision.try_get_untracked().as_ref() != Some(&requested_revision)
                    {
                        return;
                    }
                    let report = match response {
                        Ok(r) if r.id == batch as u32 + 1 && r.revision == binding => {
                            match r.result {
                                AnalysisResult::RestrictionInspection(r) => r,
                                AnalysisResult::Error(e) => {
                                    failure = Some(e);
                                    break;
                                }
                                _ => {
                                    failure = Some("Unexpected restriction scan response.".into());
                                    break;
                                }
                            }
                        }
                        Ok(_) => {
                            failure = Some(
                                "Restriction scan response does not match its request.".into(),
                            );
                            break;
                        }
                        Err(e) => {
                            failure = Some(e);
                            break;
                        }
                    };
                    if let Err(e) = report.validate_bindings_for(&role, &m, chunk, &options) {
                        failure = Some(e.to_string());
                        break;
                    }
                    complete &= report.scan_complete;
                    for e in chunk {
                        if let Some(n) = counts.get_mut(&e.name) {
                            n.complete = report.scan_complete;
                        }
                    }
                    for mut site in report.sites.iter().cloned() {
                        if let Some(n) = counts.get_mut(&site.enzyme) {
                            if site.certain {
                                n.certain += 1;
                            } else {
                                n.possible += 1;
                            }
                        }
                        site.index = sites.len();
                        site.enzyme_index += batch * 64;
                        sites.push(site);
                    }
                    reports.push(report);
                }
                w.close();
                if !active.load(Ordering::Relaxed) || generation.try_get_untracked() != Some(next) {
                    return;
                }
                worker.set_value(None);
                busy.set(false);
                if let Some(error) = failure {
                    status.set(format!("Restriction scan unavailable: {error}"));
                    return;
                }
                status.set(if complete{format!("{} recognition sites across {} enzymes. Predicted cuts; methylation and reaction conditions are not assessed.",sites.len(),enzymes.len())}else{"Restriction scan incomplete. Exact-count filters exclude incomplete or ambiguous enzyme counts.".into()});
                scan.set(Some(Arc::new(Scan {
                    context: requested_context,
                    revision: requested_revision,
                    sites: Arc::new(sites),
                    counts,
                    reports: Arc::new(reports),
                })));
            });
        });
        if let Some(window) = web_sys::window() {
            match window.set_timeout_with_callback_and_timeout_and_arguments_0(
                callback.as_ref().unchecked_ref(),
                350,
            ) {
                Ok(id) => timer.set_value(Some(Timer {
                    id,
                    _callback: callback,
                })),
                Err(_) => {
                    busy.set(false);
                    status.set("Cannot schedule restriction scan.".into());
                }
            }
        }
    });
    Restrictions {
        visible,
        names,
        count_filter,
        region_filter,
        basis,
        region_status,
        scoped_counts,
        pinned,
        source: molecule,
        context,
        selection,
        selected,
        scan,
        displayed,
        status,
        busy,
        catalog,
        refresh,
    }
}

#[derive(Clone, PartialEq, Eq)]
struct Choice {
    name: String,
    recognition: String,
}

#[component]
pub fn RestrictionPanel(
    state: Restrictions,
    on_navigate: Callback<InspectionSite>,
    saved: Memo<Vec<(String, enzyme_sets::Set)>>,
    on_save: Callback<enzyme_sets::Set, Result<(), String>>,
    on_delete: Callback<String, Result<(), String>>,
) -> impl IntoView {
    let search = RwSignal::new(String::new());
    let set_name = RwSignal::new(String::new());
    let set_status = RwSignal::new(String::new());
    let site_page = RwSignal::new(0usize);
    Effect::new(move |_| {
        let _ = state.displayed.get();
        site_page.set(0);
    });
    let choices = Memo::new(move |_| {
        let needle = search.get().trim().to_ascii_lowercase();
        state
            .catalog
            .get_value()
            .into_iter()
            .filter(|e| {
                (e.name.to_ascii_lowercase().contains(&needle)
                    || e.recognition.to_ascii_lowercase().contains(&needle))
                    && state.accepts(&e.name)
            })
            .map(|e| Choice {
                name: e.name,
                recognition: e.recognition,
            })
            .collect::<Vec<_>>()
    });
    view! {<section class="inspector-content restriction-panel"><h2>"Restriction sites"</h2>
        <label class="restriction-visibility"><input type="checkbox" aria-label="Show restriction sites" prop:checked=move||state.visible.get() on:change=move|e|state.visible.set(event_target_checked(&e))/><span>"Show in sequence and map"</span></label>
        <label>"Find an enzyme"<input aria-label="Find restriction enzyme" placeholder="Name or recognition sequence" prop:value=move||search.get() on:input=move|e|search.set(event_target_value(&e))/></label>
        <label>"Region relation"<select aria-label="Restriction region relation" prop:value=move||state.region_filter.get() on:change=move|e|{
            let value=event_target_value(&e);
            let next=match value.as_str(){"inside"=>"inside","outside"=>"outside","crossing"=>"crossing","undetermined"=>"undetermined",_=>"all"};
            if state.region_filter.get_untracked()=="all" && next!="all" {state.capture_selection();}
            state.region_filter.set(next);
        }><option value="all">"Whole source"</option><option value="inside">"Inside pinned region"</option><option value="outside">"Outside pinned region"</option><option value="crossing">"Crossing region boundary"</option><option value="undetermined">"Undetermined geometry"</option></select></label>
        <label>"Classify using"<select aria-label="Restriction region geometry" prop:value=move||{match state.basis.get(){Basis::RecognitionBases=>"recognition",Basis::StrandCuts=>"cuts"}} on:change=move|e|{state.basis.set(if event_target_value(&e)=="cuts"{Basis::StrandCuts}else{Basis::RecognitionBases});}><option value="recognition">"Recognition bases"</option><option value="cuts">"Both strand-cut boundaries"</option></select></label>
        <button class="text-tool" on:click=move|_|{state.capture_selection();}>"Use current selection"</button>
        <p class="helper" role="status">{move||state.region_status.get()}</p>
        <p class="helper">"Recognition: all bases inside, none outside, overlap crossing. Cuts: both inside, both outside, mixed crossing; endpoints count as inside. Type IIS cuts can lie outside their recognition sequence. Missing or ambiguous cuts are undetermined; exact relation counts stay unavailable when such sites could change the count."</p>
        <label>"Number of sites"<select aria-label="Restriction site count" prop:value=move||state.count_filter.get() on:change=move|e|{state.count_filter.set(match event_target_value(&e).as_str(){"0"=>"0","1"=>"1","2"=>"2","many"=>"many","uncertain"=>"uncertain",_=>"all"});}><option value="all">"All counts"</option><option value="0">"0 sites"</option><option value="1">"1 site"</option><option value="2">"2 sites"</option><option value="many">"3+ sites"</option><option value="uncertain">"Uncertain / incomplete"</option></select></label>
        <div class="restriction-picker-actions"><button class="text-tool" on:click=move|_|state.names.set(choices.with_untracked(|v|v.iter().map(|e|e.name.clone()).collect()))>"Select shown"</button><button class="text-tool" on:click=move|_|state.names.set(Vec::new())>"Clear"</button><button class="text-tool" disabled=move||state.busy.get() on:click=move|_|state.refresh.update(|v|*v=v.wrapping_add(1))>"Rescan"</button></div>
        <details><summary>"Named enzyme sets"</summary>
            <label>"Set name"<input aria-label="Enzyme set name" prop:value=move||set_name.get() on:input=move|e|{set_name.set(event_target_value(&e));}/></label>
            <button class="button secondary full" on:click=move|_|{
                let result=(||->Result<(),String>{
                    let catalog=catalog_digest::builtin_catalog().map_err(|e|e.to_string())?;
                    let selected=state.names.get_untracked();
                    let available=state.catalog.get_value();
                    let choices=selected.iter().filter_map(|name|available.iter().find(|e|&e.name==name)).map(|e|restriction_inspection::EnzymeSpec{name:e.name.clone(),recognition:e.recognition.clone(),top_cut:e.top_cut,bottom_cut:e.bottom_cut}).collect::<Vec<_>>();
                    if choices.len()!=selected.len(){return Err("Some selected enzyme names are unavailable; resolve them before saving.".into());}
                    let set=enzyme_sets::build_for_catalog(&set_name.get_untracked(),&choices,&catalog).map_err(|e|e.to_string())?;
                    on_save.run(set)
                })();
                set_status.set(match result{Ok(())=>"Saved the checked enzymes with their exact geometry. An existing set with this name is replaced.".into(),Err(e)=>e});
            }>"Save / replace checked enzymes"</button>
            {move||saved.get().into_iter().map(|(key,set)|{
                let delete_key=key.clone();
                let label=format!("{} · {} enzymes",set.name,set.choices.len());
                let name=set.name.clone();
                view!{<div class="restriction-picker-actions"><span>{label}</span>
                    <button class="text-tool" on:click=move|_|{
                        let result=catalog_digest::builtin_catalog().map_err(|e|e.to_string()).and_then(|catalog|enzyme_sets::resolve_for_catalog(&set,&catalog).map_err(|e|e.to_string()));
                        match result{
                            Ok(resolved)=>{
                                if !resolved.missing.is_empty() || !resolved.changed.is_empty(){
                                    set_status.set(format!("Set not loaded. Missing: {}. Changed geometry: {}. Your checked enzymes are unchanged.",if resolved.missing.is_empty(){"none".into()}else{resolved.missing.join(", ")},if resolved.changed.is_empty(){"none".into()}else{resolved.changed.join(", ")}));
                                }else{
                                    state.names.set(resolved.available_names);
                                    set_name.set(name.clone());
                                    set_status.set(if resolved.catalog_changed{"Loaded exact matching enzymes; the catalog edition has changed.".into()}else{"Loaded saved enzyme set.".into()});
                                }
                            },Err(e)=>set_status.set(e)
                        }
                    }>"Load"</button>
                    <button class="text-tool" on:click=move|_|{set_status.set(match on_delete.run(delete_key.clone()){Ok(())=>"Deleted enzyme set.".into(),Err(e)=>e});}>"Delete"</button>
                </div>}
            }).collect_view()}
            <p class="helper" role="status">{move||set_status.get()}</p>
        </details>
        <p class="helper restriction-status" role="status">{move||state.status.get()}</p>
        <div class="restriction-enzyme-list">{move||choices.get().into_iter().map(|e|{let name=e.name.clone();let change=e.name.clone();let count=e.name.clone();view!{<label class="restriction-enzyme"><input type="checkbox" aria-label=format!("Show enzyme {}",e.name) prop:checked=move||state.names.with(|v|v.contains(&name)) on:change=move|event|{let checked=event_target_checked(&event);state.names.update(|v|{v.retain(|n|n!=&change);if checked{v.push(change.clone());}});}/><span><strong>{e.name}</strong><small>{e.recognition}</small></span><span class="restriction-count">{move||state.count(&count).label()}</span></label>}}).collect_view()}</div>
        <Show when=move||choices.with(Vec::is_empty)><p class="helper">"No enzymes match this filter. Incomplete scans do not establish zero sites."</p></Show>
        {move||state.selected.get().and_then(|index|state.scan.with(|s|s.as_ref().and_then(|s|s.sites.iter().find(|site|site.index==index).cloned()))).map(|site|view!{<article class="restriction-site-details"><h3>{site.enzyme.clone()}</h3><p>{format!("Recognition: {}{}",site.spans.iter().map(|p|format!("{}–{}",p.start+1,p.end)).collect::<Vec<_>>().join(", "),if site.certain{""}else{" · possible match"})}</p><p>{format!("Top cut: {} · bottom cut: {}",cut_label(site.top_cut,site.top_boundary),cut_label(site.bottom_cut,site.bottom_boundary))}</p><button class="text-tool" on:click=move|_|on_navigate.run(site.clone())>"Select recognition sequence"</button></article>})}
        <button class="button secondary full" disabled=move||state.scan.with(|s|s.is_none()) on:click=move|_|{
            if let Some(scan)=state.scan.get_untracked(){
                let evidence=serde_json::json!({"schema_version":1,"kind":"restriction-recognition-scan","source_context":scan.context,"source_revision":scan.revision,"counts":scan.counts,"reports":scan.reports.as_ref(),"interpretation":"Recognition/cut geometry evidence only. Not a digest product or an experimental cleavage claim."});
                match serde_json::to_vec_pretty(&evidence){Ok(bytes)=>if let Err(e)=crate::app::download("restriction-sites.json","application/json",&bytes){state.status.set(e);},Err(e)=>state.status.set(e.to_string())}
            }
        }>"Download RE evidence"</button>
        <details class="restriction-site-list"><summary>"Shown recognition sites"</summary>{move||state.displayed.get().map(|sites|sites.iter().skip(site_page.get()*100).take(100).cloned().map(|site|{let label=format!("{} · {}{}",site.enzyme,site.spans.iter().map(|p|format!("{}–{}",p.start+1,p.end)).collect::<Vec<_>>().join(", "),if site.certain{""}else{" · possible"});view!{<button class="result-entry" on:click=move|_|{state.selected.set(Some(site.index));on_navigate.run(site.clone());}>{label}</button>}}).collect_view())}<div class="restriction-picker-actions"><button class="text-tool" disabled=move||site_page.get()==0 on:click=move|_|site_page.update(|p|*p=p.saturating_sub(1))>"Previous sites"</button><span>{move||state.displayed.with(|sites|{let total=sites.as_ref().map_or(0,|s|s.len());if total==0{"0 sites".into()}else{format!("{}–{} of {}",site_page.get()*100+1,((site_page.get()+1)*100).min(total),total)}})}</span><button class="text-tool" disabled=move||{state.displayed.with(|sites|sites.as_ref().is_none_or(|s|(site_page.get()+1)*100>=s.len()))} on:click=move|_|site_page.update(|p|*p+=1)>"Next sites"</button></div></details>
    </section>}.into_any()
}
fn cut_label(value: Option<usize>, raw: i64) -> String {
    value
        .map(|n| format!("boundary {n}"))
        .unwrap_or_else(|| format!("outside substrate ({raw})"))
}
