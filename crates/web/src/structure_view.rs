use crate::{structure_camera::Camera, structure_renderer::Renderer};
use collagene_core::{
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    structure_view::{
        self as native,
        coordinates::{
            Format,
            scene::{Alternates, Chains, Selection},
        },
    },
};
use leptos::{html, prelude::*};
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};
use wasm_bindgen_futures::{JsFuture, spawn_local};
#[derive(Clone)]
pub struct SourceFile {
    pub name: String,
    pub bytes: Vec<u8>,
    pub format: Format,
}
fn identity_label(id: &native::coordinates::Identity) -> String {
    format!(
        "{} {} · chain {} · residue {}",
        id.component.as_deref().unwrap_or("undeclared"),
        id.atom.as_deref().unwrap_or("undeclared"),
        id.chain.as_deref().unwrap_or("blank / undeclared"),
        id.sequence.as_deref().unwrap_or("undeclared")
    )
}
fn optional_number(v: Option<f64>) -> String {
    v.map(|v| v.to_string())
        .unwrap_or_else(|| "undeclared".into())
}
fn example_ids<'a>(ids: impl Iterator<Item = &'a str>) -> String {
    let mut names = std::collections::BTreeSet::new();
    for id in ids {
        names.insert(id);
        if names.len() == 16 {
            break;
        }
    }
    if names.is_empty() {
        "none declared".into()
    } else {
        names.into_iter().collect::<Vec<_>>().join(", ")
    }
}
#[component]
pub fn StructurePanel(
    entries: Memo<Vec<(String, String, Format)>>,
    get_source: Callback<String, Result<SourceFile, String>>,
    retain: Callback<(SourceFile, u64), Result<String, String>>,
    remove: Callback<(String, u64), Result<(), String>>,
    epoch: RwSignal<u64>,
    locked: Memo<bool>,
) -> impl IntoView {
    let canvas = NodeRef::<html::Canvas>::new();
    let source = RwSignal::new(None::<Arc<SourceFile>>);
    let report = RwSignal::new(None::<Arc<native::Prepared>>);
    let catalog = RwSignal::new(None::<Arc<native::Prepared>>);
    let selection = RwSignal::new(None::<Selection>);
    let camera = RwSignal::new(Camera::default());
    let picked = RwSignal::new(None::<u32>);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let alive = StoredValue::new(Arc::new(AtomicBool::new(true)));
    let job = StoredValue::new(None::<Arc<AtomicBool>>);
    let renderer = StoredValue::new_local(None::<Renderer>);
    let resize = StoredValue::new_local(None::<crate::renderer::ResizeListener>);
    let pointer = RwSignal::new(None::<(i32, f64, f64, f64, f64)>);
    let dragged = RwSignal::new(false);
    let saved = RwSignal::new(String::new());
    let format_choice = RwSignal::new("auto".to_string());
    let cancel = Callback::new(move |()| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(a) = job.get_value() {
            a.store(false, Ordering::Relaxed);
        }
        job.set_value(None);
        busy.set(false);
    });
    Effect::new(move |_| {
        let _ = epoch.get();
        cancel.run(());
        source.set(None);
        report.set(None);
        catalog.set(None);
        selection.set(None);
        picked.set(None);
        status.set(String::new());
    });
    on_cleanup(move || {
        alive.get_value().store(false, Ordering::Relaxed);
        if let Some(a) = job.get_value() {
            a.store(false, Ordering::Relaxed);
        }
        resize.set_value(None);
        renderer.set_value(None);
    });
    let run = Callback::new(
        move |(asset, options, keep): (Arc<SourceFile>, Option<Selection>, bool)| {
            if locked.get_untracked() {
                status.set("Wait for workspace restoration.".into());
                return;
            }
            cancel.run(());
            let token = generation.get_untracked();
            let captured_epoch = epoch.get_untracked();
            let revision = match native::revision(&asset.bytes, asset.format, &options) {
                Ok(v) => v,
                Err(e) => {
                    status.set(e.to_string());
                    return;
                }
            };
            source.set(Some(asset.clone()));
            report.set(None);
            picked.set(None);
            busy.set(true);
            status.set("Preparing imported coordinates…".into());
            let flag = Arc::new(AtomicBool::new(true));
            job.set_value(Some(flag.clone()));
            let owner = alive.get_value();
            spawn_local(async move {
                let response = crate::worker::analyze_cancellable(
                    Request {
                        id: token,
                        revision: revision.clone(),
                        job: AnalysisJob::PrepareStructure {
                            bytes: asset.bytes.clone(),
                            format: asset.format,
                            selection: options,
                        },
                    },
                    flag,
                )
                .await;
                if !owner.load(Ordering::Relaxed)
                    || epoch.get_untracked() != captured_epoch
                    || generation.get_untracked() != token
                {
                    return;
                }
                busy.set(false);
                job.set_value(None);
                let result = response.and_then(|r| {
                    if r.id != token || r.revision != revision {
                        Err("Structure response context changed.".into())
                    } else {
                        match r.result {
                            AnalysisResult::StructurePrepared(p) => Ok(p),
                            AnalysisResult::Error(e) => Err(e),
                            _ => Err("Unexpected structure response.".into()),
                        }
                    }
                });
                match result {
                    Ok(value) => {
                        if locked.get_untracked() {
                            status.set(
                                "Workspace changed while importing; reopen the source.".into(),
                            );
                            return;
                        }
                        let mut message="Imported coordinate view. No sequence correspondence or experimental validity is inferred.".to_string();
                        if keep {
                            match retain.run(((*asset).clone(),captured_epoch)){Ok(_)=>message="Original structure retained in the workspace. No sequence link is inferred.".into(),Err(e)=>message=format!("View prepared, but original source was not saved: {e}")}
                        }
                        selection.set(Some(value.scene.selection.clone()));
                        let value = Arc::new(value);
                        catalog.set(Some(value.clone()));
                        report.set(Some(value));
                        camera.set(Camera::default());
                        status.set(message);
                    }
                    Err(e) => status.set(e),
                }
            });
        },
    );
    Effect::new(move |_| {
        let value = report.get();
        let Some(element) = canvas.get() else {
            return;
        };
        resize.set_value(None);
        renderer.set_value(None);
        if let Some(value) = value {
            match Renderer::new(&element, &value.scene) {
                Ok(r) => {
                    if let Err(e) = r.draw(&element, camera.get_untracked(), picked.get_untracked())
                    {
                        status.set(e);
                    }
                    renderer.set_value(Some(r));
                }
                Err(e) => status.set(e),
            }
            if let Ok(observer) = crate::renderer::ResizeListener::new(&element, move || {
                if let Some(c) = canvas.get_untracked() {
                    renderer.with_value(|r| {
                        if let Some(r) = r {
                            if let Err(e) =
                                r.draw(&c, camera.get_untracked(), picked.get_untracked())
                            {
                                status.set(e);
                            }
                        }
                    });
                }
            }) {
                resize.set_value(Some(observer));
            }
        } else {
            element.set_width(1);
            element.set_height(1);
        }
    });
    Effect::new(move |_| {
        let c = camera.get();
        let _ = picked.get();
        if let Some(canvas) = canvas.get() {
            renderer.with_value(|r| {
                if let Some(r) = r {
                    if let Err(e) = r.draw(&canvas, c, picked.get_untracked()) {
                        status.set(e);
                    }
                }
            });
        }
    });
    let ingest = move |event: web_sys::Event| {
        let input = event_target::<web_sys::HtmlInputElement>(&event);
        let Some(file) = input.files().and_then(|f| f.get(0)) else {
            return;
        };
        input.set_value("");
        if locked.get_untracked() {
            return;
        }
        cancel.run(());
        catalog.set(None);
        report.set(None);
        source.set(None);
        if file.size() <= 0.0 || file.size() > native::coordinates::MAX_SOURCE_BYTES as f64 {
            status.set("Choose a nonempty coordinate file no larger than 16 MiB.".into());
            return;
        }
        let name = file.name();
        let format = match format_choice.get_untracked().as_str() {
            "pdb" => Format::Pdb,
            "mmcif" => Format::Mmcif,
            _ => {
                let lower = name.to_ascii_lowercase();
                if lower.ends_with(".pdb") {
                    Format::Pdb
                } else if lower.ends_with(".cif") || lower.ends_with(".mmcif") {
                    Format::Mmcif
                } else {
                    status.set("Choose PDB or mmCIF explicitly for this filename.".into());
                    return;
                }
            }
        };
        if name.is_empty() || name.len() > 256 || name.chars().any(char::is_control) {
            status.set("Coordinate filename is invalid or too long.".into());
            return;
        }
        let token = generation.get_untracked();
        let captured_epoch = epoch.get_untracked();
        let owner = alive.get_value();
        busy.set(true);
        status.set("Reading coordinates…".into());
        spawn_local(async move {
            let data = JsFuture::from(file.array_buffer()).await;
            if !owner.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || epoch.get_untracked() != captured_epoch
            {
                return;
            }
            busy.set(false);
            match data {
                Ok(data) => run.run((
                    Arc::new(SourceFile {
                        name,
                        bytes: js_sys::Uint8Array::new(&data).to_vec(),
                        format,
                    }),
                    None,
                    true,
                )),
                Err(_) => status.set("Could not read coordinate file.".into()),
            }
        });
    };
    let cycle_atom = Callback::new(move |direction: i32| {
        report.with_untracked(|r| {
            if let Some(r) = r {
                let points = &r.scene.points;
                if points.is_empty() {
                    return;
                }
                let current = picked
                    .get_untracked()
                    .and_then(|i| points.iter().position(|p| p.atom_index == i));
                let next = match current {
                    Some(i) => {
                        (i as i64 + direction as i64).rem_euclid(points.len() as i64) as usize
                    }
                    None => {
                        if direction < 0 {
                            points.len() - 1
                        } else {
                            0
                        }
                    }
                };
                picked.set(Some(points[next].atom_index));
            }
        });
    });
    let download = move |_: web_sys::MouseEvent| {
        if let Some(asset) = source.get_untracked() {
            let captured_epoch = epoch.get_untracked();
            let hash = blake3::hash(&asset.bytes);
            let owner = alive.get_value();
            let guard_owner = owner.clone();
            spawn_local(async move {
                let result = crate::app::download_guarded(
                    &asset.name,
                    if asset.format == Format::Pdb {
                        "chemical/x-pdb"
                    } else {
                        "chemical/x-mmcif"
                    },
                    &asset.bytes,
                    move || {
                        if !guard_owner.load(Ordering::Relaxed)
                            || epoch.get_untracked() != captured_epoch
                            || !source.with_untracked(|s| {
                                s.as_ref().is_some_and(|s| blake3::hash(&s.bytes) == hash)
                            })
                        {
                            Err("Structure source changed before export.".into())
                        } else {
                            Ok(())
                        }
                    },
                )
                .await;
                if owner.load(Ordering::Relaxed) {
                    if let Err(e) = result {
                        status.set(e);
                    }
                }
            });
        }
    };
    view! {<section class="structure-panel"><p class="helper">"Imported coordinates only. Drag to rotate, Shift-drag to pan, scroll to zoom. Select an atom for its original coordinates."</p>
      <details class="structure-source-controls" open><summary>"Structure sources"</summary>
        <label>"Input profile"<select prop:value=move||format_choice.get() on:change=move|e|{cancel.run(());format_choice.set(event_target_value(&e));}><option value="auto">"Use filename hint, validate contents"</option><option value="pdb">"PDB"</option><option value="mmcif">"mmCIF"</option></select></label>
        <label class="structure-file">"Import coordinate file"<input type="file" accept=".pdb,.cif,.mmcif" disabled=move||locked.get() on:change=ingest/></label>
        <label>"Workspace structures"<select prop:value=move||saved.get() on:change=move|e|saved.set(event_target_value(&e))><option value="">"Choose retained source"</option>{move||entries.get().into_iter().map(|(hash,name,_)|view!{<option value=hash>{name}</option>}).collect_view()}</select></label>
        <div class="structure-actions"><button class="button secondary" disabled=move||locked.get()||saved.get().is_empty() on:click=move|_|{match get_source.run(saved.get_untracked()){Ok(a)=>{catalog.set(None);run.run((Arc::new(a),None,false));},Err(e)=>status.set(e)}}>"Open"</button>
        <button class="button secondary" disabled=move||locked.get()||saved.get().is_empty() on:click=move|_|{let hash=saved.get_untracked();match remove.run((hash.clone(),epoch.get_untracked())){Ok(())=>{if source.with_untracked(|s|s.as_ref().is_some_and(|s|blake3::hash(&s.bytes).to_hex().as_str()==hash)){cancel.run(());source.set(None);report.set(None);catalog.set(None);}saved.set(String::new());status.set("Removed retained structure source.".into());},Err(e)=>status.set(e)}}>"Remove"</button></div>
      </details>
      <p class="helper" role="status">{move||status.get()}</p>
      {move||source.get().map(|s|view!{<p class="helper"><strong>"Current source: "</strong>{s.name.clone()}</p>})}
      <Show when=move||busy.get()><button class="button secondary" on:click=move|_|{cancel.run(());status.set("Structure preparation cancelled.".into());}>"Cancel"</button></Show>
      <canvas class="structure-canvas" node_ref=canvas tabindex="0" aria-label="Imported 3D coordinates. Arrow keys rotate, plus and minus zoom, Home resets."
        on:pointerdown=move|e|{dragged.set(false);if let Some(c)=canvas.get(){let _=c.set_pointer_capture(e.pointer_id());}pointer.set(Some((e.pointer_id(),e.client_x() as f64,e.client_y() as f64,e.client_x() as f64,e.client_y() as f64)));}
        on:pointermove=move|e|{if let Some((id,x,y,start_x,start_y))=pointer.get_untracked(){if id!=e.pointer_id(){return;}if let Some(c)=canvas.get(){let rect=c.get_bounding_client_rect();if (e.client_x()as f64-start_x).hypot(e.client_y()as f64-start_y)>=4.0{dragged.set(true);}camera.update(|cam|cam.drag(e.client_x()as f64-x,e.client_y()as f64-y,e.shift_key(),rect.width(),rect.height()));pointer.set(Some((id,e.client_x()as f64,e.client_y()as f64,start_x,start_y)));}}}
        on:pointerup=move|e|{if let Some((id,_,_,x,y))=pointer.get_untracked(){if id!=e.pointer_id(){return;}pointer.set(None);if !dragged.get_untracked() && (e.client_x()as f64-x).hypot(e.client_y()as f64-y)<4.0 {if let (Some(c),Some(r))=(canvas.get(),report.get_untracked()){let rect=c.get_bounding_client_rect();picked.set(camera.get_untracked().pick(&r.scene,e.client_x()as f64-rect.left(),e.client_y()as f64-rect.top(),rect.width(),rect.height()));}}}}
        on:pointercancel=move|_|pointer.set(None) on:lostpointercapture=move|_|pointer.set(None)
        on:wheel=move|e|{e.prevent_default();camera.update(|c|c.scroll(e.delta_y()));}
        on:keydown=move|e|{let key=e.key();if matches!(key.as_str(),"ArrowLeft"|"ArrowRight"|"ArrowUp"|"ArrowDown"|"+"|"="|"-"|"Home"){e.prevent_default();camera.update(|c|match key.as_str(){"ArrowLeft"=>c.drag(-8.0,0.0,false,400.0,400.0),"ArrowRight"=>c.drag(8.0,0.0,false,400.0,400.0),"ArrowUp"=>c.drag(0.0,-8.0,false,400.0,400.0),"ArrowDown"=>c.drag(0.0,8.0,false,400.0,400.0),"+"|"="=>c.scroll(-100.0),"-"=>c.scroll(100.0),"Home"=>*c=Camera::default(),_=>{}});}}
      ></canvas>
      <div class="structure-actions"><button class="button secondary" on:click=move|_|camera.set(Camera::default())>"Reset view"</button><button class="button secondary" disabled=move||source.get().is_none() on:click=download>"Download original"</button></div>
      {move||catalog.get().map(|data|{let models=data.atoms.iter().map(|a|a.model).collect::<std::collections::BTreeSet<_>>();let chain_examples=example_ids(data.atoms.iter().filter_map(|a|a.label.chain.as_deref()));let author_examples=example_ids(data.atoms.iter().map(|a|a.author.chain.as_deref().unwrap_or("(blank / undeclared)")));let alternate_examples=example_ids(data.atoms.iter().filter_map(|a|a.alternate.as_deref()));view!{<details><summary>"Display filters"</summary>
        <label>"Model"<select prop:value=move||selection.with(|s|s.as_ref().map(|s|s.model.to_string()).unwrap_or_default()) on:change=move|e|if let Ok(n)=event_target_value(&e).parse(){selection.update(|s|if let Some(s)=s{s.model=n;});}>{models.into_iter().map(|n|view!{<option value=n.to_string()>{n}</option>}).collect_view()}</select></label>
        <label>"Chain namespace"<select prop:value=move||selection.with(|s|s.as_ref().map(|s|match &s.chains{Chains::All=>"all",Chains::Label(_)=>"label",Chains::Author(_)=>"author"}).unwrap_or("all")) on:change=move|e|{let v=event_target_value(&e);selection.update(|s|if let Some(s)=s{s.chains=match v.as_str(){"label"=>Chains::Label(String::new()),"author"=>Chains::Author(None),_=>Chains::All};});}><option value="all">"All chains"</option><option value="label">"Exact label ID"</option><option value="author">"Exact author ID"</option></select></label>
        <label>"Case-sensitive chain ID"<input maxlength="256" disabled=move||selection.with(|s|s.as_ref().is_none_or(|s|matches!(s.chains,Chains::All))) prop:value=move||selection.with(|s|s.as_ref().and_then(|s|match &s.chains{Chains::Label(v)=>Some(v.clone()),Chains::Author(v)=>v.clone(),_=>None}).unwrap_or_default()) on:input=move|e|{let v=event_target_value(&e);selection.update(|s|if let Some(s)=s{match &s.chains{Chains::Label(_)=>s.chains=Chains::Label(v),Chains::Author(_)=>s.chains=Chains::Author((!v.is_empty()).then_some(v)),_=>{}}});}/></label>
        <p class="helper">{format!("Label examples: {chain_examples}. Author examples: {author_examples}. Blank author ID selects blank/undeclared author chains.")}</p>
        <label>"Alternate locations"<select prop:value=move||selection.with(|s|s.as_ref().map(|s|if matches!(s.alternates,Alternates::All){"all"}else{"named"}).unwrap_or("all")) on:change=move|e|{let v=event_target_value(&e);selection.update(|s|if let Some(s)=s{s.alternates=if v=="all"{Alternates::All}else{Alternates::Named(String::new())};});}><option value="all">"All deposited alternates"</option><option value="named">"Exact alternate ID plus unlabelled atoms"</option></select></label>
        <label>"Case-sensitive alternate ID"<input maxlength="256" disabled=move||selection.with(|s|s.as_ref().is_none_or(|s|matches!(s.alternates,Alternates::All))) prop:value=move||selection.with(|s|s.as_ref().and_then(|s|if let Alternates::Named(n)=&s.alternates{Some(n.clone())}else{None}).unwrap_or_default()) on:input=move|e|{let v=event_target_value(&e);selection.update(|s|if let Some(s)=s{s.alternates=Alternates::Named(v);});}/></label>
        <p class="helper">{format!("Alternate examples: {alternate_examples}. Examples show at most 16 distinct identifiers; any exact deposited ID can be entered.")}</p>
        <label class="checkbox-label"><input type="checkbox" prop:checked=move||selection.with(|s|s.as_ref().is_some_and(|s|s.hydrogens)) on:change=move|e|selection.update(|s|if let Some(s)=s{s.hydrogens=event_target_checked(&e);})/>"Show deposited H/D atoms"</label>
        <button class="button secondary" disabled=move||busy.get()||locked.get() on:click=move|_|{if let Some(a)=source.get_untracked(){run.run((a,selection.get_untracked(),false));}}>"Apply display filters"</button><p class="helper">"Filters apply when requested. Alternate labels are display filters, not evidence of a complete conformer."</p>
      </details>}})}
      {move||report.get().map(|r|{let s=&r.scene;view!{<p class="helper">{format!("Displayed model {} · {} / {} atoms · {} declared edges · {} unresolved · {} relations/symmetry withheld",s.selection.model,s.points.len(),r.atoms.len(),s.edges.len(),s.unresolved_connections,s.withheld_symmetry_or_relation)}</p><details><summary>"Source provenance"</summary><p class="structure-hash">{format!("{:?} · source hash {}",r.format,r.source_hash)}</p><p class="helper">{format!("{} atoms lack an element; {} have alternate labels. Original source retains unprojected records/categories.",r.notes.atoms_without_element,r.notes.atoms_with_alternates)}</p></details>}})}
      <details><summary>"Atom inspector"</summary><div class="structure-actions"><button class="button secondary" disabled=move||report.get().is_none() on:click=move|_|cycle_atom.run(-1)>"Previous atom"</button><button class="button secondary" disabled=move||report.get().is_none() on:click=move|_|cycle_atom.run(1)>"Next atom"</button></div>
      {move||picked.get().and_then(|i|report.with(|r|r.as_ref().and_then(|r|r.atoms.get(i as usize).cloned()))).map(|a|view!{<div class="structure-atom"><strong>{format!("Atom {} · model {} · {}",a.id,a.model,a.element.as_deref().unwrap_or("unknown element"))}</strong><p>{format!("Label: {} · Author: {}",identity_label(&a.label),identity_label(&a.author))}</p><p>{format!("x {} · y {} · z {} Å",a.position[0],a.position[1],a.position[2])}</p><p>{format!("Alternate {} · insertion {} · occupancy {} · B {} Å² · charge {}",a.alternate.as_deref().unwrap_or("none / undeclared"),a.insertion.as_deref().unwrap_or("none / undeclared"),optional_number(a.occupancy),optional_number(a.b_factor),a.formal_charge.map(|v|v.to_string()).unwrap_or_else(||"undeclared".into()))}</p></div>})}
      </details>
    </section>}
}
