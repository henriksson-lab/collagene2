//! Explicit target-free material review; installation uses the shared native archive transaction.
use crate::history::run_origin_model::PinnedOrigin;
use crate::platform::worker::CachedWorker;
use crate::restriction::material_derivation_model::MismatchChoice;
use crate::restriction::material_derivation_model::OverlapExtensionChoice;
use crate::restriction::material_derivation_model::RepairChoice;
use crate::restriction::material_derivation_model::SourcePin;
use crate::restriction::material_derivation_model::TerminalChoice;
use crate::sources::identity::SourceInfo;
use collagene_core::{
    chemistry::EndRepairMethod,
    cloning::pcr::mutagenesis::PcrOptions,
    documents::workspace::Workspace,
    material::derivation as native,
    protocol::{
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
};
use collagene_model::Molecule;
use leptos::prelude::*;
use serde::Serialize;
use std::{
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    MismatchPcr,
    EndRepair,
    OverlapExtensionPcr,
    TerminalProcessing,
}
impl Kind {
    fn title(self) -> &'static str {
        match self {
            Self::MismatchPcr => "Open reusable mismatch-PCR material",
            Self::EndRepair => "Repair saved material directly",
            Self::OverlapExtensionPcr => "Reconstruct reusable overlap-extension PCR material",
            Self::TerminalProcessing => "Treat or tail source material directly",
        }
    }
    fn product_name(self) -> &'static str {
        match self {
            Self::MismatchPcr => "Reusable mismatch-PCR product",
            Self::EndRepair => "Repaired material",
            Self::OverlapExtensionPcr => "Reusable overlap-extension PCR product",
            Self::TerminalProcessing => "Terminally processed material",
        }
    }
    fn preview_label(self) -> &'static str {
        match self {
            Self::MismatchPcr => "Preview reusable mismatch PCR",
            Self::EndRepair => "Preview direct end repair",
            Self::OverlapExtensionPcr => "Preview reusable overlap-extension PCR",
            Self::TerminalProcessing => "Preview terminal processing",
        }
    }
    fn evidence_name(self) -> &'static str {
        match self {
            Self::MismatchPcr => "material-mismatch-pcr.json",
            Self::EndRepair => "material-end-repair.json",
            Self::OverlapExtensionPcr => "material-overlap-extension-pcr.json",
            Self::TerminalProcessing => "material-terminal-processing.json",
        }
    }
}
#[derive(Clone, Serialize)]
pub enum Settings {
    Mismatch {
        forward: Primer,
        reverse: Primer,
        options: PcrOptions,
    },
    Repair {
        method: EndRepairMethod,
    },
    OverlapExtension {
        preparation: collagene_core::cloning::overlap::extension::pcr::Request,
    },
    TerminalProcessing {
        action: native::TerminalAction,
    },
}
enum Preview {
    Mismatch(native::pcr::MismatchPreview),
    Repair(native::end_repair::RepairPreview),
    OverlapExtension(native::overlap_extension::OverlapExtensionPreview),
    TerminalProcessing(native::terminal::TerminalProcessingPreview),
}
struct Reviewed {
    source: SourcePin,
    settings: Settings,
    settings_revision: String,
    preview: Preview,
}
fn snippet(text: &str) -> String {
    let n = 240.min(text.len());
    if text.len() > n {
        format!("{}… ({} bases)", &text[..n], text.len())
    } else {
        text.into()
    }
}
#[component]
pub fn ReusableMaterialView(
    source: Memo<Option<Molecule>>,
    source_info: Memo<Option<SourceInfo>>,
    settings: Callback<(), Result<Settings, String>>,
    settings_revision: Memo<String>,
    kind: Kind,
    origin: RwSignal<Option<PinnedOrigin>>,
    epoch: RwSignal<u64>,
    get_workspace: Callback<(), Result<Workspace, String>>,
    on_open: Callback<crate::restriction::digest_product::MaterialOpenRequest, Result<(), String>>,
) -> impl IntoView {
    let reviewed = RwSignal::new(None::<Arc<Reviewed>>);
    let selected = RwSignal::new(String::new());
    let name = RwSignal::new(kind.product_name().to_string());
    let save_primers = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let busy = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
    });
    Effect::new(move || {
        let _ = (
            source_info.get(),
            source.with(|m| m.as_ref().map(Molecule::revision)),
            settings_revision.get(),
            origin.get(),
            epoch.get(),
        );
        generation.update(|v| *v = v.wrapping_add(1));
        reviewed.set(None);
        selected.set(String::new());
        busy.set(false);
        status.set(String::new());
        if let Some(worker) = worker.get_value() {
            worker.close();
        }
        worker.set_value(None);
        if origin.with_untracked(|p| {
            p.as_ref().is_some_and(|p| {
                source_info
                    .get_untracked()
                    .as_ref()
                    .is_none_or(|s| p.event() != &s.binding)
            })
        }) {
            origin.set(None);
        }
    });
    let select_origin = Callback::new(move |next: Option<PinnedOrigin>| -> Result<(), String> {
        if let Some(next) = &next {
            let info = source_info
                .get_untracked()
                .ok_or("Choose the exact source document.")?;
            next.as_ref()
                .validate_for_workspace(&get_workspace.run(())?, &info.binding)
                .map_err(|e| e.to_string())?;
        }
        origin.set(next);
        Ok(())
    });
    let guard_alive = alive.clone();
    let guard = Callback::new(
        move |(pin, key): (SourcePin, String)| -> Result<(), String> {
            if !guard_alive.load(Ordering::Relaxed) {
                return Err("The material review was closed.".into());
            }
            if settings_revision.get_untracked() != key {
                return Err("The operation settings changed; review again.".into());
            }
            let info = source_info
                .get_untracked()
                .ok_or("The source document is unavailable.")?;
            let molecule = source
                .get_untracked()
                .ok_or("The source revision is unavailable.")?;
            pin.validate_current(
                &info.binding,
                epoch.get_untracked(),
                &molecule,
                origin.get_untracked().as_ref().map(PinnedOrigin::as_ref),
            )
        },
    );
    let run_alive = alive.clone();
    let run = move |_| {
        let prepared = (|| -> Result<_, String> {
            let molecule = source
                .get_untracked()
                .ok_or("Choose an available source.")?;
            let info = source_info
                .get_untracked()
                .ok_or("Choose an exact source event.")?;
            let frame = get_workspace.run(())?;
            if kind == Kind::EndRepair && origin.get_untracked().is_none() {
                return Err(
                    "Select the saved Run or Derived material origin before direct end repair."
                        .into(),
                );
            }
            if !frame.documents.iter().any(|doc| {
                doc.id == info.binding.document_id
                    && doc.head == info.binding.revision
                    && doc.head_event.as_ref() == Some(&info.binding.event_id)
            }) {
                return Err("The selected source head changed; choose it again.".into());
            }
            let pin = SourcePin::capture(
                info.binding,
                epoch.get_untracked(),
                &molecule,
                origin.get_untracked().map(|p| Arc::new(p.as_ref().clone())),
            )?;
            pin.attest(&frame)?;
            let configured = settings.run(())?;
            let (revision, job) = match &configured {
                Settings::Mismatch {
                    forward,
                    reverse,
                    options,
                } => (
                    native::pcr::mismatch_pcr_preview_revision(
                        &molecule,
                        &pin.native_origin(),
                        forward,
                        reverse,
                        options,
                    )
                    .map_err(|e| e.to_string())?,
                    AnalysisJob::InspectMaterialMismatchPcr {
                        source: molecule,
                        origin: pin.native_origin(),
                        forward: forward.clone(),
                        reverse: reverse.clone(),
                        options: options.clone(),
                    },
                ),
                Settings::OverlapExtension { preparation } => (
                    native::overlap_extension::overlap_extension_preview_revision(
                        &frame,
                        pin.event(),
                        &pin.native_origin(),
                        preparation,
                    )
                    .map_err(|e| e.to_string())?,
                    AnalysisJob::InspectMaterialOverlapExtension {
                        workspace: frame,
                        source: pin.event().clone(),
                        origin: pin.native_origin(),
                        preparation: preparation.clone(),
                    },
                ),
                Settings::TerminalProcessing { action } => (
                    native::terminal::terminal_processing_preview_revision(
                        &frame,
                        pin.event(),
                        &pin.native_origin(),
                        action,
                    )
                    .map_err(|e| e.to_string())?,
                    AnalysisJob::InspectMaterialTerminalProcessing {
                        workspace: frame,
                        source: pin.event().clone(),
                        origin: pin.native_origin(),
                        action: *action,
                    },
                ),
                Settings::Repair { method } => (
                    native::end_repair::repair_preview_revision(
                        &molecule,
                        &pin.native_origin(),
                        *method,
                    )
                    .map_err(|e| e.to_string())?,
                    AnalysisJob::InspectMaterialRepair {
                        source: molecule,
                        origin: pin.native_origin(),
                        method: *method,
                    },
                ),
            };
            Ok((pin, configured, revision, job))
        })();
        let (pin, configured, revision, job) = match prepared {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let key = settings_revision.get_untracked();
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        reviewed.set(None);
        busy.set(true);
        status.set("Reconstructing the selected source and operation…".into());
        if let Some(old) = worker.get_value() {
            old.close();
        }
        let active = match CachedWorker::new() {
            Ok(w) => Rc::new(w),
            Err(e) => {
                busy.set(false);
                status.set(e);
                return;
            }
        };
        worker.set_value(Some(active.clone()));
        let alive = run_alive.clone();
        spawn_local(async move {
            let response = active
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || guard.run((pin.clone(), key.clone())).is_err()
            {
                return;
            }
            busy.set(false);
            let result = (|| -> Result<Preview, String> {
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("The worker reply names another operation context.".into());
                }
                match (kind, response.result) {
                    (Kind::MismatchPcr, AnalysisResult::MaterialMismatchPcrPreview(p)) => {
                        p.validate_binding(&revision).map_err(|e| e.to_string())?;
                        if p.source_revision != pin.event().revision {
                            return Err("Mismatch-PCR preview names another source.".into());
                        }
                        Ok(Preview::Mismatch(p))
                    }
                    (Kind::EndRepair, AnalysisResult::MaterialRepairPreview(p)) => {
                        p.validate_binding(&revision).map_err(|e| e.to_string())?;
                        if p.source_revision != pin.event().revision {
                            return Err("Repair preview names another source.".into());
                        }
                        Ok(Preview::Repair(p))
                    }
                    (
                        Kind::OverlapExtensionPcr,
                        AnalysisResult::MaterialOverlapExtensionPreview(p),
                    ) => {
                        p.validate_binding(&revision).map_err(|e| e.to_string())?;
                        if &p.source != pin.event() || p.source_revision != pin.event().revision {
                            return Err(
                                "Overlap-extension PCR preview names another exact source event."
                                    .into(),
                            );
                        }
                        Ok(Preview::OverlapExtension(p))
                    }
                    (
                        Kind::TerminalProcessing,
                        AnalysisResult::MaterialTerminalProcessingPreview(p),
                    ) => {
                        p.validate_binding(&revision).map_err(|e| e.to_string())?;
                        if &p.source != pin.event() || p.source_revision != pin.event().revision {
                            return Err(
                                "Terminal-processing preview names another exact source event."
                                    .into(),
                            );
                        }
                        Ok(Preview::TerminalProcessing(p))
                    }
                    (_, AnalysisResult::Error(e)) => Err(e),
                    _ => Err("Unexpected material preview reply.".into()),
                }
            })();
            match result {
                Ok(preview) => {
                    reviewed.set(Some(Arc::new(Reviewed {
                        source: pin,
                        settings: configured,
                        settings_revision: key,
                        preview,
                    })));
                    status.set("Review the actual strands and assumptions before opening reusable material.".into());
                }
                Err(e) => status.set(e),
            }
        });
    };
    let open = move |_| {
        let result = (|| -> Result<(), String> {
            let review = reviewed
                .get_untracked()
                .ok_or("Compute and review this operation first.")?;
            guard.run((review.source.clone(), review.settings_revision.clone()))?;
            let molecule = source.get_untracked().ok_or("The source is unavailable.")?;
            let frame = get_workspace.run(())?;
            let request = match (&review.settings, &review.preview) {
                (
                    Settings::Mismatch {
                        forward,
                        reverse,
                        options,
                    },
                    Preview::Mismatch(preview),
                ) => {
                    let index = selected
                        .get_untracked()
                        .parse::<usize>()
                        .ok()
                        .and_then(|v| v.checked_sub(1))
                        .ok_or("Explicitly choose a reported geometry number.")?;
                    MismatchChoice::from_preview(
                        review.source.clone(),
                        &molecule,
                        forward.clone(),
                        reverse.clone(),
                        options.clone(),
                        preview,
                        index,
                    )?
                    .request(
                        &frame,
                        &name.get_untracked(),
                        save_primers.get_untracked(),
                    )?
                }
                (
                    Settings::OverlapExtension { preparation },
                    Preview::OverlapExtension(preview),
                ) => OverlapExtensionChoice::from_preview(
                    review.source.clone(),
                    &molecule,
                    &frame,
                    preparation.clone(),
                    preview,
                )?
                .request(
                    &frame,
                    &name.get_untracked(),
                    save_primers.get_untracked(),
                )?,
                (Settings::TerminalProcessing { action }, Preview::TerminalProcessing(preview)) => {
                    TerminalChoice::from_preview(
                        review.source.clone(),
                        &molecule,
                        &frame,
                        *action,
                        preview,
                    )?
                    .request(&frame, &name.get_untracked())?
                }
                (Settings::Repair { method }, Preview::Repair(preview)) => {
                    RepairChoice::from_preview(review.source.clone(), &molecule, *method, preview)?
                        .request(&frame, &name.get_untracked())?
                }
                _ => return Err("The operation and preview differ.".into()),
            };
            let held = review.clone();
            let held_name = name.get_untracked();
            let held_selected = selected.get_untracked();
            let held_primers = save_primers.get_untracked();
            let bound_guard = Callback::new(move |()| -> Result<(), String> {
                guard.run((held.source.clone(), held.settings_revision.clone()))?;
                if name.get_untracked() != held_name
                    || selected.get_untracked() != held_selected
                    || save_primers.get_untracked() != held_primers
                {
                    return Err("The selected output or name changed; review again.".into());
                }
                Ok(())
            });
            on_open.run(crate::restriction::digest_product::MaterialOpenRequest {
                request,
                epoch: review.source.epoch(),
                guard: bound_guard,
            })
        })();
        if let Err(e) = result {
            status.set(e);
        }
    };
    let download_alive = alive.clone();
    let download = move |_| {
        let result = (|| -> Result<_, String> {
            let review = reviewed
                .get_untracked()
                .ok_or("Compute a material preview first.")?;
            guard.run((review.source.clone(), review.settings_revision.clone()))?;
            let bytes = match &review.preview {
                Preview::Mismatch(p) => serde_json::to_vec_pretty(p),
                Preview::Repair(p) => serde_json::to_vec_pretty(p),
                Preview::OverlapExtension(p) => serde_json::to_vec_pretty(p),
                Preview::TerminalProcessing(p) => serde_json::to_vec_pretty(p),
            }
            .map_err(|e| e.to_string())?;
            Ok((
                crate::io::download::prepare_download(
                    kind.evidence_name(),
                    "application/json",
                    &bytes,
                )?,
                review,
            ))
        })();
        match result {
            Ok((artifact, review)) => {
                let alive = download_alive.clone();
                let completion_pin = review.source.clone();
                let completion_key = review.settings_revision.clone();
                artifact.complete_guarded(
                    move || {
                        if !alive.load(Ordering::Relaxed) {
                            return Err("The preview was closed.".into());
                        }
                        guard.run((review.source.clone(), review.settings_revision.clone()))
                    },
                    move |result| {
                        if status.try_get_untracked().is_some()
                            && guard.run((completion_pin, completion_key)).is_ok()
                        {
                            if let Err(e) = result {
                                status.set(e);
                            }
                        }
                    },
                );
            }
            Err(e) => status.set(e),
        }
    };
    view! {<details class="physical-evidence"><summary>{kind.title()}</summary>
        <p class="helper">"Reconstruct the exact selected source and its retained chemistry. This creates a target-free material product; it does not verify or replace an intended assembly target."</p>
        <crate::history::run_origin_view::SourceOriginPicker source=source_info epoch get_workspace selected=Memo::new(move|_|origin.get()) on_select=select_origin/>
        <Show when=move||kind==Kind::EndRepair><p class="helper">"Choose an explicit repair method above and a saved material origin here. No restriction digest, PCR, nick sealing, or activation step is inserted."</p></Show>
        <button class="button secondary full" disabled=move||busy.get() on:click=run>{kind.preview_label()}</button>
        <p class="helper" role="status">{move||status.get()}</p>
        {move||reviewed.get().map(|review|{
            let (summary,details,incomplete,total,truncated)=match &review.preview {
                Preview::Mismatch(p)=>{
                    let report=&p.preview.report;
                    let products=report.products.iter().enumerate().take(100).map(|(i,g)|format!("Geometry {} · {} bases · certain {} · overlapping domains {}\n{}",i+1,g.sequence.len(),g.certain,g.overlapping_domains,snippet(g.sequence.as_str()))).collect::<Vec<_>>().join("\n\n");
                    (format!("{} geometries · search complete {} · enumeration complete {}. Initial physical PCR products have unmodified blunt 5′ OH / 3′ OH ends. Only a single supported certain non-overlapping geometry can be reusable. Incomplete search stays conditional and never claims isolation.",report.products.len(),report.search_complete,report.enumeration_complete),format!("Source chemistry:\n{}\n\n{}\n\n{}",p.source_duplex.describe(),products,p.findings.join("\n")),p.incomplete,p.findings_total,p.findings_truncated)
                },
                Preview::OverlapExtension(p)=>{
                    let stages=p.pcr_stages.iter().map(|stage|format!("{} · {} geometries · {}",stage.stage,stage.geometries.len(),stage.product_length.map(|n|format!("{n} bases")).unwrap_or_else(||"product withheld".into()))).collect::<Vec<_>>().join("\n");
                    let extension=p.extension.as_ref().map(|g|format!("Selected terminal overlap {} bp · pairings {:?} · search complete {} · unique {}",g.options.overlap,g.terminal_pairings.overlaps,g.terminal_pairings.complete,g.terminal_pairings.unique)).unwrap_or_else(||"Reciprocal extension withheld.".into());
                    let product=p.product.as_ref().map(|m|format!("Outer-PCR product · {} bases\n{}",m.sequence.len(),snippet(m.sequence.as_str()))).unwrap_or_else(||"Reusable outer-PCR product withheld: review the stage diagnostics.".into());
                    (format!("{} PCR stages · reusable outer product {}",p.pcr_stages.len(),p.reusable()),format!("{stages}\n\n{extension}\n\n{product}\n\n{}",p.findings.join("\n")),p.incomplete,p.findings_total,p.findings_truncated)
                },
                Preview::TerminalProcessing(p)=>{
                    let action=match p.action {
                        native::TerminalAction::TreatEnds{treatment}=>match treatment {
                            collagene_core::chemistry::EndTreatment::Phosphorylate=>"Phosphorylate accessible 5′ ends".to_string(),
                            collagene_core::chemistry::EndTreatment::Dephosphorylate=>"Dephosphorylate accessible 5′ ends".to_string(),
                        },
                        native::TerminalAction::Tail3Prime{base}=>format!("Add one {} to each eligible 3′ end",base.as_str()),
                    };
                    (format!("{action} · {} top-reference bases",p.product.sequence.len()),format!("Source chemistry:\n{}\n\nResult chemistry:\n{}\n\nResult top-reference DNA:\n{}\n\n{}",p.source_duplex.describe(),p.duplex.describe(),snippet(p.product.sequence.as_str()),p.findings.join("\n")),p.incomplete,p.findings_total,p.findings_truncated)
                },
                Preview::Repair(p)=>(format!("{} bases after explicit repair",p.molecule.sequence.len()),format!("Source chemistry:\n{}\n\nResult chemistry:\n{}\n\nTop 3′ added '{}' / removed '{}'; bottom 3′ added '{}' / removed '{}'\nResult DNA:\n{}\n\n{}",p.source_duplex.describe(),p.repair.duplex.summary().map(|v|v.describe()).unwrap_or_else(|e|e.to_string()),p.repair.top_added_3prime,p.repair.top_removed_3prime,p.repair.bottom_added_3prime,p.repair.bottom_removed_3prime,snippet(p.molecule.sequence.as_str()),p.findings.join("\n")),p.incomplete,p.findings_total,p.findings_truncated)
            };
            view!{<section><h4>{summary}</h4><pre>{details}</pre><p class="helper">{format!("Inherited/operation incompleteness: {incomplete}. {total} findings; excerpt truncated: {truncated}. At most 100 geometries and 240 bases per sequence are shown; the evidence download retains the full native preview.")}</p></section>}
        })}
        <Show when=move||reviewed.get().is_some()><Show when=move||kind==Kind::MismatchPcr><label>"Geometry number"<input aria-label="Reusable mismatch PCR geometry" type="number" min="1" prop:value=move||selected.get() on:input=move|e|selected.set(event_target_value(&e))/></label></Show><Show when=move||matches!(kind,Kind::MismatchPcr|Kind::OverlapExtensionPcr)><label class="checkbox-label"><input type="checkbox" prop:checked=move||save_primers.get() on:change=move|e|save_primers.set(event_target_checked(&e))/>"Save ordered primer collection"</label></Show><label>"Material product name"<input aria-label="Reusable material name" prop:value=move||name.get() on:input=move|e|name.set(event_target_value(&e))/></label><button class="button secondary full" on:click=download.clone()>"Download full material evidence"</button><button class="button primary full" disabled=move||busy.get()||reviewed.with(|r|r.as_ref().is_none_or(|r|match &r.preview {Preview::Mismatch(p)=>!p.reusable(),Preview::Repair(_)=>false,Preview::OverlapExtension(p)=>!p.reusable(),Preview::TerminalProcessing(_)=>false})) on:click=open.clone()>"Review reusable material installation"</button></Show>
    </details>}
}
