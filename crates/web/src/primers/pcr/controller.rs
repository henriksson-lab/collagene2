use crate::platform::worker;
use collagene_core::documents::transfer;
use collagene_core::material::lineage::EventBinding;
use collagene_core::protocol;
use collagene_core::sequence::analysis::Primer;
use collagene_core::verification;
use collagene_core::{
    chemistry,
    cloning::pcr::{material as material_pcr, mutagenesis, outputs as pcr_outputs, primers},
    documents::{history, workspace},
    formats::oligos as oligo_formats,
    material::lineage as run_lineage,
    sequence::analysis,
};

use leptos::prelude::*;

use wasm_bindgen_futures::spawn_local;

use crate::workspace::session::WorkspaceSession;

use crate::workspace::session::MAX_WORKSPACE_BYTES;

use crate::workspace::session::current_event_binding;

use crate::workspace::session::document_bytes;

use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::shell_state::ShellState;

#[derive(Clone, Copy)]
pub(crate) struct PcrState {
    pub(crate) pcr_forward_name: RwSignal<String>,
    pub(crate) pcr_reverse_name: RwSignal<String>,
    pub(crate) pcr_forward: RwSignal<String>,
    pub(crate) pcr_reverse: RwSignal<String>,
    pub(crate) pcr_forward_tail: RwSignal<String>,
    pub(crate) pcr_reverse_tail: RwSignal<String>,
    pub(crate) pcr_generation: RwSignal<u32>,
    pub(crate) pcr_review: RwSignal<Option<crate::primers::pcr::pcr_outputs::Review>>,
    pub(crate) pcr_origin: RwSignal<Option<crate::history::run_origin_model::PinnedOrigin>>,
    pub(crate) published_pcr_overlay:
        RwSignal<Option<crate::primers::pcr::pcr_overlay::Publication>>,
    pub(crate) pcr_output_busy: RwSignal<bool>,
    pub(crate) pcr_output_worker:
        StoredValue<Option<std::rc::Rc<worker::CachedWorker>>, leptos::prelude::LocalStorage>,
    pub(crate) pcr_alive: StoredValue<std::sync::Arc<std::sync::atomic::AtomicBool>>,
}
impl PcrState {
    pub(crate) fn workspace_replaced(&self) {
        // Invalidate in-flight callbacks immediately, before the epoch effect
        // clears displayed evidence and closes the existing output worker.
        self.pcr_generation
            .update(|generation| *generation = generation.wrapping_add(1));
    }
    pub(crate) fn new(
        session: WorkspaceSession,
        pending_pcr_source: RwSignal<Option<crate::primers::pcr::pcr_source_view::Snapshot>>,
    ) -> Self {
        let WorkspaceSession {
            documents,
            active,
            workspace_epoch,
            workspace_restoring,
            ..
        } = session;
        let pcr_forward_name = RwSignal::new("Forward".to_string());
        let pcr_reverse_name = RwSignal::new("Reverse".to_string());
        let pcr_forward = RwSignal::new(String::new());
        let pcr_reverse = RwSignal::new(String::new());
        let pcr_forward_tail = RwSignal::new(String::new());
        let pcr_reverse_tail = RwSignal::new(String::new());
        let pcr_generation = RwSignal::new(0u32);
        let pcr_review = RwSignal::new(None::<crate::primers::pcr::pcr_outputs::Review>);
        let pcr_origin = RwSignal::new(None::<crate::history::run_origin_model::PinnedOrigin>);
        Effect::new(move |_| {
            let current = documents.with(|docs| docs.get(active.get()).map(current_event_binding));
            if pcr_origin
                .get_untracked()
                .as_ref()
                .is_some_and(|origin| current.as_ref() != Some(origin.event()))
            {
                pcr_origin.set(None);
            }
        });
        let published_pcr_overlay =
            RwSignal::new(None::<crate::primers::pcr::pcr_overlay::Publication>);
        let pcr_output_busy = RwSignal::new(false);
        let pcr_output_worker = StoredValue::new_local(None::<std::rc::Rc<worker::CachedWorker>>);
        let pcr_alive = StoredValue::new(std::sync::Arc::new(std::sync::atomic::AtomicBool::new(
            true,
        )));
        let pcr_cleanup_alive = pcr_alive.get_value();
        on_cleanup(move || {
            pcr_cleanup_alive.store(false, std::sync::atomic::Ordering::Relaxed);
            if let Some(worker) = pcr_output_worker.get_value() {
                worker.close();
            }
        });
        Effect::new(move |_| {
            let _ = (
                pcr_forward_name.get(),
                pcr_reverse_name.get(),
                pcr_forward.get(),
                pcr_reverse.get(),
                pcr_forward_tail.get(),
                pcr_reverse_tail.get(),
                workspace_epoch.get(),
                workspace_restoring.get(),
                pcr_origin.get(),
            );
            documents.with(|docs| {
                let _ = docs.get(active.get()).map(current_event_binding);
            });
            pcr_generation.update(|g| *g = g.wrapping_add(1));
            pcr_review.set(None);
            if published_pcr_overlay.with_untracked(|p| {
                p.as_ref()
                    .and_then(|p| p.overlay.as_ref())
                    .is_some_and(|overlay| {
                        overlay.kind == crate::primers::pcr::pcr_overlay::Kind::ExactProduct
                    })
            }) {
                published_pcr_overlay.set(None);
            }
            pcr_output_busy.set(false);
            if let Some(worker) = pcr_output_worker.get_value() {
                worker.close();
            }
            pcr_output_worker.set_value(None);
        });
        Effect::new(move |_| {
            let _ = workspace_epoch.get();
            published_pcr_overlay.set(None);
            pending_pcr_source.set(None);
        });
        Self {
            pcr_forward_name,
            pcr_reverse_name,
            pcr_forward,
            pcr_reverse,
            pcr_forward_tail,
            pcr_reverse_tail,
            pcr_generation,
            pcr_review,
            pcr_origin,
            published_pcr_overlay,
            pcr_output_busy,
            pcr_output_worker,
            pcr_alive,
        }
    }
}
#[derive(Clone, Copy)]
pub(crate) struct PcrController {
    pub(crate) session: WorkspaceSession,
    pub(crate) viewer: ViewerSelection,
    pub(crate) shell: ShellState,
    pub(crate) state: PcrState,
}
impl PcrController {
    pub(crate) fn do_pcr(
        self,
        library_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
    ) -> impl Fn(web_sys::MouseEvent) + Copy + 'static {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        let ShellState { message, .. } = self.shell;
        let PcrState {
            pcr_forward_name,
            pcr_reverse_name,
            pcr_forward,
            pcr_reverse,
            pcr_forward_tail,
            pcr_reverse_tail,
            pcr_generation,
            pcr_review,
            pcr_origin,
            pcr_output_busy,
            pcr_output_worker,
            pcr_alive,
            ..
        } = self.state;
        move |_: web_sys::MouseEvent| {
            if [
                pcr_forward_name.get_untracked(),
                pcr_reverse_name.get_untracked(),
            ]
            .iter()
            .any(|name| {
                name.trim().is_empty() || name.len() > 256 || name.chars().any(char::is_control)
            }) {
                message.set("Primer names must be nonempty, at most 256 bytes, and contain no control characters.".into());
                return;
            }
            let primers = Primer::new(
                pcr_forward_name.get_untracked(),
                &pcr_forward.get_untracked(),
                &pcr_forward_tail.get_untracked(),
            )
            .and_then(|f| {
                Primer::new(
                    pcr_reverse_name.get_untracked(),
                    &pcr_reverse.get_untracked(),
                    &pcr_reverse_tail.get_untracked(),
                )
                .map(|r| (f, r))
            });
            let (forward, reverse) = match primers {
                Ok(p) => p,
                Err(error) => {
                    message.set(error.to_string());
                    return;
                }
            };
            let mol = molecule.get_untracked();
            let source_revision = mol.revision();
            let origin = pcr_origin.get_untracked();
            let revision = match origin.as_ref().map_or_else(
                || Ok(source_revision.clone()),
                |origin| {
                    collagene_core::cloning::pcr::material::revision(
                        &mol,
                        origin.as_ref(),
                        &forward,
                        &reverse,
                        collagene_core::chemistry::Strand::Top,
                    )
                    .map_err(|e| e.to_string())
                },
            ) {
                Ok(revision) => revision,
                Err(error) => {
                    message.set(error);
                    return;
                }
            };
            let pcr_source = documents
                .with_untracked(|docs| current_event_binding(&docs[active.get_untracked()]));
            if let Some(origin) = &origin {
                if origin.event() != &pcr_source {
                    message.set("PCR source changed; choose its material origin again.".into());
                    return;
                }
                let checked = library_snapshot.run(()).and_then(|workspace| {
                    origin
                        .as_ref()
                        .validate_for_workspace(&workspace, &pcr_source)
                        .map_err(|e| e.to_string())
                });
                if let Err(error) = checked {
                    message.set(error);
                    return;
                }
            }
            let pcr_epoch = workspace_epoch.get_untracked();
            let pcr_inputs = (
                pcr_forward_name.get_untracked(),
                pcr_reverse_name.get_untracked(),
                pcr_forward.get_untracked(),
                pcr_reverse.get_untracked(),
                pcr_forward_tail.get_untracked(),
                pcr_reverse_tail.get_untracked(),
            );
            let generation = pcr_generation.get_untracked().wrapping_add(1);
            pcr_generation.set(generation);
            pcr_review.set(None);
            pcr_output_busy.set(false);
            if let Some(worker) = pcr_output_worker.get_value() {
                worker.close();
            }
            pcr_output_worker.set_value(None);
            let reviewed_forward = forward.clone();
            let reviewed_reverse = reverse.clone();
            let request = protocol::request::Request {
                id: generation,
                revision: revision.clone(),
                job: if let Some(origin) = &origin {
                    protocol::schema::AnalysisJob::MaterialPcrDocuments {
                        template_strand: collagene_core::chemistry::Strand::Top,
                        molecule: mol,
                        origin: origin.as_ref().clone(),
                        forward,
                        reverse,
                    }
                } else {
                    protocol::schema::AnalysisJob::PcrDocuments {
                        template_strand: collagene_core::chemistry::Strand::Top,
                        molecule: mol,
                        forward,
                        reverse,
                    }
                },
            };
            message.set("Computing PCR products in background…".into());
            let alive = pcr_alive.get_value();
            spawn_local(async move {
                let result = worker::analyze_cancellable(request, alive.clone()).await;
                if !alive.load(std::sync::atomic::Ordering::Relaxed) {
                    return;
                }
                if workspace_restoring.get_untracked()
                    || workspace_epoch.get_untracked() != pcr_epoch
                    || !documents.with_untracked(|docs| {
                        docs.get(active.get_untracked())
                            .is_some_and(|doc| current_event_binding(doc) == pcr_source)
                    })
                    || molecule.get_untracked().revision() != source_revision
                    || pcr_origin.get_untracked() != origin
                    || pcr_generation.get_untracked() != generation
                    || pcr_inputs
                        != (
                            pcr_forward_name.get_untracked(),
                            pcr_reverse_name.get_untracked(),
                            pcr_forward.get_untracked(),
                            pcr_reverse.get_untracked(),
                            pcr_forward_tail.get_untracked(),
                            pcr_reverse_tail.get_untracked(),
                        )
                {
                    return;
                }
                match (origin.is_some(), result) {
                    (
                        false,
                        Ok(protocol::response::Response {
                            result: protocol::schema::AnalysisResult::PcrDocuments(products),
                            id,
                            revision: response_revision,
                        }),
                    )
                    | (
                        true,
                        Ok(protocol::response::Response {
                            result: protocol::schema::AnalysisResult::MaterialPcrDocuments(products),
                            id,
                            revision: response_revision,
                        }),
                    ) if id == generation && response_revision == revision => {
                        if products.iter().any(|document| {
                            document.geometry.template_strand
                                != collagene_core::chemistry::Strand::Top
                        }) {
                            message.set("PCR preview returned another source coordinate frame. Recompute the current preview.".into());
                            return;
                        }
                        let count = products.len();
                        pcr_review.set(Some(crate::primers::pcr::pcr_outputs::Review {
                            source: pcr_source,
                            origin,
                            epoch: pcr_epoch,
                            generation,
                            forward: reviewed_forward,
                            reverse: reviewed_reverse,
                            products: std::sync::Arc::new(products),
                            inputs: pcr_inputs,
                        }));
                        message.set(format!("{count} exact-annealing PCR geometries. Review and choose outputs; no documents have been created."));
                    }
                    (
                        _,
                        Ok(protocol::response::Response {
                            result: protocol::schema::AnalysisResult::Error(error),
                            ..
                        }),
                    )
                    | (_, Err(error)) => message.set(error),
                    _ => message.set("Unexpected worker response".into()),
                }
            });
        }
    }
    pub(crate) fn create_pcr_outputs(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        library_snapshot: Callback<
            (),
            Result<collagene_core::documents::workspace::Workspace, String>,
        >,
        oligo_tables: RwSignal<Vec<std::sync::Arc<collagene_core::formats::oligos::OligoTable>>>,
        pcr_review_current: Callback<crate::primers::pcr::pcr_outputs::Review, bool>,
        pending_material: RwSignal<Option<crate::restriction::digest_product::MaterialOpenRequest>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::primers::pcr::pcr_outputs::Choice, ()> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_epoch,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        let ShellState {
            message,
            modal,
            modal_error,
            ..
        } = self.shell;
        let PcrState {
            pcr_forward_name,
            pcr_reverse_name,
            pcr_forward,
            pcr_reverse,
            pcr_forward_tail,
            pcr_reverse_tail,
            pcr_review,
            pcr_origin,
            pcr_output_busy,
            pcr_output_worker,
            pcr_alive,
            ..
        } = self.state;
        Callback::new(move |choice: crate::primers::pcr::pcr_outputs::Choice| {
            if pcr_output_busy.get_untracked() || !pcr_review_current.run(choice.review.clone()) {
                message.set("The PCR draft or template changed; compute and review again.".into());
                return;
            }
            let Some(selected) = choice.review.products.get(choice.product_index) else {
                message.set("Choose a reported PCR geometry.".into());
                return;
            };
            let Some(source) = documents.with_untracked(|docs| {
                docs.get(active.get_untracked())
                    .map(|doc| doc.molecule.clone())
            }) else {
                return;
            };
            if choice.save_product {
                let prepared = (|| -> Result<crate::restriction::digest_product::MaterialOpenRequest, String> {
                let pin = crate::restriction::material_derivation_model::SourcePin::capture(
                    choice.review.source.clone(),
                    choice.review.epoch,
                    &source,
                    choice
                        .review
                        .origin
                        .as_ref()
                        .map(|origin| std::sync::Arc::new(origin.as_ref().clone())),
                )?;
                let captured = crate::restriction::material_derivation_model::PcrChoice::capture(
                    pin,
                    choice.review.forward.clone(),
                    choice.review.reverse.clone(),
                    choice.product_index,
                    std::sync::Arc::new(selected.clone()),
                )?;
                let frame = library_snapshot.run(())?;
                let request = captured.request(&frame, &choice.name, choice.save_primers)?;
                let bound = choice.review.clone();
                let guard = Callback::new(move |()| {
                    if modal.get() != "material_product"
                        || workspace_epoch.get() != bound.epoch
                        || pcr_origin.get() != bound.origin
                        || bound.inputs
                            != (
                                pcr_forward_name.get(),
                                pcr_reverse_name.get(),
                                pcr_forward.get(),
                                pcr_reverse.get(),
                                pcr_forward_tail.get(),
                                pcr_reverse_tail.get(),
                            )
                        || !documents.with(|docs| {
                            docs.get(active.get())
                                .is_some_and(|doc| current_event_binding(doc) == bound.source)
                        })
                    {
                        return Err("PCR source, primers or material origin changed; recompute and select again.".into());
                    }
                    Ok(())
                });
                Ok(crate::restriction::digest_product::MaterialOpenRequest {
                    request,
                    epoch: choice.review.epoch,
                    guard,
                })
            })();
                match prepared {
                    Ok(request) => {
                        pending_material.set(Some(request));
                        modal_error.set(String::new());
                        modal.set("material_product");
                    }
                    Err(error) => message.set(error),
                }
                return;
            }
            let request = collagene_core::cloning::pcr::outputs::Request {
                template_strand: collagene_core::chemistry::Strand::Top,
                source,
                forward: choice.review.forward.clone(),
                reverse: choice.review.reverse.clone(),
                selected_index: choice.product_index,
                expected_product_revision: selected.molecule.revision(),
                options: collagene_core::cloning::pcr::outputs::Options {
                    name: choice.name.clone(),
                    save_product: choice.save_product,
                    save_primers: choice.save_primers,
                },
            };
            let preview_revision = request.expected_product_revision.clone();
            let revision = match choice.review.origin.as_ref().map_or_else(
                || {
                    collagene_core::cloning::pcr::outputs::revision(&request)
                        .map_err(|e| e.to_string())
                },
                |origin| {
                    collagene_core::cloning::pcr::material::outputs_revision(
                        origin.as_ref(),
                        &request,
                    )
                    .map_err(|e| e.to_string())
                },
            ) {
                Ok(revision) => revision,
                Err(error) => {
                    message.set(error.to_string());
                    return;
                }
            };
            let w = match worker::CachedWorker::new() {
                Ok(w) => std::rc::Rc::new(w),
                Err(error) => {
                    message.set(error);
                    return;
                }
            };
            pcr_output_worker.set_value(Some(w.clone()));
            pcr_output_busy.set(true);
            message.set("Reconstructing the selected PCR geometry before creating outputs…".into());
            let alive = pcr_alive.get_value();
            spawn_local(async move {
                let result = w
                    .ask(protocol::request::Request {
                        id: choice.review.generation,
                        revision: revision.clone(),
                        job: if let Some(origin) = &choice.review.origin {
                            protocol::schema::AnalysisJob::PrepareMaterialPcrOutputs {
                                origin: origin.as_ref().clone(),
                                request,
                            }
                        } else {
                            protocol::schema::AnalysisJob::PreparePcrOutputs { request }
                        },
                    })
                    .await;
                w.close();
                if !alive.load(std::sync::atomic::Ordering::Relaxed)
                    || !pcr_review_current.run(choice.review.clone())
                {
                    return;
                }
                pcr_output_worker.set_value(None);
                pcr_output_busy.set(false);
                let outcome = (|| -> Result<(), String> {
                    let response = result?;
                    if response.id != choice.review.generation || response.revision != revision {
                        return Err("PCR output response does not match this request.".into());
                    }
                    let prepared =
                        match response.result {
                            protocol::schema::AnalysisResult::PcrOutputs(prepared)
                                if choice.review.origin.is_none() =>
                            {
                                prepared
                            }
                            protocol::schema::AnalysisResult::MaterialPcrOutputsPrepared(
                                prepared,
                            ) if choice.review.origin.is_some() => prepared,
                            protocol::schema::AnalysisResult::Error(error) => return Err(error),
                            _ => return Err("Unexpected PCR output response.".into()),
                        };
                    prepared
                        .validate_binding(&revision)
                        .map_err(|e| e.to_string())?;
                    if prepared.context_revision != revision
                        || prepared.selected_index != choice.product_index
                        || prepared.product.is_some() != choice.save_product
                        || prepared.primers.is_some() != choice.save_primers
                    {
                        return Err(
                            "PCR output selection does not match the reviewed request.".into()
                        );
                    }
                    let output_findings = prepared.findings.join(" ");
                    let mut staged_docs = documents.get_untracked();
                    let mut staged_active = active.get_untracked();
                    let mut staged_assets = attachments.get_untracked();
                    let mut product_event = None;
                    let mut primer_asset = None;
                    if let Some(mut product) = prepared.product {
                        product.metadata.insert(
                            "pcr:template_document".into(),
                            choice.review.source.document_id.clone(),
                        );
                        product.metadata.insert(
                            "pcr:template_event".into(),
                            choice.review.source.event_id.clone(),
                        );
                        let source = staged_docs
                            .iter()
                            .find(|doc| current_event_binding(doc) == choice.review.source)
                            .ok_or("The exact PCR template event changed; recompute products.")?;
                        let prefix = format!("pcr-{}", &product.revision()[..16]);
                        let mut suffix = 1usize;
                        let id = loop {
                            let id = format!("{prefix}-{suffix}");
                            if !staged_docs.iter().any(|doc| doc.id == id) {
                                break id;
                            }
                            suffix += 1;
                        };
                        let derived = source.derive_event(&choice.review.source.event_id,id,product,
                        "Selected exact-PCR reference preview; ordered primers retained, physical chemistry and experimental success unverified")
                        .map_err(|e|e.to_string())?;
                        product_event = Some(current_event_binding(&derived));
                        staged_active = staged_docs.len();
                        staged_docs.push(derived);
                    }
                    if staged_docs.len() > 1000
                        || staged_docs
                            .iter()
                            .map(|d| document_bytes(&d.molecule))
                            .sum::<usize>()
                            > MAX_WORKSPACE_BYTES
                    {
                        return Err("Selected PCR output exceeds the workspace document budget; nothing created.".into());
                    }
                    let retained_table = prepared.primers.map(std::sync::Arc::new);
                    if let Some(table) = &retained_table {
                        table.validate().map_err(|e| e.to_string())?;
                        let bytes =
                            serde_json::to_vec(table.as_ref()).map_err(|e| e.to_string())?;
                        let asset_hash = blake3::hash(&bytes).to_hex().to_string();
                        primer_asset = Some(asset_hash.clone());
                        staged_assets.insert(
                            asset_hash,
                            workspace::WorkspaceAttachment {
                                kind: "oligo-table-v1".into(),
                                name: format!("{} primers", choice.name),
                                media_type: "application/json".into(),
                                bytes,
                            },
                        );
                    }
                    let receipt = serde_json::to_vec(&serde_json::json!({
                    "schema_version":1,"source":choice.review.source,"context_revision":revision,
                    "selected_index":choice.product_index,"original_preview_revision":preview_revision,
                    "product":product_event,"primer_attachment":primer_asset,
                    "geometry":{"forward_start":prepared.geometry.forward_start,
                        "reverse_start":prepared.geometry.reverse_start,"crosses_origin":prepared.geometry.crosses_origin,
                        "forward_primer":prepared.geometry.forward_primer,"reverse_primer":prepared.geometry.reverse_primer},
                    "forward":choice.review.forward,"reverse":choice.review.reverse,
                    "output_name":choice.name,"save_product":choice.save_product,"save_primers":choice.save_primers,
                    "model":"Exact-annealing reference preview; no independent target verification or experimental success asserted"
                })).map_err(|e|e.to_string())?;
                    staged_assets.insert(
                        blake3::hash(&receipt).to_hex().to_string(),
                        workspace::WorkspaceAttachment {
                            kind: "pcr-output-reference-v1".into(),
                            name: format!("{} PCR source", choice.name),
                            media_type: "application/json".into(),
                            bytes: receipt,
                        },
                    );
                    WorkspaceSession::snapshot(
                        &staged_docs,
                        staged_active,
                        candidate.get_untracked().as_ref(),
                        &stored_candidates.get_untracked(),
                        &staged_assets,
                        &workspace_name.get_untracked(),
                    )?
                    .validate()
                    .map_err(|e| e.to_string())?;
                    batch(|| {
                        documents.set(staged_docs);
                        active.set(staged_active);
                        attachments.set(staged_assets);
                        if let Some(table) = retained_table {
                            oligo_tables.update(|tables| {
                                if !tables.iter().any(|old| {
                                    old.source_blake3 == table.source_blake3
                                        && old.options == table.options
                                }) {
                                    tables.push(table);
                                }
                            });
                        }
                        pcr_review.set(None);
                        message.set(format!("Selected PCR outputs created. Sequence-reference previews; physical chemistry and experimental success are unverified. {output_findings}"));
                    });
                    Ok(())
                })();
                if let Err(error) = outcome {
                    message.set(error);
                }
            });
        })
    }
    pub(crate) fn open_mismatch_outputs(
        &self,
        candidate: RwSignal<Option<transfer::Candidate>>,
        oligo_tables: RwSignal<Vec<std::sync::Arc<collagene_core::formats::oligos::OligoTable>>>,
        stored_candidates: RwSignal<Vec<transfer::Candidate>>,
    ) -> Callback<crate::cloning::methods::mutagenesis::OutputOpen, Result<(), String>> {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            workspace_name,
            attachments,
            ..
        } = self.session;
        let ViewerSelection { molecule, .. } = self.viewer;
        let ShellState { message, .. } = self.shell;
        Callback::new(
            move |opened: crate::cloning::methods::mutagenesis::OutputOpen| -> Result<(), String> {
                if workspace_restoring.get_untracked()
                    || opened.epoch != workspace_epoch.get_untracked()
                {
                    return Err(
                        "The workspace changed; compute and review the current PCR draft again."
                            .into(),
                    );
                }
                let mut staged_docs = documents.get_untracked();
                let mut staged_active = active.get_untracked();
                if staged_docs
                    .get(staged_active)
                    .is_none_or(|doc| current_event_binding(doc) != opened.target)
                    || !staged_docs
                        .iter()
                        .any(|doc| current_event_binding(doc) == opened.source)
                {
                    return Err(
                        "The original PCR source or viewed reference changed; compute again."
                            .into(),
                    );
                }
                opened.options.validate().map_err(|e| e.to_string())?;
                opened
                    .prepared
                    .validate_binding(&opened.revision)
                    .map_err(|e| e.to_string())?;
                let prepared = opened.prepared;
                if prepared.search.source_revision != opened.source.revision
                    || prepared.product.is_some() != opened.options.save_product
                    || prepared.primers.is_some() != opened.options.save_primers
                    || prepared
                        .product
                        .as_ref()
                        .is_some_and(|product| product.name != opened.options.name)
                {
                    return Err("PCR output choices differ from the reviewed request.".into());
                }
                let mut assets = attachments.get_untracked();
                let mut product_event = None;
                if let Some(mut product) = prepared.product {
                    product.metadata.insert(
                        "pcr:template_document".into(),
                        opened.source.document_id.clone(),
                    );
                    product
                        .metadata
                        .insert("pcr:template_event".into(), opened.source.event_id.clone());
                    let source = staged_docs
                        .iter()
                        .find(|doc| current_event_binding(doc) == opened.source)
                        .ok_or("Original PCR source event is unavailable")?;
                    let prefix = format!("mismatch-pcr-{}", &product.revision()[..16]);
                    let mut suffix = 1usize;
                    let id = loop {
                        let id = format!("{prefix}-{suffix}");
                        if !staged_docs.iter().any(|doc| doc.id == id) {
                            break id;
                        }
                        suffix += 1;
                    };
                    let derived=source.derive_event(&opened.source.event_id,id,product,
                "Selected intentional-mismatch PCR reference preview; no independent target or physical chemistry verification asserted")
                .map_err(|e|e.to_string())?;
                    product_event = Some(current_event_binding(&derived));
                    staged_active = staged_docs.len();
                    staged_docs.push(derived);
                }
                if staged_docs.len() > 1000
                    || staged_docs
                        .iter()
                        .map(|doc| document_bytes(&doc.molecule))
                        .sum::<usize>()
                        > MAX_WORKSPACE_BYTES
                {
                    return Err("The selected PCR output exceeds the workspace document budget; nothing created.".into());
                }
                let table = prepared.primers.map(std::sync::Arc::new);
                let mut primer_asset = None;
                if let Some(table) = &table {
                    table.validate().map_err(|e| e.to_string())?;
                    let bytes = serde_json::to_vec(table.as_ref()).map_err(|e| e.to_string())?;
                    let hash = blake3::hash(&bytes).to_hex().to_string();
                    primer_asset = Some(hash.clone());
                    assets.insert(
                        hash,
                        workspace::WorkspaceAttachment {
                            kind: "oligo-table-v1".into(),
                            name: format!("{} primers", opened.options.name),
                            media_type: "application/json".into(),
                            bytes,
                        },
                    );
                }
                let geometry = prepared.geometry;
                let bytes=serde_json::to_vec(&serde_json::json!({
            "schema_version":1,"source":opened.source,"viewed_reference":opened.target,
            "context_revision":opened.revision,"selected_index":prepared.selected_index,
            "outputs":opened.options,"product":product_event,"primer_attachment":primer_asset,
            "search":prepared.search,"findings":prepared.findings,
            "geometry":{"forward_start":geometry.forward_start,"reverse_start":geometry.reverse_start,
                "template_length":geometry.template_length,"crosses_origin":geometry.crosses_origin,
                "forward_primer_index":geometry.forward_primer_index,"reverse_primer_index":geometry.reverse_primer_index,
                "forward_tail_length":geometry.forward_tail_length,"reverse_tail_length":geometry.reverse_tail_length,
                "forward_binding":geometry.forward_binding,"reverse_binding":geometry.reverse_binding,
                "substitutions":geometry.substitutions,"overlapping_domains":geometry.overlapping_domains,"certain":geometry.certain},
            "model":"Intentional-mismatch sequence-reference preview; search limits retained. No independent target, physical chemistry or experimental success verification asserted."
        })).map_err(|e|e.to_string())?;
                assets.insert(
                    blake3::hash(&bytes).to_hex().to_string(),
                    workspace::WorkspaceAttachment {
                        kind: "mismatch-pcr-output-v1".into(),
                        name: format!("{} PCR source", opened.options.name),
                        media_type: "application/json".into(),
                        bytes,
                    },
                );
                WorkspaceSession::snapshot(
                    &staged_docs,
                    staged_active,
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &assets,
                    &workspace_name.get_untracked(),
                )?
                .validate()
                .map_err(|e| e.to_string())?;
                batch(|| {
                    documents.set(staged_docs);
                    active.set(staged_active);
                    attachments.set(assets);
                    if let Some(table) = table {
                        oligo_tables.update(|tables| {
                            if !tables.iter().any(|old| {
                                old.source_blake3 == table.source_blake3
                                    && old.options == table.options
                            }) {
                                tables.push(table);
                            }
                        });
                    }
                    message.set("Selected mismatch-PCR reference outputs retained with the original source event. Search limitations and conditional substitutions remain evidence, not a successful cloning claim.".into());
                });
                Ok(())
            },
        )
    }
    pub(crate) fn navigate_pcr_overlay(
        &self,
        pending_pcr_source: RwSignal<Option<crate::primers::pcr::pcr_source_view::Snapshot>>,
    ) -> Callback<
        (
            collagene_core::material::lineage::EventBinding,
            u64,
            Vec<(usize, usize)>,
        ),
        Result<(), String>,
    > {
        let session = self.session;
        let WorkspaceSession {
            documents,
            active,
            workspace_restoring,
            workspace_epoch,
            ..
        } = self.session;
        let ViewerSelection {
            selection,
            selected_feature,
            search_regions,
            molecule,
            ..
        } = self.viewer;
        let ShellState { modal, .. } = self.shell;
        let PcrState {
            published_pcr_overlay,
            ..
        } = self.state;
        Callback::new(
            move |(source, epoch, spans): (
                collagene_core::material::lineage::EventBinding,
                u64,
                Vec<(usize, usize)>,
            )|
                  -> Result<(), String> {
                if workspace_restoring.get_untracked() || epoch != workspace_epoch.get_untracked() {
                    return Err("The PCR workspace changed; select the product again.".into());
                }
                let publication = published_pcr_overlay
                    .get_untracked()
                    .ok_or("Select a current PCR geometry first.")?;
                let overlay = publication
                    .overlay
                    .ok_or("The PCR overlay is unavailable.")?;
                if !overlay.is_current(&source, epoch)
                    || spans.is_empty()
                    || spans
                        .iter()
                        .any(|span| !overlay.spans.iter().any(|s| (s.start, s.end) == *span))
                {
                    return Err("The selected PCR geometry changed; select it again.".into());
                }
                let (index, molecule) = documents
                    .with_untracked(|docs| {
                        docs.iter()
                            .enumerate()
                            .find(|(_, d)| current_event_binding(d) == source)
                            .map(|(i, d)| (i, d.molecule.clone()))
                    })
                    .ok_or("The exact PCR template event changed; compute products again.")?;
                if index == active.get_untracked() {
                    batch(|| {
                        selected_feature.set(None);
                        selection.set(if spans.len() == 1 {
                            spans.first().copied()
                        } else {
                            None
                        });
                        search_regions.set(Some((source.revision, spans)));
                    });
                } else {
                    pending_pcr_source.set(Some(crate::primers::pcr::pcr_source_view::Snapshot {
                        owner: publication.owner,
                        overlay,
                        molecule: std::sync::Arc::new(molecule),
                        spans,
                    }));
                    modal.set("pcr_source");
                }
                Ok(())
            },
        )
    }
}

/// Editable primer-pair intent and pending primer/mutagenesis seeds.
#[derive(Clone, Copy)]
pub(crate) struct PrimerDesignState {
    pub(crate) pair_draft: RwSignal<crate::primers::primer_pair_model::PairDraft>,
    pub(crate) primer_seed: RwSignal<Option<collagene_core::sequence::analysis::Primer>>,
    pub(crate) mutation_seed:
        RwSignal<Option<crate::cloning::methods::mutagenesis::MutationSelection>>,
}
impl PrimerDesignState {
    pub(crate) fn new() -> Self {
        Self {
            pair_draft: RwSignal::new(crate::primers::primer_pair_model::PairDraft::default()),
            primer_seed: RwSignal::new(None::<Primer>),
            mutation_seed: RwSignal::new(
                None::<crate::cloning::methods::mutagenesis::MutationSelection>,
            ),
        }
    }
    pub(crate) fn restore(&self, draft: crate::primers::primer_pair_model::PairDraft) {
        self.pair_draft.set(draft);
        self.primer_seed.set(None);
        self.mutation_seed.set(None);
    }
}
