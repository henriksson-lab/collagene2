use crate::reads::controller::{ReadReviewActions, ReadReviewState};
use crate::viewer::sequence::editing::SequenceEditor;
use crate::workspace::autosave::WorkspaceAutosave;
use crate::workspace::tabs::{DocumentItems, DocumentTabs};
#[cfg(target_arch = "wasm32")]
pub(crate) mod shell;

use crate::viewer::map::renderer::MapRenderer;
use crate::workspace::session::RestoredWorkspace;
use crate::workspace::session::WorkspaceSession;
use collagene_core::documents::transfer;
use collagene_core::formats::sequence as formats;
use collagene_core::planning::plan::DesignPlan;
use collagene_core::planning::plan::Operation;
use collagene_core::planning::plan::Source;
use collagene_core::protocol;
use collagene_core::sequence::analysis::Primer;
use collagene_core::{
    chemistry::{self, methylation, tracked_duplex as modified_duplex},
    documents::workspace,
    formats::{oligos as oligo_formats, snapgene},
    review::reads::traces,
    sequence::editing,
    verification,
};
use collagene_model::Alphabet;
use collagene_model::Feature;
use collagene_model::Location;
use collagene_model::Molecule;
use collagene_model::Topology;
use leptos::html;
use leptos::prelude::*;
use wasm_bindgen::JsCast;
use web_sys::HtmlCanvasElement;
use web_sys::HtmlInputElement;

use crate::app::shell::dismiss_viewer_menus;
use crate::viewer::sequence::selection::ViewerSelection;
use crate::workspace::session::MAX_WORKSPACE_BYTES;
use crate::workspace::session::bound_primer_source;
use crate::workspace::session::current_event_binding;
use crate::workspace::session::document_bytes;
use crate::workspace::session::open_histories;
use crate::workspace::session::recombination_source_info;
use crate::workspace::shell_state::ShellState;

use crate::cloning::preparation::restriction_fields::RestrictionFields;
use crate::protocol::step_description::describe_step;
use demo::demonstration;
mod demo;
use crate::assembly::controller::AssemblyController;
use crate::biology::controller::BiologyController;
use crate::io::controller::ImportController;
use crate::library::controller::LibraryController;
use crate::primers::pcr::controller::PcrController;
use crate::primers::pcr::controller::PcrState;
use crate::protocol::controller::ProtocolController;
use crate::reads::controller::ReadController;
use crate::restriction::controller::RestrictionController;
use crate::viewer::alignment::controller::ComparisonController;
use crate::viewer::sequence::annotation_controller::AnnotationController;
use crate::workspace::documents::DocumentController;
use crate::workspace::reference_controller::ReferenceController;
#[component]
pub fn App() -> impl IntoView {
    let exports = crate::io::download::ExportService::install();
    let native_export_notice = exports.notice;
    crate::platform::transient_menus::install();
    let session = WorkspaceSession::new(demonstration());
    let WorkspaceSession {
        documents,
        active,
        workspace_restoring,
        workspace_epoch,
        workspace_name,
        attachments,
    } = session;

    let genome_links = Memo::new(move |_| {
        attachments.with(|assets| {
            assets
                .values()
                .filter(|a| a.kind == collagene_core::tracks::genome::WORKSPACE_ATTACHMENT_KIND)
                .filter_map(|a| collagene_core::tracks::genome::load_link(&a.bytes).ok())
                .collect::<Vec<_>>()
        })
    });
    let run_history_state = crate::protocol::controller::RunHistoryState::new();
    let crate::protocol::controller::RunHistoryState {
        recorded_runs,
        run_page,
        pending_run_event,
        hidden_products,
        show_prepared_history,
    } = run_history_state;
    let read_archive_state = crate::reads::controller::ReadArchiveState::new();
    let crate::reads::controller::ReadArchiveState {
        saved_assemblies,
        assembly_page,
        pending_assembly,
        inspected_assembly,
    } = read_archive_state;

    let assembly_state = crate::assembly::controller::AssemblyState::new(workspace_epoch);
    let crate::assembly::controller::AssemblyState {
        method_config,
        assembly_draft,
        previous_assembly_draft,
        assembly_visible,
        pending_finalize,
        pending_release,
        selected_release_manifest,
        pending_donor,
    } = assembly_state;
    let digest_workspace_state = crate::restriction::controller::DigestWorkspaceState::new();
    let crate::restriction::controller::DigestWorkspaceState {
        pending_digest,
        pending_material,
        digest_boundary_pick,
    } = digest_workspace_state;

    let gel_digest_seeds = RwSignal::new(Vec::<crate::restriction::gel::DigestLaneSeed>::new());

    let selected_run = RwSignal::new(String::new());
    let reference_transfer_state =
        crate::workspace::reference_controller::ReferenceTransferState::new();
    let crate::workspace::reference_controller::ReferenceTransferState {
        internal_clipboard,
        pending_reference_paste,
    } = reference_transfer_state;

    let comparison_state = crate::viewer::alignment::controller::ComparisonState::new();
    let crate::viewer::alignment::controller::ComparisonState {
        alignments,
        pending_alignment,
        pending_comparison_source,
    } = comparison_state;
    let pending_translation_source =
        RwSignal::new(None::<crate::viewer::alignment::pairwise_source_view::Snapshot>);
    let pending_pcr_source = RwSignal::new(None::<crate::primers::pcr::pcr_source_view::Snapshot>);

    let import_state = crate::io::controller::ImportState::new();
    let crate::io::controller::ImportState {
        oligo_tables,
        existing_oligo,
        pending_oligo,
        pending_fastq,
        pending_sbol,
        pending_variants,
        pending_read_alignment,
    } = import_state;
    let annotation_state = crate::viewer::sequence::annotation_controller::AnnotationState::new();
    let crate::viewer::sequence::annotation_controller::AnnotationState {
        libraries,
        annotation_source,
    } = annotation_state;

    let web_mcp_assembly = RwSignal::new(None::<crate::ai::web_mcp_assembly::bridge::Endpoint>);

    let read_state = ReadReviewState::new();
    let ReadReviewState {
        read_reference,
        pending_read_check,
        reviewed_check,
        review_seed,
        pending_read_review,
    } = read_state;

    let digest_run_origin = RwSignal::new(None::<crate::history::run_origin_model::PinnedOrigin>);
    Effect::new(move || {
        let _ = workspace_epoch.get();
        digest_run_origin.set(None);
    });

    let annotation_format =
        RwSignal::new(collagene_core::formats::annotations::AnnotationFormat::Gff3);

    let fastq_encoding = RwSignal::new(None::<oligo_formats::QualityEncoding>);
    let plan_state = crate::protocol::controller::PlanState::new();
    let crate::protocol::controller::PlanState {
        candidate,
        stored_candidates,
        candidate_events,
        plan_report,
        last_plan,
        last_plan_target,
        last_plan_events,
        plan_generation,
    } = plan_state;

    // Working tabs are presentation state, separate from retained histories and
    // the archive's last working document index (which must remain valid).
    let tabs = DocumentTabs::new(session);
    let open_tabs = tabs.open_tabs;
    let has_open_document = tabs.has_open_document;
    let molecule = Memo::new(move |_| {
        documents.with(|docs| docs[active.get().min(docs.len() - 1)].molecule.clone())
    });
    let molecule_revision = Memo::new(move |_| molecule.with(Molecule::revision));
    let selection = RwSignal::new(None::<(usize, usize)>);
    let selected_feature = RwSignal::new(None::<usize>);
    let feature_for_editor = RwSignal::new(None::<usize>);
    let feature_initial_location = RwSignal::new(None::<Location>);
    let primer_design_state = crate::primers::pcr::controller::PrimerDesignState::new();
    let crate::primers::pcr::controller::PrimerDesignState {
        pair_draft,
        primer_seed,
        mutation_seed,
    } = primer_design_state;

    let inspected_feature = Memo::new(move |_| {
        molecule.with(|m| {
            selected_feature
                .get()
                .and_then(|i| m.features.get(i).cloned())
        })
    });
    let search_regions = RwSignal::new(None::<(String, Vec<(usize, usize)>)>);
    let viewer = ViewerSelection {
        selection,
        selected_feature,
        search_regions,
        molecule,
        has_open_document,
    };
    // Navigation receipts are local UI state, never imported scientific evidence.
    // Initialize before effects so the first tick cannot erase a pointer selection.
    let SequenceEditor {
        sequence_transitions,
        sequence_on_edit,
        sequence_on_history,
    } = SequenceEditor::new(session, viewer);
    let biology_directed =
        RwSignal::new(None::<crate::biology::biology_selection::DirectedSelection>);
    let highlighted_spans = Memo::new(move |_| {
        if let Some(feature) = inspected_feature.get() {
            return feature
                .location
                .spans()
                .into_iter()
                .map(|(a, b, _)| (a, b))
                .collect::<Vec<_>>();
        }
        if let Some(span) = selection.get() {
            return vec![span];
        }
        search_regions.with(|regions| {
            regions
                .as_ref()
                .filter(|(revision, _)| *revision == molecule_revision.get())
                .map(|(_, parts)| parts.clone())
                .unwrap_or_default()
        })
    });
    let selection_status = Memo::new(move |_| {
        if let Some((a, b)) = selection.get().filter(|(a, b)| a == b) {
            let _ = b;
            return format!("Insertion boundary {a} · type to insert");
        }
        let parts = highlighted_spans.get();
        if parts.len() > 1 {
            let positions = parts
                .iter()
                .take(4)
                .map(|(a, b)| format!("{}–{}", a + 1, b))
                .collect::<Vec<_>>()
                .join(", ");
            return format!(
                "Selection {positions}{} · {} reference spans",
                if parts.len() > 4 { "…" } else { "" },
                parts.len()
            );
        }
        selection
            .get()
            .map(|(a, b)| {
                format!(
                    "Selection {}–{} · {} {}",
                    a + 1,
                    b,
                    b - a,
                    molecule.with(|m| m.sequence.alphabet().map_unit())
                )
            })
            .unwrap_or_else(|| {
                molecule.with(|m| {
                    if m.sequence.alphabet() == Alphabet::Protein {
                        "Click a residue or feature to select".into()
                    } else {
                        "Click a base or feature to select".into()
                    }
                })
            })
    });

    let hovered_restriction = RwSignal::new(None::<usize>);
    let shell = ShellState::new();
    let ShellState {
        panel,
        message,
        modal,
        modal_error,
    } = shell;

    let open_bulk_release = Callback::new(move |hash: String| {
        if workspace_restoring.get_untracked() {
            return;
        }
        pending_release.set(Some(hash));
        modal_error.set(String::new());
        modal.set("combinatorial_release");
    });
    let released_bulk = Callback::new(move |hash: String| {
        selected_release_manifest.set(Some(hash));
        assembly_visible.set(true);
        pending_release.set(None);
        modal.set("");
    });
    let import_text = RwSignal::new(String::new());
    let pending_compat_import = RwSignal::new(
        None::<(
            std::sync::Arc<Vec<u8>>,
            String,
            crate::io::sequence_compat_import::Profile,
        )>,
    );
    let edit_text = RwSignal::new(String::new());
    let trace_view = RwSignal::new(false);
    let file_input = NodeRef::<html::Input>::new();
    let gl_canvas = NodeRef::<html::Canvas>::new();
    let label_canvas = NodeRef::<html::Canvas>::new();
    let renderer = StoredValue::new_local(None::<MapRenderer>);
    let render_error = RwSignal::new(String::new());
    let map_labels = RwSignal::new(crate::viewer::map::map_labels::LabelStats::default());
    let editor_view = RwSignal::new("split");
    let show_library = RwSignal::new(true);
    let layout = crate::workspace::workspace_layout::use_layout(
        editor_view,
        show_library,
        Memo::new(move |_| workspace_restoring.get()),
    );
    let restriction_identity = Memo::new(move |_| {
        documents.with(|docs| {
            let doc = &docs[active.get()];
            format!("{}:{}:{}", workspace_epoch.get(), doc.id, doc.head_event())
        })
    });
    let annotation_hidden = Memo::new(move |_| {
        documents.with(|docs| {
            let doc = &docs[active.get()];
            attachments.with(|assets| {
                match crate::viewer::sequence::annotation_visibility::load(
                    assets,
                    &current_event_binding(doc),
                    doc.molecule.features.len(),
                ) {
                    Ok(hidden) => hidden,
                    Err(error) => {
                        message.set(format!(
                            "Annotation display preferences need attention: {error}"
                        ));
                        Default::default()
                    }
                }
            })
        })
    });
    let restriction_selection = Memo::new(move |_| {
        molecule.with(|m| {
            if let Some(index) = selected_feature.get() {
                return m.features.get(index).and_then(|feature| {
                    collagene_core::cloning::restriction::selection::footprint(
                        &feature.location,
                        m.sequence.len(),
                    )
                    .ok()
                });
            }
            if let Some((revision, spans)) = search_regions.get() {
                if revision == molecule_revision.get() && !spans.is_empty() {
                    return Some(spans);
                }
            }
            selection.get().and_then(|(start, end)| {
                (start < end && end <= m.sequence.len()).then_some(vec![(start, end)])
            })
        })
    });
    let restrictions = crate::restriction::restriction_view::use_restrictions(
        molecule,
        restriction_identity,
        Memo::new(move |_| workspace_restoring.get()),
        restriction_selection,
    );
    let navigate_restriction = RestrictionController {
        session,
        viewer,
        shell,
    }
    .navigate_restriction(
        assembly_visible,
        digest_boundary_pick,
        molecule_revision,
        restrictions,
    );
    let DocumentItems {
        document_items,
        explorer_active,
        explorer_documents,
        tab_items,
        explorer_locked,
        explorer_activate,
    } = tabs.items(show_prepared_history, hidden_products);
    let close_document_tab = DocumentController {
        session,
        viewer,
        shell,
    }
    .close_document_tab(assembly_visible, open_tabs, show_library);

    let autosave = WorkspaceAutosave::new(session, shell);
    let WorkspaceAutosave {
        db,
        native_library_context,
        loaded,
        save_status,
        save_generation,
        ..
    } = autosave;
    let pcr_state = PcrState::new(session, pending_pcr_source);
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
        published_pcr_overlay,
        pcr_output_busy,
        pcr_output_worker,
        pcr_alive,
    } = pcr_state;
    let pcr_overlay = Memo::new(move |_| {
        if workspace_restoring.get() {
            return None;
        }
        let epoch = workspace_epoch.get();
        let inputs = (
            pcr_forward_name.get(),
            pcr_reverse_name.get(),
            pcr_forward.get(),
            pcr_reverse.get(),
            pcr_forward_tail.get(),
            pcr_reverse_tail.get(),
        );
        documents.with(|docs| {
            let doc = &docs[active.get()];
            let binding = current_event_binding(doc);
            if let Some(overlay) = published_pcr_overlay.get().and_then(|p| p.overlay) {
                let fresh = overlay.is_current(&binding, epoch)
                    && (overlay.kind != crate::primers::pcr::pcr_overlay::Kind::ExactProduct
                        || pcr_review.with(|review| {
                            review.as_ref().is_some_and(|r| {
                                r.source == binding && r.epoch == epoch && r.inputs == inputs
                            })
                        }));
                if fresh {
                    return Some(overlay);
                }
            }
            let forward = Primer::new(&inputs.0, &inputs.2, &inputs.4).ok()?;
            let reverse = Primer::new(&inputs.1, &inputs.3, &inputs.5).ok()?;
            pair_draft.with(|draft| {
                crate::primers::pcr::pcr_overlay::captured(
                    &doc.molecule,
                    &binding,
                    epoch,
                    draft,
                    &forward,
                    &reverse,
                )
                .ok()
                .flatten()
                .map(std::sync::Arc::new)
            })
        })
    });
    let publish_pcr_overlay = Callback::new(
        move |publication: crate::primers::pcr::pcr_overlay::Publication| match publication
            .overlay
            .as_ref()
        {
            Some(overlay) => {
                let valid = !workspace_restoring.get_untracked()
                    && workspace_epoch.get_untracked() == overlay.epoch
                    && documents.with_untracked(|docs| {
                        docs.iter().any(|d| {
                            current_event_binding(d) == overlay.source
                                && overlay.spans.iter().all(|span| {
                                    span.start < span.end && span.end <= d.molecule.sequence.len()
                                })
                        })
                    });
                if valid {
                    published_pcr_overlay.set(Some(publication));
                }
            }
            None => {
                if published_pcr_overlay
                    .with_untracked(|p| p.as_ref().is_some_and(|p| p.owner == publication.owner))
                {
                    published_pcr_overlay.set(None);
                }
            }
        },
    );
    let primer_sources = Memo::new(move |_| {
        documents.with(|docs| {
            docs.iter()
                .filter(|doc| doc.molecule.sequence.alphabet() == Alphabet::Dna)
                .map(recombination_source_info)
                .collect::<Vec<_>>()
        })
    });
    let list_saved_primers = Callback::new(move |binding: collagene_core::material::lineage::EventBinding| -> Result<collagene_core::cloning::pcr::reuse::List,String> {
        if workspace_restoring.get_untracked() {return Err("Wait for workspace restoration to finish.".into());}
        documents.with_untracked(|docs| {
            let source=bound_primer_source(docs,&binding)?;
            let list=collagene_core::cloning::pcr::reuse::list(source,1000).map_err(|e|e.to_string())?;
            if list.source_revision!=binding.revision {return Err("The saved-primer source revision changed.".into());}
            Ok(list)
        })
    });
    let resolve_saved_primer = Callback::new(move |request: crate::primers::primer_picker::ResolveRequest| -> Result<collagene_core::cloning::pcr::reuse::Resolved,String> {
        if workspace_restoring.get_untracked() {return Err("Wait for workspace restoration to finish.".into());}
        documents.with_untracked(|docs| {
            let source=bound_primer_source(docs,&request.source)?;
            collagene_core::cloning::pcr::reuse::resolve_at_revision(source,request.feature_index,&request.source.revision).map_err(|e|e.to_string())
        })
    });

    Effect::new(move |_| {
        if panel.get() == "reads" && read_reference.get().is_none() {
            read_reference.set(Some(documents.with_untracked(|docs| {
                current_event_binding(&docs[active.get_untracked()])
            })));
        }
    });
    let read_check_source = Memo::new(move |_| {
        let pinned = read_reference.get();
        documents.with(|docs| {
            let doc = pinned
                .as_ref()
                .and_then(|pin| docs.iter().find(|doc| doc.id == pin.document_id))
                .unwrap_or_else(|| &docs[active.get()]);
            recombination_source_info(doc)
        })
    });
    let read_check_attachment_available = Callback::new(move |hash: String| {
        attachments.with(|assets| {
            assets.get(&hash).is_some_and(|asset| {
                asset.kind == collagene_core::review::reads::checks::ATTACHMENT_KIND
            })
        })
    });
    let read_check_reference = Memo::new(move |_| read_check_source.get().binding);
    let read_check_molecule = Memo::new(move |_| {
        let id = read_check_source.get().binding.document_id;
        documents.with(|docs| {
            docs.iter()
                .find(|doc| doc.id == id)
                .unwrap_or_else(|| &docs[active.get()])
                .molecule
                .clone()
        })
    });
    let primer_context =
        Memo::new(move |_| documents.with(|docs| current_event_binding(&docs[active.get()])));
    let pair_selection = Callback::new(move |()| {
        if workspace_restoring.get_untracked() {
            return Err("Wait for the workspace operation to finish.".into());
        }
        documents.with_untracked(|docs| {
            let source = &docs[active.get_untracked()];
            let (start, end) = crate::viewer::sequence::selection_toolbar::primer_range(
                &source.molecule,
                selection.get_untracked(),
                selected_feature.get_untracked(),
                &highlighted_spans.get_untracked(),
            )?;
            crate::primers::primer_pair_model::SelectionSeed::capture(
                &source.molecule,
                &current_event_binding(source),
                start,
                end,
            )
        })
    });
    let navigate_primer_binding = Callback::new(
        move |(source, spans): (
            collagene_core::material::lineage::EventBinding,
            Vec<(usize, usize)>,
        )|
              -> Result<(), String> {
            if workspace_restoring.get_untracked() || primer_context.get_untracked() != source {
                return Err(
                    "The source event changed; rerun primer analysis before navigating.".into(),
                );
            }
            let n = molecule.with_untracked(|m| m.sequence.len());
            if spans.is_empty()
                || spans.len() > 2
                || spans.iter().any(|(start, end)| start >= end || *end > n)
            {
                return Err("The primer placement has invalid source spans.".into());
            }
            batch(|| {
                selected_feature.set(None);
                selection.set(if spans.len() == 1 {
                    spans.first().copied()
                } else {
                    None
                });
                search_regions.set(Some((source.revision, spans)));
            });
            Ok(())
        },
    );
    // Optional browser tools reuse the same exact-event navigation command as
    // the editor. No document/capsule copy or scientific job runs at registration.
    crate::ai::web_mcp::install(Callback::new(move |(tool, arguments): (String, String)| {
        if workspace_restoring.get_untracked() {
            return Err("The workspace is still loading.".into());
        }
        use crate::ai::web_mcp_assembly as assembly;
        use crate::ai::web_mcp_assembly::bridge::Command;
        if tool == "collagene_open_assembly_workspace" {
            crate::ai::web_mcp::arguments::<crate::ai::web_mcp::Empty>(&arguments)?;
            assembly_visible.set(true);
            return crate::ai::web_mcp::response(
                &serde_json::json!({"opened":true,"next":"Read assembly context after the workspace mounts, then request intended-reference review with its exact pin."}),
            );
        }
        if tool == "collagene_assembly_context" {
            crate::ai::web_mcp::arguments::<crate::ai::web_mcp::Empty>(&arguments)?;
            let mut context = assembly_draft.with_untracked(|draft| {
                assembly::context(draft, workspace_epoch.get_untracked())
            })?;
            let expected: assembly::Pin =
                serde_json::from_value(context["pin"].clone()).map_err(|e| e.to_string())?;
            context["review"] = match web_mcp_assembly.get_untracked() {
                Some(endpoint) if endpoint.alive.load(std::sync::atomic::Ordering::Relaxed) => {
                    endpoint.command.run(Command::Status(expected))?
                }
                _ => {
                    serde_json::json!({"phase":"workspace_closed","message":"Open the assembly workspace to review this draft."})
                }
            };
            return crate::ai::web_mcp::response(&context);
        }
        if matches!(
            tool.as_str(),
            "collagene_review_assembly_draft"
                | "collagene_cancel_assembly_review"
                | "collagene_propose_assembly_edits"
                | "collagene_discard_assembly_proposal"
        ) {
            let endpoint = web_mcp_assembly
                .get_untracked()
                .filter(|endpoint| endpoint.alive.load(std::sync::atomic::Ordering::Relaxed))
                .ok_or("Open the assembly workspace before requesting or cancelling review.")?;
            let command = match tool.as_str() {
                "collagene_review_assembly_draft" => {
                    let request: assembly::ReviewRequest =
                        crate::ai::web_mcp::arguments(&arguments)?;
                    Command::Review(request.expected)
                }
                "collagene_cancel_assembly_review" => {
                    let request: assembly::CancelRequest =
                        crate::ai::web_mcp::arguments(&arguments)?;
                    Command::Cancel(request.expected, request.review)
                }
                "collagene_propose_assembly_edits" => {
                    Command::Propose(crate::ai::web_mcp::arguments(&arguments)?)
                }
                "collagene_discard_assembly_proposal" => {
                    Command::Discard(crate::ai::web_mcp::arguments(&arguments)?)
                }
                _ => unreachable!(),
            };
            return crate::ai::web_mcp::response(&endpoint.command.run(command)?);
        }
        let (response, navigation, source) = documents.with_untracked(|docs| {
            let doc = docs
                .get(active.get_untracked())
                .ok_or("No active sequence document.")?;
            let source = current_event_binding(doc);
            let pin = crate::ai::web_mcp::Pin {
                workspace_epoch: workspace_epoch.get_untracked(),
                source: source.clone(),
            };
            let (response, navigation) = crate::ai::web_mcp::execute(
                &tool,
                &arguments,
                &pin,
                &doc.molecule,
                selection.get_untracked(),
                selected_feature.get_untracked(),
            )?;
            Ok::<_, String>((response, navigation, source))
        })?;
        if let Some(span) = navigation {
            navigate_primer_binding.run((source, vec![span]))?;
        }
        Ok(response)
    }));
    // Index event identities on workspace changes, not full DNA on each oligo keystroke.
    let primer_event_index = Memo::new(move |_| {
        documents.with(|docs| {
            docs.iter()
                .flat_map(|doc| {
                    doc.summaries()
                        .into_iter()
                        .map(|event| (doc.id.clone(), event.event_id, event.revision))
                        .collect::<Vec<_>>()
                })
                .collect::<std::collections::BTreeSet<_>>()
        })
    });
    let mutation_selection = Callback::new(move |()| {
        if workspace_restoring.get_untracked() {
            return Err("Wait for the workspace operation to finish.".into());
        }
        documents.with_untracked(|docs| {
            let document = active.get_untracked();
            let source = &docs[document];
            let (start, end) = crate::viewer::sequence::selection_toolbar::mutation_range(
                &source.molecule,
                selection.get_untracked(),
                selected_feature.get_untracked(),
                &highlighted_spans.get_untracked(),
            )?;
            Ok(crate::cloning::methods::mutagenesis::MutationSelection {
                source: current_event_binding(source),
                document,
                start,
                end,
            })
        })
    });
    let primer_anchor_available = Callback::new(
        move |binding: collagene_core::material::lineage::EventBinding| {
            primer_event_index.with(|events| {
                events.contains(&(binding.document_id, binding.event_id, binding.revision))
            })
        },
    );

    let second_input = RwSignal::new(0usize);
    let first_input = RwSignal::new(0usize);
    let overlap = RwSignal::new("20".to_string());
    let closure_overlap = RwSignal::new("20".to_string());

    Effect::new(move |_| {
        if let Some(value) = candidate.get() {
            stored_candidates.update(|values| {
                if !values
                    .iter()
                    .any(|old| old.plan_revision == value.plan_revision)
                {
                    values.push(value);
                }
            });
        }
    });
    let plan_method = RwSignal::new("overlap");
    let first_enzymes = RwSignal::new("EcoRI".to_owned());
    let second_enzymes = RwSignal::new("EcoRI".to_owned());
    let first_fragment = RwSignal::new("1".to_owned());
    let second_fragment = RwSignal::new("1".to_owned());
    let first_treatment = RwSignal::new("none".to_owned());
    let second_treatment = RwSignal::new("none".to_owned());
    let first_reverse = RwSignal::new(false);
    let second_reverse = RwSignal::new(false);
    let restriction_policy = RwSignal::new("assumed");
    let restriction_final = RwSignal::new("ligate");
    let second_restriction = RwSignal::new(true);
    let first_repair = RwSignal::new("none".to_owned());
    let second_repair = RwSignal::new("none".to_owned());
    let first_modifications = crate::cloning::methods::checked_cloning::AssumptionSignals::new();
    let second_modifications = crate::cloning::methods::checked_cloning::AssumptionSignals::new();

    let allow_nicked = RwSignal::new(false);

    let capture_junction_plan = Callback::new(
        move |()| -> Result<(DesignPlan, Molecule), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for workspace restoration.".into());
            }
            let (bound_plan, epoch, binding) = last_plan_target
                .get_untracked()
                .ok_or("Reconstruct against the exact intended-target event first.")?;
            if epoch != workspace_epoch.get_untracked() {
                return Err("The workspace changed; reconstruct this exact plan again.".into());
            }
            last_plan.with_untracked(|last| {
            let (plan,target)=last.as_ref().ok_or("Reconstruct the current cloning plan before inspecting junctions.")?;
            collagene_core::cloning::assembly::junctions::source_preflight(plan,target).map_err(|e|e.to_string())?;
            let revision=target.revision();
            if !documents.with_untracked(|docs|crate::protocol::plan_editing_sources::resolve_target(docs,&binding,&revision).is_ok()) {
                return Err("The independently selected target changed; reconstruct the plan again.".into());
            }
            let expected=collagene_core::planning::plan::DesignPlan::revision(plan);
            if expected!=bound_plan {return Err("The reviewed plan changed; reconstruct again.".into());}
            if !plan_report.with_untracked(|report|report.as_ref().is_some_and(|r|r.plan_revision==expected&&r.target_revision==revision&&r.engine_version==verification::ENGINE_VERSION)) {
                return Err("The current reconstruction changed; verify the plan again before inspecting junctions.".into());
            }
            Ok((plan.clone(),target.clone()))
        })
        },
    );
    let quantity_binding = ProtocolController {
        session,
        viewer,
        shell,
    }
    .quantity_binding(
        candidate,
        last_plan,
        last_plan_events,
        last_plan_target,
        stored_candidates,
    );
    let quantity_identity = Memo::new(move |_| {
        let present = last_plan.with(Option::is_some);
        last_plan_target.with(|pin| {
            if present {
                pin.as_ref()
                    .map(|(revision, _, target)| (revision.clone(), target.revision.clone()))
                    .unwrap_or_default()
            } else {
                Default::default()
            }
        })
    });

    Effect::new(move |_| {
        let _ = (
            first_input.get(),
            second_input.get(),
            overlap.get(),
            closure_overlap.get(),
        );
        let _ = (
            plan_method.get(),
            first_enzymes.get(),
            second_enzymes.get(),
            first_fragment.get(),
            second_fragment.get(),
            first_treatment.get(),
            second_treatment.get(),
            first_reverse.get(),
            second_reverse.get(),
            allow_nicked.get(),
            restriction_policy.get(),
            restriction_final.get(),
            second_restriction.get(),
            first_repair.get(),
            second_repair.get(),
        );
        first_modifications.track();
        second_modifications.track();
        plan_state.invalidate_evidence();
    });
    let apply_restored = Callback::new(move |state: RestoredWorkspace| {
        // Invalidate every outstanding source receipt before any owner installs data.
        batch(|| {
            session.begin_replacement();
            pcr_state.workspace_replaced();
            plan_state.restore(state.candidate, state.candidates);
            import_state.restore(state.oligos);
            comparison_state.restore(state.alignments);
            annotation_state.restore(state.libraries);
            read_state.clear();
            read_archive_state.restore(state.assemblies);
            primer_design_state.restore(state.pair_draft);
            run_history_state.restore(state.runs, state.hidden_products);
            reference_transfer_state.restore();
            assembly_state.restore(state.assembly_draft);
            digest_workspace_state.restore();
            session.install_documents(state.name, state.documents, state.active, state.attachments);
            viewer.clear();
            shell.workspace_replaced();
        });
    });
    // Persistence callbacks outlive the dock/modal that initiated a transaction.
    // A committed write must reach App even after that child owner is disposed.
    let validate_workspace_frame =
        Callback::new(|frame: collagene_core::documents::workspace::Workspace| {
            RestoredWorkspace::from_frame(&frame).map(|_| ())
        });
    let restore_workspace_frame = Callback::new(
        move |frame: collagene_core::documents::workspace::Workspace| {
            match RestoredWorkspace::from_frame(&frame) {
                Ok(state) => batch(|| apply_restored.run(state)),
                Err(error) => message.set(error),
            }
        },
    );
    let restore_library_frame = Callback::new(
        move |frame: collagene_core::documents::workspace::Workspace| {
            match RestoredWorkspace::from_frame(&frame) {
                Ok(state) => {
                    let visible = assembly_visible.get_untracked();
                    batch(|| {
                        apply_restored.run(state);
                        assembly_visible.set(visible);
                    });
                }
                Err(error) => message.set(error),
            }
        },
    );
    let library_current = Memo::new(move |_| {
        let epoch = workspace_epoch.get();
        documents.with(|docs| {
            let document = docs.get(active.get())?;
            Some(crate::library::library_view::LocalDocument {
                id: document.id.clone(),
                name: document.molecule.name.clone(),
                revision: collagene_library_protocol::Digest::new(document.head_revision()).ok()?,
                head_event: Some(document.head_event().into()),
                workspace_epoch: epoch,
            })
        })
    });
    let library_bindings = Memo::new(move |_| {
        attachments
            .with(|assets| crate::library::library_documents::bindings(assets).unwrap_or_default())
    });
    let native_catalog_refresh = RwSignal::new(0u64);
    // Native activation status must be checked before enabling saved-library
    // writes or synchronization. Working-editor navigation/autosave stays live.
    let native_catalog_pending = RwSignal::new(crate::platform::desktop_rpc::is_desktop());
    let native_catalog_locked =
        Memo::new(move |_| workspace_restoring.get() || native_catalog_pending.get());
    let native_destination = RwSignal::new(None::<collagene_library_protocol::FolderId>);
    let native_catalog_epoch = RwSignal::new(0u64);
    let catalog_pending_changed = LibraryController {
        session,
        viewer,
        shell,
    }
    .catalog_pending_changed(
        native_catalog_epoch,
        native_catalog_pending,
        native_catalog_refresh,
        native_destination,
        native_library_context,
    );
    let native_destination_memo = Memo::new(move |_| native_destination.get());
    let native_folder_changed = Callback::new(move |folder| native_destination.set(folder));
    let native_catalog_changed = Callback::new(move |()| {
        native_catalog_refresh.update(|generation| *generation = generation.wrapping_add(1));
    });
    let native_saved = Memo::new(move |_| {
        let _ = loaded.get();
        let _ = native_catalog_refresh.get();
        let context = native_library_context.get()?;
        let current = library_current.get()?;
        let remote = crate::library::library_bindings::Remote {
            server: "collagene-desktop://local".into(),
            lab: context.lab_id,
        };
        let binding =
            library_bindings.with(|bindings| bindings.get(&remote, &current.id).cloned())?;
        Some(crate::library::remote_library_view::SavedLocalDocument {
            id: binding.remote_document,
            revision: binding.remote_revision,
            name: current.name,
            dirty: binding.local_revision != current.revision
                || binding.local_head_event != current.head_event,
        })
    });
    let library_snapshot = Callback::new(move |()| {
        documents.with_untracked(|docs| {
            WorkspaceSession::snapshot(
                docs,
                active.get_untracked(),
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &attachments.get_untracked(),
                &workspace_name.get_untracked(),
            )
        })
    });
    let library_capture = Callback::new(
        move |request: crate::library::library_documents::CaptureRequest| -> crate::library::library_view::LocalResult<
            crate::library::library_documents::CaptureOutcome,
        > {
            let prepared = (|| {
                if workspace_restoring.get_untracked() || !loaded.get_untracked() {
                    return Err("Wait for local workspace storage to finish opening.".to_string());
                }
                Ok((library_snapshot.run(())?, workspace_epoch.get_untracked()))
            })();
            Box::pin(async move {
                let (frame, epoch) = prepared?;
                let captured = if crate::platform::desktop_rpc::is_desktop() {
                    crate::library::library_documents::capture_private_local(frame, epoch, request)
                        .await?
                } else {
                    crate::library::library_documents::capture(frame, epoch, request).await?
                };
                if workspace_epoch.try_get_untracked() != Some(epoch)
                    || workspace_restoring.try_get_untracked().unwrap_or(true)
                {
                    return Err(
                        "Workspace changed while preparing the shared plasmid; nothing was sent."
                            .into(),
                    );
                }
                Ok(captured)
            })
        },
    );
    let library_record_binding = LibraryController {
        session,
        viewer,
        shell,
    }
    .library_record_binding(library_bindings);
    let library_open = LibraryController {
        session,
        viewer,
        shell,
    }
    .library_open(
        apply_restored,
        assembly_visible,
        db,
        library_bindings,
        library_snapshot,
        loaded,
        save_generation,
    );
    autosave.observe(
        assembly_draft,
        pair_draft,
        stored_candidates,
        candidate,
        apply_restored,
    );
    Effect::new(move |_| {
        if !has_open_document.get() {
            if let Some(map) = renderer.get_value() {
                map.pointer_cancel();
            }
            return;
        }
        molecule.with(|current| {
            let selected = highlighted_spans.get();
            let selected_index = selected_feature.get();
            if let (Some(canvas), Some(overlay)) = (gl_canvas.get(), label_canvas.get()) {
                if let Some(map) = renderer.get_value() {
                    map.set_document_spans(
                        current,
                        &molecule_revision.get(),
                        &selected,
                        selected_index,
                    );
                    map.set_hidden_features(&annotation_hidden.get());
                } else {
                    match MapRenderer::new(
                        canvas.unchecked_ref::<HtmlCanvasElement>().clone(),
                        overlay.unchecked_ref::<HtmlCanvasElement>().clone(),
                        current.clone(),
                        std::rc::Rc::new(move |stats| {
                            if map_labels
                                .try_get_untracked()
                                .is_some_and(|current| current != stats)
                            {
                                let _ = map_labels.try_set(stats);
                            }
                        }),
                    ) {
                        Ok(map) => {
                            map.set_restriction_hover_callback(std::rc::Rc::new(move |site| {
                                if hovered_restriction.try_get_untracked() != Some(site) {
                                    let _ = hovered_restriction.try_set(site);
                                }
                            }));
                            map.set_selection_mode(false);
                            map.set_hidden_features(&annotation_hidden.get());
                            renderer.set_value(Some(map));
                        }
                        Err(error) => render_error.set(error),
                    }
                }
            }
        });
    });
    Effect::new(move |_| {
        let _ = restriction_identity.get();
        hovered_restriction.set(None);
    });
    Effect::new(move |_| {
        let context = restriction_identity.get();
        let sites = restrictions
            .displayed
            .get()
            .unwrap_or_else(|| std::sync::Arc::new(Vec::new()));
        let selected = hovered_restriction.get().or(restrictions.selected.get());
        let revision = molecule_revision.get();
        if let Some(map) = renderer.get_value() {
            map.set_interaction_context(&context);
            map.set_restriction_sites(&revision, sites, selected);
            map.set_pcr_overlay(pcr_overlay.get());
        }
    });
    Effect::new(move |_| {
        // All document writes still invalidate multi-source scientific reports.
        // Selection follows only the exact active document/event and restore epoch.
        let epoch = workspace_epoch.get();
        let index = active.get();
        let pin = documents.with(|docs| {
            docs.get(index).map(|doc| {
                crate::viewer::sequence::sequence_edit_transition::SourcePin::of(epoch, doc)
            })
        });
        let mut decision = crate::viewer::sequence::sequence_edit_transition::Decision::Unchanged;
        sequence_transitions.update_value(|tracker| decision = tracker.observe(pin));
        match decision {
            crate::viewer::sequence::sequence_edit_transition::Decision::Unchanged => {}
            crate::viewer::sequence::sequence_edit_transition::Decision::Clear => {
                selected_feature.set(None);
                selection.set(None);
                search_regions.set(None);
            }
            crate::viewer::sequence::sequence_edit_transition::Decision::PreserveCaret(caret) => {
                // Read untracked: caret changes must not retrigger source processing.
                // A later pointer selection takes precedence over the local edit.
                if selection.get_untracked() == Some((caret, caret)) {
                    selected_feature.set(None);
                    search_regions.set(None);
                    selection.set(Some((caret, caret)));
                }
            }
        }
        plan_state.invalidate_evidence();
    });

    let open_file = ImportController {
        session,
        viewer,
        shell,
    }
    .open_file(
        annotation_format,
        annotation_source,
        candidate,
        candidate_events,
        existing_oligo,
        fastq_encoding,
        last_plan,
        last_plan_target,
        pending_alignment,
        pending_compat_import,
        pending_fastq,
        pending_oligo,
        pending_read_alignment,
        pending_sbol,
        pending_variants,
        plan_generation,
        plan_report,
    );
    let save_file = move |_: web_sys::MouseEvent| {
        let mol = molecule.get_untracked();
        match protocol::save_document(&mol) {
            Ok(bytes) => crate::io::download::download_with_completion(
                &format!("{}.cg2", mol.name),
                "application/octet-stream",
                &bytes,
                move |result| {
                    if message.try_get_untracked().is_none() {
                        return;
                    }
                    if let Err(e) = result {
                        message.set(e)
                    }
                },
            ),
            Err(e) => message.set(e.to_string()),
        }
    };
    let do_pcr = PcrController {
        session,
        viewer,
        shell,
        state: pcr_state,
    }
    .do_pcr(library_snapshot);
    let pcr_review_current =
        Callback::new(move |review: crate::primers::pcr::pcr_outputs::Review| {
            !workspace_restoring.get_untracked()
                && workspace_epoch.get_untracked() == review.epoch
                && pcr_generation.get_untracked() == review.generation
                && pcr_origin.get_untracked() == review.origin
                && review.inputs
                    == (
                        pcr_forward_name.get_untracked(),
                        pcr_reverse_name.get_untracked(),
                        pcr_forward.get_untracked(),
                        pcr_reverse.get_untracked(),
                        pcr_forward_tail.get_untracked(),
                        pcr_reverse_tail.get_untracked(),
                    )
                && documents.with_untracked(|docs| {
                    docs.get(active.get_untracked())
                        .is_some_and(|doc| current_event_binding(doc) == review.source)
                })
                && pcr_review.with_untracked(|current| {
                    current.as_ref().is_some_and(|current| {
                        current.generation == review.generation
                            && std::sync::Arc::ptr_eq(&current.products, &review.products)
                    })
                })
        });
    let open_material_derivation = Callback::new(
        move |request: crate::restriction::digest_product::MaterialOpenRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || workspace_epoch.get_untracked() != request.epoch
            {
                return Err("The material source workspace changed; review again.".into());
            }
            request.guard.run(())?;
            let original = request.guard;
            let bound_epoch = request.epoch;
            let guarded = Callback::new(move |()| -> Result<(), String> {
                if workspace_restoring.get_untracked()
                    || workspace_epoch.get_untracked() != bound_epoch
                    || modal.get_untracked() != "material_product"
                {
                    return Err("The material installation context changed.".into());
                }
                original.run(())
            });
            pending_material.set(Some(crate::restriction::digest_product::MaterialOpenRequest {
                guard: guarded,
                ..request
            }));
            modal_error.set(String::new());
            modal.set("material_product");
            Ok(())
        },
    );
    let create_pcr_outputs = PcrController {
        session,
        viewer,
        shell,
        state: pcr_state,
    }
    .create_pcr_outputs(
        candidate,
        library_snapshot,
        oligo_tables,
        pcr_review_current,
        pending_material,
        stored_candidates,
    );
    let run_plan = ProtocolController {
        session,
        viewer,
        shell,
    }
    .run_plan(plan_state);
    let verify_candidate = AssemblyController {
        session,
        viewer,
        shell,
    }
    .verify_candidate(
        candidate,
        candidate_events,
        last_plan,
        last_plan_target,
        plan_report,
        run_plan,
    );
    let export_workspace = move |_: web_sys::MouseEvent| {
        let mut unique = std::collections::BTreeMap::new();
        documents.with_untracked(|docs| {
            for d in docs {
                unique
                    .entry(d.molecule.revision())
                    .or_insert_with(|| d.molecule.clone());
            }
        });
        let snapshot = transfer::PlanningWorkspace {
            target_revision: molecule.get_untracked().revision(),
            documents: unique.into_values().collect(),
        };
        match transfer::save_workspace(&snapshot) {
            Ok(bytes) if bytes.len()<=8*1024*1024=>{crate::io::download::download_with_completion("planning-workspace.cg2design","application/octet-stream",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})},
            Ok(_)=>message.set("The initial MCP exchange supports snapshots up to 8 MiB. Open only the target and needed inputs, then export again.".into()),
            Err(e)=>message.set(e.to_string()),
        }
    };
    let do_plan = AssemblyController {
        session,
        viewer,
        shell,
    }
    .do_plan(
        allow_nicked,
        closure_overlap,
        first_enzymes,
        first_fragment,
        first_input,
        first_modifications,
        first_repair,
        first_reverse,
        first_treatment,
        overlap,
        plan_method,
        restriction_final,
        restriction_policy,
        run_plan,
        second_enzymes,
        second_fragment,
        second_input,
        second_modifications,
        second_repair,
        second_restriction,
        second_reverse,
        second_treatment,
    );

    let use_digest_in_assembly = Callback::new(
        move |fragment: crate::restriction::digest_fragment_model::AssemblyFragment| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || workspace_epoch.get_untracked() != fragment.epoch
            {
                return Err("The workspace changed; run the digest again.".into());
            }
            documents.with_untracked(|docs| -> Result<(), String> {
                let doc = docs
                    .iter()
                    .find(|d| current_event_binding(d) == fragment.source.binding)
                    .ok_or("The original digest source event changed; run the digest again.")?;
                if doc.molecule.revision() != fragment.molecule.revision() {
                    return Err("The digest source snapshot changed.".into());
                }
                Ok(())
            })?;
            let mut staged = assembly_draft.get_untracked();
            staged.set_digest_fragment(fragment)?;
            assembly_draft.set(staged);
            assembly_visible.set(true);
            message.set("Selected digest fragment placed in the active assembly slot. Its original source and pinned digest preparation are retained.".into());
            Ok(())
        },
    );
    let finalize_assembly = AssemblyController {
        session,
        viewer,
        shell,
    }
    .finalize_assembly(candidate, pending_finalize, stored_candidates);
    let finalize_reviewed_candidate = AssemblyController {
        session,
        viewer,
        shell,
    }
    .finalize_reviewed_candidate(
        candidate,
        capture_junction_plan,
        last_plan_events,
        last_plan_target,
        pending_finalize,
        stored_candidates,
    );
    let reviewed_product_current = Memo::new(move |_| {
        let _ = (workspace_epoch.get(), workspace_restoring.get());
        let Some((revision, epoch, binding)) = last_plan_target.get() else {
            return false;
        };
        epoch == workspace_epoch.get()
            && !workspace_restoring.get()
            && documents.with(|docs| {
                crate::protocol::plan_editing_sources::resolve_target(
                    docs,
                    &binding,
                    &binding.revision,
                )
                .is_ok()
            })
            && plan_report.with(|report| {
                report.as_ref().is_some_and(|report| {
                    report.plan_revision == revision
                        && report.target_revision == binding.revision
                        && report.engine_version == verification::ENGINE_VERSION
                })
            })
            && last_plan.with(Option::is_some)
    });
    let reviewed_product_guard = Callback::new(move |()| capture_junction_plan.run(()).map(|_| ()));
    let finalize_snapshot = Callback::new(
        move |()| -> Result<collagene_core::documents::workspace::Workspace, String> {
            let pending = pending_finalize
                .get_untracked()
                .ok_or("No assembly output choice is pending.")?;
            let current = documents.with_untracked(|docs| {
                WorkspaceSession::snapshot(
                    docs,
                    active.get_untracked(),
                    candidate.get_untracked().as_ref(),
                    &stored_candidates.get_untracked(),
                    &attachments.get_untracked(),
                    &workspace_name.get_untracked(),
                )
            })?;
            if collagene_core::material::lineage::workspace_fingerprint(&current)
                .map_err(|e| e.to_string())?
                != pending.source_fingerprint
            {
                return Err(
                    "The workspace changed while choosing outputs. Review this assembly again."
                        .into(),
                );
            }
            Ok(pending.workspace.as_ref().clone())
        },
    );
    let review_assembly = AssemblyController {
        session,
        viewer,
        shell,
    }
    .review_assembly(
        assembly_visible,
        candidate,
        candidate_events,
        last_plan,
        last_plan_events,
        last_plan_target,
        plan_generation,
        plan_report,
        stored_candidates,
    );
    let navigate_assembly = AssemblyController {
        session,
        viewer,
        shell,
    }
    .navigate_assembly(candidate, stored_candidates);

    let variant_asset =
        move || -> Result<(String, collagene_core::documents::workspace::WorkspaceAttachment), String> {
            let Some((source, format, name)) = pending_variants.get_untracked() else {
                return Err("No variant source is open.".into());
            };
            collagene_core::review::variants::archive::make_attachment(&source, format, &name)
                .map_err(|e| e.to_string())
        };
    let keep_variant_source = Callback::new(move |()| -> Result<(), String> {
        if workspace_restoring.get_untracked() {
            return Err("Wait for workspace restoration.".into());
        }
        let (hash, asset) = variant_asset()?;
        let mut assets = attachments.get_untracked();
        if let Some(existing) = assets.get(&hash) {
            if existing.kind != asset.kind {
                return Err(
                    "These original bytes are already retained under a different source format."
                        .into(),
                );
            }
            return Ok(());
        }
        assets.insert(hash, asset);
        let frame = documents.with_untracked(|docs| {
            WorkspaceSession::snapshot(
                docs,
                active.get_untracked(),
                candidate.get_untracked().as_ref(),
                &stored_candidates.get_untracked(),
                &assets,
                &workspace_name.get_untracked(),
            )
        })?;
        frame.validate().map_err(|e| e.to_string())?;
        attachments.set(assets);
        Ok(())
    });
    let open_variant_preview = ImportController {
        session,
        viewer,
        shell,
    }
    .open_variant_preview(candidate, stored_candidates, variant_asset);
    let keep_read_source = ImportController {
        session,
        viewer,
        shell,
    }
    .keep_read_source(candidate, pending_read_alignment, stored_candidates);
    let keep_alignment = ImportController {
        session,
        viewer,
        shell,
    }
    .keep_alignment(alignments, candidate, stored_candidates);
    let keep_library = LibraryController {
        session,
        viewer,
        shell,
    }
    .keep_library(candidate, libraries, stored_candidates);
    let reference_copy_selection = ReferenceController {
        session,
        viewer,
        shell,
    }
    .reference_copy_selection();
    let biology_selection_location = Callback::new(move |()| {
        molecule.with_untracked(|m| {
            let source = documents
                .with_untracked(|docs| current_event_binding(&docs[active.get_untracked()]));
            crate::biology::biology_selection::capture(
                m,
                &source,
                selected_feature.get_untracked(),
                selection.get_untracked(),
                search_regions.get_untracked().as_ref(),
                biology_directed.get_untracked().as_ref(),
            )
        })
    });
    let select_biology_location = Callback::new(
        move |(revision, location): (String, Location)| -> Result<(), String> {
            if workspace_restoring.get_untracked() || molecule_revision.get_untracked() != revision
            {
                return Err("The ORF source changed; rerun the analysis.".into());
            }
            let source = documents
                .with_untracked(|docs| current_event_binding(&docs[active.get_untracked()]));
            let directed = molecule.with_untracked(|m| {
                crate::biology::biology_selection::DirectedSelection::new(source, m, location)
            })?;
            batch(|| {
                selected_feature.set(None);
                selection.set(directed.contiguous());
                search_regions.set(Some((revision, directed.spans.clone())));
                biology_directed.set(Some(directed));
            });
            Ok(())
        },
    );
    Effect::new(move |_| {
        let _ = workspace_epoch.get();
        biology_directed.set(None);
        pending_translation_source.set(None);
    });
    let reference_paste_guard = ReferenceController {
        session,
        viewer,
        shell,
    }
    .reference_paste_guard(internal_clipboard);
    let preview_reference_paste = Callback::new(
        move |opened: crate::workspace::reference_transfer_view::PasteOpen| -> Result<(), String> {
            reference_paste_guard.run(opened.clone())?;
            pending_reference_paste.set(Some(opened));
            modal_error.set(String::new());
            modal.set("reference_paste");
            Ok(())
        },
    );
    let apply_reference_paste = ReferenceController {
        session,
        viewer,
        shell,
    }
    .apply_reference_paste(
        candidate,
        pending_reference_paste,
        reference_paste_guard,
        stored_candidates,
    );
    let commit_annotations = Callback::new(move |(base, edited): (String, Molecule)| {
        session.commit_document(&base, edited)?;
        selected_feature.set(None);
        selection.set(None);
        Ok::<_, String>(())
    });
    let import_oligos = ImportController {
        session,
        viewer,
        shell,
    }
    .import_oligos(candidate, oligo_tables, stored_candidates);
    let open_mismatch_outputs = PcrController {
        session,
        viewer,
        shell,
        state: pcr_state,
    }
    .open_mismatch_outputs(candidate, oligo_tables, stored_candidates);
    let restore_revision = Callback::new(move |event: String| {
        if workspace_restoring.get_untracked() {
            return Err("Wait for workspace restoration to finish.".into());
        }
        let index = active.get_untracked();
        let mut docs = documents.get_untracked();
        docs.get_mut(index)
            .ok_or("No active document")?
            .checkout(&event)
            .map_err(|e| e.to_string())?;
        if docs
            .iter()
            .map(|d| document_bytes(&d.molecule))
            .sum::<usize>()
            > MAX_WORKSPACE_BYTES
        {
            return Err("Restoring this revision exceeds the current document size limit; open it in a smaller workspace.".into());
        }
        documents.set(docs);
        Ok(())
    });
    let fork_revision = DocumentController {
        session,
        viewer,
        shell,
    }
    .fork_revision(candidate, stored_candidates);

    let contig_snapshot = Callback::new(move |()| {
        documents.with_untracked(|docs| {
            stored_candidates.with_untracked(|plans| {
                attachments.with_untracked(|assets| {
                    candidate.with_untracked(|current| {
                        WorkspaceSession::snapshot(
                            docs,
                            active.get_untracked(),
                            current.as_ref(),
                            plans,
                            assets,
                            &workspace_name.get_untracked(),
                        )
                    })
                })
            })
        })
    });
    let navigate_pcr_overlay = PcrController {
        session,
        viewer,
        shell,
        state: pcr_state,
    }
    .navigate_pcr_overlay(pending_pcr_source);
    let open_translated_protein = BiologyController {
        session,
        viewer,
        shell,
    }
    .open_translated_protein(candidate, stored_candidates);
    let retain_translation_cds = BiologyController {
        session,
        viewer,
        shell,
    }
    .retain_translation_cds();
    let open_translation_origin = BiologyController {
        session,
        viewer,
        shell,
    }
    .open_translation_origin(
        contig_snapshot,
        pending_translation_source,
        restriction_selection,
    );
    let set_annotation_visibility = AnnotationController {
        session,
        viewer,
        shell,
    }
    .set_annotation_visibility(contig_snapshot);
    let apply_annotation_batch = AnnotationController {
        session,
        viewer,
        shell,
    }
    .apply_annotation_batch(candidate, stored_candidates);
    let restriction_controller = RestrictionController {
        session,
        viewer,
        shell,
    };
    let saved_enzyme_sets = restriction_controller.saved_enzyme_sets(candidate, stored_candidates);
    let save_enzyme_set = restriction_controller.save_enzyme_set(contig_snapshot);
    let delete_enzyme_set = restriction_controller.delete_enzyme_set(contig_snapshot);
    let resolve_retained_source = Callback::new(
        move |event: collagene_core::material::lineage::EventBinding| -> Result<Molecule, String> {
            documents.with_untracked(|docs| {
                let doc = docs
                    .iter()
                    .find(|doc| doc.id == event.document_id)
                    .ok_or("The retained source document is unavailable.")?;
                let molecule = doc
                    .revision_molecule(&event.event_id)
                    .map_err(|e| e.to_string())?;
                if molecule.revision() != event.revision {
                    return Err(
                        "The exact retained source event differs from its saved revision.".into(),
                    );
                }
                Ok(molecule)
            })
        },
    );
    let saved_digests = restriction_controller.saved_digests(candidate, stored_candidates);
    let save_digest_setup = RestrictionController {
        session,
        viewer,
        shell,
    }
    .save_digest_setup(contig_snapshot);
    let load_digest_setup = restriction_controller.load_digest_setup(contig_snapshot);
    let delete_digest_setup = restriction_controller.delete_digest_setup(contig_snapshot);
    let queue_digest_gel = restriction_controller.queue_digest_gel(
        contig_snapshot,
        resolve_retained_source,
        gel_digest_seeds,
    );
    let fork_digest_source =
        restriction_controller.fork_digest_source(candidate, stored_candidates);
    let biology_controller = BiologyController {
        session,
        viewer,
        shell,
    };
    let (antibody_overlay, publish_antibody) = biology_controller.antibody_publication();
    let prepare_antibody = biology_controller.prepare_antibody();
    let select_antibody_residue = biology_controller.select_antibody_residue(biology_directed);
    let prepare_protein_liabilities = biology_controller.prepare_protein_liabilities();
    let select_protein_liability = BiologyController {
        session,
        viewer,
        shell,
    }
    .select_protein_liability(biology_directed);
    let annotate_protein_liability = BiologyController {
        session,
        viewer,
        shell,
    }
    .annotate_protein_liability(candidate);
    let saved_rna = biology_controller.saved_rna(candidate, stored_candidates);
    let resolve_rna = biology_controller.resolve_rna();
    let save_rna = biology_controller.save_rna(contig_snapshot);
    let load_rna = biology_controller.load_rna(contig_snapshot);
    let delete_rna = biology_controller.delete_rna(contig_snapshot);
    let comparison_controller = ComparisonController {
        session,
        viewer,
        shell,
    };
    let saved_comparisons = comparison_controller.saved_comparisons(candidate, stored_candidates);
    let save_comparison = comparison_controller.save_comparison(contig_snapshot);
    let reopen_comparison = comparison_controller.reopen_comparison(contig_snapshot);
    let prepare_similarity = ComparisonController {
        session,
        viewer,
        shell,
    }
    .prepare_similarity();
    let guard_similarity = comparison_controller.guard_similarity();
    let navigate_comparison = ComparisonController {
        session,
        viewer,
        shell,
    }
    .navigate_comparison(pending_comparison_source);
    let open_run_event = Callback::new(
        move |request: crate::history::run_event_model::Request| -> Result<(), String> {
            if workspace_restoring.get_untracked()
                || request.epoch != workspace_epoch.get_untracked()
            {
                return Err("The workspace changed; reopen the cloning run.".into());
            }
            let frame = contig_snapshot.run(())?;
            let resolved = crate::history::run_event_model::resolve(&frame, &request)?;
            pending_run_event.set(Some(std::sync::Arc::new(resolved)));
            modal_error.set(String::new());
            modal.set("run_event");
            Ok(())
        },
    );
    let fork_run_event = ReadController {
        session,
        viewer,
        shell,
    }
    .fork_run_event(
        contig_snapshot,
        last_plan,
        last_plan_events,
        last_plan_target,
        pending_run_event,
        plan_report,
    );
    let back_to_run = Callback::new(move |run_id: String| {
        if workspace_restoring.get_untracked() {
            return;
        }
        selected_run.set(run_id);
        pending_run_event.set(None);
        modal_error.set(String::new());
        modal.set("run_details");
    });
    let record_read_check = ReadController {
        session,
        viewer,
        shell,
    }
    .record_read_check(pending_read_check, read_check_reference);
    let ReadReviewActions {
        read_check_snapshot,
        restore_read_check,
        open_read_review,
        reopen_read_review,
        save_read_review,
        read_review_snapshot,
        restore_read_review,
    } = ReadController {
        session,
        viewer,
        shell,
    }
    .review_actions(
        read_state,
        contig_snapshot,
        read_check_reference,
        read_check_attachment_available,
        apply_restored,
    );

    let navigate_check_reference = ReadController {
        session,
        viewer,
        shell,
    }
    .navigate_check_reference(candidate, stored_candidates);

    let open_mutation_intent = AssemblyController {
        session,
        viewer,
        shell,
    }
    .open_mutation_intent(contig_snapshot, last_plan, last_plan_target, plan_report);
    let open_domestication = AssemblyController {
        session,
        viewer,
        shell,
    }
    .open_domestication(
        candidate,
        candidate_events,
        contig_snapshot,
        last_plan,
        last_plan_events,
        last_plan_target,
        plan_generation,
        plan_report,
        stored_candidates,
    );
    let contig_restore = Callback::new(
        move |frame: collagene_core::documents::workspace::Workspace| {
            match RestoredWorkspace::from_frame(&frame) {
                Ok(state) => batch(|| apply_restored.run(state)),
                Err(error) => message.set(error),
            }
        },
    );
    let record_assembly = Callback::new(
        move |request: crate::reads::contigs::ArchiveRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            if workspace_epoch.get_untracked() != request.epoch {
                return Err("The workspace changed; reopen the assembly.".into());
            }
            pending_assembly.set(Some(request));
            modal_error.set(String::new());
            modal.set("contig_archive");
            Ok(())
        },
    );
    let record_donor = Callback::new(
        move |request: crate::cloning::methods::crispr::DonorRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the native workspace update to finish.".into());
            }
            if request.epoch != workspace_epoch.get_untracked()
                || !documents.with_untracked(|docs| {
                    [&request.source, &request.target].into_iter().all(|info| {
                        docs.iter().any(|d| {
                            d.id == info.binding.document_id
                                && d.head_event() == info.binding.event_id
                                && d.molecule.revision() == info.binding.revision
                        })
                    })
                })
            {
                return Err(
                    "Source or independent target event changed. Design and check again.".into(),
                );
            }
            pending_donor.set(Some(request));
            modal_error.set(String::new());
            modal.set("crispr_product");
            Ok(())
        },
    );
    let open_digest = Callback::new(
        move |request: crate::restriction::digest_product::DigestOpenRequest| -> Result<(), String> {
            if workspace_restoring.get_untracked() {
                return Err("Wait for the workspace update to finish.".into());
            }
            if request.epoch != workspace_epoch.get_untracked()
                || !documents.with_untracked(|docs| {
                    docs.get(active.get_untracked()).is_some_and(|d| {
                        d.id == request.source.document_id
                            && d.head_event() == request.source.event_id
                            && molecule_revision.get_untracked() == request.source.revision
                    })
                })
            {
                return Err("The source document/event changed. Run the digest again before opening a fragment.".into());
            }
            pending_digest.set(Some(request));
            modal_error.set(String::new());
            modal.set("digest_product");
            Ok(())
        },
    );
    let navigate_digest = Callback::new(move |request: crate::restriction::digest::Navigation| {
        if workspace_restoring.get_untracked()
            || request.epoch != workspace_epoch.get_untracked()
            || !documents.with_untracked(|docs| {
                docs.get(active.get_untracked()).is_some_and(|d| {
                    d.id == request.source.document_id
                        && d.head_event() == request.source.event_id
                        && molecule_revision.get_untracked() == request.source.revision
                })
            })
        {
            return;
        }
        let n = molecule.with_untracked(|m| m.sequence.len());
        if request.spans.is_empty() || request.spans.iter().any(|(a, b)| a >= b || *b > n) {
            return;
        }
        batch(|| {
            selected_feature.set(None);
            selection.set(if request.spans.len() == 1 {
                request.spans.first().copied()
            } else {
                None
            });
            search_regions.set(Some((request.source.revision, request.spans)));
        });
    });
    let navigate_read_call = ReadController {
        session,
        viewer,
        shell,
    }
    .navigate_read_call(candidate, stored_candidates);

    view! {
                        <Show when=move||workspace_restoring.get()&&!matches!(modal.get(),"assembly_finalize"|"read_check_capture"|"read_review_commit")><div role="status" aria-live="polite" style="position:fixed;inset:0;z-index:10000;background:#f5f8f4ed;display:grid;place-items:center;pointer-events:all"><p>"Updating workspace…"</p></div></Show>
                        <div class="application workspace-shell" on:click=dismiss_viewer_menus><div class="workspace-contents" inert=move||workspace_restoring.get()>
                            <header class="topbar">
                                <a class="brand" href="#" aria-label="Collagene workspace"><svg viewBox="0 0 32 32" aria-hidden="true"><circle cx="16" cy="16" r="11" fill="none" stroke="currentColor" stroke-width="3"/><path d="M9 7L23 25M7 10L11 6M21 26L25 22" stroke="currentColor" stroke-width="3"/></svg><span>"collagene"</span><span class="version">"2"</span></a>
                                <div class="top-divider"></div><span class="workspace-label">{move||workspace_name.get()}</span>
                                <div class="top-actions"><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("workspace");}>"Workspace"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("import");}>"Paste sequence"</button><button class="button primary" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open files"</button></div>
                                <input node_ref=file_input type="file" multiple accept=".gb,.gbk,.gbff,.genbank,.gp,.gpt,.genpept,.gcg,.seq,.pep,.uniprot,.swiss,.swissprot,.sbol,.jsonld,.json,.ape,.fa,.fasta,.fna,.faa,.frn,.txt,.cg2,.cg2design,.cg2plan,.dna,.prot,.embl,.emb,.ab1,.abi,.scf,.fastq,.fq,.csv,.tsv,.gff3,.gff,.gtf,.bed,.xml,.collagene,.seqfile,.xdna,.xdgn,.xprt,.aln,.clustal,.sto,.stockholm,.afa,.afasta,.cg2aln,.phy,.phylip,.msf,.nex,.nexus,.pir,.nbrf,.gde,.slx,.selex,.sel,.sam,.bam,.vcf,.bcf,.gz" style="display:none" on:change=open_file/>
                            </header>
                            <div class="body-layout" class:library-hidden=move||!show_library.get()>
                                <nav class="workspace-rail" aria-label="Workspace navigation">
                                    <button class="rail-button" class:active=move||show_library.get() title="Documents" aria-label="Documents" aria-expanded=move||show_library.get().to_string() disabled=move||!layout.ready.get() on:click=move|_|show_library.update(|visible|*visible=!*visible)>
                                        <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M3 5h7l2 2h9v13H3z"/><path d="M3 10h18"/></svg>
                                    </button>
                                    <button class="rail-button" title="Import sequence files" aria-label="Import sequence files" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>
                                        <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M12 3v12m-4-4 4 4 4-4M4 16v5h16v-5"/></svg>
                                    </button>
                                    <button class="rail-button" title="Import from Addgene or NCBI" aria-label="Import from Addgene or NCBI" on:click=move|_|{modal_error.set(String::new());modal.set("url");}>
                                        <svg viewBox="0 0 24 24" aria-hidden="true"><path d="m10 13 4-4m-6 6-1 1a4 4 0 0 1-6-6l4-4a4 4 0 0 1 6 0m2 3 1-1a4 4 0 0 1 6 6l-4 4a4 4 0 0 1-6 0"/></svg>
                                    </button>
                                    <span class="rail-spacer"></span>
                                    <button class="rail-button" title="Workspace files and settings" aria-label="Workspace files and settings" on:click=move|_|{modal_error.set(String::new());modal.set("workspace");}>
                                        <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M4 6h16M4 12h16M4 18h16"/><circle cx="9" cy="6" r="2"/><circle cx="15" cy="12" r="2"/><circle cx="8" cy="18" r="2"/></svg>
                                    </button>
                                </nav>
                                <div class="explorer-frame" id="document-explorer-pane"><Show when=crate::library::library_view::enabled><crate::library::library_view::LibraryPane client=crate::library::library_client::browser() current=library_current bindings=library_bindings capture=library_capture on_open=library_open on_binding=library_record_binding locked=native_catalog_locked refresh=native_catalog_refresh on_folder=native_folder_changed/></Show><Show when=crate::platform::desktop_rpc::is_desktop><crate::library::catalog_backup_view::CatalogBackupControls client=Callback::new(|command|->crate::library::catalog_backup_view::LocalResult<collagene_desktop_protocol::CatalogBackupSuccess>{Box::pin(async move{match crate::platform::desktop_rpc::call(collagene_desktop_protocol::Command::CatalogBackup(command)).await?{collagene_desktop_protocol::Success::CatalogBackup(value)=>Ok(value),_=>Err("Unexpected native catalog backup reply.".into())}})}) context=Memo::new(move|_|native_library_context.get()) locked=explorer_locked on_catalog_pending=catalog_pending_changed/><crate::library::remote_library_view::RemoteLibraryPane client=crate::platform::desktop_rpc::sync_client() native=true saved=native_saved destination=native_destination_memo locked=native_catalog_locked on_local_change=native_catalog_changed/></Show><crate::workspace::document_explorer::DocumentExplorer documents=explorer_documents active_document=explorer_active locked=explorer_locked on_activate=explorer_activate on_open=Callback::new(move|()|if let Some(input)=file_input.get(){input.click();})/><Show when=move||!hidden_products.with(|ids|ids.is_empty())><label class="preparation-visibility"><input type="checkbox" prop:checked=move||show_prepared_history.get() on:change=move|e|show_prepared_history.set(event_target::<HtmlInputElement>(&e).checked())/>{move||format!("Show preparation history ({})",hidden_products.with(|ids|ids.len()))}</label></Show><crate::workspace::workspace_layout::PaneResize layout=layout pane=crate::workspace::workspace_layout::Pane::Documents/></div>
                                <main>
                                    <crate::workspace::document_explorer::DocumentTabs documents=tab_items active_document=explorer_active locked=explorer_locked on_activate=explorer_activate on_close=close_document_tab/>
                                    <Show when=move||!has_open_document.get()><section class="empty-document-view" id="no-open-documents" tabindex="-1" aria-label="No open documents"><h2>"No documents open"</h2><p>"Choose a sequence from Documents, or import a file. Closing a tab keeps its document and history in the workspace."</p><div class="empty-document-actions"><button class="button primary" disabled=move||workspace_restoring.get() on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open sequence…"</button><button class="button secondary" disabled=move||workspace_restoring.get() on:click=move|_|{modal_error.set(String::new());modal.set("url");}>"Import from URL…"</button></div></section></Show>
                                    <div class="active-document" hidden=move||!has_open_document.get() id="active-sequence-document" role="tabpanel" aria-label="Active sequence document" tabindex="0">
                                    <div class="document-header"><div><div class="eyebrow">"SEQUENCE"</div><h1>{move||molecule.get().name}</h1><div class="document-meta"><span>{move||format!("{} {}",molecule.get().sequence.len(),molecule.get().sequence.alphabet().short_unit())}</span><span>"·"</span><span>{move||format!("{} {}",if molecule.get().topology==Topology::Circular{"Circular"}else{"Linear"},match molecule.get().sequence.alphabet(){Alphabet::Dna=>"DNA",Alphabet::Rna=>"RNA",Alphabet::Protein=>"protein"})}</span><span>"·"</span><span>{move||molecule.get().sequence.gc_fraction().map(|v|format!("{:.1}% GC",v*100.)).unwrap_or_else(||"GC unavailable".into())}</span></div><Show when=move||molecule.get().metadata.contains_key("source:provider")><p class="source-provenance">{move||{let m=molecule.get();let mut parts=vec![format!("Source: {}",m.metadata.get("source:provider").map(String::as_str).unwrap_or("unknown")),m.metadata.get("source:accession_version").or_else(||m.metadata.get("source:identifier")).cloned().unwrap_or_else(||"record".into())];if let Some(id)=m.metadata.get("source:sequence_id"){parts.push(format!("sequence {id}"));}if let Some(category)=m.metadata.get("source:addgene_completeness"){parts.push(format!("{category} record (provider label)"));}parts.push(format!("retrieved {}",m.metadata.get("source:retrieved_at").map(String::as_str).unwrap_or("unknown")));parts.join(" · ")}}</p></Show></div>
                                        <div class="document-actions"><details class="document-file-menu" data-transient-menu=""><summary class="button secondary">"File ▾"</summary><div class="document-file-menu-content"><button class="button secondary" on:click=save_file>"Save .cg2"</button><button class="button secondary" on:click=move|_|{let mol=molecule.get_untracked();if mol.sequence.alphabet()==Alphabet::Protein {modal_error.set(String::new());modal.set("protein_genpept");return;}match formats::export_with_report(&mol,formats::TextFormat::GenBank){Ok(report)=>crate::io::download::download_with_completion(&format!("{}.gb",mol.name),"text/plain",report.text.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e)} else {message.set(report.warnings.join(" "));}}),Err(e)=>message.set(e.to_string())}}>{move||if molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein){"Export GenPept"}else{"Export GenBank"}}</button><hr/><div class="document-extra-file-actions"><Show when=move||molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein)><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("protein_uniprot");}>"Export UniProt"</button></Show><button class="button secondary" disabled=move||workspace_restoring.get()||!has_open_document.get() on:click=move|_|{modal_error.set(String::new());modal.set("map_export");}>"Export map"</button><button class="button secondary" disabled=move||workspace_restoring.get()||!has_open_document.get() on:click=move|_|{modal_error.set(String::new());modal.set("image_evidence");}>"Images…"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("compat_snapgene");}>"Export SnapGene"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("compat_xdna");}>"Export XDNA"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("compat_gcg");}>"Export GCG"</button><button class="button secondary" on:click=move|_|{let retained=molecule.with_untracked(|m|m.metadata.get("sbol:source_hash").and_then(|hash|attachments.with_untracked(|assets|assets.get(hash).filter(|a|a.kind==crate::io::sbol_exchange::ATTACHMENT_KIND&&a.media_type==crate::io::sbol_exchange::MEDIA_TYPE&&blake3::hash(&a.bytes).to_hex().as_str()==hash).map(|a|std::sync::Arc::new(a.bytes.clone())))));pending_sbol.set(retained);modal_error.set(String::new());modal.set("sbol");}>"Export SBOL"</button><button class="button secondary" on:click=move|_|{modal_error.set(String::new());modal.set("document");}>"Document settings"</button><button class="button secondary" on:click=move|_|{let mol=molecule.get_untracked();match formats::export_with_report(&mol,formats::TextFormat::Fasta){Ok(report)=>crate::io::download::download_with_completion(&format!("{}.fasta",mol.name),"text/plain",report.text.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(error)=result {message.set(error);} else {message.set(report.warnings.join(" "));}}),Err(e)=>message.set(e.to_string())}}>"FASTA"</button><Show when=move||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Protein)><button class="button secondary" on:click=move|_|match formats::export_with_report(&molecule.get_untracked(),formats::TextFormat::Embl){Ok(report)=>{crate::io::download::download_with_completion("sequence.embl","text/plain",report.text.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);} else {message.set(report.warnings.join(" "));}})},Err(e)=>message.set(e.to_string())}>"EMBL"</button></Show><Show when=move||molecule.with(|m|m.metadata.contains_key("protein:original:chunks"))><button class="button secondary" on:click=move|_|{let source=molecule.get_untracked();let guard=session.export_context_guard();match collagene_core::formats::protein::original_bytes(&source){Ok(bytes)=>crate::io::download::download_with_completion_guarded("original-protein-input.txt","text/plain",&bytes,move||guard.run(()),move|result|{if message.try_get_untracked().is_some(){match result{Ok(())=>message.set("Original protein source saved; it does not represent subsequent edits.".into()),Err(e)=>message.set(e)}}}),Err(e)=>message.set(e.to_string())}} >"Original protein input"</button></Show><Show when=move||molecule.with(|m|m.original_genbank.is_some())><button class="button secondary" title="Download the exact original GenBank input, before edits" on:click=move|_|molecule.with_untracked(|m|{if let Some(raw)=&m.original_genbank{crate::io::download::download_with_completion("original-input.gb","text/plain",raw.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})}})>"Original GenBank input"</button></Show><Show when=move||molecule.get().metadata.contains_key(snapgene::ORIGINAL_BASE64_KEY)><button class="button secondary" on:click=move|_|match snapgene::original_snapgene_bytes(&molecule.get_untracked()){Ok(bytes)=>{crate::io::download::download_with_completion("original.dna","application/octet-stream",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})},Err(e)=>message.set(e.to_string())}>"Original .dna"</button></Show><Show when=move||molecule.get().metadata.contains_key("trace_original_base64")><button class="button secondary" on:click=move|_|match traces::original_trace_bytes(&molecule.get_untracked()){Ok(bytes)=>{let extension=if bytes.starts_with(b"ABIF"){"ab1"}else{"scf"};crate::io::download::download_with_completion(&format!("original.{extension}"),"application/octet-stream",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})},Err(e)=>message.set(e.to_string())}>"Original trace"</button></Show><Show when=move||molecule.with(|doc|doc.metadata.contains_key(collagene_core::formats::original::ORIGINAL_BASE64_KEY))><button class="button secondary" on:click=move|_|match collagene_core::formats::original::original_legacy_bytes(&molecule.get_untracked()){Ok(bytes)=>{let extension=if collagene_core::formats::original::is_collagene_xml(&bytes){"xml"}else{"xdna"};crate::io::download::download_with_completion(&format!("original-legacy.{extension}"),"application/octet-stream",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(error)=result {message.set(error);}})},Err(error)=>message.set(error.to_string())}>"Original legacy file"</button></Show></div></div></details></div>
                                    </div>
                                    <crate::restriction::digest::DigestProvenance molecule=molecule/>
                                    <Show when=move||molecule.with(|m|m.metadata.contains_key("crispr:producer_json"))><div class="crispr-historical"><strong>"CRISPR producer provenance · unverified on opening"</strong><p>"This saved sequence design does not restore a current source/guide/target check or prove cellular editing. Recompute against the exact original source and independent target before reuse."</p><button class="button secondary" on:click=move|_|molecule.with_untracked(|m|{if let Some(json)=m.metadata.get("crispr:producer_json"){crate::io::download::download_with_completion("crispr-retained-provenance.json","application/json",json.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})}})>"Download retained CRISPR provenance"</button></div></Show>
                                    <div class="editor-toolbar"><div class="tool-group"><button class="icon-button" title="Show documents" disabled=move||!layout.ready.get() on:click=move|_|show_library.update(|v|*v=!*v)>"☷"</button><button class="icon-button" title="Undo" disabled=move||documents.with(|docs|!docs[active.get()].can_undo()) on:click=move|_|documents.update(|docs|{let d=&mut docs[active.get_untracked()];if let Err(error)=d.undo(){message.set(error.to_string());}})>"↶"</button><button class="icon-button" title="Redo" disabled=move||documents.with(|docs|!docs[active.get()].can_redo()) on:click=move|_|documents.update(|docs|{let d=&mut docs[active.get_untracked()];if let Err(error)=d.redo(){message.set(error.to_string());}})>"↷"</button><span class="tool-separator"></span><button class="text-tool" disabled=move||selection.get().is_none()||inspected_feature.get().is_some_and(|f|f.location.spans().len()>1) on:click=move|_|{edit_text.set(selection.get_untracked().and_then(|(a,b)|molecule.get_untracked().sequence.as_str().get(a..b).map(str::to_owned)).unwrap_or_default());modal_error.set(String::new());modal.set("edit");}>"Edit selection"</button></div>
                                        <div class="viewer-command-menus"><details class="toolbar-action-menu" data-transient-menu=""><summary>"Create ▾"</summary><div class="toolbar-action-menu-content"><button class="text-tool" on:click=move |_|{feature_for_editor.set(None);modal.set("feature");}>"Annotation…"</button><button class="text-tool" on:click=move |_|panel.set("pcr")>"PCR product…"</button><button class="text-tool" on:click=move |_|assembly_visible.set(true)>"Assembly…"</button></div></details><details class="toolbar-action-menu" data-transient-menu=""><summary>"Analyze ▾"</summary><div class="toolbar-action-menu-content"><button class="button secondary" disabled=move||workspace_restoring.get()||molecule.with(|doc|doc.sequence.alphabet()!=Alphabet::Dna) on:click=move|_|{read_reference.set(Some(primer_context.get_untracked()));panel.set("reads");}>"Check sequencing reads"</button><button class="text-tool" on:click=move |_|panel.set("biology")>"Translation and sequence analysis…"</button><button class="text-tool" on:click=move |_|panel.set("compare")>"Compare sequences…"</button><button class="text-tool" on:click=move |_|panel.set("enzymes")>"Virtual digest…"</button></div></details></div>
                                        <crate::workspace::clipboard::ClipboardTools source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) internal_clipboard source_selection=reference_copy_selection on_preview=preview_reference_paste molecule=molecule selection=selection selected_feature=selected_feature highlighted_spans locked=workspace_restoring epoch=workspace_epoch identity=Memo::new(move|_|documents.with(|docs|docs.get(active.get()).map(|doc|doc.id.clone()).unwrap_or_default())) on_commit=commit_annotations message=message/>
                                        <crate::viewer::sequence::search_view::SearchTools molecule=molecule identity=Memo::new(move|_|documents.with(|docs|docs[active.get()].id.clone())) epoch=workspace_epoch selection=selection selected_feature=selected_feature regions=search_regions message=message/>
                                    </div>
                                    <Show when=move||!native_export_notice.get().is_empty()><div class="notice" role="status"><span>{move||native_export_notice.get()}</span><button class="icon-button" aria-label="Dismiss saved export receipt" on:click=move|_|native_export_notice.set(String::new())>"×"</button></div></Show>
                                    <Show when=move||!message.get().is_empty()><div class="notice" role="status"><span>{move||message.get()}</span><button class="icon-button" aria-label="Dismiss message" on:click=move|_|message.set(String::new())>"×"</button></div></Show>
                                    <Show when=move||molecule.with(|m|m.metadata.get("collagene:example").is_some_and(|value|value=="interface-demonstration"))><div class="demo-note mvp-example"><div><strong>"Explore the example, or open your sequence"</strong><span>"Synthetic annotations illustrate the viewer. Open a GenBank file to inspect your own sequence and simulate a digest."</span></div><button class="button secondary" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open a GenBank file"</button></div></Show>
                                    <div class="editor-view-bar"><div class="view-switcher" role="group" aria-label="Editor view">
                                        <button aria-label="Split view" disabled=move||!layout.ready.get() aria-pressed=move||(editor_view.get()=="split").to_string() class:active=move||editor_view.get()=="split" on:click=move|_|editor_view.set("split")>"Split"</button>
                                        <button aria-label="Sequence view" disabled=move||!layout.ready.get() aria-pressed=move||(editor_view.get()=="sequence").to_string() class:active=move||editor_view.get()=="sequence" on:click=move|_|editor_view.set("sequence")>"Sequence"</button>
                                        <button aria-label="Map view" disabled=move||!layout.ready.get() aria-pressed=move||(editor_view.get()=="map").to_string() class:active=move||editor_view.get()=="map" on:click=move|_|editor_view.set("map")>"Map"</button>
                                    </div><div class="restriction-quick-controls"><label><input type="checkbox" aria-label="Show restriction enzymes in editor" prop:checked=move||restrictions.visible.get() disabled=move||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Dna) on:change=move|e|restrictions.visible.set(event_target_checked(&e))/><span>"RE sites"</span></label><button class="text-tool" on:click=move|_|panel.set("restriction")>"Choose enzymes"</button><span class="restriction-quick-count">{move||restrictions.displayed.with(|v|v.as_ref().map(|sites|format!("{} shown",sites.len())).unwrap_or_else(||if restrictions.busy.get(){"Scanning…".into()}else{String::new()}))}</span></div></div>
                                    <crate::viewer::sequence::selection_toolbar::SelectionToolbar molecule=molecule selection=selection selected_feature=selected_feature highlighted_spans=highlighted_spans locked=Memo::new(move|_|workspace_restoring.get())
                                        on_annotate=Callback::new(move|()|{
                                            let location=molecule.with_untracked(|m|selected_feature.get_untracked().and_then(|i|m.features.get(i).map(|f|f.location.clone()))).or_else(||{let parts=highlighted_spans.get_untracked();if parts.is_empty()||parts.len()>4096{return None;}Some(if parts.len()==1{Location::range(parts[0].0,parts[0].1)}else{Location::Join(parts.into_iter().map(|(a,b)|Location::range(a,b)).collect())})});
                                            if let Some(location)=location{feature_for_editor.set(None);feature_initial_location.set(Some(location));modal_error.set(String::new());modal.set("feature");}
                                        })
                                        on_mutation=Callback::new(move|()|match mutation_selection.run(()){Ok(seed)=>{mutation_seed.set(Some(seed));panel.set("mutagenesis");},Err(error)=>message.set(error)})
                                        on_domestication=Callback::new(move|()|panel.set("domestication"))
                                        on_primer=Callback::new(move|reverse:bool|{match pair_selection.run(()).and_then(|selection|selection.oligo(if reverse{"Reverse"}else{"Forward"},reverse,"")){Ok(oligo)=>{match oligo.primer(){Ok(primer)=>{pair_draft.update(|draft|if reverse{draft.reverse=oligo}else{draft.forward=oligo});primer_seed.set(Some(primer));panel.set("primers");},Err(error)=>message.set(error)}},Err(error)=>message.set(error)}})
                                        on_edit=Callback::new(move|()|{if let Some((a,b))=selection.get_untracked(){edit_text.set(molecule.with_untracked(|m|m.sequence.as_str().get(a..b).unwrap_or_default().to_owned()));modal_error.set(String::new());modal.set("edit");}})
                                        on_clear=Callback::new(move|()|batch(||{selection.set(None);selected_feature.set(None);search_regions.set(None);restrictions.selected.set(None);}))
                                    />
                                    <div class="editor-grid" class:sequence-only=move||editor_view.get()=="sequence" class:map-only=move||editor_view.get()=="map">
                                        <section class="sequence-panel"><div class="panel-heading"><div><h2>"Sequence"</h2></div><span class="panel-subtitle">{move||molecule.with(|doc|if doc.sequence.alphabet()==Alphabet::Protein{"Amino-acid reference"}else{"Reference & complement"})}</span></div><Show when=move||molecule.with(|doc|doc.metadata.contains_key("legacy:upper_aligned"))><details class="legacy-source-rows"><summary>"Original aligned source rows"</summary><p class="helper">"Upper is stored 5′→3′; lower is aligned 3′→5′. These original rows may describe an earlier document revision. The editor displays the reference and its calculated complement."</p><pre>{move||molecule.with(|doc|format!("{}\n{}",doc.metadata.get("legacy:upper_aligned").map(|value|value.chars().take(2000).collect::<String>()).unwrap_or_default(),doc.metadata.get("legacy:lower_aligned").map(|value|value.chars().take(2000).collect::<String>()).unwrap_or_default()))}</pre><p class="helper">"First 2,000 columns shown; original legacy-file download retains all source evidence."</p></details></Show><Show when=move||molecule.get().metadata.contains_key("trace_original_base64")><button class="text-tool" on:click=move|_|trace_view.update(|v|*v=!*v)>{move||if trace_view.get(){"Show sequence"}else{"Show chromatogram"}}</button></Show><Show when=move||trace_view.get()&&molecule.get().metadata.contains_key("trace_original_base64") fallback=move||view!{<crate::viewer::sequence::sequence_view::SequenceView antibody_overlay=Some(antibody_overlay) pcr_overlay=Some(pcr_overlay) hidden_features=Some(annotation_hidden) molecule=molecule selection=selection selected_feature=selected_feature highlighted_spans=highlighted_spans search_regions=search_regions on_show_features=Callback::new(move|()|panel.set("features")) identity=Memo::new(move|_|documents.with(|docs|crate::biology::antibody::antibody_overlay::sequence_identity(workspace_epoch.get(),&current_event_binding(&docs[active.get()])))) pcr_context=Some(restriction_identity) locked=Memo::new(move|_|workspace_restoring.get()) on_message=Callback::new(move|text:String|message.set(text)) on_history=sequence_on_history restriction_sites=restrictions.displayed restriction_visible=restrictions.visible on_show_restrictions=Callback::new(move|()|panel.set("restriction")) restriction_focus=Some(hovered_restriction) on_select_restriction=navigate_restriction on_edit=sequence_on_edit/>}><crate::reads::trace_view::TraceView molecule=molecule selection=selection source_scope=Memo::new(move |_| {
        let event = documents.with(|docs| docs.get(active.get()).map(current_event_binding));
        format!("{}:{}", workspace_epoch.get(), serde_json::to_string(&event).unwrap_or_default())
    })/></Show></section>
                                        <crate::workspace::workspace_layout::PaneResize layout=layout pane=crate::workspace::workspace_layout::Pane::Sequence/>
                                        <section class="map-panel"><div class="panel-heading"><div><h2>{move||if molecule.get().topology==Topology::Circular{"Plasmid map"}else{"Linear map"}}</h2></div><div class="map-controls"><button class="icon-button" title="Zoom out" on:click=move|_|if let Some(r)=renderer.get_value(){r.zoom(1./1.15);}>"−"</button><button class="icon-button" title="Reset map view" on:click=move|_|if let Some(r)=renderer.get_value(){r.reset();}>"Fit"</button><button class="icon-button" title="Zoom in" on:click=move|_|if let Some(r)=renderer.get_value(){r.zoom(1.15);}>"+"</button></div></div>
                                            <div class="map-stage"><canvas node_ref=gl_canvas aria-label=move||if molecule.with(|m|m.topology==Topology::Circular){"Interactive circular sequence overview; drag to rotate, scroll to zoom, click annotations to inspect"}else{"Interactive linear sequence overview; scroll to zoom, click annotations to inspect"} on:pointerdown=move|e|{if let Some(canvas)=gl_canvas.get(){let _=canvas.set_pointer_capture(e.pointer_id());}if let Some(r)=renderer.get_value(){r.pointer_down(e.client_x()as f64,e.client_y()as f64);}} on:pointermove=move|e|if let Some(r)=renderer.get_value(){r.pointer_move(e.client_x()as f64,e.client_y()as f64);} on:pointerup=move|e|if let Some(r)=renderer.get_value(){if let Some(hit)=r.pointer_up(e.client_x()as f64,e.client_y()as f64){if let Some(index)=hit.restriction_index{if let Some(site)=restrictions.scan.with_untracked(|s|s.as_ref().and_then(|s|s.sites.iter().find(|site|site.index==index).cloned())){navigate_restriction.run(site);}}else{batch(||{restrictions.selected.set(None);search_regions.set(None);selected_feature.set(hit.feature_index);selection.set(if hit.spans.len()==1{Some(hit.span)}else{None});if hit.spans.len()>1{search_regions.set(Some((molecule_revision.get_untracked(),hit.spans.clone())));}if hit.feature_index.is_some(){panel.set("features");}});}}} on:pointercancel=move|_|if let Some(r)=renderer.get_value(){r.pointer_cancel();} on:lostpointercapture=move|_|if let Some(r)=renderer.get_value(){r.pointer_cancel();} on:wheel=move|e|{e.prevent_default();if let Some(r)=renderer.get_value(){r.zoom(if e.delta_y()>0.{0.95}else{1.05});}}></canvas><canvas node_ref=label_canvas class="map-labels" aria-hidden="true"></canvas><p class="map-geometry-note" role="status">"Overview is simplified. Browse the annotation pages or export the full annotation set."</p><p class="restriction-map-note" role="status">"Some restriction labels are hidden for space. Open RE sites to inspect the selected positions."</p>
                                                <Show when=move||!render_error.get().is_empty()><div class="map-error">{move||render_error.get()}</div></Show>
                                            </div><div class="map-caption"><span>{move||if molecule.get().topology==Topology::Circular{"Drag to rotate · scroll to zoom · select exact ranges in Sequence"}else{"Scroll to zoom · click features or sites to inspect"}}</span></div>
                                            <Show when=move||{let stats=map_labels.get();stats.omitted>0||stats.abbreviated>0||stats.limited}><div class="map-label-note"><span>{move||{let stats=map_labels.get();let mut parts=vec![format!("{} of {} annotation labels shown",stats.displayed,stats.total)];if stats.abbreviated>0{parts.push(format!("{} abbreviated",stats.abbreviated));}if stats.unplaced>0{parts.push(format!("{} without local coordinates",stats.unplaced));}if stats.lane_overflow>0{parts.push(format!("{} overlapping features hidden (3 map lanes)",stats.lane_overflow));}if stats.limited{parts.push("display limit reached".into());}parts.join(" · ")}}</span><button class="text-tool" on:click=move|_|panel.set("features")>"Inspect annotations"</button></div></Show>
                                            <div class="legend"><Show when=move||restrictions.visible.get()><span><i style="background:#8b5ca4"></i>"RE recognition"</span><span><i style="background:#2563a6"></i>"Top cut"</span><span><i style="background:#b45c19"></i>"Bottom cut"</span></Show><span><i style="background:#2b8775"></i>"Coding sequence"</span><span><i style="background:#d6933b"></i>"Promoter"</span><span><i style="background:#666bb3"></i>"Origin"</span></div>
                                        </section>
                                        <aside class="inspector viewer-inspector" class:is-open=move||!panel.get().is_empty()><crate::platform::desktop_shell::InspectorNavigation panel=panel on_assemble=Callback::new(move|()|assembly_visible.set(true)) on_auxiliary=Callback::new(move|action|{use crate::platform::desktop_shell::AuxiliaryAction;modal_error.set(String::new());match action{AuxiliaryAction::GoldenGateWorkflow=>{method_config.update(|config|config.mode=crate::cloning::preparation::cloning_methods_model::Mode::GoldenGate);panel.set("methods");},AuxiliaryAction::ImportUrl=>modal.set("url"),AuxiliaryAction::Alignments=>{pending_alignment.set(None);modal.set("alignments");},AuxiliaryAction::Oligos=>{existing_oligo.set(None);pending_oligo.set(None);modal.set("oligos");}}})/><div class="inspector-body"><div class="inspector-drawer-heading"><h2>{move||crate::platform::desktop_shell::panel_title(panel.get())}</h2><button class="icon-button" aria-label="Close tool panel" title="Close tool panel" on:click=move |_|panel.set("")>"×"</button></div>
                                            <Show when=move||panel.get()=="restriction"><crate::restriction::restriction_view::RestrictionPanel state=restrictions on_navigate=navigate_restriction saved=saved_enzyme_sets on_save=save_enzyme_set on_delete=delete_enzyme_set/></Show>
                                            <Show when=move||panel.get()=="library"><div class="inspector-content"><crate::viewer::sequence::feature_tools::FeatureLibraryPanel molecule=molecule selected_feature=selected_feature libraries=libraries on_keep=keep_library on_commit=commit_annotations selection=selection/></div></Show>
                                            <Show when=move||panel.get()=="annotation_files"><div class="inspector-content"><crate::viewer::sequence::feature_tools::AnnotationFilesPanel molecule=molecule source=annotation_source format=annotation_format on_commit=commit_annotations selection=selection selected_feature=selected_feature/></div></Show>
                                            <Show when=move||panel.get()=="history"><div class="inspector-content"><crate::genomes::variant_source_library::VariantSourceLibrary entries=Memo::new(move|_|attachments.with(|assets|assets.iter().filter_map(|(hash,asset)|collagene_core::review::variants::archive::inspect(asset).ok().map(|info|(hash.clone(),info))).collect())) locked=Memo::new(move|_|workspace_restoring.get()) on_open=Callback::new(move|hash:String|{if workspace_restoring.get_untracked(){return;}let value=attachments.with_untracked(|assets|assets.get(&hash).and_then(|asset|collagene_core::review::variants::archive::inspect(asset).ok().map(|info|(std::sync::Arc::new(asset.bytes.clone()),info.format,info.name))));if let Some(value)=value{pending_variants.set(Some(value));modal_error.set(String::new());modal.set("variants");}}) on_remove=Callback::new(move|hash:String|{if !workspace_restoring.get_untracked(){attachments.update(|assets|{if assets.get(&hash).is_some_and(|asset|collagene_core::review::variants::archive::inspect(asset).is_ok()){assets.remove(&hash);}});}})/><crate::reads::read_source_library::ReadSourceLibrary entries=Memo::new(move|_|attachments.with(|assets|assets.iter().filter_map(|(hash,asset)|collagene_core::review::reads::archive::inspect(asset).ok().map(|info|(hash.clone(),info))).collect())) locked=Memo::new(move|_|workspace_restoring.get()) on_open=Callback::new(move|hash:String|{if workspace_restoring.get_untracked(){return;}let source=attachments.with_untracked(|assets|assets.get(&hash).and_then(|asset|collagene_core::review::reads::archive::inspect(asset).ok().map(|info|(std::sync::Arc::new(asset.bytes.clone()),match info.format{collagene_core::review::reads::store::Format::Sam=>collagene_core::review::reads::session::FileFormat::Sam,collagene_core::review::reads::store::Format::Bam=>collagene_core::review::reads::session::FileFormat::Bam},info.name))));if let Some(source)=source{pending_read_alignment.set(Some(source));modal_error.set(String::new());modal.set("read_alignments");}}) on_remove=Callback::new(move|hash:String|{if !workspace_restoring.get_untracked(){attachments.update(|assets|{if assets.get(&hash).is_some_and(|a|collagene_core::review::reads::archive::inspect(a).is_ok()){assets.remove(&hash);}});}})/><crate::history::history_view::HistoryPanel molecule=molecule events=Memo::new_owning(move|_|(documents.with(|docs|docs[active.get()].summaries()),true)) current_event=Memo::new(move|_|documents.with(|docs|docs[active.get()].head_event().to_string())) navigation=Memo::new(move|_|documents.with(|docs|docs[active.get()].navigation_event_ids())) external_events=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().filter(|(i,_)|*i!=active.get()).flat_map(|(_,doc)|{let label=format!("[{}] {}",doc.id,doc.molecule.name.chars().take(100).collect::<String>());doc.summaries().into_iter().map(move|e|(e.event_id,label.clone(),e.revision))}).collect())) on_restore_revision=restore_revision on_fork_revision=fork_revision candidates=Memo::new(move|_|{let current=candidate.get().map(|c|c.plan_revision);stored_candidates.with(|values|values.iter().map(|c|(c.plan_revision.clone(),c.plan.name.clone(),c.plan.target_revision.clone(),current.as_ref()==Some(&c.plan_revision))).collect())}) on_open_candidate=Callback::new(move|revision:String|{if workspace_restoring.get_untracked(){return;}if let Some(value)=stored_candidates.with_untracked(|values|values.iter().find(|c|c.plan_revision==revision).cloned()){candidate_events.set(None);candidate.set(Some(value));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");}}) on_remove_candidate=Callback::new(move|revision:String|{if workspace_restoring.get_untracked(){return;}batch(||{if candidate.with_untracked(|c|c.as_ref().is_some_and(|c|c.plan_revision==revision)){candidate_events.set(None);candidate.set(None);plan_report.set(None);last_plan.set(None);last_plan_target.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));}stored_candidates.update(|values|values.retain(|c|c.plan_revision!=revision));});}) collections=Memo::new(move|_|oligo_tables.with(|tables|tables.iter().enumerate().map(|(i,t)|(i,format!("Oligo collection {}",t.source_blake3.chars().take(12).collect::<String>()),t.records.len())).collect())) on_open_oligos=Callback::new(move|index:usize|{if let Some(table)=oligo_tables.with_untracked(|tables|tables.get(index).cloned()){existing_oligo.set(Some(table.as_ref().clone()));pending_oligo.set(None);modal_error.set(String::new());modal.set("oligos");}})/><Show when=move||!recorded_runs.with(|runs|runs.is_empty())><section class="recorded-runs"><h3>"Recorded reconstructions"</h3><p class="helper">"Saved computational history. Each record needs a fresh replay before its verification can be used."</p>{move||recorded_runs.with(|runs|runs.iter().skip(run_page.get()*20).take(20).map(|run|{let id=run.run_id.clone();view!{<article class="run-product"><strong>{run.plan_name.clone()}</strong><p class="helper">{format!("{} products · recorded {:?} · {} · unverified",run.products.len(),run.recorded_status,run.recorded_engine)}</p><button class="button secondary" on:click=move|_|{selected_run.set(id.clone());modal_error.set(String::new());modal.set("run_details");}>"View cloning run"</button></article>}}).collect_view())}<div class="pagination"><button class="button secondary" disabled=move||run_page.get()==0 on:click=move|_|run_page.update(|p|*p=p.saturating_sub(1))>"Previous records"</button><button class="button secondary" disabled=move||{(run_page.get()+1)*20>=recorded_runs.with(|r|r.len())} on:click=move|_|run_page.update(|p|*p+=1)>"Next records"</button></div></section></Show><Show when=move||!saved_assemblies.with(|rows|rows.is_empty())><section class="saved-assemblies"><h3>"Saved read assemblies"</h3><p class="helper">"Historical evidence only; replay exact retained reads before using a consensus."</p>{move||saved_assemblies.with(|rows|rows.iter().skip(assembly_page.get()*20).take(20).map(|row|{let hash=row.attachment_hash.clone();view!{<article><p>{format!("{} reads · {} contigs · recorded complete {} · {}",row.read_count,row.contig_count,row.recorded_complete,row.engine_version)}</p><button class="button secondary" on:click=move|_|{inspected_assembly.set(hash.clone());modal_error.set(String::new());modal.set("saved_contig");}>"View saved read assembly"</button></article>}}).collect_view())}<button class="button secondary" disabled=move||assembly_page.get()==0 on:click=move|_|assembly_page.update(|p|*p=p.saturating_sub(1))>"Previous saved assemblies"</button><button class="button secondary" disabled=move||{(assembly_page.get()+1)*20>=saved_assemblies.with(Vec::len)} on:click=move|_|assembly_page.update(|p|*p+=1)>"Next saved assemblies"</button></section></Show></div></Show>
                                            <Show when=move||panel.get()=="reads"><div class="inspector-content"><h2>"Read validation"</h2><Show when=move||molecule.with(|m|m.metadata.contains_key(oligo_formats::FASTQ_ORIGINAL_KEY))><section><h3>"FASTQ quality evidence"</h3><p class="helper">{move||match oligo_formats::quality_evidence(&molecule.get()){Ok(e)=>format!("{} calls · {:?} · {}. {}",e.qualities.len(),e.encoding,if e.current_valid{"Quality mapping matches current calls"}else{"Current calls changed: quality mapping is stale"},e.warnings.join(" ")),Err(e)=>e.to_string()}}</p><label>"FASTQ output encoding"<select aria-label="FASTQ output encoding" prop:value=move||match fastq_encoding.get(){None=>"",Some(oligo_formats::QualityEncoding::Phred33)=>"33",Some(oligo_formats::QualityEncoding::Phred64)=>"64"} on:change=move|e|fastq_encoding.set(match event_target_value(&e).as_str(){"33"=>Some(oligo_formats::QualityEncoding::Phred33),"64"=>Some(oligo_formats::QualityEncoding::Phred64),_=>None})><option value="">"Choose encoding"</option><option value="33">"Phred+33"</option><option value="64">"Phred+64"</option></select></label><button class="button secondary full" disabled=move||fastq_encoding.get().is_none() on:click=move|_|{let Some(encoding)=fastq_encoding.get_untracked()else{return;};match oligo_formats::export_fastq(&molecule.get_untracked(),encoding){Ok(report)=>crate::io::download::download_with_completion("current.fastq","text/plain",report.text.as_bytes(),move|result|{if message.try_get_untracked().is_none(){return;}match result{Ok(())=>message.set(report.warnings.join(" ")),Err(e)=>message.set(e)}}),Err(e)=>message.set(e.to_string())}}>"Export current FASTQ"</button><button class="button secondary full" on:click=move|_|match oligo_formats::original_fastq_bytes(&molecule.get_untracked()){Ok(bytes)=>crate::io::download::download_with_completion("original.fastq","text/plain",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e)}}),Err(e)=>message.set(e.to_string())}>"Original FASTQ record"</button><button class="button secondary full" disabled=move||selection.get().is_none() on:click=move|_|{let Some((a,b))=selection.get_untracked()else{return;};match oligo_formats::trim_fastq(&molecule.get_untracked(),a,b){Ok(doc)=>if let Err(e)=session.add_documents(vec![doc]){message.set(e)},Err(e)=>message.set(e.to_string())}}>"Open selected FASTQ calls"</button></section></Show><label>"Reference construct"<select aria-label="Read check reference construct" prop:value=move||read_check_source.get().binding.document_id on:change=move|e|{let id=event_target_value(&e);if let Some(binding)=documents.with_untracked(|docs|docs.iter().find(|doc|doc.id==id).map(current_event_binding)){read_reference.set(Some(binding));}}>{move||documents.with(|docs|docs.iter().filter(|doc|doc.molecule.sequence.alphabet()==Alphabet::Dna).map(|doc|view!{<option value=doc.id.clone()>{doc.molecule.name.clone()}</option>}).collect_view())}</select></label>
                                            <crate::reads::read_checks_view::SavedReadChecksPanel on_review=open_read_review event_available=primer_anchor_available attachment_available=read_check_attachment_available reference=read_check_reference snapshot=contig_snapshot epoch=workspace_epoch on_read=navigate_read_call on_reference=navigate_check_reference/>
                                            <crate::reads::read_review_history::ReadReviewHistory reference=read_check_reference snapshot=contig_snapshot epoch=workspace_epoch on_open=reopen_read_review on_reference=navigate_check_reference/>
                                            <crate::reads::trace_validation::TraceValidationPanel molecule=read_check_molecule reference_info=read_check_source source_info=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(recombination_source_info))) epoch=workspace_epoch on_read=navigate_read_call on_record=record_read_check on_reference=navigate_check_reference query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) selection=selection selected_feature=selected_feature on_consensus=Callback::new(move|doc:Molecule|session.add_documents(vec![doc]))/></div></Show>
                                            <Show when=move||panel.get()=="features"><div class="inspector-content"><h2>"Annotations"<span class="count">{move||molecule.get().features.len()}</span></h2><p class="helper">"Select an annotation to inspect all of its parts and biological sequence."</p><crate::viewer::sequence::annotation_table::AnnotationTable molecule context=restriction_identity locked=Memo::new(move|_|workspace_restoring.get()) selected_feature hidden=annotation_hidden on_visibility=set_annotation_visibility on_apply=apply_annotation_batch on_inspect=Callback::new(move|index:usize|{molecule.with_untracked(|m|{if let Some(feature)=m.features.get(index){batch(||{selected_feature.set(Some(index));selection.set(feature.location.spans().first().map(|(a,b,_)|(*a,*b)));search_regions.set(None);});}});}) on_edit=Callback::new(move|index:Option<usize>|{feature_for_editor.set(index);feature_initial_location.set(None);modal.set("feature");})/>{move||inspected_feature.get().map(|f|{let location=collagene_model::annotations::insdc::location_text(&f.location).unwrap_or_else(|_|"unresolved".into());let size=f.location.spans().iter().try_fold(0usize,|n,(a,b,_)|n.checked_add(b-a)).unwrap_or(usize::MAX);let exact=if size<=10_000 {f.location.extract(&molecule.get().sequence).unwrap_or_else(|e|e.to_string())}else{"Sequence exceeds this inspector's 10,000-symbol display limit; use a region/sequence export.".into()};view!{<div class="annotation-details"><h3>{f.label().to_string()}</h3><button class="text-tool" on:click=move|_|panel.set("library")>"Use selected annotation in library"</button><p>{location}</p><pre class="annotation-sequence">{exact}</pre><p class="helper">"Compound selections highlight every span. Choose contiguous bases to edit the underlying sequence."</p><button class="button secondary" on:click=move|_|{feature_for_editor.set(selected_feature.get_untracked());modal.set("feature");}>"Edit annotation"</button><button class="text-tool" on:click=move|_|if let Some(i)=selected_feature.get_untracked(){let current=molecule.get_untracked();match editing::delete_feature(&current,i){Ok(edited)=>if let Err(e)=session.commit_document(&current.revision(),edited){message.set(e);},Err(e)=>message.set(e.to_string())}}>"Delete annotation"</button><details><summary>"Annotation qualifiers"</summary><dl>{f.qualifiers.into_iter().take(200).map(|(key,value)|view!{<dt>{key}</dt><dd>{value.unwrap_or_else(||"flag".into())}</dd>}).collect_view()}</dl></details></div>}})}</div></Show>
                                            <Show when=move||panel.get()=="enzymes"><crate::restriction::digest::DigestPanel origin=digest_run_origin get_workspace=contig_snapshot saved=saved_digests on_save=save_digest_setup on_load=load_digest_setup on_delete=delete_digest_setup on_gel=queue_digest_gel on_fork=fork_digest_source source=Memo::new(move|_|documents.with(|docs|{let d=&docs[active.get()];crate::sources::identity::SourceInfo{binding:collagene_core::material::lineage::EventBinding{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:molecule_revision.get()},name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet(),topology:d.molecule.topology}})) get_molecule=Callback::new(move|()|molecule.get_untracked()) epoch=workspace_epoch on_navigate=navigate_digest on_open=open_digest on_assembly=use_digest_in_assembly boundary_pick=digest_boundary_pick initial_enzymes=restrictions.names.get_untracked() on_methylation=Callback::new(move|()|panel.set("methylation"))/></Show>
                                            <Show when=move||panel.get()=="pcr"><div class="inspector-content"><h2>"Simulate PCR"</h2><p class="helper">"Enter both primers 5′→3′. Annealing must match exactly; tails are introduced into the product."</p><crate::primers::primer_picker::PrimerPicker label="Forward".to_string() sources=primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer epoch=workspace_epoch on_copy=Callback::new(move|primer:Primer|batch(||{pcr_forward_name.set(primer.name);pcr_forward.set(primer.annealing.as_str().into());pcr_forward_tail.set(primer.tail);pcr_generation.update(|g|*g=g.wrapping_add(1));}))/>
            <crate::primers::primer_picker::PrimerPicker label="Reverse".to_string() sources=primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer epoch=workspace_epoch on_copy=Callback::new(move|primer:Primer|batch(||{pcr_reverse_name.set(primer.name);pcr_reverse.set(primer.annealing.as_str().into());pcr_reverse_tail.set(primer.tail);pcr_generation.update(|g|*g=g.wrapping_add(1));}))/>
            <label>"Forward name"<input maxlength="256" prop:value=move||pcr_forward_name.get() on:input=move|e|pcr_forward_name.set(event_target_value(&e))/></label>
            <crate::history::run_origin_view::SourceOriginPicker source=Memo::new(move |_|documents.with(|docs|docs.get(active.get()).map(|doc|crate::sources::identity::SourceInfo{binding:current_event_binding(doc),name:doc.molecule.name.clone(),length:doc.molecule.sequence.len(),alphabet:doc.molecule.sequence.alphabet(),topology:doc.molecule.topology}))) epoch=workspace_epoch get_workspace=contig_snapshot selected=Memo::new(move |_|pcr_origin.get()) on_select=Callback::new(move |origin:Option<crate::history::run_origin_model::PinnedOrigin>|{if origin.as_ref().is_some_and(|origin|documents.with_untracked(|docs|docs.get(active.get_untracked()).is_none_or(|doc|current_event_binding(doc)!=*origin.event()))){return Err("PCR source changed; select its material origin again.".into());}pcr_origin.set(origin);Ok(())})/>
            <label>"Reverse name"<input maxlength="256" prop:value=move||pcr_reverse_name.get() on:input=move|e|pcr_reverse_name.set(event_target_value(&e))/></label><label>"Forward annealing"<input prop:value=move||pcr_forward.get() on:input=move|e|pcr_forward.set(event_target_value(&e)) placeholder="ACGT…"/></label><label>"Forward 5′ tail"<input prop:value=move||pcr_forward_tail.get() on:input=move|e|pcr_forward_tail.set(event_target_value(&e)) placeholder="Optional"/></label><label>"Reverse annealing"<input prop:value=move||pcr_reverse.get() on:input=move|e|pcr_reverse.set(event_target_value(&e)) placeholder="ACGT…"/></label><label>"Reverse 5′ tail"<input prop:value=move||pcr_reverse_tail.get() on:input=move|e|pcr_reverse_tail.set(event_target_value(&e)) placeholder="Optional"/></label><button class="button primary full" on:click=do_pcr>"Compute products"</button><p class="helper">"Review reported geometries, then choose a named reference preview and optional primers. This checks sequence geometry, not laboratory success."</p>{move||pcr_review.get().map(|review|view!{<crate::primers::pcr::pcr_outputs::PcrOutputsPanel source_molecule=molecule on_overlay=publish_pcr_overlay on_navigate=navigate_pcr_overlay review=review on_create=create_pcr_outputs guard=pcr_review_current busy=pcr_output_busy/>})}</div></Show>
                                            <Show when=move||panel.get()=="plan"><div class="inspector-content"><h2>"Cloning plan"</h2><div class="plan-exchange"><button class="button secondary full" disabled=move||workspace_restoring.get()||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Dna) on:click=move|_|{modal_error.set(String::new());modal.set("plan_new");}>"New manual plan"</button><button class="button secondary full" on:click=export_workspace>"Export for AI"</button><button class="button secondary full" on:click=move|_|if let Some(input)=file_input.get(){input.click();}>"Open candidate"</button></div><p class="helper">"Select a separate intended target, then create a manual plan or open an AI candidate. Review the steps and verify before creating products."</p><Show when=move||candidate.get().is_some()><div class="candidate-review"><h3>"Candidate · review before use"</h3><p class="helper">{move||candidate_events.get().filter(|(revision,epoch,_,_)|candidate.with(|c|c.as_ref().is_some_and(|c|c.plan_revision==*revision))&&*epoch==workspace_epoch.get()).map(|(_,_,target,_)|format!("Independent retained target: document {} · event {} · revision {}. The viewer document is unchanged.",target.document_id,target.event_id,target.revision)).unwrap_or_else(||"Unbound candidate: the active document must match its independent intended target.".into())}</p><button class="button secondary full" disabled=move||workspace_restoring.get() on:click=move|_|{modal_error.set(String::new());modal.set("plan_edit");}>"Edit steps and sources"</button>{move||candidate.get().map(|c|view!{<p><strong>{c.plan.name.clone()}</strong></p><p class="helper">{format!("Target revision {}…",&c.plan.target_revision[..12])}</p><h4>"Bound inputs"</h4><ul>{c.plan.inputs.iter().map(|i|view!{<li>{format!("{}: {} · {} bases · {}…",i.id,i.molecule.name,i.molecule.sequence.len(),&i.revision[..12])}</li>}).collect_view()}</ul><h4>"Ordered steps"</h4><ol>{c.plan.steps.iter().map(|step|{let description=match step{Operation::ResectThreePrime{cuts,..}=>format!("Shorten actual 3′ suffixes: left {} bp, right {} bp",cuts.left,cuts.right),Operation::HomologyAnnealing{fragments,..}=>format!("Anneal {} fragments; retain unfilled gaps and unsealed nicks",fragments.len()),Operation::HomologyAssembly{fragments,..}=>format!("Reconstruct homology-assembly strands from {} fragments",fragments.len()),Operation::OverlapExtension{options,..}=>format!("Reciprocal overlap extension across {} bp of terminal homology",options.overlap),Operation::ReconstructRunProduct{..}=>"Reconstruct saved material from its original run".into(),Operation::ReconstructMaterial{..}=>"Reconstruct saved material from its exact native derivation".into(),Operation::CatalogDigest{selection,..}=>format!("Catalog digest with {}",selection.enzymes.iter().map(|e|e.enzyme.name.as_str()).collect::<Vec<_>>().join(", ")),Operation::DeclareSourceTermini{..}=>"Declare original duplex termini".into(),Operation::TopoActions{..}=>"Accounted flap processing and physical duplex reconstruction".into(),Operation::CheckedTopo{..}=>"TOPO graph with explicit target acceptance".into(),Operation::Topo{..}=>"TOPO physical strand reconstruction".into(),Operation::MultiSiteRecombineReference{request,..}=>format!("Multi-parent {:?} reference exchange",request.reaction.reaction),Operation::RecombineReference{request,..}=>format!("Reference-only {:?} exchange",request.reaction),Operation::DeclareSourceModifications{..}=>"Declare original source modifications".into(),Operation::CheckedDigest{selection,..}=>format!("Checked digest with {}",selection.enzymes.iter().map(|e|e.enzyme.name.as_str()).collect::<Vec<_>>().join(", ")),Operation::Tail3Prime{base,..}=>format!("Add {} to both 3′ ends",base.as_str()),Operation::ProfiledOverlap{options,..}=>format!("Profiled {:?} overlap assembly",options.method),Operation::RepairEnds{method,..}=>format!("End repair: {:?}",method),Operation::AnnealOligos{candidate_index,..}=>format!("Anneal oligos: zero-based candidate {}",candidate_index),Operation::PcrWithMismatches{..}=>"Intentional-mismatch PCR".into(),Operation::SelectedGoldenGate{options,candidate_index,..}=>format!("Golden Gate with {:?} · explicit candidate {} and exact keys",options.enzyme,candidate_index+1),Operation::GoldenGate{options,..}=>format!("Golden Gate with {:?}",options.enzyme),Operation::Digest{enzymes,..}=>format!("Digest with {}",enzymes.join(", ")),Operation::TreatEnds{treatment,..}=>format!("Terminal {treatment:?}"),Operation::Ligate{fragments,circular,..}=>format!("Ligation of {} fragments → {}",fragments.len(),if *circular {"circular"}else{"linear"}),Operation::Pcr{forward,reverse,..}=>format!("PCR with {} and {}",forward.name,reverse.name),Operation::Overlap{fragments,circular,..}=>format!("Overlap assembly of {} fragments → {}",fragments.len(),if *circular {"circular"} else {"linear"})};let details=describe_step(step);view!{<li><details><summary>{description}</summary><pre>{details}</pre></details></li>}}).collect_view()}</ol>})}<button class="button primary full" on:click=verify_candidate>"Verify candidate"</button><button class="button secondary full" on:click=move|_|{candidate_events.set(None);candidate.set(None);plan_report.set(None);last_plan.set(None);last_plan_target.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));}>"Discard candidate"</button></div></Show><Show when=move||candidate.get().is_none()><h3>"Manual cloning plan"</h3><p class="helper">"The active document is the intended target. Choose source documents and a cloning method; the core reconstructs the ordered candidate."</p><label>"First fragment"<select aria-label="First source document" prop:value=move||first_input.get().to_string() on:change=move|e|first_input.set(event_target_value(&e).parse().unwrap_or(0))>{move||documents.get().iter().enumerate().map(|(i,d)|view!{<option value=i.to_string()>{d.molecule.name.clone()}</option>}).collect_view()}</select></label><Show when=move||plan_method.get()!="restriction"||(restriction_final.get()=="ligate"&&second_restriction.get())><label>"Second fragment"<select aria-label="Second source document" prop:value=move||second_input.get().to_string() on:change=move|e|second_input.set(event_target_value(&e).parse().unwrap_or(0))>{move||documents.get().iter().enumerate().map(|(i,d)|view!{<option value=i.to_string()>{d.molecule.name.clone()}</option>}).collect_view()}</select></label></Show><label>"Method"<select aria-label="Cloning method" prop:value=move||plan_method.get() on:change=move|e|plan_method.set(if event_target_value(&e)=="restriction"{"restriction"}else{"overlap"})><option value="overlap">"Overlap assembly"</option><option value="restriction">"Restriction and ligation"</option></select></label><Show when=move||plan_method.get()=="restriction"><div class="checked-restriction-options"><label>"Restriction model"<select aria-label="Restriction model" prop:value=move||restriction_policy.get() on:change=move|e|restriction_policy.set(if event_target_value(&e)=="checked"{"checked"}else{"assumed"})><option value="assumed">"Assume complete cleavage"</option><option value="checked">"Check inherited modification state"</option></select></label><label>"Final restriction product"<select aria-label="Restriction final product" prop:value=move||restriction_final.get() on:change=move|e|restriction_final.set(if event_target_value(&e)=="fragment"{"fragment"}else{"ligate"})><option value="ligate">"Ligate selected substrates"</option><option value="fragment">"Single processed fragment (no ligation)"</option></select></label><p class="helper">"Single-fragment mode keeps the selected digest product, or its explicitly repaired/treated descendant, as the final product. A blocked circle can remain an uncut circle without claiming a ligation."</p><Show when=move||restriction_final.get()=="ligate"><label class="checkbox-label"><input type="checkbox" aria-label="Use second restriction substrate" prop:checked=move||second_restriction.get() on:change=move|e|second_restriction.set(event_target_checked(&e))/>"Use a second restriction substrate"</label></Show><Show when=move||restriction_policy.get()=="checked"><p class="helper">"Declare the original source generation under a closed Dam/Dcm/CpG model. Original eligible positions inherit the declared state; positions outside these original motifs and newly synthesized bases are unmodified under this scope. Descendant motifs are not rescanned. Other modification systems and experimental conditions are not modeled. Unknown or impaired checked cleavage withholds undetermined products; blocked sites stay uncut."</p><details class="checked-profile-identities"><summary>"Pinned positional profiles and formulation identities"</summary><p>{format!("{} · revision {}",modified_duplex::PROFILE_VERSION,modified_duplex::profile_revision())}</p>{methylation::profiles().into_iter().map(|p|view!{<p>{format!("{}: {}",p.enzyme,p.supplier_identity)}</p>}).collect_view()}<p class="helper">"This positional profile is distinct from the source-context methylation review. Generic and HF/v2 enzyme formulations are not substituted. Imported summaries cannot authorize checked state."</p></details></Show>
                <RestrictionFields side="First" enzymes=first_enzymes fragment=first_fragment treatment=first_treatment repair=first_repair reverse=first_reverse/><Show when=move||restriction_policy.get()=="checked"><p class="checked-source-binding">{move||documents.with(|docs|docs.get(first_input.get()).map(|d|format!("First source: document {} · event {}",d.id,d.head_event())).unwrap_or_else(||"Choose a first source".into()))}</p><crate::cloning::methods::checked_cloning::SourceAssumptions side="First" signals=first_modifications/></Show>
                <Show when=move||second_restriction.get()&&restriction_final.get()=="ligate"><RestrictionFields side="Second" enzymes=second_enzymes fragment=second_fragment treatment=second_treatment repair=second_repair reverse=second_reverse/><Show when=move||restriction_policy.get()=="checked"><p class="checked-source-binding">{move||documents.with(|docs|docs.get(second_input.get()).map(|d|format!("Second source: document {} · event {}",d.id,d.head_event())).unwrap_or_else(||"Choose a second source".into()))}</p><crate::cloning::methods::checked_cloning::SourceAssumptions side="Second" signals=second_modifications/></Show></Show><Show when=move||restriction_final.get()=="ligate"><label class="checkbox-label"><input type="checkbox" aria-label="Allow restriction ligation nicks" prop:checked=move||allow_nicked.get() on:change=move|e|allow_nicked.set(event_target::<HtmlInputElement>(&e).checked())/>"Allow explicitly reported nicks"</label></Show><p class="helper">"Fragment numbers follow ascending reference-strand cut order, not size. Intact circles cannot ligate. The active document remains the independent target; declared source state is not measured modification evidence."</p></div></Show><Show when=move||plan_method.get()=="overlap"><label>"Junction overlap (bp)"<input type="number" min="1" prop:value=move||overlap.get() on:input=move|e|overlap.set(event_target_value(&e))/></label><Show when=move||molecule.get().topology==Topology::Circular><label>"Circular closure overlap (bp)"<input type="number" min="1" prop:value=move||closure_overlap.get() on:input=move|e|closure_overlap.set(event_target_value(&e))/></label></Show></Show><button class="button primary full" on:click=do_plan>"Reconstruct and verify"</button></Show>{move||plan_report.with(|report|report.as_ref().map(|report|view!{<div class="verification-result"><strong>{format!("{:?}",report.status)}</strong>{report.findings.iter().cloned().map(|f|view!{<p>{f}</p>}).collect_view()}{report.physical_processing_reports.iter().map(|step|{match &step.evidence{collagene_core::verification::evidence::PhysicalProcessingEvidence::ThreePrimeResection(evidence)=>crate::assembly::review::assembly_product_view::three_prime_processing_view(step.step_index,evidence),collagene_core::verification::evidence::PhysicalProcessingEvidence::HomologyAnnealing(evidence)=>crate::assembly::review::assembly_product_view::annealing_processing_view(step.step_index,evidence),collagene_core::verification::evidence::PhysicalProcessingEvidence::HomologyAssembly(evidence)=>crate::assembly::review::assembly_product_view::homology_processing_view(step.step_index,evidence),      collagene_core::verification::evidence::PhysicalProcessingEvidence::TopoActions(evidence)=>{let graph=&evidence.projection.physical.product;view!{<details><summary>{format!("Step {} · accounted chemical actions",step.step_index+1)}</summary><p>{format!("All {} source calls retained · {} components · {} remaining nicks · {} unpaired calls",evidence.projection.physical.source_bases,graph.components.len(),graph.nicks.len(),graph.unpaired.iter().map(|r|r.end-r.start).sum::<usize>())}</p><p class="helper">"Ideal hydrolysis is a declared model, not a named nuclease or host-repair guarantee. Complete original, main and excised component graphs remain in run JSON evidence."</p>{evidence.bridge_notes.iter().map(|s|view!{<p>{s.clone()}</p>}).collect_view()}</details>}.into_any()},     collagene_core::verification::evidence::PhysicalProcessingEvidence::OverlapExtension(geometry)=>{view!{<details><summary>{format!("Step {} · reciprocal extension",step.step_index+1)}</summary><p>{format!("Homology {} bp · terminal pairings {:?} · complete {} · unique within terminal model {} · {} comparisons",geometry.options.overlap,geometry.terminal_pairings.overlaps,geometry.terminal_pairings.complete,geometry.terminal_pairings.unique,geometry.terminal_pairings.comparisons)}</p></details>}.into_any()}     collagene_core::verification::evidence::PhysicalProcessingEvidence::LigationOutcomes(evidence)=>crate::cloning::preparation::ligation_outcomes_view::ligation_outcomes_view(step.step_index,evidence), collagene_core::verification::evidence::PhysicalProcessingEvidence::DerivedOrigin(evidence)=>{view!{<details><summary>"Reusable material derivation replay"</summary><p>{format!("Exact event {} · inherited incomplete {} · {} retained findings",evidence.origin.product.event_id,evidence.material.origin_incomplete,evidence.material.origin_findings_total)}</p><p class="helper">"The retained target-free PCR/digest recipe was replayed; imported reference or end summaries were not stock."</p></details>}.into_any()}     collagene_core::verification::evidence::PhysicalProcessingEvidence::RunProduct(evidence)=>{let origin=evidence.decode_report().ok();view!{<details><summary>{format!("Step {} · saved material reconstruction",step.step_index+1)}</summary><p>{format!("Origin run {} · product {}",evidence.origin.selection.run_id,evidence.origin.selection.step_index+1)}</p><p>{origin.as_ref().map(|r|format!("Origin assessment: {:?} · {}",r.status,r.engine_version))}</p><p class="helper">"Already-prepared material is reconstructed against its original target. Nicks, termini, positional modifications and incomplete findings remain. The complete origin report is retained in reconstruction JSON."</p>{origin.into_iter().flat_map(|r|r.findings.into_iter().take(100)).map(|finding|view!{<p>{finding}</p>}).collect_view()}</details>}.into_any()}     }}).collect_view()}{report.step_reports.iter().filter_map(|step|step.topo.as_ref().and_then(|topo|topo.annotations.as_ref()).map(|evidence|(step.step_index,evidence))).map(|(step_index,evidence)|view!{<details class="annotation-evidence"><summary>{format!("Step {} · source annotation correspondence",step_index+1)}</summary><p class="helper">"Mapped coordinates follow the actual selected top-reference calls. Review PCR-derived and disrupted annotations. Original unmapped features remain evidence in the reconstruction JSON, not invented product coordinates."</p>{evidence.sources.iter().map(|source|view!{<p>{format!("{} · {} projected · {} outside contour · {} unmapped originals · projection supported {}",source.role,source.projected,source.dropped,source.unmapped.len(),source.projection_supported)}</p>}).collect_view()}</details>}).collect_view()}{report.step_reports.iter().filter_map(|step|step.topo.as_ref().map(|topo|{let outcome=match &topo.projection.outcome{collagene_core::cloning::topo::physical::projection::Outcome::Proven(contour)=>format!("Paired contour proven; {} nicks retained. No repair is inferred.",contour.nicks.len()),collagene_core::cloning::topo::physical::projection::Outcome::Withheld(_)=>"No unique paired contour; physical graph retained and reference withheld.".into()};view!{<details class="physical-evidence"><summary>{format!("Step {} · TOPO physical graph",step.step_index+1)}</summary><p>{outcome}</p><p>{format!("{} physical components · batch {}",topo.projection.physical.product.components.len(),topo.batch_revision)}</p><p class="helper">"Download the reconstruction report for the full physical graph, original source maps and projection evidence. Saved evidence needs fresh replay before use."</p></details>}})).collect_view()}{report.products.iter().enumerate().filter_map(|(i,p)|p.metadata.get("duplex_model").and_then(|s|serde_json::from_str::<chemistry::DuplexSummary>(s).ok()).map(|summary|view!{<details class="physical-evidence"><summary>{format!("Step {} · {}",i+1,p.name)}</summary><pre>{summary.describe()}</pre></details>})).collect_view()}{report.products.iter().enumerate().filter_map(|(i,p)|crate::cloning::methods::checked_cloning::modification_summary(p).map(|summary|view!{<details class="checked-modification-summary"><summary>{format!("Step {} · positional modification evidence",i+1)}</summary><p>{summary}</p></details>})).collect_view()}<details><summary>"Assumptions and limits"</summary><ul>{report.assumptions.iter().cloned().map(|a|view!{<li>{a}</li>}).collect_view()}</ul></details><p class="helper">"Alternative assemblies and experimental feasibility are not exhaustively checked."</p></div>}))}<Show when=move||reviewed_product_current.get()>{move||last_plan.get().zip(plan_report.get()).map(|((plan,target),report)|view!{<div><crate::assembly::review::assembly_product_view::AssemblyProductView plan target report epoch=workspace_epoch guard=Some(reviewed_product_guard) on_products=Callback::new(move|products:Vec<Molecule>|{reviewed_product_guard.run(())?;session.add_documents(products)}) on_finalize=Some(finalize_reviewed_candidate)/><crate::protocol::protocol_quantities_view::ReactionQuantitiesView capture=quantity_binding identity=quantity_identity epoch=workspace_epoch snapshot=contig_snapshot validate=validate_workspace_frame on_restore=restore_workspace_frame db workspace_busy=workspace_restoring save_generation message/></div>})}</Show><button class="button secondary full" disabled=move||plan_report.with(Option::is_none) on:click=move|_|{let result=plan_report.with_untracked(|report|report.as_ref().map(crate::cloning::methods::checked_cloning::prepare_report_export));if let Some(result)=result{crate::io::download::complete_download(result,move|result|{if message.try_get_untracked().is_some(){if let Err(error)=result{message.set(error);}}});}}>"Download reconstruction report"</button><button class="button secondary full" disabled=move||candidate.get().is_none()&&last_plan.get().is_none() on:click=move|_|{let value=candidate.get_untracked().or_else(||last_plan.get_untracked().map(|(plan,_)|transfer::Candidate::new(plan)));if let Some(value)=value {match transfer::save_candidate(&value){Ok(bytes)=>{crate::io::download::download_with_completion("candidate.cg2plan","application/octet-stream",&bytes,move|result|{if message.try_get_untracked().is_none(){return;}if let Err(e)=result {message.set(e);}})},Err(e)=>message.set(e.to_string())}}}>"Save candidate"</button><button class="button secondary full" disabled=move||last_plan.get().is_none() on:click=move|_|{modal_error.set(String::new());modal.set("record_run");}>"Keep reconstructed products"</button></div></Show>
                                            <Show when=move||panel.get()=="primers"><crate::primers::primers::PrimerPanel on_select_binding=navigate_primer_binding resolve_source=resolve_retained_source on_collection=Callback::new(move|table:oligo_formats::OligoTable|import_oligos.run((table,Vec::new()))) molecule=molecule selection=selection seed=primer_seed source_binding=primer_context pair_draft pair_selection anchor_available=primer_anchor_available primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer epoch=workspace_epoch on_use=Callback::new(move|(forward,reverse):(Primer,Primer)|{pcr_forward_name.set(forward.name);pcr_reverse_name.set(reverse.name);pcr_forward.set(forward.annealing.as_str().into());pcr_forward_tail.set(forward.tail);pcr_reverse.set(reverse.annealing.as_str().into());pcr_reverse_tail.set(reverse.tail);panel.set("pcr");message.set("Candidate primers loaded. Compute PCR products to inspect alternative binding geometries.".into());}) on_annotate=Callback::new(move|(binding,feature):(collagene_core::material::lineage::EventBinding,Feature)|{if primer_context.get_untracked()!=binding{return Err("The document event changed; rerun primer analysis before annotating it.".into());}let current=molecule.get_untracked();let edited=editing::set_feature(&current,None,feature).map_err(|e|e.to_string())?;session.commit_document(&binding.revision,edited)})/></Show>
                                            <div style:display=move||{if panel.get()=="gel"{"contents"}else{"none"}}><crate::restriction::gel::GelPanel seeds=gel_digest_seeds epoch=workspace_epoch resolve_event=resolve_retained_source molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone())))/></div>
                                            <Show when=move||panel.get()=="structures"><crate::biology::structure::structure_view::StructurePanel epoch=workspace_epoch locked=Memo::new(move|_|workspace_restoring.get())
                                                entries=Memo::new(move|_|attachments.with(|assets|assets.iter().filter_map(|(hash,a)|collagene_core::protein::structure::source_format(&a.kind).map(|f|(hash.clone(),a.name.clone(),f))).collect()))
                                                get_source=Callback::new(move|hash:String|->Result<crate::biology::structure::structure_view::SourceFile,String>{
                                                    if workspace_restoring.get_untracked(){return Err("Wait for workspace restoration.".into());}
                                                    attachments.with_untracked(|assets|{let a=assets.get(&hash).ok_or("Structure source is no longer retained.")?;let format=collagene_core::protein::structure::source_format(&a.kind).ok_or("Attachment is not a structure source.")?;
                                                        collagene_core::protein::structure::revision(&a.bytes,format,&None).map_err(|e|e.to_string())?;
                                                        if blake3::hash(&a.bytes).to_hex().as_str()!=hash{return Err("Structure attachment hash changed.".into());}
                                                        Ok(crate::biology::structure::structure_view::SourceFile{name:a.name.clone(),bytes:a.bytes.clone(),format})})
                                                })
                                                retain=Callback::new(move|(asset,epoch):(crate::biology::structure::structure_view::SourceFile,u64)|->Result<String,String>{
                                                    if workspace_restoring.get_untracked()||workspace_epoch.get_untracked()!=epoch{return Err("Workspace changed before retaining structure source.".into());}
                                                    attachments.try_update(|assets|collagene_core::protein::structure::retain(assets,&asset.name,asset.bytes,asset.format).map_err(|e|e.to_string())).ok_or("Workspace is unavailable.")?
                                                })
                                                remove=Callback::new(move|(hash,epoch):(String,u64)|->Result<(),String>{
                                                    if workspace_restoring.get_untracked()||workspace_epoch.get_untracked()!=epoch{return Err("Workspace changed before removing structure source.".into());}
                                                    attachments.try_update(|assets|{let a=assets.get(&hash).ok_or("Structure source is no longer retained.")?;if collagene_core::protein::structure::source_format(&a.kind).is_none(){return Err("Attachment is not a structure source.".into());}assets.remove(&hash);Ok(())}).ok_or("Workspace is unavailable.")?
                                                })
                                            /></Show>
                                            <Show when=move||panel.get()=="methods"><crate::cloning::preparation::cloning_methods::CloningMethodsPanel molecule=molecule get_workspace=Some(contig_snapshot) selected_region=Callback::new(move|():()|{if inspected_feature.get_untracked().is_some_and(|f|f.location.spans().len()>1){return None;}selection.get_untracked().filter(|(a,b)|a<b).map(|(a,b)|(active.get_untracked(),a,b))}) target_info=Memo::new(move|_|crate::cloning::preparation::review::TargetContext::Document(documents.with(|docs|recombination_source_info(&docs[active.get()])))) source_info=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(recombination_source_info))) primer_sources=primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) config=method_config epoch=workspace_epoch on_candidate=Callback::new(move|value:transfer::Candidate|{
                                            if workspace_restoring.get_untracked(){return;}
                                            let state=method_config.get_untracked();

                                            let captured=documents.with_untracked(|docs|->Result<_,String>{
                                                let indices=state.source_indices()?;
                                                let target=docs.get(active.get_untracked()).map(current_event_binding).ok_or("No independent target is active")?;
                                                if target.revision!=value.plan.target_revision||indices.len()!=value.plan.inputs.len(){return Err("Cloning source/target choices changed; reconstruct the method evidence again.".into());}
                                                let bindings=value.plan.inputs.iter().zip(indices).map(|(input,index)|{
                                                    let doc=docs.get(index).ok_or("A selected source document is unavailable")?;
                                                    let event=current_event_binding(doc);
                                                    if event.revision!=input.revision{return Err("A source event changed; reconstruct the method evidence again.".to_string());}
                                                    Ok(collagene_core::material::lineage::InputEventBinding{input_id:input.id.clone(),event})
                                                }).collect::<Result<Vec<_>,String>>()?;
                                                Ok((target,bindings))
                                            });
                                            let (target,bindings)=match captured{Ok(value)=>value,Err(error)=>{message.set(error);return;}};
                                            let preferred=(value.plan_revision.clone(),workspace_epoch.get_untracked(),target,bindings);
                                            candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("Method candidate opened with exact source-event bindings for independent target verification. No experimental success or mixture uniqueness is established.".into());
                                        }) on_products=Callback::new(move|products:Vec<Molecule>|session.add_documents(products))/></Show>
                                            <Show when=move||panel.get()=="ends"><crate::cloning::preparation::end_repair::EndRepairPanel epoch=workspace_epoch get_workspace=library_snapshot on_material=open_material_derivation source_info=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(recombination_source_info))) molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,doc)|(i,doc.molecule.name.clone())).collect())) get_query=Callback::new(move|index:usize|documents.with(|docs|docs.get(index).map(|doc|doc.molecule.clone()))) on_candidate=Callback::new(move|value:transfer::Candidate|{if workspace_restoring.get_untracked(){return;}candidate_events.set(None);candidate.set(Some(value));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);plan_generation.update(|generation|*generation=generation.wrapping_add(1));panel.set("plan");message.set("Physical candidate opened for review. Verify it against the independently selected active target; alternatives and chemistry remain explicit.".into());}) on_products=Callback::new(move|products:Vec<Molecule>|session.add_documents(products))/></Show>
                                            <Show when=move||panel.get()=="domestication"><crate::cloning::methods::domestication_view::DomesticationView source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) molecule selected_feature epoch=workspace_epoch get_workspace=contig_snapshot on_open=open_domestication on_navigate=Callback::new(move|(start,end):(usize,usize)|->Result<(),String>{
                                                let n=molecule.with_untracked(|m|m.sequence.len());
                                                if workspace_restoring.get_untracked()||start>=n||end<=start||end-start>n{return Err("The selected source range is no longer available.".into());}
                                                let spans=if end>n{vec![(start,n),(0,end-n)]}else{vec![(start,end)]};
                                                batch(||{selected_feature.set(None);selection.set(if spans.len()==1{spans.first().copied()}else{None});search_regions.set(Some((molecule.with_untracked(Molecule::revision),spans)));});
                                                Ok(())
                                            })/></Show>
                                            <Show when=move||panel.get()=="mutagenesis"><crate::cloning::methods::mutagenesis::MutagenesisPanel get_workspace=library_snapshot on_material=open_material_derivation resolve_source=resolve_retained_source on_overlay=Some(publish_pcr_overlay) on_navigate=Some(navigate_pcr_overlay) molecule=molecule mutation_selection mutation_seed on_intent=open_mutation_intent pair_draft pair_selection anchor_available=primer_anchor_available primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer epoch=workspace_epoch source_info=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(recombination_source_info))) target_binding=primer_context query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) on_candidate=Callback::new(move|opened:crate::cloning::methods::mutagenesis::CandidateOpen|{
                                                let result=(||->Result<_,String>{
                                                    if workspace_restoring.get_untracked() || opened.epoch!=workspace_epoch.get_untracked(){return Err("The mutation workspace changed; review the current design again.".into());}
                                                    let frame=documents.with_untracked(|docs|WorkspaceSession::snapshot(docs,active.get_untracked(),candidate.get_untracked().as_ref(),&stored_candidates.get_untracked(),&attachments.get_untracked(),&workspace_name.get_untracked()))?;
                                                    crate::assembly::review::assembly_review::stage_existing(&frame,&opened.candidate,&opened.target,&opened.sources)?;
                                                    Ok(())
                                                })();
                                                if let Err(error)=result{message.set(error);return;}
                                                let preferred=(opened.candidate.plan_revision.clone(),opened.epoch,opened.target,opened.sources);
                                                batch(||{candidate_events.set(None);candidate.set(Some(opened.candidate));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");});
                                                message.set("Inverse-PCR candidate opened with its original source and intended-target events. Verify, inspect the reconstructed product, then choose outputs.".into());
                                            }) on_outputs=open_mismatch_outputs/></Show>

                                            <Show when=move||panel.get()=="topo"><crate::cloning::methods::topo_view::TopoPanel
                                                target=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()])))
                                                sources=Memo::new(move|_|documents.with(|docs|docs.iter().map(recombination_source_info).collect()))
                                                capture=Callback::new(move|capture:crate::cloning::methods::topo_preparation::CaptureRequest|->Result<crate::cloning::methods::topo_preparation::Captured,String>{
                                                    if workspace_restoring.get_untracked(){return Err("Wait for the workspace operation to finish.".into());}
                                                    let frame=library_snapshot.run(())?;
                                                    crate::cloning::methods::topo_preparation::capture(&frame,&capture)
                                                }) get_workspace=library_snapshot epoch=workspace_epoch locked=Memo::new(move|_|workspace_restoring.get())
                                                validate=validate_workspace_frame on_restore=restore_workspace_frame db workspace_busy=workspace_restoring save_generation message
                                                on_candidate=Callback::new(move|(value,bindings):(transfer::Candidate,Vec<collagene_core::material::lineage::InputEventBinding>)|{
                                                    if workspace_restoring.get_untracked(){return;}
                                                    let target=documents.with_untracked(|docs|docs.get(active.get_untracked()).map(current_event_binding));let Some(target)=target else{return;};
                                                    if value.plan.target_revision!=target.revision||bindings.len()!=value.plan.inputs.len()||!documents.with_untracked(|docs|bindings.iter().all(|binding|value.plan.inputs.iter().any(|input|input.id==binding.input_id&&input.revision==binding.event.revision)&&docs.iter().any(|doc|current_event_binding(doc)==binding.event))){message.set("TOPO source/target event bindings changed. Inspect again.".into());return;}
                                                    let preferred=(value.plan_revision.clone(),workspace_epoch.get_untracked(),target,bindings);
                                                    candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("TOPO candidate opened for source-bound review and independent verification. Full physical evidence remains distinct from reference DNA.".into());
                                                })/></Show>
                                            <Show when=move||panel.get()=="recombination"><details><summary>"Multi-parent reference exchange"</summary><crate::cloning::methods::multisite::MultiSitePanel
                                                target=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()])))
                                                sources=Memo::new(move|_|documents.with(|docs|docs.iter().map(recombination_source_info).collect()))
                                                capture=Callback::new(move|capture:crate::cloning::methods::multisite::CaptureRequest| {
                                                    if workspace_restoring.get_untracked(){return Err("Wait for the workspace operation to finish.".into());}
                                                    documents.with_untracked(|docs| {
                                                        let resolve=|binding:&collagene_core::material::lineage::EventBinding|docs.iter().find(|d|d.id==binding.document_id&&current_event_binding(d)==*binding).ok_or_else(||"The selected document event changed. Inspect again.".to_string());
                                                        if !(2..=5).contains(&capture.sources.len()){return Err("Choose two to five original parents.".into());}
                                                        let parents=capture.sources.iter().map(resolve).collect::<Result<Vec<_>,_>>()?;
                                                        let target=resolve(&capture.target)?;
                                                        if docs.get(active.get_untracked()).map(current_event_binding).as_ref()!=Some(&capture.target){return Err("The active independent target changed.".into());}
                                                        let molecules=parents.iter().map(|d|&d.molecule).collect::<Vec<_>>();
                                                        collagene_core::cloning::recombination::multisite::source_preflight(&molecules,&capture.request).map_err(|e|e.to_string())?;
                                                        if target.molecule.sequence.len()>2_000_000{return Err("Independent target exceeds 2,000,000 bases.".into());}
                                                        collagene_model::bounds::model_preflight(&[&target.molecule]).map_err(|e|e.to_string())?;
                                                        target.molecule.validate().map_err(|e|e.to_string())?;
                                                        let inputs=parents.iter().enumerate().map(|(i,d)|collagene_core::cloning::recombination::reference::Input{role:format!("parent_{}",i+1),molecule:d.molecule.clone()}).collect();
                                                        Ok(crate::cloning::methods::multisite::Captured{inputs,target:target.molecule.clone()})
                                                    })
                                                }) epoch=workspace_epoch locked=Memo::new(move|_|workspace_restoring.get())
                                                on_candidate=Callback::new(move|(value,bindings):(transfer::Candidate,Vec<collagene_core::material::lineage::InputEventBinding>)| {
                                                    if workspace_restoring.get_untracked(){return;}
                                                    let Some(target)=documents.with_untracked(|docs|docs.get(active.get_untracked()).map(current_event_binding)) else{return;};
                                                    if value.plan.target_revision!=target.revision||bindings.len()!=value.plan.inputs.len()||!documents.with_untracked(|docs|bindings.iter().all(|binding|value.plan.inputs.iter().any(|input|input.id==binding.input_id&&input.revision==binding.event.revision)&&docs.iter().any(|doc|current_event_binding(doc)==binding.event))){message.set("Multi-parent source or target event bindings changed. Inspect again.".into());return;}
                                                    let preferred=(value.plan_revision.clone(),workspace_epoch.get_untracked(),target,bindings);
                                                    candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("Multi-parent reference candidate opened. Verify against the independent target; reference ancestry remains chemically restricted.".into());
                                                })/></details></Show>
                                            <Show when=move||panel.get()=="recombination"><crate::cloning::methods::recombination::RecombinationPanel
                                                target=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()])))
                                                sources=Memo::new(move|_|documents.with(|docs|docs.iter().map(recombination_source_info).collect()))
                                                capture=Callback::new(move|capture:crate::cloning::methods::recombination::CaptureRequest| {
                                                    if workspace_restoring.get_untracked(){return Err("Wait for the workspace operation to finish.".into());}
                                                    documents.with_untracked(|docs| {
                                                        let resolve=|binding:&collagene_core::material::lineage::EventBinding|docs.iter().find(|d|d.id==binding.document_id && current_event_binding(d)==*binding).ok_or_else(||"The selected document event changed. Inspect again.".to_string());
                                                        let a=resolve(&capture.sources[0])?;
                                                        let b=resolve(&capture.sources[1])?;
                                                        let target=resolve(&capture.target)?;
                                                        if docs.get(active.get_untracked()).map(current_event_binding).as_ref()!=Some(&capture.target){return Err("The active independent target changed.".into());}
                                                        collagene_core::cloning::recombination::source_target_preflight([&a.molecule,&b.molecule],&capture.request,&target.molecule).map_err(|e|e.to_string())?;
                                                        let inputs=verification::sources::recombination_inputs(&[Source::Input("source_1".into()),Source::Input("source_2".into())],[&a.molecule,&b.molecule],&capture.request).map_err(|e|e.to_string())?;
                                                        Ok(crate::cloning::methods::recombination::Captured{inputs,target:target.molecule.clone()})
                                                    })
                                                }) epoch=workspace_epoch locked=Memo::new(move|_|workspace_restoring.get())
                                                on_candidate=Callback::new(move|(value,bindings):(transfer::Candidate,Vec<collagene_core::material::lineage::InputEventBinding>)| {
                                                    if workspace_restoring.get_untracked(){return;}
                                                    let target=documents.with_untracked(|docs|docs.get(active.get_untracked()).map(current_event_binding));
                                                    let Some(target)=target else{return;};
                                                    if value.plan.target_revision!=target.revision || bindings.len()!=2 || !documents.with_untracked(|docs|bindings.iter().all(|b|docs.iter().any(|d|current_event_binding(d)==b.event))){message.set("Reference-exchange event bindings changed. Inspect again.".into());return;}
                                                    let preferred=(value.plan_revision.clone(),workspace_epoch.get_untracked(),target,bindings);
                                                    candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");message.set("Reference-only candidate opened for review. Verify against the independent target before recording source-event ancestry.".into());
                                                })/></Show>
                                            <Show when=move||panel.get()=="methylation"><crate::biology::methylation::MethylationPanel source=Memo::new(move|_|documents.with(|docs|{let d=&docs[active.get()];crate::biology::methylation::SourceInfo{binding:crate::reads::contigs::ReadBinding{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:d.molecule.revision()},name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet()}})) get_molecule=Callback::new(move|()|molecule.get_untracked()) epoch=workspace_epoch on_navigate=navigate_read_call on_digest=Callback::new(move|()|panel.set("enzymes"))/></Show>
                                            <Show when=move||panel.get()=="genomes"><crate::genomes::genomes::GenomePanel links=genome_links epoch=workspace_epoch on_link=Callback::new(move|link:collagene_core::tracks::genome::ReferenceLink|->Result<(),String>{
                                                if workspace_restoring.get_untracked(){return Err("Wait for workspace restoration.".into());}
                                                let(hash,attachment)=collagene_core::tracks::genome::make_attachment(&link).map_err(|e|e.to_string())?;
                                                let mut assets=attachments.get_untracked();assets.retain(|_,a|a.kind!=collagene_core::tracks::genome::WORKSPACE_ATTACHMENT_KIND||collagene_core::tracks::genome::load_link(&a.bytes).is_ok_and(|old|old.role!=link.role));assets.insert(hash,attachment);
                                                let frame=documents.with_untracked(|docs|WorkspaceSession::snapshot(docs,active.get_untracked(),candidate.get_untracked().as_ref(),&stored_candidates.get_untracked(),&assets,&workspace_name.get_untracked()))?;frame.validate().map_err(|e|e.to_string())?;attachments.set(assets);Ok(())
                                            }) on_unlink=Callback::new(move|role:String|->Result<(),String>{if workspace_restoring.get_untracked(){return Err("Wait for workspace restoration.".into());}attachments.update(|assets|assets.retain(|_,a|a.kind!=collagene_core::tracks::genome::WORKSPACE_ATTACHMENT_KIND||collagene_core::tracks::genome::load_link(&a.bytes).is_ok_and(|l|l.role!=role)));Ok(())})/></Show>
                                            <Show when=move||panel.get()=="crispr"><crate::cloning::methods::crispr::CrisprPanel source=Memo::new(move|_|documents.with(|docs|{let d=&docs[active.get()];crate::cloning::methods::crispr::DocumentInfo{binding:crate::reads::contigs::ReadBinding{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:d.molecule.revision()},name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet(),topology:d.molecule.topology}})) documents=Memo::new(move|_|documents.with(|docs|docs.iter().map(|d|crate::cloning::methods::crispr::DocumentInfo{binding:crate::reads::contigs::ReadBinding{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:d.molecule.revision()},name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet(),topology:d.molecule.topology}).collect())) get_document=Callback::new(move|id:String|documents.with_untracked(|docs|docs.iter().find(|d|d.id==id).map(|d|d.molecule.clone()))) epoch=workspace_epoch on_navigate=navigate_read_call on_record=record_donor/></Show>
                                            <Show when=move||panel.get()=="contigs"><crate::reads::contigs::ContigPanel read_options=Memo::new(move|_|documents.with(|docs|docs.iter().map(|d|crate::reads::contigs::ReadInfo{document_id:d.id.clone(),event_id:d.head_event().to_string(),revision:d.molecule.revision(),name:d.molecule.name.clone(),length:d.molecule.sequence.len(),eligible:d.molecule.sequence.alphabet()==Alphabet::Dna&&d.molecule.topology==Topology::Linear,quality_label:if d.molecule.metadata.contains_key("trace_original_base64"){"Sanger trace · native caller evidence".into()}else if d.molecule.metadata.contains_key(collagene_core::formats::oligos::FASTQ_ORIGINAL_KEY){format!("FASTQ · {}",d.molecule.metadata.get("fastq:encoding").map(String::as_str).unwrap_or("declared encoding"))}else{"No quality evidence".into()}}).collect())) get_read=Callback::new(move|id:String|documents.with_untracked(|docs|docs.iter().find(|d|d.id==id).map(|d|d.molecule.clone()))) epoch=workspace_epoch on_record=record_assembly on_navigate=navigate_read_call/></Show>
                                            <Show when=move||panel.get()=="codons"><crate::biology::codon_design::CodonDesignPanel molecule=molecule source_identity=Memo::new(move|_|documents.with(|docs|format!("{}:{}",docs[active.get()].id,docs[active.get()].head_event()))) epoch=workspace_epoch on_open=Callback::new(move|product:Molecule|{
                                                if workspace_restoring.get_untracked(){return Err("Wait for the workspace update to finish.".into());}
                                                let mut docs=documents.get_untracked();let first=docs.len();
                                                if first>=1000||docs.iter().map(|d|document_bytes(&d.molecule)).sum::<usize>()+document_bytes(&product)>MAX_WORKSPACE_BYTES{return Err("This derived document exceeds the workspace document or memory limit.".into());}
                                                let source=docs.get(active.get_untracked()).ok_or("No active source document")?;
                                                let prefix=format!("document-{}",&product.revision()[..16]);let mut suffix=1;
                                                let id=loop{let id=format!("{prefix}-{suffix}");if docs.iter().all(|d|d.id!=id){break id;}suffix+=1;};
                                                let derived=source.derive_event(source.head_event(),id,product,"Opened codon candidate after independent translation and constraint checks").map_err(|e|e.to_string())?;
                                                docs.push(derived);
                                                let frame=WorkspaceSession::snapshot(&docs,first,candidate.get_untracked().as_ref(),&stored_candidates.get_untracked(),&attachments.get_untracked(),&workspace_name.get_untracked())?;
                                                frame.validate().map_err(|e|e.to_string())?;
                                                batch(||{documents.set(docs);active.set(first);});Ok(())
                                            })/></Show>
                                            <Show when=move||panel.get()=="biology"><Show when=move||molecule.with(|m|m.sequence.alphabet()!=Alphabet::Protein)><crate::biology::biology::BiologyPanel source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) epoch=workspace_epoch selection_location=biology_selection_location on_select_location=select_biology_location on_open_protein=open_translated_protein on_retain_cds=retain_translation_cds molecule=molecule selected_feature=selected_feature selection=selection on_annotate=Callback::new(move|(base,feature):(String,Feature)|{let current=molecule.get_untracked();if current.revision()!=base{return Err("The document changed; rerun ORF analysis before annotating it.".into());}let edited=editing::set_feature(&current,None,feature).map_err(|e|e.to_string())?;session.commit_document(&base,edited)})/></Show><crate::biology::biology_source_view::TranslationSourceLink molecule on_open=open_translation_origin/><Show when=move||molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein)><crate::biology::protein_properties::ProteinPanel molecule=molecule/><crate::biology::antibody::antibody_view::AntibodyPanel source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) epoch=workspace_epoch workspace_busy=workspace_restoring prepare=prepare_antibody on_select=select_antibody_residue on_overlay=publish_antibody/><crate::biology::protein_liability_view::ProteinLiabilityPanel source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) epoch=workspace_epoch workspace_busy=workspace_restoring prepare=prepare_protein_liabilities on_select=select_protein_liability on_annotate=annotate_protein_liability selection=Memo::new(move|_|selection.get())/></Show><Show when=move||molecule.with(|m|m.sequence.alphabet()==Alphabet::Rna)><crate::biology::rna_view::RnaPanel active_source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) epoch=workspace_epoch workspace_busy=workspace_restoring resolve=resolve_rna saved=saved_rna on_save=save_rna on_load=load_rna on_delete=delete_rna/></Show></Show>
                                            <Show when=move||panel.get()=="dotplot"><crate::viewer::alignment::dotplot::DotPlotPanel molecule=molecule active_index=Memo::new(move|_|active.get()) source_identity=Memo::new(move|_|documents.with(|docs|format!("{}:{}",docs[active.get()].id,docs[active.get()].head_event()))) query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,format!("{} [{}]",d.molecule.name,d.id))).collect())) query_info=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|crate::viewer::alignment::dotplot::QueryInfo{identity:format!("{}:{}",d.id,d.head_event()),revision:d.molecule.revision(),name:d.molecule.name.clone(),length:d.molecule.sequence.len(),alphabet:d.molecule.sequence.alphabet()}))) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) epoch=workspace_epoch on_navigate=Callback::new(move|(query,index,spans):(bool,usize,Vec<(usize,usize)>)|{
                                                if workspace_restoring.get_untracked(){return;}
                                                let index=if query{index}else{active.get_untracked()};
                                                let Some((document_id,event_id,revision))=documents.with_untracked(|docs|docs.get(index).map(|d|(d.id.clone(),d.head_event().to_string(),d.molecule.revision())))else{return;};
                                                let expected_epoch=workspace_epoch.get_untracked();
                                                let apply=move|deferred:bool|{
                                                    if workspace_epoch.try_get_untracked()!=Some(expected_epoch)||workspace_restoring.try_get_untracked().unwrap_or(true){return;}
                                                    let Some(current_index)=active.try_get_untracked()else{return;};
                                                    if documents.with_untracked(|docs|docs.iter().position(|d|d.id==document_id&&d.head_event()==event_id))!=Some(current_index){return;}
                                                    if deferred&&(selected_feature.get_untracked().is_some()||selection.get_untracked().is_some()||search_regions.get_untracked().is_some()){return;}
                                                    batch(||{selected_feature.set(None);selection.set(if spans.len()==1{spans.first().copied()}else{None});search_regions.set(Some((revision,spans)));});
                                                };
                                                if active.get_untracked()==index&&has_open_document.get_untracked(){apply(false);}else{
                                                    active.set(index);
                                                    // App clears old selections when the active source changes.
                                                    // Install this source-bound selection after those effects settle.
                                                    let frame=wasm_bindgen::closure::Closure::once_into_js(move||apply(true));
                                                    if let Some(window)=web_sys::window(){let _=window.request_animation_frame(frame.unchecked_ref());}
                                                }
                                                if query{panel.set("features");}
                                            })/></Show><Show when=move||panel.get()=="compare"><crate::viewer::alignment::compare::ComparePanel molecule=molecule query_options=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().map(|(i,d)|(i,d.molecule.name.clone())).collect())) get_query=Callback::new(move|i:usize|documents.with(|docs|docs.get(i).map(|d|d.molecule.clone()))) active_source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()]))) source_info=Callback::new(move|index:usize|documents.with(|docs|docs.get(index).map(recombination_source_info))) resolve=resolve_retained_source epoch=workspace_epoch prepare_similarity=prepare_similarity guard_similarity=guard_similarity saved=saved_comparisons on_save=save_comparison on_reopen=reopen_comparison on_navigate=navigate_comparison/></Show>
                                        </div></aside>
                                    </div>
                                    <Show when=move||assembly_visible.get()><section class="assembly-dock" aria-label="Assembly workspace"><div class="assembly-dock-bar"><span>"Assembly design"</span><button class="text-tool" on:click=move|_|assembly_visible.set(false)>"Hide assembly"</button></div>
                                        <fieldset aria-label="Assembly source and preparation controls" disabled=move||workspace_restoring.get() style="border:0;padding:0;margin:0;min-width:0">
                                        <crate::assembly::draft::assembly_workspace::AssemblyWorkspace web_mcp=web_mcp_assembly draft=assembly_draft previous_draft=previous_assembly_draft epoch=workspace_epoch snapshot=contig_snapshot workspace_busy=workspace_restoring
                                            apply_saved_attachments=Callback::new(move|(before,after):(workspace::Workspace,workspace::Workspace)|{
                                                if workspace_restoring.get_untracked(){return Err("Wait for the workspace update to finish.".to_owned());}
                                                let current=contig_snapshot.run(())?;
                                                crate::assembly::library::assembly_saved_designs_commit::validate(&before,&after,&current)?;
                                                attachments.set(after.attachments);
                                                Ok(())
                                            })
                                            sources=Memo::new(move|_|documents.with(|docs|docs.iter().map(recombination_source_info).collect()))
                                            active_source=Memo::new(move|_|documents.with(|docs|recombination_source_info(&docs[active.get()])))
                                            selected_region=Callback::new(move|():()|{if inspected_feature.get_untracked().is_some_and(|f|f.location.spans().len()>1){return None;}selection.get_untracked()})
                                            get_source=Callback::new(move|binding:collagene_core::material::lineage::EventBinding|documents.with_untracked(|docs|{
                                                let doc=docs.iter().find(|d|d.id==binding.document_id).ok_or("The selected source document is unavailable.")?;
                                                let exact=doc.revision_molecule(&binding.event_id).map_err(|e|e.to_string())?;
                                                if exact.revision()!=binding.revision{return Err("The selected source event changed.".into());}
                                                Ok(exact)
                                            }))
                                            on_navigate=navigate_assembly on_review=review_assembly on_finalize=finalize_assembly
                                            on_products=Callback::new(move|products:Vec<Molecule>|session.add_documents(products))
                                            primer_sources=primer_sources list_primers=list_saved_primers resolve_primer=resolve_saved_primer/>
                                        </fieldset>
                                        <crate::assembly::library::assembly_library_view::AssemblyLibraryView draft=assembly_draft epoch=workspace_epoch snapshot=contig_snapshot
                                            sources=Memo::new(move|_|documents.with(|docs|docs.iter().map(recombination_source_info).collect()))
                                            get_source=Callback::new(move|binding:collagene_core::material::lineage::EventBinding|documents.with_untracked(|docs|{
                                                let doc=docs.iter().find(|doc|doc.id==binding.document_id).ok_or("The source document is unavailable.")?;
                                                let source=doc.revision_molecule(&binding.event_id).map_err(|error|error.to_string())?;
                                                if source.revision()!=binding.revision{return Err("The selected source event changed.".into());}
                                                Ok(source)
                                            }))
                                            validate=validate_workspace_frame
                                            on_restore=restore_library_frame
                                            db workspace_busy=workspace_restoring loaded save_generation message/>
                                        <crate::assembly::draft::assembly_pool_view::AssemblyPoolView draft=assembly_draft snapshot=contig_snapshot epoch=workspace_epoch busy=workspace_restoring on_review=review_assembly on_finalize=finalize_assembly/>
                                        <crate::assembly::library::assembly_batch_view::AssemblyBatchView draft=assembly_draft snapshot=contig_snapshot
                                            validate=validate_workspace_frame
                                            on_restore=restore_workspace_frame
                                            attachment_available=Callback::new(move|hash:String|attachments.with_untracked(|assets|assets.get(&hash).is_some_and(|asset|asset.kind==collagene_core::planning::combinatorial::batch::KIND)))
                                            db workspace_busy=workspace_restoring loaded save_generation epoch=workspace_epoch message on_release=open_bulk_release/>
                                        <crate::assembly::review::assembly_archive_view::AssemblyArchiveView snapshot=contig_snapshot epoch=workspace_epoch selected=selected_release_manifest busy=workspace_restoring
                                            attachment_available=Callback::new(move|hash:String|attachments.with_untracked(|assets|assets.get(&hash).is_some_and(|asset|asset.kind==collagene_core::planning::combinatorial::release::KIND)))/>
                                    </section></Show>
                                    <footer class="statusbar"><span>{move||save_status.get()}</span><Show when=move||!layout.status.get().is_empty()><span role="status">{move||layout.status.get()}</span></Show><span>{move||selection_status.get()}</span><span>{move||documents.with(|docs|format!("{} {}",docs.len(),if docs.len()==1{"document"}else{"documents"}))}</span></footer>
                                    </div>
                                </main>
                            </div>
                            </div>
                            <Show when=move||modal.get()!="">{view!{<div class="modal-backdrop"><section class="modal" class:alignment-modal=move||modal.get()=="alignments" class:read-review-modal=move||matches!(modal.get(),"plan_edit"|"plan_new")||modal.get()=="combinatorial_release"||modal.get()=="read_review"||modal.get()=="run_event"||modal.get()=="comparison_source"||modal.get()=="translation_source"||modal.get()=="pcr_source" role="dialog" aria-modal="true" aria-labelledby="dialog-title"><div class="modal-heading"><h2 id="dialog-title">{move||if matches!(modal.get(),"plan_edit"|"plan_new"){"Edit cloning plan"}else if modal.get()=="variants"{"Variant source"}else if modal.get()=="read_alignments"{"Read alignments"}else if modal.get()=="material_product"{"Open reusable PCR product"}else if modal.get()=="digest_product"{"Open virtual digest fragment"}else if modal.get()=="crispr_product"{"Open checked sequence design"}else if modal.get()=="contig_archive"{"Record read assembly"}else if modal.get()=="saved_contig"{"Saved read assembly"}else if modal.get()=="read_review"{"Review sequencing differences"}else if modal.get()=="read_review_commit"{"Save read-review decisions"}else if modal.get()=="read_check_capture"{"Save sequencing check"}else if modal.get()=="combinatorial_release"{"Create selected construct outputs"}else if modal.get()=="assembly_finalize"{"Create assembly outputs"}else if modal.get()=="record_run"{"Keep reconstructed products"}else if modal.get()=="run_details"{"Recorded reconstruction"}else if modal.get()=="run_event"{"Recorded sequence version"}else if modal.get()=="comparison_source"{"Comparison source version"}else if modal.get()=="translation_source"{"Original translation source"}else if modal.get()=="pcr_source"{"PCR template placements"}else if modal.get()=="reference_paste"{"Paste annotated reference"}else if modal.get()=="alignments"{"Alignments"}else if modal.get()=="image_evidence"{"Documentary images"}else if modal.get()=="workspace"{"Workspace"}else if modal.get()=="oligos"{"Oligo table"}else if modal.get()=="fastq"{"Import FASTQ reads"}else if modal.get()=="import"{"Open a sequence"}else if modal.get()=="url"{"Import from a sequence source"}else if modal.get()=="feature"{"Edit annotation"}else if modal.get()=="compat_import"{"Review sequence import"}else if modal.get()=="compat_snapgene"{"Export SnapGene sequence"}else if modal.get()=="compat_xdna"{"Export XDNA sequence"}else if modal.get()=="compat_gcg"{"Export GCG sequence"}else if modal.get()=="protein_genpept"{"Export GenPept protein"}else if modal.get()=="protein_uniprot"{"Export UniProt protein"}else if modal.get()=="sbol"{"SBOL design exchange"}else if modal.get()=="document"{"Document settings"}else if molecule.with(|m|m.sequence.alphabet()==Alphabet::Protein){"Edit selected residues"}else{"Edit selected bases"}}</h2><button class="icon-button" aria-label="Close dialog" disabled=move||workspace_restoring.get() on:click=move|_|{if !workspace_restoring.get_untracked(){modal.set("");pending_finalize.set(None);pending_release.set(None);pending_run_event.set(None);pending_comparison_source.set(None);pending_reference_paste.set(None);pending_read_check.set(None);pending_read_review.set(None);reviewed_check.set(None);review_seed.set(None);}}>"×"</button></div><Show when=move||!modal_error.get().is_empty()><p class="dialog-error" role="alert">{move||modal_error.get()}</p></Show><Show when=move||modal.get()=="material_product">{move||pending_material.get_untracked().map(|request|view!{<crate::restriction::digest_product::MaterialOpenPanel request snapshot=contig_snapshot validate=validate_workspace_frame on_restore=restore_workspace_frame db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="digest_product">{move||pending_digest.get_untracked().map(|request|view!{<crate::restriction::digest_product::DigestOpenPanel request=request snapshot=contig_snapshot validate=validate_workspace_frame on_restore=contig_restore db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="crispr_product">{move||pending_donor.get_untracked().map(|request|view!{<crate::cloning::methods::crispr::DonorOpenPanel request=request snapshot=contig_snapshot validate=validate_workspace_frame on_restore=contig_restore db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="contig_archive">{move||pending_assembly.get_untracked().map(|request|view!{<crate::reads::contigs::ContigArchivePanel request=request snapshot=contig_snapshot validate=validate_workspace_frame on_restore=contig_restore db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||modal.get()=="saved_contig"><crate::reads::contigs::SavedAssemblyPanel attachment_hash=inspected_assembly.get_untracked() snapshot=contig_snapshot epoch=workspace_epoch on_record=record_assembly on_navigate=navigate_read_call/></Show><Show when=move||modal.get()=="record_run">{move||last_plan.get_untracked().map(|(plan,_)|view!{<crate::history::run_capture::RunCapturePanel preferred_target=last_plan_target.get_untracked().filter(|(revision,epoch,_)|*epoch==workspace_epoch.get_untracked()&&*revision==collagene_core::planning::plan::DesignPlan::revision(&plan)).map(|(_,_,target)|target) plan=plan preferred_bindings=last_plan_events.get_untracked().filter(|(revision,_)|last_plan.get_untracked().is_some_and(|(p,_)|collagene_core::planning::plan::DesignPlan::revision(&p)==*revision)).map(|(_,bindings)|bindings).unwrap_or_default() snapshot=Callback::new(move|()|documents.with_untracked(|docs|stored_candidates.with_untracked(|plans|attachments.with_untracked(|assets|candidate.with_untracked(|current|WorkspaceSession::snapshot(docs,active.get_untracked(),current.as_ref(),plans,assets,&workspace_name.get_untracked())))))) validate=validate_workspace_frame on_restore=restore_workspace_frame db=db busy=workspace_restoring loaded=loaded save_generation=save_generation epoch=workspace_epoch message=message/>})}</Show><Show when=move||matches!(modal.get(),"plan_edit"|"plan_new")><crate::protocol::plan_editing_host::ManualPlanEditor candidate preferences=candidate_events documents active epoch=workspace_epoch locked=Memo::new(move|_|workspace_restoring.get()) creating=modal.get_untracked()=="plan_new" on_close=Callback::new(move|()|modal.set("")) on_apply=Callback::new(move|(value,target,bindings):(transfer::Candidate,collagene_core::material::lineage::EventBinding,Vec<collagene_core::material::lineage::InputEventBinding>)|->Result<(),String>{
                                    if workspace_restoring.get_untracked(){return Err("Wait for workspace restoration.".into());}
                                    let preferred=(value.plan_revision.clone(),workspace_epoch.get_untracked(),target,bindings);
                                    batch(||{
                                        plan_generation.update(|g|*g=g.wrapping_add(1));
                                        plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);
                                        pending_finalize.set(None);candidate_events.set(None);candidate.set(Some(value));candidate_events.set(Some(preferred));
                                        modal.set("");message.set("Edited candidate saved for review. Verify it again before creating products or protocols.".into());
                                    });Ok(())
                                })/></Show><Show when=move||modal.get()=="reference_paste">{move||pending_reference_paste.get_untracked().map(|opened|view!{<crate::workspace::reference_transfer_view::ReferencePastePanel opened guard=reference_paste_guard on_apply=apply_reference_paste on_cancel=Callback::new(move|()|{pending_reference_paste.set(None);modal.set("");})/>})}</Show><Show when=move||modal.get()=="pcr_source">{move||pending_pcr_source.get().map(|snapshot|{
                                let owner=snapshot.owner.clone();let bound=snapshot.overlay.clone();
                                let current=Memo::new(move|_|workspace_epoch.get()==bound.epoch&&!workspace_restoring.get()
                                    &&published_pcr_overlay.with(|p|p.as_ref().is_some_and(|p|p.owner==owner&&p.overlay.as_ref().is_some_and(|o|std::sync::Arc::ptr_eq(o,&bound))))
                                    &&documents.with(|docs|docs.iter().any(|d|current_event_binding(d)==bound.source)));
                                view!{<crate::primers::pcr::pcr_source_view::PcrSourceView snapshot current on_back=Callback::new(move|()|{pending_pcr_source.set(None);modal.set("");})/>}
                            })}</Show><Show when=move||modal.get()=="translation_source">{move||pending_translation_source.get().map(|snapshot|view!{<crate::biology::biology_source_view::TranslationSourceView snapshot epoch=workspace_epoch on_back=Callback::new(move|()|{pending_translation_source.set(None);modal.set("");})/>})}</Show><Show when=move||modal.get()=="comparison_source">{move||pending_comparison_source.get().map(|snapshot|view!{<crate::viewer::alignment::pairwise_source_view::PairwiseSourceView snapshot epoch=workspace_epoch on_back=Callback::new(move|()|{pending_comparison_source.set(None);modal.set("");})/>})}</Show><Show when=move||modal.get()=="run_event">{move||pending_run_event.get().map(|snapshot|view!{<crate::history::run_event_view::RunEventPanel snapshot epoch=workspace_epoch on_fork=fork_run_event on_back=back_to_run/>})}</Show><Show when=move||modal.get()=="run_details"><crate::history::run_capture::RunDetailsPanel validate=validate_workspace_frame on_restore=restore_workspace_frame db workspace_busy=workspace_restoring save_generation message snapshot=Callback::new(move|()|documents.with_untracked(|docs|stored_candidates.with_untracked(|plans|attachments.with_untracked(|assets|candidate.with_untracked(|current|WorkspaceSession::snapshot(docs,active.get_untracked(),current.as_ref(),plans,assets,&workspace_name.get_untracked())))))) run_id=selected_run.get_untracked() epoch=workspace_epoch on_open_event=open_run_event on_edit_plan=Callback::new(move|prepared:std::sync::Arc<collagene_core::material::review::PreparedReview>|->Result<(),String>{
                if workspace_restoring.get_untracked(){return Err("Wait for workspace restoration.".into());}
                let frame=library_snapshot.run(())?;
                prepared.validate_for_workspace(&frame).map_err(|e|e.to_string())?;
                let preferred=(prepared.candidate.plan_revision.clone(),workspace_epoch.get_untracked(),prepared.target_event.clone(),prepared.inputs.clone());
                batch(||{candidate.set(Some(prepared.candidate.clone()));candidate_events.set(Some(preferred));plan_report.set(None);last_plan.set(None);last_plan_target.set(None);last_plan_events.set(None);plan_generation.update(|g|*g=g.wrapping_add(1));panel.set("plan");modal_error.set(String::new());modal.set("plan_edit");});
                message.set("Editing the retained plan against its exact independent target and source events. Viewer heads are unchanged; verify edits afresh.".into());
                Ok(())
            })/></Show><Show when=move||modal.get()=="workspace"><crate::workspace::workspace_view::WorkspacePanel snapshot=Callback::new(move|()|documents.with_untracked(|docs|stored_candidates.with_untracked(|plans|attachments.with_untracked(|assets|candidate.with_untracked(|current|WorkspaceSession::snapshot(docs,active.get_untracked(),current.as_ref(),plans,assets,&workspace_name.get_untracked())))))) validate=validate_workspace_frame on_restore=restore_workspace_frame db=db busy=workspace_restoring loaded=loaded save_generation=save_generation message=message/></Show>
                                <Show when=move||modal.get()=="read_review">{move||reviewed_check.get().map(|hash|view!{<crate::reads::read_review_view::ReadReviewPanel attachment_hash=hash reference=read_check_reference snapshot=contig_snapshot epoch=workspace_epoch on_commit=save_read_review on_read=navigate_read_call on_reference=navigate_check_reference event_available=primer_anchor_available attachment_available=read_check_attachment_available initial_request=review_seed.get_untracked()/>})}</Show>
                                <Show when=move||modal.get()=="read_review_commit">{move||pending_read_review.get().map(|pending|view!{<crate::reads::read_review_capture::ReadReviewCapturePanel request=pending.request snapshot=read_review_snapshot validate=validate_workspace_frame on_restore=restore_read_review db busy=workspace_restoring loaded save_generation epoch=workspace_epoch message/>})}</Show>
                                <Show when=move||modal.get()=="read_check_capture">{move||pending_read_check.get().map(|request|view!{<crate::reads::read_checks_view::ReadCheckCapturePanel request snapshot=read_check_snapshot validate=validate_workspace_frame on_restore=restore_read_check db busy=workspace_restoring loaded save_generation epoch=workspace_epoch message/>})}</Show>
                                <Show when=move||modal.get()=="combinatorial_release">{move||pending_release.get().map(|hash|view!{<crate::assembly::review::assembly_release_view::AssemblyReleaseView batch_hash=hash draft=assembly_draft snapshot=contig_snapshot
                                    validate=validate_workspace_frame
                                    on_restore=restore_workspace_frame
                                    db workspace_busy=workspace_restoring loaded save_generation epoch=workspace_epoch message on_success=released_bulk/>})}</Show>
                                <Show when=move||modal.get()=="assembly_finalize">{move||pending_finalize.get().map(|pending|view!{<crate::assembly::review::assembly_finalize::AssemblyFinalizePanel target=pending.target.clone() plan=pending.plan.clone() preferred_bindings=pending.sources.clone() snapshot=finalize_snapshot recovery_snapshot=contig_snapshot validate=validate_workspace_frame on_restore=restore_workspace_frame db busy=workspace_restoring loaded save_generation epoch=workspace_epoch message/>})}</Show>
                                <Show when=move||modal.get()=="alignments"><crate::viewer::alignment::alignment_view::AlignmentPanel alignments=alignments pending=pending_alignment.get_untracked() epoch=workspace_epoch on_keep={let opened_epoch=workspace_epoch.get_untracked();Callback::new(move|item:crate::viewer::alignment::alignment_view::AlignmentRevision|{if opened_epoch!=workspace_epoch.get_untracked(){return Err("The workspace changed during alignment import. Open the alignment again.".into());}keep_alignment.run(item)})} on_extract={let opened_epoch=workspace_epoch.get_untracked();Callback::new(move|molecules:Vec<Molecule>|{if opened_epoch!=workspace_epoch.get_untracked(){return Err("The workspace changed during row extraction.".into());}session.add_documents(molecules)})}/></Show>
                                <Show when=move||modal.get()=="oligos"><crate::io::oligo_import::OligoImport on_import=import_oligos existing=existing_oligo.get_untracked() initial_text=pending_oligo.get_untracked()/></Show>
                                <Show when=move||modal.get()=="fastq"><p>"Select the quality encoding supplied by the sequencing source. Phred+33 and Phred+64 overlap; it cannot be inferred from printable characters."</p><label>"Supplied FASTQ encoding"<select aria-label="Supplied FASTQ encoding" prop:value=move||match fastq_encoding.get(){None=>"",Some(oligo_formats::QualityEncoding::Phred33)=>"33",Some(oligo_formats::QualityEncoding::Phred64)=>"64"} on:change=move|e|fastq_encoding.set(match event_target_value(&e).as_str(){"33"=>Some(oligo_formats::QualityEncoding::Phred33),"64"=>Some(oligo_formats::QualityEncoding::Phred64),_=>None})><option value="">"Choose supplied encoding"</option><option value="33">"Phred+33"</option><option value="64">"Phred+64"</option></select></label><button class="button primary" disabled=move||fastq_encoding.get().is_none() on:click=move|_|{let(Some(bytes),Some(encoding))=(pending_fastq.get_untracked(),fastq_encoding.get_untracked())else{return;};match oligo_formats::import_fastq(&bytes,encoding){Ok(report)=>match session.add_documents(report.molecules){Ok(())=>{message.set(report.warnings.join(" "));pending_fastq.set(None);modal.set("");panel.set("reads");},Err(e)=>modal_error.set(e)},Err(e)=>modal_error.set(e.to_string())}}>"Import FASTQ with selected encoding"</button></Show>
                                <Show when=move||modal.get()=="variants"><p class="helper">"Declared file alleles and genotype columns are source evidence. Choose an actual DNA reference and a one-based alternative explicitly before checking or creating a design."</p><button class="button secondary" disabled=move||workspace_restoring.get() on:click=move|_|if let Err(error)=keep_variant_source.run(()){modal_error.set(error)}>"Keep original source in workspace"</button>{move||pending_variants.get().map(|(source,format,name)|view!{<crate::genomes::variants_view::VariantPanel source=source format=format name=name documents=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().filter(|(_,doc)|doc.molecule.sequence.alphabet()==Alphabet::Dna).map(|(index,doc)|(index,doc.molecule.name.clone(),format!("{}:{}:{}",workspace_epoch.get(),doc.id,doc.head_event()))).collect())) get_document=Callback::new(move|index:usize|documents.with_untracked(|docs|docs.get(index).and_then(|doc|collagene_model::bounds::model_preflight(&[&doc.molecule]).ok().map(|_|doc.molecule.clone())))) on_navigate=Callback::new(move|(index,parts):(usize,Vec<(usize,usize)>)|{if workspace_restoring.get_untracked(){return;}let Some((id,event,revision))=documents.with_untracked(|docs|docs.get(index).map(|doc|(doc.id.clone(),doc.head_event().to_string(),doc.molecule.revision())))else{return;};let epoch=workspace_epoch.get_untracked();let apply=move|deferred:bool|{if workspace_epoch.try_get_untracked()!=Some(epoch)||workspace_restoring.try_get_untracked().unwrap_or(true){return;}if !documents.with_untracked(|docs|docs.get(active.get_untracked()).is_some_and(|doc|doc.id==id&&doc.head_event()==event)){return;}if deferred&&(selection.get_untracked().is_some()||selected_feature.get_untracked().is_some()||search_regions.get_untracked().is_some()){return;}batch(||{selected_feature.set(None);selection.set(if parts.len()==1{parts.first().copied()}else{None});search_regions.set(Some((revision,parts)));});};modal.set("");if active.get_untracked()==index&&has_open_document.get_untracked(){apply(false);}else{active.set(index);let frame=wasm_bindgen::closure::Closure::once_into_js(move||apply(true));if let Some(window)=web_sys::window(){let _=window.request_animation_frame(frame.unchecked_ref());}}}) on_preview=open_variant_preview/>})}</Show><Show when=move||modal.get()=="read_alignments"><p class="helper">"Keep the original source in this workspace to reopen it later. Mapping claims remain unverified until checked against an explicitly selected DNA document."</p><button class="button secondary" disabled=move||workspace_restoring.get() on:click=move|_|if let Err(error)=keep_read_source.run(()){modal_error.set(error)}>"Keep original source in workspace"</button>{move||pending_read_alignment.get().map(|(source,format,name)|view!{
                                    <crate::reads::read_alignments_view::ReadAlignmentPanel source=source format=format name=name documents=Memo::new(move|_|documents.with(|docs|docs.iter().enumerate().filter(|(_,doc)|doc.molecule.sequence.alphabet()==Alphabet::Dna).map(|(index,doc)|(index,doc.molecule.name.clone(),format!("{}:{}:{}",workspace_epoch.get(),doc.id,doc.head_event()))).collect())) get_document=Callback::new(move|index:usize|documents.with_untracked(|docs|docs.get(index).and_then(|doc|collagene_model::bounds::model_preflight(&[&doc.molecule]).ok().map(|_|doc.molecule.clone())))) on_navigate=Callback::new(move|(index,parts):(usize,Vec<(usize,usize)>)|{if workspace_restoring.get_untracked(){return;}let Some((id,event,revision))=documents.with_untracked(|docs|docs.get(index).map(|doc|(doc.id.clone(),doc.head_event().to_string(),doc.molecule.revision())))else{return;};
                                        let epoch=workspace_epoch.get_untracked();
                                        let apply=move|deferred:bool|{if workspace_epoch.try_get_untracked()!=Some(epoch)||workspace_restoring.try_get_untracked().unwrap_or(true){return;}if !documents.with_untracked(|docs|docs.get(active.get_untracked()).is_some_and(|doc|doc.id==id&&doc.head_event()==event)){return;}if deferred&&(selection.get_untracked().is_some()||selected_feature.get_untracked().is_some()||search_regions.get_untracked().is_some()){return;}batch(||{selected_feature.set(None);selection.set(if parts.len()==1{parts.first().copied()}else{None});search_regions.set(Some((revision,parts)));});};
                                        modal.set("");if active.get_untracked()==index&&has_open_document.get_untracked(){apply(false);}else{active.set(index);let frame=wasm_bindgen::closure::Closure::once_into_js(move||apply(true));if let Some(window)=web_sys::window(){let _=window.request_animation_frame(frame.unchecked_ref());}}})/>
                                })}</Show><Show when=move||modal.get()=="import"><p>"Paste GenBank, GenPept, UniProt, EMBL, FASTA, or plain DNA. No project setup needed."</p><button class="button secondary" disabled=move||workspace_restoring.get() on:click=move|_|{let text=import_text.get_untracked();if !text.trim_start_matches('\u{feff}').trim_start().starts_with('>'){modal_error.set("Paste FASTA records beginning with > before choosing their alphabet.".into());return;}pending_compat_import.set(Some((std::sync::Arc::new(text.into_bytes()),"Pasted FASTA".into(),crate::io::sequence_compat_import::Profile::Fasta(None))));modal_error.set(String::new());modal.set("compat_import");}>"Choose FASTA alphabet"</button><textarea aria-label="Sequence input" rows="12" prop:value=move||import_text.get() on:input=move|e|import_text.set(event_target_value(&e)) placeholder=">my_sequence\nACGT…"></textarea><button class="button primary" on:click=move|_|{let text=import_text.get_untracked();if text.starts_with("##fileformat=VCF"){if text.len()>16*1024*1024{modal_error.set("Variant source exceeds the 16 MiB browser profile.".into());return;}pending_variants.set(Some((std::sync::Arc::new(text.into_bytes()),collagene_core::review::variants::session::FileFormat::Vcf,"pasted.vcf".into())));modal_error.set(String::new());modal.set("variants");return;}if [b"@HD\t".as_slice(),b"@SQ\t".as_slice(),b"@RG\t".as_slice(),b"@PG\t".as_slice(),b"@CO\t".as_slice()].iter().any(|prefix|text.as_bytes().starts_with(prefix)){if text.len()>collagene_core::review::reads::session::MAX_SOURCE_BYTES{modal_error.set("This browser read-alignment profile accepts sources up to 16 MiB.".into());return;}pending_read_alignment.set(Some((std::sync::Arc::new(text.into_bytes()),collagene_core::review::reads::session::FileFormat::Sam,"pasted.sam".into())));modal_error.set(String::new());modal.set("read_alignments");return;}if text.as_bytes().first()==Some(&b'@'){pending_fastq.set(Some(std::sync::Arc::new(text.into_bytes())));fastq_encoding.set(None);modal.set("fastq");return;}match formats::import_text(&text){Ok(report)=>{if let Err(error)=session.add_documents(report.molecules){modal_error.set(error);return;}message.set(report.warnings.join(" "));modal.set("");import_text.set(String::new());},Err(e)=>modal_error.set(e.to_string())}}>"Open sequence"</button></Show><Show when=move||modal.get()=="url"><crate::io::source_import::SourceImport on_import={let epoch=workspace_epoch.get_untracked();Callback::new(move |report:formats::ImportReport|{if epoch!=workspace_epoch.get_untracked(){return Err("The workspace changed during retrieval. Import the source again in this workspace.".into());}session.add_documents(report.molecules)?;message.set(report.warnings.join(" "));modal.set("");Ok(())})}/></Show><Show when=move||modal.get()=="feature"><crate::viewer::sequence::editors::FeatureEditor molecule=molecule.get_untracked() index=feature_for_editor.get_untracked() selection=selection.get_untracked() initial_location=feature_initial_location.get_untracked() on_commit=Callback::new(move|(base,edited):(String,Molecule)|{session.commit_document(&base,edited)?;modal.set("");Ok(())})/></Show><Show when=move||matches!(modal.get(),"protein_genpept"|"protein_uniprot")>{move||{let profile=if modal.get()=="protein_genpept"{collagene_core::formats::protein::Profile::GenPept}else{collagene_core::formats::protein::Profile::UniProt};let guard=session.export_context_guard();view!{<crate::biology::protein_export::ProteinExportDialog source=molecule.get_untracked() profile guard on_saved=Callback::new(move|notice:String|{if message.try_get_untracked().is_some(){message.set(notice);}})/>}}}</Show>
        <Show when=move||modal.get()=="compat_import">{move||pending_compat_import.get().map(|(source,name,profile)|{
            let selected=source.clone();let context=session.export_context_guard();
            let guard=Callback::new(move|()|{context.run(())?;if modal.get_untracked()!="compat_import"{return Err("The import review closed.".into());}if !pending_compat_import.with_untracked(|p|p.as_ref().is_some_and(|(current,_,_)|std::sync::Arc::ptr_eq(current,&selected))){return Err("The selected import source changed; reopen its review.".into());}Ok(())});
            view!{<crate::io::sequence_compat_import::SequenceCompatImport source name profile guard on_import=Callback::new(move|report:formats::ImportReport|{let notice=report.warnings.join(" ");session.add_documents(report.molecules)?;message.set(notice);modal.set("");pending_compat_import.set(None);Ok(())})/>}
        })}</Show>
        <Show when=move||modal.get()=="image_evidence">{move||{
            let guard=session.export_context_guard();
            let snapshot=Callback::new(move|()|->Result<crate::io::image_evidence_model::Snapshot,String>{
                guard.run(())?;
                if modal.try_get_untracked()!=Some("image_evidence")||!loaded.get_untracked(){return Err("The image view closed or its workspace is unavailable.".into());}
                crate::io::image_evidence_model::Snapshot::new(library_snapshot.run(())?,workspace_epoch.get_untracked())
            });
            let on_attach=Callback::new(move|request:crate::io::image_evidence_model::Install|->crate::io::image_evidence::AttachResult {
                let proposed=(||{
                    request.ensure_live()?;
                    let current=snapshot.run(())?;
                    request.context.matches(&current.context)?;
                    let summary=request.prepared.summary().map_err(|e|e.to_string())?;
                    if summary.record.association!=current.context.association{return Err("The prepared image belongs to another selected document.".into());}
                    let attached=request.prepared.apply_to(&current.workspace).map_err(|e|e.to_string())?;
                    if attached.source_workspace_fingerprint!=request.context.workspace_fingerprint{return Err("The workspace changed before image attachment.".into());}
                    Ok::<_,String>(attached.workspace)
                })();
                Box::pin(async move {
                    request.ensure_live()?;
                    let frame=proposed?;
                    // The archive's 64MiB encoded limit is independent of its expanded/model
                    // limits. Check representability in the private persistence worker first.
                    let encoded=crate::workspace::session::encode_workspace(frame,save_generation.get_untracked() as u32).await?;
                    drop(encoded);
                    request.ensure_live()?;
                    let current=snapshot.run(())?;
                    request.context.matches(&current.context)?;
                    // Rebuild only from the unchanged current frame. No retained attachment
                    // clone can overwrite another edit while archive encoding was pending.
                    let attached=request.prepared.apply_to(&current.workspace).map_err(|e|e.to_string())?;
                    if attached.source_workspace_fingerprint!=request.context.workspace_fingerprint{return Err("The workspace changed before image attachment admission.".into());}
                    attachments.set(attached.workspace.attachments);
                    // attachments is observed by the existing generation-guarded autosave lane.
                    // No document, feature/origin, history or candidate is replaced here.
                    modal.set("");
                    Ok(attached.summary)
                })
            });
            let on_saved:std::sync::Arc<dyn Fn(String)+Send+Sync>=std::sync::Arc::new(move|notice:String|{if message.try_get_untracked().is_some(){message.set(notice);}});
            view!{<crate::io::image_evidence::ImageEvidenceDialog snapshot on_attach on_saved/>}
        }}</Show>
        <Show when=move||modal.get()=="map_export">{move||{
            let context=session.export_context_guard();
            let snapshot=Callback::new(move|()|->Result<crate::io::map_export_model::Snapshot,String>{
                context.run(())?;
                if modal.try_get_untracked()!=Some("map_export") || !loaded.get_untracked() {
                    return Err("The map export view closed or its workspace is unavailable.".into());
                }
                let map=renderer.get_value().ok_or("Wait for the source map to finish opening.")?.export_view()?;
                let (source,event)=documents.with_untracked(|docs|docs.get(active.get_untracked()).map(|doc|(doc.molecule.clone(),current_event_binding(doc)))).ok_or("The selected map document is unavailable.")?;
                let hidden=annotation_hidden.get_untracked().into_iter().collect::<Vec<_>>();
                if map.source_revision!=event.revision || map.interaction_context!=restriction_identity.get_untracked()
                    || map.hidden_features!=hidden || map.selected_feature!=selected_feature.get_untracked()
                    || map.highlighted_spans!=highlighted_spans.get_untracked() {
                    return Err("Map source or display settings are still updating; prepare the export again.".into());
                }
                let restriction_layer=if restrictions.visible.get_untracked(){
                    if restrictions.busy.get_untracked(){return Err("Wait for the displayed restriction scan before exporting it.".into());}
                    let scan=restrictions.scan.get_untracked().ok_or("No exact restriction scan is ready; scan or hide restriction sites before exporting.")?;
                    if scan.revision!=event.revision || scan.context!=restriction_identity.get_untracked(){return Err("Restriction scan belongs to a different source event.".into());}
                    let sites=restrictions.displayed.get_untracked().ok_or("Displayed restriction sites are still updating.")?;
                    if sites.as_ref()!=map.restriction_sites.as_ref(){return Err("Restriction map filters are still updating; prepare the export again.".into());}
                    Some(crate::io::map_export_model::RestrictionLayer{
                        source_revision:event.revision.clone(),
                        display_context:{let settings=serde_json::to_vec(&(scan.context.as_str(),restrictions.names.get_untracked(),restrictions.count_filter.get_untracked(),restrictions.region_filter.get_untracked(),format!("{:?}",restrictions.basis.get_untracked()))).map_err(|error|error.to_string())?;blake3::hash(&settings).to_hex().to_string()},
                        sites:sites.as_ref().clone(),selected_site:map.selected_site,
                    })
                }else{
                    if !map.restriction_sites.is_empty(){return Err("Restriction visibility is still updating; prepare the export again.".into());}
                    None
                };
                let view=crate::io::map_export_model::ViewSettings{
                    selection:selection.get_untracked().map(|(start,end)|crate::io::map_export_model::Selection{start,end,wrap:start>end}),
                    rotation_radians:map.rotation_radians,zoom:map.zoom,selected_feature:map.selected_feature,
                    highlighted_spans:map.highlighted_spans,hidden_features:map.hidden_features,
                    display_context:serde_json::to_string(&(map.interaction_context,editor_view.get_untracked(),map.overlay_context)).map_err(|error|error.to_string())?,
                };
                Ok(crate::io::map_export_model::Snapshot{source,event,epoch:workspace_epoch.get_untracked(),view,restrictions:restriction_layer})
            });
            view!{<crate::io::map_export::MapExportDialog snapshot render=crate::io::map_export::worker_renderer() on_saved=Callback::new(move|notice:String|{if message.try_get_untracked().is_some(){message.set(notice);}})/>}
        }}</Show>
        <Show when=move||matches!(modal.get(),"compat_snapgene"|"compat_xdna"|"compat_gcg")>{move||{let expected_modal=modal.get();let mode=match expected_modal{"compat_snapgene"=>crate::io::sequence_export_model::Mode::SnapGene,"compat_xdna"=>crate::io::sequence_export_model::Mode::Xdna,_=>crate::io::sequence_export_model::Mode::Gcg};let context=session.export_context_guard();let guard=Callback::new(move|()|{context.run(())?;if modal.try_get_untracked().as_ref()!=Some(&expected_modal){return Err("The export view closed or changed format.".into());}Ok(())});view!{<crate::io::sequence_compat_export::SequenceCompatExport source=molecule.get_untracked() mode guard on_saved=Callback::new(move|notice:String|{if message.try_get_untracked().is_some(){message.set(notice);}})/>}}}</Show>
        <Show when=move||modal.get()=="sbol">{move||{
            let source=pending_sbol.get();
            let expected_source=source.clone();
            let context_guard=session.export_context_guard();
            let guard=Callback::new(move|()|{context_guard.run(())?;let current=pending_sbol.try_get_untracked().ok_or("The SBOL source view closed.")?;let same=match(&expected_source,&current){(None,None)=>true,(Some(a),Some(b))=>std::sync::Arc::ptr_eq(a,b),_=>false};if !same{return Err("The selected SBOL source changed; reopen its review.".into());}Ok(())});
            let document_id=documents.with_untracked(|docs|docs.get(active.get_untracked()).map(|d|d.id.clone()).unwrap_or_else(||"no-active-document".into()));
            view!{<crate::io::sbol_exchange::SbolExchangeDialog bytes=source current=std::sync::Arc::new(molecule.get_untracked()) document_id epoch=workspace_epoch guard on_close=Callback::new(move|()|{pending_sbol.set(None);modal.set("");}) on_import=Callback::new(move|(product,attachment):(Molecule,Option<collagene_core::documents::workspace::WorkspaceAttachment>)|->Result<(),String>{
                guard.run(())?;
                let mut docs=documents.get_untracked();
                if docs.len()>=1000 || docs.iter().map(|d|document_bytes(&d.molecule)).sum::<usize>()+document_bytes(&product)>MAX_WORKSPACE_BYTES {return Err("This SBOL projection exceeds the workspace document/memory limit.".into());}
                let index=docs.len();let additions=open_histories(vec![product],&docs)?;docs.extend(additions);
                let mut assets=attachments.get_untracked();
                if let Some(asset)=attachment {
                    let hash=blake3::hash(&asset.bytes).to_hex().to_string();
                    if let Some(existing)=assets.get(&hash) {if existing.kind!=asset.kind || existing.media_type!=asset.media_type || existing.bytes!=asset.bytes {return Err("The SBOL source collides with an existing attachment; no workspace data changed.".into());}}
                    else {assets.insert(hash,asset);}
                }
                let frame=WorkspaceSession::snapshot(&docs,index,candidate.get_untracked().as_ref(),&stored_candidates.get_untracked(),&assets,&workspace_name.get_untracked())?;
                frame.validate().map_err(|e|e.to_string())?;
                batch(||{documents.set(docs);active.set(index);attachments.set(assets);pending_sbol.set(None);});
                message.set("Reviewed SBOL component imported with its original graph when available; source annotations remain documentary.".into());modal.set("");Ok(())
            })/>}
        }}</Show><Show when=move||modal.get()=="document"><crate::viewer::sequence::editors::DocumentEditor molecule=molecule.get_untracked() selection=selection.get_untracked() on_commit=Callback::new(move|(base,edited):(String,Molecule)|{session.commit_document(&base,edited)?;modal.set("");Ok(())})/></Show><Show when=move||modal.get()=="edit"><p>{move||selection.get().map(|(a,b)|if a==b{format!("Insert at boundary {a}. Affected annotations will be marked for review.")}else{format!("Replace {} {}–{}; clear the field to delete. Affected annotations will be marked for review.",molecule.with(|m|m.sequence.alphabet().selected_symbols()),a+1,b)}).unwrap_or_else(||format!("Select {} in the sequence view first.",molecule.with(|m|m.sequence.alphabet().selected_symbols())))}</p><textarea aria-label="Replacement sequence" rows="6" prop:value=move||edit_text.get() on:input=move|e|edit_text.set(event_target_value(&e))></textarea><button class="button primary" disabled=move||selection.get().is_none() on:click=move|_|if let Some((a,b))=selection.get_untracked(){match molecule.get_untracked().replace(a,b,&edit_text.get_untracked()){Ok(edited)=>{let base=molecule.get_untracked().revision();if let Err(error)=session.commit_document(&base,edited){modal_error.set(error);return;}selection.set(None);plan_report.set(None);last_plan.set(None);last_plan_target.set(None);modal.set("");},Err(e)=>modal_error.set(e.to_string())}}>"Apply edit"</button></Show></section></div>}.into_any()}</Show>
                        </div>
                    }
}
