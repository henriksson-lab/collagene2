#[cfg(target_arch = "wasm32")]
mod alignment_view;
#[cfg(target_arch = "wasm32")]
mod anneal_preparation_view;
#[cfg(any(target_arch = "wasm32", test))]
mod annotation_geometry;
#[cfg(any(target_arch = "wasm32", test))]
mod annotation_layout;
#[cfg(target_arch = "wasm32")]
mod annotation_table;
#[cfg(any(target_arch = "wasm32", test))]
mod annotation_visibility;
#[cfg(target_arch = "wasm32")]
mod app;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_alternatives;
#[cfg(target_arch = "wasm32")]
mod assembly_alternatives_view;
#[cfg(target_arch = "wasm32")]
mod assembly_archive_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_batch_projection;
#[cfg(target_arch = "wasm32")]
mod assembly_batch_view;
#[cfg(target_arch = "wasm32")]
mod assembly_draft_actions;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_draft_model;
#[cfg(target_arch = "wasm32")]
mod assembly_finalize;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_fusion_design;
#[cfg(target_arch = "wasm32")]
mod assembly_fusion_view;
#[cfg(target_arch = "wasm32")]
mod assembly_library_add_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_library_commit;
#[cfg(target_arch = "wasm32")]
mod assembly_library_history_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_library_topology;
#[cfg(target_arch = "wasm32")]
mod assembly_library_view;
#[cfg(target_arch = "wasm32")]
mod assembly_multisite_view;
#[cfg(target_arch = "wasm32")]
mod assembly_outputs;
#[cfg(target_arch = "wasm32")]
mod assembly_pool_view;
#[cfg(target_arch = "wasm32")]
mod assembly_product_view;
#[cfg(target_arch = "wasm32")]
mod assembly_profile_view;
#[cfg(target_arch = "wasm32")]
mod assembly_release_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_review;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_role_ids;
#[cfg(target_arch = "wasm32")]
mod assembly_strategy_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_substitution_model;
#[cfg(target_arch = "wasm32")]
mod assembly_substitution_view;
#[cfg(target_arch = "wasm32")]
mod assembly_workspace;
#[cfg(target_arch = "wasm32")]
mod checked_cloning;
#[cfg(target_arch = "wasm32")]
mod clipboard;
#[cfg(any(target_arch = "wasm32", test))]
mod clipboard_model;
#[cfg(target_arch = "wasm32")]
mod cloning_methods;
#[cfg(any(target_arch = "wasm32", test))]
mod cloning_methods_model;
#[cfg(any(target_arch = "wasm32", test))]
mod cloning_preparations;
#[cfg(target_arch = "wasm32")]
mod codon_design;
#[cfg(any(target_arch = "wasm32", test))]
mod commit_lease;
#[cfg(target_arch = "wasm32")]
mod contigs;
#[cfg(target_arch = "wasm32")]
mod crispr;
#[cfg(target_arch = "wasm32")]
mod desktop_shell;
#[cfg(target_arch = "wasm32")]
mod digest;
#[cfg(target_arch = "wasm32")]
mod digest_product;
#[cfg(target_arch = "wasm32")]
mod document_explorer;
#[cfg(any(target_arch = "wasm32", test))]
mod document_tabs_model;
#[cfg(any(target_arch = "wasm32", test))]
mod domestication_form;
#[cfg(any(target_arch = "wasm32", test))]
mod domestication_review;
#[cfg(target_arch = "wasm32")]
mod domestication_view;
#[cfg(target_arch = "wasm32")]
mod dotplot;
#[cfg(any(target_arch = "wasm32", test))]
mod end_geometry;
#[cfg(target_arch = "wasm32")]
mod end_repair;
#[cfg(target_arch = "wasm32")]
mod end_view;
#[cfg(target_arch = "wasm32")]
mod feature_tools;
#[cfg(target_arch = "wasm32")]
mod genome_store;
#[cfg(target_arch = "wasm32")]
mod genome_view;
#[cfg(target_arch = "wasm32")]
mod genomes;
#[cfg(target_arch = "wasm32")]
mod glyph_renderer;
#[cfg(target_arch = "wasm32")]
mod golden_gate_candidate_view;
#[cfg(any(target_arch = "wasm32", test))]
mod golden_gate_pool_form;
#[cfg(target_arch = "wasm32")]
mod history_view;
#[cfg(target_arch = "wasm32")]
mod literal_intent_view;
#[cfg(target_arch = "wasm32")]
mod map_labels;
#[cfg(target_arch = "wasm32")]
mod mismatch_pcr_outputs;
#[cfg(any(target_arch = "wasm32", test))]
mod multisite_form;
#[cfg(target_arch = "wasm32")]
mod oligo_import;
#[cfg(any(target_arch = "wasm32", test))]
mod pcr_digest_model;
#[cfg(target_arch = "wasm32")]
mod pcr_digest_view;
#[cfg(target_arch = "wasm32")]
mod protein_properties;
#[cfg(target_arch = "wasm32")]
mod read_alignments_view;
#[cfg(target_arch = "wasm32")]
mod read_checks_view;
#[cfg(target_arch = "wasm32")]
mod read_review_capture;
#[cfg(target_arch = "wasm32")]
mod read_review_history;
#[cfg(target_arch = "wasm32")]
mod read_review_view;
#[cfg(target_arch = "wasm32")]
mod read_source_library;
#[cfg(target_arch = "wasm32")]
mod reference_transfer_view;
#[cfg(target_arch = "wasm32")]
mod restriction_map;
#[cfg(target_arch = "wasm32")]
mod restriction_view;
#[cfg(target_arch = "wasm32")]
mod run_capture;
#[cfg(any(target_arch = "wasm32", test))]
mod run_event_model;
#[cfg(target_arch = "wasm32")]
mod run_event_view;
#[cfg(target_arch = "wasm32")]
mod search_view;
#[cfg(any(target_arch = "wasm32", test))]
mod selection_toolbar;
#[cfg(target_arch = "wasm32")]
mod sequence_view;
#[cfg(any(target_arch = "wasm32", test))]
mod source_identity;
#[cfg(target_arch = "wasm32")]
mod topo_actions_view;
#[cfg(target_arch = "wasm32")]
mod topo_nick_view;
#[cfg(any(target_arch = "wasm32", test))]
mod topo_preparation;
#[cfg(target_arch = "wasm32")]
mod topo_view;
#[cfg(target_arch = "wasm32")]
mod trace_validation;
#[cfg(target_arch = "wasm32")]
mod variant_source_library;
#[cfg(target_arch = "wasm32")]
mod variants_view;
#[cfg(target_arch = "wasm32")]
mod workspace_commit;
#[cfg(target_arch = "wasm32")]
mod workspace_layout;
#[cfg(any(target_arch = "wasm32", test))]
mod workspace_preferences;
#[cfg(target_arch = "wasm32")]
mod workspace_view;

#[cfg(target_arch = "wasm32")]
mod biology;
#[cfg(any(target_arch = "wasm32", test))]
mod biology_selection;
#[cfg(target_arch = "wasm32")]
mod biology_source_view;
#[cfg(target_arch = "wasm32")]
mod compare;
#[cfg(target_arch = "wasm32")]
mod editors;
#[cfg(target_arch = "wasm32")]
mod gel;
#[cfg(target_arch = "wasm32")]
mod golden_gate;
#[cfg(target_arch = "wasm32")]
mod methylation;
#[cfg(target_arch = "wasm32")]
mod multisite;
#[cfg(target_arch = "wasm32")]
mod mutagenesis;
#[cfg(any(target_arch = "wasm32", test))]
mod mutation_review;
#[cfg(target_arch = "wasm32")]
mod pairwise_context_view;
#[cfg(target_arch = "wasm32")]
mod pairwise_export_view;
#[cfg(any(target_arch = "wasm32", test))]
mod pairwise_pages;
#[cfg(target_arch = "wasm32")]
mod pairwise_result_view;
#[cfg(target_arch = "wasm32")]
mod pairwise_search_view;
#[cfg(target_arch = "wasm32")]
mod pairwise_source_view;
#[cfg(target_arch = "wasm32")]
mod pcr_outputs;
#[cfg(any(target_arch = "wasm32", test))]
mod pcr_overlay;
#[cfg(target_arch = "wasm32")]
mod pcr_overlay_draw;
#[cfg(target_arch = "wasm32")]
mod pcr_source_view;
#[cfg(target_arch = "wasm32")]
mod primer_binding_view;
#[cfg(any(target_arch = "wasm32", test))]
mod primer_domain_edit;
#[cfg(target_arch = "wasm32")]
mod primer_domain_edit_view;
#[cfg(target_arch = "wasm32")]
mod primer_orders;
#[cfg(target_arch = "wasm32")]
mod primer_pair_editor;
#[cfg(any(target_arch = "wasm32", test))]
mod primer_pair_model;
#[cfg(any(target_arch = "wasm32", test))]
mod primer_pair_storage;
#[cfg(target_arch = "wasm32")]
mod primer_picker;
#[cfg(target_arch = "wasm32")]
mod primer_reference_capture_view;
#[cfg(target_arch = "wasm32")]
mod primers;
#[cfg(target_arch = "wasm32")]
mod recombination;
#[cfg(target_arch = "wasm32")]
mod renderer;
#[cfg(target_arch = "wasm32")]
mod source_import;
#[cfg(target_arch = "wasm32")]
mod storage;
#[cfg(target_arch = "wasm32")]
mod trace_view;
#[cfg(target_arch = "wasm32")]
mod worker;

#[cfg(target_arch = "wasm32")]
fn main() {
    leptos::mount::mount_to_body(app::App);
}

#[cfg(not(target_arch = "wasm32"))]
fn main() {
    println!("Build the web application with: trunk serve crates/web/index.html");
}

#[cfg(target_arch = "wasm32")]
mod assembly_junctions_view;

#[cfg(target_arch = "wasm32")]
mod plan_editing_host;
#[cfg(target_arch = "wasm32")]
mod plan_editing_view;

#[cfg(any(target_arch = "wasm32", test))]
mod plan_editing_sources;

#[cfg(all(test, not(target_arch = "wasm32")))]
mod sequence_geometry;

#[cfg(any(target_arch = "wasm32", test))]
mod run_origin_model;
#[cfg(target_arch = "wasm32")]
mod run_origin_view;

#[cfg(any(target_arch = "wasm32", test))]
mod digest_fragment_model;
