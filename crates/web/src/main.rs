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
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_comparison_model;
#[cfg(target_arch = "wasm32")]
mod assembly_comparison_view;
#[cfg(target_arch = "wasm32")]
mod assembly_draft_actions;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_draft_model;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_edit_model;
#[cfg(target_arch = "wasm32")]
mod assembly_edit_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_edited_primer_commit;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_edited_primer_model;
#[cfg(target_arch = "wasm32")]
mod assembly_edited_primer_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_extension_primer_commit;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_extension_primer_model;
#[cfg(target_arch = "wasm32")]
mod assembly_extension_primer_view;
#[cfg(target_arch = "wasm32")]
mod assembly_finalize;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_fusion_design;
#[cfg(target_arch = "wasm32")]
mod assembly_fusion_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_inverse_pcr_model;
#[cfg(target_arch = "wasm32")]
mod assembly_inverse_pcr_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_junction_navigation_model;
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
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_output_model;
#[cfg(target_arch = "wasm32")]
mod assembly_outputs;
#[cfg(target_arch = "wasm32")]
mod assembly_pool_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_product_export_model;
#[cfg(target_arch = "wasm32")]
mod assembly_product_exports;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_product_model;
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
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_saved_designs_commit;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_saved_designs_model;
#[cfg(target_arch = "wasm32")]
mod assembly_saved_designs_view;
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
mod ligation_outcomes_view;
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
#[cfg(target_arch = "wasm32")]
mod overlap_extension_preparation_view;
#[cfg(any(target_arch = "wasm32", test))]
mod pcr_digest_model;
#[cfg(target_arch = "wasm32")]
mod pcr_digest_view;
#[cfg(target_arch = "wasm32")]
mod protein_export;
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
#[cfg(any(target_arch = "wasm32", test))]
mod restriction_label_groups;
#[cfg(target_arch = "wasm32")]
mod restriction_map;
#[cfg(target_arch = "wasm32")]
mod restriction_view;
#[cfg(any(target_arch = "wasm32", test))]
mod rna_model;
#[cfg(target_arch = "wasm32")]
mod rna_view;
#[cfg(target_arch = "wasm32")]
mod run_capture;
#[cfg(any(target_arch = "wasm32", test))]
mod run_event_model;
#[cfg(target_arch = "wasm32")]
mod run_event_view;
#[cfg(target_arch = "wasm32")]
mod sbol_exchange;
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
#[cfg(any(target_arch = "wasm32", test))]
mod topo_graph_model;
#[cfg(target_arch = "wasm32")]
mod topo_graph_view;
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
#[cfg(any(target_arch = "wasm32", test))]
mod web_mcp;
#[cfg(any(target_arch = "wasm32", test))]
mod web_mcp_assembly;
#[cfg(any(target_arch = "wasm32", test))]
mod web_mcp_design;
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
#[cfg(any(target_arch = "wasm32", test))]
mod catalog_backup_model;
#[cfg(target_arch = "wasm32")]
mod catalog_backup_view;
#[cfg(target_arch = "wasm32")]
mod compare;
#[cfg(target_arch = "wasm32")]
mod desktop_rpc;
#[cfg(target_arch = "wasm32")]
mod desktop_sync_client;
#[cfg(target_arch = "wasm32")]
mod editors;
#[cfg(target_arch = "wasm32")]
mod gel;
#[cfg(target_arch = "wasm32")]
mod methylation;
#[cfg(target_arch = "wasm32")]
mod multisite;
#[cfg(target_arch = "wasm32")]
mod mutagenesis;
#[cfg(any(target_arch = "wasm32", test))]
mod mutation_review;
#[cfg(any(target_arch = "wasm32", test))]
mod native_publication_model;
#[cfg(target_arch = "wasm32")]
mod native_publication_view;
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
mod remote_library_view;
#[cfg(any(target_arch = "wasm32", test))]
mod remote_removal_model;
#[cfg(target_arch = "wasm32")]
mod remote_removal_view;
#[cfg(target_arch = "wasm32")]
mod renderer;
#[cfg(any(target_arch = "wasm32", test))]
mod scf_export_model;
#[cfg(target_arch = "wasm32")]
mod scf_export_view;
#[cfg(target_arch = "wasm32")]
mod source_gateway_client;
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
mod plan_processing_model;

#[cfg(any(target_arch = "wasm32", test))]
mod plan_editing_sources;

#[cfg(any(target_arch = "wasm32", test))]
mod sequence_edit_transition;

#[cfg(all(test, not(target_arch = "wasm32")))]
mod sequence_geometry;

#[cfg(any(target_arch = "wasm32", test))]
mod run_origin_model;
#[cfg(target_arch = "wasm32")]
mod run_origin_view;

#[cfg(any(target_arch = "wasm32", test))]
mod digest_fragment_model;

#[cfg(target_arch = "wasm32")]
mod transient_menus;

#[cfg(target_arch = "wasm32")]
mod library_bindings;
#[cfg(any(target_arch = "wasm32", test))]
mod library_capture;
#[cfg(target_arch = "wasm32")]
mod library_client;
#[cfg(target_arch = "wasm32")]
mod library_documents;
#[cfg(target_arch = "wasm32")]
mod library_view;

#[cfg(any(target_arch = "wasm32", test))]
mod material_derivation_model;
#[cfg(target_arch = "wasm32")]
mod reusable_material_view;
#[cfg(target_arch = "wasm32")]
mod sequence_compat_export;
#[cfg(target_arch = "wasm32")]
mod sequence_compat_import;
#[cfg(any(target_arch = "wasm32", test))]
mod sequence_export_model;

#[cfg(target_arch = "wasm32")]
mod map_export;
#[cfg(any(target_arch = "wasm32", test))]
mod map_export_model;
#[cfg(any(target_arch = "wasm32", test))]
mod saved_run_review_model;
#[cfg(target_arch = "wasm32")]
mod saved_run_review_view;

#[cfg(any(target_arch = "wasm32", test))]
mod protocol_layout_model;
#[cfg(target_arch = "wasm32")]
mod protocol_layout_view;

#[cfg(any(target_arch = "wasm32", test))]
mod protocol_quantities_model;
#[cfg(target_arch = "wasm32")]
mod protocol_quantities_view;

#[cfg(any(target_arch = "wasm32", test))]
mod library_discussion_model;
#[cfg(target_arch = "wasm32")]
mod library_discussion_view;

#[cfg(target_arch = "wasm32")]
mod library_discussion_client;

#[cfg(any(target_arch = "wasm32", test))]
mod sequence_similarity_model;
#[cfg(target_arch = "wasm32")]
mod sequence_similarity_view;

#[cfg(any(target_arch = "wasm32", test))]
mod structure_camera;
#[cfg(target_arch = "wasm32")]
mod structure_renderer;
#[cfg(target_arch = "wasm32")]
mod structure_view;
