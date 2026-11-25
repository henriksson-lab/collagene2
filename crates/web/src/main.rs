#[cfg(target_arch = "wasm32")]
mod alignment_view;
#[cfg(any(target_arch = "wasm32", test))]
mod annotation_geometry;
#[cfg(any(target_arch = "wasm32", test))]
mod annotation_layout;
#[cfg(target_arch = "wasm32")]
mod app;
#[cfg(target_arch = "wasm32")]
mod assembly_finalize;
#[cfg(target_arch = "wasm32")]
mod assembly_outputs;
#[cfg(target_arch = "wasm32")]
mod assembly_product_view;
#[cfg(any(target_arch = "wasm32", test))]
mod assembly_review;
#[cfg(target_arch = "wasm32")]
mod assembly_workspace;
#[cfg(target_arch = "wasm32")]
mod checked_cloning;
#[cfg(target_arch = "wasm32")]
mod clipboard;
#[cfg(target_arch = "wasm32")]
mod cloning_methods;
#[cfg(target_arch = "wasm32")]
mod codon_design;
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
#[cfg(target_arch = "wasm32")]
mod dotplot;
#[cfg(target_arch = "wasm32")]
mod end_repair;
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
mod history_view;
#[cfg(target_arch = "wasm32")]
mod map_labels;
#[cfg(target_arch = "wasm32")]
mod oligo_import;
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
mod restriction_map;
#[cfg(target_arch = "wasm32")]
mod restriction_view;
#[cfg(target_arch = "wasm32")]
mod run_capture;
#[cfg(target_arch = "wasm32")]
mod search_view;
#[cfg(target_arch = "wasm32")]
mod selection_toolbar;
#[cfg(target_arch = "wasm32")]
mod sequence_view;
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
mod workspace_layout;
#[cfg(any(target_arch = "wasm32", test))]
mod workspace_preferences;
#[cfg(target_arch = "wasm32")]
mod workspace_view;

#[cfg(target_arch = "wasm32")]
mod biology;
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
mod mutagenesis;
#[cfg(any(target_arch = "wasm32", test))]
mod mutation_review;
#[cfg(target_arch = "wasm32")]
mod primer_pair_editor;
#[cfg(any(target_arch = "wasm32", test))]
mod primer_pair_model;
#[cfg(any(target_arch = "wasm32", test))]
mod primer_pair_storage;
#[cfg(target_arch = "wasm32")]
mod primer_picker;
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
