#[cfg(target_arch = "wasm32")]
mod alignment_view;
#[cfg(target_arch = "wasm32")]
mod app;
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
mod history_view;
#[cfg(target_arch = "wasm32")]
mod oligo_import;
#[cfg(target_arch = "wasm32")]
mod protein_properties;
#[cfg(target_arch = "wasm32")]
mod run_capture;
#[cfg(target_arch = "wasm32")]
mod search_view;
#[cfg(target_arch = "wasm32")]
mod sequence_view;
#[cfg(target_arch = "wasm32")]
mod trace_validation;
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
#[cfg(target_arch = "wasm32")]
mod primers;
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
