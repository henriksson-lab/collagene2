//! biology ownership boundary.
pub(crate) mod antibody;
#[cfg(target_arch = "wasm32")]
pub(crate) mod biology;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod biology_selection;
#[cfg(target_arch = "wasm32")]
pub(crate) mod biology_source_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod codon_design;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod codon_preferences;
#[cfg(target_arch = "wasm32")]
pub(crate) mod construct_qc_elements;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod construct_qc_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod construct_qc_scans;
#[cfg(target_arch = "wasm32")]
pub(crate) mod construct_qc_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod controller;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod ivt;
#[cfg(target_arch = "wasm32")]
pub(crate) mod methylation;
#[cfg(target_arch = "wasm32")]
pub(crate) mod protein_export;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod protein_liability_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod protein_liability_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod protein_properties;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod rna_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod rna_view;
pub(crate) mod structure;
#[cfg(all(test, not(target_arch = "wasm32")))]
#[path = "construct_qc_view/vector_class/mod.rs"]
mod vector_class_tests;

#[cfg(all(test, not(target_arch = "wasm32")))]
#[path = "construct_qc_view/causal/mod.rs"]
mod causal_qc_tests;
