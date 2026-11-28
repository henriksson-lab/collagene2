//! reads ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod contigs;
#[cfg(target_arch = "wasm32")]
pub(crate) mod controller;
#[cfg(target_arch = "wasm32")]
pub(crate) mod read_alignments_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod read_checks_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod read_review_capture;
#[cfg(target_arch = "wasm32")]
pub(crate) mod read_review_history;
#[cfg(target_arch = "wasm32")]
pub(crate) mod read_review_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod read_source_library;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod scf_export_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod scf_export_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod trace_validation;
#[cfg(target_arch = "wasm32")]
pub(crate) mod trace_view;
