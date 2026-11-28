//! assembly / primer_design ownership boundary.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_edited_primer_commit;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_edited_primer_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_edited_primer_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_extension_primer_commit;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_extension_primer_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_extension_primer_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod assembly_inverse_pcr_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod assembly_inverse_pcr_view;
