//! primers ownership boundary.
pub(crate) mod pcr;
#[cfg(target_arch = "wasm32")]
pub(crate) mod primer_binding_view;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod primer_domain_edit;
#[cfg(target_arch = "wasm32")]
pub(crate) mod primer_domain_edit_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod primer_orders;
#[cfg(target_arch = "wasm32")]
pub(crate) mod primer_pair_editor;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod primer_pair_model;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod primer_pair_storage;
#[cfg(target_arch = "wasm32")]
pub(crate) mod primer_picker;
#[cfg(target_arch = "wasm32")]
pub(crate) mod primer_reference_capture_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod primers;

#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod collection_binding_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod collection_binding_view;
