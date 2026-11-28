//! restriction ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod controller;
#[cfg(target_arch = "wasm32")]
pub(crate) mod digest;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod digest_fragment_model;
#[cfg(target_arch = "wasm32")]
pub(crate) mod digest_product;
#[cfg(target_arch = "wasm32")]
pub(crate) mod gel;
#[cfg(target_arch = "wasm32")]
pub(crate) mod gel_export;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod material_derivation_model;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod restriction_label_groups;
#[cfg(target_arch = "wasm32")]
pub(crate) mod restriction_map;
#[cfg(target_arch = "wasm32")]
pub(crate) mod restriction_view;
#[cfg(target_arch = "wasm32")]
pub(crate) mod reusable_material_view;
