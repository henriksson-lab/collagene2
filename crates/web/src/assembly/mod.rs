//! assembly ownership boundary.
#[cfg(target_arch = "wasm32")]
pub(crate) mod controller;
pub(crate) mod draft;
pub(crate) mod edit;
pub(crate) mod library;
pub(crate) mod primer_design;
pub(crate) mod review;
