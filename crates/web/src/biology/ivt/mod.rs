//! Source-bound IVT template inspection controls.
#[cfg(any(target_arch = "wasm32", test))]
mod linearization_model;
#[cfg(target_arch = "wasm32")]
mod linearization_view;
pub(crate) mod model;
#[cfg(target_arch = "wasm32")]
mod view;
#[cfg(target_arch = "wasm32")]
pub(crate) use view::IvtPanel;
#[cfg(test)]
mod tests;
