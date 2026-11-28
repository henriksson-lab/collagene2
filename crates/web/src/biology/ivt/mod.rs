//! Source-bound IVT template inspection controls.
pub(crate) mod model;
#[cfg(target_arch = "wasm32")]
mod view;
#[cfg(target_arch = "wasm32")]
pub(crate) use view::IvtPanel;
#[cfg(test)]
mod tests;
