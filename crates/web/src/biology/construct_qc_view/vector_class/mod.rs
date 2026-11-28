//! Declared vector graph controls share the parent QC premises and references.
mod model;
#[cfg(target_arch = "wasm32")]
mod view;
#[cfg(target_arch = "wasm32")]
pub(crate) use view::VectorClassPanel;
#[cfg(test)]
mod tests;
