//! Explicit named-product screening analysis controls; no product/history writes.
#[cfg(target_arch = "wasm32")]
mod editors;
mod model;
#[cfg(target_arch = "wasm32")]
mod rendering;
#[cfg(target_arch = "wasm32")]
mod view;
#[cfg(target_arch = "wasm32")]
pub(crate) use view::ScreeningPanel;
#[cfg(test)]
mod tests;
