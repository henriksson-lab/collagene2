//! Design/consensus mode of the existing exact-version comparison workflow.
#[cfg(any(target_arch = "wasm32", test))]
mod model;
#[cfg(target_arch = "wasm32")]
pub(super) mod rendering;
#[cfg(target_arch = "wasm32")]
mod view;

#[cfg(target_arch = "wasm32")]
pub(crate) use view::DesignConsensusPanel;
