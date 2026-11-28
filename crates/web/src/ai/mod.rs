//! ai ownership boundary.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod web_mcp;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod web_mcp_assembly;
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) mod web_mcp_design;
