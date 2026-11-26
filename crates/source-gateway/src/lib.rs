//! Independent MIT source gateway. Scientific parsing and network policy remain
//! in `collagene-sources`; this crate grants no filesystem or executable access.
#[cfg(not(target_arch = "wasm32"))]
mod http;
#[cfg(not(target_arch = "wasm32"))]
mod service;
pub mod wire;
#[cfg(not(target_arch = "wasm32"))]
pub use http::{CONTENT_TYPE, PATH, serve};
#[cfg(not(target_arch = "wasm32"))]
pub use service::{Backend, Cancellation, Config, ConfigError, Gateway, Limits, NativeBackend};
