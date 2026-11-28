//! Fresh verification over current typed plans, with explicit evidence and output owners.
pub(crate) mod ancestry;
pub mod checks;
pub mod engine;
pub mod evidence;
pub mod output;
pub mod pcr;
pub mod projection;
pub(crate) mod replay_budget;
pub mod sources;
pub mod terminal;
pub const ENGINE_VERSION: &str = concat!(env!("CARGO_PKG_VERSION"), "/science-38");
