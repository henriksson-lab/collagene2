//! Session-scoped native MCP access to deterministic cloning tools.

pub mod server;
pub mod session;
pub mod tools;

/// Number of tools registered by the current native server.
pub const CURRENT_TOOL_COUNT: usize = 219;
