use crate::{Backend, Result, Scope, Secret};
#[cfg(target_os = "linux")]
#[path = "linux.rs"]
mod native;
#[cfg(target_os = "macos")]
#[path = "macos.rs"]
mod native;
#[cfg(target_os = "windows")]
#[path = "windows.rs"]
mod native;

/// Lazy OS backend: constructing it never connects to a vault, unlocks it or
/// prompts. Only explicit credential operations reach platform services.
#[derive(Default)]
pub struct OsVault;
impl OsVault {
    pub fn new() -> Self {
        Self
    }
}
impl Backend for OsVault {
    fn read(&mut self, scope: &Scope) -> Result<Secret> {
        native::read(scope)
    }
    fn store(&mut self, scope: &Scope, secret: &Secret) -> Result<()> {
        native::store(scope, secret)
    }
    fn delete(&mut self, scope: &Scope) -> Result<()> {
        native::delete(scope)
    }
}
#[cfg(not(any(target_os = "linux", target_os = "macos", target_os = "windows")))]
mod native {
    use super::*;
    pub fn read(_: &Scope) -> Result<Secret> {
        Err(crate::Error::Unavailable)
    }
    pub fn store(_: &Scope, _: &Secret) -> Result<()> {
        Err(crate::Error::Unavailable)
    }
    pub fn delete(_: &Scope) -> Result<()> {
        Err(crate::Error::Unavailable)
    }
}
