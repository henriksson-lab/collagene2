//! Windows generic credentials in the current user's local credential set.
use crate::{Error, MAX_TOKEN_BYTES, Result, SERVICE, Scope, Secret};
use windows_sys::Win32::{
    Foundation::{ERROR_ACCESS_DENIED, ERROR_NOT_FOUND, GetLastError},
    Security::Credentials::{
        CRED_PERSIST_LOCAL_MACHINE, CRED_TYPE_GENERIC, CREDENTIALW, CredDeleteW, CredFree,
        CredReadW, CredWriteW,
    },
};
use zeroize::Zeroize;
fn error() -> Error {
    match unsafe { GetLastError() } {
        ERROR_NOT_FOUND => Error::Missing,
        ERROR_ACCESS_DENIED => Error::Denied,
        _ => Error::Unavailable,
    }
}
fn target(scope: &Scope) -> Vec<u16> {
    format!("{SERVICE}:{}", scope.account())
        .encode_utf16()
        .chain(std::iter::once(0))
        .collect()
}
struct Allocation(*mut CREDENTIALW);
impl Drop for Allocation {
    fn drop(&mut self) {
        // SAFETY: only created from successful CredReadW, which owns this native
        // allocation and its blob. Clear its byte data before matching CredFree.
        unsafe {
            if !self.0.is_null() {
                let credential = &mut *self.0;
                if !credential.CredentialBlob.is_null() && credential.CredentialBlobSize > 0 {
                    std::slice::from_raw_parts_mut(
                        credential.CredentialBlob,
                        credential.CredentialBlobSize as usize,
                    )
                    .zeroize();
                }
                CredFree(self.0.cast());
            }
        }
    }
}
pub fn read(scope: &Scope) -> Result<Secret> {
    let name = target(scope);
    let mut pointer = std::ptr::null_mut();
    // SAFETY: terminated UTF-16 input and valid output pointer. Windows returns
    // one owned credential; its lifetime is bounded by Allocation below.
    if unsafe { CredReadW(name.as_ptr(), CRED_TYPE_GENERIC, 0, &mut pointer) } == 0 {
        return Err(error());
    }
    let allocation = Allocation(pointer);
    if allocation.0.is_null() {
        return Err(Error::Unavailable);
    }
    let credential = unsafe { &*allocation.0 };
    let length = credential.CredentialBlobSize as usize;
    if credential.Type != CRED_TYPE_GENERIC
        || length == 0
        || length > MAX_TOKEN_BYTES
        || credential.CredentialBlob.is_null()
    {
        return Err(Error::InvalidSecret);
    }
    // SAFETY: validated size and nonnull blob returned by successful CredRead.
    Secret::new(unsafe { std::slice::from_raw_parts(credential.CredentialBlob, length) }.to_vec())
}
pub fn store(scope: &Scope, secret: &Secret) -> Result<()> {
    let mut name = target(scope);
    let mut account = scope
        .account()
        .encode_utf16()
        .chain(std::iter::once(0))
        .collect::<Vec<_>>();
    let credential = CREDENTIALW {
        Type: CRED_TYPE_GENERIC,
        TargetName: name.as_mut_ptr(),
        CredentialBlobSize: secret.as_bytes().len() as u32,
        CredentialBlob: secret.as_bytes().as_ptr().cast_mut(),
        Persist: CRED_PERSIST_LOCAL_MACHINE,
        UserName: account.as_mut_ptr(),
        ..Default::default()
    };
    // SAFETY: all buffers outlive this synchronous call; CredWrite copies the
    // input blob and does not retain or modify our borrowed secret data.
    if unsafe { CredWriteW(&credential, 0) } == 0 {
        Err(error())
    } else {
        Ok(())
    }
}
pub fn delete(scope: &Scope) -> Result<()> {
    let name = target(scope);
    // SAFETY: valid terminated exact target; generic credential type and flags.
    if unsafe { CredDeleteW(name.as_ptr(), CRED_TYPE_GENERIC, 0) } != 0 {
        return Ok(());
    }
    match error() {
        Error::Missing => Ok(()),
        other => Err(other),
    }
}
