//! Generic Keychain items, device-local/non-synchronizable, with UI forbidden.
use crate::{Error, MAX_TOKEN_BYTES, Result, SERVICE, Scope, Secret};
use core_foundation::{
    base::{CFType, CFTypeRef, TCFType},
    boolean::CFBoolean,
    data::CFData,
    dictionary::CFDictionary,
    string::{CFString, CFStringRef},
};
use security_framework_sys::{
    base::{errSecAuthFailed, errSecDuplicateItem, errSecItemNotFound},
    item::*,
    keychain_item::{SecItemAdd, SecItemCopyMatching, SecItemDelete, SecItemUpdate},
};
#[link(name = "Security", kind = "framework")]
unsafe extern "C" {
    static kSecUseAuthenticationUIFail: CFStringRef;
}
fn status(code: i32) -> Error {
    // Convert only the numeric status; never format Apple's possibly revealing
    // error descriptions into logs or public IPC failures.
    let code = security_framework::base::Error::from_code(code).code();
    match code {
        errSecItemNotFound => Error::Missing,
        -25308 | -128 => Error::Locked,
        errSecAuthFailed => Error::Denied,
        _ => Error::Unavailable,
    }
}
unsafe fn key(value: CFStringRef) -> CFString {
    unsafe { CFString::wrap_under_get_rule(value) }
}
fn pairs(scope: &Scope) -> Vec<(CFString, CFType)> {
    // SAFETY: Security exports process-lifetime CFString constants. Wrapping
    // under get-rule retains each value before its safe dictionary takes it.
    unsafe {
        vec![
            (key(kSecClass), key(kSecClassGenericPassword).into_CFType()),
            (key(kSecAttrService), CFString::new(SERVICE).into_CFType()),
            (
                key(kSecAttrAccount),
                CFString::new(scope.account()).into_CFType(),
            ),
            (
                key(kSecAttrSynchronizable),
                CFBoolean::false_value().into_CFType(),
            ),
            (
                key(kSecUseAuthenticationUI),
                key(kSecUseAuthenticationUIFail).into_CFType(),
            ),
        ]
    }
}
pub fn read(scope: &Scope) -> Result<Secret> {
    let mut pairs = pairs(scope);
    unsafe {
        pairs.push((key(kSecReturnData), CFBoolean::true_value().into_CFType()));
    }
    let query = CFDictionary::from_CFType_pairs(&pairs);
    let mut output: CFTypeRef = std::ptr::null();
    // SAFETY: valid owned CF dictionary; Security transfers an owned result on
    // success. We check type and bounds before copying its byte slice.
    let code = unsafe { SecItemCopyMatching(query.as_concrete_TypeRef(), &mut output) };
    if code != 0 {
        return Err(status(code));
    }
    if output.is_null() {
        return Err(Error::Unavailable);
    }
    let result = unsafe { CFType::wrap_under_create_rule(output) };
    let data = result.downcast::<CFData>().ok_or(Error::InvalidSecret)?;
    if data.bytes().len() > MAX_TOKEN_BYTES {
        return Err(Error::InvalidSecret);
    }
    Secret::new(data.bytes().to_vec())
}
pub fn store(scope: &Scope, secret: &Secret) -> Result<()> {
    let pairs = pairs(scope);
    let query = CFDictionary::from_CFType_pairs(&pairs);
    let data = CFData::from_buffer(secret.as_bytes());
    let updates =
        unsafe { CFDictionary::from_CFType_pairs(&[(key(kSecValueData), data.into_CFType())]) };
    // SAFETY: both dictionaries own valid CF objects for the full call; the
    // system copies data into its Keychain and no pointers escape this function.
    let mut code =
        unsafe { SecItemUpdate(query.as_concrete_TypeRef(), updates.as_concrete_TypeRef()) };
    if code == errSecItemNotFound {
        let mut creating = pairs;
        unsafe {
            creating.push((
                key(kSecValueData),
                CFData::from_buffer(secret.as_bytes()).into_CFType(),
            ));
        }
        let creating = CFDictionary::from_CFType_pairs(&creating);
        code = unsafe { SecItemAdd(creating.as_concrete_TypeRef(), std::ptr::null_mut()) };
        if code == errSecDuplicateItem {
            code = unsafe {
                SecItemUpdate(query.as_concrete_TypeRef(), updates.as_concrete_TypeRef())
            };
        }
    }
    if code == 0 { Ok(()) } else { Err(status(code)) }
}
pub fn delete(scope: &Scope) -> Result<()> {
    let query = CFDictionary::from_CFType_pairs(&pairs(scope));
    // SAFETY: valid owned query and exact service/account, no wildcard search.
    let code = unsafe { SecItemDelete(query.as_concrete_TypeRef()) };
    if code == 0 || code == errSecItemNotFound {
        Ok(())
    } else {
        Err(status(code))
    }
}
