//! Injected backend fixtures only. No user vault, D-Bus connection or OS prompt.
use collagene_desktop_credentials::*;
use collagene_library_protocol::LabId;
use collagene_repository::{RemoteConfig, http::CredentialProvider};
use std::{
    collections::BTreeMap,
    sync::{Arc, Mutex},
};
use zeroize::Zeroizing;
#[derive(Default)]
struct State {
    items: BTreeMap<String, Zeroizing<Vec<u8>>>,
    calls: usize,
    failure: Option<Error>,
}
#[derive(Clone, Default)]
struct Vault(Arc<Mutex<State>>);
impl Backend for Vault {
    fn read(&mut self, scope: &Scope) -> Result<Secret> {
        let mut state = self.0.lock().unwrap();
        state.calls += 1;
        if let Some(error) = state.failure {
            return Err(error);
        }
        let bytes = state.items.get(scope.account()).ok_or(Error::Missing)?;
        Secret::new(bytes.to_vec())
    }
    fn store(&mut self, scope: &Scope, secret: &Secret) -> Result<()> {
        let mut state = self.0.lock().unwrap();
        state.calls += 1;
        if let Some(error) = state.failure {
            return Err(error);
        }
        state.items.insert(
            scope.account().into(),
            Zeroizing::new(secret.as_bytes().to_vec()),
        );
        Ok(())
    }
    fn delete(&mut self, scope: &Scope) -> Result<()> {
        let mut state = self.0.lock().unwrap();
        state.calls += 1;
        if let Some(error) = state.failure {
            return Err(error);
        }
        state.items.remove(scope.account());
        Ok(())
    }
}
fn local(n: char) -> LabId {
    LabId::new(format!("lab-{}", n.to_string().repeat(32))).unwrap()
}
fn remote() -> RemoteConfig {
    RemoteConfig {
        id: format!("rem-{}", "1".repeat(32)),
        name: "Lab display name".into(),
        endpoint: "https://lab.example".into(),
        lab: local('2'),
    }
}
fn secret(text: &str) -> Secret {
    Secret::new(text.as_bytes().to_vec()).unwrap()
}
#[test]
fn every_security_scope_pin_changes_account_but_cosmetic_name_and_canonical_origin_do_not() {
    let a = local('a');
    let original = remote();
    let key = Scope::new(&a, &original).unwrap().account().to_owned();
    let mut renamed = original.clone();
    renamed.name = "Another display name".into();
    renamed.endpoint = "https://LAB.EXAMPLE:443/".into();
    assert_eq!(Scope::new(&a, &renamed).unwrap().account(), key);
    let mut alternatives = Vec::new();
    let mut changed = original.clone();
    changed.id = format!("rem-{}", "3".repeat(32));
    alternatives.push(changed);
    let mut changed = original.clone();
    changed.endpoint = "https://other.example/".into();
    alternatives.push(changed);
    let mut changed = original.clone();
    changed.endpoint = "https://lab.example:444/".into();
    alternatives.push(changed);
    let mut changed = original.clone();
    changed.endpoint = "http://lab.example/".into();
    alternatives.push(changed);
    let mut changed = original.clone();
    changed.lab = local('4');
    alternatives.push(changed);
    for changed in alternatives {
        assert_ne!(Scope::new(&a, &changed).unwrap().account(), key);
    }
    assert_ne!(Scope::new(&local('b'), &original).unwrap().account(), key);
    assert!(!key.contains("lab.example"));
    assert!(!key.contains("Lab display name"));
}
#[test]
fn session_mode_never_calls_a_locked_or_missing_os_service() {
    let vault = Vault::default();
    vault.0.lock().unwrap().failure = Some(Error::Locked);
    let mut credentials = Credentials::with_backend(local('a'), vault.clone());
    let remote = remote();
    assert_eq!(vault.0.lock().unwrap().calls, 0);
    credentials
        .store(&remote, secret("session-token"), Persistence::Session)
        .unwrap();
    assert_eq!(
        credentials.status(&remote).unwrap(),
        CredentialStatus::Session
    );
    assert!(credentials.bearer_for(&remote).is_ok());
    assert_eq!(vault.0.lock().unwrap().calls, 0);
    let mut other = remote.clone();
    other.lab = local('b');
    assert!(matches!(credentials.bearer(&other), Err(Error::Locked)));
    credentials.clear_session();
    assert!(matches!(credentials.bearer(&remote), Err(Error::Locked)));
}
#[test]
fn remembered_credentials_survive_provider_recreation_and_delete_only_exact_scope() {
    let vault = Vault::default();
    let a = local('a');
    let first = remote();
    let mut second = first.clone();
    second.lab = local('3');
    let mut credentials = Credentials::with_backend(a.clone(), vault.clone());
    credentials
        .store(&first, secret("first-token"), Persistence::Remember)
        .unwrap();
    credentials
        .store(&second, secret("second-token"), Persistence::Remember)
        .unwrap();
    drop(credentials);
    let mut reopened = Credentials::with_backend(a.clone(), vault.clone());
    assert_eq!(
        reopened.status(&first).unwrap(),
        CredentialStatus::Remembered
    );
    assert!(reopened.bearer_for(&first).is_ok());
    let mut other_repository = Credentials::with_backend(local('b'), vault.clone());
    assert!(matches!(
        other_repository.bearer(&first),
        Err(Error::Missing)
    ));
    reopened
        .store(&first, secret("temporary-override"), Persistence::Session)
        .unwrap();
    reopened.delete(&first).unwrap();
    assert_eq!(reopened.status(&first).unwrap(), CredentialStatus::Missing);
    assert_eq!(
        reopened.status(&second).unwrap(),
        CredentialStatus::Remembered
    );
    let state = vault.0.lock().unwrap();
    assert!(
        !state
            .items
            .contains_key(Scope::new(&a, &first).unwrap().account())
    );
    assert_eq!(
        state
            .items
            .get(Scope::new(&a, &second).unwrap().account())
            .unwrap()
            .as_slice(),
        b"second-token"
    );
}
#[test]
fn failed_remember_or_delete_has_explicit_error_without_hidden_session_fallback() {
    let vault = Vault::default();
    let mut credentials = Credentials::with_backend(local('a'), vault.clone());
    let remote = remote();
    credentials
        .store(&remote, secret("old-session"), Persistence::Session)
        .unwrap();
    vault.0.lock().unwrap().failure = Some(Error::Locked);
    assert_eq!(
        credentials.store(
            &remote,
            secret("new-secret-must-not-fall-back"),
            Persistence::Remember
        ),
        Err(Error::Locked)
    );
    assert_eq!(credentials.status(&remote), Err(Error::Locked));
    assert!(matches!(credentials.bearer(&remote), Err(Error::Locked)));
    assert!(vault.0.lock().unwrap().items.is_empty());
    credentials
        .store(
            &remote,
            secret("explicit-new-session"),
            Persistence::Session,
        )
        .unwrap();
    assert_eq!(credentials.delete(&remote), Err(Error::Locked));
    assert!(matches!(credentials.bearer(&remote), Err(Error::Locked)));
    let error = credentials.bearer_for(&remote).err().unwrap();
    assert_eq!(error.to_string(), "credential provider unavailable");
    assert!(!format!("{:?}", Error::Locked).contains("secret"));
}
#[test]
fn invalid_scope_and_token_material_are_rejected_without_entering_backend() {
    let vault = Vault::default();
    let mut credentials = Credentials::with_backend(local('a'), vault.clone());
    for endpoint in [
        "https://user:password@lab.example/",
        "https://lab.example/other",
        "https://lab.example/?key=secret",
        "https://lab.example/#secret",
        "file:///tmp/credentials",
    ] {
        let mut remote = remote();
        remote.endpoint = endpoint.into();
        assert_eq!(
            credentials.store(
                &remote,
                secret("valid-but-never-stored"),
                Persistence::Remember
            ),
            Err(Error::InvalidScope)
        );
    }
    for invalid in [
        vec![],
        b"contains space".to_vec(),
        b"header\r\ninjection".to_vec(),
        vec![0],
        vec![b'x'; MAX_TOKEN_BYTES + 1],
        "nonascii-é".as_bytes().to_vec(),
    ] {
        assert!(matches!(Secret::new(invalid), Err(Error::InvalidSecret)));
    }
    assert!(Secret::new(vec![b'x'; MAX_TOKEN_BYTES]).is_ok());
    assert_eq!(vault.0.lock().unwrap().calls, 0);
}
