//! Trusted native credential boundary. No credential is a field in a repository,
//! configuration, archive or synchronization message. The caller owns approval
//! to use/remember a token and supplies the actual local repository lab scope.
use collagene_library_protocol::LabId;
use collagene_repository::{
    RemoteConfig,
    http::{BearerToken, CredentialError, CredentialProvider},
};
use std::collections::BTreeMap;
use zeroize::Zeroizing;
mod platform;
pub use platform::OsVault;

pub const MAX_TOKEN_BYTES: usize = 2048; // portable below Windows generic blob bound
pub const SERVICE: &str = "org.henrikssonlab.collagene2.remote-access.v1";
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum Error {
    #[error("invalid trusted credential scope")]
    InvalidScope,
    #[error("invalid access token")]
    InvalidSecret,
    #[error("no credential is stored for this repository and remote")]
    Missing,
    #[error("the native credential vault is locked or requires interaction")]
    Locked,
    #[error("native credential access was denied")]
    Denied,
    #[error("the native credential vault is unavailable")]
    Unavailable,
    #[error("multiple native credentials match this exact scope")]
    Ambiguous,
}
pub type Result<T> = std::result::Result<T, Error>;

/// Deliberately no Clone, Debug, Display or serde implementations. Owned buffers
/// are zeroized on drop; OS/HTTP libraries still manage their own native copies.
pub struct Secret(Zeroizing<Vec<u8>>);
impl Secret {
    pub fn new(bytes: Vec<u8>) -> Result<Self> {
        let bytes = Zeroizing::new(bytes);
        if bytes.is_empty()
            || bytes.len() > MAX_TOKEN_BYTES
            || !bytes.iter().all(|b| (0x21..=0x7e).contains(b))
        {
            return Err(Error::InvalidSecret);
        }
        Ok(Self(bytes))
    }
    /// Explicit borrowed access for the native backend, never a log formatter.
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }
    fn bearer(&self) -> Result<BearerToken> {
        BearerToken::new(std::str::from_utf8(&self.0).map_err(|_| Error::InvalidSecret)?)
            .map_err(|_| Error::InvalidSecret)
    }
}
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Scope {
    account: String,
}
impl Scope {
    /// local_repository_lab is trusted native context, not a remote message field.
    /// Canonicalizing an HTTP origin grants no permission to transmit a token:
    /// HttpRemoteTransport separately requires HTTPS or explicit literal loopback.
    pub fn new(local_repository_lab: &LabId, remote: &RemoteConfig) -> Result<Self> {
        if remote.endpoint.is_empty()
            || remote.endpoint.len() > 2048
            || remote.endpoint.chars().any(char::is_control)
            || remote.endpoint.trim() != remote.endpoint
            || remote.id.len() != 36
            || !remote.id.starts_with("rem-")
            || !remote.id[4..]
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err(Error::InvalidScope);
        }
        let endpoint = url::Url::parse(&remote.endpoint).map_err(|_| Error::InvalidScope)?;
        if !matches!(endpoint.scheme(), "https" | "http")
            || endpoint.host_str().is_none()
            || !endpoint.username().is_empty()
            || endpoint.password().is_some()
            || endpoint.path() != "/"
            || endpoint.query().is_some()
            || endpoint.fragment().is_some()
        {
            return Err(Error::InvalidScope);
        }
        let origin = endpoint.origin().ascii_serialization();
        let mut hash = blake3::Hasher::new();
        hash.update(b"collagene-native-credential-scope-v1\0");
        for field in [
            local_repository_lab.as_str(),
            remote.id.as_str(),
            origin.as_str(),
            remote.lab.as_str(),
        ] {
            hash.update(&(field.len() as u64).to_le_bytes());
            hash.update(field.as_bytes());
        }
        Ok(Self {
            account: format!("scope-v1-{}", hash.finalize().to_hex()),
        })
    }
    /// Opaque nonsecret account identifier; no raw URL, member name or token.
    pub fn account(&self) -> &str {
        &self.account
    }
}
/// Injectable native backend. This interface is not serialized or exposed to IPC.
/// Backends must use the fixed service and exact scope without key enumeration.
pub trait Backend {
    fn read(&mut self, scope: &Scope) -> Result<Secret>;
    fn store(&mut self, scope: &Scope, secret: &Secret) -> Result<()>;
    fn delete(&mut self, scope: &Scope) -> Result<()>;
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CredentialStatus {
    Session,
    Remembered,
    Missing,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Persistence {
    Session,
    Remember,
}
pub struct Credentials<B = OsVault> {
    local_lab: LabId,
    backend: B,
    session: BTreeMap<Scope, Secret>,
}
impl Credentials<OsVault> {
    pub fn new(trusted_local_repository_lab: LabId) -> Self {
        Self::with_backend(trusted_local_repository_lab, OsVault::new())
    }
}
impl<B: Backend> Credentials<B> {
    pub fn with_backend(trusted_local_repository_lab: LabId, backend: B) -> Self {
        Self {
            local_lab: trusted_local_repository_lab,
            backend,
            session: BTreeMap::new(),
        }
    }
    /// Nonsecret trusted local namespace, for native service scope checks.
    pub fn local_repository_lab(&self) -> &LabId {
        &self.local_lab
    }
    /// Consumes the submitted secret. Session mode never calls the OS. Remember
    /// errors are explicit and never install a plaintext/session fallback.
    pub fn store(
        &mut self,
        remote: &RemoteConfig,
        secret: Secret,
        persistence: Persistence,
    ) -> Result<()> {
        let scope = Scope::new(&self.local_lab, remote)?;
        match persistence {
            Persistence::Session => {
                self.session.insert(scope, secret);
                Ok(())
            }
            Persistence::Remember => {
                self.session.remove(&scope);
                self.backend.store(&scope, &secret)
            }
        }
    }
    /// Removes active memory before attempting exact scoped OS deletion. A vault
    /// error means persisted deletion was not confirmed; never claim otherwise.
    pub fn delete(&mut self, remote: &RemoteConfig) -> Result<()> {
        let scope = Scope::new(&self.local_lab, remote)?;
        self.session.remove(&scope);
        match self.backend.delete(&scope) {
            Err(Error::Missing) => Ok(()),
            result => result,
        }
    }
    pub fn status(&mut self, remote: &RemoteConfig) -> Result<CredentialStatus> {
        let scope = Scope::new(&self.local_lab, remote)?;
        if self.session.contains_key(&scope) {
            return Ok(CredentialStatus::Session);
        }
        match self.backend.read(&scope) {
            Ok(_) => Ok(CredentialStatus::Remembered),
            Err(Error::Missing) => Ok(CredentialStatus::Missing),
            Err(error) => Err(error),
        }
    }
    pub fn clear_session(&mut self) {
        self.session.clear();
    }
    pub fn forget_session(&mut self, remote: &RemoteConfig) -> Result<()> {
        let scope = Scope::new(&self.local_lab, remote)?;
        self.session.remove(&scope);
        Ok(())
    }
    pub fn bearer(&mut self, remote: &RemoteConfig) -> Result<BearerToken> {
        let scope = Scope::new(&self.local_lab, remote)?;
        if let Some(secret) = self.session.get(&scope) {
            return secret.bearer();
        }
        self.backend.read(&scope)?.bearer()
    }
}
impl<B: Backend> CredentialProvider for Credentials<B> {
    fn bearer_for(
        &mut self,
        remote: &RemoteConfig,
    ) -> std::result::Result<BearerToken, CredentialError> {
        self.bearer(remote).map_err(|_| CredentialError)
    }
}
