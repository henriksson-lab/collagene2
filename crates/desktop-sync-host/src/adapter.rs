use collagene_desktop_credentials as credentials;
use collagene_library_protocol::LabId;
use collagene_repository::{RemoteConfig, RemoteTransport, http::HttpRemoteTransport};
/// Native-only controls. The adapter owns authentication/transport material;
/// neither this trait nor its credential values are serialized as metadata.
pub trait SyncAdapter: RemoteTransport {
    fn local_repository_lab(&mut self) -> LabId;
    fn credential_status(
        &mut self,
        remote: &RemoteConfig,
    ) -> credentials::Result<credentials::CredentialStatus>;
    fn store_credential(
        &mut self,
        remote: &RemoteConfig,
        secret: credentials::Secret,
        persistence: credentials::Persistence,
    ) -> credentials::Result<()>;
    fn delete_credential(&mut self, remote: &RemoteConfig) -> credentials::Result<()>;
    fn clear_session(&mut self, remote: Option<&RemoteConfig>) -> credentials::Result<()>;
}
impl<B: credentials::Backend> SyncAdapter for HttpRemoteTransport<credentials::Credentials<B>> {
    fn local_repository_lab(&mut self) -> LabId {
        self.credentials_mut().local_repository_lab().clone()
    }
    fn credential_status(
        &mut self,
        remote: &RemoteConfig,
    ) -> credentials::Result<credentials::CredentialStatus> {
        self.credentials_mut().status(remote)
    }
    fn store_credential(
        &mut self,
        remote: &RemoteConfig,
        secret: credentials::Secret,
        persistence: credentials::Persistence,
    ) -> credentials::Result<()> {
        self.credentials_mut().store(remote, secret, persistence)
    }
    fn delete_credential(&mut self, remote: &RemoteConfig) -> credentials::Result<()> {
        self.credentials_mut().delete(remote)
    }
    fn clear_session(&mut self, remote: Option<&RemoteConfig>) -> credentials::Result<()> {
        match remote {
            Some(remote) => self.credentials_mut().forget_session(remote),
            None => {
                self.credentials_mut().clear_session();
                Ok(())
            }
        }
    }
}
pub type ProductionAdapter = HttpRemoteTransport<credentials::Credentials<credentials::OsVault>>;
