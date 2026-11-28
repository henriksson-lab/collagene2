//! Same discussion UI, with one explicit authenticated transport per deployment.
use crate::library::library_client;
use crate::library::library_client::LibraryClient;
use crate::platform::desktop_sync_client;
use collagene_desktop_sync_protocol as native;
use collagene_library_protocol::{LabId, discussion as api};
use std::rc::Rc;
#[derive(Clone)]
pub enum Client {
    Browser(Rc<library_client::Session>),
    Remote {
        client: desktop_sync_client::Client,
        remote: native::Remote,
    },
}
impl Client {
    pub fn browser() -> Self {
        Self::Browser(library_client::browser_session())
    }
    pub fn lab_id(&self) -> Option<LabId> {
        match self {
            Self::Browser(s) => s.context().map(|c| c.lab_id),
            Self::Remote { remote, .. } => Some(remote.expected_lab.clone()),
        }
    }
    pub fn publication_epoch(&self) -> u64 {
        match self {
            Self::Browser(s) => s.publication_epoch(),
            Self::Remote { .. } => 0,
        }
    }
    pub async fn request(&self, request: api::Request) -> Result<api::ResponseEnvelope, String> {
        match self {
            Self::Browser(s) => s.discussion(request).await,
            Self::Remote { client, remote } => {
                match client
                    .request(
                        native::Command::Discussion {
                            remote: remote.id.clone(),
                            request,
                        },
                        Some(remote.clone()),
                    )
                    .await?
                {
                    native::Success::Discussion(reply) => Ok(reply),
                    _ => Err("Native discussion reply has the wrong operation type.".into()),
                }
            }
        }
    }
}
