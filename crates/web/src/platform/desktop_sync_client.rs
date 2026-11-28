//! Native-only DSQ boundary. Root injects the existing collagene_rpc lane;
//! this module never invokes Tauri, HTTP or credential storage itself.
use collagene_desktop_sync_protocol as wire;
use std::{future::Future, pin::Pin, rc::Rc};
pub type LocalFuture<T> = Pin<Box<dyn Future<Output = Result<T, String>>>>;
pub type RequestFactory = Rc<dyn Fn(wire::Command) -> Result<wire::RequestEnvelope, String>>;
pub type RawRpc = Rc<dyn Fn(Vec<u8>) -> LocalFuture<Vec<u8>>>;
#[derive(Clone)]
pub struct Client {
    factory: RequestFactory,
    rpc: RawRpc,
}
impl Client {
    /// Factory must allocate unique per-window operation IDs in the SAME native
    /// lane as DTQ; it must not clone/replay requests containing credentials.
    pub fn new(factory: RequestFactory, rpc: RawRpc) -> Self {
        Self { factory, rpc }
    }
    pub async fn request(
        &self,
        command: wire::Command,
        expected: Option<wire::Remote>,
    ) -> Result<wire::Success, String> {
        match self.request_outcome(command, expected).await? {
            wire::Outcome::Success(value) => Ok(value),
            wire::Outcome::Failure(failure) => Err(failure_message(&failure)),
        }
    }
    /// Preserve typed conflicts for an explicit, pinned review. Transport and
    /// binding failures remain errors and must never enable a checkout action.
    pub async fn request_outcome(
        &self,
        command: wire::Command,
        expected: Option<wire::Remote>,
    ) -> Result<wire::Outcome, String> {
        let request = (self.factory)(command)
            .map_err(|_| "Cannot allocate native database operation.".to_owned())?;
        let bytes = wire::encode_request(&request)
            .map_err(|_| "Invalid native database request.".to_owned())?;
        let discussion = matches!(&request.command, wire::Command::Discussion { .. });
        let bytes = (self.rpc)(bytes).await.map_err(|_| {
            if discussion {
                "Remote discussion interrupted. Reload its history before choosing to resubmit."
                    .to_owned()
            } else {
                "Database operation interrupted. Review its durable job before retrying.".to_owned()
            }
        })?;
        let response = wire::decode_response_for(&bytes, &request)
            .map_err(|_| "Native database response failed its request checks.".to_owned())?;
        // Some local-only acknowledgements deliberately have no remote context;
        // whenever native execution returns one, bind the complete known scope.
        if let Some(expected) = expected {
            if response.context.is_some() {
                response
                    .validate_in_remote(&request, &expected)
                    .map_err(|_| {
                        "Database response came from a different configured remote.".to_owned()
                    })?;
            }
        }
        Ok(response.outcome)
    }
}
pub fn failure_message(failure: &wire::Failure) -> String {
    match failure {
        wire::Failure::LocalConflict{..}=>"The saved local plasmid changed. Refresh the local library and review it before retrying.".into(),
        wire::Failure::Remote(f)=>format!("Remote {:?}: {}",f.code,f.message.as_str()),
        wire::Failure::Error{code,message}=>format!("{:?}: {}",code,message.as_str()),
    }
}
