//! Bounded stateless worker entrypoint, typed domain routing and current response serialization.
mod assembly;
mod context;
mod documents;
mod library;
mod material;
mod output;
mod preparation;
mod restriction;
mod review;
mod router;
mod sequence;
mod topology;
mod tracks;

use super::{
    envelope::{WORKSPACE_REQUEST_SELECTORS, decode_with_limit, worker_limit},
    request::Request,
    response::Response,
};
use collagene_model::Error;

pub fn process_request(bytes: &[u8]) -> Result<Vec<u8>, Error> {
    let maximum = worker_limit(b"CGQ\0", bytes, bytes.len(), WORKSPACE_REQUEST_SELECTORS)?;
    let Request { id, revision, job } = decode_with_limit(b"CGQ\0", bytes, maximum)?;
    let context = context::RequestContext {
        revision: &revision,
    };
    let result = router::dispatch(job, &context)?;
    Response {
        id,
        revision,
        result,
    }
    .encode()
}
