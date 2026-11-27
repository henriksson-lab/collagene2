//! Human discussion dispatch accepts only a trusted individual member capability.
use super::*;
use store::collaboration as native;
use wire::discussion as api;

pub fn dispatch(
    db: &mut store::Library,
    member: &native::MemberAccess,
    request: &api::RequestEnvelope,
) -> Result<api::ResponseEnvelope, Error> {
    request.validate().map_err(|_| Error::Wire)?;
    let access = db
        .trusted_host_access(member.lab_id())
        .map_err(|_| Error::Unavailable)?;
    let lab = db.lab(&access).map_err(|_| Error::Unavailable)?;
    let target = |t: &api::Target| native::Target {
        document: document(&t.document),
        revision: revision(&t.revision),
    };
    let result: Result<api::Success, native::DiscussionError> = match &request.request {
        api::Request::Identify => Ok(api::Success::Identified),
        api::Request::List {
            target: t,
            offset,
            limit,
            expected_generation,
        } => db
            .discussion_page(member, &target(t), *offset, *limit, *expected_generation)
            .and_then(|p| {
                convert(&p)
                    .map(api::Success::Listed)
                    .map_err(|_| store::Error::Corrupt("discussion conversion").into())
            }),
        api::Request::Append {
            target: t,
            expected_generation,
            action,
        } => {
            let action: native::Action = convert(action)?;
            db.append_discussion(member, &target(t), *expected_generation, action)
                .and_then(|e| {
                    convert(&e)
                        .map(|event| api::Success::Appended {
                            target: t.clone(),
                            event,
                        })
                        .map_err(|_| store::Error::Corrupt("discussion conversion").into())
                })
        }
    };
    let outcome = match result {
        Ok(s) => api::Outcome::Success(s),
        Err(e) => {
            use native::DiscussionError as D;
            let (code, actual_generation) = match e {
                D::Forbidden => (wire::ErrorCode::Forbidden, None),
                D::Conflict { actual, .. } => (wire::ErrorCode::Conflict, Some(actual)),
                D::Storage(store::Error::NotFound) => (wire::ErrorCode::NotFound, None),
                D::Storage(store::Error::Invalid(_)) => (wire::ErrorCode::Invalid, None),
                D::Storage(store::Error::Limit(_)) => (wire::ErrorCode::TooLarge, None),
                _ => (wire::ErrorCode::Unavailable, None),
            };
            api::Outcome::Failure(api::Failure {
                code,
                actual_generation,
            })
        }
    };
    let response = api::ResponseEnvelope {
        request_id: request.request_id,
        context: wire::LabContext {
            lab_id: wire::LabId::new(lab.id.as_str()).map_err(|_| Error::Unavailable)?,
            name: wire::Text::new(lab.name).map_err(|_| Error::Unavailable)?,
        },
        member: api::Member {
            author: convert(member.author())?,
            role: convert(&member.role())?,
        },
        outcome,
    };
    response
        .validate_for(request, &response.context.lab_id)
        .map_err(|_| Error::Unavailable)?;
    Ok(response)
}
