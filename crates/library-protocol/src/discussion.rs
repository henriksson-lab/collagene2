//! Bounded human collaboration messages; no client-supplied author or role.
use super::*;
pub const ENDPOINT: &str = "/api/library/discussion";
pub const MEDIA_TYPE: &str = "application/vnd.collagene.discussion+postcard";
pub const MAX_REQUEST_BYTES: usize = 16 * 1024;
pub const MAX_RESPONSE_BYTES: usize = 2 * 1024 * 1024;
pub const MAX_EVENTS: u64 = 4096;
pub const MAX_BODY_BYTES: usize = 8192;
const REQUEST_MAGIC: &[u8; 4] = b"CDQ\0";
const RESPONSE_MAGIC: &[u8; 4] = b"CDR\0";
dto! {pub struct Target {pub document:DocumentId,pub revision:RevisionId}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Role {
    Reader,
    Contributor,
    Maintainer,
}
dto! {pub struct Author {pub id:Text<128>,pub display_name:Text<256>}}
impl Author {
    pub fn validate(&self) -> Result<(), Error> {
        if self.id.is_empty()
            || !self
                .id
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"-_.:@".contains(&b))
            || self.display_name.trim().is_empty()
            || self.display_name.chars().any(char::is_control)
        {
            return Err(Error::Invalid("discussion author"));
        }
        Ok(())
    }
}
dto! {pub struct Member {pub author:Author,pub role:Role}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReviewDecision {
    Approved,
    ChangesRequested,
    Withdrawn,
}
dto! {pub enum Action {
    OpenThread{body:Text<MAX_BODY_BYTES>},Reply{thread:u64,body:Text<MAX_BODY_BYTES>},
    SetResolved{thread:u64,resolved:bool},Review{decision:ReviewDecision,body:Text<MAX_BODY_BYTES>}
}}
impl Action {
    pub fn validate(&self) -> Result<(), Error> {
        let (body, required) = match self {
            Self::OpenThread { body } | Self::Reply { body, .. } => (Some(body), true),
            Self::Review { body, .. } => (Some(body), false),
            Self::SetResolved { .. } => (None, false),
        };
        if body.is_some_and(|b| {
            (required && b.trim().is_empty())
                || b.chars().any(|c| c.is_control() && c != '\n' && c != '\t')
        }) {
            return Err(Error::Invalid("discussion body"));
        }
        if matches!(self,Self::Reply{thread,..}|Self::SetResolved{thread,..} if *thread==0||*thread>MAX_EVENTS)
        {
            return Err(Error::Invalid("discussion thread"));
        }
        Ok(())
    }
}
dto! {pub struct Event {pub sequence:u64,pub author:Author,pub role:Role,pub created_at:u64,pub action:Action}}
impl Event {
    fn validate(&self) -> Result<(), Error> {
        self.author.validate()?;
        self.action.validate()?;
        if self.sequence == 0 || self.sequence > MAX_EVENTS || self.role == Role::Reader {
            return Err(Error::Invalid("discussion event"));
        }
        if matches!(&self.action,Action::Reply{thread,..}|Action::SetResolved{thread,..} if *thread>=self.sequence)
        {
            return Err(Error::Invalid("discussion thread order"));
        }
        Ok(())
    }
}
dto! {pub enum Request {
    Identify,
    List{target:Target,offset:u32,limit:u16,expected_generation:Option<u64>},
    Append{target:Target,expected_generation:u64,action:Action}
}}
impl Request {
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Identify => Ok(()),
            Self::List {
                offset,
                limit,
                expected_generation,
                ..
            } => {
                if u64::from(*offset) > MAX_EVENTS
                    || *limit == 0
                    || usize::from(*limit) > MAX_PAGE_ROWS
                    || (*offset > 0 && expected_generation.is_none())
                    || expected_generation.is_some_and(|v| v > MAX_EVENTS)
                {
                    return Err(Error::Invalid("discussion page"));
                }
                Ok(())
            }
            Self::Append {
                expected_generation,
                action,
                ..
            } => {
                if *expected_generation >= MAX_EVENTS {
                    return Err(Error::Invalid("discussion generation"));
                }
                action.validate()
            }
        }
    }
    pub fn is_write(&self) -> bool {
        matches!(self, Self::Append { .. })
    }
}
dto! {pub struct Page {pub target:Target,pub generation:u64,pub events:Rows<Event,MAX_PAGE_ROWS>,pub next_offset:Option<u32>}}
dto! {pub enum Success {Identified,Listed(Page),Appended{target:Target,event:Event}}}
dto! {pub struct Failure {pub code:ErrorCode,pub actual_generation:Option<u64>}}
dto! {pub enum Outcome {Success(Success),Failure(Failure)}}
dto! {pub struct RequestEnvelope {pub request_id:u64,pub request:Request}}
dto! {pub struct ResponseEnvelope {pub request_id:u64,pub context:LabContext,pub member:Member,pub outcome:Outcome}}
impl RequestEnvelope {
    pub fn validate(&self) -> Result<(), Error> {
        if self.request_id == 0 {
            return Err(Error::Invalid("discussion correlation"));
        }
        self.request.validate()
    }
}
impl ResponseEnvelope {
    pub fn validate(&self) -> Result<(), Error> {
        if self.request_id == 0 {
            return Err(Error::Invalid("discussion correlation"));
        }
        validate_name(&self.context.name)?;
        self.member.author.validate()?;
        match &self.outcome {
            Outcome::Failure(f) => {
                if (f.code == ErrorCode::Conflict) != f.actual_generation.is_some()
                    || f.actual_generation.is_some_and(|n| n > MAX_EVENTS)
                {
                    return Err(Error::Invalid("discussion conflict"));
                }
            }
            Outcome::Success(Success::Identified) => (),
            Outcome::Success(Success::Appended { event, .. }) => event.validate()?,
            Outcome::Success(Success::Listed(p)) => {
                if p.generation > MAX_EVENTS
                    || p.next_offset
                        .is_some_and(|o| o == 0 || u64::from(o) >= p.generation)
                {
                    return Err(Error::Invalid("discussion page"));
                }
                for e in p.events.iter() {
                    e.validate()?;
                    if e.sequence > p.generation {
                        return Err(Error::Invalid("discussion sequence"));
                    }
                }
                if p.events
                    .windows(2)
                    .any(|w| w[0].sequence + 1 != w[1].sequence)
                {
                    return Err(Error::Invalid("discussion page order"));
                }
            }
        }
        Ok(())
    }
    pub fn validate_for(
        &self,
        request: &RequestEnvelope,
        expected_lab: &LabId,
    ) -> Result<(), Error> {
        self.validate()?;
        request.validate()?;
        if self.request_id != request.request_id || &self.context.lab_id != expected_lab {
            return Err(Error::ResponseMismatch);
        }
        let matches = match (&request.request, &self.outcome) {
            (_, Outcome::Failure(f)) => match f.actual_generation {
                None => true,
                Some(actual) => match &request.request {
                    Request::List {
                        expected_generation: Some(expected),
                        ..
                    }
                    | Request::Append {
                        expected_generation: expected,
                        ..
                    } => *expected != actual,
                    _ => false,
                },
            },
            (Request::Identify, Outcome::Success(Success::Identified)) => true,
            (
                Request::Append {
                    target,
                    expected_generation,
                    action,
                },
                Outcome::Success(Success::Appended { target: t, event }),
            ) => {
                target == t
                    && event.sequence == expected_generation + 1
                    && event.action == *action
                    && event.author == self.member.author
                    && event.role == self.member.role
                    && self.member.role != Role::Reader
            }
            (
                Request::List {
                    target,
                    offset,
                    limit,
                    expected_generation,
                },
                Outcome::Success(Success::Listed(p)),
            ) => {
                let start = u64::from(*offset).min(p.generation);
                let end = (u64::from(*offset) + u64::from(*limit)).min(p.generation);
                target == &p.target
                    && expected_generation.is_none_or(|g| g == p.generation)
                    && p.events.len() as u64 == end - start
                    && p.events
                        .iter()
                        .enumerate()
                        .all(|(i, e)| e.sequence == start + i as u64 + 1)
                    && p.next_offset == (end < p.generation).then_some(end as u32)
            }
            _ => false,
        };
        if !matches {
            return Err(Error::ResponseMismatch);
        }
        Ok(())
    }
}
pub fn encode_request(v: &RequestEnvelope) -> Result<Vec<u8>, Error> {
    v.validate()?;
    let b = super::encode(REQUEST_MAGIC, v)?;
    if b.len() > MAX_REQUEST_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(b)
}
pub fn decode_request(b: &[u8]) -> Result<RequestEnvelope, Error> {
    if b.len() > MAX_REQUEST_BYTES {
        return Err(Error::TooLarge);
    }
    let v: RequestEnvelope = super::decode(REQUEST_MAGIC, b)?;
    v.validate()?;
    Ok(v)
}
pub fn encode_response(v: &ResponseEnvelope) -> Result<Vec<u8>, Error> {
    v.validate()?;
    let b = super::encode(RESPONSE_MAGIC, v)?;
    if b.len() > MAX_RESPONSE_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(b)
}
pub fn decode_response(b: &[u8]) -> Result<ResponseEnvelope, Error> {
    if b.len() > MAX_RESPONSE_BYTES {
        return Err(Error::TooLarge);
    }
    let v: ResponseEnvelope = super::decode(RESPONSE_MAGIC, b)?;
    v.validate()?;
    Ok(v)
}
