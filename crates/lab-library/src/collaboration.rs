//! Revision-bound human discussion. These records never carry scientific approval.
//! Member capabilities come only from the host's authenticated identity registry.
use super::*;

pub const MAX_BODY_BYTES: usize = 8192;
pub const MAX_REVISION_EVENTS: u64 = 4096;
pub const MAX_LAB_EVENTS: u64 = 100_000;
pub const MAX_LAB_DISCUSSION_BYTES: u64 = 16 * 1024 * 1024;

#[derive(Debug, thiserror::Error)]
pub enum DiscussionError {
    #[error("member is not permitted to perform this discussion action")]
    Forbidden,
    #[error("discussion changed: expected generation {expected}, current {actual}")]
    Conflict { expected: u64, actual: u64 },
    #[error(transparent)]
    Storage(#[from] super::Error),
    #[error(transparent)]
    Database(#[from] rusqlite::Error),
}
pub type DiscussionResult<T> = std::result::Result<T, DiscussionError>;

/// Stable identifier in the operator's identity registry, never request-supplied.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(try_from = "String", into = "String")]
pub struct ActorId(String);
impl ActorId {
    pub fn parse(value: impl Into<String>) -> super::Result<Self> {
        let value = value.into();
        if value.is_empty()
            || value.len() > 128
            || !value
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"-_.:@".contains(&b))
        {
            return Err(Error::Invalid("actor ID"));
        }
        Ok(Self(value))
    }
    pub fn as_str(&self) -> &str {
        &self.0
    }
}
impl TryFrom<String> for ActorId {
    type Error = Error;
    fn try_from(value: String) -> super::Result<Self> {
        Self::parse(value)
    }
}
impl From<ActorId> for String {
    fn from(value: ActorId) -> Self {
        value.0
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Author {
    pub id: ActorId,
    /// Historical display name at the time of posting, not an identity key.
    pub display_name: String,
}
impl Author {
    pub fn new(id: ActorId, display_name: impl Into<String>) -> super::Result<Self> {
        let display_name = display_name.into();
        if display_name.trim().is_empty()
            || display_name.len() > 256
            || display_name.chars().any(char::is_control)
        {
            return Err(Error::Invalid("author display name"));
        }
        Ok(Self { id, display_name })
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Role {
    Reader,
    Contributor,
    Maintainer,
}

/// Non-wire capability. It intentionally cannot be deserialized from a request.
#[derive(Debug, Clone)]
pub struct MemberAccess {
    lab: LabId,
    author: Author,
    role: Role,
}
impl MemberAccess {
    pub fn lab_id(&self) -> &LabId {
        &self.lab
    }
    pub fn author(&self) -> &Author {
        &self.author
    }
    pub fn role(&self) -> Role {
        self.role
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Target {
    pub document: DocumentId,
    pub revision: RevisionId,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReviewDecision {
    Approved,
    ChangesRequested,
    Withdrawn,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub enum Action {
    /// This event's sequence number is the new thread ID.
    OpenThread {
        body: String,
    },
    Reply {
        thread: u64,
        body: String,
    },
    SetResolved {
        thread: u64,
        resolved: bool,
    },
    /// Append-only; a later decision supersedes only this author's previous
    /// decision on this exact library document/revision, never other revisions.
    Review {
        decision: ReviewDecision,
        body: String,
    },
}
impl Action {
    fn validate(&self) -> super::Result<()> {
        match self {
            Self::OpenThread { body } | Self::Reply { body, .. } => body_valid(body, false)?,
            Self::Review { body, .. } => body_valid(body, true)?,
            Self::SetResolved { .. } => (),
        }
        if matches!(
            self,
            Self::Reply { thread: 0, .. } | Self::SetResolved { thread: 0, .. }
        ) {
            return Err(Error::Invalid("discussion thread"));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Event {
    pub sequence: u64,
    pub author: Author,
    pub role: Role,
    pub created_at: u64,
    pub action: Action,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Page {
    pub target: Target,
    pub generation: u64,
    pub events: Vec<Event>,
    pub next_offset: Option<u32>,
}
fn body_valid(body: &str, allow_empty: bool) -> super::Result<()> {
    if body.len() > MAX_BODY_BYTES
        || (!allow_empty && body.trim().is_empty())
        || body
            .chars()
            .any(|c| c.is_control() && c != '\n' && c != '\t')
    {
        return Err(Error::Invalid("discussion body"));
    }
    Ok(())
}
fn target_exists(c: &Connection, lab: &LabId, target: &Target) -> super::Result<()> {
    get_document(c, lab, &target.document)?;
    get_revision(c, &target.document, &target.revision)?;
    Ok(())
}
fn generation(c: &Connection, target: &Target) -> super::Result<u64> {
    Ok(c.query_row(
        "SELECT generation FROM discussion_heads WHERE document_id=?1 AND revision_id=?2",
        params![target.document.as_str(), target.revision.as_str()],
        |r| r.get(0),
    )
    .optional()?
    .unwrap_or(0))
}
fn event(c: &Connection, target: &Target, sequence: u64) -> super::Result<Event> {
    let (author, display_name, role, created_at, action): (String, String, String, u64, String) = c.query_row(
        "SELECT actor_id,display_name,role,created_at,action FROM discussion_events WHERE document_id=?1 AND revision_id=?2 AND sequence=?3",
        params![target.document.as_str(), target.revision.as_str(), sequence],
        |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get(3)?, r.get(4)?)))
        .optional()?.ok_or(Error::NotFound)?;
    let role: Role = serde_json::from_str(&role).map_err(|_| Error::Corrupt("discussion role"))?;
    let action: Action =
        serde_json::from_str(&action).map_err(|_| Error::Corrupt("discussion action"))?;
    action
        .validate()
        .map_err(|_| Error::Corrupt("discussion action"))?;
    let author = Author::new(ActorId::parse(author)?, display_name)?;
    if role == Role::Reader || sequence == 0 || sequence > MAX_REVISION_EVENTS {
        return Err(Error::Corrupt("discussion event"));
    }
    Ok(Event {
        sequence,
        author,
        role,
        created_at,
        action,
    })
}
fn check_action(
    c: &Connection,
    target: &Target,
    member: &MemberAccess,
    action: &Action,
) -> DiscussionResult<()> {
    if member.role == Role::Reader {
        return Err(DiscussionError::Forbidden);
    }
    let thread = match action {
        Action::Reply { thread, .. } | Action::SetResolved { thread, .. } => *thread,
        _ => return Ok(()),
    };
    let root = event(c, target, thread)?;
    if !matches!(root.action, Action::OpenThread { .. }) {
        return Err(Error::Invalid("discussion thread root").into());
    }
    let mut resolved = false;
    let mut query = c.prepare("SELECT sequence FROM discussion_events WHERE document_id=?1 AND revision_id=?2 AND sequence>?3 ORDER BY sequence")?;
    for row in query.query_map(
        params![target.document.as_str(), target.revision.as_str(), thread],
        |r| r.get::<_, u64>(0),
    )? {
        if let Action::SetResolved {
            thread: t,
            resolved: value,
        } = event(c, target, row?)?.action
        {
            if t == thread {
                resolved = value;
            }
        }
    }
    match action {
        Action::Reply { .. } if resolved => {
            return Err(Error::Invalid("resolved discussion thread").into());
        }
        Action::SetResolved {
            resolved: value, ..
        } => {
            if member.role != Role::Maintainer && member.author.id != root.author.id {
                return Err(DiscussionError::Forbidden);
            }
            if *value == resolved {
                return Err(Error::Invalid("unchanged thread resolution").into());
            }
        }
        _ => (),
    }
    Ok(())
}
impl Library {
    /// Host administration only. The caller must authenticate an individual and
    /// resolve their current role from trusted configuration before calling.
    /// LabAccess remains the separate operator/storage capability; it is not a
    /// substitute for a member's identity and cannot author discussion events.
    pub fn trusted_member_access(
        &self,
        access: &LabAccess,
        author: Author,
        role: Role,
    ) -> super::Result<MemberAccess> {
        ensure_lab(&self.connection, &access.lab)?;
        let author = Author::new(ActorId::parse(author.id.as_str())?, author.display_name)?;
        Ok(MemberAccess {
            lab: access.lab.clone(),
            author,
            role,
        })
    }
    pub fn discussion_page(
        &mut self,
        member: &MemberAccess,
        target: &Target,
        offset: u32,
        limit: u16,
        expected_generation: Option<u64>,
    ) -> DiscussionResult<Page> {
        pagination(offset, limit)?;
        continuation(offset, expected_generation)?;
        if u64::from(offset) > MAX_REVISION_EVENTS {
            return Err(Error::Invalid("discussion offset").into());
        }
        let tx = self.connection.transaction()?;
        target_exists(&tx, &member.lab, target)?;
        let current = generation(&tx, target)?;
        if let Some(expected) = expected_generation {
            if current != expected {
                return Err(DiscussionError::Conflict {
                    expected,
                    actual: current,
                });
            }
        }
        let end = (u64::from(offset) + u64::from(limit)).min(current);
        let mut events = Vec::new();
        for sequence in (u64::from(offset) + 1)..=end {
            events.push(event(&tx, target, sequence)?);
        }
        tx.commit()?;
        Ok(Page {
            target: target.clone(),
            generation: current,
            events,
            next_offset: (end < current).then_some(end as u32),
        })
    }
    /// Compare-and-append in one immediate transaction. Retry after a lost reply
    /// produces a conflict, never a duplicate comment/approval. Read the exact
    /// generation again before deciding to resubmit.
    pub fn append_discussion(
        &mut self,
        member: &MemberAccess,
        target: &Target,
        expected_generation: u64,
        action: Action,
    ) -> DiscussionResult<Event> {
        action.validate()?;
        let encoded =
            serde_json::to_string(&action).map_err(|_| Error::Invalid("discussion action"))?;
        let tx = self
            .connection
            .transaction_with_behavior(TransactionBehavior::Immediate)?;
        target_exists(&tx, &member.lab, target)?;
        if member.role == Role::Reader {
            return Err(DiscussionError::Forbidden);
        }
        let current = generation(&tx, target)?;
        if current != expected_generation {
            return Err(DiscussionError::Conflict {
                expected: expected_generation,
                actual: current,
            });
        }
        if current >= MAX_REVISION_EVENTS {
            return Err(Error::Limit("revision discussion events").into());
        }
        check_action(&tx, target, member, &action)?;
        let (count, bytes): (u64, u64) = tx.query_row(
            "SELECT discussion_events,discussion_bytes FROM labs WHERE id=?1",
            [member.lab.as_str()],
            |r| Ok((r.get(0)?, r.get(1)?)),
        )?;
        let role =
            serde_json::to_string(&member.role).map_err(|_| Error::Invalid("discussion role"))?;
        let charge = encoded.len() as u64
            + member.author.id.as_str().len() as u64
            + member.author.display_name.len() as u64
            + role.len() as u64;
        if count >= MAX_LAB_EVENTS || bytes.saturating_add(charge) > MAX_LAB_DISCUSSION_BYTES {
            return Err(Error::Limit("lab discussion quota").into());
        }
        tx.execute("INSERT INTO discussion_heads(document_id,revision_id,generation) VALUES(?1,?2,?3) ON CONFLICT(document_id,revision_id) DO UPDATE SET generation=excluded.generation",
            params![target.document.as_str(), target.revision.as_str(), current + 1])?;
        tx.execute("INSERT INTO discussion_events(document_id,revision_id,sequence,actor_id,display_name,role,action,charged_bytes) VALUES(?1,?2,?3,?4,?5,?6,?7,?8)",
            params![target.document.as_str(), target.revision.as_str(), current + 1, member.author.id.as_str(), &member.author.display_name, role, encoded, charge])?;
        tx.execute("UPDATE labs SET discussion_events=discussion_events+1,discussion_bytes=discussion_bytes+?2 WHERE id=?1", params![member.lab.as_str(), charge])?;
        let result = event(&tx, target, current + 1)?;
        tx.commit()?;
        Ok(result)
    }
}

/// Backup verification includes the event grammar, authorship permissions,
/// append order and quota accounting; foreign keys alone do not establish these.
pub(crate) fn check_integrity(c: &Connection) -> super::Result<()> {
    let mut labs = c.prepare("SELECT id,discussion_events,discussion_bytes FROM labs")?;
    for row in labs.query_map([], |r| {
        Ok((
            r.get::<_, String>(0)?,
            r.get::<_, u64>(1)?,
            r.get::<_, u64>(2)?,
        ))
    })? {
        let (lab, count, bytes) = row?;
        let (actual_count, actual_bytes): (u64, u64) = c.query_row(
            "SELECT count(*),coalesce(sum(e.charged_bytes),0) FROM discussion_events e JOIN documents d ON d.id=e.document_id WHERE d.lab_id=?1", [&lab], |r| Ok((r.get(0)?, r.get(1)?)))?;
        if count != actual_count
            || bytes != actual_bytes
            || count > MAX_LAB_EVENTS
            || bytes > MAX_LAB_DISCUSSION_BYTES
        {
            return Err(Error::Corrupt("discussion quota accounting"));
        }
    }
    let mut heads = c.prepare("SELECT h.document_id,h.revision_id,h.generation,d.lab_id FROM discussion_heads h JOIN documents d ON d.id=h.document_id")?;
    for row in heads.query_map([], |r| {
        Ok((
            r.get::<_, String>(0)?,
            r.get::<_, String>(1)?,
            r.get::<_, u64>(2)?,
            r.get::<_, String>(3)?,
        ))
    })? {
        let (document, revision, generation, lab) = row?;
        let target = Target {
            document: DocumentId::parse(document)?,
            revision: RevisionId::parse(revision)?,
        };
        let lab = LabId::parse(lab)?;
        target_exists(c, &lab, &target)?;
        let (count, last): (u64, u64) = c.query_row("SELECT count(*),coalesce(max(sequence),0) FROM discussion_events WHERE document_id=?1 AND revision_id=?2",
            params![target.document.as_str(), target.revision.as_str()], |r| Ok((r.get(0)?, r.get(1)?)))?;
        if generation == 0
            || generation != count
            || generation != last
            || generation > MAX_REVISION_EVENTS
        {
            return Err(Error::Corrupt("discussion generation"));
        }
        // Replay only earlier events when checking whether a resolution/reply
        // was authorized at that time (later reopen events cannot change it).
        let mut threads = std::collections::BTreeMap::<u64, (ActorId, bool)>::new();
        for sequence in 1..=generation {
            let e = event(c, &target, sequence)?;
            let encoded = serde_json::to_string(&e.action)
                .map_err(|_| Error::Corrupt("discussion action"))?;
            let role =
                serde_json::to_string(&e.role).map_err(|_| Error::Corrupt("discussion role"))?;
            let expected = encoded.len()
                + e.author.id.as_str().len()
                + e.author.display_name.len()
                + role.len();
            let charge: u64 = c.query_row("SELECT charged_bytes FROM discussion_events WHERE document_id=?1 AND revision_id=?2 AND sequence=?3",
                params![target.document.as_str(), target.revision.as_str(), sequence], |r| r.get(0))?;
            if charge != expected as u64 {
                return Err(Error::Corrupt("discussion event charge"));
            }
            match e.action {
                Action::OpenThread { .. } => {
                    threads.insert(sequence, (e.author.id, false));
                }
                Action::Reply { thread, .. } => {
                    if !matches!(threads.get(&thread), Some((_, false))) {
                        return Err(Error::Corrupt("discussion reply"));
                    }
                }
                Action::SetResolved { thread, resolved } => {
                    let (author, state) = threads
                        .get_mut(&thread)
                        .ok_or(Error::Corrupt("discussion thread"))?;
                    if *state == resolved || (e.role != Role::Maintainer && *author != e.author.id)
                    {
                        return Err(Error::Corrupt("discussion resolution"));
                    }
                    *state = resolved;
                }
                Action::Review { .. } => (),
            }
        }
    }
    Ok(())
}

pub(crate) const SCHEMA: &str = r#"
CREATE TABLE discussion_heads(
 document_id TEXT NOT NULL,
 revision_id TEXT NOT NULL,
 generation INTEGER NOT NULL CHECK(typeof(generation)='integer' AND generation BETWEEN 1 AND 4096),
 PRIMARY KEY(document_id,revision_id),
 FOREIGN KEY(document_id,revision_id) REFERENCES revisions(document_id,id)
);
CREATE TABLE discussion_events(
 document_id TEXT NOT NULL,
 revision_id TEXT NOT NULL,
 sequence INTEGER NOT NULL CHECK(typeof(sequence)='integer' AND sequence BETWEEN 1 AND 4096),
 actor_id TEXT NOT NULL CHECK(length(actor_id) BETWEEN 1 AND 128),
 display_name TEXT NOT NULL CHECK(length(display_name) BETWEEN 1 AND 256),
 role TEXT NOT NULL CHECK(length(role)<=32),
 action TEXT NOT NULL CHECK(length(CAST(action AS BLOB))<=50000),
 charged_bytes INTEGER NOT NULL CHECK(typeof(charged_bytes)='integer' AND charged_bytes BETWEEN 1 AND 51000),
 created_at INTEGER NOT NULL DEFAULT (unixepoch()) CHECK(created_at>=0),
 PRIMARY KEY(document_id,revision_id,sequence),
 FOREIGN KEY(document_id,revision_id) REFERENCES discussion_heads(document_id,revision_id)
);
CREATE TRIGGER immutable_discussion_update BEFORE UPDATE ON discussion_events BEGIN SELECT RAISE(ABORT,'immutable discussion event'); END;
CREATE TRIGGER immutable_discussion_delete BEFORE DELETE ON discussion_events BEGIN SELECT RAISE(ABORT,'immutable discussion event'); END;
"#;
