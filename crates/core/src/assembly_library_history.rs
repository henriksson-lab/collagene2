//! Immutable recovery bytes for explicitly reviewed source-library edits.
//! Draft semantics belong to the host. Declared source pins are checked against
//! retained native events; neither a record nor its hash restores scientific evidence.
use crate::{
    Error,
    run_lineage::{self, EventBinding},
    workspace::{Workspace, WorkspaceAttachment},
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub const SCHEMA_VERSION: u16 = 2;
pub const KIND: &str = "assembly-library-history-v2";
pub const ACTIVE_DRAFT_VERSION: u32 = 19;
pub const ACTIVE_DRAFT_KIND: &str = "assembly-draft-v18";
pub const MEDIA_TYPE: &str = "application/x-collagene-assembly-library-history+postcard";
pub const MAX_RECORD_BYTES: usize = 32 << 20;
pub const MAX_HISTORY_BYTES: usize = 32 << 20;
pub const MAX_RECORDS: usize = 128;
pub const MAX_SOURCES: usize = 536;
pub const MAX_IMPACTS: usize = 512;
const MAX_BINS: usize = 12;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Impact {
    pub construct_id:String,pub name:String,
    #[serde(deserialize_with="decode_bins")]pub primer_review_bins:Vec<String>,
    #[serde(deserialize_with="decode_bins")]pub removed_preparation_bins:Vec<String>,
    pub retained_internal_edits:usize,pub recovered_internal_edits:usize,
    pub retained_junctions:usize,pub recovered_junctions:usize,
    #[serde(deserialize_with="decode_junctions")]pub new_junctions:Vec<usize>,
    pub overlap_text_changed:bool,pub role_order_changed:bool,pub candidate_cleared:bool,
    pub construct_removed:bool,pub construct_added:bool,
}}
dto! {pub struct Record {
    pub schema_version:u16,
    /// Persisted active bytes and raw working bytes can intentionally differ.
    pub before_active_hash:String,pub before_hash:String,
    #[serde(deserialize_with="decode_before")]pub before_draft:Vec<u8>,
    pub after_kind:String,pub after_hash:String,
    pub library_id:String,pub before_next_id:u64,pub after_next_id:u64,pub allocation_highwater:u64,
    pub proposal_revision:String,pub source_workspace_fingerprint:String,
    /// Retained declared pins only: the host must separately attest every capture
    /// in both decoded drafts, including removed/inactive choices and working slots.
    #[serde(deserialize_with="decode_sources")]pub sources:Vec<EventBinding>,
    #[serde(deserialize_with="decode_impacts")]pub impacts:Vec<Impact>,
}}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Summary {
    pub attachment_hash: String,
    pub proposal_revision: String,
    pub library_id: String,
    pub before_hash: String,
    pub after_hash: String,
    pub before_next_id: u64,
    pub after_next_id: u64,
    pub allocation_highwater: u64,
    pub impact_count: usize,
    pub source_count: usize,
}
dto! {pub struct Commit {
    #[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub workspace:Workspace,
    pub record_hash:String,pub after_hash:String,pub source_workspace_fingerprint:String,
}}
#[derive(Deserialize)]
struct Header {
    version: u32,
    alternatives: LibraryHeader,
}
#[derive(Deserialize)]
struct LibraryHeader {
    id: String,
    next_id: u64,
}
fn bad(message: &str) -> Error {
    Error::Format(format!("assembly library recovery: {message}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn id(s: &str, max: usize) -> bool {
    !s.is_empty()
        && s.len() <= max
        && s.bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.'))
}
fn header(bytes: &[u8]) -> Result<LibraryHeader, Error> {
    if bytes.is_empty() || bytes.len() > MAX_RECORD_BYTES {
        return Err(Error::TooLarge);
    }
    // Ignored values are skipped by serde; no WebDraft/Molecule is materialized.
    let h: Header = serde_json::from_slice(bytes)
        .map_err(|_| bad("raw draft requires a bounded current-version library header"))?;
    if h.version != ACTIVE_DRAFT_VERSION
        || !id(&h.alternatives.id, 48)
        || h.alternatives.next_id == 0
    {
        return Err(bad(
            "invalid current-version library identity/allocation header",
        ));
    }
    Ok(h.alternatives)
}
fn bounded<'de, D, T, const N: usize>(d: D) -> Result<Vec<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    struct V<T, const N: usize>(std::marker::PhantomData<T>);
    impl<'de, T: Deserialize<'de>, const N: usize> serde::de::Visitor<'de> for V<T, N> {
        type Value = Vec<T>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "at most {N} recovery values")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(self, mut a: A) -> Result<Vec<T>, A::Error> {
            if a.size_hint().is_some_and(|n| n > N) {
                return Err(serde::de::Error::custom("recovery value limit"));
            }
            let mut v = Vec::with_capacity(a.size_hint().unwrap_or(0).min(N));
            while v.len() < N {
                match a.next_element()? {
                    Some(value) => v.push(value),
                    None => return Ok(v),
                }
            }
            if a.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom("recovery value limit"));
            }
            Ok(v)
        }
    }
    d.deserialize_seq(V::<T, N>(std::marker::PhantomData))
}
fn decode_before<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<u8>, D::Error> {
    bounded::<D, u8, MAX_RECORD_BYTES>(d)
}
fn decode_sources<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<EventBinding>, D::Error> {
    bounded::<D, EventBinding, MAX_SOURCES>(d)
}
fn decode_impacts<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<Impact>, D::Error> {
    bounded::<D, Impact, MAX_IMPACTS>(d)
}
fn decode_bins<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<String>, D::Error> {
    bounded::<D, String, MAX_BINS>(d)
}
fn decode_junctions<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<usize>, D::Error> {
    bounded::<D, usize, MAX_BINS>(d)
}
fn hash(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}
impl Impact {
    fn validate(&self) -> Result<(), Error> {
        if !id(&self.construct_id, 48)
            || self.name.trim().is_empty()
            || self.name.len() > 256
            || self.name.chars().any(char::is_control)
            || self.construct_removed && self.construct_added
            || self.retained_internal_edits > 512
            || self.recovered_internal_edits > 512
            || self
                .retained_internal_edits
                .checked_add(self.recovered_internal_edits)
                .is_none_or(|total| total > 512)
            || self.retained_junctions > MAX_BINS
            || self.recovered_junctions > MAX_BINS
            || self.new_junctions.len() > MAX_BINS
            || self.new_junctions.iter().any(|n| *n >= MAX_BINS)
            || self.new_junctions.iter().collect::<BTreeSet<_>>().len() != self.new_junctions.len()
        {
            return Err(bad("invalid bounded row impact"));
        }
        for bins in [&self.primer_review_bins, &self.removed_preparation_bins] {
            if bins.len() > MAX_BINS
                || bins.iter().any(|b| !id(b, 48))
                || bins.iter().collect::<BTreeSet<_>>().len() != bins.len()
            {
                return Err(bad("invalid impact bin identity"));
            }
        }
        Ok(())
    }
}
impl Record {
    /// Structural historical data only. This cannot validate the host's complete
    /// draft schema, exhaustiveness of declared source pins, or scientific claims.
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != SCHEMA_VERSION
            || self.after_kind != ACTIVE_DRAFT_KIND
            || !id(&self.library_id, 48)
            || [
                &self.before_active_hash,
                &self.before_hash,
                &self.after_hash,
                &self.proposal_revision,
                &self.source_workspace_fingerprint,
            ]
            .iter()
            .any(|s| !digest(s))
            || self.sources.len() > MAX_SOURCES
            || self.impacts.len() > MAX_IMPACTS
            || self.before_next_id == 0
            || self.after_next_id == 0
            || self.after_next_id < self.before_next_id
            || self.allocation_highwater != self.before_next_id.max(self.after_next_id)
        {
            return Err(bad("invalid recovery record identity/dimensions"));
        }
        let h = header(&self.before_draft)?;
        if hash(&self.before_draft) != self.before_hash
            || h.id != self.library_id
            || h.next_id != self.before_next_id
        {
            return Err(bad("raw before draft hash/header changed"));
        }
        let mut pins = BTreeSet::new();
        for pin in &self.sources {
            if pin.document_id.is_empty()
                || pin.document_id.len() > 128
                || !pin
                    .document_id
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b':'))
                || !pin.event_id.strip_prefix("event-").is_some_and(digest)
                || !digest(&pin.revision)
                || !pins.insert((&pin.document_id, &pin.event_id, &pin.revision))
            {
                return Err(bad("invalid or duplicate retained source pin"));
            }
        }
        let mut rows = BTreeSet::new();
        for impact in &self.impacts {
            impact.validate()?;
            if !rows.insert(&impact.construct_id) {
                return Err(bad("duplicate impacted construct"));
            }
        }
        // Bounded native fields limit the allocation independently of input envelopes.
        if postcard::to_allocvec(self)
            .map_err(|_| Error::TooLarge)?
            .len()
            > MAX_RECORD_BYTES
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
pub fn encode_record(record: &Record) -> Result<Vec<u8>, Error> {
    record.validate()?;
    postcard::to_allocvec(record).map_err(|_| bad("could not encode recovery record"))
}
pub fn load_record(bytes: &[u8]) -> Result<Record, Error> {
    if bytes.len() > MAX_RECORD_BYTES {
        return Err(Error::TooLarge);
    }
    let (record, tail): (Record, &[u8]) =
        postcard::take_from_bytes(bytes).map_err(|_| bad("invalid recovery record encoding"))?;
    if !tail.is_empty() {
        return Err(bad("trailing recovery fields are not accepted"));
    }
    record.validate()?;
    Ok(record)
}
fn attest(w: &Workspace, record: &Record) -> Result<(), Error> {
    for pin in &record.sources {
        crate::combinatorial_design::resolve_event(w, pin)?;
    }
    Ok(())
}
fn active(w: &Workspace) -> Result<(&str, &WorkspaceAttachment), Error> {
    let mut matches = w
        .attachments
        .iter()
        .filter(|(_, a)| a.kind.starts_with("assembly-draft-v"));
    let (hash, a) = matches
        .next()
        .ok_or_else(|| bad("one persisted active current-version draft is required"))?;
    if matches.next().is_some() || a.kind != ACTIVE_DRAFT_KIND || a.media_type != "application/json"
    {
        return Err(bad(
            "exactly one active current-version JSON draft is required",
        ));
    }
    Ok((hash.as_str(), a))
}
fn records(w: &Workspace) -> Result<Vec<(String, Record)>, Error> {
    crate::combinatorial_design::workspace_guard(w)?;
    let mut total = 0usize;
    let mut result = Vec::new();
    for (hash, a) in &w.attachments {
        if a.kind != KIND {
            continue;
        }
        if result.len() >= MAX_RECORDS || a.media_type != MEDIA_TYPE {
            return Err(bad("recovery history count/media limit"));
        }
        total = total
            .checked_add(a.bytes.len())
            .filter(|n| *n <= MAX_HISTORY_BYTES)
            .ok_or(Error::TooLarge)?;
        let record = load_record(&a.bytes)?;
        attest(w, &record)?;
        result.push((hash.clone(), record));
    }
    Ok(result)
}
/// Loading attests declared original source events, independently of live heads.
/// The host must decode and re-attest all actual captures before recovery use.
pub fn load(w: &Workspace, attachment_hash: &str) -> Result<Record, Error> {
    if !digest(attachment_hash) {
        return Err(bad("invalid recovery attachment hash"));
    }
    records(w)?
        .into_iter()
        .find(|(hash, _)| hash == attachment_hash)
        .map(|(_, r)| r)
        .ok_or_else(|| bad("recovery attachment unavailable"))
}
pub fn list(w: &Workspace) -> Result<Vec<Summary>, Error> {
    Ok(records(w)?
        .into_iter()
        .map(|(attachment_hash, r)| Summary {
            attachment_hash,
            proposal_revision: r.proposal_revision,
            library_id: r.library_id,
            before_hash: r.before_hash,
            after_hash: r.after_hash,
            before_next_id: r.before_next_id,
            after_next_id: r.after_next_id,
            allocation_highwater: r.allocation_highwater,
            impact_count: r.impacts.len(),
            source_count: r.sources.len(),
        })
        .collect())
}
/// Preserve this lower bound when restoring raw historical drafts; restoration
/// is a new host-validated edit, never a rewind of identities or native science.
pub fn allocation_highwater(w: &Workspace, library_id: &str) -> Result<u64, Error> {
    if !id(library_id, 48) {
        return Err(bad("invalid library identity"));
    }
    let mut high = records(w)?
        .into_iter()
        .filter(|(_, r)| r.library_id == library_id)
        .map(|(_, r)| r.allocation_highwater)
        .max()
        .unwrap_or(1);
    let active_count = w
        .attachments
        .values()
        .filter(|a| a.kind.starts_with("assembly-draft-v"))
        .count();
    if active_count > 0 {
        let (_, a) = active(w)?;
        let h = header(&a.bytes)?;
        if h.id == library_id {
            high = high.max(h.next_id);
        }
    }
    Ok(high)
}
/// Host precondition: decode/validate BOTH full drafts and attest every capture
/// before calling. The native probe below only establishes identity/allocation.
pub fn create(
    w: &Workspace,
    before_raw: Vec<u8>,
    after: &[u8],
    proposal_revision: &str,
    mut sources: Vec<EventBinding>,
    impacts: Vec<Impact>,
) -> Result<Record, Error> {
    if sources.len() > MAX_SOURCES || impacts.len() > MAX_IMPACTS {
        return Err(Error::TooLarge);
    }
    if !digest(proposal_revision) {
        return Err(bad("invalid reviewed proposal revision"));
    }
    for impact in &impacts {
        impact.validate()?;
    }
    crate::combinatorial_design::workspace_guard(w)?;
    for source in &sources {
        crate::combinatorial_design::resolve_event(w, source)?;
    }
    let (before_active_hash, _) = active(w)?;
    let before = header(&before_raw)?;
    let next = header(after)?;
    if before.id != next.id
        || next.next_id < before.next_id
        || next.next_id < allocation_highwater(w, &before.id)?
    {
        return Err(bad(
            "library identity changed or allocation counter regressed",
        ));
    }
    sources.sort_by(|a, b| {
        (&a.document_id, &a.event_id, &a.revision).cmp(&(&b.document_id, &b.event_id, &b.revision))
    });
    sources.dedup_by(|a, b| {
        a.document_id == b.document_id && a.event_id == b.event_id && a.revision == b.revision
    });
    let record = Record {
        schema_version: SCHEMA_VERSION,
        before_active_hash: before_active_hash.into(),
        before_hash: hash(&before_raw),
        before_draft: before_raw,
        after_kind: ACTIVE_DRAFT_KIND.into(),
        after_hash: hash(after),
        library_id: before.id,
        before_next_id: before.next_id,
        after_next_id: next.next_id,
        allocation_highwater: before.next_id.max(next.next_id),
        proposal_revision: proposal_revision.into(),
        source_workspace_fingerprint: run_lineage::workspace_fingerprint(w)?,
        sources,
        impacts,
    };
    record.validate()?;
    attest(w, &record)?;
    Ok(record)
}
/// Consumes a private frame and stages both attachments together. The host must
/// install this complete result in one transaction after final live-context checks.
pub fn install(mut w: Workspace, record: Record, after: Vec<u8>) -> Result<Commit, Error> {
    crate::combinatorial_design::workspace_guard(&w)?;
    record.validate()?;
    attest(&w, &record)?;
    if run_lineage::workspace_fingerprint(&w)? != record.source_workspace_fingerprint
        || active(&w)?.0 != record.before_active_hash
    {
        return Err(bad(
            "workspace or persisted active draft changed before recovery installation",
        ));
    }
    let h = header(&after)?;
    if hash(&after) != record.after_hash
        || h.id != record.library_id
        || h.next_id != record.after_next_id
        || h.next_id < allocation_highwater(&w, &h.id)?
    {
        return Err(bad("replacement draft/header differs or counter regressed"));
    }
    let existing = records(&w)?;
    let bytes = encode_record(&record)?;
    let record_hash = hash(&bytes);
    let total = w
        .attachments
        .values()
        .filter(|a| a.kind == KIND)
        .try_fold(bytes.len(), |n, a| {
            n.checked_add(a.bytes.len()).ok_or(Error::TooLarge)
        })?;
    if existing.len() >= MAX_RECORDS
        || total > MAX_HISTORY_BYTES
        || w.attachments.len() >= crate::workspace::MAX_ATTACHMENTS
    {
        return Err(Error::TooLarge);
    }
    if w.attachments
        .get(&record.after_hash)
        .is_some_and(|a| a.kind != ACTIVE_DRAFT_KIND)
        || w.attachments.contains_key(&record_hash)
    {
        return Err(bad(
            "recovery/replacement hash already belongs to a retained attachment",
        ));
    }
    // No operation above mutates the frame. History retention cannot be pruned to
    // make an edit succeed; any failure publishes neither replacement nor history.
    w.attachments.retain(|_, a| a.kind != ACTIVE_DRAFT_KIND);
    w.attachments.insert(
        record.after_hash.clone(),
        WorkspaceAttachment {
            kind: ACTIVE_DRAFT_KIND.into(),
            name: "assembly-draft.json".into(),
            media_type: "application/json".into(),
            bytes: after,
        },
    );
    w.attachments.insert(
        record_hash.clone(),
        WorkspaceAttachment {
            kind: KIND.into(),
            name: format!("assembly-library-{record_hash}.postcard"),
            media_type: MEDIA_TYPE.into(),
            bytes,
        },
    );
    crate::combinatorial_design::workspace_guard(&w)?;
    Ok(Commit {
        workspace: w,
        record_hash,
        after_hash: record.after_hash,
        source_workspace_fingerprint: record.source_workspace_fingerprint,
    })
}
