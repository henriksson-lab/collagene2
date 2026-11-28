//! Explicit caller decisions over fresh read evidence. No original call or source
//! document is edited, and a reviewed reference is never physical sample proof.
use crate::documents::history::transition::{AnalysisOrigin, AnalysisProducer, Transition};
use crate::{
    documents::{
        history::DocumentHistory,
        workspace::{Workspace, WorkspaceAttachment},
    },
    material::lineage::{self as run_lineage, EventBinding},
    review::reads::{
        checks as read_checks,
        trace_validation::{ColumnKind, SiteState, TraceValidationReport},
    },
};
use collagene_model::{Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ATTACHMENT_KIND: &str = "read-review-v1";
pub const MAX_REVIEWS: usize = 128;
const MAX_BYTES: usize = 32 << 20;
const MAX_DECISIONS: usize = 256;
const REFERENCE_MODEL: &str =
    "Reference-only human read-review decision; physical material state unavailable";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq,PartialOrd,Ord)] pub enum EvidenceKey {Site{position:usize},Insertion{index:usize}}}
dto! {#[derive(PartialEq,Eq)] pub enum Choice {RetainReference,Unresolved,UseBase{read_index:usize,call_index:usize,base:char},DeleteBase{read_index:usize},UseInsertion{read_index:usize,boundary:usize,sequence:String}}}
dto! {#[derive(PartialEq,Eq)] pub struct Decision {pub key:EvidenceKey,pub choice:Choice,pub reason:String}}
dto! {#[derive(PartialEq,Eq)] pub struct Request {pub attachment_hash:String,pub name:String,pub decisions:Vec<Decision>}}
dto! {#[derive(PartialEq,Eq)] pub struct Edit {pub start:usize,pub end:usize,pub replacement:String,pub reason:String}}
dto! {pub struct Preview {pub context_revision:String,pub request:Request,pub reference:EventBinding,pub reads:Vec<EventBinding>,pub target:Option<Molecule>,pub edits:Vec<Edit>,pub warnings:Vec<String>}}
dto! {pub struct SavedReview {pub schema_version:u32,pub context_revision:String,pub request:Request,pub reference:EventBinding,pub reads:Vec<EventBinding>,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))] pub report:TraceValidationReport,pub edits:Vec<Edit>,pub warnings:Vec<String>,pub target:Option<EventBinding>}}
dto! {pub struct Commit {#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))] pub workspace:Workspace,pub context_revision:String,pub attachment_hash:String,pub target:Option<EventBinding>}}
dto! {pub struct Summary {pub attachment_hash:String,pub reference:EventBinding,pub decision_count:usize,pub change_count:usize,pub target:Option<EventBinding>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("read review: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn request_guard(r: &Request) -> Result<(), Error> {
    if !hex(&r.attachment_hash)
        || r.name.trim().is_empty()
        || r.name.len() > 256
        || r.name.chars().any(char::is_control)
        || r.decisions.is_empty()
        || r.decisions.len() > MAX_DECISIONS
    {
        return Err(bad(
            "bounded check hash, name and1..256 explicit decisions are required",
        ));
    }
    let mut seen = BTreeSet::new();
    let mut inserted = 0usize;
    for d in &r.decisions {
        if d.reason.trim().is_empty()
            || d.reason.len() > 1024
            || d.reason.chars().any(char::is_control)
            || !seen.insert(&d.key)
        {
            return Err(bad(
                "each evidence key requires one explicit decision and a bounded reason",
            ));
        }
        match &d.choice {
            Choice::UseBase { base, .. } if !matches!(base, 'A' | 'C' | 'G' | 'T') => {
                return Err(bad("expected allele must be concrete uppercase DNA"));
            }
            Choice::UseInsertion { sequence, .. } => {
                inserted = inserted
                    .checked_add(sequence.len())
                    .filter(|n| *n <= 100_000)
                    .ok_or(Error::TooLarge)?;
                if sequence.is_empty() || !sequence.bytes().all(|b| b"ACGT".contains(&b)) {
                    return Err(bad("expected insertion must be concrete uppercase DNA"));
                }
            }
            _ => {}
        }
    }
    Ok(())
}
struct Bytes {
    bytes: Vec<u8>,
}
impl Write for Bytes {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > MAX_BYTES)
        {
            return Err(std::io::Error::other("review output cap"));
        }
        self.bytes.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encode(v: &impl Serialize) -> Result<Vec<u8>, Error> {
    let mut w = Bytes { bytes: Vec::new() };
    serde_json::to_writer(&mut w, v).map_err(|_| Error::TooLarge)?;
    Ok(w.bytes)
}
fn json(v: &impl Serialize) -> Result<String, Error> {
    String::from_utf8(encode(v)?).map_err(|_| bad("invalid evidence encoding"))
}
struct Count(usize);
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| std::io::Error::other("source cap"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
/// Borrowed workspace, source, annotation-work and caller-text guards. No read
/// alignment is performed here; worker capture/commit must replay before use.
pub fn source_preflight(w: &Workspace, r: &Request) -> Result<(), Error> {
    request_guard(r)?;
    read_checks::workspace_guard(w)?;
    if w.attachments
        .values()
        .filter(|a| a.kind == ATTACHMENT_KIND)
        .count()
        > MAX_REVIEWS
    {
        return Err(Error::TooLarge);
    }
    let saved = read_checks::load(w, &r.attachment_hash)?;
    let source = read_checks::resolve(w, &saved.reference)?;
    if source.sequence.len() > 8 << 20 {
        return Err(Error::TooLarge);
    }
    collagene_model::bounds::model_preflight(&[source])?;
    let mut bytes = Count(0);
    serde_json::to_writer(&mut bytes, source).map_err(|_| Error::TooLarge)?;
    if bytes
        .0
        .checked_mul(
            r.decisions
                .iter()
                .filter(|d| !matches!(d.choice, Choice::RetainReference | Choice::Unresolved))
                .count()
                + 1,
        )
        .is_none_or(|n| n > 128 << 20)
    {
        return Err(Error::TooLarge);
    }
    let mut replaced = 0usize;
    for d in &r.decisions {
        if let EvidenceKey::Site { position } = d.key {
            if position >= source.sequence.len() {
                return Err(bad("site outside original reference"));
            }
        }
        match &d.choice {
            Choice::UseBase { read_index, .. } | Choice::DeleteBase { read_index }
                if *read_index >= saved.reads.len() =>
            {
                return Err(bad("unknown nominated read"));
            }
            Choice::UseInsertion {
                read_index,
                boundary,
                ..
            } if *read_index >= saved.reads.len() || *boundary > source.sequence.len() => {
                return Err(bad("invalid insertion source coordinates"));
            }
            _ => {}
        }

        if let (EvidenceKey::Insertion { index }, Choice::UseInsertion { .. }) = (&d.key, &d.choice)
        {
            let i = saved
                .report
                .insertions
                .get(*index)
                .ok_or_else(|| bad("unknown insertion evidence"))?;
            replaced = replaced
                .checked_add(i.sequence.len())
                .filter(|n| *n <= 100_000)
                .ok_or(Error::TooLarge)?;
        }
    }
    Ok(())
}
fn revision(w: &Workspace, r: &Request, domain: &str) -> Result<String, Error> {
    source_preflight(w, r)?;
    crate::encoding::bounded_json::hash(&(
        domain,
        run_lineage::workspace_fingerprint(w)?,
        r,
        crate::verification::ENGINE_VERSION,
        crate::review::reads::trace_validation::MODEL,
    ))
}
pub fn preview_revision(w: &Workspace, r: &Request) -> Result<String, Error> {
    revision(w, r, "preview-read-review-1")
}
pub fn commit_revision(w: &Workspace, r: &Request) -> Result<String, Error> {
    revision(w, r, "commit-read-review-1")
}
fn stable_revision(
    r: &Request,
    reference: &EventBinding,
    reads: &[EventBinding],
    report: &TraceValidationReport,
) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(
        "read-review-evidence-1",
        r,
        reference,
        reads,
        crate::encoding::bounded_json::hash(report)?,
    ))
}
fn warnings_guard(warnings: &[String]) -> Result<(), Error> {
    if warnings.len() > 1024 || warnings.iter().any(|w| w.len() > 4096) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn choose(
    source: &Molecule,
    r: &Request,
    report: &TraceValidationReport,
) -> Result<(Vec<Edit>, Vec<String>), Error> {
    let n = source.sequence.len();
    let mut edits = Vec::new();
    let mut boundaries: BTreeMap<usize, bool> = BTreeMap::new();
    let mut total = 0usize;
    let mut insertion_calls = BTreeMap::new();
    for (read_index, read) in report.reads.iter().enumerate() {
        for c in &read.columns {
            if c.kind == ColumnKind::Insertion {
                if let Some(call) = c.call_index {
                    if insertion_calls
                        .insert((read_index, c.reference_boundary, call), c)
                        .is_some()
                    {
                        return Err(bad("insertion call identity is repeated"));
                    }
                }
            }
        }
    }
    let mut warnings=vec!["Explicit caller decisions over freshly recomputed evidence; these do not establish consensus, sample identity, physical material state or experimental success.".into()];
    for d in &r.decisions {
        let edit = match (&d.key, &d.choice) {
            (EvidenceKey::Site { position }, choice) => {
                if *position >= n {
                    return Err(bad("site outside original reference"));
                }
                match choice {
                    Choice::RetainReference | Choice::Unresolved => None,
                    Choice::UseBase {
                        read_index,
                        call_index,
                        base: expected_base,
                    } => {
                        let read = report
                            .reads
                            .get(*read_index)
                            .ok_or_else(|| bad("unknown nominated read"))?;
                        let mut candidates = read.columns.iter().filter(|c| {
                            c.kind == ColumnKind::Pair
                                && c.reference_position == Some(*position)
                                && c.call_index == Some(*call_index)
                                && c.eligible
                        });
                        let call = candidates.next().ok_or_else(|| {
                            bad("nominated eligible current call does not support this site")
                        })?;
                        if candidates.next().is_some() {
                            return Err(bad("nominated call has ambiguous coordinate identity"));
                        }
                        let base = call
                            .observed_base
                            .filter(|b| matches!(b, 'A' | 'C' | 'G' | 'T'))
                            .ok_or_else(|| {
                                bad("a selected allele must be one concrete observed DNA call")
                            })?;
                        if base != *expected_base {
                            return Err(bad(
                                "fresh observed allele differs from the explicitly requested base",
                            ));
                        }
                        if report.sites.iter().any(|s| {
                            s.position == *position && s.state == SiteState::ConflictingCalls
                        }) {
                            warnings.push(format!("Site {} selects one eligible conflicting allele; the other observed evidence is retained and this choice is not consensus.",position));
                        }
                        (source.sequence.as_str().as_bytes()[*position] != base as u8).then(|| {
                            Edit {
                                start: *position,
                                end: *position + 1,
                                replacement: base.to_string(),
                                reason: d.reason.clone(),
                            }
                        })
                    }
                    Choice::DeleteBase { read_index } => {
                        let read = report
                            .reads
                            .get(*read_index)
                            .ok_or_else(|| bad("unknown nominated read"))?;
                        if !read.columns.iter().any(|c| {
                            c.kind == ColumnKind::Deletion
                                && c.reference_position == Some(*position)
                                && c.eligible
                        }) {
                            return Err(bad(
                                "no eligible deletion observation at the nominated site",
                            ));
                        }
                        Some(Edit {
                            start: *position,
                            end: *position + 1,
                            replacement: String::new(),
                            reason: d.reason.clone(),
                        })
                    }
                    Choice::UseInsertion { .. } => {
                        return Err(bad("insertion choice requires insertion evidence"));
                    }
                }
            }
            (EvidenceKey::Insertion { index }, choice) => {
                let i = report
                    .insertions
                    .get(*index)
                    .ok_or_else(|| bad("unknown insertion evidence"))?;
                let canonical = if source.topology == Topology::Circular {
                    i.reference_boundary % n
                } else {
                    i.reference_boundary
                };
                let uses = matches!(choice, Choice::UseInsertion { .. });
                if let Some(previous) = boundaries.insert(canonical, uses) {
                    if previous || uses {
                        return Err(bad(
                            "contradictory insertion decisions share one boundary, including circular0/length",
                        ));
                    }
                }
                match choice {
                    Choice::RetainReference | Choice::Unresolved => None,
                    Choice::UseInsertion {
                        read_index,
                        boundary,
                        sequence,
                    } => {
                        if *read_index != i.read_index
                            || *boundary != i.reference_boundary
                            || sequence != &i.sequence
                        {
                            return Err(bad(
                                "fresh insertion read, boundary or sequence differs from the explicit request",
                            ));
                        }
                        if !i.eligible
                            || i.sequence.is_empty()
                            || !i.sequence.bytes().all(|b| b"ACGT".contains(&b))
                            || i.sequence.len() != i.call_indices.len()
                        {
                            return Err(bad(
                                "insertion requires its complete eligible concrete observed sequence",
                            ));
                        }
                        for (&index, base) in i.call_indices.iter().zip(i.sequence.bytes()) {
                            if !insertion_calls
                                .get(&(i.read_index, i.reference_boundary, index))
                                .is_some_and(|c| {
                                    c.observed_base == Some(base as char) && c.eligible
                                })
                            {
                                return Err(bad(
                                    "insertion is not supported by all original eligible calls",
                                ));
                            }
                        }
                        Some(Edit {
                            start: i.reference_boundary,
                            end: i.reference_boundary,
                            replacement: i.sequence.clone(),
                            reason: d.reason.clone(),
                        })
                    }
                    _ => {
                        return Err(bad(
                            "insertion evidence accepts only insertion, retain or unresolved decisions",
                        ));
                    }
                }
            }
        };
        if let Some(e) = edit {
            total = total
                .checked_add(e.replacement.len())
                .filter(|n| *n <= 100_000)
                .ok_or(Error::TooLarge)?;
            edits.push(e);
        }
    }
    // Original coordinates, descending. At the same boundary remove/replace the
    // original base first, then insert before the resulting base/path.
    edits.sort_by(|a, b| b.start.cmp(&a.start).then_with(|| b.end.cmp(&a.end)));
    for pair in edits.windows(2) {
        let later = &pair[0];
        let earlier = &pair[1];
        if earlier.end > later.start
            && !(earlier.start == later.start && earlier.start == earlier.end)
        {
            return Err(bad("overlapping original-coordinate edits"));
        }
    }
    if edits.is_empty() {
        warnings.push(
            "Decisions were retained as a review record; no new target reference is created."
                .into(),
        );
    }
    warnings_guard(&warnings)?;
    Ok((edits, warnings))
}
fn apply(
    source: &Molecule,
    r: &Request,
    edits: &[Edit],
    stable: &str,
    reference: &EventBinding,
    reads: &[EventBinding],
) -> Result<Option<Molecule>, Error> {
    if edits.is_empty() {
        return Ok(None);
    }
    let mut target = source.clone();
    for e in edits {
        target = target.replace(e.start, e.end, &e.replacement)?;
    }
    if let Some(original) = target.metadata.remove("duplex_model") {
        let key = format!("read_review:source_duplex_model:{}", reference.revision);
        if target.metadata.get(&key).is_some_and(|v| v != &original) {
            return Err(bad("source descriptor evidence collision"));
        }
        target.metadata.insert(key, original);
    }
    target.name = r.name.clone();
    for (key, value) in [
        ("duplex_model", REFERENCE_MODEL.into()),
        ("read_review:reference_only", "true".into()),
        ("read_review:context_revision", stable.into()),
        ("read_review:check_hash", r.attachment_hash.clone()),
        (
            "read_review:source_length",
            source.sequence.len().to_string(),
        ),
        ("read_review:reference", json(reference)?),
        ("read_review:reads", json(&reads)?),
        ("read_review:edits", json(&edits)?),
        (
            "read_review:sequence_blake3",
            blake3::hash(target.sequence.as_str().as_bytes())
                .to_hex()
                .to_string(),
        ),
    ] {
        target.metadata.insert(key.into(), value);
    }
    collagene_model::bounds::model_preflight(&[&target])?;
    target.validate()?;
    Ok(Some(target))
}
fn fresh(
    w: &Workspace,
    r: &Request,
    domain: &str,
) -> Result<(Preview, TraceValidationReport, String), Error> {
    let context_revision = revision(w, r, domain)?;
    let checked = read_checks::replay(w, &r.attachment_hash)?;
    let reference = checked.saved.reference;
    let reads = checked.saved.reads;
    let source = read_checks::resolve(w, &reference)?;
    let (edits, warnings) = choose(source, r, &checked.report)?;
    let stable = stable_revision(r, &reference, &reads, &checked.report)?;
    let target = apply(source, r, &edits, &stable, &reference, &reads)?;
    let p = Preview {
        context_revision: context_revision.clone(),
        request: r.clone(),
        reference,
        reads,
        target,
        edits,
        warnings,
    };
    p.validate_binding(&context_revision)?;
    Ok((p, checked.report, stable))
}
pub fn preview(w: &Workspace, r: &Request) -> Result<Preview, Error> {
    Ok(fresh(w, r, "preview-read-review-1")?.0)
}
impl Preview {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        request_guard(&self.request)?;
        if !hex(expected)
            || self.context_revision != expected
            || self.edits.len() > self.request.decisions.len()
            || self.reads.is_empty()
            || self.reads.len() > 32
            || self.target.is_some() != !self.edits.is_empty()
        {
            return Err(bad("preview context or shape mismatch"));
        }
        warnings_guard(&self.warnings)?;
        let mut total = 0usize;
        for e in &self.edits {
            if e.start > e.end
                || e.end - e.start > 1
                || e.reason.len() > 1024
                || !e.replacement.bytes().all(|b| b"ACGT".contains(&b))
            {
                return Err(bad("invalid readonly edit"));
            }
            total = total
                .checked_add(e.replacement.len())
                .filter(|n| *n <= 100_000)
                .ok_or(Error::TooLarge)?;
        }
        if let Some(m) = &self.target {
            target_guard(
                m,
                &self.request,
                &self.reference,
                &self.reads,
                &self.edits,
                None,
            )?;
        }
        crate::encoding::bounded_json::hash(self)?;
        Ok(())
    }
}
fn load_inner(w: &Workspace, hash: &str) -> Result<SavedReview, Error> {
    if !hex(hash) {
        return Err(bad("invalid review hash"));
    }
    let a = w
        .attachments
        .get(hash)
        .ok_or_else(|| bad("missing review attachment"))?;
    if a.kind != ATTACHMENT_KIND
        || a.media_type != "application/json"
        || a.bytes.len() > MAX_BYTES
        || blake3::hash(&a.bytes).to_hex().as_str() != hash
    {
        return Err(bad("invalid review asset"));
    }
    let s: SavedReview =
        serde_json::from_slice(&a.bytes).map_err(|_| bad("malformed review record"))?;
    request_guard(&s.request)?;
    if s.schema_version != 1 || !hex(&s.context_revision) {
        return Err(bad("invalid review schema"));
    }
    read_checks::report_guard(&s.report, &s.reference, &s.reads)?;
    let check = read_checks::load(w, &s.request.attachment_hash)?;
    if check.reference != s.reference
        || check.reads != s.reads
        || check.report.options != s.report.options
    {
        return Err(bad("review does not bind its saved check"));
    }
    if stable_revision(&s.request, &s.reference, &s.reads, &s.report)? != s.context_revision {
        return Err(bad("review context changed"));
    }
    let source = read_checks::resolve(w, &s.reference)?;
    let (edits, _) = choose(source, &s.request, &s.report)?;
    if edits != s.edits || s.target.is_some() != !s.edits.is_empty() {
        return Err(bad("recorded choices disagree with their evidence"));
    }
    warnings_guard(&s.warnings)?;
    if let Some(binding) = &s.target {
        let m = read_checks::resolve(w, binding)?;
        if m.name != s.request.name
            || m.metadata.get("read_review:context_revision") != Some(&s.context_revision)
        {
            return Err(bad("recorded target identity changed"));
        }
        target_guard(
            m,
            &s.request,
            &s.reference,
            &s.reads,
            &s.edits,
            Some(&s.context_revision),
        )?;
        let event = w
            .history
            .iter()
            .find(|e| e.id == binding.event_id)
            .ok_or_else(|| bad("missing review target event"))?;
        if event.transition
            != (Transition::Analysis {
                origin: AnalysisOrigin::context(
                    AnalysisProducer::ReadReview,
                    s.context_revision.clone(),
                ),
            })
        {
            return Err(bad(
                "review target typed producer/context differs from its retained record",
            ));
        }
        let required = s
            .reads
            .iter()
            .chain(std::iter::once(&s.reference))
            .map(|b| b.event_id.as_str())
            .collect::<BTreeSet<_>>();
        if event
            .parents
            .iter()
            .map(String::as_str)
            .collect::<BTreeSet<_>>()
            != required
        {
            return Err(bad(
                "review target does not retain all source event parents",
            ));
        }
    }
    encode(&s)?;
    Ok(s)
}
pub fn load(w: &Workspace, hash: &str) -> Result<SavedReview, Error> {
    read_checks::workspace_guard(w)?;
    if w.attachments
        .values()
        .filter(|a| a.kind == ATTACHMENT_KIND)
        .count()
        > MAX_REVIEWS
    {
        return Err(Error::TooLarge);
    }
    load_inner(w, hash)
}
pub fn list(w: &Workspace) -> Result<Vec<Summary>, Error> {
    read_checks::workspace_guard(w)?;
    let hashes = w
        .attachments
        .iter()
        .filter(|(_, a)| a.kind == ATTACHMENT_KIND)
        .map(|(h, _)| h)
        .collect::<Vec<_>>();
    if hashes.len() > MAX_REVIEWS {
        return Err(Error::TooLarge);
    }
    hashes
        .into_iter()
        .map(|h| {
            let s = load_inner(w, h)?;
            Ok(Summary {
                attachment_hash: h.clone(),
                reference: s.reference,
                decision_count: s.request.decisions.len(),
                change_count: s.edits.len(),
                target: s.target,
            })
        })
        .collect()
}
pub fn commit(w: &Workspace, r: &Request) -> Result<Commit, Error> {
    let (p, report, stable) = fresh(w, r, "commit-read-review-1")?;
    let mut target = None;
    let mut existing_count = 0usize;
    for (hash, _) in w
        .attachments
        .iter()
        .filter(|(_, a)| a.kind == ATTACHMENT_KIND)
    {
        existing_count += 1;
        let s = load_inner(w, hash)?;
        if s.context_revision != stable {
            continue;
        }
        let reusable = match (&s.target, &p.target) {
            (None, None) => true,
            (Some(b), Some(m)) => {
                w.documents.iter().any(|d| {
                    d.id == b.document_id
                        && d.head_event.as_deref() == Some(b.event_id.as_str())
                        && d.head == b.revision
                }) && read_checks::resolve(w, b)? == m
            }
            _ => false,
        };
        if reusable {
            let c = Commit {
                workspace: w.clone(),
                context_revision: p.context_revision.clone(),
                attachment_hash: hash.clone(),
                target: s.target.clone(),
            };
            c.validate_binding(&p.context_revision)?;
            return Ok(c);
        }
    }
    let mut next = w.clone();
    if let Some(m) = &p.target {
        let histories = DocumentHistory::restore_all(&next)?;
        let bindings = std::iter::once(&p.reference)
            .chain(p.reads.iter())
            .collect::<Vec<_>>();
        let mut sources = Vec::with_capacity(bindings.len());
        for b in bindings {
            let h = histories
                .iter()
                .find(|h| h.id == b.document_id)
                .ok_or_else(|| bad("source history missing"))?;
            sources.push((h, b.event_id.as_str()));
        }
        let base = format!("read-review-{stable}");
        let mut id = base.clone();
        let mut suffix = 1usize;
        while next.documents.iter().any(|d| d.id == id) {
            id = format!("{base}-{suffix}");
            suffix += 1;
        }
        let derived = DocumentHistory::derive_from_events(
            &sources,
            id.clone(),
            m.clone(),
            "Explicit reviewed reference from original construct and all read events",
            AnalysisOrigin::context(AnalysisProducer::ReadReview, stable.clone()),
        )?;
        let binding = EventBinding {
            document_id: id,
            event_id: derived.head_event().into(),
            revision: derived.head_revision().into(),
        };
        derived.contribute(&mut next)?;
        target = Some(binding);
    }
    let saved = SavedReview {
        schema_version: 1,
        context_revision: stable,
        request: r.clone(),
        reference: p.reference,
        reads: p.reads,
        report,
        edits: p.edits,
        warnings: p.warnings,
        target: target.clone(),
    };
    let bytes = encode(&saved)?;
    let hash = blake3::hash(&bytes).to_hex().to_string();
    if !next.attachments.contains_key(&hash) && existing_count >= MAX_REVIEWS {
        return Err(Error::TooLarge);
    }
    if next
        .attachments
        .get(&hash)
        .is_some_and(|a| a.kind != ATTACHMENT_KIND || a.bytes != bytes)
    {
        return Err(bad("review attachment identity collision"));
    }
    next.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: ATTACHMENT_KIND.into(),
            name: format!("Read review {}", &saved.context_revision[..12]),
            media_type: "application/json".into(),
            bytes,
        },
    );
    let c = Commit {
        workspace: next,
        context_revision: p.context_revision.clone(),
        attachment_hash: hash,
        target,
    };
    c.validate_binding(&p.context_revision)?;
    Ok(c)
}
impl Commit {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !hex(expected) || self.context_revision != expected {
            return Err(bad("commit request changed"));
        }
        let s = load(&self.workspace, &self.attachment_hash)?;
        if self.target != s.target {
            return Err(bad("commit target disagrees with historical review"));
        }
        if let Some(b) = &self.target {
            if !self.workspace.documents.iter().any(|d| {
                d.id == b.document_id
                    && d.head_event.as_deref() == Some(b.event_id.as_str())
                    && d.head == b.revision
            }) {
                return Err(bad("reviewed target is not its current document head"));
            }
        }
        Ok(())
    }
}

fn target_guard(
    m: &Molecule,
    r: &Request,
    reference: &EventBinding,
    reads: &[EventBinding],
    edits: &[Edit],
    stable: Option<&str>,
) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[m])?;
    m.validate()?;
    let n = m
        .metadata
        .get("read_review:source_length")
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|n| *n > 0 && *n <= 8 << 20)
        .ok_or_else(|| bad("missing bounded original reference length"))?;
    let mut removed = 0usize;
    let mut added = 0usize;
    for e in edits {
        if e.start > e.end || e.end > n || e.end - e.start > 1 {
            return Err(bad("review edits exceed their source coordinates"));
        }
        removed = removed
            .checked_add(e.end - e.start)
            .ok_or(Error::TooLarge)?;
        added = added
            .checked_add(e.replacement.len())
            .ok_or(Error::TooLarge)?;
    }
    if n.checked_sub(removed).and_then(|n| n.checked_add(added)) != Some(m.sequence.len()) {
        return Err(bad("reviewed target length disagrees with edits"));
    }
    if m.name != r.name
        || m.metadata.get("duplex_model").map(String::as_str) != Some(REFERENCE_MODEL)
        || m.metadata
            .get("read_review:reference_only")
            .map(String::as_str)
            != Some("true")
        || m.metadata.get("read_review:check_hash") != Some(&r.attachment_hash)
        || m.metadata.get("read_review:reference") != Some(&json(reference)?)
        || m.metadata.get("read_review:reads") != Some(&json(&reads)?)
        || m.metadata.get("read_review:edits") != Some(&json(&edits)?)
        || m.metadata
            .get("read_review:sequence_blake3")
            .map(String::as_str)
            != Some(
                blake3::hash(m.sequence.as_str().as_bytes())
                    .to_hex()
                    .as_str(),
            )
        || !m
            .metadata
            .get("read_review:context_revision")
            .is_some_and(|s| hex(s) && stable.is_none_or(|v| v == s))
    {
        return Err(bad("review target evidence mismatch"));
    }
    Ok(())
}
