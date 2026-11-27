//! Documentary well assignments to immutable native construct declarations.
//! No scientific approval, quantities, inventory state or instrument instructions.
use crate::{
    Error, combinatorial_batch as batch, combinatorial_release as release,
    run_lineage::{self, EventBinding},
    verification::Status,
    workspace::{Workspace, WorkspaceAttachment},
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::{self, Write},
};
pub const KIND: &str = "protocol-layout-v1";
pub const ENGINE_VERSION: &str = "protocol-layout-1";
pub const MAX_WELLS: usize = 384;
pub const MAX_LAYOUTS: usize = 128;
pub const MAX_BYTES: usize = 4 * 1024 * 1024;
pub const MAX_REFERENCED_BYTES: usize = 32 * 1024 * 1024;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Plate {
    Wells96,
    Wells384,
}
impl Plate {
    pub fn dimensions(self) -> (u8, u8) {
        match self {
            Self::Wells96 => (8, 12),
            Self::Wells384 => (16, 24),
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Well {
    pub row: u8,
    pub column: u8,
}
impl Well {
    pub fn label(self) -> String {
        format!(
            "{}{}",
            char::from_u32(u32::from(b'A') + u32::from(self.row)).unwrap_or('?'),
            u16::from(self.column) + 1
        )
    }
}
dto! {pub enum Source{Planned{batch_hash:String,construct_id:String,declaration_revision:String},Released{manifest_hash:String,construct_id:String,run_id:String,product:EventBinding}}}
dto! {pub struct Assignment{pub id:String,pub well:Well,pub label:String,pub source:Source}}
dto! {pub struct Request{pub name:String,pub plate:Plate,#[serde(deserialize_with="assignments")]pub assignments:Vec<Assignment>}}
fn assignments<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<Assignment>, D::Error> {
    struct V;
    impl<'de> serde::de::Visitor<'de> for V {
        type Value = Vec<Assignment>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "at most384 explicit well assignments")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            if seq.size_hint().is_some_and(|n| n > MAX_WELLS) {
                return Err(serde::de::Error::custom("well bound"));
            }
            let mut out = vec![];
            while out.len() < MAX_WELLS {
                match seq.next_element()? {
                    Some(x) => out.push(x),
                    None => return Ok(out),
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom("well bound"));
            }
            Ok(out)
        }
    }
    d.deserialize_seq(V)
}
dto! {pub struct ResolvedAssignment{pub assignment_id:String,pub name:String,pub recorded_engine:String,pub recorded_status:Option<Status>,pub declaration_revision:String,pub target:Option<EventBinding>,pub run_id:Option<String>,pub product:Option<EventBinding>,pub missing_quantities:bool}}
dto! {pub struct Layout{pub schema_version:u16,pub engine_version:String,pub source_workspace_fingerprint:String,pub context_revision:String,pub request:Request,#[serde(deserialize_with="resolved_rows")]pub resolved:Vec<ResolvedAssignment>,#[serde(deserialize_with="warning_rows")]pub warnings:Vec<String>}}
dto! {pub struct Capture{#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub workspace:Workspace,pub context_revision:String,pub source_workspace_fingerprint:String,pub attachment_hash:String,pub layout:Layout}}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ExportFormat {
    Csv,
    Html,
    Pdf,
}
dto! {pub struct Exported{pub context_revision:String,pub layout_hash:String,pub format:ExportFormat,pub filename:String,pub media_type:String,#[serde(deserialize_with="export_bytes")]pub bytes:Vec<u8>,pub content_hash:String,#[serde(deserialize_with="warning_rows")]pub warnings:Vec<String>}}
dto! {pub struct Summary{pub attachment_hash:String,pub name:String,pub plate:Plate,pub assignment_count:usize,pub context_revision:String}}
#[derive(Debug, Clone)]
pub struct Dependencies {
    pub attachments: Vec<String>,
    pub runs: Vec<String>,
    pub events: Vec<EventBinding>,
}
fn bounded_rows<'de, D, T, const N: usize>(d: D) -> Result<Vec<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    struct V<T, const N: usize>(std::marker::PhantomData<T>);
    impl<'de, T: Deserialize<'de>, const N: usize> serde::de::Visitor<'de> for V<T, N> {
        type Value = Vec<T>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "at most{N} bounded rows")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            if seq.size_hint().is_some_and(|n| n > N) {
                return Err(serde::de::Error::custom("row bound"));
            }
            let mut out = vec![];
            while out.len() < N {
                match seq.next_element()? {
                    Some(x) => out.push(x),
                    None => return Ok(out),
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom("row bound"));
            }
            Ok(out)
        }
    }
    d.deserialize_seq(V::<T, N>(std::marker::PhantomData))
}
fn resolved_rows<'de, D: serde::Deserializer<'de>>(
    d: D,
) -> Result<Vec<ResolvedAssignment>, D::Error> {
    bounded_rows::<D, ResolvedAssignment, MAX_WELLS>(d)
}
fn warning_rows<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<String>, D::Error> {
    bounded_rows::<D, String, 16>(d)
}
fn export_bytes<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<u8>, D::Error> {
    bounded_rows::<D, u8, MAX_BYTES>(d)
}
fn bad(s: &str) -> Error {
    Error::Format(format!("Protocol layout: {s}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn text(s: &str, n: usize, empty: bool) -> bool {
    (empty || !s.trim().is_empty())
        && s.chars().count() <= n
        && s.len() <= n * 4
        && !s.chars().any(char::is_control)
}
fn id(s: &str, n: usize) -> bool {
    !s.is_empty()
        && s.len() <= n
        && s.bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_-.:".contains(&b))
}
fn event(e: &EventBinding) -> Result<(), Error> {
    if !id(&e.document_id, 128) || !id(&e.event_id, 128) || !digest(&e.revision) {
        return Err(bad("invalid exact event pin"));
    }
    Ok(())
}
struct Bytes(Vec<u8>);
impl Write for Bytes {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        if b.len() > MAX_BYTES.saturating_sub(self.0.len()) {
            return Err(io::Error::other("layout byte bound"));
        }
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn encode(v: &impl Serialize) -> Result<Vec<u8>, Error> {
    let mut b = Bytes(vec![]);
    serde_json::to_writer(&mut b, v).map_err(|_| Error::TooLarge)?;
    Ok(b.0)
}
fn hash(v: &impl Serialize) -> Result<String, Error> {
    Ok(blake3::hash(&encode(v)?).to_hex().to_string())
}
fn warnings() -> Vec<String> {
    vec!["Documentary well assignment only. Recorded construct status is historical evidence; assignment does not approve a candidate, establish physical material, experimental success or execution.".into(),"Quantities, concentrations, reaction volumes and reagent recipes are missing for every assignment. No zero, standard recipe, well capacity or recommended amount was inferred.".into()]
}
impl Request {
    pub fn validate(&self) -> Result<(), Error> {
        if !text(&self.name, 64, false) || self.assignments.len() > MAX_WELLS {
            return Err(bad("invalid layout label or assignment bound"));
        }
        let (rows, cols) = self.plate.dimensions();
        let mut ids = BTreeSet::new();
        let mut wells = BTreeSet::new();
        for a in &self.assignments {
            if !id(&a.id, 64)
                || !ids.insert(&a.id)
                || !wells.insert(a.well)
                || a.well.row >= rows
                || a.well.column >= cols
                || !text(&a.label, 64, true)
            {
                return Err(bad(
                    "duplicate assignment/well, out-of-range coordinate or invalid label",
                ));
            }
            match &a.source {
                Source::Planned {
                    batch_hash,
                    construct_id,
                    declaration_revision,
                } => {
                    if !digest(batch_hash) || !id(construct_id, 48) || !digest(declaration_revision)
                    {
                        return Err(bad("invalid planned identity"));
                    }
                }
                Source::Released {
                    manifest_hash,
                    construct_id,
                    run_id,
                    product,
                } => {
                    if !digest(manifest_hash) || !id(construct_id, 48) || !digest(run_id) {
                        return Err(bad("invalid released identity"));
                    }
                    event(product)?
                }
            }
        }
        encode(self)?;
        Ok(())
    }
}
fn layouts_bound(w: &Workspace) -> Result<(), Error> {
    let mut count = 0;
    for a in w.attachments.values().filter(|a| a.kind == KIND) {
        count += 1;
        if count > MAX_LAYOUTS || a.bytes.len() > MAX_BYTES {
            return Err(Error::TooLarge);
        }
    }
    Ok(())
}
fn context(fp: &str, r: &Request) -> Result<String, Error> {
    hash(&(ENGINE_VERSION, fp, r))
}
pub fn preview_revision(w: &Workspace, r: &Request) -> Result<String, Error> {
    r.validate()?;
    layouts_bound(w)?;
    context(&run_lineage::workspace_fingerprint(w)?, r)
}
struct Cache<'a> {
    w: &'a Workspace,
    charged: BTreeSet<String>,
    bytes: usize,
    batches: BTreeMap<String, batch::Record>,
    releases: BTreeMap<String, release::Manifest>,
}
impl<'a> Cache<'a> {
    fn new(w: &'a Workspace) -> Self {
        Self {
            w,
            charged: BTreeSet::new(),
            bytes: 0,
            batches: BTreeMap::new(),
            releases: BTreeMap::new(),
        }
    }
    fn charge(&mut self, hash: &str, kind: &str) -> Result<(), Error> {
        let a = self
            .w
            .attachments
            .get(hash)
            .ok_or_else(|| bad("referenced immutable attachment missing"))?;
        if !digest(hash)
            || a.kind != kind
            || a.media_type != "application/json"
            || a.bytes.len() > MAX_REFERENCED_BYTES
        {
            return Err(bad("referenced attachment identity/type/hash differs"));
        }
        if self.charged.insert(hash.into()) {
            self.bytes = self
                .bytes
                .checked_add(a.bytes.len())
                .filter(|n| *n <= MAX_REFERENCED_BYTES)
                .ok_or(Error::TooLarge)?
        }
        if blake3::hash(&a.bytes).to_hex().as_str() != hash {
            return Err(bad("referenced attachment hash differs"));
        }
        Ok(())
    }
    fn batch(&mut self, h: &str) -> Result<&batch::Record, Error> {
        if !self.batches.contains_key(h) {
            self.charge(h, batch::KIND)?;
            let r = batch::load(self.w, h)?;
            self.batches.insert(h.into(), r);
        }
        Ok(&self.batches[h])
    }
    fn release(&mut self, h: &str) -> Result<&release::Manifest, Error> {
        if !self.releases.contains_key(h) {
            self.charge(h, release::KIND)?;
            #[derive(Deserialize)]
            struct Header {
                batch_hash: String,
            }
            let header: Header = serde_json::from_slice(&self.w.attachments[h].bytes)
                .map_err(|e| bad(&e.to_string()))?;
            self.batch(&header.batch_hash)?;
            let r = release::load(self.w, h)?;
            self.releases.insert(h.into(), r);
        }
        Ok(&self.releases[h])
    }
}
fn resolve<'a>(
    w: &'a Workspace,
    r: &Request,
) -> Result<(Vec<ResolvedAssignment>, Cache<'a>), Error> {
    let mut cache = Cache::new(w);
    let mut out = Vec::with_capacity(r.assignments.len());
    for a in &r.assignments {
        let row = match &a.source {
            Source::Planned {
                batch_hash,
                construct_id,
                declaration_revision,
            } => {
                let b = cache.batch(batch_hash)?;
                let mut found = b.rows.iter().filter(|r| &r.construct_id == construct_id);
                let row = found
                    .next()
                    .ok_or_else(|| bad("planned construct ID missing"))?;
                if found.next().is_some()
                    || row.declaration_revision.as_ref() != Some(declaration_revision)
                {
                    return Err(bad("exact planned declaration pin changed or unresolved"));
                }
                ResolvedAssignment {
                    assignment_id: a.id.clone(),
                    name: row.name.clone(),
                    recorded_engine: b.engine_version.clone(),
                    recorded_status: row.status,
                    declaration_revision: declaration_revision.clone(),
                    target: row.target_binding.clone(),
                    run_id: None,
                    product: None,
                    missing_quantities: true,
                }
            }
            Source::Released {
                manifest_hash,
                construct_id,
                run_id,
                product,
            } => {
                let m = cache.release(manifest_hash)?;
                let mut found = m.rows.iter().filter(|r| &r.construct_id == construct_id);
                let row = found
                    .next()
                    .ok_or_else(|| bad("released construct ID missing"))?;
                if found.next().is_some()
                    || row.output.run.run_id != *run_id
                    || row.output.construct != *product
                {
                    return Err(bad("released run/product event pin differs"));
                }
                ResolvedAssignment {
                    assignment_id: a.id.clone(),
                    name: row.options.name.clone(),
                    recorded_engine: m.engine_version.clone(),
                    recorded_status: Some(row.output.run.recorded_status),
                    declaration_revision: row.declaration_revision.clone(),
                    target: Some(row.target.clone()),
                    run_id: Some(run_id.clone()),
                    product: Some(product.clone()),
                    missing_quantities: true,
                }
            }
        };
        out.push(row)
    }
    Ok((out, cache))
}
pub fn preview(w: &Workspace, r: &Request) -> Result<Layout, Error> {
    let revision = preview_revision(w, r)?;
    let (resolved, _) = resolve(w, r)?;
    let result = Layout {
        schema_version: 1,
        engine_version: ENGINE_VERSION.into(),
        source_workspace_fingerprint: run_lineage::workspace_fingerprint(w)?,
        context_revision: revision.clone(),
        request: r.clone(),
        resolved,
        warnings: warnings(),
    };
    result.validate_binding(&revision)?;
    Ok(result)
}
impl Layout {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.request.validate()?;
        if self.schema_version != 1
            || self.engine_version != ENGINE_VERSION
            || !digest(&self.source_workspace_fingerprint)
            || !digest(expected)
            || self.context_revision != expected
            || context(&self.source_workspace_fingerprint, &self.request)? != expected
            || self.resolved.len() != self.request.assignments.len()
            || self.warnings != warnings()
        {
            return Err(bad("layout request/context/shape changed"));
        }
        for (a, r) in self.request.assignments.iter().zip(&self.resolved) {
            if r.assignment_id != a.id
                || !text(&r.name, 256, false)
                || !text(&r.recorded_engine, 256, false)
                || !digest(&r.declaration_revision)
                || !r.missing_quantities
            {
                return Err(bad("resolved documentary row changed"));
            }
            if let Some(e) = &r.target {
                event(e)?
            }
            match &a.source {
                Source::Planned {
                    declaration_revision,
                    ..
                } => {
                    if r.declaration_revision != *declaration_revision
                        || r.run_id.is_some()
                        || r.product.is_some()
                    {
                        return Err(bad("planned row acquired released authority"));
                    }
                }
                Source::Released {
                    run_id, product, ..
                } => {
                    if r.run_id.as_ref() != Some(run_id)
                        || r.product.as_ref() != Some(product)
                        || r.target.is_none()
                        || r.recorded_status.is_none()
                    {
                        return Err(bad("released row pin changed"));
                    }
                }
            }
        }
        encode(self)?;
        Ok(())
    }
    pub fn validate_for_workspace(&self, w: &Workspace) -> Result<(), Error> {
        self.validate_binding(&self.context_revision)?;
        layouts_bound(w)?;
        w.validate()?;
        let (actual, _) = resolve(w, &self.request)?;
        if encode(&actual)? != encode(&self.resolved)? {
            return Err(bad(
                "layout source bindings no longer match immutable native archives",
            ));
        }
        Ok(())
    }
}
pub fn load_record(bytes: &[u8]) -> Result<Layout, Error> {
    if bytes.len() > MAX_BYTES {
        return Err(Error::TooLarge);
    }
    let r: Layout = serde_json::from_slice(bytes).map_err(|e| bad(&e.to_string()))?;
    r.validate_binding(&r.context_revision)?;
    Ok(r)
}
pub fn load(w: &Workspace, h: &str) -> Result<Layout, Error> {
    layouts_bound(w)?;
    let a = w
        .attachments
        .get(h)
        .ok_or_else(|| bad("layout attachment missing"))?;
    if !digest(h)
        || a.kind != KIND
        || a.media_type != "application/json"
        || a.bytes.len() > MAX_BYTES
        || blake3::hash(&a.bytes).to_hex().as_str() != h
    {
        return Err(bad("layout hash/type differs"));
    }
    let r = load_record(&a.bytes)?;
    r.validate_for_workspace(w)?;
    Ok(r)
}
pub fn list(w: &Workspace) -> Result<Vec<Summary>, Error> {
    layouts_bound(w)?;
    w.validate()?;
    w.attachments
        .iter()
        .filter(|(_, a)| a.kind == KIND)
        .map(|(h, a)| {
            let r = load_record(&a.bytes)?;
            Ok(Summary {
                attachment_hash: h.clone(),
                name: r.request.name,
                plate: r.request.plate,
                assignment_count: r.request.assignments.len(),
                context_revision: r.context_revision,
            })
        })
        .collect()
}
pub fn capture(w: &Workspace, r: &Request) -> Result<Capture, Error> {
    if w.attachments.values().filter(|a| a.kind == KIND).count() >= MAX_LAYOUTS {
        return Err(Error::TooLarge);
    }
    let layout = preview(w, r)?;
    let bytes = encode(&layout)?;
    let attachment_hash = blake3::hash(&bytes).to_hex().to_string();
    let mut workspace = w.clone();
    workspace.attachments.insert(
        attachment_hash.clone(),
        WorkspaceAttachment {
            kind: KIND.into(),
            name: format!("layout-{}.json", &attachment_hash[..16]),
            media_type: "application/json".into(),
            bytes,
        },
    );
    layouts_bound(&workspace)?;
    workspace.validate()?;
    let result = Capture {
        context_revision: layout.context_revision.clone(),
        source_workspace_fingerprint: layout.source_workspace_fingerprint.clone(),
        layout,
        attachment_hash,
        workspace,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl Capture {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.layout.validate_binding(expected)?;
        if self.context_revision != expected
            || self.source_workspace_fingerprint != self.layout.source_workspace_fingerprint
            || !digest(&self.attachment_hash)
        {
            return Err(bad("capture context changed"));
        }
        let a = self
            .workspace
            .attachments
            .get(&self.attachment_hash)
            .ok_or_else(|| bad("capture attachment missing"))?;
        if a.kind != KIND
            || a.media_type != "application/json"
            || a.bytes != encode(&self.layout)?
            || blake3::hash(&a.bytes).to_hex().as_str() != self.attachment_hash
        {
            return Err(bad("capture attachment differs"));
        }
        self.layout.validate_for_workspace(&self.workspace)?;
        let mut original = self.workspace.clone();
        original.attachments.remove(&self.attachment_hash);
        if run_lineage::workspace_fingerprint(&original)? != self.source_workspace_fingerprint {
            return Err(bad(
                "capture changed source workspace beyond its single layout attachment",
            ));
        }
        Ok(())
    }
}
/// Full actual dependencies of indivisible included records, not only displayed
/// assigned rows. Explicit layout publication therefore includes these whole
/// batch/release archives. Ordinary plasmid publication must omit this sidecar.
pub fn source_bindings(w: &Workspace, layout: &Layout) -> Result<Dependencies, Error> {
    layout.validate_binding(&layout.context_revision)?;
    let (rows, cache) = resolve(w, &layout.request)?;
    if encode(&rows)? != encode(&layout.resolved)? {
        return Err(bad("layout rows differ from archives"));
    }
    let mut attachments = cache.charged.clone();
    let mut runs = BTreeSet::new();
    let mut events = BTreeMap::new();
    let mut add = |e: &EventBinding| {
        events.insert(
            (
                e.document_id.clone(),
                e.event_id.clone(),
                e.revision.clone(),
            ),
            e.clone(),
        );
    };
    for b in cache.batches.values() {
        if let Some(d) = &b.request.design {
            for bin in &d.bins {
                for choice in &bin.choices {
                    add(&choice.source)
                }
            }
        }
        for r in &b.request.rows {
            if let batch::RowRequest::Ready(c) = r {
                if let crate::combinatorial_design::TargetDeclaration::Existing(e) = &c.target {
                    add(e)
                }
            }
        }
        for r in &b.rows {
            for i in &r.input_bindings {
                add(&i.event)
            }
            if let Some(e) = &r.target_binding {
                add(e)
            }
        }
    }
    for m in cache.releases.values() {
        attachments.insert(m.batch_hash.clone());
        for row in &m.rows {
            runs.insert(row.output.run.run_id.clone());
            attachments.insert(row.output.run.attachment_hash.clone());
            if let Some(h) = &row.output.primer_attachment {
                attachments.insert(h.clone());
            }
            for i in &row.inputs {
                add(&i.event)
            }
            add(&row.target);
            add(&row.output.construct);
            for p in &row.output.run.products {
                add(&EventBinding {
                    document_id: p.document_id.clone(),
                    event_id: p.event_id.clone(),
                    revision: p.revision.clone(),
                });
            }
        }
    }
    Ok(Dependencies {
        attachments: attachments.into_iter().collect(),
        runs: runs.into_iter().collect(),
        events: events.into_values().collect(),
    })
}
fn export_context(h: &str, f: ExportFormat) -> Result<String, Error> {
    hash(&(ENGINE_VERSION, "export", h, f))
}
pub fn export_revision(w: &Workspace, h: &str, f: ExportFormat) -> Result<String, Error> {
    load(w, h)?;
    export_context(h, f)
}
fn identity(source: &Source) -> String {
    match source {
        Source::Planned {
            batch_hash,
            construct_id,
            declaration_revision,
        } => format!(
            "Planned batch {batch_hash}; construct {construct_id}; declaration {declaration_revision}"
        ),
        Source::Released {
            manifest_hash,
            construct_id,
            run_id,
            product,
        } => format!(
            "Released manifest {manifest_hash}; construct {construct_id}; run {run_id}; document {}; event {}; revision {}",
            product.document_id, product.event_id, product.revision
        ),
    }
}
fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}
fn csv(s: &str) -> String {
    format!("\"{}\"", s.replace('"', "\"\""))
}
fn profile(f: ExportFormat) -> (&'static str, &'static str) {
    match f {
        ExportFormat::Csv => ("csv", "text/csv;charset=utf-8"),
        ExportFormat::Html => ("html", "text/html;charset=utf-8"),
        ExportFormat::Pdf => ("pdf", "application/pdf"),
    }
}
pub fn export(w: &Workspace, h: &str, f: ExportFormat) -> Result<Exported, Error> {
    let l = load(w, h)?;
    let mut order = (0..l.request.assignments.len()).collect::<Vec<_>>();
    order.sort_by_key(|i| l.request.assignments[*i].well);
    let mut extra = vec![];
    let bytes = match f {
        ExportFormat::Csv => {
            let mut out = Bytes(vec![]);
            writeln!(
                out,
                "{}",
                [
                    "plate",
                    "well",
                    "row_index",
                    "column_index",
                    "assignment_id",
                    "label",
                    "source_json",
                    "source_name",
                    "declaration_revision",
                    "target_json",
                    "recorded_engine",
                    "recorded_status",
                    "missing_quantities",
                    "layout_hash",
                    "layout_context"
                ]
                .map(csv)
                .join(",")
            )
            .map_err(|_| Error::TooLarge)?;
            for i in &order {
                let a = &l.request.assignments[*i];
                let r = &l.resolved[*i];
                let values = [
                    format!("{:?}", l.request.plate),
                    a.well.label(),
                    a.well.row.to_string(),
                    a.well.column.to_string(),
                    a.id.clone(),
                    a.label.clone(),
                    serde_json::to_string(&a.source).map_err(|_| Error::TooLarge)?,
                    r.name.clone(),
                    r.declaration_revision.clone(),
                    serde_json::to_string(&r.target).map_err(|_| Error::TooLarge)?,
                    r.recorded_engine.clone(),
                    r.recorded_status
                        .map(|s| format!("{s:?}"))
                        .unwrap_or_else(|| "Not evaluated".into()),
                    "Missing; no quantities inferred".into(),
                    h.into(),
                    l.context_revision.clone(),
                ];
                writeln!(
                    out,
                    "{}",
                    values.iter().map(|x| csv(x)).collect::<Vec<_>>().join(",")
                )
                .map_err(|_| Error::TooLarge)?;
            }
            out.0
        }
        ExportFormat::Html => {
            let mut out = Bytes(vec![]);
            write!(out,"<!doctype html><html><meta charset=\"utf-8\"><title>{}</title><body><h1>{}</h1><p>Documentary plate layout. Layout hash <code>{h}</code>; context <code>{}</code>.</p>",escape(&l.request.name),escape(&l.request.name),l.context_revision).map_err(|_|Error::TooLarge)?;
            for warning in &l.warnings {
                write!(out, "<p>{}</p>", escape(warning)).map_err(|_| Error::TooLarge)?
            }
            let (rows, cols) = l.request.plate.dimensions();
            write!(
                out,
                "<table border=\"1\"><caption>{:?}</caption>",
                l.request.plate
            )
            .map_err(|_| Error::TooLarge)?;
            for row in 0..rows {
                write!(out, "<tr>").map_err(|_| Error::TooLarge)?;
                for column in 0..cols {
                    let well = Well { row, column };
                    let found = l.request.assignments.iter().position(|a| a.well == well);
                    write!(out, "<td><strong>{}</strong>", well.label())
                        .map_err(|_| Error::TooLarge)?;
                    if let Some(i) = found {
                        write!(
                            out,
                            "<br>{}<br>{}",
                            escape(&l.request.assignments[i].label),
                            escape(&l.resolved[i].name)
                        )
                        .map_err(|_| Error::TooLarge)?
                    }
                    write!(out, "</td>").map_err(|_| Error::TooLarge)?
                }
                write!(out, "</tr>").map_err(|_| Error::TooLarge)?
            }
            write!(out, "</table><ol>").map_err(|_| Error::TooLarge)?;
            for i in &order {
                let a = &l.request.assignments[*i];
                let r = &l.resolved[*i];
                write!(out,"<li>{} / {} / {}: {}. <code>{}</code>. Declaration {}; target <code>{}</code>; recorded engine {}; historical status {:?}; quantities missing.</li>",a.well.label(),escape(&a.id),escape(&a.label),escape(&r.name),escape(&identity(&a.source)),r.declaration_revision,escape(&serde_json::to_string(&r.target).map_err(|_|Error::TooLarge)?),escape(&r.recorded_engine),r.recorded_status).map_err(|_|Error::TooLarge)?
            }
            write!(out, "</ol></body></html>").map_err(|_| Error::TooLarge)?;
            out.0
        }
        ExportFormat::Pdf => {
            use crate::pdf::{DocumentaryKind, ReviewBlock};
            let mut blocks = vec![
                ReviewBlock::Heading(l.request.name.clone()),
                ReviewBlock::Paragraph(format!(
                    "Layout hash {h}; context {}. {:?}.",
                    l.context_revision, l.request.plate
                )),
            ];
            blocks.extend(l.warnings.iter().cloned().map(ReviewBlock::Paragraph));
            for i in &order {
                let a = &l.request.assignments[*i];
                let r = &l.resolved[*i];
                blocks.push(ReviewBlock::Heading(format!(
                    "{} / {} / {}",
                    a.well.label(),
                    a.id,
                    a.label
                )));
                blocks.push(ReviewBlock::Paragraph(format!("{}; {}. Declaration {}; target {}. Recorded engine {}; historical status {:?}. Quantities/concentrations/recipe missing.",r.name,identity(&a.source),r.declaration_revision,serde_json::to_string(&r.target).map_err(|_|Error::TooLarge)?,r.recorded_engine,r.recorded_status)));
            }
            let (bytes, _, warnings) = crate::pdf::export_documentary_pdf(
                DocumentaryKind::PlateLayout,
                &l.request.name,
                "Documentary well assignments; quantities missing",
                h,
                &l.context_revision,
                &blocks,
            )?;
            extra = warnings;
            bytes
        }
    };
    if bytes.len() > MAX_BYTES {
        return Err(Error::TooLarge);
    }
    let (ext, media) = profile(f);
    let result = Exported {
        context_revision: export_context(h, f)?,
        layout_hash: h.into(),
        format: f,
        filename: format!("protocol-layout-{}.{ext}", &h[..16]),
        media_type: media.into(),
        content_hash: blake3::hash(&bytes).to_hex().to_string(),
        bytes,
        warnings: extra,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl Exported {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let (ext, media) = profile(self.format);
        if !digest(&self.layout_hash)
            || !digest(&self.content_hash)
            || self.context_revision != expected
            || export_context(&self.layout_hash, self.format)? != expected
            || self.filename != format!("protocol-layout-{}.{ext}", &self.layout_hash[..16])
            || self.media_type != media
            || self.bytes.is_empty()
            || self.bytes.len() > MAX_BYTES
            || blake3::hash(&self.bytes).to_hex().as_str() != self.content_hash
            || self.warnings.len() > 16
            || self.warnings.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("export receipt changed"));
        }
        match self.format {
            ExportFormat::Pdf if !self.bytes.starts_with(b"%PDF-") => {
                return Err(bad("invalid native PDF"));
            }
            ExportFormat::Html if !self.bytes.starts_with(b"<!doctype html>") => {
                return Err(bad("invalid HTML"));
            }
            ExportFormat::Csv if !self.bytes.starts_with(b"\"plate\",\"well\"") => {
                return Err(bad("invalid CSV"));
            }
            _ => {}
        }
        Ok(())
    }
}
