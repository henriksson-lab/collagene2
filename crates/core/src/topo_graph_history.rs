//! Exact native elementary TOPO graphs. Documentary records never restore a
//! supplied ProductGraph or ordinary duplex; actual calls/transfers replay fresh.
use crate::{
    Error, assembly_recipe,
    run_lineage::{self, EventBinding},
    run_material::{self, MaterialOrigin},
    topo, verification,
    workspace::{Workspace, WorkspaceAttachment},
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const ENGINE_VERSION: &str = "topo-ordered-graph-2";
pub const WORKSPACE_ATTACHMENT_KIND: &str = "topo-graph-v2";
pub const MEDIA_TYPE: &str = "application/vnd.collagene.topo-graph";
pub const MAX_RECORD_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_CAPSULE_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_RECORDS: usize = 128;
const MAGIC: &[u8; 8] = b"TGRF\0\0\0\x02";
macro_rules! dto{($($t:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)]pub struct Selection{pub record_hash:String,pub state_revision:String}}
dto! {pub struct SourceBinding{pub role:String,pub event:EventBinding,pub origin:Option<MaterialOrigin>}}
dto! {pub struct CaptureRequest{pub name:String,pub bindings:Vec<SourceBinding>,pub request:topo::Request}}
dto! {pub struct Record{pub schema_version:u16,pub engine_version:String,pub scientific_engine:String,pub source_workspace_fingerprint:String,pub context_revision:String,pub request:CaptureRequest,pub state_revision:String,pub incomplete:bool,pub state:StateSummary}}
dto! {pub struct Capsule{pub schema_version:u16,pub selection:Selection,#[serde(with="opaque_bytes")]#[cfg_attr(feature="json-schema",schemars(with="String"))]pub record_bytes:Vec<u8>,pub source_snapshots:Vec<topo::Input>}}
dto! {pub struct StateSummary{pub source_bases:usize,pub chains:usize,pub components:usize,pub remaining_activations:Vec<String>,pub released_activations:Vec<String>,pub findings:Vec<String>,pub assumptions:Vec<String>}}
dto! {pub struct Summary{pub selection:Selection,pub name:String,pub source_events:Vec<EventBinding>,pub source_bases:usize,pub chains:usize,pub components:usize,pub remaining_activations:Vec<String>,pub released_activations:Vec<String>,pub incomplete:bool,pub findings:Vec<String>,pub assumptions:Vec<String>}}
dto! {pub struct Receipt{pub context_revision:String,pub source_workspace_fingerprint:String,pub record_hash:String,pub selection:Selection,pub summary:Summary}}
dto! {pub struct PreparedCapture{#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub workspace:Workspace,pub receipt:Receipt}}
dto! {pub struct Listed{pub context_revision:String,pub source_workspace_fingerprint:String,pub entries:Vec<Summary>}}
dto! {pub struct Inspection{pub context_revision:String,pub source_workspace_fingerprint:String,pub selection:Selection,pub summary:Summary,pub record:Record,pub report:topo::TopoReport}}
/// Constructed only by a current native source/transcript execution.
pub struct FreshGraph {
    pub(crate) report: topo::TopoReport,
    pub(crate) incomplete: bool,
}
impl FreshGraph {
    pub fn report(&self) -> &topo::TopoReport {
        &self.report
    }
    pub fn incomplete(&self) -> bool {
        self.incomplete
    }
}
mod opaque_bytes {
    use super::*;
    use base64::{Engine, engine::general_purpose::STANDARD};
    pub fn serialize<S: serde::Serializer>(b: &[u8], s: S) -> Result<S::Ok, S::Error> {
        if b.len() > MAX_RECORD_BYTES {
            return Err(serde::ser::Error::custom("graph record bytes bound"));
        }
        if s.is_human_readable() {
            s.serialize_str(&STANDARD.encode(b))
        } else {
            s.serialize_bytes(b)
        }
    }
    pub fn deserialize<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<u8>, D::Error> {
        struct V;
        impl<'de> serde::de::Visitor<'de> for V {
            type Value = Vec<u8>;
            fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                f.write_str("bounded opaque native graph record")
            }
            fn visit_str<E: serde::de::Error>(self, s: &str) -> Result<Vec<u8>, E> {
                if s.len() > MAX_RECORD_BYTES.div_ceil(3) * 4 {
                    return Err(E::custom("graph record encoded bound"));
                }
                let b = STANDARD.decode(s).map_err(E::custom)?;
                if b.len() > MAX_RECORD_BYTES || STANDARD.encode(&b) != s {
                    return Err(E::custom("graph record bytes/canonical bound"));
                }
                Ok(b)
            }
            fn visit_bytes<E: serde::de::Error>(self, b: &[u8]) -> Result<Vec<u8>, E> {
                if b.len() > MAX_RECORD_BYTES {
                    return Err(E::custom("graph record bytes bound"));
                }
                Ok(b.to_vec())
            }
            fn visit_byte_buf<E: serde::de::Error>(self, b: Vec<u8>) -> Result<Vec<u8>, E> {
                if b.len() > MAX_RECORD_BYTES {
                    return Err(E::custom("graph record bytes bound"));
                }
                Ok(b)
            }
        }
        if d.is_human_readable() {
            d.deserialize_str(V)
        } else {
            d.deserialize_bytes(V)
        }
    }
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("TOPO graph history: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
}
fn text(s: &str, n: usize) -> bool {
    !s.trim().is_empty() && s.len() <= n && !s.chars().any(char::is_control)
}
fn bound(v: &impl Serialize, max: usize) -> Result<(), Error> {
    struct Count(usize, usize);
    impl Write for Count {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.0 = self
                .0
                .checked_add(b.len())
                .filter(|n| *n <= self.1)
                .ok_or_else(|| std::io::Error::other("TOPO graph byte limit"))?;
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(Count(0, max), v).map_err(|_| Error::TooLarge)
}
fn shape(r: &CaptureRequest) -> Result<(), Error> {
    topo::physical::steps_preflight(&r.request.steps, &r.request.engine_limits)
        .map_err(|e| bad(&e.to_string()))?;
    if !text(&r.name, 256) || r.bindings.len() > 32 {
        return Err(Error::TooLarge);
    }
    let mut seen = BTreeSet::new();
    for b in &r.bindings {
        if !text(&b.role, 256) || !seen.insert(&b.role) {
            return Err(bad("distinct bounded source roles required"));
        }
        run_material::event_shape(&b.event)?;
        if let Some(o) = &b.origin {
            if o.event() != &b.event {
                return Err(bad("ordinary PCR origin names another event"));
            }
        }
    }
    let mut used = BTreeSet::new();
    for m in &r.request.recipe.materials {
        match m {
            topo::MaterialRecipe::Prepared { strands, .. } => {
                for p in strands.iter().flat_map(|s| &s.segments) {
                    used.insert(p.source_role.as_str());
                }
            }
            topo::MaterialRecipe::ExactPcr { template_role, .. } => {
                used.insert(template_role.as_str());
            }
            topo::MaterialRecipe::SavedGraph { .. } => {}
        }
    }
    if used.len() != r.bindings.len() || r.bindings.iter().any(|b| !used.contains(b.role.as_str()))
    {
        return Err(bad(
            "graph capture bindings must exactly match used external roles",
        ));
    }
    Ok(())
}
impl Record {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        shape(&self.request)?;
        state_shape(&self.state)?;
        if self.schema_version != 2
            || !text(&self.engine_version, 128)
            || !text(&self.scientific_engine, 128)
            || ![
                &self.context_revision,
                &self.source_workspace_fingerprint,
                &self.state_revision,
            ]
            .iter()
            .all(|s| hex(s))
        {
            return Err(bad("unsupported graph record schema/binding"));
        }
        bound(self, MAX_RECORD_BYTES - 44)?;
        let context = assembly_recipe::hash(&(
            &self.engine_version,
            &self.scientific_engine,
            &self.source_workspace_fingerprint,
            &self.request,
        ))?;
        if context != self.context_revision {
            return Err(bad("graph record source/request context differs"));
        }
        Ok(())
    }
}
pub fn save_record(r: &Record) -> Result<Vec<u8>, Error> {
    r.validate_documentary()?;
    let body = serde_json::to_vec(r).map_err(|e| bad(&e.to_string()))?;
    let mut out = Vec::with_capacity(44 + body.len());
    out.extend_from_slice(MAGIC);
    out.extend_from_slice(&(body.len() as u32).to_le_bytes());
    out.extend_from_slice(blake3::hash(&body).as_bytes());
    out.extend(body);
    Ok(out)
}
pub fn load_record(bytes: &[u8]) -> Result<Record, Error> {
    if bytes.len() < 44
        || bytes.len() > MAX_RECORD_BYTES
        || &bytes[..8] != MAGIC
        || u32::from_le_bytes(bytes[8..12].try_into().unwrap()) as usize != bytes.len() - 44
        || blake3::hash(&bytes[44..]).as_bytes() != &bytes[12..44]
    {
        return Err(bad("graph record frame/hash differs"));
    }
    let r: Record = serde_json::from_slice(&bytes[44..]).map_err(|e| bad(&e.to_string()))?;
    r.validate_documentary()?;
    if save_record(&r)? != bytes {
        return Err(bad("graph record is not canonical"));
    }
    Ok(r)
}
impl Capsule {
    pub fn validate_shallow(&self) -> Result<Record, Error> {
        if self.schema_version != 2
            || !hex(&self.selection.record_hash)
            || !hex(&self.selection.state_revision)
            || self.record_bytes.len() > MAX_RECORD_BYTES
            || self.source_snapshots.len() > 32
        {
            return Err(bad("graph capsule shape/version"));
        }
        crate::model::bounds::model_preflight(
            &self
                .source_snapshots
                .iter()
                .map(|i| &i.molecule)
                .collect::<Vec<_>>(),
        )?;
        bound(self, MAX_CAPSULE_BYTES)?;
        let r = load_record(&self.record_bytes)?;
        if blake3::hash(&self.record_bytes).to_hex().as_str() != self.selection.record_hash
            || r.state_revision != self.selection.state_revision
        {
            return Err(bad("graph capsule selection differs"));
        }
        snapshot_binding(&r, &self.source_snapshots)?;
        Ok(r)
    }
}
impl Capsule {
    pub fn strand_count(&self) -> Result<usize, Error> {
        run_material::native_origins_preflight([run_material::NativeOriginRef::TopoGraph(self)])?;
        fn count(c: &Capsule, depth: usize) -> Result<usize, Error> {
            if depth > 32 {
                return Err(Error::TooLarge);
            }
            let r = c.validate_shallow()?;
            let mut n = 0usize;
            for m in &r.request.request.recipe.materials {
                n = n
                    .checked_add(match m {
                        topo::MaterialRecipe::Prepared { strands, .. } => strands.len(),
                        topo::MaterialRecipe::ExactPcr { .. } => 2,
                        topo::MaterialRecipe::SavedGraph { origin, .. } => {
                            count(origin, depth + 1)?
                        }
                    })
                    .filter(|n| *n <= 64)
                    .ok_or(Error::TooLarge)?;
            }
            Ok(n)
        }
        count(self, 0)
    }
}
fn snapshot_binding(r: &Record, inputs: &[topo::Input]) -> Result<(), Error> {
    if r.request.bindings.len() != inputs.len() {
        return Err(bad("graph snapshots/bindings differ"));
    }
    for (b, i) in r.request.bindings.iter().zip(inputs) {
        if b.role != i.role || b.event.revision != i.molecule.revision() {
            return Err(bad("graph source snapshot differs from its exact event"));
        }
    }
    Ok(())
}
pub(crate) fn attest_record_shallow(
    w: &Workspace,
    r: &Record,
    snapshots: Option<&[topo::Input]>,
) -> Result<(), Error> {
    r.validate_documentary()?;
    let bytes = save_record(r)?;
    let hash = blake3::hash(&bytes).to_hex().to_string();
    let a = w
        .attachments
        .get(&hash)
        .ok_or_else(|| bad("graph record is not retained in this native workspace"))?;
    if a.kind != WORKSPACE_ATTACHMENT_KIND || a.media_type != MEDIA_TYPE || a.bytes != bytes {
        return Err(bad("native graph record owner differs"));
    }
    for b in &r.request.bindings {
        let source = run_material::resolve(w, &b.event)?;
        if let Some(inputs) = snapshots {
            let i = inputs
                .iter()
                .find(|i| i.role == b.role)
                .ok_or_else(|| bad("missing closed snapshot"))?;
            if &i.molecule != source {
                return Err(bad("retained graph snapshot/source event differs"));
            }
        }
    }
    if let Some(s) = snapshots {
        snapshot_binding(r, s)?
    }
    Ok(())
}
pub(crate) fn documentary_capsule_parts(w: &Workspace, s: &Selection) -> Result<Record, Error> {
    if !hex(&s.record_hash) || !hex(&s.state_revision) {
        return Err(bad("invalid graph selection"));
    }
    let r = load(w, &s.record_hash)?;
    if r.state_revision != s.state_revision {
        return Err(bad("graph selection state differs"));
    }
    attest_record_shallow(w, &r, None)?;
    Ok(r)
}
pub fn load(w: &Workspace, hash: &str) -> Result<Record, Error> {
    if !hex(hash) {
        return Err(bad("invalid graph record hash"));
    }
    let a = w
        .attachments
        .get(hash)
        .ok_or_else(|| bad("graph record absent"))?;
    if a.kind != WORKSPACE_ATTACHMENT_KIND
        || a.media_type != MEDIA_TYPE
        || blake3::hash(&a.bytes).to_hex().as_str() != hash
    {
        return Err(bad("graph record owner/hash differs"));
    }
    load_record(&a.bytes)
}
pub fn source_bindings(r: &Record) -> Result<Vec<EventBinding>, Error> {
    let mut runs = vec![];
    let mut derived = vec![];
    let mut graphs = vec![];
    crate::run_material_capsule::graph_ancestry_walk(
        r,
        |o, _| {
            runs.push(o.target_event.clone());
            runs.push(o.selection.event.clone());
            runs.extend(o.inputs.iter().map(|i| i.event.clone()));
            Ok(())
        },
        |d, _| {
            derived.extend(d.request.sources.iter().map(|source| source.event.clone()));
            derived.push(d.product.clone());
            Ok(())
        },
        |g, _| {
            graphs.extend(g.request.bindings.iter().map(|b| b.event.clone()));
            Ok(())
        },
    )?;
    let mut unique = std::collections::BTreeMap::new();
    for e in runs.into_iter().chain(derived).chain(graphs) {
        unique.insert(
            (
                e.document_id.clone(),
                e.event_id.clone(),
                e.revision.clone(),
            ),
            e,
        );
    }
    Ok(unique.into_values().collect())
}
pub fn required_graph_records(r: &Record) -> Result<Vec<Selection>, Error> {
    let own = blake3::hash(&save_record(r)?).to_hex().to_string();
    let mut selections = std::collections::BTreeMap::new();
    crate::run_material_capsule::graph_ancestry_walk(
        r,
        |_, _| Ok(()),
        |_, _| Ok(()),
        |g, _| {
            let hash = blake3::hash(&save_record(g)?).to_hex().to_string();
            if hash != own {
                selections.insert(
                    hash.clone(),
                    Selection {
                        record_hash: hash,
                        state_revision: g.state_revision.clone(),
                    },
                );
            }
            Ok(())
        },
    )?;
    Ok(selections.into_values().collect())
}
pub fn required_run_ids(r: &Record) -> Result<Vec<String>, Error> {
    let mut ids = BTreeSet::new();
    crate::run_material_capsule::graph_ancestry_walk(
        r,
        |o, _| {
            ids.insert(o.selection.run_id.clone());
            Ok(())
        },
        |_, _| Ok(()),
        |_, _| Ok(()),
    )?;
    Ok(ids.into_iter().collect())
}
fn state(report: &topo::TopoReport) -> Result<(String, StateSummary), Error> {
    let p = &report.physical.product;
    let mut remaining = BTreeSet::new();
    for c in &p.chains {
        if let topo::physical::ChainEnds::Linear {
            three_prime: topo::physical::ThreePrime::Phosphotyrosyl { activation_id, .. },
            ..
        } = &c.ends
        {
            remaining.insert(activation_id.clone());
        }
    }
    let s = StateSummary {
        source_bases: report.physical.source_bases,
        chains: p.chains.len(),
        components: p.components.len(),
        remaining_activations: remaining.into_iter().collect(),
        released_activations: p.released_activations.clone(),
        findings: report.findings.clone(),
        assumptions: report
            .preparation
            .assumptions
            .iter()
            .chain(&report.physical.assumptions)
            .take(128)
            .map(|s| s.chars().take(4096).collect())
            .collect(),
    };
    Ok((assembly_recipe::hash(&(ENGINE_VERSION, p))?, s))
}
fn summary(hash: &str, r: &Record) -> Result<Summary, Error> {
    Ok(Summary {
        selection: Selection {
            record_hash: hash.into(),
            state_revision: r.state_revision.clone(),
        },
        name: r.request.name.clone(),
        source_events: source_bindings(r)?,
        source_bases: r.state.source_bases,
        chains: r.state.chains,
        components: r.state.components,
        remaining_activations: r.state.remaining_activations.clone(),
        released_activations: r.state.released_activations.clone(),
        incomplete: r.incomplete,
        findings: r.state.findings.clone(),
        assumptions: r.state.assumptions.clone(),
    })
}
fn sources(w: &Workspace, r: &CaptureRequest) -> Result<Vec<topo::Input>, Error> {
    r.bindings
        .iter()
        .map(|b| {
            Ok(topo::Input {
                role: b.role.clone(),
                molecule: run_material::resolve(w, &b.event)?.clone(),
            })
        })
        .collect()
}
pub fn preflight(w: &Workspace, r: &CaptureRequest) -> Result<(), Error> {
    shape(r)?;
    let roots = r
        .bindings
        .iter()
        .filter_map(|b| b.origin.as_ref())
        .map(|o| match o {
            MaterialOrigin::Run(c) => run_material::NativeOriginRef::Run(c),
            MaterialOrigin::Derived(c) => run_material::NativeOriginRef::Derived(c),
        })
        .chain(r.request.recipe.materials.iter().filter_map(|m| {
            if let topo::MaterialRecipe::SavedGraph { origin, .. } = m {
                Some(run_material::NativeOriginRef::TopoGraph(origin))
            } else {
                None
            }
        }));
    run_material::native_origins_preflight(roots)?;
    w.validate()?;
    for b in &r.bindings {
        let m = run_material::resolve(w, &b.event)?;
        if let Some(o) = &b.origin {
            if o.event() != &b.event {
                return Err(bad("source origin differs"));
            }
            o.validate_for_workspace(w, &b.event)?;
        }
        if m.sequence.alphabet() != crate::Alphabet::Dna {
            return Err(bad("DNA source required"));
        }
    }
    for m in &r.request.recipe.materials {
        if let topo::MaterialRecipe::SavedGraph { origin, .. } = m {
            validate_capsule_for_workspace(w, origin)?
        }
    }
    bound(r, MAX_RECORD_BYTES - 44)?;
    Ok(())
}
pub fn revision(w: &Workspace, r: &CaptureRequest) -> Result<String, Error> {
    preflight(w, r)?;
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        run_lineage::workspace_fingerprint(w)?,
        r,
    ))
}
pub fn prepare_capture(w: &Workspace, r: &CaptureRequest) -> Result<PreparedCapture, Error> {
    let context_revision = revision(w, r)?;
    let input = sources(w, r)?;
    let report = topo::simulate_bound(&input, &r.request, &r.bindings)?;
    let (state_revision, state) = state(&report)?;
    let fingerprint = run_lineage::workspace_fingerprint(w)?;
    let record = Record {
        schema_version: 2,
        engine_version: ENGINE_VERSION.into(),
        scientific_engine: verification::ENGINE_VERSION.into(),
        source_workspace_fingerprint: fingerprint.clone(),
        context_revision: context_revision.clone(),
        request: r.clone(),
        state_revision,
        incomplete: report.incomplete,
        state,
    };
    let bytes = save_record(&record)?;
    let hash = blake3::hash(&bytes).to_hex().to_string();
    if w.attachments
        .values()
        .filter(|a| a.kind == WORKSPACE_ATTACHMENT_KIND)
        .count()
        >= MAX_RECORDS
        && !w.attachments.contains_key(&hash)
    {
        return Err(Error::TooLarge);
    }
    let mut next = w.clone();
    if let Some(old) = next.attachments.get(&hash) {
        if old.kind != WORKSPACE_ATTACHMENT_KIND
            || old.media_type != MEDIA_TYPE
            || old.bytes != bytes
        {
            return Err(bad("graph attachment collision"));
        }
    } else {
        next.attachments.insert(
            hash.clone(),
            WorkspaceAttachment {
                kind: WORKSPACE_ATTACHMENT_KIND.into(),
                name: r.name.clone(),
                media_type: MEDIA_TYPE.into(),
                bytes,
            },
        );
    }
    next.validate()?;
    let receipt = Receipt {
        context_revision,
        source_workspace_fingerprint: fingerprint,
        record_hash: hash.clone(),
        selection: Selection {
            record_hash: hash.clone(),
            state_revision: record.state_revision.clone(),
        },
        summary: summary(&hash, &record)?,
    };
    receipt.validate_binding(&receipt.context_revision)?;
    Ok(PreparedCapture {
        workspace: next,
        receipt,
    })
}
pub fn documentary_capsule(w: &Workspace, s: &Selection) -> Result<Capsule, Error> {
    let mut origins = run_material::documentary_origins(w, &[], std::slice::from_ref(s))?;
    if origins.graphs.len() != 1 {
        return Err(bad("documentary graph factory count differs"));
    }
    Ok(origins.graphs.remove(0))
}
pub fn validate_capsule_for_workspace(w: &Workspace, c: &Capsule) -> Result<(), Error> {
    crate::run_material::validate_native_origins(w, [run_material::NativeOriginRef::TopoGraph(c)])
}
pub fn replay(c: &Capsule) -> Result<FreshGraph, Error> {
    run_material::native_origins_preflight([run_material::NativeOriginRef::TopoGraph(c)])?;
    let r = c.validate_shallow()?;
    if r.engine_version != ENGINE_VERSION || r.scientific_engine != verification::ENGINE_VERSION {
        return Err(bad(
            "historical graph model cannot establish current source material",
        ));
    }
    let report =
        topo::simulate_bound(&c.source_snapshots, &r.request.request, &r.request.bindings)?;
    let (revision, found) = state(&report)?;
    if revision != c.selection.state_revision
        || assembly_recipe::hash(&found)? != assembly_recipe::hash(&r.state)?
        || report.incomplete != r.incomplete
    {
        return Err(bad(
            "fresh graph differs from exact retained source/transcript",
        ));
    }
    Ok(FreshGraph {
        incomplete: report.incomplete,
        report,
    })
}
pub fn list_revision(w: &Workspace) -> Result<String, Error> {
    w.validate()?;
    let hashes = w
        .attachments
        .iter()
        .filter(|(_, a)| a.kind == WORKSPACE_ATTACHMENT_KIND)
        .map(|(h, _)| h)
        .collect::<Vec<_>>();
    if hashes.len() > MAX_RECORDS {
        return Err(Error::TooLarge);
    }
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        "list",
        run_lineage::workspace_fingerprint(w)?,
        hashes,
    ))
}
pub fn list(w: &Workspace) -> Result<Listed, Error> {
    let context_revision = list_revision(w)?;
    let mut entries = vec![];
    for (h, a) in &w.attachments {
        if a.kind == WORKSPACE_ATTACHMENT_KIND {
            let r = load(w, h)?;
            attest_record_shallow(w, &r, None)?;
            entries.push(summary(h, &r)?);
        }
    }
    let out = Listed {
        context_revision,
        source_workspace_fingerprint: run_lineage::workspace_fingerprint(w)?,
        entries,
    };
    out.validate_binding(&out.context_revision)?;
    Ok(out)
}
pub fn inspect_revision(w: &Workspace, s: &Selection) -> Result<String, Error> {
    documentary_capsule_parts(w, s)?;
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        "inspect",
        verification::ENGINE_VERSION,
        run_lineage::workspace_fingerprint(w)?,
        s,
    ))
}
pub fn inspect(w: &Workspace, s: &Selection) -> Result<Inspection, Error> {
    let context_revision = inspect_revision(w, s)?;
    let c = documentary_capsule(w, s)?;
    let record = c.validate_shallow()?;
    let fresh = replay(&c)?;
    let out = Inspection {
        context_revision,
        source_workspace_fingerprint: run_lineage::workspace_fingerprint(w)?,
        selection: s.clone(),
        summary: summary(&s.record_hash, &record)?,
        record,
        report: fresh.report,
    };
    out.validate_binding(&out.context_revision)?;
    Ok(out)
}
fn state_shape(s: &StateSummary) -> Result<(), Error> {
    if s.source_bases > 500_000
        || s.chains > 128
        || s.components > 128
        || s.remaining_activations.len() > 64
        || s.released_activations.len() > 64
        || s.remaining_activations
            .iter()
            .chain(&s.released_activations)
            .any(|s| !text(s, 256))
        || s.findings.len() > 128
        || s.assumptions.len() > 256
        || s.findings
            .iter()
            .chain(&s.assumptions)
            .any(|s| s.len() > 16384)
    {
        return Err(bad("graph state summary bounds"));
    }
    Ok(())
}
pub fn validate_workspace_records(w: &Workspace) -> Result<(), Error> {
    let mut count = 0;
    for (hash, a) in &w.attachments {
        if a.kind == WORKSPACE_ATTACHMENT_KIND {
            count += 1;
            if count > MAX_RECORDS {
                return Err(Error::TooLarge);
            }
            let record = load(w, hash)?;
            attest_record_shallow(w, &record, None)?;
            let roots = record
                .request
                .bindings
                .iter()
                .filter_map(|b| b.origin.as_ref())
                .map(|o| match o {
                    MaterialOrigin::Run(c) => run_material::NativeOriginRef::Run(c),
                    MaterialOrigin::Derived(c) => run_material::NativeOriginRef::Derived(c),
                })
                .chain(
                    record
                        .request
                        .request
                        .recipe
                        .materials
                        .iter()
                        .filter_map(|m| {
                            if let topo::MaterialRecipe::SavedGraph { origin, .. } = m {
                                Some(run_material::NativeOriginRef::TopoGraph(origin))
                            } else {
                                None
                            }
                        }),
                );
            run_material::validate_native_origins(w, roots)?;
        }
    }
    Ok(())
}
fn summary_shape(s: &Summary) -> Result<(), Error> {
    if !hex(&s.selection.record_hash)
        || !hex(&s.selection.state_revision)
        || !text(&s.name, 256)
        || s.source_events.len() > 4096
        || s.chains > 128
        || s.components > 128
        || s.source_bases > 500_000
        || s.remaining_activations.len() > 64
        || s.released_activations.len() > 64
        || s.findings.len() > 128
        || s.assumptions.len() > 256
        || s.findings
            .iter()
            .chain(&s.assumptions)
            .any(|s| s.len() > 16384)
    {
        return Err(bad("graph summary bounds"));
    }
    Ok(())
}
impl Receipt {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        summary_shape(&self.summary)?;
        if self.context_revision != expected
            || !hex(expected)
            || !hex(&self.source_workspace_fingerprint)
            || self.record_hash != self.selection.record_hash
            || self.summary.selection != self.selection
        {
            return Err(bad("graph capture receipt changed"));
        }
        Ok(())
    }
}
impl Listed {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.context_revision != expected
            || !hex(expected)
            || !hex(&self.source_workspace_fingerprint)
            || self.entries.len() > MAX_RECORDS
        {
            return Err(bad("graph list context/bounds"));
        }
        for e in &self.entries {
            summary_shape(e)?
        }
        bound(self, MAX_CAPSULE_BYTES)
    }
}
impl Inspection {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.record.validate_documentary()?;
        summary_shape(&self.summary)?;
        let (state, summary) = state(&self.report)?;
        if self.context_revision != expected
            || !hex(expected)
            || !hex(&self.source_workspace_fingerprint)
            || self.selection != self.summary.selection
            || self.record.state_revision != self.selection.state_revision
            || state != self.selection.state_revision
            || assembly_recipe::hash(&summary)? != assembly_recipe::hash(&self.record.state)?
        {
            return Err(bad("graph inspection correspondence changed"));
        }
        bound(self, MAX_CAPSULE_BYTES)
    }
}

impl PreparedCapture {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.receipt.validate_binding(expected)?;
        let r = load(&self.workspace, &self.receipt.record_hash)?;
        if r.context_revision != expected
            || r.source_workspace_fingerprint != self.receipt.source_workspace_fingerprint
            || summary(&self.receipt.record_hash, &r)?.selection != self.receipt.selection
        {
            return Err(bad("captured graph workspace/receipt differs"));
        }
        let mut original = self.workspace.clone();
        original.attachments.remove(&self.receipt.record_hash);
        if run_lineage::workspace_fingerprint(&original)?
            != self.receipt.source_workspace_fingerprint
        {
            return Err(bad("graph capture changed native biological/history scope"));
        }
        Ok(())
    }
}
