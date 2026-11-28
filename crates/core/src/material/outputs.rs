//! Atomic, presentation-aware output capture. Historical run products are never pruned.
use crate::{
    documents::{
        history::DocumentHistory,
        workspace::{self, Workspace, WorkspaceAttachment},
    },
    material::lineage::{self as run_lineage, EventBinding, InputEventBinding, RunReceipt},
    planning::{oligos as plan_oligos, plan::DesignPlan},
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const ENGINE_VERSION: &str = "cloning-outputs-1";
pub const VISIBILITY_KIND: &str = "cloning-output-visibility-v1";
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Options {pub name:String,pub show_prepared:bool,pub save_primers:bool,pub allow_incomplete:bool}}
impl Default for Options {
    fn default() -> Self {
        Self {
            name: "Assembly".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        }
    }
}
impl Options {
    pub fn validate(&self) -> Result<(), Error> {
        if self.name.is_empty()
            || self.name.trim() != self.name
            || self.name.len() > 256
            || self.name.chars().any(char::is_control)
        {
            return Err(bad("name must be nonempty, trimmed and at most 256 bytes"));
        }
        Ok(())
    }
}
dto! {pub struct Receipt {pub context_revision:String,pub run:RunReceipt,pub construct:EventBinding,pub primer_attachment:Option<String>,pub visible_products:Vec<String>,pub hidden_products:Vec<String>,pub source_workspace_fingerprint:String}}
dto! {pub struct Commit {#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))] pub workspace:Workspace,pub receipt:Receipt}}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Visibility {
    version: u8,
    hidden: Vec<String>,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("cloning outputs: {s}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
}
fn workspace_guard(w: &Workspace) -> Result<(), Error> {
    if w.documents.len() > workspace::MAX_DOCUMENTS
        || w.revisions.len() > workspace::MAX_REVISIONS
        || w.history.len() > workspace::MAX_HISTORY_EVENTS
        || w.attachments.len() > workspace::MAX_ATTACHMENTS
        || w.candidates.len() > workspace::MAX_CANDIDATES
    {
        return Err(Error::TooLarge);
    }
    for m in w.revisions.values() {
        // Unrelated documents retain native annotation limits. The actual
        // requested plan and target receive scientific guards below.
        m.validate()?;
    }
    for c in &w.candidates {
        if c.plan.inputs.len() > 1000 || c.plan.steps.len() > 1000 {
            return Err(Error::TooLarge);
        }
        for i in &c.plan.inputs {
            i.molecule.validate()?;
        }
    }
    w.validate()
}
fn resolve<'a>(w: &'a Workspace, b: &EventBinding) -> Result<&'a Molecule, Error> {
    if b.document_id.len() > 128
        || !digest(&b.revision)
        || !b.event_id.strip_prefix("event-").is_some_and(digest)
    {
        return Err(bad("invalid event binding"));
    }
    let event = w
        .history
        .iter()
        .find(|e| e.id == b.event_id)
        .ok_or_else(|| bad("missing exact event"))?;
    if event.document_id != b.document_id || event.revision != b.revision {
        return Err(bad("event binding mismatch"));
    }
    w.revisions
        .get(&b.revision)
        .ok_or_else(|| bad("missing event revision"))
}
/// Presentation state only; no scientific or material authority is restored.
pub fn load_hidden(w: &Workspace) -> Result<BTreeSet<String>, Error> {
    if w.attachments.len() > workspace::MAX_ATTACHMENTS
        || w.documents.len() > workspace::MAX_DOCUMENTS
    {
        return Err(Error::TooLarge);
    }
    let mut assets = w
        .attachments
        .iter()
        .filter(|(_, a)| a.kind == VISIBILITY_KIND);
    let Some((hash, a)) = assets.next() else {
        return Ok(BTreeSet::new());
    };
    if assets.next().is_some() {
        return Err(bad("multiple visibility attachments"));
    }
    if a.bytes.len() > 256 * 1024 {
        return Err(Error::TooLarge);
    }
    if blake3::hash(&a.bytes).to_hex().as_str() != hash || a.media_type != "application/json" {
        return Err(bad("invalid visibility attachment identity"));
    }
    let v: Visibility =
        serde_json::from_slice(&a.bytes).map_err(|_| bad("invalid visibility attachment"))?;
    if v.version != 1 || v.hidden.len() > workspace::MAX_DOCUMENTS {
        return Err(bad("invalid visibility version or count"));
    }
    let mut result = BTreeSet::new();
    for id in v.hidden {
        if id.len() > 128 || !w.documents.iter().any(|d| d.id == id) || !result.insert(id) {
            return Err(bad("invalid or duplicate hidden document"));
        }
    }
    Ok(result)
}
struct Bounded {
    bytes: Vec<u8>,
}
impl Write for Bounded {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > workspace::MAX_ENTRY_BYTES)
        {
            return Err(std::io::Error::other("output cap"));
        }
        self.bytes.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encoded(v: &impl Serialize) -> Result<Vec<u8>, Error> {
    let mut w = Bounded { bytes: Vec::new() };
    serde_json::to_writer(&mut w, v).map_err(|_| Error::TooLarge)?;
    Ok(w.bytes)
}
fn attach(w: &mut Workspace, kind: &str, name: String, bytes: Vec<u8>) -> Result<String, Error> {
    let hash = blake3::hash(&bytes).to_hex().to_string();
    if let Some(old) = w.attachments.get(&hash) {
        if old.kind != kind || old.media_type != "application/json" || old.bytes != bytes {
            return Err(bad("attachment identity conflict"));
        }
    } else {
        w.attachments.insert(
            hash.clone(),
            WorkspaceAttachment {
                kind: kind.into(),
                name,
                media_type: "application/json".into(),
                bytes,
            },
        );
    }
    Ok(hash)
}
/// Hash current workspace and exact historical bindings; no biological reconstruction.
pub fn finalize_revision(
    w: &Workspace,
    plan: &DesignPlan,
    target: EventBinding,
    inputs: Vec<InputEventBinding>,
    options: &Options,
) -> Result<String, Error> {
    options.validate()?;
    workspace_guard(w)?;
    if inputs.len() > 128 {
        return Err(Error::TooLarge);
    }
    let actual = resolve(w, &target)?;
    crate::cloning::assembly::junctions::source_preflight(plan, actual)?;
    crate::documents::transfer::Candidate::new(plan.clone()).validate()?;
    if actual.revision() != plan.target_revision || inputs.len() != plan.inputs.len() {
        return Err(bad("target revision or input binding count mismatch"));
    }
    for source in &plan.inputs {
        let mut matching = inputs.iter().filter(|b| b.input_id == source.id);
        let binding = matching
            .next()
            .ok_or_else(|| bad("missing input binding"))?;
        if matching.next().is_some() || resolve(w, &binding.event)?.revision() != source.revision {
            return Err(bad("duplicate or stale input binding"));
        }
    }
    if options.save_primers {
        plan_oligos::preflight(plan, actual, &plan_oligos::Options::default())?;
    }
    for input in &inputs {
        if input.input_id.len() > 128 {
            return Err(Error::TooLarge);
        }
        resolve(w, &input.event)?;
    }
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        run_lineage::workspace_fingerprint(w)?,
        crate::planning::plan::DesignPlan::revision(plan),
        target,
        inputs,
        options,
    ))
}
/// Freshly reconstruct and capture all steps, then publish presentation changes atomically.
pub fn finalize(
    w: &Workspace,
    plan: &DesignPlan,
    target: EventBinding,
    inputs: Vec<InputEventBinding>,
    options: &Options,
) -> Result<Commit, Error> {
    let context_revision = finalize_revision(w, plan, target.clone(), inputs.clone(), options)?;
    let mut hidden = load_hidden(w)?;
    let source_workspace_fingerprint = run_lineage::workspace_fingerprint(w)?;
    let actual_target = resolve(w, &target)?;
    let captured = run_lineage::record_run(w, plan, target, inputs, options.allow_incomplete)?;
    let mut next = captured.workspace;
    let final_product = captured
        .receipt
        .products
        .last()
        .ok_or_else(|| bad("run has no final product"))?;
    let recorded_final_id = final_product.document_id.clone();
    let original = next
        .revisions
        .get(&final_product.revision)
        .ok_or_else(|| bad("missing original product"))?;
    let current = next.head(&recorded_final_id)?;
    let edited = !same_reference(current, original);
    let mut final_id = recorded_final_id.clone();
    if edited {
        // Reuse only an exact named descendant of this recorded product event.
        if let Some(d) = next.documents.iter().find(|d| {
            next.revisions
                .get(&d.head)
                .is_some_and(|m| m.name == options.name && same_reference(m, original))
                && d.head_event
                    .as_ref()
                    .is_some_and(|e| ancestry(&next, e, &final_product.event_id))
        }) {
            final_id = d.id.clone();
        } else {
            let mut histories = DocumentHistory::restore_all(&next)?;
            let h = histories
                .iter()
                .find(|h| h.id == recorded_final_id)
                .ok_or_else(|| bad("missing product history"))?;
            let base = format!("construct-{}", captured.receipt.run_id);
            final_id = base.clone();
            let mut suffix = 1usize;
            while next.documents.iter().any(|d| d.id == final_id) {
                final_id = format!("{base}-{suffix}");
                suffix += 1;
            }
            let mut fork = h.fork_event(&final_product.event_id, final_id.clone())?;
            let mut named = fork.molecule.clone();
            named.name = options.name.clone();
            fork.commit(named, "Name finalized assembly construct")?;
            histories.push(fork);
            next.documents.clear();
            next.history.clear();
            for h in &histories {
                h.contribute(&mut next)?;
            }
        }
    } else if current.name != options.name {
        let mut histories = DocumentHistory::restore_all(&next)?;
        let h = histories
            .iter_mut()
            .find(|h| h.id == final_id)
            .ok_or_else(|| bad("missing product history"))?;
        let mut named = h.molecule.clone();
        named.name = options.name.clone();
        h.commit(named, "Name finalized assembly construct")?;
        next.documents.clear();
        next.history.clear();
        for h in &histories {
            h.contribute(&mut next)?;
        }
    }
    let primer_attachment = if options.save_primers {
        let report = plan_oligos::build(plan, actual_target, &plan_oligos::Options::default())?;
        if !report.complete || report.oligos.records.is_empty() {
            return Err(bad(
                "primer collection requires complete, nonempty declared oligos",
            ));
        }
        report.oligos.validate()?;
        Some(attach(
            &mut next,
            "oligo-table-v1",
            format!("Assembly primers {}", captured.receipt.run_id),
            encoded(&report.oligos)?,
        )?)
    } else {
        None
    };
    for p in &captured.receipt.products {
        if options.show_prepared || p.document_id == recorded_final_id {
            hidden.remove(&p.document_id);
        } else {
            hidden.insert(p.document_id.clone());
        }
    }
    hidden.remove(&final_id);
    let mut visible_products: Vec<String> = captured
        .receipt
        .products
        .iter()
        .filter(|p| !hidden.contains(&p.document_id))
        .map(|p| p.document_id.clone())
        .collect();
    if final_id != recorded_final_id {
        visible_products.push(final_id.clone());
    }
    let hidden_products = captured
        .receipt
        .products
        .iter()
        .filter(|p| hidden.contains(&p.document_id))
        .map(|p| p.document_id.clone())
        .collect();
    next.attachments.retain(|_, a| a.kind != VISIBILITY_KIND);
    attach(
        &mut next,
        VISIBILITY_KIND,
        "Assembly document visibility".into(),
        encoded(&Visibility {
            version: 1,
            hidden: hidden.into_iter().collect(),
        })?,
    )?;
    next.active_document = final_id.clone();
    let d = next
        .documents
        .iter()
        .find(|d| d.id == final_id)
        .ok_or_else(|| bad("missing final document"))?;
    let construct = EventBinding {
        document_id: final_id,
        event_id: d
            .head_event
            .clone()
            .ok_or_else(|| bad("missing head event"))?,
        revision: d.head.clone(),
    };
    let commit = Commit {
        workspace: next,
        receipt: Receipt {
            context_revision: context_revision.clone(),
            run: captured.receipt,
            construct,
            primer_attachment,
            visible_products,
            hidden_products,
            source_workspace_fingerprint,
        },
    };
    commit.validate_binding(&context_revision)?;
    Ok(commit)
}
impl Commit {
    /// Readonly structure and historical identity check, never a fresh biological verdict.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let r = &self.receipt;
        let w = &self.workspace;
        if !digest(expected)
            || r.context_revision != expected
            || !digest(&r.source_workspace_fingerprint)
            || r.run.source_workspace_fingerprint != r.source_workspace_fingerprint
        {
            return Err(bad("receipt context mismatch"));
        }
        workspace_guard(w)?;
        let hidden = load_hidden(w)?;
        let run = run_lineage::load_run(w, &r.run.run_id)?;
        let asset = w
            .attachments
            .get(&r.run.attachment_hash)
            .ok_or_else(|| bad("missing run attachment"))?;
        if asset.kind != run_lineage::WORKSPACE_ATTACHMENT_KIND
            || run_lineage::load_run_record(&asset.bytes)?.run_id != run.run_id
            || crate::encoding::bounded_json::hash(&run.products)?
                != crate::encoding::bounded_json::hash(&r.run.products)?
            || run.assessment.status != r.run.recorded_status
            || run.assessment.engine_version != r.run.recorded_engine
        {
            return Err(bad("run receipt mismatch"));
        }
        let last = run.products.last().ok_or_else(|| bad("empty run"))?;
        if w.active_document != r.construct.document_id {
            return Err(bad("construct identity mismatch"));
        }
        let actual = resolve(w, &r.construct)?;
        if !ancestry(w, &r.construct.event_id, &last.event_id) {
            return Err(bad("construct has no recorded-product ancestry"));
        }
        let d = w
            .documents
            .iter()
            .find(|d| d.id == r.construct.document_id)
            .ok_or_else(|| bad("missing construct"))?;
        if d.head != r.construct.revision
            || d.head_event.as_deref() != Some(r.construct.event_id.as_str())
        {
            return Err(bad("construct is not current head"));
        }
        let original = w
            .revisions
            .get(&last.revision)
            .ok_or_else(|| bad("missing recorded product"))?;
        if !same_reference(actual, original) {
            return Err(bad("construct differs from recorded reference beyond name"));
        }
        let mut seen = BTreeSet::new();
        for (ids, is_hidden) in [(&r.visible_products, false), (&r.hidden_products, true)] {
            if ids.len() > 1000 {
                return Err(Error::TooLarge);
            }
            for id in ids {
                if !seen.insert(id)
                    || (!run.products.iter().any(|p| &p.document_id == id)
                        && id != &r.construct.document_id)
                    || hidden.contains(id) != is_hidden
                {
                    return Err(bad("visibility receipt mismatch"));
                }
            }
        }
        if seen.len()
            != run.products.len() + usize::from(last.document_id != r.construct.document_id)
            || !seen.contains(&r.construct.document_id)
            || hidden.contains(&last.document_id)
            || hidden.contains(&r.construct.document_id)
        {
            return Err(bad("incomplete visibility receipt"));
        }
        if let Some(hash) = &r.primer_attachment {
            let a = w
                .attachments
                .get(hash)
                .ok_or_else(|| bad("missing primer attachment"))?;
            if a.kind != "oligo-table-v1" || a.bytes.len() > workspace::MAX_ENTRY_BYTES {
                return Err(bad("invalid primer attachment"));
            }
            let table: crate::formats::oligos::OligoTable =
                serde_json::from_slice(&a.bytes).map_err(|_| bad("invalid primer table"))?;
            table.validate()?;
            if table.records.is_empty() {
                return Err(bad("empty primer table"));
            }
        }
        Ok(())
    }
}

pub(crate) fn ancestry(w: &Workspace, event: &str, ancestor: &str) -> bool {
    let by_id: std::collections::BTreeMap<_, _> =
        w.history.iter().map(|e| (e.id.as_str(), e)).collect();
    let mut seen = BTreeSet::new();
    let mut pending = vec![event];
    while let Some(id) = pending.pop() {
        if id == ancestor {
            return true;
        }
        if !seen.insert(id) {
            continue;
        }
        let Some(e) = by_id.get(id) else {
            return false;
        };
        for p in &e.parents {
            if !seen.contains(p.as_str()) {
                pending.push(p);
            }
        }
    }
    false
}

pub(crate) fn same_reference(a: &Molecule, b: &Molecule) -> bool {
    a.sequence == b.sequence
        && a.topology == b.topology
        && a.features == b.features
        && a.metadata == b.metadata
        && a.original_genbank == b.original_genbank
}
