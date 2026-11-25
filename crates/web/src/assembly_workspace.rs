//! Persistent declared assembly slots. Frozen intent/scientific evidence is
//! component-private and is never restored from a serialized draft.
use crate::{
    cloning_methods::{
        CloningMethodsPanel, MethodsState, Mode, OverlapInput, OverlapPreparation, TargetContext,
    },
    digest::SourceInfo,
    primer_picker::ResolveRequest,
    worker::CachedWorker,
};
use collagene_core::{
    Alphabet, Molecule, Topology, assembly_intent, primer_reuse,
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_lineage::EventBinding,
    transfer::Candidate,
    verification::PlanInput,
};
use leptos::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    io::Write,
    rc::Rc,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
};
use wasm_bindgen_futures::spawn_local;
const MAX_BYTES: usize = 32 << 20;
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CapturedSource {
    pub source: SourceInfo,
    pub molecule: Molecule,
}
mod optional_snapshot {
    use super::*;
    pub fn serialize<S: serde::Serializer>(
        value: &Option<Arc<CapturedSource>>,
        serializer: S,
    ) -> Result<S::Ok, S::Error> {
        value.as_deref().serialize(serializer)
    }
    pub fn deserialize<'de, D: serde::Deserializer<'de>>(
        deserializer: D,
    ) -> Result<Option<Arc<CapturedSource>>, D::Error> {
        Option::<CapturedSource>::deserialize(deserializer).map(|v| v.map(Arc::new))
    }
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DigestGeometry {
    pub top_start: usize,
    pub top_length: usize,
    pub bottom_start: i64,
    pub bottom_length: usize,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Slot {
    pub label: String,
    #[serde(with = "optional_snapshot")]
    pub captured: Option<Arc<CapturedSource>>,
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
    #[serde(default)]
    pub digest_geometry: Option<DigestGeometry>,
}
impl Slot {
    fn empty(label: &str) -> Self {
        Self {
            label: label.into(),
            captured: None,
            start: 0,
            end: 0,
            reverse: false,
            digest_geometry: None,
        }
    }
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Draft {
    pub version: u32,
    pub name: String,
    pub circular: bool,
    #[serde(deserialize_with = "deserialize_slots")]
    pub slots: Vec<Slot>,
    pub active_slot: usize,
    pub methods: MethodsState,
}
fn deserialize_slots<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<Vec<Slot>, D::Error> {
    struct Visitor;
    impl<'de> serde::de::Visitor<'de> for Visitor {
        type Value = Vec<Slot>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("at most12 assembly slots")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            if seq.size_hint().is_some_and(|n| n > 12) {
                return Err(serde::de::Error::custom("at most12 assembly slots"));
            }
            let mut slots = Vec::with_capacity(2);
            while slots.len() < 12 {
                match seq.next_element()? {
                    Some(slot) => slots.push(slot),
                    None => return Ok(slots),
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom("at most12 assembly slots"));
            }
            Ok(slots)
        }
    }
    deserializer.deserialize_seq(Visitor)
}
impl Default for Draft {
    fn default() -> Self {
        let mut value = Self {
            version: 1,
            name: "New assembly".into(),
            circular: true,
            slots: vec![Slot::empty("Backbone"), Slot::empty("Insert")],
            active_slot: 0,
            methods: MethodsState::default(),
        };
        value.sync_rows();
        for row in &mut value.methods.rows {
            row.preparation = OverlapPreparation::ExactPcr;
        }
        value
    }
}
fn token(text: &str, max: usize) -> bool {
    !text.trim().is_empty() && text.len() <= max && !text.chars().any(char::is_control)
}
fn hex(text: &str) -> bool {
    text.len() == 64
        && text
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
struct Bounded {
    bytes: Vec<u8>,
}
struct Count(usize);
impl Write for Count {
    fn write(&mut self, data: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(data.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| std::io::Error::other("Assembly draft exceeds32MiB"))?;
        Ok(data.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
impl Write for Bounded {
    fn write(&mut self, data: &[u8]) -> std::io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(data.len())
            .is_none_or(|n| n > MAX_BYTES)
        {
            return Err(std::io::Error::other("Assembly draft exceeds32MiB"));
        }
        self.bytes.extend_from_slice(data);
        Ok(data.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
impl Draft {
    /// Structural/bounded source consistency only. Historical event IDs require
    /// the host workspace to attest them before installing any candidate.
    pub fn validate(&self) -> Result<(), String> {
        if self.version != 1
            || !token(&self.name, 256)
            || self.slots.is_empty()
            || self.slots.len() > 12
            || self.active_slot >= self.slots.len()
            || self.methods.rows.len() != self.slots.len()
            || self.methods.mode == Mode::Tailing
            || self.methods.circular != self.circular
        {
            return Err("Invalid bounded assembly draft or unsupported draft mode".into());
        }
        let models = self
            .slots
            .iter()
            .filter_map(|slot| slot.captured.as_deref().map(|c| &c.molecule))
            .collect::<Vec<_>>();
        collagene_core::model::bounds::model_preflight(&models).map_err(|e| e.to_string())?;
        for row in &self.methods.rows {
            if let Some(selection) = &row.catalog_selection {
                selection.validate().map_err(|e| e.to_string())?;
            }
            for text in [
                &row.enzymes,
                &row.fragment,
                &row.forward_name,
                &row.reverse_name,
                &row.forward,
                &row.reverse_primer,
                &row.forward_tail,
                &row.reverse_tail,
                &row.interval_start,
                &row.interval_end,
                &row.right_fusion,
            ] {
                if text.len() > 65536 {
                    return Err("Assembly option field exceeds64KiB".into());
                }
            }
        }
        let gate = &self.methods.golden_gate;
        if gate.clamp.len() > 12
            || gate.spacer.len() > 2
            || gate.forward_length.len() > 4
            || gate.reverse_length.len() > 4
        {
            return Err("Type IIS primer settings exceed their field bounds".into());
        }
        for row in [&self.methods.vector, &self.methods.insert] {
            for text in [
                &row.enzymes,
                &row.fragment,
                &row.forward_name,
                &row.reverse_name,
                &row.forward,
                &row.reverse,
                &row.forward_tail,
                &row.reverse_tail,
            ] {
                if text.len() > 65536 {
                    return Err("Assembly option field exceeds64KiB".into());
                }
            }
        }
        for text in [
            &self.methods.overlaps,
            &self.methods.comparisons,
            &self.methods.occurrences,
        ] {
            if text.len() > 65536 {
                return Err("Assembly option field exceeds64KiB".into());
            }
        }
        serde_json::to_writer(&mut Count(0), self).map_err(|e| e.to_string())?;
        let mut bases = 0usize;
        for (i, slot) in self.slots.iter().enumerate() {
            if !token(&slot.label, 80) {
                return Err("Slot labels must be nonblank and at most80 bytes".into());
            }
            let row = &self.methods.rows[i];
            if row.document != i
                || row.reverse != slot.reverse
                || slot.start.checked_add(1).map(|n| n.to_string()).as_deref()
                    != Some(row.interval_start.as_str())
                || row.interval_end
                    != if slot.captured.is_some() {
                        slot.end.to_string()
                    } else {
                        String::new()
                    }
            {
                return Err("Draft preparation source/order differs from its pinned slots".into());
            }
            for text in [
                &row.enzymes,
                &row.fragment,
                &row.forward_name,
                &row.reverse_name,
                &row.forward,
                &row.reverse_primer,
                &row.forward_tail,
                &row.reverse_tail,
                &self.methods.overlaps,
                &self.methods.comparisons,
                &self.methods.occurrences,
            ] {
                if text.len() > 65536 {
                    return Err("Assembly option field exceeds64KiB".into());
                }
            }
            let Some(captured) = slot.captured.as_deref() else {
                if slot.start != 0 || slot.end != 0 || slot.digest_geometry.is_some() {
                    return Err("An empty slot cannot retain a source interval".into());
                }
                continue;
            };
            let source = &captured.source;
            let molecule = &captured.molecule;
            let n = molecule.sequence.len();
            bases = bases
                .checked_add(n)
                .filter(|n| *n <= 8 << 20)
                .ok_or("Assembly pinned source bases exceed8MiB")?;
            if !token(&source.binding.document_id, 512)
                || !token(&source.binding.event_id, 512)
                || !hex(&source.binding.revision)
                || source.alphabet != Alphabet::Dna
                || source.alphabet != molecule.sequence.alphabet()
                || source.topology != molecule.topology
                || source.length != n
                || source.name != molecule.name
                || n == 0
                || slot.start >= n
                || slot.end <= slot.start
                || slot
                    .end
                    .checked_sub(slot.start)
                    .is_none_or(|length| length > n)
                || (molecule.topology == Topology::Linear && slot.end > n)
            {
                return Err("Invalid pinned DNA source or one-turn interval".into());
            }
            if let Some(geometry) = &slot.digest_geometry {
                digest_interval(source, geometry, false)?;
                digest_interval(source, geometry, true)?;
            }
            if molecule.revision() != source.binding.revision {
                return Err("Pinned source content differs from its captured revision".into());
            }
        }
        Ok(())
    }
    fn sync_rows(&mut self) {
        self.methods.circular = self.circular;
        self.methods
            .rows
            .resize_with(self.slots.len(), OverlapInput::default);
        for (i, (slot, row)) in self.slots.iter().zip(&mut self.methods.rows).enumerate() {
            row.document = i;
            row.reverse = slot.reverse;
            row.interval_start = slot.start.saturating_add(1).to_string();
            row.interval_end = if slot.captured.is_some() {
                slot.end.to_string()
            } else {
                String::new()
            };
        }
    }
    fn identity(&self) -> IntentIdentity {
        IntentIdentity {
            name: self.name.clone(),
            circular: self.circular,
            slots: self
                .slots
                .iter()
                .map(|s| {
                    (
                        s.label.clone(),
                        s.captured.as_ref().map(|c| c.source.clone()),
                        s.start,
                        s.end,
                        s.reverse,
                        s.digest_geometry.clone(),
                    )
                })
                .collect(),
        }
    }
}
pub fn encode(draft: &Draft) -> Result<Vec<u8>, String> {
    draft.validate()?;
    let mut writer = Bounded { bytes: Vec::new() };
    serde_json::to_writer(&mut writer, draft).map_err(|e| e.to_string())?;
    Ok(writer.bytes)
}
pub fn decode(bytes: &[u8]) -> Result<Draft, String> {
    if bytes.len() > MAX_BYTES {
        return Err("Assembly draft exceeds32MiB".into());
    }
    let draft: Draft = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
    draft.validate()?;
    Ok(draft)
}
#[derive(Clone, PartialEq, Eq)]
struct IntentIdentity {
    name: String,
    circular: bool,
    slots: Vec<(
        String,
        Option<SourceInfo>,
        usize,
        usize,
        bool,
        Option<DigestGeometry>,
    )>,
}
#[derive(Clone)]
pub struct SelectionNavigation {
    pub source: EventBinding,
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
}
#[derive(Clone)]
pub struct ReviewRequest {
    pub candidate: Candidate,
    pub sources: Vec<EventBinding>,
    pub intent: assembly_intent::Report,
}
struct Frozen {
    identity: IntentIdentity,
    epoch: u64,
    intent: assembly_intent::Report,
    sources: Vec<Arc<CapturedSource>>,
    target_revision: String,
}

fn inputs_request(
    draft: &Draft,
) -> Result<
    (
        Vec<PlanInput>,
        assembly_intent::Request,
        Vec<Arc<CapturedSource>>,
    ),
    String,
> {
    draft.validate()?;
    let captures = draft
        .slots
        .iter()
        .map(|s| {
            s.captured.clone().ok_or_else(|| {
                format!(
                    "Set the {} fragment before reviewing the intended construct.",
                    s.label
                )
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let ids = (0..captures.len())
        .map(|i| format!("source-{}", i + 1))
        .collect::<Vec<_>>();
    let request = assembly_intent::Request {
        name: draft.name.clone(),
        circular: draft.circular,
        fragments: draft
            .slots
            .iter()
            .enumerate()
            .map(|(i, s)| collagene_core::overlap_primer_design::Fragment {
                role: format!("fragment-{}", i + 1),
                input_id: ids[i].clone(),
                start: s.start,
                end: s.end,
                reverse: s.reverse,
            })
            .collect(),
    };
    let borrowed = ids
        .iter()
        .zip(&captures)
        .map(|(id, c)| (id.as_str(), &c.molecule))
        .collect::<Vec<_>>();
    assembly_intent::preflight(&borrowed, &request).map_err(|e| e.to_string())?;
    let inputs = ids
        .into_iter()
        .zip(&captures)
        .map(|(id, c)| PlanInput {
            id,
            revision: c.source.binding.revision.clone(),
            molecule: c.molecule.clone(),
        })
        .collect();
    Ok((inputs, request, captures))
}

#[component]
pub fn AssemblyWorkspace(
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    sources: Memo<Vec<SourceInfo>>,
    active_source: Memo<SourceInfo>,
    selected_region: Callback<(), Option<(usize, usize)>>,
    get_source: Callback<EventBinding, Result<Molecule, String>>,
    on_navigate: Callback<SelectionNavigation>,
    on_review: Callback<ReviewRequest, Result<(), String>>,
    on_finalize: Callback<ReviewRequest, Result<(), String>>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
    primer_sources: Memo<Vec<SourceInfo>>,
    list_primers: Callback<EventBinding, Result<primer_reuse::List, String>>,
    resolve_primer: Callback<ResolveRequest, Result<primer_reuse::Resolved, String>>,
) -> impl IntoView {
    let identity = Memo::new(move |_| draft.with(Draft::identity));
    let methods = RwSignal::new(draft.with_untracked(|d| d.methods.clone()));
    Effect::new(move |_| {
        let value = draft.with(|d| d.methods.clone());
        if methods.get_untracked() != value {
            methods.set(value);
        }
    });
    Effect::new(move |_| {
        let value = methods.get();
        if draft.with_untracked(|d| d.methods != value) {
            draft.update(|d| d.methods = value);
        }
    });
    let frozen = RwSignal::new(None::<Arc<Frozen>>);
    let busy = RwSignal::new(false);
    let status = RwSignal::new(String::new());
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
    });
    Effect::new(move |_| {
        let _ = (identity.get(), epoch.get());
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        frozen.set(None);
        busy.set(false);
        status.set(String::new());
    });
    let pin = Callback::new(move |kind: u8| {
        let source = active_source.get_untracked();
        let n = source.length;
        if source.alphabet != Alphabet::Dna || n == 0 {
            status.set("Choose a nonempty DNA source in the sequence viewer.".into());
            return;
        }
        let (start, end) = if kind == 0 {
            (0, n)
        } else {
            let Some((a, b)) = selected_region.run(()) else {
                status
                    .set("Select a contiguous source region in the sequence viewer first.".into());
                return;
            };
            if a >= b || b > n {
                status.set("Selection is outside the active source. Select it again.".into());
                return;
            }
            if kind == 2 {
                if source.topology != Topology::Circular {
                    status.set("Keeping outside a selection requires a circular source; choose an explicit region for linear DNA.".into());
                    return;
                }
                let length = n - (b - a);
                if length == 0 {
                    status.set("The outside of a whole-source selection is empty.".into());
                    return;
                }
                let start = b % n;
                let Some(end) = start.checked_add(length) else {
                    status.set("Selection range overflow.".into());
                    return;
                };
                (start, end)
            } else {
                (a, b)
            }
        };
        let molecule = match get_source.run(source.binding.clone()) {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if active_source.get_untracked().binding != source.binding {
            status.set("The active source changed before capture. Select the region again.".into());
            return;
        }
        let mut next = draft.get_untracked();
        let index = next.active_slot;
        if index >= next.slots.len() {
            status.set("Choose an assembly slot first.".into());
            return;
        }
        next.slots[index].captured = Some(Arc::new(CapturedSource { source, molecule }));
        next.slots[index].start = start;
        next.slots[index].end = end;
        next.slots[index].digest_geometry = None;
        let orientation = next.slots[index].reverse;
        next.methods.rows[index] = OverlapInput {
            preparation: OverlapPreparation::ExactPcr,
            reverse: orientation,
            ..Default::default()
        };
        next.sync_rows();
        match next.validate() {
            Ok(()) => {
                draft.set(next);
                status.set("Source region pinned to this slot. Review the intended construct, then choose or design its PCR primers/preparation.".into());
            }
            Err(e) => status.set(e),
        }
    });
    let review_intent = Callback::new(move |()| {
        let declared = draft.get_untracked();
        let state = identity.get_untracked();
        let bound_epoch = epoch.get_untracked();
        let (inputs, request, captures) = match inputs_request(&declared) {
            Ok(v) => v,
            Err(e) => {
                status.set(e);
                return;
            }
        };
        let revision = match assembly_intent::intent_revision(&inputs, &request) {
            Ok(v) => v,
            Err(e) => {
                status.set(e.to_string());
                return;
            }
        };
        let current = match CachedWorker::new() {
            Ok(v) => Rc::new(v),
            Err(e) => {
                status.set(e);
                return;
            }
        };
        if let Some(old) = worker.get_value() {
            old.close();
        }
        worker.set_value(Some(current.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        busy.set(true);
        frozen.set(None);
        status.set("Building the intended reference from your pinned source regions…".into());
        let alive = alive.clone();
        spawn_local(async move {
            let response = current
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job: AnalysisJob::BuildAssemblyIntent { inputs, request },
                })
                .await;
            current.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || identity.get_untracked() != state
                || epoch.get_untracked() != bound_epoch
            {
                return;
            }
            worker.set_value(None);
            busy.set(false);
            match response {
                Ok(response) => match response.result {
                    AnalysisResult::AssemblyIntent(intent) => {
                        match intent.validate_binding(&revision) {
                            Ok(()) => {
                                let target_revision = intent.target.revision();
                                frozen.set(Some(Arc::new(Frozen {
                                    identity: state,
                                    epoch: bound_epoch,
                                    intent,
                                    sources: captures,
                                    target_revision,
                                })));
                                status.set("Intended construct ready for review. It describes the sequence you want; no cloning method has been verified.".into());
                            }
                            Err(e) => status.set(e.to_string()),
                        }
                    }
                    AnalysisResult::Error(error) => status.set(error),
                    _ => status.set("Unexpected intended-construct response.".into()),
                },
                Err(error) => status.set(error),
            }
        });
    });
    let cancel = move |_| {
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        busy.set(false);
        frozen.set(None);
        status.set("Intended-construct review cancelled.".into());
    };
    view! {<section class="assembly-workspace"><header class="assembly-workspace-header"><h2>"Assembly workspace"</h2><p class="helper">"Choose a slot, select its source in the viewer above, then set the fragment. Switching source tabs preserves every pinned fragment. Preparation is reviewed against a separate intended construct."</p></header>
        <div class="assembly-draft-settings"><label>"Assembly name"<input maxlength="256" aria-label="Assembly draft name" prop:value=move||draft.with(|d|d.name.clone()) on:input=move|e|draft.update(|d|d.name=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Intended assembly circular" prop:checked=move||draft.with(|d|d.circular) on:change=move|e|draft.update(|d|{d.circular=event_target_checked(&e);d.sync_rows();})/>"Circular intended construct"</label></div>
        <h3>"1. Set original source fragments"</h3><div class="assembly-slot-list" role="group" aria-label="Ordered assembly fragment slots">
        <For each=move||draft.with(|d|(0..d.slots.len()).collect::<Vec<_>>()) key=|i|*i children=move|i|view!{<article class="assembly-slot" class:active=move||draft.with(|d|d.active_slot==i)><button class="button secondary" aria-pressed=move||draft.with(|d|d.active_slot==i) on:click=move|_|draft.update(|d|d.active_slot=i)>{move||draft.with(|d|format!("{} · {}",i+1,d.slots.get(i).map(|s|s.label.as_str()).unwrap_or("slot")))}</button><label>"Slot label"<input maxlength="80" aria-label=format!("Assembly slot {} label",i+1) prop:value=move||draft.with(|d|d.slots.get(i).map(|s|s.label.clone()).unwrap_or_default()) on:input=move|e|draft.update(|d|if let Some(s)=d.slots.get_mut(i){s.label=event_target_value(&e)})/></label>
            <p>{move||draft.with(|d|d.slots.get(i).map(|s|match s.captured.as_deref(){Some(c)=>format!("{} · {}–{} (1-based inclusive, unrolled) · {} bp · {}",c.source.name,s.start+1,s.end,s.end-s.start,if s.reverse{"reverse"}else{"forward"}),None=>"No source fragment set".into()}).unwrap_or_default())}</p>
            <details class="assembly-source-details"><summary>"Source details"</summary><p class="helper">{move||draft.with(|d|d.slots.get(i).and_then(|s|s.captured.as_ref()).map(|c|format!("Pinned source event {} · revision {}",c.source.binding.event_id,c.source.binding.revision)).unwrap_or_default())}</p></details>
            <p class="helper">{move||{let captured=draft.with(|d|d.slots.get(i).and_then(|s|s.captured.clone()));captured.map(|c|sources.with(|current|if current.iter().any(|s|s.binding==c.source.binding){"Current source version is pinned"}else{"Using the captured source version; current tabs do not replace it"})).unwrap_or("")}}</p>
            <label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse assembly slot {}",i+1) prop:checked=move||draft.with(|d|d.slots.get(i).is_some_and(|s|s.reverse)) on:change=move|e|draft.update(|d|{let reversed=event_target_checked(&e);if i<d.slots.len(){let adjusted=d.slots[i].captured.as_ref().zip(d.slots[i].digest_geometry.as_ref()).filter(|_|d.methods.rows[i].preparation==OverlapPreparation::CatalogDigest&&d.methods.rows[i].repair.is_none()).map(|(capture,geometry)|digest_interval(&capture.source,geometry,reversed));match adjusted{Some(Ok((start,end)))=>{d.slots[i].start=start;d.slots[i].end=end;},Some(Err(error))=>{status.set(error);return;},None=>{}}d.slots[i].reverse=reversed;}d.sync_rows();})/>"Reverse fragment in construct"</label>
            <div class="method-row-actions"><button class="text-tool" on:click=move|_|{if let Some(nav)=draft.with_untracked(|d|d.slots.get(i).and_then(|s|s.captured.as_ref().map(|c|SelectionNavigation{source:c.source.binding.clone(),start:s.start,end:s.end,reverse:s.reverse}))){on_navigate.run(nav);}}>"View source region"</button><button class="text-tool" disabled=i==0 on:click=move|_|draft.update(|d|{if i>0&&i<d.slots.len(){d.slots.swap(i,i-1);d.methods.rows.swap(i,i-1);d.active_slot=i-1;d.sync_rows();}})>"Move earlier"</button><button class="text-tool" disabled=move||draft.with(|d|i+1>=d.slots.len()) on:click=move|_|draft.update(|d|{if i+1<d.slots.len(){d.slots.swap(i,i+1);d.methods.rows.swap(i,i+1);d.active_slot=i+1;d.sync_rows();}})>"Move later"</button><button class="text-tool" disabled=move||draft.with(|d|d.slots.len()<=1) on:click=move|_|draft.update(|d|{if d.slots.len()>1&&i<d.slots.len(){d.slots.remove(i);d.methods.rows.remove(i);d.active_slot=d.active_slot.min(d.slots.len()-1);d.sync_rows();}})>"Remove slot"</button></div>
        </article>}/></div>
        <button class="button secondary" disabled=move||draft.with(|d|d.slots.len()>=12) on:click=move|_|draft.update(|d|{if d.slots.len()<12{d.slots.push(Slot::empty(&format!("Fragment {}",d.slots.len()+1)));d.methods.rows.push(OverlapInput{preparation:OverlapPreparation::ExactPcr,..Default::default()});d.active_slot=d.slots.len()-1;d.sync_rows();}})>"Add fragment slot"</button>
        <section class="assembly-set-fragment"><p>{move||format!("Active viewer source: {} · set {}",active_source.get().name,draft.with(|d|d.slots.get(d.active_slot).map(|s|s.label.as_str()).unwrap_or("selected slot").to_owned()))}</p><button class="button primary" on:click=move|_|pin.run(1)>"Set from sequence selection"</button><button class="button secondary" on:click=move|_|pin.run(0)>"Use whole source"</button><button class="button secondary" on:click=move|_|pin.run(2)>"Keep outside selection"</button><p class="helper">"Keep outside selection retains the complementary path on circular DNA, including an origin crossing. This sets intended source coordinates; it does not digest, amplify or prove a physical fragment."</p></section>
        <button class="button primary" disabled=move||busy.get() on:click=move|_|review_intent.run(())>"Review intended construct"</button><button class="text-tool" disabled=move||!busy.get() on:click=cancel>"Cancel intended construct"</button><p class="assembly-workspace-status" role="status">{move||status.get()}</p>
        {move||frozen.get().map(|bound|view!{<FrozenWorkspace bound draft methods epoch on_navigate on_review on_finalize on_products primer_sources list_primers resolve_primer/>})}
    </section>}
}

fn intended_window(molecule: &Molecule, boundary: usize) -> String {
    let dna = molecule.sequence.as_str().as_bytes();
    let n = dna.len();
    if n == 0 {
        return "Intended sequence unavailable".into();
    }
    let (left, right) = if molecule.topology == Topology::Circular {
        let b = boundary % n;
        let count = 20.min(n);
        (
            (0..count)
                .map(|i| dna[(b + n - count + i) % n] as char)
                .collect::<String>(),
            (0..count)
                .map(|i| dna[(b + i) % n] as char)
                .collect::<String>(),
        )
    } else {
        let b = boundary.min(n);
        (
            String::from_utf8_lossy(&dna[b.saturating_sub(20)..b]).into_owned(),
            String::from_utf8_lossy(&dna[b..b.saturating_add(20).min(n)]).into_owned(),
        )
    };
    format!("{left}|{right}")
}
#[component]
fn FrozenWorkspace(
    bound: Arc<Frozen>,
    draft: RwSignal<Draft>,
    methods: RwSignal<MethodsState>,
    epoch: RwSignal<u64>,
    on_navigate: Callback<SelectionNavigation>,
    on_review: Callback<ReviewRequest, Result<(), String>>,
    on_finalize: Callback<ReviewRequest, Result<(), String>>,
    on_products: Callback<Vec<Molecule>, Result<(), String>>,
    primer_sources: Memo<Vec<SourceInfo>>,
    list_primers: Callback<EventBinding, Result<primer_reuse::List, String>>,
    resolve_primer: Callback<ResolveRequest, Result<primer_reuse::Resolved, String>>,
) -> impl IntoView {
    let bound = StoredValue::new(bound);
    let junction = RwSignal::new(0usize);
    let status = RwSignal::new(String::new());
    let target = Memo::new(move |_| bound.with_value(|f| f.intent.target.clone()));
    let target_info = Memo::new(move |_| {
        bound.with_value(|f| TargetContext::Intent {
            context_revision: f.intent.context_revision.clone(),
            revision: f.target_revision.clone(),
        })
    });
    let query_options = Memo::new(move |_| {
        bound.with_value(|f| {
            f.sources
                .iter()
                .enumerate()
                .map(|(i, c)| (i, format!("{} · {}", f.identity.slots[i].0, c.source.name)))
                .collect::<Vec<_>>()
        })
    });
    let get_query = Callback::new(move |i: usize| {
        bound.with_value(|f| {
            f.sources.get(i).and_then(|c| {
                collagene_core::model::bounds::model_preflight(&[&c.molecule])
                    .ok()
                    .map(|_| c.molecule.clone())
            })
        })
    });
    let source_info = Callback::new(move |i: usize| {
        bound.with_value(|f| f.sources.get(i).map(|c| c.source.clone()))
    });
    let prepare_review = Callback::new(move |candidate: Candidate| {
        let outcome = bound.with_value(|f| {
            if epoch.get_untracked() != f.epoch
                || draft.with_untracked(Draft::identity) != f.identity
            {
                return Err("Source slots changed; review the intended construct again.".into());
            }
            if candidate.plan.target_revision != f.target_revision
                || candidate.plan.inputs.len() != f.sources.len()
                || candidate.plan.inputs.iter().enumerate().any(|(i, p)| {
                    p.id != format!("source-{}", i + 1)
                        || p.revision != f.sources[i].source.binding.revision
                })
            {
                return Err(
                    "Candidate inputs differ from the pinned source order; no candidate opened."
                        .into(),
                );
            }
            Ok(ReviewRequest {
                candidate,
                sources: f.sources.iter().map(|c| c.source.binding.clone()).collect(),
                intent: f.intent.clone(),
            })
        });
        outcome
    });
    let review = Callback::new(move |candidate: Candidate| {
        let result = prepare_review
            .run(candidate)
            .and_then(|request| on_review.run(request));
        match result { Ok(()) => status.set("Candidate opened for fresh verification against the intended reference. No experimental success is established.".into()), Err(error) => status.set(error) }
    });
    let finalize = Callback::new(move |candidate: Candidate| {
        prepare_review
            .run(candidate)
            .and_then(|request| on_finalize.run(request))
    });
    let navigate = Callback::new(move |i: usize| {
        if let Some(nav) = bound.with_value(|f| {
            f.intent.request.fragments.get(i).and_then(|region| {
                f.sources.get(i).map(|source| SelectionNavigation {
                    source: source.source.binding.clone(),
                    start: region.start,
                    end: region.end,
                    reverse: region.reverse,
                })
            })
        }) {
            on_navigate.run(nav);
        }
    });
    view! {<section class="assembly-intended-review"><h3>"2. Review the intended construct"</h3><p>{bound.with_value(|f|format!("{} · {} bp · {}",f.intent.target.name,f.intent.target.sequence.len(),if f.intent.target.topology==Topology::Circular{"circular"}else{"linear"}))}</p><p class="helper">"This sequence is your declared design requirement, assembled only from the chosen source regions. It is separate from any simulated cloning outcome. Editing a slot clears this review."</p>
        <div class="assembly-fragment-ribbon" role="group" aria-label="Intended construct fragments">{bound.with_value(|f|f.intent.fragments.iter().enumerate().map(|(i,fragment)|view!{<button class="assembly-ribbon-fragment" on:click=move|_|navigate.run(i)>{format!("{}{} · {} bp · {}–{}",f.identity.slots[i].0,if fragment.reverse{" ←"}else{" →"},fragment.output_end-fragment.output_start,fragment.output_start+1,fragment.output_end)}</button>}).collect_view())}</div>
        {bound.with_value(|f|{let count=f.intent.fragments.len().saturating_sub(1)+usize::from(f.intent.target.topology==Topology::Circular);if count==0{return view!{<p>"One linear fragment has no internal junction."</p>}.into_any();}view!{<div class="assembly-intended-junction">{move||bound.with_value(|f|{let i=junction.get().min(count-1);let next=(i+1)%f.sources.len();let boundary=f.intent.fragments[i].output_end;view!{<h4>{format!("Intended junction {} · {} → {}{}",i+1,f.identity.slots[i].0,f.identity.slots[next].0,if i+1==f.sources.len(){" · circle closure"}else{""})}</h4><pre>{intended_window(&f.intent.target,boundary)}</pre><p class="helper">{format!("| marks intended product boundary {} (zero-based). These letters are not a digest/ligation/homology validation.",if boundary==f.intent.target.sequence.len()&&f.intent.target.topology==Topology::Circular{0}else{boundary})}</p><button class="text-tool" on:click=move|_|navigate.run(i)>"View upstream source"</button><button class="text-tool" on:click=move|_|navigate.run(next)>"View downstream source"</button>}})}<div class="method-row-actions"><button class="text-tool" disabled=move||junction.get()==0 on:click=move|_|junction.update(|i|*i=i.saturating_sub(1))>"Previous intended junction"</button><button class="text-tool" disabled=move||{junction.get()+1>=count} on:click=move|_|junction.update(|i|*i=(*i+1).min(count-1))>"Next intended junction"</button></div></div>}.into_any()})}
        {bound.with_value(|f|f.intent.warnings.iter().map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view())}
        <h3>"3. Prepare fragments and verify the method"</h3><p class="helper">"PCR rows begin with no primer sequences. Suggest and explicitly apply primers, or choose another declared preparation. The pinned original sources remain the inputs; the intended reference is never used as physical source material."</p>
        <CloningMethodsPanel molecule=target query_options get_query source_info target_info selected_region=Callback::new(|()|None) config=methods epoch on_candidate=review on_finalize=Some(finalize) on_products primer_sources list_primers resolve_primer fixed_sources=true/>
        <p role="status">{move||status.get()}</p>
    </section>}
}

fn digest_interval(
    source: &SourceInfo,
    geometry: &DigestGeometry,
    reverse: bool,
) -> Result<(usize, usize), String> {
    let n = source.length;
    if n == 0 {
        return Err("Digest source is empty".into());
    }
    let (start, length) = if reverse {
        geometry
            .bottom_start
            .checked_add(
                i64::try_from(geometry.bottom_length)
                    .map_err(|_| "Bottom interval length overflow")?,
            )
            .ok_or("Bottom interval coordinate overflow")?;
        if geometry.bottom_length == 0 || geometry.bottom_length > n {
            return Err("Digest bottom-strand interval is empty or exceeds one source turn".into());
        }
        let start = if source.topology == Topology::Circular {
            let size =
                i64::try_from(n).map_err(|_| "Source length exceeds signed coordinate range")?;
            usize::try_from(geometry.bottom_start.rem_euclid(size))
                .map_err(|_| "Invalid bottom-strand coordinate")?
        } else {
            usize::try_from(geometry.bottom_start)
                .map_err(|_| "Linear digest bottom interval starts outside its source")?
        };
        (start, geometry.bottom_length)
    } else {
        (geometry.top_start, geometry.top_length)
    };
    if start >= n || length == 0 || length > n {
        return Err("Digest interval is outside its original source".into());
    }
    let end = start
        .checked_add(length)
        .ok_or("Digest interval overflow")?;
    if source.topology == Topology::Linear && end > n {
        return Err("Linear digest interval extends beyond its source".into());
    }
    Ok((start, end))
}
/// Pin a fresh private digest selection, retaining its original source and full
/// catalog identity. This checks bindings/bounds, not cleavage; the compiled
/// CatalogDigest operation must reconstruct the actual fragment in the worker.
/// The host must check payload epoch and actual retained source-event identity.
/// Error leaves the original draft untouched.
pub fn set_digest_fragment(
    draft: &mut Draft,
    fragment: crate::digest::AssemblyFragment,
) -> Result<(), String> {
    use collagene_core::{catalog_digest, restriction_history};
    draft.validate()?;
    if !token(&fragment.source.binding.document_id, 512)
        || !token(&fragment.source.binding.event_id, 512)
        || !hex(&fragment.source.binding.revision)
    {
        return Err("Invalid original-source event binding".into());
    }
    let index = draft.active_slot;
    let catalog = catalog_digest::builtin_catalog().map_err(|e| e.to_string())?;
    let role =
        restriction_history::source_role(&fragment.source.binding).map_err(|e| e.to_string())?;
    catalog_digest::source_preflight(
        &catalog,
        &fragment.selection,
        &role,
        &fragment.molecule,
        &fragment.options,
    )
    .map_err(|e| e.to_string())?;
    let n = fragment.molecule.sequence.len();
    if fragment.source.length != n
        || fragment.source.name != fragment.molecule.name
        || fragment.source.topology != fragment.molecule.topology
        || fragment.source.alphabet != fragment.molecule.sequence.alphabet()
        || fragment.source.alphabet != Alphabet::Dna
        || fragment.molecule.revision() != fragment.source.binding.revision
    {
        return Err("Digest original source differs from its captured document revision".into());
    }
    if fragment.fragment_index >= fragment.options.max_fragments
        || !hex(&fragment.batch_revision)
        || catalog_digest::batch_revision(
            &catalog,
            &fragment.selection,
            &role,
            &fragment.molecule,
            &fragment.options,
        )
        .map_err(|e| e.to_string())?
            != fragment.batch_revision
    {
        return Err("Digest selection/source/options changed; run the digest again".into());
    }
    let geometry = DigestGeometry {
        top_start: fragment.top_start,
        top_length: fragment.top_length,
        bottom_start: fragment.bottom_start,
        bottom_length: fragment.bottom_length,
    };
    digest_interval(&fragment.source, &geometry, false)?;
    digest_interval(&fragment.source, &geometry, true)?;
    let (start, end) = digest_interval(&fragment.source, &geometry, draft.slots[index].reverse)?;
    let mut borrowed = draft
        .slots
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != index)
        .filter_map(|(_, slot)| slot.captured.as_ref().map(|c| &c.molecule))
        .collect::<Vec<_>>();
    borrowed.push(&fragment.molecule);
    collagene_core::model::bounds::model_preflight(&borrowed).map_err(|e| e.to_string())?;
    let mut bases = 0usize;
    for source in &borrowed {
        bases = bases
            .checked_add(source.sequence.len())
            .filter(|n| *n <= 8 << 20)
            .ok_or("Pinned assembly sources exceed8MiB")?;
    }
    let display_index = fragment
        .fragment_index
        .checked_add(1)
        .ok_or("Digest fragment index overflow")?;
    let mut next = draft.clone();
    next.slots[index].captured = Some(Arc::new(CapturedSource {
        source: fragment.source,
        molecule: fragment.molecule.as_ref().clone(),
    }));
    next.slots[index].start = start;
    next.slots[index].end = end;
    next.slots[index].digest_geometry = Some(geometry);
    next.methods.rows[index] = OverlapInput {
        preparation: OverlapPreparation::CatalogDigest,
        catalog_selection: Some(fragment.selection),
        fragment: display_index.to_string(),
        repair: None,
        reverse: next.slots[index].reverse,
        ..Default::default()
    };
    next.sync_rows();
    next.validate()?;
    *draft = next;
    Ok(())
}
