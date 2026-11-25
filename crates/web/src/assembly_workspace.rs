//! Persistent declared assembly slots. Frozen intent/scientific evidence is
//! component-private and is never restored from a serialized draft.
use crate::{
    assembly_review::IntendedReference,
    cloning_methods::{
        CloningMethodsPanel, MethodsState, Mode, OverlapInput, OverlapPreparation, TargetContext,
    },
    digest::SourceInfo,
    primer_picker::ResolveRequest,
    worker::CachedWorker,
};
use collagene_core::{
    Alphabet, Molecule, Topology, assembly_intent, assembly_junction_edit,
    golden_gate_junction_edit, primer_reuse,
    protocol::{AnalysisJob, AnalysisResult, Request},
    run_lineage::EventBinding,
    selective_junction_edit,
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
    #[serde(default, deserialize_with = "deserialize_junction_edits")]
    pub junction_edits: Vec<assembly_junction_edit::Edit>,
}
fn deserialize_junction_edits<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<Vec<assembly_junction_edit::Edit>, D::Error> {
    struct Visitor;
    impl<'de> serde::de::Visitor<'de> for Visitor {
        type Value = Vec<assembly_junction_edit::Edit>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("at most 12 explicit junction edits")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            if seq.size_hint().is_some_and(|n| n > 12) {
                return Err(serde::de::Error::custom("at most 12 junction edits"));
            }
            let mut values = Vec::new();
            while values.len() < 12 {
                match seq.next_element()? {
                    Some(value) => values.push(value),
                    None => return Ok(values),
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom("at most 12 junction edits"));
            }
            Ok(values)
        }
    }
    deserializer.deserialize_seq(Visitor)
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
            junction_edits: Vec::new(),
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
        if let Some(pool) = &self.methods.golden_gate_pool {
            pool.validate()?;
            if pool.choices.len() != self.slots.len() {
                return Err("Reaction substrate choices must match the assembly slots.".into());
            }
        }
        let junction_count = self.slots.len().saturating_sub(1) + usize::from(self.circular);
        if self.junction_edits.len() > 12
            || (!self.junction_edits.is_empty()
                && !matches!(self.methods.mode, Mode::Overlap | Mode::GoldenGate))
        {
            return Err(
                "Junction edits require overlap or Golden Gate assembly and at most 12 edits"
                    .into(),
            );
        }
        let mut junctions = std::collections::BTreeSet::new();
        for edit in &self.junction_edits {
            if edit.junction >= junction_count
                || !junctions.insert(edit.junction)
                || edit.inserted.len() > 200
            {
                return Err(
                    "Invalid or duplicate junction edit; inserted text is limited to 200 bytes"
                        .into(),
                );
            }
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
            junction_edits: self.junction_edits.clone(),
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
fn same_declaration(current: &Draft, applied: &Draft) -> bool {
    // The active slot is view state, not a change to the declared design.
    current.version == applied.version
        && current.name == applied.name
        && current.circular == applied.circular
        && current.slots == applied.slots
        && current.methods == applied.methods
        && current.junction_edits == applied.junction_edits
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
    junction_edits: Vec<assembly_junction_edit::Edit>,
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
    pub intent: IntendedReference,
}
fn context(intent: &IntendedReference) -> &str {
    match intent {
        IntendedReference::Regions(i) => &i.context_revision,
        IntendedReference::Junctions(i) => &i.context_revision,
    }
}
fn warnings(intent: &IntendedReference) -> &[String] {
    match intent {
        IntendedReference::Regions(i) => &i.warnings,
        IntendedReference::Junctions(i) => &i.warnings,
    }
}
struct Frozen {
    identity: IntentIdentity,
    epoch: u64,
    intent: IntendedReference,
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
    if !draft.junction_edits.is_empty() {
        assembly_junction_edit::intent_preflight(&borrowed, &request, &draft.junction_edits)
            .map_err(|e| e.to_string())?;
    }
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
    let undo_junction = RwSignal::new(None::<(Draft, Draft)>);
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
        let (current, current_epoch) = (identity.get(), epoch.get());
        if frozen.with_untracked(|value| {
            value
                .as_ref()
                .is_some_and(|f| f.identity == current && f.epoch == current_epoch)
        }) {
            return;
        }
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
        next.junction_edits.clear();
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
        let revision = match if declared.junction_edits.is_empty() {
            assembly_intent::intent_revision(&inputs, &request)
        } else {
            assembly_junction_edit::intent_revision(&inputs, &request, &declared.junction_edits)
        } {
            Ok(value) => value,
            Err(error) => {
                status.set(error.to_string());
                return;
            }
        };
        let job = if declared.junction_edits.is_empty() {
            AnalysisJob::BuildAssemblyIntent { inputs, request }
        } else {
            AnalysisJob::BuildEditedAssemblyIntent {
                inputs,
                base: request,
                edits: declared.junction_edits.clone(),
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
                    job,
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
                Ok(response) if response.id == token && response.revision == revision => {
                    let result = match response.result {
                        AnalysisResult::AssemblyIntent(intent) => intent
                            .validate_binding(&revision)
                            .map(|()| IntendedReference::Regions(intent)),
                        AnalysisResult::EditedAssemblyIntent(intent) => intent
                            .validate_binding(&revision)
                            .map(|()| IntendedReference::Junctions(intent)),
                        AnalysisResult::Error(error) => {
                            status.set(error);
                            return;
                        }
                        _ => {
                            status.set("Unexpected intended-construct response.".into());
                            return;
                        }
                    };
                    match result {
                        Ok(intent) => {
                            let target_revision = intent.target().revision();
                            frozen.set(Some(Arc::new(Frozen {
                                identity: state,
                                epoch: bound_epoch,
                                intent,
                                sources: captures,
                                target_revision,
                            })));
                            status.set("Intended construct ready for review. No cloning method has been verified.".into());
                        }
                        Err(error) => status.set(error.to_string()),
                    }
                }
                Ok(_) => status.set("Intended-construct response identity changed.".into()),
                Err(error) => status.set(error),
            }
        });
    });
    let apply_junction = Callback::new(
        move |proposal: Arc<JunctionProposal>| -> Result<(), String> {
            if epoch.get_untracked() != proposal.epoch || draft.get_untracked() != proposal.before {
                return Err(
                    "Source slots or settings changed. Preview the junction edit again.".into(),
                );
            }
            let mut next = proposal.before.clone();
            if proposal.report.len() != next.slots.len() {
                return Err(
                    "Proposed PCR preparation does not cover every pinned source slot.".into(),
                );
            }
            for i in 0..proposal.report.len() {
                let fragment = proposal
                    .report
                    .fragment(i)
                    .ok_or("Proposed fragment is missing.")?;
                if fragment.input_id != format!("source-{}", i + 1)
                    || fragment.role != format!("fragment-{}", i + 1)
                    || fragment.reverse != next.slots[i].reverse
                {
                    return Err("Proposed source order or orientation changed.".into());
                }
                match proposal.report.retained(i) {
                    Some(true) => continue,
                    Some(false) => {}
                    None => return Err("Proposed row preservation choices are incomplete.".into()),
                }
                let (forward, reverse, selection, fragment_index) = match &fragment.preparation {
                    collagene_core::assembly_recipe::Preparation::ExactPcr { forward, reverse }
                        if matches!(&proposal.report, JunctionReport::Overlap(_)) =>
                    {
                        (forward, reverse, None, None)
                    }
                    collagene_core::assembly_recipe::Preparation::PcrCatalogDigest {
                        forward,
                        reverse,
                        selection,
                        fragment_index,
                        repair: None,
                    } if matches!(&proposal.report, JunctionReport::GoldenGate(_)) => (
                        forward,
                        reverse,
                        Some(selection.clone()),
                        Some(*fragment_index),
                    ),
                    _ => {
                        return Err(
                            "Selected row does not contain its proposed preparation chain.".into(),
                        );
                    }
                };
                let row = &mut next.methods.rows[i];
                row.preparation = if selection.is_some() {
                    OverlapPreparation::PcrCatalogDigest
                } else {
                    OverlapPreparation::ExactPcr
                };
                row.catalog_selection = selection;
                if let Some(index) = fragment_index {
                    row.fragment = (index + 1).to_string();
                }
                row.repair = None;
                row.treatment = crate::cloning_methods::Treatment::Keep;
                row.forward_name = forward.name.clone();
                row.forward = forward.annealing.as_str().into();
                row.forward_tail = forward.tail.clone();
                row.reverse_name = reverse.name.clone();
                row.reverse_primer = reverse.annealing.as_str().into();
                row.reverse_tail = reverse.tail.clone();
            }
            next.junction_edits = proposal.report.edits().to_vec();
            if let JunctionReport::Overlap(report) = &proposal.report {
                next.methods.overlaps = report
                    .recipe
                    .overlaps
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",");
            }
            next.sync_rows();
            next.validate()?;
            let identity = next.identity();
            let intent = IntendedReference::Junctions(proposal.report.intent().clone());
            let target_revision = intent.target().revision();
            undo_junction.set(Some((proposal.before.clone(), next.clone())));
            methods.set(next.methods.clone());
            draft.set(next);
            frozen.set(Some(Arc::new(Frozen {
                identity,
                epoch: proposal.epoch,
                intent,
                sources: proposal.sources.clone(),
                target_revision,
            })));
            status.set("Junction edits and selected PCR primers applied together; other preparations were preserved. Compute fresh method evidence against the revised intended reference.".into());
            Ok(())
        },
    );
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
        <div class="assembly-draft-settings"><label>"Assembly name"<input maxlength="256" aria-label="Assembly draft name" prop:value=move||draft.with(|d|d.name.clone()) on:input=move|e|draft.update(|d|d.name=event_target_value(&e))/></label><label class="checkbox-label"><input type="checkbox" aria-label="Intended assembly circular" prop:checked=move||draft.with(|d|d.circular) on:change=move|e|draft.update(|d|{d.circular=event_target_checked(&e);d.junction_edits.clear();d.sync_rows();})/>"Circular intended construct"</label></div>
        <h3>"1. Set original source fragments"</h3><div class="assembly-slot-list" role="group" aria-label="Ordered assembly fragment slots">
        <For each=move||draft.with(|d|(0..d.slots.len()).collect::<Vec<_>>()) key=|i|*i children=move|i|view!{<article class="assembly-slot" class:active=move||draft.with(|d|d.active_slot==i)><button class="button secondary" aria-pressed=move||draft.with(|d|d.active_slot==i) on:click=move|_|draft.update(|d|d.active_slot=i)>{move||draft.with(|d|format!("{} · {}",i+1,d.slots.get(i).map(|s|s.label.as_str()).unwrap_or("slot")))}</button><label>"Slot label"<input maxlength="80" aria-label=format!("Assembly slot {} label",i+1) prop:value=move||draft.with(|d|d.slots.get(i).map(|s|s.label.clone()).unwrap_or_default()) on:input=move|e|draft.update(|d|if let Some(s)=d.slots.get_mut(i){s.label=event_target_value(&e)})/></label>
            <p>{move||draft.with(|d|d.slots.get(i).map(|s|match s.captured.as_deref(){Some(c)=>format!("{} · {}–{} (1-based inclusive, unrolled) · {} bp · {}",c.source.name,s.start+1,s.end,s.end-s.start,if s.reverse{"reverse"}else{"forward"}),None=>"No source fragment set".into()}).unwrap_or_default())}</p>
            <details class="assembly-source-details"><summary>"Source details"</summary><p class="helper">{move||draft.with(|d|d.slots.get(i).and_then(|s|s.captured.as_ref()).map(|c|format!("Pinned source event {} · revision {}",c.source.binding.event_id,c.source.binding.revision)).unwrap_or_default())}</p></details>
            <p class="helper">{move||{let captured=draft.with(|d|d.slots.get(i).and_then(|s|s.captured.clone()));captured.map(|c|sources.with(|current|if current.iter().any(|s|s.binding==c.source.binding){"Current source version is pinned"}else{"Using the captured source version; current tabs do not replace it"})).unwrap_or("")}}</p>
            <label class="checkbox-label"><input type="checkbox" aria-label=format!("Reverse assembly slot {}",i+1) prop:checked=move||draft.with(|d|d.slots.get(i).is_some_and(|s|s.reverse)) on:change=move|e|draft.update(|d|{let reversed=event_target_checked(&e);if i<d.slots.len(){let adjusted=d.slots[i].captured.as_ref().zip(d.slots[i].digest_geometry.as_ref()).filter(|_|d.methods.rows[i].preparation==OverlapPreparation::CatalogDigest&&d.methods.rows[i].repair.is_none()).map(|(capture,geometry)|digest_interval(&capture.source,geometry,reversed));match adjusted{Some(Ok((start,end)))=>{d.slots[i].start=start;d.slots[i].end=end;},Some(Err(error))=>{status.set(error);return;},None=>{}}d.slots[i].reverse=reversed;d.junction_edits.clear();}d.sync_rows();})/>"Reverse fragment in construct"</label>
            <div class="method-row-actions"><button class="text-tool" on:click=move|_|{if let Some(nav)=draft.with_untracked(|d|d.slots.get(i).and_then(|s|s.captured.as_ref().map(|c|SelectionNavigation{source:c.source.binding.clone(),start:s.start,end:s.end,reverse:s.reverse}))){on_navigate.run(nav);}}>"View source region"</button><button class="text-tool" disabled=i==0 on:click=move|_|draft.update(|d|{if i>0&&i<d.slots.len(){d.slots.swap(i,i-1);d.methods.swap_rows(i,i-1);d.junction_edits.clear();d.active_slot=i-1;d.sync_rows();}})>"Move earlier"</button><button class="text-tool" disabled=move||draft.with(|d|i+1>=d.slots.len()) on:click=move|_|draft.update(|d|{if i+1<d.slots.len(){d.slots.swap(i,i+1);d.methods.swap_rows(i,i+1);d.junction_edits.clear();d.active_slot=i+1;d.sync_rows();}})>"Move later"</button><button class="text-tool" disabled=move||draft.with(|d|d.slots.len()<=1) on:click=move|_|draft.update(|d|{if d.slots.len()>1&&i<d.slots.len(){d.slots.remove(i);d.methods.remove_row(i);d.junction_edits.clear();d.active_slot=d.active_slot.min(d.slots.len()-1);d.sync_rows();}})>"Remove slot"</button></div>
        </article>}/></div>
        <button class="button secondary" disabled=move||draft.with(|d|d.slots.len()>=12) on:click=move|_|draft.update(|d|{if d.slots.len()<12{d.slots.push(Slot::empty(&format!("Fragment {}",d.slots.len()+1)));d.methods.push_row(OverlapInput{preparation:OverlapPreparation::ExactPcr,..Default::default()});d.active_slot=d.slots.len()-1;d.junction_edits.clear();d.sync_rows();}})>"Add fragment slot"</button>
        <section class="assembly-set-fragment"><p>{move||format!("Active viewer source: {} · set {}",active_source.get().name,draft.with(|d|d.slots.get(d.active_slot).map(|s|s.label.as_str()).unwrap_or("selected slot").to_owned()))}</p><button class="button primary" on:click=move|_|pin.run(1)>"Set from sequence selection"</button><button class="button secondary" on:click=move|_|pin.run(0)>"Use whole source"</button><button class="button secondary" on:click=move|_|pin.run(2)>"Keep outside selection"</button><p class="helper">"Keep outside selection retains the complementary path on circular DNA, including an origin crossing. This sets intended source coordinates; it does not digest, amplify or prove a physical fragment."</p></section>
        <button class="button primary" disabled=move||busy.get() on:click=move|_|review_intent.run(())>"Review intended construct"</button><button class="text-tool" disabled=move||!busy.get() on:click=cancel>"Cancel intended construct"</button><p class="assembly-workspace-status" role="status">{move||status.get()}</p>
        <JunctionEditor draft epoch on_apply=apply_junction/>
        <button class="text-tool" disabled=move||undo_junction.with(|undo|undo.as_ref().is_none_or(|(_,applied)|!same_declaration(&draft.get(),applied))) on:click=move|_|{
            if let Some((mut previous,applied))=undo_junction.get_untracked(){
                let current=draft.get_untracked();
                if !same_declaration(&current,&applied){status.set("Undo is available only before further draft changes. Your later changes were preserved.".into());return;}
                if current.active_slot<previous.slots.len(){previous.active_slot=current.active_slot;}
                methods.set(previous.methods.clone());draft.set(previous);frozen.set(None);undo_junction.set(None);status.set("Previous assembly draft restored. Review its intent and recompute evidence; no positive result was restored.".into());
            }
        }>"Undo applied junction edit"</button>
        <Show when=move||undo_junction.with(|undo|undo.as_ref().is_some_and(|(_,applied)|!same_declaration(&draft.get(),applied)))><p class="helper">"Undo available before further draft changes. Later source, primer or design edits are preserved."</p></Show>
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
    let target = Memo::new(move |_| bound.with_value(|f| f.intent.target().clone()));
    let target_info = Memo::new(move |_| {
        bound.with_value(|f| TargetContext::Intent {
            context_revision: context(&f.intent).to_string(),
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
            f.intent.fragments().get(i).and_then(|region| {
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
    let diagnostic_source = Callback::new(move |(document, row): (usize, usize)| {
        let nav = bound.with_value(|f| {
            if epoch.get_untracked() != f.epoch
                || draft.with_untracked(Draft::identity) != f.identity
                || document != row
            {
                return Err(
                    "The pinned source slots changed; review the intended construct again."
                        .to_string(),
                );
            }
            let region = f
                .intent
                .fragments()
                .get(row)
                .ok_or("The preparation row is unavailable.")?;
            let source = f
                .sources
                .get(row)
                .ok_or("The original source is unavailable.")?;
            Ok(SelectionNavigation {
                source: source.source.binding.clone(),
                start: region.start,
                end: region.end,
                reverse: region.reverse,
            })
        })?;
        draft.update(|d| d.active_slot = row);
        on_navigate.run(nav);
        Ok(())
    });
    view! {<section class="assembly-intended-review"><h3>"2. Review the intended construct"</h3><p>{bound.with_value(|f|format!("{} · {} bp · {}",f.intent.target().name,f.intent.target().sequence.len(),if f.intent.target().topology==Topology::Circular{"circular"}else{"linear"}))}</p><p class="helper">"This sequence is your declared design requirement, assembled from chosen source regions and any explicit junction edits. It is separate from any simulated cloning outcome. Editing a slot clears this review."</p>
        {bound.with_value(|f|match &f.intent{IntendedReference::Junctions(intent)=>Some(view!{<p class="helper">{format!("{} explicit junction edits are included. Ribbon output spans include each adjusted source core and its outgoing synthetic insert; source navigation selects only the adjusted original-source core.",intent.edits.len())}</p>}),_=>None})}
        <div class="assembly-fragment-ribbon" role="group" aria-label="Intended construct fragments">{bound.with_value(|f|f.intent.fragments().iter().enumerate().map(|(i,fragment)|view!{<button class="assembly-ribbon-fragment" on:click=move|_|navigate.run(i)>{format!("{}{} · {} bp · {}–{}",f.identity.slots[i].0,if fragment.reverse{" ←"}else{" →"},fragment.output_end-fragment.output_start,fragment.output_start+1,fragment.output_end)}</button>}).collect_view())}</div>
        {bound.with_value(|f|{let count=f.intent.fragments().len().saturating_sub(1)+usize::from(f.intent.target().topology==Topology::Circular);if count==0{return view!{<p>"One linear fragment has no internal junction."</p>}.into_any();}view!{<div class="assembly-intended-junction">{move||bound.with_value(|f|{let i=junction.get().min(count-1);let next=(i+1)%f.sources.len();let boundary=f.intent.fragments()[i].output_end;view!{<h4>{format!("Intended junction {} · {} → {}{}",i+1,f.identity.slots[i].0,f.identity.slots[next].0,if i+1==f.sources.len(){" · circle closure"}else{""})}</h4><pre>{intended_window(f.intent.target(),boundary)}</pre><p class="helper">{format!("| marks intended product boundary {} (zero-based). These letters are not a digest/ligation/homology validation.",if boundary==f.intent.target().sequence.len()&&f.intent.target().topology==Topology::Circular{0}else{boundary})}</p><button class="text-tool" on:click=move|_|navigate.run(i)>"View upstream source"</button><button class="text-tool" on:click=move|_|navigate.run(next)>"View downstream source"</button>}})}<div class="method-row-actions"><button class="text-tool" disabled=move||junction.get()==0 on:click=move|_|junction.update(|i|*i=i.saturating_sub(1))>"Previous intended junction"</button><button class="text-tool" disabled=move||{junction.get()+1>=count} on:click=move|_|junction.update(|i|*i=(*i+1).min(count-1))>"Next intended junction"</button></div></div>}.into_any()})}
        {bound.with_value(|f|warnings(&f.intent).iter().map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view())}
        <h3>"3. Prepare fragments and verify the method"</h3><p class="helper">"PCR rows begin with no primer sequences. Suggest and explicitly apply primers, or choose another declared preparation. The pinned original sources remain the inputs; the intended reference is never used as physical source material."</p>
        <CloningMethodsPanel molecule=target query_options get_query source_info target_info selected_region=Callback::new(|()|None) config=methods epoch on_candidate=review on_finalize=Some(finalize) on_products primer_sources list_primers resolve_primer fixed_sources=true disable_plain_primer_design=bound.with_value(|f|!f.identity.junction_edits.is_empty()) on_source=Some(diagnostic_source)/>
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
    next.junction_edits.clear();
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

fn default_regeneration(draft: &Draft, proposed: &[assembly_junction_edit::Edit]) -> Vec<bool> {
    let n = draft.slots.len();
    if n == 0 || n > 12 {
        return Vec::new();
    }
    let mut selected = (0..n)
        .map(|i| {
            draft.methods.rows.get(i).is_some_and(|row| {
                matches!(
                    row.preparation,
                    OverlapPreparation::ExactPcr | OverlapPreparation::PcrCatalogDigest
                ) && (row.forward.trim().is_empty() || row.reverse_primer.trim().is_empty())
            })
        })
        .collect::<Vec<_>>();
    let junctions = draft
        .junction_edits
        .iter()
        .chain(proposed)
        .map(|edit| edit.junction)
        .collect::<std::collections::BTreeSet<_>>();
    for junction in junctions {
        if junction >= n {
            continue;
        }
        let old = draft
            .junction_edits
            .iter()
            .find(|edit| edit.junction == junction);
        let new = proposed.iter().find(|edit| edit.junction == junction);
        if old != new {
            selected[junction] = true;
            if old.is_some_and(|edit| edit.remove_right > 0)
                || new.is_some_and(|edit| edit.remove_right > 0)
            {
                selected[(junction + 1) % n] = true;
            }
        }
    }
    selected
}
#[derive(Clone, PartialEq, Eq)]
struct JunctionSettings {
    overlap: String,
    minimum: String,
    maximum: String,
    tm: String,
}
impl Default for JunctionSettings {
    fn default() -> Self {
        Self {
            overlap: "20".into(),
            minimum: "18".into(),
            maximum: "30".into(),
            tm: "60".into(),
        }
    }
}
enum JunctionReport {
    Overlap(selective_junction_edit::Report),
    GoldenGate(golden_gate_junction_edit::Report),
}
impl JunctionReport {
    fn intent(&self) -> &assembly_junction_edit::Intent {
        match self {
            Self::Overlap(r) => &r.intent,
            Self::GoldenGate(r) => &r.intent,
        }
    }
    fn edits(&self) -> &[assembly_junction_edit::Edit] {
        match self {
            Self::Overlap(r) => &r.request.edits,
            Self::GoldenGate(r) => &r.request.edits,
        }
    }
    fn len(&self) -> usize {
        match self {
            Self::Overlap(r) => r.recipe.fragments.len(),
            Self::GoldenGate(r) => r.recipe.fragments.len(),
        }
    }
    fn fragment(&self, index: usize) -> Option<&collagene_core::assembly_recipe::FragmentSpec> {
        match self {
            Self::Overlap(r) => r.recipe.fragments.get(index),
            Self::GoldenGate(r) => r.recipe.fragments.get(index).map(|f| &f.source),
        }
    }
    fn retained(&self, index: usize) -> Option<bool> {
        match self {
            Self::Overlap(r) => r
                .request
                .rows
                .get(index)
                .map(|row| matches!(row, selective_junction_edit::Row::Retain(_))),
            Self::GoldenGate(r) => r
                .request
                .rows
                .get(index)
                .map(|row| matches!(row, golden_gate_junction_edit::Row::Retain(_))),
        }
    }
    fn designed_count(&self) -> usize {
        match self {
            Self::Overlap(r) => r.designs.iter().filter(|d| d.is_some()).count(),
            Self::GoldenGate(r) => r.designs.iter().filter(|d| d.is_some()).count(),
        }
    }
    fn retained_rows(&self) -> &[usize] {
        match self {
            Self::Overlap(r) => &r.retained_rows,
            Self::GoldenGate(r) => &r.retained_rows,
        }
    }
    fn warnings(&self) -> &[String] {
        match self {
            Self::Overlap(r) => &r.warnings,
            Self::GoldenGate(r) => &r.warnings,
        }
    }
    fn validate_binding(&self, revision: &str) -> Result<(), collagene_core::Error> {
        match self {
            Self::Overlap(r) => r.validate_binding(revision),
            Self::GoldenGate(r) => r.validate_binding(revision),
        }
    }
}
struct JunctionProposal {
    before: Draft,
    report: JunctionReport,
    sources: Vec<Arc<CapturedSource>>,
    settings: JunctionSettings,
    regenerate: Vec<bool>,
    epoch: u64,
}
#[component]
fn JunctionEditor(
    draft: RwSignal<Draft>,
    epoch: RwSignal<u64>,
    on_apply: Callback<Arc<JunctionProposal>, Result<(), String>>,
) -> impl IntoView {
    let proposed = RwSignal::new(draft.with_untracked(|d| d.junction_edits.clone()));
    let junction = RwSignal::new(0usize);
    let remove_left = RwSignal::new("0".to_string());
    let remove_right = RwSignal::new("0".to_string());
    let inserted = RwSignal::new(String::new());
    let settings = RwSignal::new(JunctionSettings {
        overlap: draft.with_untracked(|d| d.methods.overlaps.clone()),
        ..Default::default()
    });
    let regenerate = RwSignal::new(default_regeneration(
        &draft.get_untracked(),
        &proposed.get_untracked(),
    ));
    let manual_regeneration = RwSignal::new(false);
    let result = RwSignal::new(None::<Arc<JunctionProposal>>);
    let status = RwSignal::new(String::new());
    let pending = RwSignal::new(false);
    let generation = RwSignal::new(0u32);
    let worker = StoredValue::new_local(None::<Rc<CachedWorker>>);
    let alive = Arc::new(AtomicBool::new(true));
    let cleanup = alive.clone();
    on_cleanup(move || {
        cleanup.store(false, Ordering::Relaxed);
        if let Some(w) = worker.get_value() {
            w.close();
        }
    });
    let source_context = Memo::new(move |_| draft.with(Draft::identity));
    Effect::new(move |_| {
        let _ = (source_context.get(), epoch.get());
        let edits = draft.with_untracked(|d| d.junction_edits.clone());
        let first = edits.iter().find(|edit| edit.junction == 0);
        remove_left.set(
            first
                .map(|e| e.remove_left.to_string())
                .unwrap_or_else(|| "0".into()),
        );
        remove_right.set(
            first
                .map(|e| e.remove_right.to_string())
                .unwrap_or_else(|| "0".into()),
        );
        inserted.set(first.map(|e| e.inserted.clone()).unwrap_or_default());
        regenerate.set(default_regeneration(&draft.get_untracked(), &edits));
        manual_regeneration.set(false);
        settings.update(|settings| {
            settings.overlap = draft.with_untracked(|d| d.methods.overlaps.clone())
        });
        proposed.set(edits);
        junction.set(0);
    });
    Effect::new(move |_| {
        let before = draft.get();
        let edits = proposed.get();
        if !manual_regeneration.get() {
            let defaults = default_regeneration(&before, &edits);
            if regenerate.get_untracked() != defaults {
                regenerate.set(defaults);
            }
        }
    });
    Effect::new(move |_| {
        let _ = (
            draft.get(),
            proposed.get(),
            settings.get(),
            regenerate.get(),
            epoch.get(),
            remove_left.get(),
            remove_right.get(),
            inserted.get(),
            junction.get(),
        );
        generation.update(|g| *g = g.wrapping_add(1));
        if let Some(w) = worker.get_value() {
            w.close();
        }
        worker.set_value(None);
        pending.set(false);
        result.set(None);
    });
    let load_edit = Callback::new(move |index: usize| {
        junction.set(index);
        let edit =
            proposed.with_untracked(|edits| edits.iter().find(|e| e.junction == index).cloned());
        remove_left.set(
            edit.as_ref()
                .map(|e| e.remove_left.to_string())
                .unwrap_or_else(|| "0".into()),
        );
        remove_right.set(
            edit.as_ref()
                .map(|e| e.remove_right.to_string())
                .unwrap_or_else(|| "0".into()),
        );
        inserted.set(edit.map(|e| e.inserted).unwrap_or_default());
    });
    let unstaged = Memo::new(move |_| {
        let current = junction.get();
        proposed.with(|edits| {
            let edit = edits.iter().find(|edit| edit.junction == current);
            remove_left.get().parse::<usize>().ok()
                != Some(edit.map(|e| e.remove_left).unwrap_or(0))
                || remove_right.get().parse::<usize>().ok()
                    != Some(edit.map(|e| e.remove_right).unwrap_or(0))
                || inserted.get() != edit.map(|e| e.inserted.as_str()).unwrap_or("")
        })
    });
    let set_edit = move |_| {
        let outcome = (|| -> Result<_, String> {
            let junction = junction.get_untracked();
            let count =
                draft.with_untracked(|d| d.slots.len().saturating_sub(1) + usize::from(d.circular));
            if junction >= count {
                return Err("Choose an available junction.".into());
            }
            let left = remove_left.get_untracked();
            let right = remove_right.get_untracked();
            let inserted = inserted.get_untracked();
            if left.len() > 20 || right.len() > 20 || inserted.len() > 200 {
                return Err(
                    "Trim fields or inserted text exceed their bounds (200 raw inserted bytes)."
                        .into(),
                );
            }
            Ok(assembly_junction_edit::Edit {
                junction,
                remove_left: left.parse().map_err(|_| "Enter an upstream trim length.")?,
                remove_right: right
                    .parse()
                    .map_err(|_| "Enter a downstream trim length.")?,
                inserted,
            })
        })();
        match outcome {
            Ok(edit) => {
                proposed.update(|items| {
                    if let Some(old) = items.iter_mut().find(|e| e.junction == edit.junction) {
                        *old = edit;
                    } else {
                        items.push(edit);
                        items.sort_by_key(|e| e.junction);
                    }
                });
                status.set("Proposed edit updated. Preview all junction edits and PCR primers before applying.".into());
            }
            Err(error) => status.set(error),
        }
    };
    let design_alive = alive.clone();
    let preview = Callback::new(move |()| {
        let before = draft.get_untracked();
        let edits = proposed.get_untracked();
        let form = settings.get_untracked();
        let selected_rows = regenerate.get_untracked();
        let captured_epoch = epoch.get_untracked();
        let bound_edits = edits.clone();
        let prepared = (|| -> Result<_, String> {
            if !matches!(before.methods.mode, Mode::Overlap | Mode::GoldenGate) {
                return Err(
                    "Choose overlap or Golden Gate mode before preparing junction edits.".into(),
                );
            }
            if edits.is_empty() && before.junction_edits.is_empty() {
                return Err("Set at least one proposed junction edit first.".into());
            }
            let mut candidate = before.clone();
            candidate.junction_edits = edits.clone();
            let (inputs, base, sources) = inputs_request(&candidate)?;
            if selected_rows.len() != before.slots.len() {
                return Err("Select a preparation choice for every slot.".into());
            }
            let mut retained = before.methods.clone();
            for (i, selected) in selected_rows.iter().enumerate() {
                if *selected {
                    retained.rows[i] = OverlapInput {
                        document: i,
                        reverse: before.slots[i].reverse,
                        ..Default::default()
                    };
                }
            }
            let (_, fragments) = crate::cloning_methods::prepared_rows(&retained, inputs.len())?;
            let (revision, job) = if before.methods.mode == Mode::GoldenGate {
                let gate = &before.methods.golden_gate;
                let (forward_length, reverse_length) =
                    if selected_rows.iter().any(|selected| *selected) {
                        (
                            gate.forward_length
                                .parse()
                                .map_err(|_| "Enter a Golden Gate forward annealing length.")?,
                            gate.reverse_length
                                .parse()
                                .map_err(|_| "Enter a Golden Gate reverse annealing length.")?,
                        )
                    } else {
                        (0, 0)
                    }; // Not serialized: Retain rows carry no PCR length fields.

                let rows = fragments
                    .into_iter()
                    .enumerate()
                    .map(|(i, fragment)| {
                        if selected_rows[i] {
                            golden_gate_junction_edit::Row::DesignPcr {
                                role: fragment.role,
                                input_id: fragment.input_id,
                                reverse: fragment.reverse,
                                forward_length,
                                reverse_length,
                            }
                        } else {
                            golden_gate_junction_edit::Row::Retain(
                                collagene_core::ligation_recipe::Fragment {
                                    source: fragment,
                                    treatment: before.methods.rows[i].treatment.operation(),
                                },
                            )
                        }
                    })
                    .collect();
                let terminal_right_overhang =
                    if !before.circular && selected_rows.last() == Some(&true) {
                        Some(
                            before
                                .methods
                                .rows
                                .last()
                                .ok_or("Missing terminal fragment.")?
                                .right_fusion
                                .clone(),
                        )
                    } else {
                        None
                    };
                let request = golden_gate_junction_edit::Request {
                    base,
                    edits,
                    rows,
                    enzyme: gate.enzyme,
                    clamp: gate.clamp.clone(),
                    spacer: gate.spacer.clone(),
                    terminal_right_overhang,
                    allow_nicked: before.methods.allow_nicked,
                };
                let revision = golden_gate_junction_edit::design_revision(&inputs, &request)
                    .map_err(|e| e.to_string())?;
                (
                    revision,
                    AnalysisJob::DesignGoldenGateJunctionPrimers { inputs, request },
                )
            } else {
                if form.overlap.len() > 256
                    || [&form.minimum, &form.maximum, &form.tm]
                        .iter()
                        .any(|text| text.len() > 32)
                {
                    return Err("Primer settings exceed their field bounds.".into());
                }
                let rows = fragments
                    .into_iter()
                    .enumerate()
                    .map(|(i, fragment)| {
                        if selected_rows[i] {
                            selective_junction_edit::Row::DesignPcr {
                                role: fragment.role,
                                input_id: fragment.input_id,
                                reverse: fragment.reverse,
                            }
                        } else {
                            selective_junction_edit::Row::Retain(fragment)
                        }
                    })
                    .collect();
                let request = selective_junction_edit::Request {
                    base,
                    edits,
                    rows,
                    overlaps: form
                        .overlap
                        .split(',')
                        .map(|part| {
                            part.trim().parse::<usize>().map_err(
                                |_| "Enter one comma-separated homology length per junction.",
                            )
                        })
                        .collect::<Result<Vec<_>, _>>()?,
                    primer_options: collagene_core::primers::DesignOptions {
                        min_length: form
                            .minimum
                            .parse()
                            .map_err(|_| "Enter a minimum annealing length.")?,
                        max_length: form
                            .maximum
                            .parse()
                            .map_err(|_| "Enter a maximum annealing length.")?,
                        target_tm_c: form
                            .tm
                            .parse()
                            .map_err(|_| "Enter an annealing target Tm.")?,
                        max_pairs: 1,
                        ..Default::default()
                    },
                    profile_options: before.methods.profile_options()?,
                };
                let revision = selective_junction_edit::design_revision(&inputs, &request)
                    .map_err(|e| e.to_string())?;
                (
                    revision,
                    AnalysisJob::DesignSelectiveAssemblyPrimers { inputs, request },
                )
            };
            let w = Rc::new(CachedWorker::new()?);
            Ok((job, sources, revision, w))
        })();
        let (job, sources, revision, w) = match prepared {
            Ok(value) => value,
            Err(error) => {
                status.set(error);
                return;
            }
        };
        if let Some(old) = worker.get_value() {
            old.close();
        }
        worker.set_value(Some(w.clone()));
        let token = generation.get_untracked().wrapping_add(1);
        generation.set(token);
        pending.set(true);
        result.set(None);
        status.set(
            "Building the edited intended reference and proposing original-source PCR primers…"
                .into(),
        );
        let alive = design_alive.clone();
        spawn_local(async move {
            let response = w
                .ask(Request {
                    id: token,
                    revision: revision.clone(),
                    job,
                })
                .await;
            w.close();
            if !alive.load(Ordering::Relaxed)
                || generation.get_untracked() != token
                || epoch.get_untracked() != captured_epoch
                || draft.get_untracked() != before
                || settings.get_untracked() != form
                || proposed.get_untracked() != bound_edits
                || regenerate.get_untracked() != selected_rows
            {
                return;
            }
            pending.set(false);
            worker.set_value(None);
            let outcome = (|| -> Result<_, String> {
                let response = response?;
                if response.id != token || response.revision != revision {
                    return Err("Junction proposal response identity changed.".into());
                }
                let report = match response.result {
                    AnalysisResult::SelectiveAssemblyPrimers(report)
                        if before.methods.mode == Mode::Overlap =>
                    {
                        JunctionReport::Overlap(report)
                    }
                    AnalysisResult::GoldenGateJunctionPrimers(report)
                        if before.methods.mode == Mode::GoldenGate =>
                    {
                        JunctionReport::GoldenGate(report)
                    }
                    AnalysisResult::Error(error) => return Err(error),
                    _ => return Err("Unexpected junction primer proposal result.".into()),
                };
                report
                    .validate_binding(&revision)
                    .map_err(|e| e.to_string())?;
                if report.edits() != bound_edits.as_slice() {
                    return Err("Junction proposal edits differ from the requested changes.".into());
                }
                Ok(Arc::new(JunctionProposal {
                    before,
                    report,
                    sources,
                    settings: form,
                    regenerate: selected_rows,
                    epoch: captured_epoch,
                }))
            })();
            match outcome {
                Ok(proposal) => {
                    result.set(Some(proposal));
                    status.set("Proposal ready. Only checked rows will receive new original-source PCR primers; retained preparations and oligos stay unchanged.".into());
                }
                Err(error) => status.set(error),
            }
        });
    });
    let apply = move |_| {
        let outcome = (|| -> Result<(), String> {
            let proof = result
                .get_untracked()
                .ok_or("Preview the proposed junction edits first.")?;
            if proof.before != draft.get_untracked()
                || proof.epoch != epoch.get_untracked()
                || proof.settings != settings.get_untracked()
                || proof.report.edits() != proposed.get_untracked().as_slice()
                || proof.regenerate != regenerate.get_untracked()
            {
                return Err(
                    "Sources, proposed edits or primer settings changed. Preview again.".into(),
                );
            }
            on_apply.run(proof)
        })();
        if let Err(error) = outcome {
            status.set(error);
        }
    };
    view! {<details class="assembly-junction-editor"><summary>"Edit a junction and redesign its PCR preparation"</summary>
        <p class="helper">"Trim bases from the upstream core's assembly-facing end and/or the downstream core's beginning, then insert literal DNA between them. Original source regions stay pinned. Choose which rows receive newly designed PCR primers. Apply replaces only those checked rows; retained manual primers and digest preparations stay unchanged. This does not establish cloning success."</p>
        <label>"Junction"<select aria-label="Assembly junction to edit" prop:value=move||junction.get().to_string() on:change=move|e|if let Ok(index)=event_target_value(&e).parse(){load_edit.run(index)}>{move||draft.with(|d|{let count=d.slots.len().saturating_sub(1)+usize::from(d.circular);(0..count).map(|i|view!{<option value=i.to_string()>{format!("{} → {}{}",d.slots[i].label,d.slots[(i+1)%d.slots.len()].label,if i+1==d.slots.len(){" (circle closure)"}else{""})}</option>}).collect_view()})}</select></label>
        <div class="primer-grid"><label>"Remove upstream bases"<input type="number" min="0" aria-label="Junction upstream trim" prop:value=move||remove_left.get() on:input=move|e|remove_left.set(event_target_value(&e))/></label><label>"Remove downstream bases"<input type="number" min="0" aria-label="Junction downstream trim" prop:value=move||remove_right.get() on:input=move|e|remove_right.set(event_target_value(&e))/></label></div>
        <label>"Insert literal DNA at this junction"<textarea maxlength="200" aria-label="Junction inserted DNA" prop:value=move||inserted.get() on:input=move|e|inserted.set(event_target_value(&e))></textarea></label>
        <button class="button secondary" on:click=set_edit>"Set proposed junction edit"</button>
        {move||proposed.get().into_iter().map(|edit|{let index=edit.junction;view!{<p>{format!("Junction {}: remove {} upstream / {} downstream; insert {}",index+1,edit.remove_left,edit.remove_right,if edit.inserted.is_empty(){"(none)"}else{&edit.inserted})}<button class="text-tool" on:click=move|_|{proposed.update(|items|items.retain(|e|e.junction!=index));if junction.get_untracked()==index{remove_left.set("0".into());remove_right.set("0".into());inserted.set(String::new());}}>"Remove proposed edit"</button></p>}}).collect_view()}
        <h4>"Preparation choices"</h4><p class="helper">"Recommended boxes start with changed upstream fragments, downstream trims and unconfigured PCR rows. After you change a box, later edit changes do not override that choice. Retaining an incompatible preparation can prevent reconstruction; it is never silently replaced."</p>
        {move||draft.with(|d|d.slots.iter().enumerate().map(|(i,slot)|view!{<label class="checkbox-label"><input type="checkbox" aria-label=format!("Regenerate PCR primers for slot {}",i+1) prop:checked=move||regenerate.with(|rows|rows.get(i).copied().unwrap_or(false)) on:change=move|e|{manual_regeneration.set(true);regenerate.update(|rows|if let Some(row)=rows.get_mut(i){*row=event_target_checked(&e)});}/>{format!("{} — regenerate PCR primers",slot.label)}</label><p class="helper">{move||draft.with(|current|current.methods.rows.get(i).map(|row|format!("{}: {}",if regenerate.with(|rows|rows.get(i).copied().unwrap_or(false)){"Replace after explicit Apply"}else{"Retain exact current preparation"},row.summary())).unwrap_or_default())}</p>}).collect_view())}
        <Show when=move||draft.with(|d|d.methods.mode==Mode::Overlap)>
        <p class="helper">{move||draft.with(|d|format!("Enter {} explicit homology lengths in junction order, including circle closure when present. Existing values are preserved; no uniform length is inferred.",d.slots.len().saturating_sub(1)+usize::from(d.circular)))}</p>
        <div class="primer-grid"><label>"Per-junction homology lengths (comma-separated)"<input aria-label="Edited junction homology" prop:value=move||settings.with(|s|s.overlap.clone()) on:input=move|e|settings.update(|s|s.overlap=event_target_value(&e))/></label><label>"Minimum annealing length"<input type="number" min="1" max="200" aria-label="Edited junction primer minimum" prop:value=move||settings.with(|s|s.minimum.clone()) on:input=move|e|settings.update(|s|s.minimum=event_target_value(&e))/></label><label>"Maximum annealing length"<input type="number" min="1" max="200" aria-label="Edited junction primer maximum" prop:value=move||settings.with(|s|s.maximum.clone()) on:input=move|e|settings.update(|s|s.maximum=event_target_value(&e))/></label><label>"Target annealing Tm (°C)"<input type="number" aria-label="Edited junction primer target Tm" prop:value=move||settings.with(|s|s.tm.clone()) on:input=move|e|settings.update(|s|s.tm=event_target_value(&e))/></label></div>
        <p class="helper">"The selected Gibson/HiFi/In-Fusion profile is retained. Primer metrics are model calculations, not reagent conditions or demonstrated binding. Concrete source and inserted DNA are required for these primer proposals; ambiguous calls are not silently resolved."</p>
        </Show>
        <Show when=move||draft.with(|d|d.methods.mode==Mode::GoldenGate)>
            <h4>"Golden Gate junction preparation"</h4><p class="helper">"Checked rows receive PCR-added Type IIS sites and the literal junction insert. Internal outgoing fusion bases come from the next adjusted source core; no overlap Tm or homology-length setting is used. Retained rows keep their exact digest/primer/end-treatment choices."</p>
            <div class="primer-grid"><label>"Type IIS enzyme"<select aria-label="Junction Type IIS enzyme" prop:value=move||format!("{:?}",draft.with(|d|d.methods.golden_gate.enzyme)) on:change=move|e|draft.update(|d|d.methods.golden_gate.enzyme=match event_target_value(&e).as_str(){"BsmBI"=>collagene_core::golden_gate::GoldenGateEnzyme::BsmBI,"BbsI"=>collagene_core::golden_gate::GoldenGateEnzyme::BbsI,"SapI"=>collagene_core::golden_gate::GoldenGateEnzyme::SapI,_=>collagene_core::golden_gate::GoldenGateEnzyme::BsaI})><option value="BsaI">"BsaI · 4-base fusion"</option><option value="BsmBI">"BsmBI · 4-base fusion"</option><option value="BbsI">"BbsI · 4-base fusion"</option><option value="SapI">"SapI · 3-base fusion"</option></select></label>
            <label>"Pre-recognition clamp"<input maxlength="12" aria-label="Junction Golden Gate clamp" prop:value=move||draft.with(|d|d.methods.golden_gate.clamp.clone()) on:input=move|e|draft.update(|d|d.methods.golden_gate.clamp=event_target_value(&e))/></label><label>"Cleavage spacer"<input maxlength="2" aria-label="Junction Golden Gate spacer" prop:value=move||draft.with(|d|d.methods.golden_gate.spacer.clone()) on:input=move|e|draft.update(|d|d.methods.golden_gate.spacer=event_target_value(&e))/></label>
            <label>"Forward annealing length"<input type="number" min="1" max="200" aria-label="Junction Golden Gate forward length" prop:value=move||draft.with(|d|d.methods.golden_gate.forward_length.clone()) on:input=move|e|draft.update(|d|d.methods.golden_gate.forward_length=event_target_value(&e))/></label><label>"Reverse annealing length"<input type="number" min="1" max="200" aria-label="Junction Golden Gate reverse length" prop:value=move||draft.with(|d|d.methods.golden_gate.reverse_length.clone()) on:input=move|e|draft.update(|d|d.methods.golden_gate.reverse_length=event_target_value(&e))/></label></div>
            <Show when=move||draft.with(|d|!d.circular)&&regenerate.with(|selected|selected.last()==Some(&true))><label>"Terminal right fusion (assembly direction, 5′→3′)"<input maxlength="4" aria-label="Junction Golden Gate terminal fusion" prop:value=move||draft.with(|d|d.methods.rows.last().map(|row|row.right_fusion.clone()).unwrap_or_default()) on:input=move|e|draft.update(|d|if let Some(row)=d.methods.rows.last_mut(){row.right_fusion=event_target_value(&e)})/></label><p class="helper">"A selected terminal linear fragment requires this explicit fusion DNA (3 bases for SapI; 4 for the other listed enzymes). It is not filled from the intended target."</p></Show>
            <label class="checkbox-label"><input type="checkbox" aria-label="Junction Golden Gate allow retained nicks" prop:checked=move||draft.with(|d|d.methods.allow_nicked) on:change=move|e|draft.update(|d|d.methods.allow_nicked=event_target_checked(&e))/>"Allow explicitly retained nicks"</label><p class="helper">"Clamp/spacer and terminal fusion are explicit design choices, not supplier conditions. Preview checks the actual PCR/digest ends and the full ordered ligation; it does not establish reaction efficiency or mixture uniqueness."</p>
        </Show>
        <Show when=move||unstaged.get()><p class="helper">"Set the proposed junction edit before previewing, so the changed trim and inserted bases are included."</p></Show>
        <button class="button primary" disabled=move||pending.get()||unstaged.get()||draft.with(|d|!matches!(d.methods.mode,Mode::Overlap|Mode::GoldenGate)) on:click=move|_|preview.run(())>"Preview edited construct and PCR primers"</button>
        <button class="text-tool" disabled=move||!pending.get() on:click=move|_|{generation.update(|g|*g=g.wrapping_add(1));if let Some(w)=worker.get_value(){w.close();}worker.set_value(None);pending.set(false);result.set(None);status.set("Junction proposal cancelled.".into());}>"Cancel junction proposal"</button>
        <p role="status">{move||status.get()}</p>
        {move||result.get().map(|proof|view!{<section class="assembly-junction-proposal"><h3>"Edited construct and proposed PCR preparation"</h3><p>{format!("{} · {} bases",proof.report.intent().target.name,proof.report.intent().target.sequence.len())}</p><pre>{let dna=proof.report.intent().target.sequence.as_str();format!("{}{}",dna.chars().take(240).collect::<String>(),if dna.len()>240{"… (bounded preview)"}else{""})}</pre>
            {proof.report.intent().fragments.iter().enumerate().map(|(i,fragment)|view!{<p>{format!("{}: original core {}–{} → adjusted core {}–{}{}; output span {}–{} includes any outgoing synthetic insert",proof.before.slots[i].label,proof.before.slots[i].start+1,proof.before.slots[i].end,fragment.start+1,fragment.end,if fragment.reverse{" reverse"}else{" forward"},fragment.output_start+1,fragment.output_end)}</p>}).collect_view()}
            <p>{format!("Retained preparations: {} · regenerated PCR rows: {}",proof.report.retained_rows().len(),proof.report.designed_count())}</p>
            {proof.report.retained_rows().iter().filter_map(|&i|proof.before.methods.rows.get(i).map(|row|(i,row))).map(|(i,row)|view!{<p class="helper">{format!("Retained {}: {}",proof.before.slots[i].label,row.summary())}</p>}).collect_view()}
            {junction_primers(&proof.report)}
            {proof.report.warnings().iter().chain(proof.report.intent().warnings.iter()).map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}
            <button class="button primary" on:click=apply>"Apply junction edits and selected PCR primers"</button>
        </section>})}
    </details>}
}

fn junction_primers(report: &JunctionReport) -> AnyView {
    match report {
        JunctionReport::Overlap(report)=>report.designs.iter().enumerate().filter_map(|(i,design)|design.as_ref().map(|design|(i,design))).map(|(i,design)|view!{<article><h4>{format!("Fragment {} proposed overlap PCR primers",i+1)}</h4>{design.pairs.first().map(|pair|view!{<pre>{format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\nAdded tail: {}\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\nAdded tail: {}",pair.forward.primer.name,pair.forward.primer.tail,pair.forward.primer.annealing.as_str(),pair.forward.primer.annealing.as_str(),pair.forward.primer.tail,pair.reverse.primer.name,pair.reverse.primer.tail,pair.reverse.primer.annealing.as_str(),pair.reverse.primer.annealing.as_str(),pair.reverse.primer.tail)}</pre><p>{format!("Annealing Tm {} / {} · pair ΔTm {:.1}°C",pair.forward.analysis.tm.as_ref().map(|t|format!("{:.1}°C",t.tm_c)).unwrap_or_else(||"unavailable".into()),pair.reverse.analysis.tm.as_ref().map(|t|format!("{:.1}°C",t.tm_c)).unwrap_or_else(||"unavailable".into()),pair.tm_difference_c)}</p>})}{design.warnings.iter().map(|warning|view!{<p class="helper">{warning.clone()}</p>}).collect_view()}</article>}).collect_view().into_any(),
        JunctionReport::GoldenGate(report)=>report.designs.iter().enumerate().filter_map(|(i,design)|design.as_ref().map(|design|(i,design))).map(|(i,design)|view!{<article><h4>{format!("Fragment {} proposed Type IIS PCR/digest chain",i+1)}</h4><pre>{format!("{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\nTail: clamp [{}] · motif [{}] · spacer [{}] · fusion [{}] · inserted [{}]\n\n{}\nOrdered 5′→3′: {}{}\nAnnealing: {}\nTail: clamp [{}] · motif [{}] · spacer [{}] · fusion [{}] · inserted [{}]",design.forward.name,design.forward.tail,design.forward.annealing.as_str(),design.forward.annealing.as_str(),design.forward_tail.clamp,design.forward_tail.motif,design.forward_tail.spacer,design.forward_tail.fusion,design.forward_tail.inserted,design.reverse.name,design.reverse.tail,design.reverse.annealing.as_str(),design.reverse.annealing.as_str(),design.reverse_tail.clamp,design.reverse_tail.motif,design.reverse_tail.spacer,design.reverse_tail.fusion,design.reverse_tail.inserted)}</pre><p>{format!("Actual digested fragment {} · assembly-oriented top DNA: {}{}",design.fragment_index+1,design.oriented_top.chars().take(240).collect::<String>(),if design.oriented_top.len()>240{"… (bounded preview)"}else{""})}</p><p>{format!("Left end: {}",crate::cloning_methods::method_end(Some(&design.left)))}</p><p>{format!("Right end: {}",crate::cloning_methods::method_end(Some(&design.right)))}</p><p class="helper">"Tail segments are shown in each actual oligo's own 5′→3′ direction. Literal insert and fusion segments may therefore be reverse-complemented; the displayed digested top DNA is in assembly direction."</p></article>}).collect_view().into_any(),
    }
}
