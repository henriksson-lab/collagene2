//! Editable assembly data and versioned persistence, without GUI state or receipts.
use crate::{
    cloning_methods_model::{MethodsState, Mode, OverlapInput, OverlapPreparation},
    source_identity::SourceInfo,
};
use collagene_core::{
    Alphabet, Molecule, Topology, assembly_edit_intent, assembly_intent, assembly_junction_edit,
    assembly_substitution_intent, verification::PlanInput,
};
use serde::{Deserialize, Serialize};
use std::{io::Write, sync::Arc};

pub const CURRENT_VERSION: u32 = collagene_core::assembly_library_history::ACTIVE_DRAFT_VERSION;
pub const ATTACHMENT_KIND: &str = collagene_core::assembly_library_history::ACTIVE_DRAFT_KIND;
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
    #[serde(deserialize_with = "required_option")]
    pub material: Option<MaterialDomain>,
}
pub use collagene_core::assembly_intent::MaterialDomain;
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Slot {
    pub label: String,
    #[serde(with = "optional_snapshot")]
    pub captured: Option<Arc<CapturedSource>>,
    pub selector: assembly_intent::Selector,
    #[serde(deserialize_with = "required_option")]
    pub digest_geometry: Option<DigestGeometry>,
    #[serde(with = "optional_snapshot")]
    pub anneal_bottom: Option<Arc<CapturedSource>>,
}
impl Slot {
    pub(crate) fn empty(label: &str) -> Self {
        Self {
            label: label.into(),
            captured: None,
            selector: assembly_intent::Selector::ReferenceRange {
                start: 0,
                end: 0,
                reverse: false,
            },
            digest_geometry: None,
            anneal_bottom: None,
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
    #[serde(deserialize_with = "deserialize_junction_edits")]
    pub junction_edits: Vec<assembly_junction_edit::Edit>,
    #[serde(deserialize_with = "deserialize_substitutions")]
    pub substitutions: Vec<collagene_core::assembly_substitution_intent::Substitution>,
    #[serde(deserialize_with = "deserialize_internal_edits")]
    pub edits: Vec<collagene_core::assembly_edit_intent::Edit>,
    #[serde(deserialize_with = "required_option")]
    pub alternatives: Option<crate::assembly_alternatives::Library>,
}
fn required_option<'de, D: serde::Deserializer<'de>, T: Deserialize<'de>>(
    deserializer: D,
) -> Result<Option<T>, D::Error> {
    Option::deserialize(deserializer)
}
fn deserialize_internal_edits<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<Vec<collagene_core::assembly_edit_intent::Edit>, D::Error> {
    struct Visitor;
    impl<'de> serde::de::Visitor<'de> for Visitor {
        type Value = Vec<collagene_core::assembly_edit_intent::Edit>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("at most256 literal contribution edits")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            if seq.size_hint().is_some_and(|n| n > 256) {
                return Err(serde::de::Error::custom("at most256 contribution edits"));
            }
            let mut result = Vec::new();
            while let Some(value) = seq.next_element()? {
                if result.len() == 256 {
                    return Err(serde::de::Error::custom("at most256 contribution edits"));
                }
                result.push(value);
            }
            Ok(result)
        }
    }
    deserializer.deserialize_seq(Visitor)
}
fn deserialize_substitutions<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<Vec<collagene_core::assembly_substitution_intent::Substitution>, D::Error> {
    struct Visitor;
    impl<'de> serde::de::Visitor<'de> for Visitor {
        type Value = Vec<collagene_core::assembly_substitution_intent::Substitution>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("at most256 literal substitutions")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            if seq.size_hint().is_some_and(|n| n > 256) {
                return Err(serde::de::Error::custom("at most256 substitutions"));
            }
            let mut result = Vec::new();
            while let Some(value) = seq.next_element()? {
                if result.len() == 256 {
                    return Err(serde::de::Error::custom("at most256 substitutions"));
                }
                result.push(value);
            }
            Ok(result)
        }
    }
    deserializer.deserialize_seq(Visitor)
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
            version: CURRENT_VERSION,
            name: "New assembly".into(),
            circular: true,
            slots: vec![Slot::empty("Backbone"), Slot::empty("Insert")],
            active_slot: 0,
            methods: MethodsState::default(),
            junction_edits: Vec::new(),
            substitutions: Vec::new(),
            edits: Vec::new(),
            alternatives: None,
        };
        value.sync_rows();
        for row in &mut value.methods.rows {
            row.preparation = OverlapPreparation::ExactPcr;
        }
        value
    }
}
pub(crate) fn token(text: &str, max: usize) -> bool {
    !text.trim().is_empty() && text.len() <= max && !text.chars().any(char::is_control)
}
pub(crate) fn hex(text: &str) -> bool {
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
pub(crate) fn bounded_size(value: &impl Serialize) -> Result<(), String> {
    serde_json::to_writer(&mut Count(0), value).map_err(|e| e.to_string())
}
pub(crate) fn validate_preparation_input(row: &OverlapInput) -> Result<(), String> {
    if let Some(raw) = &row.source_generation {
        raw.validate_raw()?;
    }

    if let Some(raw) = &row.overlap_extension_pcr {
        raw.validate_raw()?;
    }
    if [
        &row.mismatch_max_mismatches,
        &row.mismatch_three_prime_exact,
        &row.mismatch_max_hits,
        &row.mismatch_max_work,
        &row.mismatch_max_pairs,
        &row.mismatch_max_products,
        &row.mismatch_max_product_bases,
    ]
    .iter()
    .any(|value| value.len() > 20)
    {
        return Err("Mismatch-PCR policy fields exceed20 bytes.".into());
    }
    let options = &row.anneal_options;
    if options.min_paired_bases == 0
        || options.min_paired_bases > collagene_core::chemistry::MAX_ANNEAL_OLIGO_BASES
        || options.max_comparisons == 0
        || options.max_comparisons > collagene_core::chemistry::MAX_ANNEAL_COMPARISONS
        || options.max_candidates == 0
        || options.max_candidates > collagene_core::chemistry::MAX_ANNEAL_CANDIDATES
        || row
            .anneal_candidate
            .is_some_and(|index| index >= options.max_candidates)
        || row
            .anneal_bottom_document
            .is_some_and(|index| index >= 1000)
        || row
            .anneal_source_context
            .as_ref()
            .is_some_and(|context| !hex(context))
    {
        return Err("Paired-oligo settings exceed their bounded native profile.".into());
    }
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
    Ok(())
}
pub(crate) fn validate_controls(
    name: &str,
    circular: bool,
    active_slot: usize,
    count: usize,
    methods: &MethodsState,
    junction_edits: &[assembly_junction_edit::Edit],
) -> Result<(), String> {
    if !token(&name, 256)
        || count == 0
        || count > 12
        || active_slot >= count
        || methods.rows.len() != count
        || methods.circular != circular
    {
        return Err("Invalid bounded assembly draft or unsupported draft mode".into());
    }
    if methods.mode == Mode::Tailing {
        if count != 2 {
            return Err("Ordinary TA/GC cloning requires exactly two source slots.".into());
        }
    }
    if methods.mode != Mode::Ligation
        && methods
            .rows
            .iter()
            .any(|row| row.preparation == OverlapPreparation::AnnealedOligos)
    {
        return Err(
            "Paired-oligo preparation currently requires ordinary ordered ligation.".into(),
        );
    }
    if let Some(pool) = &methods.golden_gate_pool {
        pool.validate()?;
        if !pool.choices.is_empty() && pool.choices.len() != count {
            return Err("Reaction substrate choices must match the assembly slots.".into());
        }
    }
    if let Some(form) = &methods.multisite {
        form.validate_shape()?;
    }
    let junction_count = count.saturating_sub(1) + usize::from(circular);
    if junction_edits.len() > 12
        || (!junction_edits.is_empty()
            && !matches!(
                methods.mode,
                Mode::Overlap
                    | Mode::GoldenGate
                    | Mode::MultiSiteReference
                    | Mode::Tailing
                    | Mode::Ligation
            ))
    {
        return Err(
            "Junction edits require declared reference intent, overlap, Golden Gate, TA/GC or ordered ligation and at most12 edits".into(),
        );
    }
    let mut junctions = std::collections::BTreeSet::new();
    for edit in junction_edits {
        let terminal = !circular
            && matches!(methods.mode, Mode::Tailing | Mode::Ligation)
            && edit.junction == count - 1;
        if (!terminal && edit.junction >= junction_count)
            || (terminal && edit.remove_right != 0)
            || !junctions.insert(edit.junction)
            || edit.inserted.len() > 200
        {
            return Err(
                "Invalid or duplicate junction edit; inserted text is limited to 200 bytes".into(),
            );
        }
    }
    for row in &methods.rows {
        validate_preparation_input(row)?;
    }
    let gate = &methods.golden_gate;
    if gate.clamp.len() > 12
        || gate.spacer.len() > 2
        || gate.forward_length.len() > 4
        || gate.reverse_length.len() > 4
    {
        return Err("Type IIS primer settings exceed their field bounds".into());
    }

    for text in [
        &methods.overlaps,
        &methods.comparisons,
        &methods.occurrences,
    ] {
        if text.len() > 65536 {
            return Err("Assembly option field exceeds64KiB".into());
        }
    }
    if let Some(form) = &methods.homology {
        form.validate_raw()?;
    }
    if let Some(form) = &methods.annealing {
        form.validate_raw()?;
    }
    methods.all_input_ids()?;
    Ok(())
}
pub(crate) fn validate_slot_link(
    slot: &Slot,
    row: &OverlapInput,
    index: usize,
) -> Result<(), String> {
    let (start, end, reverse) = match &slot.selector {
        assembly_intent::Selector::ReferenceRange {
            start,
            end,
            reverse,
        } => (
            start
                .checked_add(1)
                .ok_or("Reference interval coordinate overflow")?
                .to_string(),
            if slot.captured.is_some() {
                end.to_string()
            } else {
                String::new()
            },
            Some(*reverse),
        ),
        assembly_intent::Selector::MaterialStrandRange { .. } => {
            (String::new(), String::new(), None)
        }
    };
    if row.document != index
        || reverse.is_some_and(|reverse| row.reverse != reverse)
        || row.interval_start != start
        || row.interval_end != end
    {
        return Err("Draft preparation source/order differs from its pinned slots".into());
    }
    Ok(())
}
pub(crate) fn validate_slot(slot: &Slot, row: &OverlapInput, index: usize) -> Result<(), String> {
    if !token(&slot.label, 80) {
        return Err("Slot labels must be nonblank and at most80 bytes".into());
    }
    validate_preparation_input(row)?;
    validate_slot_link(slot, row, index)?;
    if let Some(origin) = &row.run_origin {
        origin.validate_documentary().map_err(|e| e.to_string())?;
        let capture = slot
            .captured
            .as_ref()
            .ok_or("Saved-run material needs its exact captured source.")?;
        if *origin.event() != capture.source.binding
            || capture.molecule.revision() != origin.event().revision
        {
            return Err("Saved-run material differs from the exact captured source event.".into());
        }
        if matches!(
            row.preparation,
            OverlapPreparation::DeclaredDuplex | OverlapPreparation::AnnealedOligos
        ) {
            return Err("Saved-run material cannot be replaced by declared duplex or paired-oligo preparation.".into());
        }
    }
    if let Some(bottom) = slot.anneal_bottom.as_deref() {
        validate_oligo_capture(bottom)?;
    }
    let Some(captured) = slot.captured.as_deref() else {
        if !matches!(
            &slot.selector,
            assembly_intent::Selector::ReferenceRange {
                start: 0,
                end: 0,
                ..
            }
        ) || slot.digest_geometry.is_some()
        {
            return Err("An empty slot cannot retain a source interval".into());
        }
        return Ok(());
    };
    captured.molecule.validate().map_err(|e| e.to_string())?;
    let source = &captured.source;
    let molecule = &captured.molecule;
    let n = molecule.sequence.len();
    if !token(&source.binding.document_id, 512)
        || !token(&source.binding.event_id, 512)
        || !hex(&source.binding.revision)
        || source.alphabet != Alphabet::Dna
        || source.alphabet != molecule.sequence.alphabet()
        || source.topology != molecule.topology
        || source.length != n
        || source.name != molecule.name
        || n == 0
    {
        return Err("Invalid pinned DNA source or one-turn interval".into());
    }
    match &slot.selector {
        assembly_intent::Selector::ReferenceRange { start, end, .. } => {
            if *start >= n
                || end
                    .checked_sub(*start)
                    .is_none_or(|length| length == 0 || length > n)
                || (molecule.topology == Topology::Linear && *end > n)
            {
                return Err("Invalid pinned reference interval".into());
            }
        }
        assembly_intent::Selector::MaterialStrandRange {
            expected_state_revision,
            ..
        } => {
            if row.run_origin.is_none() || !hex(expected_state_revision) {
                return Err("Actual strand intent needs its exact saved material origin".into());
            }
            let geometry = slot
                .digest_geometry
                .as_ref()
                .ok_or("Actual strand intent needs its captured source domain")?;
            if digest_selector(source, geometry, row.reverse)? != slot.selector {
                return Err(
                    "Actual strand selection differs from its pinned digest/source state".into(),
                );
            }
        }
    }
    if let Some(geometry) = &slot.digest_geometry {
        if digest_selector(source, geometry, row.reverse)? != slot.selector {
            return Err("Digest contribution differs from its selected actual strand".into());
        }
    }
    if molecule.revision() != source.binding.revision {
        return Err("Pinned source content differs from its captured revision".into());
    }
    if row.preparation == OverlapPreparation::AnnealedOligos {
        validate_oligo_capture(captured)?;
    }
    Ok(())
}
fn validate_oligo_capture(captured: &CapturedSource) -> Result<(), String> {
    let source = &captured.source;
    let molecule = &captured.molecule;
    let n = molecule.sequence.len();
    if n == 0
        || n > collagene_core::chemistry::MAX_ANNEAL_OLIGO_BASES
        || source.topology != Topology::Linear
        || molecule.topology != Topology::Linear
        || source.alphabet != Alphabet::Dna
        || molecule.sequence.alphabet() != Alphabet::Dna
        || source.length != n
        || source.name != molecule.name
        || !token(&source.binding.document_id, 512)
        || !token(&source.binding.event_id, 512)
        || !hex(&source.binding.revision)
    {
        return Err("Paired oligos require bounded nonempty linear original DNA sources and exact event bindings.".into());
    }
    molecule.validate().map_err(|error| error.to_string())?;
    if molecule.revision() != source.binding.revision {
        return Err("Captured oligo content differs from its original source revision.".into());
    }
    Ok(())
}
impl Draft {
    pub(crate) fn clear_slot_substitutions(&mut self, index: usize) -> Result<(), String> {
        let roles = self.methods.role_ids(self.slots.len())?;
        let role = roles.get(index).ok_or("The source slot is unavailable.")?;
        self.substitutions.retain(|edit| edit.role != role.role);
        self.edits.retain(|edit| edit.role != role.role);
        Ok(())
    }
    fn validate_internal_edits(&self) -> Result<(), String> {
        if self.edits.len() > 256 {
            return Err("At most256 explicit contribution edits are supported.".into());
        }
        if self.edits.is_empty() {
            return Ok(());
        }
        if self.alternatives.is_some() || self.methods.source_roles.is_none() {
            return Err("Internal contribution edits require an ordinary construct and stable source roles.".into());
        }
        let roles = self.methods.role_ids(self.slots.len())?;
        let mut total = 0usize;
        for (index, edit) in self.edits.iter().enumerate() {
            let length = edit
                .end
                .checked_sub(edit.start)
                .ok_or("Contribution edit end precedes its start.")?;
            total = total
                .checked_add(length)
                .and_then(|n| n.checked_add(edit.replacement.len()))
                .filter(|n| *n <= 100_000)
                .ok_or("Contribution edits exceed100000 affected/inserted bases.")?;
            if !token(&edit.role, 128)
                || !roles.iter().any(|r| r.role == edit.role)
                || edit.end > collagene_core::model::MAX_SEQUENCE_LEN
                || edit.expected.len() != length
                || edit.expected == edit.replacement
                || !edit
                    .expected
                    .bytes()
                    .chain(edit.replacement.bytes())
                    .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
            {
                return Err("Contribution edits require original-role coordinates and distinct uppercase DNA literals.".into());
            }
            // Insertions sharing a boundary with another edit are intentionally
            // rejected: their relative order is not declared by these literals.
            let conflict = |a: usize, b: usize, c: usize, d: usize| {
                if a == b {
                    c <= a && a <= d
                } else if c == d {
                    a <= c && c <= b
                } else {
                    a < d && c < b
                }
            };
            if self.edits[..index].iter().any(|other| {
                other.role == edit.role && conflict(edit.start, edit.end, other.start, other.end)
            }) || self.substitutions.iter().any(|other| {
                other.role == edit.role
                    && conflict(edit.start, edit.end, other.offset, other.offset + 1)
            }) {
                return Err("Contribution edits must not overlap each other or substitutions, or share insertion boundaries.".into());
            }
        }
        Ok(())
    }
    /// Explicitly reverse the declared pair by swapping its ordered originals.
    /// This changes independent intent; it does not choose an annealing result.
    /// Only the swapped row's incoming/outgoing boundary declarations are reset;
    /// unrelated intent remains literal. Every mutation is transactional.
    pub(crate) fn swap_oligo_roles(&self, index: usize) -> Result<Draft, String> {
        self.validate()?;
        if self.methods.mode != Mode::Ligation || self.alternatives.is_some() {
            return Err("Swap ordered oligos requires an ordinary ligation draft.".into());
        }
        let row = self
            .methods
            .rows
            .get(index)
            .ok_or("The paired preparation row is unavailable.")?;
        let slot = self
            .slots
            .get(index)
            .ok_or("The paired source slot is unavailable.")?;
        if row.preparation != OverlapPreparation::AnnealedOligos
            || slot.captured.is_none()
            || slot.anneal_bottom.is_none()
        {
            return Err(
                "Capture both original strands of an active paired-oligo preparation first.".into(),
            );
        }
        let mut next = self.clone();
        next.clear_slot_substitutions(index)?;
        let slot = &mut next.slots[index];
        std::mem::swap(&mut slot.captured, &mut slot.anneal_bottom);
        slot.selector = assembly_intent::Selector::ReferenceRange {
            start: 0,
            end: slot
                .captured
                .as_ref()
                .ok_or("The new top oligo is unavailable.")?
                .molecule
                .sequence
                .len(),
            reverse: false,
        };
        slot.digest_geometry = None;
        let row = &mut next.methods.rows[index];
        std::mem::swap(
            &mut row.anneal_options.top_chemistry,
            &mut row.anneal_options.bottom_chemistry,
        );
        row.anneal_candidate = None;
        row.anneal_source_context = None;
        row.anneal_bottom_document = None;
        let incoming = if index > 0 {
            Some(index - 1)
        } else if next.circular {
            Some(next.slots.len() - 1)
        } else {
            None
        };
        next.junction_edits
            .retain(|edit| edit.junction != index && Some(edit.junction) != incoming);
        next.sync_rows();
        next.validate()?;
        Ok(next)
    }

    fn uses_paired_fields(&self) -> bool {
        self.slots.iter().any(|slot| slot.anneal_bottom.is_some())
            || self.methods.rows.iter().any(|row| {
                row.preparation == OverlapPreparation::AnnealedOligos
                    || row.anneal_options != Default::default()
                    || row.anneal_candidate.is_some()
                    || row.anneal_bottom_document.is_some()
                    || row.anneal_source_context.is_some()
            })
            || self.alternatives.as_ref().is_some_and(|library| {
                library.bins.iter().any(|bin| {
                    bin.choices
                        .iter()
                        .any(|choice| choice.slot.anneal_bottom.is_some())
                }) || library.rows.iter().any(|saved| {
                    saved.edits.methods.rows.iter().any(|row| {
                        row.preparation == OverlapPreparation::AnnealedOligos
                            || row.anneal_options != Default::default()
                            || row.anneal_candidate.is_some()
                            || row.anneal_bottom_document.is_some()
                            || row.anneal_source_context.is_some()
                    })
                })
            })
    }
    /// A documentary pin for an explicit candidate index, never physical evidence.
    pub(crate) fn anneal_context(&self, index: usize) -> Option<String> {
        let row = self.methods.rows.get(index)?;
        if row.preparation != OverlapPreparation::AnnealedOligos {
            return None;
        }
        let slot = self.slots.get(index)?;
        let top = &slot.captured.as_ref()?.source.binding;
        let bottom = &slot.anneal_bottom.as_ref()?.source.binding;
        let input_id = self
            .methods
            .role_ids(self.slots.len())
            .ok()?
            .get(index)?
            .input_id
            .clone();
        let mut hasher = blake3::Hasher::new();
        for value in [
            collagene_core::chemistry::ANNEAL_ENGINE_VERSION,
            input_id.as_str(),
            &top.document_id,
            &top.event_id,
            &top.revision,
            &bottom.document_id,
            &bottom.event_id,
            &bottom.revision,
        ] {
            if value.len() > 512 {
                return None;
            }
            hasher.update(&(value.len() as u64).to_le_bytes());
            hasher.update(value.as_bytes());
        }
        let options = serde_json::to_vec(&row.anneal_options).ok()?;
        hasher.update(&options);
        Some(hasher.finalize().to_hex().to_string())
    }
    /// Structural/bounded source consistency only. Historical event IDs require
    /// the host workspace to attest them before installing any candidate.
    pub fn validate(&self) -> Result<(), String> {
        if self.version != CURRENT_VERSION {
            return Err("Unsupported assembly draft version.".into());
        }
        if self.substitutions.len() > 256 {
            return Err("At most256 explicit substitutions are supported.".into());
        }
        if !self.substitutions.is_empty() {
            if self.alternatives.is_some() {
                return Err(
                    "Internal substitutions are currently supported only in ordinary assemblies."
                        .into(),
                );
            }
            if self.methods.source_roles.is_none() {
                return Err(
                    "Pin stable source roles before declaring internal substitutions.".into(),
                );
            }
            let roles = self.methods.role_ids(self.slots.len())?;
            let mut seen = std::collections::BTreeSet::new();
            for edit in &self.substitutions {
                if !roles.iter().any(|role| role.role == edit.role)
                    || !token(&edit.role, 128)
                    || edit.offset >= collagene_core::model::MAX_SEQUENCE_LEN
                    || edit.expected.len() != 1
                    || edit.replacement.len() != 1
                    || !edit
                        .expected
                        .bytes()
                        .chain(edit.replacement.bytes())
                        .all(|base| b"ACGT".contains(&base))
                    || edit.expected == edit.replacement
                    || !seen.insert((&edit.role, edit.offset))
                {
                    return Err("Internal substitutions require distinct concrete single-base edits on stable source roles.".into());
                }
            }
        }
        self.validate_internal_edits()?;

        let has_origin = self.methods.rows.iter().any(|row| row.run_origin.is_some());

        if self.alternatives.as_ref().is_some_and(|library| {
            has_origin
                || library
                    .bins
                    .iter()
                    .flat_map(|bin| &bin.choices)
                    .any(|choice| choice.row.run_origin.is_some())
                || library.rows.iter().any(|row| {
                    row.edits
                        .methods
                        .rows
                        .iter()
                        .any(|row| row.run_origin.is_some())
                })
        }) {
            return Err(
                "Saved-run material is currently supported only in ordinary assemblies.".into(),
            );
        }

        if self.alternatives.is_some() && self.uses_paired_fields() {
            return Err(
                "Paired-oligo sources/settings are not supported in assembly libraries.".into(),
            );
        }
        if self.alternatives.as_ref().is_some_and(|library| {
            self.methods.mode == Mode::Tailing
                || library
                    .rows
                    .iter()
                    .any(|row| row.edits.methods.mode == Mode::Tailing)
        }) {
            return Err("TA/GC cloning is an ordinary two-source workflow; library preparations are not supported.".into());
        }

        validate_controls(
            &self.name,
            self.circular,
            self.active_slot,
            self.slots.len(),
            &self.methods,
            &self.junction_edits,
        )?;
        bounded_size(self)?;
        let models = self
            .slots
            .iter()
            .flat_map(|slot| [slot.captured.as_deref(), slot.anneal_bottom.as_deref()])
            .flatten()
            .map(|captured| &captured.molecule)
            .collect::<Vec<_>>();
        collagene_core::model::bounds::model_preflight(&models).map_err(|e| e.to_string())?;
        let mut bases = 0usize;
        let mut next_bottom = self.slots.len();
        for (index, (slot, row)) in self.slots.iter().zip(&self.methods.rows).enumerate() {
            for captured in [&slot.captured, &slot.anneal_bottom].into_iter().flatten() {
                bases = bases
                    .checked_add(captured.molecule.sequence.len())
                    .filter(|n| *n <= 8 << 20)
                    .ok_or("Assembly pinned source bases exceed8MiB")?;
            }
            validate_slot(slot, row, index)?;
            if row.preparation == OverlapPreparation::AnnealedOligos {
                let expected = slot.anneal_bottom.as_ref().map(|_| {
                    let index = next_bottom;
                    next_bottom += 1;
                    index
                });
                if row.anneal_bottom_document != expected
                    || (row.anneal_candidate.is_some()
                        && (self.anneal_context(index).is_none()
                            || row.anneal_source_context != self.anneal_context(index)))
                {
                    return Err("Annealing candidate/source binding changed; capture and inspect the oligos again.".into());
                }
            }
        }
        if let Some(library) = &self.alternatives {
            if library.bins.len() != self.slots.len() {
                return Err(
                    "Working assembly slots differ from the retained alternative bins.".into(),
                );
            }
            library.validate()?;
            if self.methods.source_roles != crate::assembly_alternatives::roles_for_active(self)? {
                return Err("Working preparation identities differ from the active construct and source choices.".into());
            }
        }
        Ok(())
    }
    pub(crate) fn sync_rows(&mut self) {
        self.methods.circular = self.circular;
        self.methods
            .rows
            .resize_with(self.slots.len(), OverlapInput::default);
        for (i, (slot, row)) in self.slots.iter().zip(&mut self.methods.rows).enumerate() {
            if row.run_origin.as_ref().is_some_and(|origin| {
                slot.captured.as_ref().is_none_or(|capture| {
                    *origin.event() != capture.source.binding
                        || origin.event().revision != capture.molecule.revision()
                })
            }) {
                row.run_origin = None;
            }
            row.document = i;
            match &slot.selector {
                assembly_intent::Selector::ReferenceRange {
                    start,
                    end,
                    reverse,
                } => {
                    row.reverse = *reverse;
                    row.interval_start = start.saturating_add(1).to_string();
                    row.interval_end = if slot.captured.is_some() {
                        end.to_string()
                    } else {
                        String::new()
                    };
                }
                assembly_intent::Selector::MaterialStrandRange { .. } => {
                    row.interval_start.clear();
                    row.interval_end.clear();
                }
            }
        }

        let mut next_bottom = self.slots.len();
        for index in 0..self.slots.len() {
            let context = self.anneal_context(index);
            let row = &mut self.methods.rows[index];
            if row.anneal_source_context != context {
                row.anneal_candidate = None;
            }
            row.anneal_source_context = context;
            row.anneal_bottom_document = if row.preparation == OverlapPreparation::AnnealedOligos
                && self.slots[index].anneal_bottom.is_some()
            {
                let at = next_bottom;
                next_bottom += 1;
                Some(at)
            } else {
                None
            };
            if row.preparation != OverlapPreparation::AnnealedOligos {
                row.anneal_candidate = None;
            }
        }
        self.version = CURRENT_VERSION;
    }
}

#[cfg(test)]
pub fn encode(draft: &Draft) -> Result<Vec<u8>, String> {
    let stored = crate::assembly_alternatives::for_storage(draft)?;
    encode_raw(&stored)
}
/// Encode current editable state without promoting working fields into saved rows.
pub fn encode_raw(draft: &Draft) -> Result<Vec<u8>, String> {
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

pub(crate) fn digest_selector(
    source: &SourceInfo,
    geometry: &DigestGeometry,
    reverse: bool,
) -> Result<assembly_intent::Selector, String> {
    let top_length = source.length;
    if top_length == 0 {
        return Err("Digest source is empty".into());
    }
    if let Some(domain) = &geometry.material {
        if domain.top_length != source.length || domain.topology != source.topology {
            return Err("Actual source domain differs from the pinned reference".into());
        }
        return assembly_intent::digest_material_selector(
            domain,
            &assembly_intent::DigestGeometry {
                top_start: geometry.top_start,
                top_length: geometry.top_length,
                bottom_start: geometry.bottom_start,
                bottom_length: geometry.bottom_length,
            },
            reverse,
        )
        .map_err(|e| e.to_string());
    }
    let (start, length) = if reverse {
        let start = if source.topology == Topology::Circular {
            usize::try_from(geometry.bottom_start.rem_euclid(top_length as i64))
                .map_err(|_| "Invalid circular coordinate")?
        } else {
            usize::try_from(geometry.bottom_start)
                .map_err(|_| "Bottom interval starts outside reference")?
        };
        (start, geometry.bottom_length)
    } else {
        (geometry.top_start, geometry.top_length)
    };
    let end = start
        .checked_add(length)
        .ok_or("Reference interval overflow")?;
    if start >= top_length
        || length == 0
        || length > top_length
        || (source.topology == Topology::Linear && end > top_length)
    {
        return Err("Digest interval is outside the reference".into());
    }
    Ok(assembly_intent::Selector::ReferenceRange {
        start,
        end,
        reverse,
    })
}

#[cfg(test)]
mod tailing_tests {
    use super::*;
    use crate::{assembly_role_ids::Role, cloning_methods_model::Treatment};
    use collagene_core::{tailing::TailingBase, verification::PlanInput, workspace::Workspace};

    fn captured(name: &str, dna: &str) -> Arc<CapturedSource> {
        let molecule = Molecule::dna(name, dna, Topology::Linear).unwrap();
        let workspace = Workspace::from_documents(name, vec![molecule.clone()], 0).unwrap();
        let document = &workspace.documents[0];
        Arc::new(CapturedSource {
            source: SourceInfo {
                binding: collagene_core::run_lineage::EventBinding {
                    document_id: document.id.clone(),
                    event_id: document.head_event.clone().unwrap(),
                    revision: document.head.clone(),
                },
                name: molecule.name.clone(),
                length: molecule.sequence.len(),
                alphabet: molecule.sequence.alphabet(),
                topology: molecule.topology,
            },
            molecule,
        })
    }
    fn draft() -> Draft {
        let mut draft = Draft::default();
        draft.methods.mode = Mode::Tailing;
        draft.methods.source_roles = Some(vec![
            Role {
                input_id: "pinned-vector".into(),
                role: "vector-core".into(),
            },
            Role {
                input_id: "pinned-insert".into(),
                role: "insert-core".into(),
            },
        ]);
        for (slot, source) in draft.slots.iter_mut().zip([
            captured("Original vector", "AACC"),
            captured("Original insert", "GGTT"),
        ]) {
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    end: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = source.molecule.sequence.len();
            }
            slot.captured = Some(source);
        }
        draft.sync_rows();
        draft
    }

    #[test]
    fn two_source_tailing_roundtrip_preserves_pins_roles_preparations_and_literal_intent() {
        let mut draft = draft();
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange {
                reverse: value, ..
            } = &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = true;
        }
        draft.methods.rows[0].preparation = OverlapPreparation::Digest;
        draft.methods.rows[0].enzymes = "EcoRV".into();
        draft.methods.rows[0].fragment = "2".into();
        draft.methods.rows[1].forward = "GG".into();
        draft.methods.rows[1].reverse_primer = "AA".into();
        draft.methods.tailing_profile = collagene_core::tailing_recipe::Profile::Gc;
        draft.methods.vector.action = collagene_core::tailing_recipe::Action::AddTail {
            base: TailingBase::C,
        };
        draft.methods.insert.action = collagene_core::tailing_recipe::Action::AddTail {
            base: TailingBase::G,
        };
        draft.methods.insert.treatment = Treatment::Phosphorylate;
        draft.junction_edits.push(assembly_junction_edit::Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "T".into(),
        });
        draft.sync_rows();
        let bytes = encode_raw(&draft).unwrap();
        let restored = decode(&bytes).unwrap();
        assert!(restored == draft);
        assert_eq!(restored.methods.rows[0].document, 0);
        assert_eq!(restored.methods.rows[1].document, 1);
        assert!(restored.methods.rows[1].reverse);
        assert_eq!(restored.junction_edits[0].inserted, "T");
        assert_eq!(restored.methods.rows[0].enzymes, "EcoRV");
    }

    #[test]
    fn linear_terminal_intent_is_tailing_only_and_has_no_downstream_trim() {
        let mut draft = draft();
        draft.circular = false;
        draft.junction_edits = vec![assembly_junction_edit::Edit {
            junction: 1,
            remove_left: 1,
            remove_right: 0,
            inserted: "A".into(),
        }];
        draft.sync_rows();
        draft.validate().unwrap();
        assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
        draft.junction_edits[0].remove_right = 1;
        assert!(draft.validate().is_err());
        draft.junction_edits[0].remove_right = 0;
        for mode in [Mode::Overlap, Mode::GoldenGate, Mode::MultiSiteReference] {
            draft.methods.mode = mode;
            assert!(draft.validate().is_err());
        }
    }

    #[test]
    fn moving_tailing_sources_keeps_each_sources_tail_chemistry_and_preparation() {
        let mut draft = draft();
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange {
                reverse: value, ..
            } = &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = true;
        }
        draft.methods.rows[0].preparation = OverlapPreparation::Digest;
        draft.methods.rows[0].enzymes = "EcoRV".into();
        draft.methods.rows[1].forward = "GG".into();
        draft.methods.tailing_profile = collagene_core::tailing_recipe::Profile::Gc;
        draft.methods.vector.action = collagene_core::tailing_recipe::Action::AddTail {
            base: TailingBase::C,
        };
        draft.methods.vector.treatment = Treatment::Dephosphorylate;
        draft.methods.insert.action = collagene_core::tailing_recipe::Action::AddTail {
            base: TailingBase::G,
        };
        draft.methods.insert.treatment = Treatment::Phosphorylate;
        draft.sync_rows();
        let original = draft.clone();
        draft.slots.swap(0, 1);
        draft.methods.swap_rows(0, 1);
        draft.sync_rows();
        draft.validate().unwrap();
        assert!(draft.slots[0].captured == original.slots[1].captured);
        assert!(draft.slots[1].captured == original.slots[0].captured);
        assert!(
            draft.methods.vector.action
                == (collagene_core::tailing_recipe::Action::AddTail {
                    base: TailingBase::G
                })
        );
        assert!(draft.methods.vector.treatment == Treatment::Phosphorylate);
        assert!(
            draft.methods.insert.action
                == (collagene_core::tailing_recipe::Action::AddTail {
                    base: TailingBase::C
                })
        );
        assert!(draft.methods.insert.treatment == Treatment::Dephosphorylate);
        assert_eq!(draft.methods.rows[0].forward, "GG");
        assert_eq!(draft.methods.rows[1].enzymes, "EcoRV");
        assert!(!draft.methods.rows[0].reverse && draft.methods.rows[1].reverse);
        assert!(
            draft.methods.source_roles.as_ref().unwrap()[0]
                == original.methods.source_roles.as_ref().unwrap()[1]
        );
        draft.methods.swap_rows(0, 0);
        assert!(
            draft.methods.vector.action
                == (collagene_core::tailing_recipe::Action::AddTail {
                    base: TailingBase::G
                })
        );
        draft.slots.swap(1, 0);
        draft.methods.swap_rows(1, 0);
        draft.sync_rows();
        assert!(draft == original);
    }

    #[test]
    fn kept_ends_and_added_tail_choices_follow_exact_sources_on_reorder() {
        use collagene_core::tailing_recipe::{Action, Profile};
        let mut draft = draft();
        draft.methods.tailing_profile = Profile::Ta;
        draft.methods.vector.action = Action::KeepExisting;
        draft.methods.vector.treatment = Treatment::Dephosphorylate;
        draft.methods.insert.action = Action::AddTail {
            base: TailingBase::A,
        };
        draft.methods.insert.treatment = Treatment::Phosphorylate;
        let before = draft.clone();
        draft.slots.swap(0, 1);
        draft.methods.swap_rows(0, 1);
        draft.sync_rows();
        assert!(draft.methods.vector.action == before.methods.insert.action);
        assert!(draft.methods.insert.action == Action::KeepExisting);
        assert!(draft.methods.vector.treatment == Treatment::Phosphorylate);
        assert!(draft.methods.insert.treatment == Treatment::Dephosphorylate);
        assert!(draft.methods.tailing_profile == Profile::Ta);
        assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
        draft.slots.swap(0, 1);
        draft.methods.swap_rows(0, 1);
        draft.sync_rows();
        assert!(draft == before);
    }

    #[test]
    fn tailing_rejects_wrong_count_source_drift_unsupported_schema_and_library_rows() {
        let mut invalid = draft();
        invalid.methods.rows[1].document = 0;
        assert!(invalid.validate().is_err());
        invalid.sync_rows();
        invalid.methods.rows[1].reverse = true;
        assert!(invalid.validate().is_err());
        invalid.sync_rows();
        invalid.version = CURRENT_VERSION - 1;
        assert!(invalid.validate().is_err());
        invalid.version = CURRENT_VERSION;
        invalid.slots.push(Slot::empty("Third source"));
        invalid.sync_rows();
        assert!(invalid.validate().is_err());

        let mut library =
            crate::assembly_alternatives::enable(&Draft::default(), "ta-exclusion").unwrap();
        let saved = &mut library.alternatives.as_mut().unwrap().rows[0].edits.methods;
        saved.mode = Mode::Tailing;
        assert!(library.validate().unwrap_err().contains("library"));
    }

    #[test]
    fn circular_target_uses_only_operator_declared_junction_bases() {
        let mut draft = draft();
        draft.junction_edits = vec![
            assembly_junction_edit::Edit {
                junction: 0,
                remove_left: 0,
                remove_right: 0,
                inserted: "T".into(),
            },
            assembly_junction_edit::Edit {
                junction: 1,
                remove_left: 0,
                remove_right: 0,
                inserted: "A".into(),
            },
        ];
        draft.validate().unwrap();
        let roles = draft.methods.role_ids(2).unwrap();
        let inputs = draft
            .slots
            .iter()
            .zip(&roles)
            .map(|(slot, role)| {
                let captured = slot.captured.as_ref().unwrap();
                PlanInput {
                    id: role.input_id.clone(),
                    revision: captured.source.binding.revision.clone(),
                    molecule: captured.molecule.clone(),
                }
            })
            .collect::<Vec<_>>();
        let base = collagene_core::assembly_intent::Request {
            material_bindings: Vec::new(),
            name: draft.name.clone(),
            circular: true,
            fragments: draft
                .slots
                .iter()
                .zip(&roles)
                .map(
                    |(slot, role)| collagene_core::assembly_intent::Contribution {
                        role: role.role.clone(),
                        input_id: role.input_id.clone(),
                        selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                            start: slot.selector.reference_interval().unwrap().0,
                            end: slot.selector.reference_interval().unwrap().1,
                            reverse: slot.selector.reference_interval().unwrap().2,
                        },
                    },
                )
                .collect(),
        };
        let target =
            assembly_junction_edit::build_intent(&inputs, &base, &draft.junction_edits, None)
                .unwrap();
        assert_eq!(target.target.sequence.as_str(), "AACCTGGTTA");
        draft.methods.tailing_profile = collagene_core::tailing_recipe::Profile::Gc;
        draft.methods.vector.action = collagene_core::tailing_recipe::Action::AddTail {
            base: TailingBase::C,
        };
        draft.methods.insert.action = collagene_core::tailing_recipe::Action::AddTail {
            base: TailingBase::G,
        };
        draft.sync_rows();
        let unchanged =
            assembly_junction_edit::build_intent(&inputs, &base, &draft.junction_edits, None)
                .unwrap();
        assert_eq!(target.context_revision, unchanged.context_revision);
        assert_eq!(target.target, unchanged.target);
    }
}

#[cfg(test)]
mod paired_oligo_tests {
    use super::*;
    use crate::assembly_role_ids::Role;
    use collagene_core::workspace::Workspace;

    fn captured(name: &str, dna: &str) -> Arc<CapturedSource> {
        captured_topology(name, dna, Topology::Linear)
    }
    fn captured_topology(name: &str, dna: &str, topology: Topology) -> Arc<CapturedSource> {
        let molecule = Molecule::dna(name, dna, topology).unwrap();
        let workspace = Workspace::from_documents(name, vec![molecule.clone()], 0).unwrap();
        let document = &workspace.documents[0];
        Arc::new(CapturedSource {
            source: SourceInfo {
                binding: collagene_core::run_lineage::EventBinding {
                    document_id: document.id.clone(),
                    event_id: document.head_event.clone().unwrap(),
                    revision: document.head.clone(),
                },
                name: molecule.name.clone(),
                length: molecule.sequence.len(),
                alphabet: Alphabet::Dna,
                topology,
            },
            molecule,
        })
    }
    fn paired() -> Draft {
        let mut draft = Draft::default();
        draft.methods.mode = Mode::Ligation;
        draft.methods.source_roles = Some(vec![
            Role {
                input_id: "backbone-parent".into(),
                role: "backbone".into(),
            },
            Role {
                input_id: "adaptor-top".into(),
                role: "adaptor".into(),
            },
        ]);
        for (slot, source) in draft.slots.iter_mut().zip([
            captured("Backbone", "AAACCC"),
            captured("Top oligo", "AAGGCC"),
        ]) {
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    end: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = source.molecule.sequence.len();
            }
            slot.captured = Some(source);
        }
        draft.methods.rows[1].preparation = OverlapPreparation::AnnealedOligos;
        draft.methods.rows[1].anneal_options.min_paired_bases = 4;
        draft.slots[1].anneal_bottom = Some(captured("Bottom oligo", "GGCCTT"));
        draft.sync_rows();
        draft.methods.rows[1].anneal_candidate = Some(0);
        draft.validate().unwrap();
        draft
    }

    #[test]
    fn ligation_intent_literals_and_linear_suffixes_are_independent_of_preparation() {
        let mut draft = paired();
        let preparation = draft.methods.clone();
        draft.junction_edits = vec![
            assembly_junction_edit::Edit {
                junction: 0,
                remove_left: 0,
                remove_right: 0,
                inserted: "G".into(),
            },
            assembly_junction_edit::Edit {
                junction: 1,
                remove_left: 0,
                remove_right: 0,
                inserted: "C".into(),
            },
        ];
        draft.validate().unwrap();
        assert!(draft.methods == preparation);
        assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
        draft.circular = false;
        draft.sync_rows();
        draft.validate().unwrap();
        draft.junction_edits[1].remove_right = 1;
        assert!(draft.validate().is_err());
        draft.junction_edits[1].remove_right = 0;
        draft.junction_edits[1].inserted = "A".repeat(201);
        assert!(draft.validate().is_err());
    }

    #[test]
    fn swapped_ordered_roles_reset_intent_to_the_literal_old_bottom_and_clear_stale_pins() {
        use collagene_core::{
            chemistry::{FivePrime, Terminus, ThreePrime},
            verification::PlanInput,
        };
        let mut draft = paired();
        draft.slots[0].captured = Some(captured_topology(
            "Backbone",
            "GGAATTCTTACG",
            Topology::Circular,
        ));
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { start: value, .. } =
                &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 2;
        }
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 14;
        }
        draft.slots[1].captured = Some(captured("Top", "AATTGCTACG"));
        draft.slots[1].anneal_bottom = Some(captured("Bottom", "AATTCGTAGC"));
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { start: value, .. } =
                &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 1;
        }
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 9;
        }
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange {
                reverse: value, ..
            } = &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = true;
        }
        draft.slots[1].digest_geometry = Some(DigestGeometry {
            material: None,
            top_start: 1,
            top_length: 8,
            bottom_start: 1,
            bottom_length: 8,
        });
        draft.methods.rows[1].anneal_options.min_paired_bases = 6;
        draft.methods.rows[1].anneal_options.top_chemistry = Terminus {
            five_prime: FivePrime::Hydroxyl,
            three_prime: ThreePrime::Blocked,
        };
        draft.methods.rows[1].anneal_options.bottom_chemistry = Terminus::CLEAVED;
        draft.methods.rows[1].treatment = crate::cloning_methods_model::Treatment::Phosphorylate;
        draft.junction_edits = vec![assembly_junction_edit::Edit {
            junction: 0,
            remove_left: 1,
            remove_right: 0,
            inserted: "T".into(),
        }];
        draft.sync_rows();
        draft.methods.rows[1].anneal_candidate = Some(0);
        draft.validate().unwrap();
        let before = draft.clone();
        let swapped = draft.swap_oligo_roles(1).unwrap();
        assert!(draft == before);
        assert!(swapped.slots[1].captured == before.slots[1].anneal_bottom);
        assert!(swapped.slots[1].anneal_bottom == before.slots[1].captured);
        assert_eq!(
            (
                swapped.slots[1].selector.reference_interval().unwrap().0,
                swapped.slots[1].selector.reference_interval().unwrap().1,
                swapped.slots[1].selector.reference_interval().unwrap().2
            ),
            (0, 10, false)
        );
        assert_eq!(
            swapped.methods.rows[1].anneal_options.top_chemistry,
            Terminus::CLEAVED
        );
        assert_eq!(
            swapped.methods.rows[1].anneal_options.bottom_chemistry,
            before.methods.rows[1].anneal_options.top_chemistry
        );
        assert!(swapped.methods.rows[1].treatment == before.methods.rows[1].treatment);
        assert!(swapped.methods.rows[0] == before.methods.rows[0]);
        assert_eq!(swapped.methods.source_roles, before.methods.source_roles);
        assert_eq!(swapped.methods.rows[1].anneal_candidate, None);
        assert_ne!(
            swapped.methods.rows[1].anneal_source_context,
            before.methods.rows[1].anneal_source_context
        );
        assert!(swapped.slots[1].digest_geometry.is_none() && swapped.junction_edits.is_empty());
        assert!(decode(&encode_raw(&swapped).unwrap()).unwrap() == swapped);
        let ids = swapped.methods.all_input_ids().unwrap();
        let sources = [
            swapped.slots[0].captured.as_ref().unwrap(),
            swapped.slots[1].captured.as_ref().unwrap(),
            swapped.slots[1].anneal_bottom.as_ref().unwrap(),
        ];
        let inputs = sources
            .into_iter()
            .zip(&ids)
            .map(|(capture, id)| PlanInput {
                id: id.clone(),
                revision: capture.source.binding.revision.clone(),
                molecule: capture.molecule.clone(),
            })
            .collect::<Vec<_>>();
        let roles = swapped.methods.role_ids(2).unwrap();
        let request = collagene_core::assembly_intent::Request {
            material_bindings: Vec::new(),
            name: "Independently reversed adaptor".into(),
            circular: true,
            fragments: swapped
                .slots
                .iter()
                .zip(&roles)
                .map(
                    |(slot, role)| collagene_core::assembly_intent::Contribution {
                        role: role.role.clone(),
                        input_id: role.input_id.clone(),
                        selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                            start: slot.selector.reference_interval().unwrap().0,
                            end: slot.selector.reference_interval().unwrap().1,
                            reverse: slot.selector.reference_interval().unwrap().2,
                        },
                    },
                )
                .collect(),
        };
        let intent = collagene_core::assembly_intent::build(&inputs, &request, None).unwrap();
        assert_eq!(intent.target.sequence.as_str(), "AATTCTTACGGGAATTCGTAGC");
    }

    #[test]
    fn swapping_one_pair_preserves_unrelated_four_slot_junction_intent() {
        let mut draft = paired();
        for (label, dna, id) in [
            ("Third", "GTAC", "third-parent"),
            ("Fourth", "CAAT", "fourth-parent"),
        ] {
            let source = captured(label, dna);
            let mut slot = Slot::empty(label);
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    end: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = source.molecule.sequence.len();
            }
            slot.captured = Some(source);
            draft.slots.push(slot);
            draft.methods.rows.push(OverlapInput::default());
            draft.methods.source_roles.as_mut().unwrap().push(Role {
                input_id: id.into(),
                role: label.to_lowercase(),
            });
        }
        draft.methods.rows[2].preparation = OverlapPreparation::ExactPcr;
        draft.methods.rows[2].forward = "GT".into();
        draft.methods.rows[2].reverse_primer = "GT".into();
        draft.junction_edits = vec![
            assembly_junction_edit::Edit {
                junction: 0,
                remove_left: 0,
                remove_right: 0,
                inserted: "A".into(),
            },
            assembly_junction_edit::Edit {
                junction: 1,
                remove_left: 0,
                remove_right: 0,
                inserted: "C".into(),
            },
            assembly_junction_edit::Edit {
                junction: 2,
                remove_left: 1,
                remove_right: 1,
                inserted: "GG".into(),
            },
            assembly_junction_edit::Edit {
                junction: 3,
                remove_left: 1,
                remove_right: 0,
                inserted: "T".into(),
            },
        ];
        draft.sync_rows();
        draft.methods.rows[1].anneal_candidate = Some(0);
        let before = draft.clone();
        let swapped = draft.swap_oligo_roles(1).unwrap();
        assert_eq!(swapped.junction_edits, before.junction_edits[2..]);
        assert!(swapped.slots[0] == before.slots[0] && swapped.slots[2..] == before.slots[2..]);
        assert!(
            swapped.methods.rows[0] == before.methods.rows[0]
                && swapped.methods.rows[2..] == before.methods.rows[2..]
        );
        assert_eq!(swapped.methods.source_roles, before.methods.source_roles);
        assert_eq!(swapped.methods.rows[1].anneal_candidate, None);
        assert!(swapped.slots[1].captured == before.slots[1].anneal_bottom);
        assert!(decode(&encode_raw(&swapped).unwrap()).unwrap() == swapped);
        assert!(draft == before);
        // The first circular contribution has the closing boundary as its
        // incoming neighbor; it must not erase the middle two declarations.
        draft.slots.swap(0, 1);
        draft.methods.swap_rows(0, 1);
        draft.sync_rows();
        let first = draft.swap_oligo_roles(0).unwrap();
        assert_eq!(first.junction_edits, draft.junction_edits[1..3]);
    }

    #[test]
    fn swapping_originals_and_full_termini_matches_asymmetric_physical_reversal() {
        use collagene_core::chemistry::{self, FivePrime, Terminus, ThreePrime};
        let mut draft = paired();
        draft.slots[1].captured = Some(captured("Top", "AAGCTACG"));
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 8;
        }
        draft.slots[1].anneal_bottom = Some(captured("Bottom", "TTTCGTAGC"));
        draft.methods.rows[1].anneal_options.min_paired_bases = 6;
        draft.methods.rows[1].anneal_options.top_chemistry = Terminus {
            five_prime: FivePrime::Unknown,
            three_prime: ThreePrime::Blocked,
        };
        draft.methods.rows[1].anneal_options.bottom_chemistry = Terminus::CLEAVED;
        draft.sync_rows();
        let original = chemistry::anneal_oligos(
            &draft.slots[1].captured.as_ref().unwrap().molecule.sequence,
            &draft.slots[1]
                .anneal_bottom
                .as_ref()
                .unwrap()
                .molecule
                .sequence,
            &draft.methods.rows[1].anneal_options,
        )
        .unwrap();
        let swapped = draft.swap_oligo_roles(1).unwrap();
        let fresh = chemistry::anneal_oligos(
            &swapped.slots[1]
                .captured
                .as_ref()
                .unwrap()
                .molecule
                .sequence,
            &swapped.slots[1]
                .anneal_bottom
                .as_ref()
                .unwrap()
                .molecule
                .sequence,
            &swapped.methods.rows[1].anneal_options,
        )
        .unwrap();
        assert!(original.complete && fresh.complete && !original.candidates.is_empty());
        assert_eq!(original.candidates.len(), fresh.candidates.len());
        for candidate in &original.candidates {
            let reversed = candidate.duplex.reverse().unwrap();
            assert!(
                fresh
                    .candidates
                    .iter()
                    .any(|candidate| candidate.duplex == reversed)
            );
        }
        let asymmetric = original
            .candidates
            .iter()
            .find(|candidate| candidate.duplex.bottom_start == 2)
            .unwrap();
        let (left, right) = asymmetric.duplex.ends().unwrap();
        assert_eq!(left.overhang.len(), 2);
        assert_eq!(right.overhang.len(), 3);
        assert!(asymmetric.duplex.nicks.is_empty());
        let before = serde_json::to_value(&draft).unwrap();
        assert!(draft.swap_oligo_roles(0).is_err());
        assert!(draft.swap_oligo_roles(99).is_err());
        assert_eq!(serde_json::to_value(&draft).unwrap(), before);
    }

    #[test]
    fn paired_draft_flattens_three_exact_originals_and_roundtrips_current_schema() {
        let draft = paired();
        assert_eq!(draft.version, CURRENT_VERSION);
        let expected = format!("anneal-bottom-{}", blake3::hash(b"adaptor-top").to_hex());
        assert_eq!(draft.methods.bottom_input_id(1).unwrap(), expected);
        assert_eq!(
            draft.methods.all_input_ids().unwrap(),
            vec![
                "backbone-parent".to_string(),
                "adaptor-top".into(),
                expected
            ]
        );
        assert_eq!(draft.methods.source_indices().unwrap(), vec![0, 1, 2]);
        assert_eq!(
            crate::assembly_alternatives::all_captures(&draft).count(),
            3,
            "Workspace restore must validate the bottom event as well as primary contributions"
        );
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        assert!(restored == draft);
        assert_eq!(restored.methods.rows[1].anneal_candidate, Some(0));
        let mut legacy = draft.clone();
        legacy.version = CURRENT_VERSION - 1;
        assert!(legacy.validate().is_err());
        assert!(
            crate::assembly_alternatives::for_storage(&draft)
                .unwrap()
                .version
                == CURRENT_VERSION
        );
    }

    #[test]
    fn moving_source_slots_keeps_pair_events_options_candidate_and_stable_identity() {
        let mut draft = paired();
        let before = draft.clone();
        let bottom_id = draft.methods.bottom_input_id(1).unwrap();
        draft.slots.swap(0, 1);
        draft.methods.swap_rows(0, 1);
        draft.sync_rows();
        draft.validate().unwrap();
        assert!(draft.slots[0].anneal_bottom == before.slots[1].anneal_bottom);
        assert_eq!(draft.methods.bottom_input_id(0).unwrap(), bottom_id);
        assert_eq!(draft.methods.rows[0].anneal_candidate, Some(0));
        assert_eq!(
            draft.methods.rows[0].anneal_source_context,
            before.methods.rows[1].anneal_source_context
        );
        assert_eq!(draft.methods.source_indices().unwrap(), vec![0, 1, 2]);
        draft.slots.swap(1, 0);
        draft.methods.swap_rows(1, 0);
        draft.sync_rows();
        assert!(draft == before);
    }

    #[test]
    fn equal_dna_different_bottom_event_and_changed_options_clear_candidate() {
        let mut draft = paired();
        let prior = draft.methods.rows[1].anneal_source_context.clone();
        let bottom = Arc::make_mut(draft.slots[1].anneal_bottom.as_mut().unwrap());
        bottom.source.binding.event_id = "another-original-bottom-event".into();
        assert!(draft.validate().is_err());
        draft.sync_rows();
        assert_eq!(draft.methods.rows[1].anneal_candidate, None);
        assert_ne!(draft.methods.rows[1].anneal_source_context, prior);
        draft.methods.rows[1].anneal_candidate = Some(0);
        draft.methods.rows[1].anneal_options.min_paired_bases = 5;
        draft.sync_rows();
        assert_eq!(draft.methods.rows[1].anneal_candidate, None);
        draft.slots[1].anneal_bottom = None;
        draft.sync_rows();
        assert!(draft.methods.source_indices().is_err());
        assert_eq!(draft.methods.rows[1].anneal_bottom_document, None);
    }

    #[test]
    fn collisions_non_ligation_libraries_and_unbounded_oligos_are_rejected() {
        let mut draft = paired();
        let secondary = draft.methods.bottom_input_id(1).unwrap();
        draft.methods.source_roles.as_mut().unwrap()[0].input_id = secondary;
        assert!(draft.methods.all_input_ids().is_err());
        let mut wrong = paired();
        wrong.methods.mode = Mode::Overlap;
        assert!(wrong.validate().is_err());
        let mut library =
            crate::assembly_alternatives::enable(&Draft::default(), "pair-ban").unwrap();
        library.slots[0].anneal_bottom = paired().slots[1].anneal_bottom.clone();
        library.sync_rows();
        assert!(library.validate().unwrap_err().contains("libraries"));
        let mut too_long = paired();
        too_long.slots[1].anneal_bottom = Some(captured(
            "Oversized oligo",
            &"A".repeat(collagene_core::chemistry::MAX_ANNEAL_OLIGO_BASES + 1),
        ));
        too_long.sync_rows();
        assert!(too_long.validate().is_err());
    }
}

#[cfg(test)]
mod run_origin_tests {
    use super::*;
    use collagene_core::{
        analysis::Primer,
        cloning_outputs,
        run_lineage::{EventBinding, InputEventBinding},
        run_material,
        verification::{DesignPlan, Operation, PlanInput, Source},
        workspace::Workspace,
    };

    fn saved_pcr() -> Draft {
        saved_pcr_workspace().0
    }
    fn saved_pcr_workspace() -> (Draft, Workspace) {
        let source = Molecule::dna("Original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
        let target = Molecule::dna("Independent", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
        let workspace =
            Workspace::from_documents("origin", vec![source.clone(), target.clone()], 0).unwrap();
        let event = |index: usize| {
            let doc = &workspace.documents[index];
            EventBinding {
                document_id: doc.id.clone(),
                event_id: doc.head_event.clone().unwrap(),
                revision: doc.head.clone(),
            }
        };
        let plan = DesignPlan {
            name: "Original PCR".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "original".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps: vec![Operation::Pcr {
                template: Source::Input("original".into()),
                forward: Primer::new("Original forward", "ACGACCCC", "").unwrap(),
                reverse: Primer::new("Original reverse", "GTTACCCC", "").unwrap(),
            }],
        };
        let commit = cloning_outputs::finalize(
            &workspace,
            &plan,
            event(1),
            vec![InputEventBinding {
                input_id: "original".into(),
                event: event(0),
            }],
            &cloning_outputs::Options {
                name: "Actual saved PCR".into(),
                save_primers: false,
                show_prepared: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let selection = run_material::Selection {
            run_id: commit.receipt.run.run_id.clone(),
            step_index: 0,
            event: commit.receipt.construct.clone(),
        };
        let molecule = commit
            .workspace
            .head(&selection.event.document_id)
            .unwrap()
            .clone();
        let captured = Arc::new(CapturedSource {
            source: SourceInfo {
                binding: selection.event.clone(),
                name: molecule.name.clone(),
                length: molecule.sequence.len(),
                alphabet: Alphabet::Dna,
                topology: molecule.topology,
            },
            molecule,
        });
        let mut draft = Draft::default();
        draft.methods.mode = Mode::Ligation;
        draft.circular = false;
        draft.slots.truncate(1);
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = captured.molecule.sequence.len();
        }
        draft.slots[0].captured = Some(captured);
        draft.methods.rows.truncate(1);
        draft.methods.rows[0].preparation = OverlapPreparation::WholeLinear;
        draft.methods.rows[0].run_origin = Some(
            crate::run_origin_model::PinnedOrigin::from_material(
                collagene_core::run_material::MaterialOrigin::Run(
                    run_material::capsule(&commit.workspace, &selection).unwrap(),
                ),
            )
            .unwrap(),
        );
        draft.sync_rows();
        draft.validate().unwrap();
        (draft, commit.workspace)
    }

    #[test]
    fn saved_origin_preserves_inactive_original_generation_text_without_redeclaration() {
        let mut draft = saved_pcr();
        draft.methods.rows[0].source_generation =
            Some(crate::cloning_methods_model::SourceGenerationInput {
                dam: "unfinished".into(),
                dcm: "fully".into(),
                cpg: "unknown".into(),
            });
        draft.validate().unwrap();
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        assert!(restored == draft);
        let (_, rows) = crate::cloning_preparations::prepared_rows(&restored.methods, 1).unwrap();
        assert!(rows[0].generation.is_none());
        assert!(matches!(
            rows[0].preparation,
            collagene_core::assembly_recipe::Preparation::FromRun { .. }
        ));
        let mut original = restored;
        original.methods.rows[0].run_origin = None;
        assert!(crate::cloning_preparations::prepared_rows(&original.methods, 1).is_err());
    }

    #[test]
    fn derived_digest_roundtrip_retains_actual_preparation_and_refuses_run_only_routes() {
        let (original, mut workspace) = saved_pcr_workspace();
        let run = original.methods.rows[0]
            .run_origin
            .as_ref()
            .unwrap()
            .as_run()
            .unwrap()
            .clone();
        let source = run.selection.event.clone();
        workspace.active_document = source.document_id.clone();
        let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
        let enzyme = catalog
            .edition()
            .enzymes
            .iter()
            .find(|e| e.name == "SmaI")
            .unwrap();
        let selection =
            collagene_core::catalog_digest::select(&catalog, &[enzyme.id.clone()]).unwrap();
        let installed = collagene_core::material_derivation::prepare_install(
            &workspace,
            &collagene_core::material_derivation::Request {
                source: source.clone(),
                origin: collagene_core::material_derivation::SourceOrigin::Saved(
                    collagene_core::run_material::MaterialOrigin::Run(run),
                ),
                derivation: collagene_core::material_derivation::Derivation::CatalogDigest {
                    selection,
                    options: Default::default(),
                    fragment_index: 0,
                },
                name: "Saved SmaI fragment".into(),
                save_primers: false,
            },
        )
        .unwrap();
        assert_eq!(installed.receipt.source, source);
        assert_eq!(
            (
                installed.receipt.duplex.top_length,
                installed.receipt.duplex.bottom_length,
                installed.receipt.duplex.bottom_start
            ),
            (8, 8, 0)
        );
        let (left, right) = installed.receipt.duplex.ends.as_ref().unwrap();
        assert!(left.polarity == collagene_core::chemistry::Polarity::Blunt);
        assert!(right.polarity == collagene_core::chemistry::Polarity::Blunt);
        assert!(left.chemistry.five_prime == collagene_core::chemistry::FivePrime::Hydroxyl);
        assert!(right.chemistry.five_prime == collagene_core::chemistry::FivePrime::Phosphate);
        assert!(installed.receipt.duplex.nicks.is_empty());
        let event = installed.receipt.product.clone();
        let capsule =
            collagene_core::material_derivation::capsule(&installed.workspace, &event).unwrap();
        let before = serde_json::to_vec(&installed.workspace).unwrap();
        let record =
            collagene_core::material_derivation::load_record(&capsule.record_bytes).unwrap();
        assert_eq!(record.request.source, source);
        assert!(matches!(
            record.request.origin,
            collagene_core::material_derivation::SourceOrigin::Saved(
                collagene_core::run_material::MaterialOrigin::Run(_)
            )
        ));
        assert!(matches!(
            record.request.derivation,
            collagene_core::material_derivation::Derivation::CatalogDigest {
                fragment_index: 0,
                ..
            }
        ));
        let material = collagene_core::run_material::MaterialOrigin::Derived(capsule);
        let pin = crate::run_origin_model::PinnedOrigin::from_material(material).unwrap();
        assert!(pin.as_run().is_err());
        let molecule = installed.workspace.revisions[&event.revision].clone();
        assert_eq!(molecule.sequence.as_str(), "ACGACCCC");
        let mut draft = original;
        draft.slots[0].captured = Some(Arc::new(CapturedSource {
            source: SourceInfo {
                binding: event,
                name: molecule.name.clone(),
                length: molecule.sequence.len(),
                alphabet: Alphabet::Dna,
                topology: molecule.topology,
            },
            molecule,
        }));
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 8;
        }
        draft.methods.rows[0].run_origin = Some(pin);
        draft.methods.rows[0].preparation = OverlapPreparation::ExactPcr;
        draft.methods.rows[0].forward = "ACGACCCC".into();
        draft.methods.rows[0].reverse_primer = "GGGGTCGT".into();
        draft.sync_rows();
        assert_eq!(draft.version, CURRENT_VERSION);
        let bytes = encode_raw(&draft).unwrap();
        let restored = decode(&bytes).unwrap();
        assert!(restored == draft);
        let mut obsolete = serde_json::to_value(&draft).unwrap();
        let current = obsolete["methods"]["rows"][0]["run_origin"]
            .as_object_mut()
            .unwrap();
        let body = current.remove("Derived").unwrap();
        current.insert("Digest".into(), body);
        assert!(
            decode(&serde_json::to_vec(&obsolete).unwrap()).is_err(),
            "Removed origin tag must never be reinterpreted"
        );

        let (_, rows) = crate::cloning_preparations::prepared_rows(&restored.methods, 1).unwrap();
        assert!(matches!(
            &rows[0].preparation,
            collagene_core::assembly_recipe::Preparation::FromDerived {
                preparation: collagene_core::assembly_recipe::RunPreparation::ExactPcr { .. },
                ..
            }
        ));
        let mut downgraded = draft.clone();
        downgraded.version = CURRENT_VERSION - 1;
        assert!(downgraded.validate().is_err());
        let mut future = draft;
        future.version = CURRENT_VERSION + 1;
        assert!(future.validate().is_err());
        assert_eq!(serde_json::to_vec(&installed.workspace).unwrap(), before);
    }

    #[test]
    fn reusable_pcr_derived_origin_keeps_native_record_and_primers() {
        let (mut draft, mut workspace) = saved_pcr_workspace();
        let source_event = {
            let document = &workspace.documents[0];
            EventBinding {
                document_id: document.id.clone(),
                event_id: document.head_event.clone().unwrap(),
                revision: document.head.clone(),
            }
        };
        workspace.active_document = source_event.document_id.clone();
        let source = &workspace.revisions[&source_event.revision];
        let forward = Primer::new("F", "ACGACCCC", "").unwrap();
        let reverse = Primer::new("R", "GTTACCCC", "").unwrap();
        let preview =
            collagene_core::verification::simulate_pcr_documents(source, &forward, &reverse)
                .unwrap()
                .remove(0);
        let request = collagene_core::material_derivation::Request {
            source: source_event,
            origin: collagene_core::material_derivation::SourceOrigin::OriginalReference,
            derivation: collagene_core::material_derivation::Derivation::ExactPcr {
                forward,
                reverse,
                selection: collagene_core::material_derivation::PcrSelection {
                    selected_index: 0,
                    expected_preview_revision: preview.molecule.revision(),
                    geometry: preview.geometry,
                },
            },
            name: "Reusable exact PCR".into(),
            save_primers: true,
        };
        let installed =
            collagene_core::material_derivation::prepare_install(&workspace, &request).unwrap();
        assert!(installed.receipt.primer_attachment.is_some());
        let product = installed.receipt.product.clone();
        let molecule = installed.workspace.revisions[&product.revision].clone();
        let capsule =
            collagene_core::material_derivation::capsule(&installed.workspace, &product).unwrap();
        draft.slots[0].captured = Some(Arc::new(CapturedSource {
            source: SourceInfo {
                binding: product,
                name: molecule.name.clone(),
                length: molecule.sequence.len(),
                alphabet: Alphabet::Dna,
                topology: molecule.topology,
            },
            molecule,
        }));
        draft.methods.rows[0].run_origin = Some(
            crate::run_origin_model::PinnedOrigin::from_material(
                run_material::MaterialOrigin::Derived(capsule),
            )
            .unwrap(),
        );
        draft.sync_rows();
        assert_eq!(draft.version, CURRENT_VERSION);
        let bytes = encode_raw(&draft).unwrap();
        assert!(decode(&bytes).unwrap() == draft);

        let (_, rows) = crate::cloning_preparations::prepared_rows(&draft.methods, 1).unwrap();
        assert!(matches!(
            &rows[0].preparation,
            collagene_core::assembly_recipe::Preparation::FromDerived { .. }
        ));
        draft.version = CURRENT_VERSION - 1;
        assert!(draft.validate().is_err());
        draft.version = CURRENT_VERSION + 1;
        assert!(draft.validate().is_err());
    }

    #[test]
    fn saved_origin_mismatch_preparation_keeps_original_event_and_policy() {
        let mut draft = saved_pcr();
        let row = &mut draft.methods.rows[0];
        row.preparation = OverlapPreparation::MismatchPcr;
        row.forward = "ACGACCCC".into();
        row.reverse_primer = "GTTACCCC".into();
        row.mismatch_max_mismatches = "2".into();
        row.mismatch_three_prime_exact = "3".into();
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        let (_, rows) = crate::cloning_preparations::prepared_rows(&restored.methods, 1).unwrap();
        let collagene_core::assembly_recipe::Preparation::FromRun {
            origin,
            preparation:
                collagene_core::assembly_recipe::RunPreparation::MismatchPcr { options, .. },
        } = &rows[0].preparation
        else {
            panic!("saved mismatch preparation")
        };
        assert_eq!(
            origin.selection.event,
            restored.slots[0].captured.as_ref().unwrap().source.binding
        );
        assert_eq!(options.binding.max_mismatches, 2);
        assert_eq!(options.binding.min_three_prime_exact, 3);
    }

    #[test]
    fn exact_saved_source_survives_v5_roundtrip_and_compiles_explicit_reconstruction() {
        let mut draft = saved_pcr();
        assert_eq!(draft.version, CURRENT_VERSION);
        let pin = draft.methods.rows[0].run_origin.clone();
        draft.sync_rows();
        assert!(draft.methods.rows[0].run_origin == pin);
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        assert!(restored == draft);
        let (_, rows) = crate::cloning_preparations::prepared_rows(&restored.methods, 1).unwrap();
        let collagene_core::assembly_recipe::Preparation::FromRun {
            origin,
            preparation,
        } = &rows[0].preparation
        else {
            panic!("saved material must remain explicit")
        };
        assert_eq!(
            origin.selection.event,
            restored.slots[0].captured.as_ref().unwrap().source.binding
        );
        assert!(matches!(
            preparation,
            collagene_core::assembly_recipe::RunPreparation::WholeLinear
        ));
        let mut old = draft.clone();
        old.version = CURRENT_VERSION - 1;
        assert!(old.validate().is_err());
        draft.methods.rows[0].preparation = OverlapPreparation::DeclaredDuplex;
        assert!(crate::cloning_preparations::prepared_rows(&draft.methods, 1).is_err());
    }

    #[test]
    fn mismatched_persisted_origin_is_refused_while_explicit_source_edit_clears_it() {
        let mut draft = saved_pcr();
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { start: value, .. } =
                &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 2;
        }
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 14;
        }
        draft.methods.rows[0].treatment = crate::cloning_methods_model::Treatment::Phosphorylate;
        draft.sync_rows();
        let mut changed = draft.clone();
        Arc::make_mut(changed.slots[0].captured.as_mut().unwrap())
            .source
            .binding
            .document_id = "unrelated-source".into();
        assert!(decode(&serde_json::to_vec(&changed).unwrap()).is_err());
        changed.sync_rows();
        assert!(changed.methods.rows[0].run_origin.is_none());
        assert_eq!(
            (
                changed.slots[0].selector.reference_interval().unwrap().0,
                changed.slots[0].selector.reference_interval().unwrap().1
            ),
            (2, 14)
        );
        assert!(changed.methods.rows[0].treatment == draft.methods.rows[0].treatment);
        assert_eq!(changed.name, draft.name);
        assert_eq!(changed.circular, draft.circular);
        assert!(changed.junction_edits == draft.junction_edits);
    }
}

/// Pin a fresh private digest selection, retaining its original source and full
/// catalog identity. This checks bindings/bounds, not cleavage; the compiled
/// CatalogDigest operation must reconstruct the actual fragment in the worker.
/// The host must check payload epoch and actual retained source-event identity.
/// Error leaves the original draft untouched.
pub fn set_digest_fragment(
    draft: &mut Draft,
    fragment: crate::digest_fragment_model::AssemblyFragment,
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
    if let Some(origin) = &fragment.origin {
        if *origin.event() != fragment.source.binding {
            return Err("Saved material belongs to a different exact source event".into());
        }
        collagene_core::material_derivation::digest_inspection_revision(
            &role,
            &fragment.molecule,
            &collagene_core::material_derivation::SourceOrigin::Saved(origin.as_ref().clone()),
            &fragment.selection,
            &fragment.options,
        )
        .map_err(|e| e.to_string())?;
    } else {
        catalog_digest::source_preflight(
            &catalog,
            &fragment.selection,
            &role,
            &fragment.molecule,
            &fragment.options,
        )
        .map_err(|e| e.to_string())?;
    }
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
        || fragment
            .origin
            .as_ref()
            .map_or_else(
                || {
                    catalog_digest::batch_revision(
                        &catalog,
                        &fragment.selection,
                        &role,
                        &fragment.molecule,
                        &fragment.options,
                    )
                },
                |origin| {
                    collagene_core::material_derivation::digest_inspection_revision(
                        &role,
                        &fragment.molecule,
                        &collagene_core::material_derivation::SourceOrigin::Saved(
                            origin.as_ref().clone(),
                        ),
                        &fragment.selection,
                        &fragment.options,
                    )
                },
            )
            .map_err(|e| e.to_string())?
            != fragment.batch_revision
    {
        return Err("Digest selection/source/options changed; run the digest again".into());
    }
    let geometry = DigestGeometry {
        material: fragment.material_domain.clone(),
        top_start: fragment.top_start,
        top_length: fragment.top_length,
        bottom_start: fragment.bottom_start,
        bottom_length: fragment.bottom_length,
    };
    if geometry.material.is_some() != fragment.origin.is_some() {
        return Err("Actual source geometry and origin must be pinned together".into());
    }
    let selector = digest_selector(
        &fragment.source,
        &geometry,
        draft.methods.rows[index].reverse,
    )?;
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
    next.clear_slot_substitutions(index)?;
    next.slots[index].captured = Some(Arc::new(CapturedSource {
        source: fragment.source,
        molecule: fragment.molecule.as_ref().clone(),
    }));
    next.slots[index].selector = selector;
    next.slots[index].digest_geometry = Some(geometry);
    next.junction_edits.clear();
    next.methods.rows[index] = OverlapInput {
        preparation: OverlapPreparation::CatalogDigest,
        run_origin: fragment.origin,
        catalog_selection: Some(fragment.selection),
        fragment: display_index.to_string(),
        repair: None,
        reverse: next.methods.rows[index].reverse,
        ..Default::default()
    };
    next.sync_rows();
    next.validate()?;
    *draft = next;
    Ok(())
}

#[cfg(test)]
mod mismatch_draft_tests {
    use super::*;
    #[test]
    fn unfinished_mismatch_primers_and_policy_are_persisted_without_claiming_valid_preparation() {
        let mut draft = Draft::default();
        draft.slots.truncate(1);
        draft.sync_rows();
        let row = &mut draft.methods.rows[0];
        row.preparation = OverlapPreparation::MismatchPcr;
        row.forward_name = "Unfinished mutation F".into();
        row.forward = String::new();
        row.reverse_primer = "GTTACCCC".into();
        row.mismatch_max_mismatches = "2e".into();
        row.mismatch_three_prime_exact = String::new();
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        assert!(restored == draft);
        assert!(crate::cloning_preparations::prepared_rows(&restored.methods, 1).is_err());
        draft.methods.rows[0].forward = "ACGACCCC".into();
        // Even with complete oligos, unfinished policy text remains an editable
        // draft while fresh preparation compilation reports the local error.
        assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
        let error = crate::cloning_preparations::prepared_rows(&draft.methods, 1)
            .err()
            .unwrap();
        assert!(error.contains("allowed mismatches"));
    }
    #[test]
    fn mismatch_policy_roundtrip_preserves_explicit_primers_and_unfinished_values() {
        let mut draft = Draft::default();
        let row = &mut draft.methods.rows[0];
        row.preparation = OverlapPreparation::MismatchPcr;
        row.forward = "ACGACCCC".into();
        row.reverse_primer = "GTTACCCC".into();
        row.mismatch_max_mismatches = "2".into();
        row.mismatch_three_prime_exact = "3".into();
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        assert!(restored == draft);
        draft.methods.rows[0].mismatch_three_prime_exact = "9".into();
        assert!(draft.validate().is_ok());
        assert!(draft.methods.rows[0].mismatch_options().is_err());
        assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
        draft.methods.rows[0].preparation = OverlapPreparation::ExactPcr;
        draft.methods.rows[0].mismatch_max_mismatches = "x".repeat(21);
        assert!(draft.validate().is_err());
    }
}

#[cfg(test)]
mod substitution_tests {
    use super::*;
    use collagene_core::assembly_substitution_intent::Substitution;
    #[test]
    fn literal_substitutions_require_stable_roles_and_survive_navigation_and_roundtrip() {
        let mut draft = Draft::default();
        let roles = draft.methods.role_ids(draft.slots.len()).unwrap();
        draft.methods.source_roles = Some(roles.clone());
        draft.substitutions = vec![Substitution {
            role: roles[0].role.clone(),
            offset: 2,
            expected: "G".into(),
            replacement: "T".into(),
        }];
        draft.sync_rows();
        assert_eq!(draft.version, CURRENT_VERSION);
        draft.active_slot = 1;
        draft.sync_rows();
        assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
        assert!(crate::assembly_alternatives::enable(&draft, "patch-library").is_err());
        let mut old = draft.clone();
        old.version = CURRENT_VERSION - 1;
        assert!(old.validate().is_err());
        old = draft.clone();
        old.methods.source_roles = None;
        assert!(old.validate().is_err());
        let original = draft.substitutions.clone();
        draft.slots.swap(0, 1);
        draft.methods.swap_rows(0, 1);
        draft.sync_rows();
        assert!(draft.substitutions == original);
        assert_eq!(draft.methods.role_ids(2).unwrap()[1].role, original[0].role);
        draft.clear_slot_substitutions(0).unwrap();
        assert!(draft.substitutions == original);
        draft.clear_slot_substitutions(1).unwrap();
        assert!(draft.substitutions.is_empty());
    }
    #[test]
    fn bounded_substitution_deserialization_preserves_documentary_shape() {
        let mut draft = Draft::default();
        let roles = draft.methods.role_ids(2).unwrap();
        draft.methods.source_roles = Some(roles.clone());
        let edit = Substitution {
            role: roles[0].role.clone(),
            offset: 0,
            expected: "A".into(),
            replacement: "T".into(),
        };
        draft.substitutions = vec![edit.clone()];
        draft.sync_rows();
        draft.substitutions.push(edit.clone());
        assert!(draft.validate().is_err());
        let mut oversized = serde_json::to_value(Draft::default()).unwrap();
        oversized["version"] = CURRENT_VERSION.into();
        oversized["substitutions"] = serde_json::to_value(vec![edit; 257]).unwrap();
        assert!(decode(&serde_json::to_vec(&oversized).unwrap()).is_err());
    }
}

#[cfg(test)]
mod mismatch_search_limit_tests {
    use super::*;
    #[test]
    fn search_limits_and_unfinished_numeric_text_roundtrip_without_changing_intent() {
        let mut draft = Draft::default();
        let row = &mut draft.methods.rows[1];
        row.preparation = OverlapPreparation::MismatchPcrCatalogDigest;
        row.forward = "AGACCGAT".into();
        row.reverse_primer = "TCCAGGTC".into();
        row.mismatch_max_products = "1".into();
        row.mismatch_max_work = "1e".into();
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        assert!(restored == draft);
        assert!(restored.methods.rows[1].mismatch_options().is_err());
        draft.methods.rows[1].mismatch_max_work = "1000".into();
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        let options = restored.methods.rows[1].mismatch_options().unwrap();
        assert_eq!(options.max_products, 1);
        assert_eq!(options.binding.max_work, 1000);
        assert!(
            restored.slots == draft.slots
                && restored.junction_edits == draft.junction_edits
                && restored.substitutions == draft.substitutions
        );
        draft.methods.rows[1].preparation = OverlapPreparation::WholeLinear;
        draft.methods.rows[1].mismatch_max_pairs = "x".repeat(21);
        assert!(draft.validate().is_err());
    }
}

#[cfg(test)]
mod contribution_edit_tests {
    use super::*;
    use collagene_core::{assembly_edit_intent::Edit, assembly_substitution_intent::Substitution};
    fn edited() -> Draft {
        let mut d = Draft::default();
        let roles = d.methods.role_ids(d.slots.len()).unwrap();
        d.edits.push(Edit {
            role: roles[0].role.clone(),
            start: 2,
            end: 4,
            expected: "CC".into(),
            replacement: "N".into(),
        });
        d.methods.source_roles = Some(roles);
        d.sync_rows();
        d
    }
    #[test]
    fn persisted_general_edits_require_stable_roles_and_bounded_literals() {
        let d = edited();
        assert_eq!(d.version, CURRENT_VERSION);
        assert!(d.validate().is_ok());
        assert!(decode(&encode_raw(&d).unwrap()).unwrap() == d);
        let mut wrong = d.clone();
        wrong.version = CURRENT_VERSION - 1;
        assert!(wrong.validate().is_err());
        wrong = d.clone();
        wrong.methods.source_roles = None;
        assert!(wrong.validate().is_err());
        wrong = d.clone();
        wrong.edits[0].expected = "C".into();
        assert!(wrong.validate().is_err());
        wrong = d.clone();
        wrong.edits[0].replacement = "Z".into();
        assert!(wrong.validate().is_err());
        wrong = d.clone();
        wrong.edits[0].replacement = "A".repeat(100_001);
        assert!(wrong.validate().is_err());
        let mut oversized = serde_json::to_value(&d).unwrap();
        oversized["edits"] = serde_json::to_value(vec![d.edits[0].clone(); 257]).unwrap();
        assert!(decode(&serde_json::to_vec(&oversized).unwrap()).is_err());
    }
    #[test]
    fn general_edits_reject_ambiguous_order_and_cross_kind_overlap() {
        let d = edited();
        let role = d.edits[0].role.clone();
        let mut next = d.clone();
        next.edits.push(Edit {
            role: role.clone(),
            start: 4,
            end: 5,
            expected: "G".into(),
            replacement: "T".into(),
        });
        assert!(next.validate().is_ok()); // Adjacent nonempty intervals have an unambiguous order.
        next.edits[1] = Edit {
            role: role.clone(),
            start: 4,
            end: 4,
            expected: "".into(),
            replacement: "T".into(),
        };
        assert!(next.validate().is_err());
        next = d.clone();
        next.substitutions.push(Substitution {
            role: role.clone(),
            offset: 3,
            expected: "C".into(),
            replacement: "G".into(),
        });
        assert!(next.validate().is_err());
        next = d.clone();
        next.substitutions.push(Substitution {
            role: role.clone(),
            offset: 4,
            expected: "G".into(),
            replacement: "T".into(),
        });
        assert!(next.validate().is_ok());
        next.edits[0] = Edit {
            role,
            start: 4,
            end: 4,
            expected: "".into(),
            replacement: "C".into(),
        };
        assert!(next.validate().is_err());
    }
}

#[cfg(test)]
mod extension_draft_tests {
    use super::*;
    use crate::cloning_methods_model::OverlapExtensionPcrInput;

    #[test]
    fn six_independent_unfinished_oligos_and_policy_roundtrip() {
        let mut draft = Draft::default();
        let mut raw = OverlapExtensionPcrInput::default();
        raw.left.forward.annealing = "ag tc unfinished".into();
        raw.left.reverse.tail = "GGTCTGAAAA".into();
        raw.right.forward.tail = "GCTACGTATT".into();
        raw.right.reverse.name = "Right outer stock primer".into();
        raw.outer.forward.tail = " separate final assembly tail ".into();
        raw.outer.reverse.annealing = "TGACGATC".into();
        raw.overlap = "1e".into();
        raw.max_candidates = String::new();
        draft.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcr;
        draft.methods.rows[1].overlap_extension_pcr = Some(raw);
        draft.sync_rows();
        assert_eq!(draft.version, CURRENT_VERSION);
        let bytes = encode_raw(&draft).unwrap();
        let restored = decode(&bytes).unwrap();
        assert!(restored == draft);
        assert!(crate::cloning_preparations::prepared_rows(&restored.methods, 2).is_err());

        let mut old = restored.clone();
        old.version = CURRENT_VERSION - 1;
        assert!(old.validate().is_err());
        // Changing preparation preserves the inactive editable fields.
        old.methods.rows[1].preparation = OverlapPreparation::WholeLinear;
        assert!(old.validate().is_err());
        old.sync_rows();
        assert_eq!(old.version, CURRENT_VERSION);
        assert!(decode(&encode_raw(&old).unwrap()).unwrap() == old);
        let mut oversized = old;
        oversized.methods.rows[1]
            .overlap_extension_pcr
            .as_mut()
            .unwrap()
            .outer
            .reverse
            .tail = "A".repeat(1001);
        assert!(encode_raw(&oversized).is_err());
    }

    #[test]
    fn physical_homology_settings_roundtrip_unfinished_inactive_values() {
        let mut draft = Draft::default();
        assert!(draft.methods.homology.is_none());
        let mut form = crate::cloning_methods_model::HomologyForm::default();
        form.enabled = false;
        form.exact_homology = false;
        form.cuts = "20:unfinished, 3:".into();
        form.outer_left = "partial".into();
        draft.methods.homology = Some(form);
        draft.sync_rows();
        assert_eq!(draft.version, CURRENT_VERSION);
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        assert!(restored == draft);
        let mut downgraded = restored.clone();
        downgraded.version = CURRENT_VERSION - 1;
        assert!(downgraded.validate().is_err());

        draft.methods.homology.as_mut().unwrap().cuts = "x".repeat(1025);
        assert!(encode_raw(&draft).is_err());
    }
    #[test]
    fn three_prime_settings_preserve_independent_five_prime_settings() {
        let mut draft = Draft::default();
        assert!(draft.methods.annealing.is_none());
        draft.methods.homology = Some(crate::cloning_methods_model::HomologyForm::default());
        draft.sync_rows();
        assert_eq!(draft.version, CURRENT_VERSION);
        let current = decode(&encode_raw(&draft).unwrap()).unwrap();
        assert!(current.methods.annealing.is_none());
        assert!(current.methods.homology.as_ref().unwrap().enabled);
        let mut form = crate::cloning_methods_model::AnnealingForm::default();
        form.enabled = false;
        form.cuts = "0:unfinished, 20:".into();
        form.max_candidates = "partial".into();
        draft.methods.annealing = Some(form);
        draft.sync_rows();
        assert_eq!(draft.version, CURRENT_VERSION);
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        assert!(restored == draft);
        assert!(restored.methods.homology.as_ref().unwrap().enabled);
        assert!(!restored.methods.annealing.as_ref().unwrap().enabled);
        let mut downgraded = restored.clone();
        downgraded.version = CURRENT_VERSION - 1;
        assert!(downgraded.validate().is_err());

        draft.methods.annealing.as_mut().unwrap().cuts = "x".repeat(1025);
        assert!(encode_raw(&draft).is_err());
    }
    #[test]
    fn extension_catalog_selector_preserves_unfinished_raw6() {
        let mut draft = Draft::default();
        draft.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcrCatalogDigest;
        draft.sync_rows();
        assert_eq!(draft.version, CURRENT_VERSION);
        assert!(draft.methods.rows[1].overlap_extension_pcr.is_none());
        let mut raw = OverlapExtensionPcrInput::default();
        raw.outer.forward.tail = "CCGAATTC".into();
        raw.right.forward.annealing = "unfinished".into();
        draft.methods.rows[1].overlap_extension_pcr = Some(raw);
        draft.sync_rows();
        let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
        assert!(restored == draft);

        let mut downgraded = draft.clone();
        downgraded.version = CURRENT_VERSION - 1;
        assert!(downgraded.validate().is_err());
        downgraded.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcr;
        downgraded.version = CURRENT_VERSION;
        assert!(downgraded.validate().is_ok());
        assert!(decode(&encode_raw(&downgraded).unwrap()).unwrap() == downgraded);
    }
    #[test]
    fn selected_empty_extension_preparation_roundtrips_explicit_null() {
        let mut draft = Draft::default();
        draft.methods.rows[0].preparation = OverlapPreparation::OverlapExtensionPcr;
        assert!(draft.validate().is_ok());
        draft.sync_rows();
        assert_eq!(draft.version, CURRENT_VERSION);
        assert!(draft.methods.rows[0].overlap_extension_pcr.is_none());
        assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
    }
}

// Shared original-source declarations for review and preparation proposals.
// Both callers use the same borrowed preflights before copying source records.
pub(crate) fn inputs_request(
    draft: &Draft,
    workspace: Option<&collagene_core::workspace::Workspace>,
) -> Result<
    (
        Vec<PlanInput>,
        assembly_intent::Request,
        Vec<Arc<CapturedSource>>,
    ),
    String,
> {
    draft.validate()?;
    crate::assembly_alternatives::require_active_choice_match(draft)?;
    let mut captures = draft
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
    let roles = draft.methods.role_ids(captures.len())?;
    let ids = draft.methods.all_input_ids()?;
    for (index, row) in draft.methods.rows.iter().enumerate() {
        if row.preparation == OverlapPreparation::AnnealedOligos {
            captures.push(draft.slots[index].anneal_bottom.clone().ok_or_else(|| {
                format!(
                    "Pin the bottom oligo for {} before reviewing intended DNA.",
                    draft.slots[index].label
                )
            })?);
        }
    }
    let mut material_bindings = Vec::new();
    let fragments = draft
        .slots
        .iter()
        .enumerate()
        .map(|(i, slot)| {
            if matches!(
                &slot.selector,
                assembly_intent::Selector::MaterialStrandRange { .. }
            ) {
                let origin = draft.methods.rows[i]
                    .run_origin
                    .as_ref()
                    .ok_or("Actual strand needs its saved source origin")?;
                material_bindings.push(assembly_intent::MaterialBinding {
                    input_id: ids[i].clone(),
                    event: captures[i].source.binding.clone(),
                    origin: origin.as_ref().clone(),
                });
            }
            Ok(assembly_intent::Contribution {
                role: roles[i].role.clone(),
                input_id: ids[i].clone(),
                selector: slot.selector.clone(),
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let request = assembly_intent::Request {
        name: draft.name.clone(),
        circular: draft.circular,
        fragments,
        material_bindings,
    };
    let borrowed = ids
        .iter()
        .zip(&captures)
        .map(|(id, c)| (id.as_str(), &c.molecule))
        .collect::<Vec<_>>();
    assembly_intent::preflight(&borrowed, &request, workspace).map_err(|e| e.to_string())?;
    if !draft.junction_edits.is_empty() {
        assembly_junction_edit::intent_preflight(
            &borrowed,
            &request,
            &draft.junction_edits,
            workspace,
        )
        .map_err(|e| e.to_string())?;
    }
    if !draft.substitutions.is_empty() {
        assembly_substitution_intent::source_preflight(
            &borrowed,
            &assembly_substitution_intent::Request {
                base: request.clone(),
                junction_edits: draft.junction_edits.clone(),
                substitutions: draft.substitutions.clone(),
            },
            workspace,
        )
        .map_err(|e| e.to_string())?;
    }
    if !draft.edits.is_empty() {
        assembly_edit_intent::source_preflight(
            &borrowed,
            &edit_request(draft, request.clone()),
            workspace,
        )
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

pub(crate) fn edit_request(
    draft: &Draft,
    base: assembly_intent::Request,
) -> assembly_edit_intent::Request {
    assembly_edit_intent::Request {
        base: assembly_substitution_intent::Request {
            base,
            junction_edits: draft.junction_edits.clone(),
            substitutions: draft.substitutions.clone(),
        },
        edits: draft.edits.clone(),
    }
}

#[cfg(test)]
mod current_codec_tests {
    use super::*;
    #[test]
    fn current_schema_requires_all_fields_including_explicit_nulls_and_defaults() {
        let draft = Draft::default();
        let value = serde_json::to_value(&draft).unwrap();
        assert!(decode(&serde_json::to_vec(&value).unwrap()).unwrap() == draft);
        for pointer in [
            "/junction_edits",
            "/substitutions",
            "/edits",
            "/alternatives",
            "/slots/0/selector",
            "/slots/0/captured",
            "/slots/0/digest_geometry",
            "/slots/0/anneal_bottom",
            "/methods/tailing_profile",
            "/methods/vector/action",
            "/methods/golden_gate",
            "/methods/golden_gate_pool",
            "/methods/multisite",
            "/methods/source_roles",
            "/methods/homology",
            "/methods/annealing",
            "/methods/rows/0/run_origin",
            "/methods/rows/0/source_generation",
            "/methods/rows/0/overlap_extension_pcr",
            "/methods/rows/0/catalog_selection",
            "/methods/rows/0/mismatch_max_mismatches",
            "/methods/rows/0/anneal_options",
            "/methods/rows/0/anneal_candidate",
            "/methods/rows/0/anneal_bottom_document",
            "/methods/rows/0/anneal_source_context",
        ] {
            let mut missing = value.clone();
            let (parent, field) = pointer.rsplit_once('/').unwrap();
            missing
                .pointer_mut(parent)
                .unwrap()
                .as_object_mut()
                .unwrap()
                .remove(field);
            assert!(
                decode(&serde_json::to_vec(&missing).unwrap()).is_err(),
                "Missing field accepted: {pointer}"
            );
        }
        for version in [CURRENT_VERSION - 1, CURRENT_VERSION + 1] {
            let mut unsupported = value.clone();
            unsupported["version"] = version.into();
            assert!(decode(&serde_json::to_vec(&unsupported).unwrap()).is_err());
        }
        let mut unknown = value;
        unknown["methods"]["obsolete_control"] = true.into();
        assert!(decode(&serde_json::to_vec(&unknown).unwrap()).is_err());
    }
}

pub(crate) fn selector_label(selector: &assembly_intent::Selector) -> String {
    match selector {
        assembly_intent::Selector::ReferenceRange {
            start,
            end,
            reverse,
        } => format!(
            "Reference {}–{} · {}",
            start.saturating_add(1),
            end,
            if *reverse { "reverse" } else { "forward" }
        ),
        assembly_intent::Selector::MaterialStrandRange {
            strand, start, end, ..
        } => format!(
            "Actual {} strand {}–{} · 5′→3′",
            match strand {
                collagene_core::chemistry::Strand::Top => "top",
                collagene_core::chemistry::Strand::Bottom => "bottom",
            },
            start + 1,
            end
        ),
    }
}
pub(crate) fn orient_fragment(draft: &Draft, index: usize, reverse: bool) -> Result<Draft, String> {
    let mut next = draft.clone();
    let slot = next
        .slots
        .get(index)
        .ok_or("Choose a source contribution first")?;
    let selector = if let Some((capture, geometry)) =
        slot.captured.as_ref().zip(slot.digest_geometry.as_ref())
    {
        digest_selector(&capture.source, geometry, reverse)?
    } else {
        match &slot.selector {
            assembly_intent::Selector::ReferenceRange { start, end, .. } => {
                assembly_intent::Selector::ReferenceRange {
                    start: *start,
                    end: *end,
                    reverse,
                }
            }
            assembly_intent::Selector::MaterialStrandRange { .. } => {
                return Err("Actual strand orientation needs its pinned fragment geometry".into());
            }
        }
    };
    next.clear_slot_substitutions(index)?;
    next.slots[index].selector = selector;
    next.methods.rows[index].reverse = reverse;
    next.junction_edits.clear();
    next.sync_rows();
    next.validate()?;
    Ok(next)
}
pub(crate) fn require_reference_selectors(draft: &Draft) -> Result<(), String> {
    if draft.slots.iter().any(|slot| {
        matches!(
            &slot.selector,
            assembly_intent::Selector::MaterialStrandRange { .. }
        )
    }) {
        return Err("This primer/library workflow requires original reference ranges. Actual-strand contributions remain available for direct assembly.".into());
    }
    Ok(())
}
