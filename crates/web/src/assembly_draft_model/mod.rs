//! Editable assembly data and versioned persistence, without GUI state or receipts.
use crate::{
    cloning_methods_model::{MethodsState, Mode, OverlapInput, OverlapPreparation},
    source_identity::SourceInfo,
};
use collagene_core::{
    assembly_edit_intent, assembly_intent, assembly_junction_edit, assembly_substitution_intent,
    verification::PlanInput,
};
use collagene_model::{Alphabet, Molecule, Topology};
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
    pub material_origin: Option<crate::run_origin_model::PinnedOrigin>,
    #[serde(deserialize_with = "required_option")]
    pub digest_geometry: Option<DigestGeometry>,
    #[serde(with = "optional_snapshot")]
    pub anneal_bottom: Option<Arc<CapturedSource>>,
}
impl Slot {
    pub(crate) fn empty(label: &str) -> Self {
        Self {
            label: label.into(),
            material_origin: None,
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
    if row.run_origin != slot.material_origin {
        return Err(
            "Physical preparation must retain the selected contribution material origin.".into(),
        );
    }
    if let Some(origin) = &slot.material_origin {
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
    let selected_reverse = match &slot.selector {
        assembly_intent::Selector::ReferenceRange { reverse, .. } => *reverse,
        assembly_intent::Selector::MaterialStrandRange { strand, .. } => {
            *strand == collagene_core::chemistry::Strand::Bottom
        }
    };
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
            if slot.material_origin.is_none() || !hex(expected_state_revision) {
                return Err("Actual strand intent needs its exact saved material origin".into());
            }
            let geometry = slot
                .digest_geometry
                .as_ref()
                .ok_or("Actual strand intent needs its captured source domain")?;
            if digest_selector(source, geometry, selected_reverse)? != slot.selector {
                return Err(
                    "Actual strand selection differs from its pinned digest/source state".into(),
                );
            }
        }
    }
    if let Some(geometry) = &slot.digest_geometry {
        if digest_selector(source, geometry, selected_reverse)? != slot.selector {
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
pub(crate) fn validate_internal_declarations(
    methods: &MethodsState,
    slots: usize,
    substitutions: &[collagene_core::assembly_substitution_intent::Substitution],
    edits: &[collagene_core::assembly_edit_intent::Edit],
) -> Result<(), String> {
    if substitutions.len() > 256 {
        return Err("At most256 explicit substitutions are supported.".into());
    }
    if !substitutions.is_empty() {
        if methods.source_roles.is_none() {
            return Err("Pin stable source roles before declaring internal substitutions.".into());
        }
        let roles = methods.role_ids(slots)?;
        let mut seen = std::collections::BTreeSet::new();
        for edit in substitutions {
            if !roles.iter().any(|role| role.role == edit.role)
                || !token(&edit.role, 128)
                || edit.offset >= collagene_model::MAX_SEQUENCE_LEN
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

    if edits.len() > 256 {
        return Err("At most256 explicit contribution edits are supported.".into());
    }
    if edits.is_empty() {
        return Ok(());
    }
    if methods.source_roles.is_none() {
        return Err("Internal contribution edits require stable source roles.".into());
    }
    let roles = methods.role_ids(slots)?;
    let mut total = 0usize;
    for (index, edit) in edits.iter().enumerate() {
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
            || edit.end > collagene_model::MAX_SEQUENCE_LEN
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
        if edits[..index].iter().any(|other| {
            other.role == edit.role && conflict(edit.start, edit.end, other.start, other.end)
        }) || substitutions.iter().any(|other| {
            other.role == edit.role
                && conflict(edit.start, edit.end, other.offset, other.offset + 1)
        }) {
            return Err("Contribution edits must not overlap each other or substitutions, or share insertion boundaries.".into());
        }
    }
    Ok(())
}
#[derive(Clone, PartialEq, Eq)]
pub(crate) struct IntentIdentity {
    name: String,
    circular: bool,
    roles: Option<Vec<crate::assembly_role_ids::Role>>,
    combination: Option<(String, String, Vec<String>)>,
    junction_edits: Vec<assembly_junction_edit::Edit>,
    substitutions: Vec<assembly_substitution_intent::Substitution>,
    edits: Vec<assembly_edit_intent::Edit>,
    material_origins: Vec<Option<crate::run_origin_model::PinnedOrigin>>,
    paired_rows: Vec<usize>,
    secondary_sources: Vec<Option<SourceInfo>>,
    slots: Vec<(
        String,
        Option<SourceInfo>,
        collagene_core::assembly_intent::Selector,
        Option<DigestGeometry>,
    )>,
}
impl IntentIdentity {
    pub(crate) fn slot_label(&self, index: usize) -> Option<&str> {
        self.slots.get(index).map(|slot| slot.0.as_str())
    }

    pub(crate) fn has_material_sources(&self) -> bool {
        self.slots.iter().any(|slot| {
            matches!(
                &slot.2,
                assembly_intent::Selector::MaterialStrandRange { .. }
            )
        })
    }

    pub(crate) fn has_sequence_edits(&self) -> bool {
        !self.junction_edits.is_empty() || !self.substitutions.is_empty() || !self.edits.is_empty()
    }
}
impl Draft {
    pub(crate) fn identity(&self) -> IntentIdentity {
        IntentIdentity {
            name: self.name.clone(),
            circular: self.circular,
            roles: self.methods.source_roles.clone().or_else(|| {
                (self.methods.mode == Mode::MultiSiteReference)
                    .then(|| self.methods.role_ids(self.slots.len()).ok())
                    .flatten()
            }),
            combination: crate::assembly_alternatives::identity(self),
            junction_edits: self.junction_edits.clone(),
            substitutions: self.substitutions.clone(),
            edits: self.edits.clone(),
            material_origins: self
                .slots
                .iter()
                .map(|slot| slot.material_origin.clone())
                .collect(),
            paired_rows: self
                .methods
                .rows
                .iter()
                .enumerate()
                .filter_map(|(i, row)| {
                    (row.preparation == OverlapPreparation::AnnealedOligos).then_some(i)
                })
                .collect(),
            secondary_sources: self
                .slots
                .iter()
                .map(|slot| slot.anneal_bottom.as_ref().map(|s| s.source.clone()))
                .collect(),
            slots: self
                .slots
                .iter()
                .map(|s| {
                    (
                        s.label.clone(),
                        s.captured.as_ref().map(|c| c.source.clone()),
                        s.selector.clone(),
                        s.digest_geometry.clone(),
                    )
                })
                .collect(),
        }
    }

    pub(crate) fn same_declaration(&self, applied: &Draft) -> bool {
        // The active slot is view state, not a change to the declared design.
        self.version == applied.version
            && self.name == applied.name
            && self.circular == applied.circular
            && self.slots == applied.slots
            && self.methods == applied.methods
            && self.junction_edits == applied.junction_edits
            && self.substitutions == applied.substitutions
            && self.edits == applied.edits
            && self.alternatives == applied.alternatives
    }
    pub(crate) fn clear_slot_substitutions(&mut self, index: usize) -> Result<(), String> {
        let roles = self.methods.role_ids(self.slots.len())?;
        let role = roles.get(index).ok_or("The source slot is unavailable.")?;
        self.substitutions.retain(|edit| edit.role != role.role);
        self.edits.retain(|edit| edit.role != role.role);
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
        validate_internal_declarations(
            &self.methods,
            self.slots.len(),
            &self.substitutions,
            &self.edits,
        )?;

        if self.alternatives.is_some() && self.uses_paired_fields() {
            return Err(
                "Paired-oligo sources/settings are not supported in assembly libraries.".into(),
            );
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
        collagene_model::bounds::model_preflight(&models).map_err(|e| e.to_string())?;
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
        for (i, (slot, row)) in self
            .slots
            .iter_mut()
            .zip(&mut self.methods.rows)
            .enumerate()
        {
            if slot.material_origin.as_ref().is_some_and(|origin| {
                slot.captured
                    .as_ref()
                    .is_none_or(|capture| *origin.event() != capture.source.binding)
            }) {
                slot.material_origin = None;
            }
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
        if geometry.material.is_some() {
            false
        } else {
            draft.methods.rows[index].reverse
        },
    )?;
    let mut borrowed = draft
        .slots
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != index)
        .filter_map(|(_, slot)| slot.captured.as_ref().map(|c| &c.molecule))
        .collect::<Vec<_>>();
    borrowed.push(&fragment.molecule);
    collagene_model::bounds::model_preflight(&borrowed).map_err(|e| e.to_string())?;
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
    next.slots[index].material_origin = fragment.origin.clone();
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
    collagene_core::run_material::material_capsules_preflight(
        draft
            .slots
            .iter()
            .filter_map(|slot| slot.material_origin.as_ref())
            .map(|origin| match origin.as_ref() {
                collagene_core::run_material::MaterialOrigin::Run(c) => {
                    collagene_core::run_material::MaterialOriginRef::Run(c)
                }
                collagene_core::run_material::MaterialOrigin::Derived(c) => {
                    collagene_core::run_material::MaterialOriginRef::Derived(c)
                }
            }),
    )
    .map_err(|e| e.to_string())?;
    let mut material_bindings = Vec::new();
    let fragments = draft
        .slots
        .iter()
        .enumerate()
        .map(|(i, slot)| {
            if let Some(origin) = &slot.material_origin {
                material_bindings.push(assembly_intent::MaterialBinding {
                    input_id: ids[i].clone(),
                    event: captures[i].source.binding.clone(),
                    origin: origin.as_ref().clone(),
                });
            } else if matches!(
                &slot.selector,
                assembly_intent::Selector::MaterialStrandRange { .. }
            ) {
                return Err("Actual strand needs its selected source origin".into());
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
    if matches!(
        &next.slots[index].selector,
        assembly_intent::Selector::ReferenceRange { .. }
    ) {
        next.methods.rows[index].reverse = reverse;
    }
    next.junction_edits.clear();
    next.sync_rows();
    next.validate()?;
    Ok(next)
}
pub(crate) fn set_material_origin(
    draft: &Draft,
    index: usize,
    pin: Option<crate::run_origin_model::PinnedOrigin>,
) -> Result<Draft, String> {
    let slot = draft
        .slots
        .get(index)
        .ok_or("The source contribution is unavailable")?;
    let capture = slot
        .captured
        .as_ref()
        .ok_or("Pin the source contribution first")?;
    if let Some(pin) = &pin {
        pin.validate_documentary()?;
        if pin.event() != &capture.source.binding {
            return Err("Saved material belongs to another exact source event".into());
        }
    }
    if matches!(
        &slot.selector,
        assembly_intent::Selector::MaterialStrandRange { .. }
    ) && pin != slot.material_origin
    {
        return Err(
            "Select a fresh digest fragment before changing an actual-strand origin".into(),
        );
    }
    let mut next = draft.clone();
    next.slots[index].material_origin = pin.clone();
    next.methods.rows[index].run_origin = pin;
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
        return Err("This single-contribution primer workflow requires original reference ranges. Use the whole-construct primer proposal for actual-strand contributions.".into());
    }
    Ok(())
}

#[cfg(test)]
#[path = "tests/contribution_edit_tests.rs"]
mod contribution_edit_tests;
#[cfg(test)]
#[path = "tests/current_codec_tests.rs"]
mod current_codec_tests;
#[cfg(test)]
#[path = "tests/extension_draft_tests.rs"]
mod extension_draft_tests;
#[cfg(test)]
#[path = "tests/mismatch_draft_tests.rs"]
mod mismatch_draft_tests;
#[cfg(test)]
#[path = "tests/mismatch_search_limit_tests.rs"]
mod mismatch_search_limit_tests;
#[cfg(test)]
#[path = "tests/paired_oligo_tests.rs"]
mod paired_oligo_tests;
#[cfg(test)]
#[path = "tests/run_origin_tests.rs"]
mod run_origin_tests;
#[cfg(test)]
#[path = "tests/substitution_tests.rs"]
mod substitution_tests;
#[cfg(test)]
#[path = "tests/tailing_tests.rs"]
mod tailing_tests;
