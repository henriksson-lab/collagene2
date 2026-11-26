//! Editable assembly data and versioned persistence, without GUI state or receipts.
use crate::{
    cloning_methods_model::{MethodsState, Mode, OverlapInput, OverlapPreparation},
    source_identity::SourceInfo,
};
use collagene_core::{Alphabet, Molecule, Topology, assembly_junction_edit};
use serde::{Deserialize, Serialize};
use std::{io::Write, sync::Arc};

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
    pub(crate) fn empty(label: &str) -> Self {
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
    #[serde(default)]
    pub alternatives: Option<crate::assembly_alternatives::Library>,
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
            version: 3,
            name: "New assembly".into(),
            circular: true,
            slots: vec![Slot::empty("Backbone"), Slot::empty("Insert")],
            active_slot: 0,
            methods: MethodsState::default(),
            junction_edits: Vec::new(),
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
        // Source preparation stays in rows. Legacy tail controls contribute only
        // the explicitly chosen base/treatment, never another source selection.
        for (index, tail) in [&methods.vector, &methods.insert].into_iter().enumerate() {
            if tail.document != index || tail.reversed != methods.rows[index].reverse {
                return Err("TA/GC tail controls differ from the pinned source slots.".into());
            }
        }
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
                Mode::Overlap | Mode::GoldenGate | Mode::MultiSiteReference | Mode::Tailing
            ))
    {
        return Err(
            "Junction edits require declared reference intent, overlap, Golden Gate or TA/GC assembly and at most12 edits".into(),
        );
    }
    let mut junctions = std::collections::BTreeSet::new();
    for edit in junction_edits {
        let terminal = !circular && methods.mode == Mode::Tailing && edit.junction == count - 1;
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
    for row in [&methods.vector, &methods.insert] {
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
        &methods.overlaps,
        &methods.comparisons,
        &methods.occurrences,
    ] {
        if text.len() > 65536 {
            return Err("Assembly option field exceeds64KiB".into());
        }
    }
    methods.role_ids(count)?;
    Ok(())
}
pub(crate) fn validate_slot_link(
    slot: &Slot,
    row: &OverlapInput,
    index: usize,
) -> Result<(), String> {
    if row.document != index
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
    Ok(())
}
pub(crate) fn validate_slot(slot: &Slot, row: &OverlapInput, index: usize) -> Result<(), String> {
    if !token(&slot.label, 80) {
        return Err("Slot labels must be nonblank and at most80 bytes".into());
    }
    validate_preparation_input(row)?;
    validate_slot_link(slot, row, index)?;
    let Some(captured) = slot.captured.as_deref() else {
        if slot.start != 0 || slot.end != 0 || slot.digest_geometry.is_some() {
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
    Ok(())
}
impl Draft {
    /// Structural/bounded source consistency only. Historical event IDs require
    /// the host workspace to attest them before installing any candidate.
    pub fn validate(&self) -> Result<(), String> {
        if !matches!(self.version, 1 | 2 | 3) || self.version == 1 && self.alternatives.is_some() {
            return Err(
                "Unsupported assembly draft version or alternatives in a legacy draft.".into(),
            );
        }
        if self.methods.mode == Mode::Tailing && self.version < 3 {
            return Err("Ordinary TA/GC cloning requires assembly draft version 3.".into());
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
        if self.version < 3
            && (self
                .methods
                .golden_gate_pool
                .as_ref()
                .is_some_and(|pool| pool.candidate.is_some())
                || self.alternatives.as_ref().is_some_and(|library| {
                    library.rows.iter().any(|row| {
                        row.edits
                            .methods
                            .golden_gate_pool
                            .as_ref()
                            .is_some_and(|pool| pool.candidate.is_some())
                    })
                }))
        {
            return Err("Pool candidate declarations require assembly draft version 3.".into());
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
            .filter_map(|slot| slot.captured.as_deref().map(|c| &c.molecule))
            .collect::<Vec<_>>();
        collagene_core::model::bounds::model_preflight(&models).map_err(|e| e.to_string())?;
        let mut bases = 0usize;
        for (index, (slot, row)) in self.slots.iter().zip(&self.methods.rows).enumerate() {
            if let Some(captured) = &slot.captured {
                bases = bases
                    .checked_add(captured.molecule.sequence.len())
                    .filter(|n| *n <= 8 << 20)
                    .ok_or("Assembly pinned source bases exceed8MiB")?;
            }
            validate_slot(slot, row, index)?;
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
            row.document = i;
            row.reverse = slot.reverse;
            row.interval_start = slot.start.saturating_add(1).to_string();
            row.interval_end = if slot.captured.is_some() {
                slot.end.to_string()
            } else {
                String::new()
            };
        }
        if self.methods.mode == Mode::Tailing && self.slots.len() == 2 {
            self.methods.vector.document = 0;
            self.methods.vector.reversed = self.slots[0].reverse;
            self.methods.insert.document = 1;
            self.methods.insert.reversed = self.slots[1].reverse;
        }
    }
}

#[cfg(test)]
pub fn encode(draft: &Draft) -> Result<Vec<u8>, String> {
    let stored = crate::assembly_alternatives::for_storage(draft)?;
    encode_raw(&stored)
}
/// Encode current editable state without promoting working fields into saved rows.
pub fn encode_raw(draft: &Draft) -> Result<Vec<u8>, String> {
    if draft.version != 3 {
        return Err("Raw assembly persistence requires draft version 3.".into());
    }
    draft.validate()?;
    let mut writer = Bounded { bytes: Vec::new() };
    serde_json::to_writer(&mut writer, draft).map_err(|e| e.to_string())?;
    Ok(writer.bytes)
}
pub fn decode(bytes: &[u8]) -> Result<Draft, String> {
    decode_expected(bytes, None)
}
pub fn decode_version(bytes: &[u8], expected: u32) -> Result<Draft, String> {
    decode_expected(bytes, Some(expected))
}
fn decode_expected(bytes: &[u8], expected: Option<u32>) -> Result<Draft, String> {
    if bytes.len() > MAX_BYTES {
        return Err("Assembly draft exceeds32MiB".into());
    }
    let mut draft: Draft = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
    if expected.is_some_and(|version| version != draft.version) {
        return Err("Assembly attachment kind disagrees with its draft version.".into());
    }
    draft.validate()?;
    draft.version = 3;
    Ok(draft)
}

pub(crate) fn digest_interval(
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
            slot.end = source.molecule.sequence.len();
            slot.captured = Some(source);
        }
        draft.sync_rows();
        draft
    }

    #[test]
    fn two_source_tailing_roundtrip_preserves_pins_roles_preparations_and_literal_intent() {
        let mut draft = draft();
        draft.slots[1].reverse = true;
        draft.methods.rows[0].preparation = OverlapPreparation::Digest;
        draft.methods.rows[0].enzymes = "EcoRV".into();
        draft.methods.rows[0].fragment = "2".into();
        draft.methods.rows[1].forward = "GG".into();
        draft.methods.rows[1].reverse_primer = "AA".into();
        draft.methods.vector.base = TailingBase::C;
        draft.methods.insert.base = TailingBase::G;
        draft.methods.insert.treatment = Treatment::Phosphorylate;
        draft.junction_edits.push(assembly_junction_edit::Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "T".into(),
        });
        draft.sync_rows();
        let bytes = encode_raw(&draft).unwrap();
        let restored = decode_version(&bytes, 3).unwrap();
        assert!(restored == draft);
        assert_eq!(restored.methods.vector.document, 0);
        assert_eq!(restored.methods.insert.document, 1);
        assert!(restored.methods.insert.reversed);
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
        assert!(decode_version(&encode_raw(&draft).unwrap(), 3).unwrap() == draft);
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
        draft.slots[0].reverse = true;
        draft.methods.rows[0].preparation = OverlapPreparation::Digest;
        draft.methods.rows[0].enzymes = "EcoRV".into();
        draft.methods.rows[1].forward = "GG".into();
        draft.methods.vector.base = TailingBase::C;
        draft.methods.vector.treatment = Treatment::Dephosphorylate;
        draft.methods.insert.base = TailingBase::G;
        draft.methods.insert.treatment = Treatment::Phosphorylate;
        draft.sync_rows();
        let original = draft.clone();
        draft.slots.swap(0, 1);
        draft.methods.swap_rows(0, 1);
        draft.sync_rows();
        draft.validate().unwrap();
        assert!(draft.slots[0].captured == original.slots[1].captured);
        assert!(draft.slots[1].captured == original.slots[0].captured);
        assert!(draft.methods.vector.base == TailingBase::G);
        assert!(draft.methods.vector.treatment == Treatment::Phosphorylate);
        assert!(draft.methods.insert.base == TailingBase::C);
        assert!(draft.methods.insert.treatment == Treatment::Dephosphorylate);
        assert_eq!(draft.methods.rows[0].forward, "GG");
        assert_eq!(draft.methods.rows[1].enzymes, "EcoRV");
        assert!(!draft.methods.vector.reversed && draft.methods.insert.reversed);
        assert!(
            draft.methods.source_roles.as_ref().unwrap()[0]
                == original.methods.source_roles.as_ref().unwrap()[1]
        );
        draft.methods.swap_rows(0, 0);
        assert!(draft.methods.vector.base == TailingBase::G);
        draft.slots.swap(1, 0);
        draft.methods.swap_rows(1, 0);
        draft.sync_rows();
        assert!(draft == original);
    }

    #[test]
    fn tailing_rejects_wrong_count_source_drift_legacy_versions_and_library_rows() {
        let mut invalid = draft();
        invalid.methods.insert.document = 0;
        assert!(invalid.validate().is_err());
        invalid.sync_rows();
        invalid.methods.insert.reversed = true;
        assert!(invalid.validate().is_err());
        invalid.sync_rows();
        invalid.version = 2;
        assert!(invalid.validate().is_err());
        invalid.version = 3;
        invalid.slots.push(Slot::empty("Third source"));
        invalid.sync_rows();
        assert!(invalid.validate().is_err());

        let mut library =
            crate::assembly_alternatives::enable(&Draft::default(), "ta-exclusion").unwrap();
        let saved = &mut library.alternatives.as_mut().unwrap().rows[0].edits.methods;
        saved.mode = Mode::Tailing;
        saved.insert.document = 1;
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
            name: draft.name.clone(),
            circular: true,
            fragments: draft
                .slots
                .iter()
                .zip(&roles)
                .map(
                    |(slot, role)| collagene_core::overlap_primer_design::Fragment {
                        role: role.role.clone(),
                        input_id: role.input_id.clone(),
                        start: slot.start,
                        end: slot.end,
                        reverse: slot.reverse,
                    },
                )
                .collect(),
        };
        let target =
            assembly_junction_edit::build_intent(&inputs, &base, &draft.junction_edits).unwrap();
        assert_eq!(target.target.sequence.as_str(), "AACCTGGTTA");
        draft.methods.vector.base = TailingBase::C;
        draft.methods.insert.base = TailingBase::G;
        draft.sync_rows();
        let unchanged =
            assembly_junction_edit::build_intent(&inputs, &base, &draft.junction_edits).unwrap();
        assert_eq!(target.context_revision, unchanged.context_revision);
        assert_eq!(target.target, unchanged.target);
    }
}
