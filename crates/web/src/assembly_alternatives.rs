//! Immutable source choices plus editable per-construct assembly settings.
//! Working source edits persist independently and need explicit choice capture
//! before they can acquire a construct identity or enter scientific preparation.
use crate::{
    assembly_draft_model::{self as assembly_workspace, CapturedSource, Draft, Slot},
    assembly_role_ids::Role,
    cloning_methods_model::{MethodsState, OverlapInput},
};
use collagene_core::assembly_junction_edit;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub const MAX_BINS: usize = 12;
pub const MAX_CHOICES_PER_BIN: usize = 64;
pub const MAX_CHOICES: usize = 256;
pub const MAX_ROWS: usize = 256;

#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Choice {
    pub id: String,
    pub slot: Slot,
    pub row: OverlapInput,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bin {
    pub id: String,
    pub label: String,
    #[serde(deserialize_with = "deserialize_choices")]
    pub choices: Vec<Choice>,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Edits {
    pub name: String,
    pub circular: bool,
    pub methods: MethodsState,
    #[serde(deserialize_with = "deserialize_junction_edits")]
    pub junction_edits: Vec<assembly_junction_edit::Edit>,
    #[serde(deserialize_with = "deserialize_substitutions")]
    pub substitutions: Vec<collagene_core::assembly_substitution_intent::Substitution>,
    #[serde(deserialize_with = "deserialize_internal_edits")]
    pub edits: Vec<collagene_core::assembly_edit_intent::Edit>,
    pub active_slot: usize,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Construct {
    pub id: String,
    #[serde(deserialize_with = "deserialize_choice_ids")]
    pub choices: Vec<String>,
    pub edits: Edits,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Library {
    /// Short design identifier; entity counters remain bounded native IDs.
    pub id: String,
    /// All entity kinds share this monotonic, checked allocation counter.
    pub next_id: u64,
    #[serde(deserialize_with = "deserialize_bins")]
    pub bins: Vec<Bin>,
    #[serde(deserialize_with = "deserialize_constructs")]
    pub rows: Vec<Construct>,
    pub active_row: String,
}

fn bounded<'de, D, T, const N: usize>(deserializer: D) -> Result<Vec<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    struct Visitor<T, const N: usize>(std::marker::PhantomData<T>);
    impl<'de, T: Deserialize<'de>, const N: usize> serde::de::Visitor<'de> for Visitor<T, N> {
        type Value = Vec<T>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "at most {N} assembly alternatives")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(self, mut seq: A) -> Result<Vec<T>, A::Error> {
            if seq.size_hint().is_some_and(|n| n > N) {
                return Err(serde::de::Error::custom("assembly alternative count limit"));
            }
            let mut result = Vec::new();
            while result.len() < N {
                match seq.next_element()? {
                    Some(value) => result.push(value),
                    None => return Ok(result),
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom("assembly alternative count limit"));
            }
            Ok(result)
        }
    }
    deserializer.deserialize_seq(Visitor::<T, N>(std::marker::PhantomData))
}
macro_rules! decoder {
    ($name:ident,$ty:ty,$n:expr) => {
        fn $name<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<$ty>, D::Error> {
            bounded::<D, $ty, $n>(d)
        }
    };
}
decoder!(deserialize_choices, Choice, MAX_CHOICES_PER_BIN);
decoder!(deserialize_bins, Bin, MAX_BINS);
decoder!(deserialize_constructs, Construct, MAX_ROWS);
decoder!(deserialize_choice_ids, String, MAX_BINS);
decoder!(
    deserialize_junction_edits,
    assembly_junction_edit::Edit,
    MAX_BINS
);

decoder!(
    deserialize_substitutions,
    collagene_core::assembly_substitution_intent::Substitution,
    256
);
decoder!(
    deserialize_internal_edits,
    collagene_core::assembly_edit_intent::Edit,
    256
);

fn token(text: &str, max: usize) -> bool {
    !text.trim().is_empty() && text.len() <= max && !text.chars().any(char::is_control)
}
fn design_id(id: &str) -> bool {
    !id.is_empty()
        && id.len() <= 16
        && id
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.'))
}
fn entity_id(library: &Library, id: &str, kind: &str) -> bool {
    id.len() <= 48
        && id
            .strip_prefix(&format!("{}-{kind}-", library.id))
            .and_then(|n| n.parse::<u64>().ok())
            .is_some_and(|n| n > 0 && n < library.next_id)
}
fn allocate(library: &mut Library, kind: &str) -> Result<String, String> {
    let value = library.next_id;
    library.next_id = value
        .checked_add(1)
        .ok_or("The assembly alternative identifier counter is exhausted.")?;
    Ok(format!("{}-{kind}-{value}", library.id))
}
fn active(library: &Library) -> Result<&Construct, String> {
    library
        .rows
        .iter()
        .find(|r| r.id == library.active_row)
        .ok_or_else(|| "The active construct row is unavailable.".into())
}
fn active_mut(library: &mut Library) -> Result<&mut Construct, String> {
    library
        .rows
        .iter_mut()
        .find(|r| r.id == library.active_row)
        .ok_or_else(|| "The active construct row is unavailable.".into())
}
fn chosen<'a>(library: &'a Library, choices: &[String]) -> Result<Vec<&'a Choice>, String> {
    if choices.len() != library.bins.len() {
        return Err("Choose exactly one source alternative from every bin.".into());
    }
    library
        .bins
        .iter()
        .zip(choices)
        .map(|(bin, id)| {
            bin.choices
                .iter()
                .find(|c| c.id == *id)
                .ok_or_else(|| "A selected alternative does not belong to its ordered bin.".into())
        })
        .collect()
}
fn roles(library: &Library, row: &Construct) -> Result<Vec<Role>, String> {
    chosen(library, &row.choices)?;
    library
        .bins
        .iter()
        .zip(&row.choices)
        .map(|(bin, choice)| {
            let (input_id, role) = collagene_core::combinatorial_design::source_role_ids(
                &library.id,
                &row.id,
                &bin.id,
                choice,
            )
            .map_err(|e| e.to_string())?;
            Ok(Role { input_id, role })
        })
        .collect()
}
pub fn roles_for_active(draft: &Draft) -> Result<Option<Vec<Role>>, String> {
    draft
        .alternatives
        .as_ref()
        .map(|library| roles(library, active(library)?))
        .transpose()
}
pub fn identity(draft: &Draft) -> Option<(String, String, Vec<String>)> {
    let library = draft.alternatives.as_ref()?;
    let row = active(library).ok()?;
    Some((library.id.clone(), row.id.clone(), row.choices.clone()))
}

fn same_geometry(a: &Slot, b: &Slot) -> bool {
    a.captured == b.captured
        && a.selector == b.selector
        && a.digest_geometry == b.digest_geometry
        && a.material_origin == b.material_origin
}
/// Source choice authority is explicit; an unfinished working source edit may be
/// persisted, but must not silently inherit the previously chosen source ID.
pub fn require_active_choice_match(draft: &Draft) -> Result<(), String> {
    let Some(library) = &draft.alternatives else {
        return Ok(());
    };
    let row = active(library)?;
    let choices = chosen(library, &row.choices)?;
    if choices.len() != draft.slots.len()
        || draft
            .slots
            .iter()
            .zip(choices)
            .any(|(slot, choice)| !same_geometry(slot, &choice.slot))
    {
        return Err("Working source, interval, strand or digest geometry differs from the selected alternatives. Explicitly capture each changed source as a new choice before saving or switching constructs or running preparation.".into());
    }
    if draft.methods.source_roles.as_ref() != Some(&roles(library, row)?) {
        return Err("Working preparation identities differ from the active construct and selected alternatives.".into());
    }
    Ok(())
}

impl Library {
    pub fn validate(&self) -> Result<(), String> {
        if !design_id(&self.id)
            || self.next_id == 0
            || self.bins.is_empty()
            || self.bins.len() > MAX_BINS
            || self.rows.is_empty()
            || self.rows.len() > MAX_ROWS
        {
            return Err("Invalid bounded assembly alternative library.".into());
        }
        assembly_workspace::bounded_size(self)?;
        let mut ids = BTreeSet::new();
        let mut total = 0usize;
        for (index, bin) in self.bins.iter().enumerate() {
            if !entity_id(self, &bin.id, "bin")
                || !ids.insert(&bin.id)
                || !token(&bin.label, 80)
                || bin.choices.is_empty()
                || bin.choices.len() > MAX_CHOICES_PER_BIN
            {
                return Err("Invalid assembly bin identity, label or choice count.".into());
            }
            total = total
                .checked_add(bin.choices.len())
                .filter(|n| *n <= MAX_CHOICES)
                .ok_or("Assembly alternatives exceed 256 retained choices.")?;
            for choice in &bin.choices {
                if !entity_id(self, &choice.id, "choice") || !ids.insert(&choice.id) {
                    return Err("Invalid or duplicate assembly choice identity.".into());
                }
                assembly_workspace::validate_slot(&choice.slot, &choice.row, index)?;
            }
        }
        for row in &self.rows {
            if !entity_id(self, &row.id, "row") || !ids.insert(&row.id) {
                return Err("Invalid or duplicate construct row identity.".into());
            }
            let selected = chosen(self, &row.choices)?;
            assembly_workspace::validate_controls(
                &row.edits.name,
                row.edits.circular,
                row.edits.active_slot,
                self.bins.len(),
                &row.edits.methods,
                &row.edits.junction_edits,
            )?;
            assembly_workspace::validate_internal_declarations(
                &row.edits.methods,
                self.bins.len(),
                &row.edits.substitutions,
                &row.edits.edits,
            )?;
            for (index, (choice, input)) in selected
                .into_iter()
                .zip(&row.edits.methods.rows)
                .enumerate()
            {
                assembly_workspace::validate_slot_link(&choice.slot, input, index)?;
            }
            if row.edits.methods.source_roles.as_ref() != Some(&roles(self, row)?) {
                return Err("Saved construct preparation roles differ from its stable selected alternatives.".into());
            }
        }
        active(self)?;
        Ok(())
    }
}

fn edits(draft: &Draft) -> Edits {
    Edits {
        name: draft.name.clone(),
        circular: draft.circular,
        methods: draft.methods.clone(),
        junction_edits: draft.junction_edits.clone(),
        substitutions: draft.substitutions.clone(),
        edits: draft.edits.clone(),
        active_slot: draft.active_slot,
    }
}
fn stage(draft: &Draft) -> Result<Draft, String> {
    draft.validate()?;
    Ok(draft.clone())
}
fn finish(mut draft: Draft) -> Result<Draft, String> {
    if draft.alternatives.is_some() {
        draft.methods.source_roles = roles_for_active(&draft)?;
    }
    draft.validate()?;
    Ok(draft)
}
fn save_inner(draft: &mut Draft) -> Result<(), String> {
    require_active_choice_match(draft)?;
    if draft.alternatives.is_some() {
        draft.methods.source_roles = roles_for_active(draft)?;
    }
    let saved = edits(draft);
    let labels = draft
        .slots
        .iter()
        .map(|s| s.label.clone())
        .collect::<Vec<_>>();
    if let Some(library) = &mut draft.alternatives {
        for (bin, label) in library.bins.iter_mut().zip(labels) {
            bin.label = label;
        }
        active_mut(library)?.edits = saved;
    }
    Ok(())
}
pub fn enable(draft: &Draft, id: &str) -> Result<Draft, String> {
    if draft.alternatives.is_some() {
        return Err("This assembly already has alternative bins.".into());
    }
    if !design_id(id) {
        return Err("Use a nonempty ASCII assembly design identifier of at most 16 bytes.".into());
    }
    let mut next = stage(draft)?;
    let mut library = Library {
        id: id.into(),
        next_id: 1,
        bins: vec![],
        rows: vec![],
        active_row: String::new(),
    };
    for (slot, row) in next.slots.iter().zip(&next.methods.rows) {
        let bin = allocate(&mut library, "bin")?;
        let choice = allocate(&mut library, "choice")?;
        library.bins.push(Bin {
            id: bin,
            label: slot.label.clone(),
            choices: vec![Choice {
                id: choice,
                slot: slot.clone(),
                row: row.clone(),
            }],
        });
    }
    let row_id = allocate(&mut library, "row")?;
    let mut row = Construct {
        id: row_id.clone(),
        choices: library
            .bins
            .iter()
            .map(|b| b.choices[0].id.clone())
            .collect(),
        edits: edits(&next),
    };
    let old_roles = next.methods.role_ids(next.slots.len())?;
    let new_roles = roles(&library, &row)?;
    for (old, new) in old_roles.iter().zip(&new_roles) {
        for edit in &mut row.edits.substitutions {
            if edit.role == old.role {
                edit.role = new.role.clone();
            }
        }
        for edit in &mut row.edits.edits {
            if edit.role == old.role {
                edit.role = new.role.clone();
            }
        }
    }
    next.substitutions = row.edits.substitutions.clone();
    next.edits = row.edits.edits.clone();
    row.edits.methods.source_roles = Some(new_roles);
    library.rows.push(row);
    library.active_row = row_id;
    next.version = crate::assembly_draft_model::CURRENT_VERSION;
    next.alternatives = Some(library);
    finish(next)
}
pub fn save_active(draft: &Draft) -> Result<Draft, String> {
    let mut next = stage(draft)?;
    save_inner(&mut next)?;
    finish(next)
}
/// Save an explicit candidate declaration for the active construct. The caller
/// binds its asynchronous inspection to the unchanged draft before invoking
/// this helper; native compilation must reconstruct and check the choice again.
pub fn set_active_pool_candidate(
    draft: &Draft,
    candidate: Option<crate::golden_gate_pool_form::CandidateIntent>,
) -> Result<Draft, String> {
    let mut next = stage(draft)?;
    require_active_choice_match(&next)?;
    if next.alternatives.is_none()
        || next.methods.mode != crate::cloning_methods_model::Mode::GoldenGate
    {
        return Err("Choose a Golden Gate construct before declaring a pool candidate.".into());
    }
    if let Some(intent) = &candidate {
        intent.validate()?;
    }
    let pool = next
        .methods
        .golden_gate_pool
        .as_mut()
        .ok_or("Inspect explicit pool settings before choosing a candidate.")?;
    pool.candidate = candidate;
    save_inner(&mut next)?;
    next.version = crate::assembly_draft_model::CURRENT_VERSION;
    finish(next)
}
/// Persist raw working edits even while a source edit awaits explicit choice
/// capture. Existing construct edits are refreshed only when their sources match.
pub fn for_storage(draft: &Draft) -> Result<Draft, String> {
    let mut next = stage(draft)?;
    if require_active_choice_match(&next).is_ok() {
        save_inner(&mut next)?;
    }
    next.version = crate::assembly_draft_model::CURRENT_VERSION;
    finish(next)
}
pub fn capture_active_choice(draft: &Draft) -> Result<Draft, String> {
    let mut next = stage(draft)?;
    let index = next.active_slot;
    // A new choice has a new role; old source literals must never move onto it.
    next.clear_slot_substitutions(index)?;
    let choice = Choice {
        id: String::new(),
        slot: next.slots[index].clone(),
        row: next.methods.rows[index].clone(),
    };
    let library = next
        .alternatives
        .as_mut()
        .ok_or("Enable source alternatives first.")?;
    if library.bins[index].choices.len() >= MAX_CHOICES_PER_BIN
        || library.bins.iter().map(|b| b.choices.len()).sum::<usize>() >= MAX_CHOICES
    {
        return Err("Assembly alternative choice limit reached.".into());
    }
    let id = allocate(library, "choice")?;
    library.bins[index].choices.push(Choice {
        id: id.clone(),
        ..choice
    });
    let saved = active_mut(library)?;
    saved.choices[index] = id;
    saved.edits.methods.rows[index] = next.methods.rows[index].clone();
    let expected = roles(library, active(library)?)?;
    active_mut(library)?.edits.methods.source_roles = Some(expected);
    next.methods.source_roles = roles_for_active(&next)?;
    if require_active_choice_match(&next).is_ok() {
        save_inner(&mut next)?;
    }
    finish(next)
}
fn apply_row(draft: &mut Draft, id: &str) -> Result<(), String> {
    let library = draft
        .alternatives
        .as_ref()
        .ok_or("Enable source alternatives first.")?;
    let row = library
        .rows
        .iter()
        .find(|r| r.id == id)
        .ok_or("The selected construct row is unavailable.")?;
    let slots = chosen(library, &row.choices)?
        .into_iter()
        .zip(&library.bins)
        .map(|(choice, bin)| {
            let mut slot = choice.slot.clone();
            slot.label = bin.label.clone();
            slot
        })
        .collect();
    let saved = row.edits.clone();
    draft.name = saved.name;
    draft.circular = saved.circular;
    draft.methods = saved.methods;
    draft.junction_edits = saved.junction_edits;
    draft.substitutions = saved.substitutions;
    draft.edits = saved.edits;
    draft.active_slot = saved.active_slot;
    draft.slots = slots;
    draft
        .alternatives
        .as_mut()
        .expect("checked above")
        .active_row = id.into();
    Ok(())
}
pub fn open_row(draft: &Draft, id: &str) -> Result<Draft, String> {
    let mut next = save_active(draft)?;
    apply_row(&mut next, id)?;
    finish(next)
}
pub fn add_row(draft: &Draft, selection: Vec<String>) -> Result<Draft, String> {
    let mut next = save_active(draft)?;
    let mut saved = edits(&next);
    saved.name = format!(
        "Construct {}",
        next.alternatives
            .as_ref()
            .ok_or("Enable source alternatives first.")?
            .rows
            .len()
            + 1
    );
    saved.junction_edits.clear();
    saved.substitutions.clear();
    saved.edits.clear();
    let library = next.alternatives.as_mut().expect("checked above");
    if library.rows.len() >= MAX_ROWS {
        return Err("The selected construct matrix is limited to 256 rows.".into());
    }
    let selected = chosen(library, &selection)?;
    saved.methods.rows = selected
        .into_iter()
        .map(|choice| choice.row.clone())
        .collect();
    let id = allocate(library, "row")?;
    let mut row = Construct {
        id: id.clone(),
        choices: selection,
        edits: saved,
    };
    row.edits.methods.source_roles = Some(roles(library, &row)?);
    library.rows.push(row);
    apply_row(&mut next, &id)?;
    finish(next)
}
pub fn remove_choice(draft: &Draft, id: &str) -> Result<Draft, String> {
    let mut next = stage(draft)?;
    let library = next
        .alternatives
        .as_mut()
        .ok_or("Enable source alternatives first.")?;
    if library
        .rows
        .iter()
        .any(|row| row.choices.iter().any(|choice| choice == id))
    {
        return Err(
            "A construct still uses this source choice; remove or revise that row first.".into(),
        );
    }
    let bin = library
        .bins
        .iter_mut()
        .find(|bin| bin.choices.iter().any(|choice| choice.id == id))
        .ok_or("The source choice is unavailable.")?;
    if bin.choices.len() == 1 {
        return Err("Every alternative bin needs at least one source choice.".into());
    }
    bin.choices.retain(|choice| choice.id != id);
    finish(next)
}
pub fn remove_row(draft: &Draft, id: &str) -> Result<Draft, String> {
    let mut next = stage(draft)?;
    let library = next
        .alternatives
        .as_ref()
        .ok_or("Enable source alternatives first.")?;
    if library.rows.len() == 1 {
        return Err("Keep at least one construct row.".into());
    }
    if !library.rows.iter().any(|row| row.id == id) {
        return Err("The construct row is unavailable.".into());
    }
    let is_active = library.active_row == id;
    if is_active {
        save_inner(&mut next)?;
    }
    let library = next.alternatives.as_mut().expect("checked above");
    library.rows.retain(|row| row.id != id);
    if is_active {
        let replacement = library.rows[0].id.clone();
        apply_row(&mut next, &replacement)?;
    }
    finish(next)
}
// Matrix display carries only bounded stable IDs. It neither copies nor attests
// source models; mutation/preparation still uses the complete draft validation.
fn matrix_counts(library: &Library) -> Result<Vec<usize>, String> {
    if !design_id(&library.id)
        || library.next_id == 0
        || library.bins.is_empty()
        || library.bins.len() > MAX_BINS
    {
        return Err("Invalid bounded alternative matrix.".into());
    }
    let mut ids = BTreeSet::new();
    let mut counts = Vec::new();
    let mut total = 0usize;
    for bin in &library.bins {
        if !entity_id(library, &bin.id, "bin")
            || !ids.insert(&bin.id)
            || bin.choices.is_empty()
            || bin.choices.len() > MAX_CHOICES_PER_BIN
        {
            return Err("Invalid alternative matrix bin.".into());
        }
        total = total
            .checked_add(bin.choices.len())
            .filter(|n| *n <= MAX_CHOICES)
            .ok_or("Alternative matrix choice limit exceeded.")?;
        for choice in &bin.choices {
            if !entity_id(library, &choice.id, "choice") || !ids.insert(&choice.id) {
                return Err("Invalid alternative matrix choice identity.".into());
            }
        }
        counts.push(bin.choices.len());
    }
    Ok(counts)
}
pub fn matrix_size(library: &Library) -> Result<u128, String> {
    collagene_core::combinatorial_design::cartesian_size(&matrix_counts(library)?)
        .map_err(|e| e.to_string())
}
pub fn combinations_page(
    library: &Library,
    offset: u128,
    limit: usize,
) -> Result<Vec<Vec<String>>, String> {
    collagene_core::combinatorial_design::cartesian_page(&matrix_counts(library)?, offset, limit)
        .map_err(|e| e.to_string())
        .map(|rows| {
            rows.into_iter()
                .map(|row| {
                    row.into_iter()
                        .enumerate()
                        .map(|(i, choice)| library.bins[i].choices[choice].id.clone())
                        .collect()
                })
                .collect()
        })
}
pub fn all_captures(draft: &Draft) -> impl Iterator<Item = &CapturedSource> {
    draft
        .slots
        .iter()
        .flat_map(|slot| {
            [slot.captured.as_deref(), slot.anneal_bottom.as_deref()]
                .into_iter()
                .flatten()
        })
        .chain(
            draft
                .alternatives
                .iter()
                .flat_map(|l| &l.bins)
                .flat_map(|bin| &bin.choices)
                .flat_map(|choice| {
                    [
                        choice.slot.captured.as_deref(),
                        choice.slot.anneal_bottom.as_deref(),
                    ]
                    .into_iter()
                    .flatten()
                }),
        )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        assembly_draft_model::{self, Slot},
        source_identity::SourceInfo,
    };
    use collagene_core::{Molecule, Topology, run_lineage::EventBinding, workspace::Workspace};
    use std::sync::Arc;

    fn source(name: &str, dna: &str) -> Arc<CapturedSource> {
        let molecule = Molecule::dna(name, dna, Topology::Linear).unwrap();
        let frame = Workspace::from_documents("source", vec![molecule.clone()], 0).unwrap();
        let doc = &frame.documents[0];
        Arc::new(CapturedSource {
            source: SourceInfo {
                binding: EventBinding {
                    document_id: doc.id.clone(),
                    event_id: doc.head_event.clone().unwrap(),
                    revision: doc.head.clone(),
                },
                name: molecule.name.clone(),
                length: molecule.sequence.len(),
                alphabet: molecule.sequence.alphabet(),
                topology: molecule.topology,
            },
            molecule,
        })
    }
    fn ordinary() -> Draft {
        let mut draft = Draft::default();
        draft.name = "First construct".into();
        draft.circular = false;
        for (slot, captured) in draft
            .slots
            .iter_mut()
            .zip([source("backbone", "AACCGGTT"), source("insert", "TTGCAACC")])
        {
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    start: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = 0;
            }
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    end: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = captured.molecule.sequence.len();
            }
            slot.captured = Some(captured);
        }
        draft.sync_rows();
        draft.validate().unwrap();
        draft
    }
    fn replace_insert(draft: &mut Draft) {
        let captured = source("alternative", "CCGTAAGG");
        draft.active_slot = 1;
        draft.slots[1] = Slot {
            material_origin: None,
            selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                start: 0,
                end: captured.molecule.sequence.len(),
                reverse: false,
            },
            label: "Insert".into(),
            digest_geometry: None,
            anneal_bottom: None,
            captured: Some(captured),
        };
        draft.sync_rows();
    }

    #[test]
    fn row_local_substitutions_and_indels_survive_switch_restart_reorder_and_recovery() {
        let first = enable(&ordinary(), "literal-rows").unwrap();
        let first_id = first.alternatives.as_ref().unwrap().active_row.clone();
        let choices = first.alternatives.as_ref().unwrap().rows[0].choices.clone();
        let mut second = add_row(&first, choices).unwrap();
        let second_id = second.alternatives.as_ref().unwrap().active_row.clone();
        second = crate::assembly_edit_model::declare(&second, 0, "3", "5", "TAG", None).unwrap();
        second = crate::assembly_substitution_model::declare(&second, 1, "1", "A", None).unwrap();
        let target = |d: &Draft| {
            let (inputs, base, _) = assembly_draft_model::inputs_request(d, None).unwrap();
            collagene_core::assembly_edit_intent::build(
                &inputs,
                &assembly_draft_model::edit_request(d, base),
                None,
            )
            .unwrap()
            .target
            .sequence
            .as_str()
            .to_owned()
        };
        assert_eq!(target(&second), "AATAGGGTTATGCAACC");
        let first_again = open_row(&second, &first_id).unwrap();
        assert!(first_again.edits.is_empty() && first_again.substitutions.is_empty());
        assert_eq!(target(&first_again), "AACCGGTTTTGCAACC");
        let restored =
            assembly_draft_model::decode(&assembly_draft_model::encode(&first_again).unwrap())
                .unwrap();
        let second_again = open_row(&restored, &second_id).unwrap();
        assert_eq!(target(&second_again), "AATAGGGTTATGCAACC");
        let bins = second_again
            .alternatives
            .as_ref()
            .unwrap()
            .bins
            .iter()
            .rev()
            .map(|b| b.id.clone())
            .collect();
        let reorder = crate::assembly_library_topology::propose(
            &second_again,
            crate::assembly_library_topology::Change::Reorder { bin_ids: bins },
        )
        .unwrap();
        assert!(
            reorder
                .impacts()
                .iter()
                .all(|i| i.recovered_internal_edits == 0)
        );
        let reordered =
            crate::assembly_library_topology::apply(&second_again, &reorder, reorder.revision())
                .unwrap();
        assert_eq!(target(&reordered.draft), "ATGCAACCAATAGGGTT");
        let removed_id = reordered.draft.alternatives.as_ref().unwrap().bins[1]
            .id
            .clone();
        let removal = crate::assembly_library_topology::propose(
            &reordered.draft,
            crate::assembly_library_topology::Change::Remove { bin_id: removed_id },
        )
        .unwrap();
        assert_eq!(
            removal
                .impacts()
                .iter()
                .find(|r| r.construct_id == second_id)
                .unwrap()
                .recovered_internal_edits,
            1
        );
        let applied =
            crate::assembly_library_topology::apply(&reordered.draft, &removal, removal.revision())
                .unwrap();
        assert_eq!(applied.recovery.edits, second.edits);
        assert!(applied.draft.edits.is_empty());
        assert_eq!(applied.draft.substitutions.len(), 1);
    }
    #[test]
    fn source_replacement_and_new_rows_never_rebind_internal_literals() {
        let first = enable(&ordinary(), "new-source").unwrap();
        let first = crate::assembly_edit_model::declare(&first, 1, "3", "5", "TAG", None).unwrap();
        let choices = first.alternatives.as_ref().unwrap().rows[0].choices.clone();
        let second = add_row(&first, choices).unwrap();
        assert!(second.edits.is_empty() && second.substitutions.is_empty());
        let first_id = first.alternatives.as_ref().unwrap().active_row.clone();
        let mut working = open_row(&second, &first_id).unwrap();
        replace_insert(&mut working);
        assert!(crate::assembly_edit_model::declare(&working, 1, "1", "2", "A", None).is_err());
        let replaced = capture_active_choice(&working).unwrap();
        assert!(replaced.edits.is_empty());
        let recovered =
            open_row(&replaced, &second.alternatives.as_ref().unwrap().active_row).unwrap();
        assert!(recovered.edits.is_empty());
        let mut missing = serde_json::to_value(&replaced).unwrap();
        missing["alternatives"]["rows"][0]["edits"]
            .as_object_mut()
            .unwrap()
            .remove("edits");
        assert!(assembly_draft_model::decode(&serde_json::to_vec(&missing).unwrap()).is_err());
    }
    #[test]
    fn source_aliases_and_separate_constructs_receive_stable_distinct_recipe_identities() {
        let first = enable(&ordinary(), "test-library").unwrap();
        let original = first.methods.source_roles.clone().unwrap();
        let library = first.alternatives.as_ref().unwrap();
        let selection = library.rows[0].choices.clone();
        let second = add_row(&first, selection).unwrap();
        assert_ne!(identity(&first), identity(&second));
        assert_ne!(original, second.methods.source_roles.clone().unwrap());
        assert!(first.slots == second.slots);
        assert_eq!(
            matrix_size(second.alternatives.as_ref().unwrap()).unwrap(),
            1
        );
        require_active_choice_match(&second).unwrap();
    }

    #[test]
    fn per_construct_unfinished_primers_and_junction_edits_survive_switching_and_storage() {
        let mut first = enable(&ordinary(), "per-row").unwrap();
        first.methods.rows[0].forward = "ag N ? unfinished".into();
        let first_id = first.alternatives.as_ref().unwrap().active_row.clone();
        let choices = first.alternatives.as_ref().unwrap().rows[0].choices.clone();
        let mut second = add_row(&first, choices).unwrap();
        second.methods.rows[0].forward = "different manual primer".into();
        second.junction_edits = vec![assembly_junction_edit::Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "AC".into(),
        }];
        let second_id = second.alternatives.as_ref().unwrap().active_row.clone();
        let restored =
            assembly_draft_model::decode(&assembly_draft_model::encode(&second).unwrap()).unwrap();
        let first_again = open_row(&restored, &first_id).unwrap();
        assert_eq!(first_again.methods.rows[0].forward, "ag N ? unfinished");
        assert!(first_again.junction_edits.is_empty());
        let second_again = open_row(&first_again, &second_id).unwrap();
        assert_eq!(
            second_again.methods.rows[0].forward,
            "different manual primer"
        );
        assert_eq!(second_again.junction_edits[0].inserted, "AC");
    }

    #[test]
    fn pending_source_edits_persist_but_cannot_inherit_an_old_choice_identity() {
        let original = enable(&ordinary(), "pending").unwrap();
        let mut working = original.clone();
        replace_insert(&mut working);
        working.methods.rows[1].forward = "keep this pending text".into();
        working.validate().unwrap();
        assert!(require_active_choice_match(&working).is_err());
        assert!(save_active(&working).is_err());
        let bytes = assembly_draft_model::encode(&working).unwrap();
        let restored = assembly_draft_model::decode(&bytes).unwrap();
        assert_eq!(
            restored.slots[1].captured.as_ref().unwrap().source.name,
            "alternative"
        );
        assert_eq!(restored.methods.rows[1].forward, "keep this pending text");
        assert!(require_active_choice_match(&restored).is_err());
        let captured = capture_active_choice(&restored).unwrap();
        require_active_choice_match(&captured).unwrap();
        assert_ne!(identity(&original), identity(&captured));
        assert_eq!(
            captured.alternatives.as_ref().unwrap().bins[1]
                .choices
                .len(),
            2
        );
        assert_eq!(captured.methods.rows[1].forward, "keep this pending text");
    }

    #[test]
    fn working_identity_tampering_and_unknown_choice_rows_fail_without_rebinding() {
        let first = enable(&ordinary(), "identity").unwrap();
        let mut bad = first.clone();
        bad.methods.source_roles.as_mut().unwrap()[0].input_id =
            "valid-but-unrelated-source-id".into();
        assert!(bad.validate().is_err());
        assert!(require_active_choice_match(&bad).is_err());
        assert!(assembly_draft_model::encode(&bad).is_err());
        let mut bad = first.clone();
        bad.alternatives.as_mut().unwrap().rows[0].choices[0] = "missing-choice".into();
        assert!(bad.validate().is_err());
        let mut bad = first.clone();
        let duplicate = bad.alternatives.as_ref().unwrap().rows[0].clone();
        bad.alternatives.as_mut().unwrap().rows.push(duplicate);
        assert!(bad.validate().is_err());
        assert!(first.alternatives.as_ref().unwrap().rows.len() == 1);
    }

    #[test]
    fn lazy_matrix_and_removal_keep_choices_stable_and_protect_referenced_sources() {
        let mut draft = enable(&ordinary(), "matrix").unwrap();
        let initial_id = draft.alternatives.as_ref().unwrap().bins[1].choices[0]
            .id
            .clone();
        assert!(remove_choice(&draft, &initial_id).is_err());
        replace_insert(&mut draft);
        draft = capture_active_choice(&draft).unwrap();
        let library = draft.alternatives.as_ref().unwrap();
        assert_eq!(matrix_size(library).unwrap(), 2);
        let first = combinations_page(library, 0, 1).unwrap();
        let second = combinations_page(library, 1, 1).unwrap();
        assert_eq!(first[0][1], initial_id);
        assert_ne!(first, second);
        assert!(combinations_page(library, 2, 1).unwrap().is_empty());
        let old_counter = library.next_id;
        let removed = remove_choice(&draft, &initial_id).unwrap();
        assert_eq!(removed.alternatives.as_ref().unwrap().next_id, old_counter);
        assert!(require_active_choice_match(&removed).is_ok());
        assert!(removed.slots == draft.slots);
    }
}
