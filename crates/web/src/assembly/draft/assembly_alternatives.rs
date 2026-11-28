//! Immutable source choices plus editable per-construct assembly settings.
//! Working source edits persist independently and need explicit choice capture
//! before they can acquire a construct identity or enter scientific preparation.
use crate::assembly::draft::assembly_role_ids::Role;
use crate::assembly::draft::model;
use crate::assembly::draft::model::CapturedSource;
use crate::assembly::draft::model::Draft;
use crate::assembly::draft::model::Slot;
use crate::cloning::preparation::cloning_methods_model::MethodsState;
use crate::cloning::preparation::cloning_methods_model::OverlapInput;
use crate::cloning::preparation::cloning_methods_model::OverlapPreparation;
use collagene_core::cloning::assembly::junction_edit as assembly_junction_edit;
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
    pub substitutions: Vec<collagene_core::cloning::assembly::substitution_intent::Substitution>,
    #[serde(deserialize_with = "deserialize_internal_edits")]
    pub edits: Vec<collagene_core::cloning::assembly::edit_intent::Edit>,
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
    collagene_core::cloning::assembly::substitution_intent::Substitution,
    256
);
decoder!(
    deserialize_internal_edits,
    collagene_core::cloning::assembly::edit_intent::Edit,
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

fn chosen<'a>(bins: &'a [Bin], choices: &[String]) -> Result<Vec<&'a Choice>, String> {
    if choices.len() != bins.len() {
        return Err("Choose exactly one source alternative from every bin.".into());
    }
    bins.iter()
        .zip(choices)
        .map(|(bin, id)| {
            bin.choices
                .iter()
                .find(|c| c.id == *id)
                .ok_or_else(|| "A selected alternative does not belong to its ordered bin.".into())
        })
        .collect()
}

pub fn roles_for_active(draft: &Draft) -> Result<Option<Vec<Role>>, String> {
    draft
        .alternatives
        .as_ref()
        .map(|library| library.roles(library.active()?))
        .transpose()
}
pub fn identity(draft: &Draft) -> Option<(String, String, Vec<String>)> {
    let library = draft.alternatives.as_ref()?;
    let row = library.active().ok()?;
    Some((library.id.clone(), row.id.clone(), row.choices.clone()))
}

/// Source choice authority is explicit; an unfinished working source edit may be
/// persisted, but must not silently inherit the previously chosen source ID.
pub fn require_active_choice_match(draft: &Draft) -> Result<(), String> {
    let Some(library) = &draft.alternatives else {
        return Ok(());
    };
    let row = library.active()?;
    let choices = chosen(&library.bins, &row.choices)?;
    if choices.len() != draft.slots.len()
        || draft
            .slots
            .iter()
            .zip(choices)
            .any(|(slot, choice)| !slot.same_source_geometry(&choice.slot))
    {
        return Err("Working source, interval, strand or digest geometry differs from the selected alternatives. Explicitly capture each changed source as a new choice before saving or switching constructs or running preparation.".into());
    }
    if draft.methods.source_roles.as_ref() != Some(&library.roles(row)?) {
        return Err("Working preparation identities differ from the active construct and selected alternatives.".into());
    }
    Ok(())
}

impl Construct {
    fn sync_source_order(&mut self, bins: &[Bin]) -> Result<(), String> {
        let selected = chosen(bins, &self.choices)?;
        if selected.len() != self.edits.methods.rows.len() {
            return Err("Construct preparation rows differ from selected source bins.".into());
        }
        let mut bottom = selected.len();
        for (index, (choice, raw)) in selected
            .iter()
            .zip(&mut self.edits.methods.rows)
            .enumerate()
        {
            raw.document = index;
            raw.anneal_bottom_document = if raw.preparation == OverlapPreparation::AnnealedOligos
                && choice.slot.anneal_bottom.is_some()
            {
                let index = bottom;
                bottom += 1;
                Some(index)
            } else {
                None
            };
        }
        Ok(())
    }
    fn validate_selected_sources(&self, bins: &[Bin]) -> Result<(), String> {
        let selected = chosen(bins, &self.choices)?;
        let mut bottom = selected.len();
        for (index, (choice, raw)) in selected.iter().zip(&self.edits.methods.rows).enumerate() {
            model::validate_slot(&choice.slot, raw, index)?;
            if raw.preparation == OverlapPreparation::AnnealedOligos {
                let expected = choice.slot.anneal_bottom.as_ref().map(|_| {
                    let at = bottom;
                    bottom += 1;
                    at
                });
                if raw.anneal_bottom_document != expected
                    || (raw.anneal_candidate.is_some()
                        && (choice.slot.anneal_context(raw).is_none()
                            || raw.anneal_source_context != choice.slot.anneal_context(raw)))
                {
                    return Err("Construct annealing selection differs from its ordered original source events or options.".into());
                }
            }
        }
        Ok(())
    }
}
impl Library {
    fn entity_id(&self, id: &str, kind: &str) -> bool {
        id.len() <= 48
            && id
                .strip_prefix(&format!("{}-{kind}-", self.id))
                .and_then(|n| n.parse::<u64>().ok())
                .is_some_and(|n| n > 0 && n < self.next_id)
    }
    fn allocate(&mut self, kind: &str) -> Result<String, String> {
        let value = self.next_id;
        self.next_id = value
            .checked_add(1)
            .ok_or("The assembly alternative identifier counter is exhausted.")?;
        Ok(format!("{}-{kind}-{value}", self.id))
    }
    fn active(&self) -> Result<&Construct, String> {
        self.rows
            .iter()
            .find(|r| r.id == self.active_row)
            .ok_or_else(|| "The active construct row is unavailable.".into())
    }
    fn active_mut(&mut self) -> Result<&mut Construct, String> {
        self.rows
            .iter_mut()
            .find(|r| r.id == self.active_row)
            .ok_or_else(|| "The active construct row is unavailable.".into())
    }
    fn roles(&self, row: &Construct) -> Result<Vec<Role>, String> {
        chosen(&self.bins, &row.choices)?;
        self.bins
            .iter()
            .zip(&row.choices)
            .map(|(bin, choice)| {
                let (input_id, role) = collagene_core::planning::combinatorial::source_role_ids(
                    &self.id, &row.id, &bin.id, choice,
                )
                .map_err(|e| e.to_string())?;
                Ok(Role { input_id, role })
            })
            .collect()
    }
    // Matrix display carries only bounded stable IDs. It neither copies nor attests
    // source models; mutation/preparation still uses the complete draft validation.
    fn matrix_counts(&self) -> Result<Vec<usize>, String> {
        if !design_id(&self.id)
            || self.next_id == 0
            || self.bins.is_empty()
            || self.bins.len() > MAX_BINS
        {
            return Err("Invalid bounded alternative matrix.".into());
        }
        let mut ids = BTreeSet::new();
        let mut counts = Vec::new();
        let mut total = 0usize;
        for bin in &self.bins {
            if !self.entity_id(&bin.id, "bin")
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
                if !self.entity_id(&choice.id, "choice") || !ids.insert(&choice.id) {
                    return Err("Invalid alternative matrix choice identity.".into());
                }
            }
            counts.push(bin.choices.len());
        }
        Ok(counts)
    }
    pub fn matrix_size(&self) -> Result<u128, String> {
        collagene_core::planning::combinatorial::cartesian_size(&self.matrix_counts()?)
            .map_err(|e| e.to_string())
    }
    pub fn combinations_page(
        &self,
        offset: u128,
        limit: usize,
    ) -> Result<Vec<Vec<String>>, String> {
        collagene_core::planning::combinatorial::cartesian_page(
            &self.matrix_counts()?,
            offset,
            limit,
        )
        .map_err(|e| e.to_string())
        .map(|rows| {
            rows.into_iter()
                .map(|row| {
                    row.into_iter()
                        .enumerate()
                        .map(|(i, choice)| self.bins[i].choices[choice].id.clone())
                        .collect()
                })
                .collect()
        })
    }

    pub(crate) fn sync_source_order(&mut self) -> Result<(), String> {
        // Choice templates own actual captured events, never construct-local bottom ordinals.
        for bin in &mut self.bins {
            for choice in &mut bin.choices {
                choice.row.anneal_bottom_document = None;
            }
        }
        for row in &mut self.rows {
            row.sync_source_order(&self.bins)?;
        }
        Ok(())
    }
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
        model::bounded_size(self)?;
        let mut ids = BTreeSet::new();
        let mut total = 0usize;
        for (index, bin) in self.bins.iter().enumerate() {
            if !self.entity_id(&bin.id, "bin")
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
                if !self.entity_id(&choice.id, "choice") || !ids.insert(&choice.id) {
                    return Err("Invalid or duplicate assembly choice identity.".into());
                }
                model::validate_slot(&choice.slot, &choice.row, index)?;
                if choice.row.anneal_bottom_document.is_some() {
                    return Err("Source alternatives must not retain construct-local bottom input ordinals.".into());
                }
                if choice.row.anneal_candidate.is_some()
                    && (choice.slot.anneal_context(&choice.row).is_none()
                        || choice.row.anneal_source_context
                            != choice.slot.anneal_context(&choice.row))
                {
                    return Err("Captured alternative annealing selection differs from its ordered source events or options.".into());
                }
            }
        }
        for row in &self.rows {
            if !self.entity_id(&row.id, "row") || !ids.insert(&row.id) {
                return Err("Invalid or duplicate construct row identity.".into());
            }
            model::validate_controls(
                &row.edits.name,
                row.edits.circular,
                row.edits.active_slot,
                self.bins.len(),
                &row.edits.methods,
                &row.edits.junction_edits,
            )?;
            model::validate_internal_declarations(
                &row.edits.methods,
                self.bins.len(),
                &row.edits.substitutions,
                &row.edits.edits,
            )?;
            row.validate_selected_sources(&self.bins)?;
            if row.edits.methods.source_roles.as_ref() != Some(&self.roles(row)?) {
                return Err("Saved construct preparation roles differ from its stable selected alternatives.".into());
            }
        }
        self.active()?;
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
    if let Some(library) = &mut draft.alternatives {
        library.sync_source_order()?;
    }
    if draft.alternatives.is_some() {
        draft.methods.source_roles = roles_for_active(&draft)?;
    }
    draft.sync_rows();
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
        library.active_mut()?.edits = saved;
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
        let bin = library.allocate("bin")?;
        let choice = library.allocate("choice")?;
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
    let row_id = library.allocate("row")?;
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
    let new_roles = library.roles(&row)?;
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
    next.version = crate::assembly::draft::model::CURRENT_VERSION;
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
    candidate: Option<crate::cloning::preparation::golden_gate_pool_form::CandidateIntent>,
) -> Result<Draft, String> {
    let mut next = stage(draft)?;
    require_active_choice_match(&next)?;
    if next.alternatives.is_none()
        || next.methods.mode != crate::cloning::preparation::cloning_methods_model::Mode::GoldenGate
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
    next.version = crate::assembly::draft::model::CURRENT_VERSION;
    finish(next)
}
/// Persist raw working edits even while a source edit awaits explicit choice
/// capture. Existing construct edits are refreshed only when their sources match.
pub fn for_storage(draft: &Draft) -> Result<Draft, String> {
    let mut next = stage(draft)?;
    if require_active_choice_match(&next).is_ok() {
        save_inner(&mut next)?;
    }
    next.version = crate::assembly::draft::model::CURRENT_VERSION;
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
    let id = library.allocate("choice")?;
    library.bins[index].choices.push(Choice {
        id: id.clone(),
        ..choice
    });
    let saved = library.active_mut()?;
    saved.choices[index] = id;
    saved.edits.methods.rows[index] = next.methods.rows[index].clone();
    let expected = library.roles(library.active()?)?;
    library.active_mut()?.edits.methods.source_roles = Some(expected);
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
    let slots = chosen(&library.bins, &row.choices)?
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
    let selected = chosen(&library.bins, &selection)?;
    saved.methods.rows = selected
        .into_iter()
        .map(|choice| choice.row.clone())
        .collect();
    let id = library.allocate("row")?;
    let mut row = Construct {
        id: id.clone(),
        choices: selection,
        edits: saved,
    };
    row.edits.methods.source_roles = Some(library.roles(&row)?);
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
    use crate::assembly::draft::model;
    use crate::assembly::draft::model::Slot;
    use crate::sources::identity::SourceInfo;
    use collagene_core::{documents::workspace::Workspace, material::lineage::EventBinding};
    use collagene_model::{Molecule, Topology};
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
                let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    start: value,
                    ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = 0;
            }
            {
                let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    end: value,
                    ..
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
            selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
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
        second =
            crate::assembly::edit::assembly_edit_model::declare(&second, 0, "3", "5", "TAG", None)
                .unwrap();
        second =
            crate::assembly::edit::assembly_substitution_model::declare(&second, 1, "1", "A", None)
                .unwrap();
        let target = |d: &Draft| {
            let (inputs, base, _) = model::inputs_request(d, None).unwrap();
            collagene_core::cloning::assembly::edit_intent::build(
                &inputs,
                &model::edit_request(d, base),
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
        let restored = model::decode(&model::encode(&first_again).unwrap()).unwrap();
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
        let reorder = crate::assembly::library::assembly_library_topology::propose(
            &second_again,
            crate::assembly::library::assembly_library_topology::Change::Reorder { bin_ids: bins },
        )
        .unwrap();
        assert!(
            reorder
                .impacts()
                .iter()
                .all(|i| i.recovered_internal_edits == 0)
        );
        let reordered = crate::assembly::library::assembly_library_topology::apply(
            &second_again,
            &reorder,
            reorder.revision(),
        )
        .unwrap();
        assert_eq!(target(&reordered.draft), "ATGCAACCAATAGGGTT");
        let removed_id = reordered.draft.alternatives.as_ref().unwrap().bins[1]
            .id
            .clone();
        let removal = crate::assembly::library::assembly_library_topology::propose(
            &reordered.draft,
            crate::assembly::library::assembly_library_topology::Change::Remove {
                bin_id: removed_id,
            },
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
        let applied = crate::assembly::library::assembly_library_topology::apply(
            &reordered.draft,
            &removal,
            removal.revision(),
        )
        .unwrap();
        assert_eq!(applied.recovery.edits, second.edits);
        assert!(applied.draft.edits.is_empty());
        assert_eq!(applied.draft.substitutions.len(), 1);
    }
    #[test]
    fn source_replacement_and_new_rows_never_rebind_internal_literals() {
        let first = enable(&ordinary(), "new-source").unwrap();
        let first =
            crate::assembly::edit::assembly_edit_model::declare(&first, 1, "3", "5", "TAG", None)
                .unwrap();
        let choices = first.alternatives.as_ref().unwrap().rows[0].choices.clone();
        let second = add_row(&first, choices).unwrap();
        assert!(second.edits.is_empty() && second.substitutions.is_empty());
        let first_id = first.alternatives.as_ref().unwrap().active_row.clone();
        let mut working = open_row(&second, &first_id).unwrap();
        replace_insert(&mut working);
        assert!(
            crate::assembly::edit::assembly_edit_model::declare(&working, 1, "1", "2", "A", None)
                .is_err()
        );
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
        assert!(model::decode(&serde_json::to_vec(&missing).unwrap()).is_err());
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
            second.alternatives.as_ref().unwrap().matrix_size().unwrap(),
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
        let restored = model::decode(&model::encode(&second).unwrap()).unwrap();
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
        let bytes = model::encode(&working).unwrap();
        let restored = model::decode(&bytes).unwrap();
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
        assert!(model::encode(&bad).is_err());
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
        assert_eq!(library.matrix_size().unwrap(), 2);
        let first = library.combinations_page(0, 1).unwrap();
        let second = library.combinations_page(1, 1).unwrap();
        assert_eq!(first[0][1], initial_id);
        assert_ne!(first, second);
        assert!(library.combinations_page(2, 1).unwrap().is_empty());
        let old_counter = library.next_id;
        let removed = remove_choice(&draft, &initial_id).unwrap();
        assert_eq!(removed.alternatives.as_ref().unwrap().next_id, old_counter);
        assert!(require_active_choice_match(&removed).is_ok());
        assert!(removed.slots == draft.slots);
    }
}
