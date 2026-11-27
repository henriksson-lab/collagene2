//! Reviewed transactional edits to ordered source-library topology.
//! Install `Applied.recovery` durably with the new draft, never discard it.
use crate::{
    assembly_alternatives::{self as alternatives, Bin, Choice, Library},
    assembly_draft_model::{Draft, Slot},
    assembly_role_ids::Role,
    cloning_methods_model::{Mode, OverlapInput},
};
use serde::Serialize;
use std::{collections::BTreeSet, io, sync::Arc};

#[derive(Clone, Serialize)]
pub struct ChoiceSeed {
    pub slot: Slot,
    pub preparation: OverlapInput,
}
#[derive(Clone, Serialize)]
pub struct AddedRow {
    pub construct_id: String,
    /// Explicit index into Add.choices, never an implicit first choice.
    pub choice: usize,
    pub preparation: OverlapInput,
    /// 0 prepared, 1 whole source, 2 whole amplicon; only used by an existing pool.
    pub pool_substrate: u8,
}
#[derive(Clone, Serialize)]
pub enum Change {
    Add {
        position: usize,
        label: String,
        choices: Vec<ChoiceSeed>,
        rows: Vec<AddedRow>,
    },
    Remove {
        bin_id: String,
    },
    Reorder {
        bin_ids: Vec<String>,
    },
}
#[derive(Clone, PartialEq, Eq, Serialize)]
pub struct RowImpact {
    pub construct_id: String,
    pub name: String,
    pub construct_removed: bool,
    pub construct_added: bool,
    pub primer_review_bins: Vec<String>,
    pub removed_preparation_bins: Vec<String>,
    pub retained_internal_edits: usize,
    pub recovered_internal_edits: usize,
    pub retained_junctions: usize,
    pub recovered_junctions: usize,
    pub new_junctions: Vec<usize>,
    pub overlap_text_changed: bool,
    pub role_order_changed: bool,
    pub candidate_cleared: bool,
}
/// Unforgeable in-process proposal: callers can inspect but cannot replace the
/// validated before/after or impact rows after receiving its review token.
pub struct Proposal {
    before: Arc<Draft>,
    after: Arc<Draft>,
    impacts: Vec<RowImpact>,
    revision: String,
}
impl Proposal {
    pub fn before(&self) -> &Draft {
        &self.before
    }
    pub fn after(&self) -> &Draft {
        &self.after
    }
    pub fn impacts(&self) -> &[RowImpact] {
        &self.impacts
    }
    pub fn revision(&self) -> &str {
        &self.revision
    }
}
pub struct Applied {
    pub draft: Draft,
    /// Mandatory recovery contains exact working raw text, removed source bins,
    /// old primers, rejected overlap text, source pins and junction declarations.
    pub recovery: Arc<Draft>,
    pub impacts: Vec<RowImpact>,
    pub proposal_revision: String,
}
struct BoundedBytes(usize);
impl io::Write for BoundedBytes {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|size| *size <= 32 << 20)
            .ok_or_else(|| io::Error::other("topology input budget"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
struct Hash(blake3::Hasher);
impl io::Write for Hash {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn edges(ids: &[String], circular: bool) -> Vec<(String, String)> {
    let mut result = ids
        .windows(2)
        .map(|pair| (pair[0].clone(), pair[1].clone()))
        .collect::<Vec<_>>();
    if circular && !ids.is_empty() {
        result.push((ids[ids.len() - 1].clone(), ids[0].clone()));
    }
    result
}
fn neighbors<'a>(
    ids: &'a [String],
    id: &str,
    circular: bool,
) -> (Option<&'a str>, Option<&'a str>) {
    let Some(position) = ids.iter().position(|value| value == id) else {
        return (None, None);
    };
    (
        position
            .checked_sub(1)
            .map(|index| ids[index].as_str())
            .or_else(|| circular.then(|| ids[ids.len() - 1].as_str())),
        ids.get(position + 1)
            .map(String::as_str)
            .or_else(|| circular.then(|| ids[0].as_str())),
    )
}
fn allocate(library: &mut Library, kind: &str) -> Result<String, String> {
    let value = library.next_id;
    library.next_id = value
        .checked_add(1)
        .ok_or("The library identifier counter is exhausted.")?;
    let id = format!("{}-{kind}-{value}", library.id);
    if id.len() > 48 {
        return Err("Library identifier exceeds its bound.".into());
    }
    Ok(id)
}
fn roles(library: &Library, row_id: &str, choices: &[String]) -> Result<Vec<Role>, String> {
    library
        .bins
        .iter()
        .zip(choices)
        .map(|(bin, choice)| {
            let (input_id, role) = collagene_core::combinatorial_design::source_role_ids(
                &library.id,
                row_id,
                &bin.id,
                choice,
            )
            .map_err(|error| error.to_string())?;
            Ok(Role { input_id, role })
        })
        .collect()
}
fn install_active(draft: &mut Draft) -> Result<(), String> {
    let library = draft
        .alternatives
        .as_ref()
        .ok_or("A source library is required.")?;
    let row = library
        .rows
        .iter()
        .find(|row| row.id == library.active_row)
        .ok_or("Active row is missing.")?;
    let mut slots = Vec::with_capacity(library.bins.len());
    for (bin, choice_id) in library.bins.iter().zip(&row.choices) {
        let mut slot = bin
            .choices
            .iter()
            .find(|choice| &choice.id == choice_id)
            .ok_or("Selected source is missing.")?
            .slot
            .clone();
        slot.label = bin.label.clone();
        slots.push(slot);
    }
    draft.slots = slots;
    draft.name = row.edits.name.clone();
    draft.circular = row.edits.circular;
    draft.methods = row.edits.methods.clone();
    draft.junction_edits = row.edits.junction_edits.clone();
    draft.substitutions = row.edits.substitutions.clone();
    draft.edits = row.edits.edits.clone();
    draft.active_slot = row.edits.active_slot;
    Ok(())
}

pub fn propose(current: &Draft, change: Change) -> Result<Proposal, String> {
    current.validate()?;
    // Reject cheap shape violations before walking any caller-supplied source
    // molecules or raw fields for the serialized input budget.
    match &change {
        Change::Add {
            choices,
            rows,
            label,
            ..
        } if choices.is_empty()
            || choices.len() > alternatives::MAX_CHOICES_PER_BIN
            || rows.len() > alternatives::MAX_ROWS
            || label.len() > 80 =>
        {
            return Err("Topology addition exceeds source, construct or label bounds.".into());
        }
        Change::Reorder { bin_ids }
            if bin_ids.len() > alternatives::MAX_BINS || bin_ids.iter().any(|id| id.len() > 48) =>
        {
            return Err("Topology ordering exceeds stable bin identifier bounds.".into());
        }
        Change::Remove { bin_id } if bin_id.len() > 48 => {
            return Err("Topology removal exceeds stable bin identifier bounds.".into());
        }
        _ => {}
    }
    let mut counter = BoundedBytes(0);
    serde_json::to_writer(&mut counter, &change)
        .map_err(|_| "Topology edit exceeds its 32 MiB input budget.")?;
    // No overwrite of pending source geometry or unsaved active raw settings.
    alternatives::require_active_choice_match(current)?;
    let mut after = alternatives::save_active(current)?;
    let library = after
        .alternatives
        .as_mut()
        .ok_or("Enable a source library before editing its bins.")?;
    let old_ids = library
        .bins
        .iter()
        .map(|bin| bin.id.clone())
        .collect::<Vec<_>>();
    let old_rows = library.rows.clone();
    let old_bins = library.bins.clone();
    let mut addition: Option<(usize, Vec<AddedRow>)> = None;
    match change {
        Change::Add {
            position,
            label,
            choices,
            rows,
        } => {
            if library.bins.len() >= alternatives::MAX_BINS || position > library.bins.len() {
                return Err("New bin position or bin count is out of bounds.".into());
            }
            if label.trim().is_empty() || label.len() > 80 || label.chars().any(char::is_control) {
                return Err("Enter a bounded bin label.".into());
            }
            if choices.is_empty()
                || choices.len() > alternatives::MAX_CHOICES_PER_BIN
                || library
                    .bins
                    .iter()
                    .map(|bin| bin.choices.len())
                    .sum::<usize>()
                    + choices.len()
                    > alternatives::MAX_CHOICES
            {
                return Err("New choices exceed the library limits.".into());
            }
            if rows.len() != library.rows.len() {
                return Err(
                    "Every construct requires an explicit new-bin choice and preparation.".into(),
                );
            }
            let mut selected = BTreeSet::new();
            for selection in &rows {
                if !selected.insert(selection.construct_id.as_str())
                    || !library
                        .rows
                        .iter()
                        .any(|row| row.id == selection.construct_id)
                    || selection.choice >= choices.len()
                    || selection.pool_substrate > 2
                {
                    return Err(
                        "New-bin construct choices are missing, duplicated or invalid.".into(),
                    );
                }
            }
            let id = allocate(library, "bin")?;
            let mut values = Vec::with_capacity(choices.len());
            for seed in choices {
                values.push(Choice {
                    id: allocate(library, "choice")?,
                    slot: seed.slot,
                    row: seed.preparation,
                });
            }
            library.bins.insert(
                position,
                Bin {
                    id,
                    label,
                    choices: values,
                },
            );
            addition = Some((position, rows));
        }
        Change::Remove { bin_id } => {
            if library.bins.len() <= 1 {
                return Err("Keep at least one source bin.".into());
            }
            let position = library
                .bins
                .iter()
                .position(|bin| bin.id == bin_id)
                .ok_or("The source bin is unavailable.")?;
            library.bins.remove(position);
        }
        Change::Reorder { bin_ids } => {
            let unique = bin_ids.iter().collect::<BTreeSet<_>>();
            if bin_ids.len() != old_ids.len()
                || unique.len() != old_ids.len()
                || bin_ids.iter().any(|id| !old_ids.contains(id))
            {
                return Err("Reordering requires each existing stable bin ID exactly once.".into());
            }
            if bin_ids == old_ids {
                return Err("The source bin order is unchanged.".into());
            }
            library.bins = bin_ids
                .iter()
                .map(|id| {
                    old_bins
                        .iter()
                        .find(|bin| &bin.id == id)
                        .expect("validated permutation")
                        .clone()
                })
                .collect();
        }
    }
    let new_ids = library
        .bins
        .iter()
        .map(|bin| bin.id.clone())
        .collect::<Vec<_>>();
    for (index, bin) in library.bins.iter_mut().enumerate() {
        for choice in &mut bin.choices {
            choice.row.document = index;
        }
    }
    let mut impacts = Vec::with_capacity(old_rows.len());
    for old in old_rows {
        let mut row = old.clone();
        let old_edges = edges(&old_ids, old.edits.circular);
        let new_edges = edges(&new_ids, old.edits.circular);
        row.choices.clear();
        row.edits.methods.rows.clear();
        let old_pool = old.edits.methods.golden_gate_pool.as_ref();
        let mut pool_choices = Vec::with_capacity(new_ids.len());
        for (index, bin) in library.bins.iter().enumerate() {
            if let Some(previous) = old_ids.iter().position(|id| id == &bin.id) {
                row.choices.push(old.choices[previous].clone());
                let mut raw = old.edits.methods.rows[previous].clone();
                raw.document = index;
                row.edits.methods.rows.push(raw);
                pool_choices.push(
                    old_pool
                        .and_then(|pool| pool.choices.get(previous))
                        .copied()
                        .unwrap_or(0),
                );
            } else {
                let (_, selections) = addition
                    .as_ref()
                    .ok_or("Missing explicit new-bin choices.")?;
                let selection = selections
                    .iter()
                    .find(|selection| selection.construct_id == row.id)
                    .ok_or("Missing construct's new-bin choice.")?;
                row.choices.push(bin.choices[selection.choice].id.clone());
                let mut raw = selection.preparation.clone();
                raw.document = index;
                row.edits.methods.rows.push(raw);
                pool_choices.push(selection.pool_substrate);
            }
        }
        let mut retained = Vec::new();
        for edit in &old.edits.junction_edits {
            if let Some(position) = old_edges
                .get(edit.junction)
                .and_then(|edge| new_edges.iter().position(|value| value == edge))
            {
                let mut value = edit.clone();
                value.junction = position;
                retained.push(value);
            }
        }
        retained.sort_by_key(|edit| edit.junction);
        row.edits.junction_edits = retained;
        let new_junctions = new_edges
            .iter()
            .enumerate()
            .filter_map(|(index, edge)| (!old_edges.contains(edge)).then_some(index))
            .collect::<Vec<_>>();
        // No inference from raw text that cannot be associated with old edges.
        // New edges deliberately leave the complete overlap profile uncomputed.
        if row.edits.methods.mode == Mode::Overlap {
            let values = old
                .edits
                .methods
                .overlaps
                .split(',')
                .map(str::trim)
                .collect::<Vec<_>>();
            row.edits.methods.overlaps =
                if new_junctions.is_empty() && values.len() == old_edges.len() {
                    new_edges
                        .iter()
                        .map(|edge| {
                            values[old_edges
                                .iter()
                                .position(|old| old == edge)
                                .expect("retained edge")]
                        })
                        .collect::<Vec<_>>()
                        .join(",")
                } else {
                    String::new()
                };
        }
        let candidate_cleared = old_pool.is_some_and(|pool| pool.candidate.is_some());
        if let Some(pool) = &mut row.edits.methods.golden_gate_pool {
            if !pool.reaction && pool_choices.iter().any(|choice| *choice != 0) {
                return Err(
                    "New whole-original/amplicon choices require this construct's reaction mode."
                        .into(),
                );
            }
            pool.choices = pool_choices;
            pool.candidate = None;
        }
        let previous_active = old_ids.get(old.edits.active_slot);
        row.edits.active_slot = previous_active
            .and_then(|id| new_ids.iter().position(|value| value == id))
            .unwrap_or(old.edits.active_slot.min(new_ids.len() - 1));
        let new_roles = roles(library, &row.id, &row.choices)?;
        row.edits
            .substitutions
            .retain(|edit| new_roles.iter().any(|r| r.role == edit.role));
        row.edits
            .edits
            .retain(|edit| new_roles.iter().any(|r| r.role == edit.role));
        let internal_retained = row.edits.substitutions.len() + row.edits.edits.len();
        row.edits.methods.source_roles = Some(new_roles);
        impacts.push(RowImpact {
            construct_id: row.id.clone(),
            name: row.edits.name.clone(),
            construct_removed: false,
            construct_added: false,
            primer_review_bins: new_ids
                .iter()
                .filter(|id| {
                    !old_ids.contains(id)
                        || neighbors(&old_ids, id, old.edits.circular)
                            != neighbors(&new_ids, id, old.edits.circular)
                })
                .cloned()
                .collect(),
            removed_preparation_bins: old_ids
                .iter()
                .filter(|id| !new_ids.contains(id))
                .cloned()
                .collect(),
            retained_internal_edits: internal_retained,
            recovered_internal_edits: old.edits.substitutions.len() + old.edits.edits.len()
                - internal_retained,
            retained_junctions: row.edits.junction_edits.len(),
            recovered_junctions: old.edits.junction_edits.len() - row.edits.junction_edits.len(),
            new_junctions,
            overlap_text_changed: old.edits.methods.overlaps != row.edits.methods.overlaps,
            role_order_changed: old.edits.methods.source_roles != row.edits.methods.source_roles,
            candidate_cleared,
        });
        let position = library
            .rows
            .iter()
            .position(|value| value.id == row.id)
            .expect("existing row");
        library.rows[position] = row;
    }
    install_active(&mut after)?;
    after.version = crate::assembly_draft_model::CURRENT_VERSION;
    after.validate()?;
    let mut writer = Hash(blake3::Hasher::new());
    serde_json::to_writer(
        &mut writer,
        &("assembly-library-topology-1", current, &after, &impacts),
    )
    .map_err(|error| error.to_string())?;
    Ok(Proposal {
        before: Arc::new(current.clone()),
        after: Arc::new(after),
        impacts,
        revision: writer.0.finalize().to_hex().to_string(),
    })
}

#[derive(Clone, Copy)]
pub struct EffectiveEdits<'a> {
    pub name: &'a str,
    pub circular: bool,
    pub methods: &'a crate::cloning_methods_model::MethodsState,
    pub junctions: &'a [collagene_core::assembly_junction_edit::Edit],
    pub substitutions: &'a [collagene_core::assembly_substitution_intent::Substitution],
    pub edits: &'a [collagene_core::assembly_edit_intent::Edit],
}
pub fn effective_edits<'a>(
    draft: &'a Draft,
    row: &'a alternatives::Construct,
) -> EffectiveEdits<'a> {
    if draft
        .alternatives
        .as_ref()
        .is_some_and(|library| library.active_row == row.id)
    {
        EffectiveEdits {
            name: &draft.name,
            circular: draft.circular,
            methods: &draft.methods,
            junctions: &draft.junction_edits,
            substitutions: &draft.substitutions,
            edits: &draft.edits,
        }
    } else {
        EffectiveEdits {
            name: &row.edits.name,
            circular: row.edits.circular,
            methods: &row.edits.methods,
            junctions: &row.edits.junction_edits,
            substitutions: &row.edits.substitutions,
            edits: &row.edits.edits,
        }
    }
}

/// Build explicit new-bin assignments against the same effective raw row
/// settings that propose() saves before transforming the library.
pub fn addition_rows(
    draft: &Draft,
    choices: &[ChoiceSeed],
    assigned: &[Option<usize>],
    substrates: &[Option<u8>],
) -> Result<Vec<AddedRow>, String> {
    draft.validate()?;
    let library = draft
        .alternatives
        .as_ref()
        .ok_or("Enable a source library first.")?;
    if choices.is_empty()
        || choices.len() > alternatives::MAX_CHOICES_PER_BIN
        || assigned.len() != library.rows.len()
        || substrates.len() != library.rows.len()
    {
        return Err("Every construct requires a bounded explicit new-bin assignment.".into());
    }
    let mut declarations = Vec::with_capacity(library.rows.len());
    let mut budget = BoundedBytes(0);
    for (index, row) in library.rows.iter().enumerate() {
        let effective = effective_edits(draft, row);
        let choice = assigned[index].ok_or_else(|| {
            format!(
                "Choose a source and copied preparation for {}.",
                effective.name
            )
        })?;
        let seed = choices
            .get(choice)
            .ok_or("An assigned choice was removed.")?;
        let substrate = if effective.methods.golden_gate_pool.is_some() {
            substrates[index]
                .ok_or_else(|| format!("Choose a reaction substrate for {}.", effective.name))?
        } else {
            0
        };
        if substrate > 2 {
            return Err("Unknown reaction substrate choice.".into());
        }
        crate::assembly_draft_model::validate_preparation_input(&seed.preparation)?;
        serde_json::to_writer(
            &mut budget,
            &(&row.id, choice, &seed.preparation, substrate),
        )
        .map_err(|_| "New-bin construct preparations exceed32MiB.")?;
        declarations.push((row, choice, seed, substrate));
    }
    Ok(declarations
        .into_iter()
        .map(|(row, choice, seed, pool_substrate)| AddedRow {
            construct_id: row.id.clone(),
            choice,
            preparation: seed.preparation.clone(),
            pool_substrate,
        })
        .collect())
}

/// Propose an exact raw recovery frame. The host independently attests retained
/// source snapshots and supplies the allocation high-water mark from its native
/// recovery history. Neither argument is scientific evidence or an ID alias.
pub fn propose_restore(
    current: &Draft,
    recovered: &Draft,
    allocation_highwater: u64,
) -> Result<Proposal, String> {
    current.validate()?;
    recovered.validate()?;
    let previous = current
        .alternatives
        .as_ref()
        .ok_or("The current source library is unavailable.")?;
    let restored = recovered
        .alternatives
        .as_ref()
        .ok_or("The recovered source library is unavailable.")?;
    if previous.id != restored.id {
        return Err("A topology recovery must belong to the same exact source library.".into());
    }
    let mut after = recovered.clone();
    after.version = crate::assembly_draft_model::CURRENT_VERSION;
    after
        .alternatives
        .as_mut()
        .expect("checked recovery library")
        .next_id = previous
        .next_id
        .max(restored.next_id)
        .max(allocation_highwater);
    after.validate()?;
    if current == &after {
        return Err("The recovered source library is unchanged.".into());
    }
    let old_bins = previous
        .bins
        .iter()
        .map(|bin| bin.id.clone())
        .collect::<Vec<_>>();
    let new_bins = restored
        .bins
        .iter()
        .map(|bin| bin.id.clone())
        .collect::<Vec<_>>();
    let union = previous
        .rows
        .iter()
        .map(|row| row.id.as_str())
        .chain(restored.rows.iter().map(|row| row.id.as_str()))
        .collect::<BTreeSet<_>>();
    if union.len() > alternatives::MAX_ROWS * 2 {
        return Err("Recovery impact rows exceed the bounded review limit.".into());
    }
    let mut impacts = Vec::with_capacity(union.len());
    for id in union {
        let old = previous
            .rows
            .iter()
            .find(|row| row.id == id)
            .map(|row| effective_edits(current, row));
        let next = restored
            .rows
            .iter()
            .find(|row| row.id == id)
            .map(|row| effective_edits(recovered, row));
        let old_edges = old
            .map(|row| edges(&old_bins, row.circular))
            .unwrap_or_default();
        let new_edges = next
            .map(|row| edges(&new_bins, row.circular))
            .unwrap_or_default();
        let retained = old
            .map(|row| {
                row.junctions
                    .iter()
                    .filter(|edit| {
                        let Some(edge) = old_edges.get(edit.junction) else {
                            return false;
                        };
                        next.is_some_and(|row| {
                            row.junctions.iter().any(|candidate| {
                                new_edges.get(candidate.junction) == Some(edge)
                                    && candidate.remove_left == edit.remove_left
                                    && candidate.remove_right == edit.remove_right
                                    && candidate.inserted == edit.inserted
                            })
                        })
                    })
                    .count()
            })
            .unwrap_or(0);
        let old_candidate = old
            .and_then(|row| row.methods.golden_gate_pool.as_ref())
            .and_then(|pool| pool.candidate.as_ref());
        let new_candidate = next
            .and_then(|row| row.methods.golden_gate_pool.as_ref())
            .and_then(|pool| pool.candidate.as_ref());
        let internal_retained = old
            .zip(next)
            .map(|(old, next)| {
                old.substitutions
                    .iter()
                    .filter(|e| next.substitutions.contains(e))
                    .count()
                    + old.edits.iter().filter(|e| next.edits.contains(e)).count()
            })
            .unwrap_or(0);
        impacts.push(RowImpact {
            construct_id: id.into(),
            name: next.or(old).expect("union contains a row").name.to_string(),
            construct_removed: next.is_none(),
            construct_added: old.is_none(),
            primer_review_bins: if next.is_some() {
                new_bins.clone()
            } else {
                old_bins.clone()
            },
            removed_preparation_bins: if next.is_none() {
                old_bins.clone()
            } else {
                old_bins
                    .iter()
                    .filter(|id| !new_bins.contains(id))
                    .cloned()
                    .collect()
            },
            retained_internal_edits: internal_retained,
            recovered_internal_edits: old
                .map(|row| row.substitutions.len() + row.edits.len())
                .unwrap_or(0)
                - internal_retained,
            retained_junctions: retained,
            recovered_junctions: old.map(|row| row.junctions.len()).unwrap_or(0) - retained,
            new_junctions: new_edges
                .iter()
                .enumerate()
                .filter_map(|(index, edge)| (!old_edges.contains(edge)).then_some(index))
                .collect(),
            overlap_text_changed: old.map(|row| &row.methods.overlaps)
                != next.map(|row| &row.methods.overlaps),
            role_order_changed: old.map(|row| &row.methods.source_roles)
                != next.map(|row| &row.methods.source_roles),
            candidate_cleared: old_candidate.is_some() && old_candidate != new_candidate,
        });
    }
    let mut writer = Hash(blake3::Hasher::new());
    serde_json::to_writer(
        &mut writer,
        &(
            "assembly-library-topology-recovery-1",
            current,
            &after,
            &impacts,
        ),
    )
    .map_err(|error| error.to_string())?;
    Ok(Proposal {
        before: Arc::new(current.clone()),
        after: Arc::new(after),
        impacts,
        revision: writer.0.finalize().to_hex().to_string(),
    })
}
/// Explicit review covers every impact row and returns the mandatory recovery
/// frame. No source mutation or scientific result is installed by this function.
pub fn apply(
    current: &Draft,
    proposal: &Proposal,
    acknowledged_revision: &str,
) -> Result<Applied, String> {
    current.validate()?;
    if current != proposal.before.as_ref() || acknowledged_revision != proposal.revision {
        return Err(
            "Topology proposal is stale or its exact impacts were not acknowledged.".into(),
        );
    }
    Ok(Applied {
        draft: proposal.after.as_ref().clone(),
        recovery: proposal.before.clone(),
        impacts: proposal.impacts.clone(),
        proposal_revision: proposal.revision.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        assembly_draft_model::CapturedSource,
        cloning_methods_model::{Mode, OverlapPreparation},
        golden_gate_pool_form::{CandidateIntent, PoolForm},
        source_identity::SourceInfo,
    };
    use collagene_core::{Molecule, Topology, run_lineage::EventBinding, workspace::Workspace};
    fn fixture() -> Draft {
        let mut draft = Draft::default();
        draft.circular = false;
        draft.methods.mode = Mode::Overlap;
        draft.slots.push(Slot::empty("Third"));
        draft.sync_rows();
        draft.methods.overlaps = "4,5".into();
        for (index, slot) in draft.slots.iter_mut().enumerate() {
            let molecule =
                Molecule::dna(format!("Source {index}"), "AACCGGTT", Topology::Linear).unwrap();
            let workspace =
                Workspace::from_documents("fixture", vec![molecule.clone()], 0).unwrap();
            let doc = &workspace.documents[0];
            {
                let collagene_core::assembly_intent::Selector::ReferenceRange {
                    end: value, ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = molecule.sequence.len();
            }
            slot.captured = Some(Arc::new(CapturedSource {
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
            }));
        }
        draft.sync_rows();
        for row in &mut draft.methods.rows {
            row.preparation = OverlapPreparation::WholeLinear;
            row.forward_tail = "AAA".into();
        }
        draft.junction_edits = vec![collagene_core::assembly_junction_edit::Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "AC".into(),
        }];
        let first = alternatives::enable(&draft, "topology").unwrap();
        let choices = first.alternatives.as_ref().unwrap().rows[0].choices.clone();
        let mut second = alternatives::add_row(&first, choices).unwrap();
        second.name = "Golden Gate construct".into();
        second.methods.mode = Mode::GoldenGate;
        second.methods.rows[0].forward_tail = "CCC".into();
        second.methods.golden_gate_pool = Some(PoolForm {
            choices: vec![0; 3],
            circular: false,
            candidate: Some(CandidateIntent {
                source_revision: "d".repeat(64),
                context_revision: "a".repeat(64),
                candidate_index: 0,
                directed_sequence_key: "b".repeat(64),
                physical_sequence_key: "c".repeat(64),
            }),
            ..Default::default()
        });
        alternatives::save_active(&second).unwrap()
    }
    fn ids(draft: &Draft) -> Vec<String> {
        draft
            .alternatives
            .as_ref()
            .unwrap()
            .bins
            .iter()
            .map(|bin| bin.id.clone())
            .collect()
    }
    #[test]
    fn reorder_preserves_ids_sources_and_manual_primers_but_proposes_edge_invalidations() {
        let draft = fixture();
        let old = ids(&draft);
        let order = vec![old[1].clone(), old[0].clone(), old[2].clone()];
        let proposal = propose(
            &draft,
            Change::Reorder {
                bin_ids: order.clone(),
            },
        )
        .unwrap();
        assert_eq!(ids(proposal.before()), old);
        assert_eq!(ids(proposal.after()), order);
        let before = draft.alternatives.as_ref().unwrap();
        let after = proposal.after().alternatives.as_ref().unwrap();
        for (previous, next) in before.rows.iter().zip(&after.rows) {
            assert_eq!(previous.id, next.id);
            assert_eq!(next.choices[1], previous.choices[0]);
            assert_eq!(
                next.edits.methods.rows[1].forward_tail,
                previous.edits.methods.rows[0].forward_tail
            );
            assert_eq!(
                next.edits.methods.source_roles.as_ref().unwrap()[1],
                previous.edits.methods.source_roles.as_ref().unwrap()[0]
            );
        }
        assert!(after.rows[0].edits.methods.overlaps.is_empty());
        assert!(
            after.rows[1]
                .edits
                .methods
                .golden_gate_pool
                .as_ref()
                .unwrap()
                .candidate
                .is_none()
        );
        assert!(
            proposal
                .impacts()
                .iter()
                .all(|impact| impact.role_order_changed && !impact.primer_review_bins.is_empty())
        );
        let applied = apply(&draft, &proposal, proposal.revision()).unwrap();
        assert!(applied.recovery.as_ref() == &draft);
        assert!(
            applied.recovery.alternatives.as_ref().unwrap().rows[0]
                .edits
                .junction_edits
                .len()
                == 1
        );
    }
    #[test]
    fn remove_requires_exact_current_frame_and_keeps_removed_sources_in_recovery() {
        let draft = fixture();
        let removed = ids(&draft)[1].clone();
        let proposal = propose(
            &draft,
            Change::Remove {
                bin_id: removed.clone(),
            },
        )
        .unwrap();
        assert!(
            proposal
                .impacts()
                .iter()
                .all(|impact| impact.removed_preparation_bins == vec![removed.clone()])
        );
        assert!(apply(&draft, &proposal, "unreviewed").is_err());
        let mut stale = draft.clone();
        stale.methods.rows[0].forward_tail = "Uncommitted manual tail".into();
        assert!(apply(&stale, &proposal, proposal.revision()).is_err());
        let applied = apply(&draft, &proposal, proposal.revision()).unwrap();
        assert_eq!(applied.draft.slots.len(), 2);
        assert_eq!(applied.recovery.slots.len(), 3);
        assert!(
            applied
                .recovery
                .alternatives
                .as_ref()
                .unwrap()
                .bins
                .iter()
                .any(|bin| bin.id == removed)
        );
        assert!(
            applied.draft.alternatives.as_ref().unwrap().next_id
                == draft.alternatives.as_ref().unwrap().next_id
        );
    }
    #[test]
    fn add_requires_all_construct_choices_and_allocates_only_fresh_ids() {
        let draft = fixture();
        let before_next = draft.alternatives.as_ref().unwrap().next_id;
        let seed = ChoiceSeed {
            slot: Slot::empty("New source"),
            preparation: OverlapInput::default(),
        };
        let rows = draft
            .alternatives
            .as_ref()
            .unwrap()
            .rows
            .iter()
            .map(|row| AddedRow {
                construct_id: row.id.clone(),
                choice: 0,
                preparation: OverlapInput::default(),
                pool_substrate: 0,
            })
            .collect::<Vec<_>>();
        assert!(
            propose(
                &draft,
                Change::Add {
                    position: 1,
                    label: "New bin".into(),
                    choices: vec![seed.clone()],
                    rows: rows[..1].to_vec()
                }
            )
            .is_err()
        );
        let proposal = propose(
            &draft,
            Change::Add {
                position: 1,
                label: "New bin".into(),
                choices: vec![seed],
                rows,
            },
        )
        .unwrap();
        let library = proposal.after().alternatives.as_ref().unwrap();
        assert_eq!(library.next_id, before_next + 2);
        assert_eq!(
            library.bins[1].id,
            format!("{}-bin-{before_next}", library.id)
        );
        assert!(
            library
                .rows
                .iter()
                .all(|row| row.choices[1] == library.bins[1].choices[0].id)
        );
        assert!(
            proposal
                .impacts()
                .iter()
                .all(|impact| !impact.new_junctions.is_empty())
        );
        let original_ids = ids(&draft);
        assert_eq!(library.bins[2].id, original_ids[1]);
    }
    #[test]
    fn exact_order_noop_rejects_without_clearing_candidate_or_normalizing_raw_text() {
        let mut draft = fixture();
        draft.alternatives.as_mut().unwrap().rows[0]
            .edits
            .methods
            .overlaps = " 004 , 5 ".into();
        draft.validate().unwrap();
        let original = draft.clone();
        let result = propose(
            &draft,
            Change::Reorder {
                bin_ids: ids(&draft),
            },
        );
        let error = match result {
            Ok(_) => panic!("identical bin order must not create a topology proposal"),
            Err(error) => error,
        };
        assert!(error.contains("unchanged"));
        assert!(draft == original);
        assert_eq!(
            draft.alternatives.as_ref().unwrap().rows[0]
                .edits
                .methods
                .overlaps,
            " 004 , 5 "
        );
        assert!(
            draft
                .methods
                .golden_gate_pool
                .as_ref()
                .unwrap()
                .candidate
                .is_some()
        );
    }

    #[test]
    fn invalid_permutations_bounds_and_pending_source_drafts_do_not_mutate_inputs() {
        let draft = fixture();
        let original = draft.clone();
        let order = ids(&draft);
        assert!(
            propose(
                &draft,
                Change::Reorder {
                    bin_ids: vec![order[0].clone(); 3]
                }
            )
            .is_err()
        );
        assert!(
            propose(
                &draft,
                Change::Add {
                    position: 99,
                    label: "bad".into(),
                    choices: vec![],
                    rows: vec![]
                }
            )
            .is_err()
        );
        let mut pending = draft.clone();
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { start: value, .. } =
                &mut pending.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 1;
        }
        pending.sync_rows();
        assert!(
            propose(
                &pending,
                Change::Remove {
                    bin_id: order[1].clone()
                }
            )
            .is_err()
        );
        assert!(draft == original);
        assert_eq!(pending.slots[0].selector.reference_interval().unwrap().0, 1);
    }

    #[test]
    fn independent_bin_choice_and_counter_limits_fail_without_allocating_ids() {
        let draft = fixture();
        let seed = ChoiceSeed {
            slot: Slot::empty("New source"),
            preparation: OverlapInput::default(),
        };
        let original = draft.clone();
        assert!(
            propose(
                &draft,
                Change::Add {
                    position: 0,
                    label: "Excess choices".into(),
                    choices: vec![seed.clone(); alternatives::MAX_CHOICES_PER_BIN + 1],
                    rows: Vec::new()
                }
            )
            .is_err()
        );
        let mut ordinary = draft.clone();
        ordinary.alternatives = None;
        ordinary.methods.source_roles = None;
        ordinary.methods.golden_gate_pool = None;
        let slot = ordinary.slots[0].clone();
        ordinary.slots.resize(alternatives::MAX_BINS, slot);
        ordinary.sync_rows();
        let full = alternatives::enable(&ordinary, "full-bins").unwrap();
        let rows = full
            .alternatives
            .as_ref()
            .unwrap()
            .rows
            .iter()
            .map(|row| AddedRow {
                construct_id: row.id.clone(),
                choice: 0,
                preparation: OverlapInput::default(),
                pool_substrate: 0,
            })
            .collect();
        assert!(
            propose(
                &full,
                Change::Add {
                    position: 0,
                    label: "Excess bin".into(),
                    choices: vec![seed.clone()],
                    rows
                }
            )
            .is_err()
        );
        let mut exhausted = draft.clone();
        exhausted.alternatives.as_mut().unwrap().next_id = u64::MAX;
        let rows = exhausted
            .alternatives
            .as_ref()
            .unwrap()
            .rows
            .iter()
            .map(|row| AddedRow {
                construct_id: row.id.clone(),
                choice: 0,
                preparation: OverlapInput::default(),
                pool_substrate: 0,
            })
            .collect();
        assert!(
            propose(
                &exhausted,
                Change::Add {
                    position: 0,
                    label: "No ID reuse".into(),
                    choices: vec![seed],
                    rows
                }
            )
            .is_err()
        );
        assert!(draft == original);
        assert_eq!(exhausted.alternatives.as_ref().unwrap().next_id, u64::MAX);
    }

    #[test]
    fn restore_preserves_exact_unsaved_frames_and_raises_allocation_highwater() {
        let original = fixture();
        let mut recovered = original.clone();
        recovered.name = "Unsaved recovery name".into();
        recovered.methods.rows[0].forward_tail = "RAW RECOVERY TEXT".into();
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { start: value, .. } =
                &mut recovered.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 1;
        }
        recovered.sync_rows();
        // Recovery may contain a working/saved discrepancy and a stale candidate
        // declaration. Both are raw intent; no save_active normalization occurs.
        recovered.validate().unwrap();
        let bin_id = ids(&original)[1].clone();
        let removal = propose(&original, Change::Remove { bin_id }).unwrap();
        let mut current = apply(&original, &removal, removal.revision())
            .unwrap()
            .draft;
        let removed_row = current.alternatives.as_ref().unwrap().rows[0].id.clone();
        current = alternatives::remove_row(&current, &removed_row).unwrap();
        let choices = current.alternatives.as_ref().unwrap().rows[0]
            .choices
            .clone();
        current = alternatives::add_row(&current, choices).unwrap();
        current.name = "Unsaved current name".into();
        current.methods.rows[0].forward_tail = "CURRENT RAW TEXT".into();
        let current_highwater = current.alternatives.as_ref().unwrap().next_id;
        let proposal = propose_restore(&current, &recovered, current_highwater + 50).unwrap();
        let mut expected = recovered.clone();
        expected.alternatives.as_mut().unwrap().next_id = current_highwater + 50;
        assert!(proposal.after() == &expected);
        assert!(proposal.before() == &current);
        assert!(
            proposal
                .impacts()
                .iter()
                .any(|impact| impact.construct_removed)
        );
        assert!(
            proposal
                .impacts()
                .iter()
                .any(|impact| impact.construct_added)
        );
        let active = proposal
            .impacts()
            .iter()
            .find(|impact| {
                impact.construct_id == recovered.alternatives.as_ref().unwrap().active_row
            })
            .unwrap();
        assert_eq!(active.name, "Unsaved recovery name");
        let applied = apply(&current, &proposal, proposal.revision()).unwrap();
        assert!(applied.recovery.as_ref() == &current);
        assert!(applied.draft == expected);
        assert_eq!(
            applied.draft.methods.rows[0].forward_tail,
            "RAW RECOVERY TEXT"
        );
        assert!(
            applied
                .draft
                .methods
                .golden_gate_pool
                .as_ref()
                .unwrap()
                .candidate
                .is_some()
        );
    }

    #[test]
    fn restore_rejects_cross_library_noops_and_stale_review_without_rewriting_raw_text() {
        let draft = fixture();
        assert!(propose_restore(&draft, &draft, 0).is_err());
        let mut ordinary = draft.clone();
        ordinary.alternatives = None;
        ordinary.methods.source_roles = None;
        let foreign = alternatives::enable(&ordinary, "foreign").unwrap();
        assert!(propose_restore(&draft, &foreign, u64::MAX).is_err());
        let mut recovered = draft.clone();
        recovered.methods.overlaps = " 004 , INVALID RAW ".into();
        let proposal = propose_restore(&draft, &recovered, 0).unwrap();
        assert_eq!(proposal.after().methods.overlaps, " 004 , INVALID RAW ");
        let mut stale = draft.clone();
        stale.name = "Changed after review".into();
        assert!(apply(&stale, &proposal, proposal.revision()).is_err());
        assert!(proposal.before() == &draft);
        assert_eq!(
            proposal.after().alternatives.as_ref().unwrap().next_id,
            draft
                .alternatives
                .as_ref()
                .unwrap()
                .next_id
                .max(recovered.alternatives.as_ref().unwrap().next_id)
        );
    }

    #[test]
    fn new_bin_assignments_require_the_active_working_pool_and_name_before_row_save() {
        let mut draft = fixture();
        draft.name = "Working unsaved construct".into();
        draft.methods.golden_gate_pool.as_mut().unwrap().reaction = true;
        let active_id = draft.alternatives.as_ref().unwrap().active_row.clone();
        let library = draft.alternatives.as_mut().unwrap();
        let index = library
            .rows
            .iter()
            .position(|row| row.id == active_id)
            .unwrap();
        library.rows[index].edits.methods.golden_gate_pool = None;
        library.rows[index].edits.name = "Old saved name".into();
        draft.validate().unwrap();
        let choices = vec![ChoiceSeed {
            slot: Slot::empty("New source"),
            preparation: OverlapInput::default(),
        }];
        let assigned = vec![Some(0); draft.alternatives.as_ref().unwrap().rows.len()];
        let mut substrates = vec![None; assigned.len()];
        let error = match addition_rows(&draft, &choices, &assigned, &substrates) {
            Ok(_) => panic!("working pool must require an explicit reaction substrate"),
            Err(error) => error,
        };
        assert!(
            error.contains("reaction substrate") && error.contains("Working unsaved construct")
        );
        assert!(!error.contains("Old saved name"));
        let row = &draft.alternatives.as_ref().unwrap().rows[index];
        assert_eq!(
            effective_edits(&draft, row).name,
            "Working unsaved construct"
        );
        substrates[index] = Some(1);
        let rows = addition_rows(&draft, &choices, &assigned, &substrates).unwrap();
        assert_eq!(rows[index].pool_substrate, 1);
        assert_eq!(rows[1 - index].pool_substrate, 0);
        // The reverse discrepancy must not demand a choice for an inactive pool.
        let saved_pool = draft.methods.golden_gate_pool.take();
        draft.alternatives.as_mut().unwrap().rows[index]
            .edits
            .methods
            .golden_gate_pool = saved_pool;
        substrates[index] = None;
        let rows = addition_rows(&draft, &choices, &assigned, &substrates).unwrap();
        assert_eq!(rows[index].pool_substrate, 0);
    }
}
