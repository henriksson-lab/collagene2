//! Exact retained-event choices. Display names and DNA equality never select ancestry.
use collagene_core::{documents::workspace::Workspace, material::lineage::EventBinding};
use collagene_model::{Alphabet, Molecule};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct EventChoice {
    pub binding: EventBinding,
    pub label: String,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CdsChoice {
    pub index: usize,
    pub label: String,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Form {
    pub source_event: String,
    pub source_feature: String,
    pub product_feature: String,
    pub rationale: String,
    pub genetic_code: String,
}
impl Default for Form {
    fn default() -> Self {
        Self {
            source_event: String::new(),
            source_feature: String::new(),
            product_feature: String::new(),
            rationale: String::new(),
            genetic_code: "1".into(),
        }
    }
}
pub(super) fn molecule<'a>(
    workspace: &'a Workspace,
    binding: &EventBinding,
) -> Result<&'a Molecule, String> {
    let event = workspace
        .history
        .iter()
        .find(|event| event.id == binding.event_id)
        .ok_or("The exact retained event is unavailable.")?;
    if event.document_id != binding.document_id || event.revision != binding.revision {
        return Err("The retained event has a different document or revision.".into());
    }
    let molecule = workspace
        .revisions
        .get(&binding.revision)
        .ok_or("The retained event's sequence snapshot is unavailable.")?;
    Ok(molecule)
}
pub(super) fn events(workspace: &Workspace) -> Vec<EventChoice> {
    workspace
        .history
        .iter()
        .filter_map(|event| {
            let molecule = workspace.revisions.get(&event.revision)?;
            if molecule.sequence.alphabet() != Alphabet::Dna {
                return None;
            }
            Some(EventChoice {
                binding: EventBinding {
                    document_id: event.document_id.clone(),
                    event_id: event.id.clone(),
                    revision: event.revision.clone(),
                },
                label: format!(
                    "{} · {} · {}",
                    molecule.name.chars().take(64).collect::<String>(),
                    event.document_id.chars().take(64).collect::<String>(),
                    event.id
                ),
            })
        })
        .collect()
}
pub(super) fn cds(molecule: &Molecule) -> Vec<CdsChoice> {
    molecule
        .features
        .iter()
        .enumerate()
        .filter(|(_, feature)| feature.kind.eq_ignore_ascii_case("CDS"))
        .map(|(index, feature)| CdsChoice {
            index,
            label: format!(
                "{} · {}{}",
                index + 1,
                feature.label().chars().take(96).collect::<String>(),
                if feature.needs_review {
                    " · needs review"
                } else {
                    ""
                }
            ),
        })
        .collect()
}
impl Form {
    pub fn source<'a>(&self, choices: &'a [EventChoice]) -> Result<&'a EventBinding, String> {
        choices
            .iter()
            .find(|choice| choice.binding.event_id == self.source_event)
            .map(|choice| &choice.binding)
            .ok_or_else(|| "Choose an exact source event from retained history.".into())
    }
}

const MENU_ROWS: usize = 100;
/// Keep the exact selected event visible without exceeding the menu bound.
pub(super) fn event_menu<'a>(
    choices: &'a [EventChoice],
    query: &str,
    selected: &str,
) -> Vec<&'a EventChoice> {
    let query = query.to_lowercase();
    let chosen = choices
        .iter()
        .find(|choice| choice.binding.event_id == selected);
    let remaining = MENU_ROWS - usize::from(chosen.is_some());
    chosen
        .into_iter()
        .chain(
            choices
                .iter()
                .filter(|choice| {
                    choice.binding.event_id != selected
                        && (choice.label.to_lowercase().contains(&query)
                            || choice.binding.document_id.to_lowercase().contains(&query)
                            || choice.binding.event_id.to_lowercase().contains(&query))
                })
                .take(remaining),
        )
        .collect()
}
pub(super) fn cds_menu<'a>(
    choices: &'a [CdsChoice],
    query: &str,
    selected: &str,
) -> Vec<&'a CdsChoice> {
    let query = query.to_lowercase();
    let index = selected.parse::<usize>().ok();
    let chosen = choices.iter().find(|choice| Some(choice.index) == index);
    let remaining = MENU_ROWS - usize::from(chosen.is_some());
    chosen
        .into_iter()
        .chain(
            choices
                .iter()
                .filter(|choice| {
                    Some(choice.index) != index && choice.label.to_lowercase().contains(&query)
                })
                .take(remaining),
        )
        .collect()
}
pub(super) fn navigation_ranges_valid(ranges: &[(usize, usize)], length: usize) -> bool {
    !ranges.is_empty()
        && ranges
            .iter()
            .all(|(start, end)| start < end && *end <= length)
}
