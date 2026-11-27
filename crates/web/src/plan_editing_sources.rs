//! Exact-event source choices for manual plan editing. Equal content supplies
//! display metadata only; it never combines or chooses distinct event identities.
use crate::source_identity::SourceInfo;
use collagene_core::{
    Alphabet, Molecule, Topology, history::DocumentHistory, run_lineage::EventBinding,
};
use std::collections::{BTreeMap, BTreeSet};
pub const MAX_DOCUMENTS: usize = 1000;
pub const MAX_EVENTS: usize = 50_000;
pub const MAX_HISTORICAL: usize = 4096;
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Catalog {
    pub choices: Vec<SourceInfo>,
    pub omitted: usize,
}
#[derive(Clone)]
struct Metadata {
    name: String,
    length: usize,
    alphabet: Alphabet,
    topology: Topology,
}
impl Metadata {
    fn from(m: &Molecule) -> Self {
        Self {
            name: m.name.clone(),
            length: m.sequence.len(),
            alphabet: m.sequence.alphabet(),
            topology: m.topology,
        }
    }
    fn choice(&self, binding: EventBinding) -> SourceInfo {
        SourceInfo {
            binding,
            name: self.name.clone(),
            length: self.length,
            alphabet: self.alphabet,
            topology: self.topology,
        }
    }
}
fn bounds(docs: &[DocumentHistory]) -> Result<(), String> {
    if docs.len() > MAX_DOCUMENTS {
        return Err("Source catalog exceeds 1,000 documents.".into());
    }
    let mut total = 0usize;
    let mut ids = BTreeSet::new();
    for doc in docs {
        if doc.id.trim().is_empty() || doc.id.len() > 128 || !ids.insert(doc.id.as_str()) {
            return Err(
                "Source catalog document identities are empty, duplicated or too long.".into(),
            );
        }
        total = total
            .checked_add(doc.event_count())
            .filter(|n| *n <= MAX_EVENTS)
            .ok_or("Source catalog exceeds 50,000 retained events.")?;
    }
    Ok(())
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
/// All current DNA heads are listed first. Historical DNA events are included
/// only for required original input revisions, in retained document/event order.
/// The historical cap never removes current heads or merges same-content events.
pub fn catalog(
    docs: &[DocumentHistory],
    required: &BTreeSet<String>,
    max_historical: usize,
) -> Result<Catalog, String> {
    bounds(docs)?;
    if max_historical > MAX_HISTORICAL
        || required.len() > 1000
        || required.iter().any(|s| !digest(s))
    {
        return Err("Source catalog historical/revision bounds exceeded.".into());
    }
    let mut choices = Vec::with_capacity(docs.len());
    // Repeated equal-content events reuse tiny display metadata, not source DNA.
    // Each choice below still receives its own exact document/event binding.
    let mut metadata = BTreeMap::<String, Metadata>::new();
    for doc in docs {
        let m = doc
            .revision_molecule_ref(doc.head_event())
            .map_err(|e| e.to_string())?;
        let info = Metadata::from(m);
        if info.alphabet == Alphabet::Dna {
            choices.push(info.choice(EventBinding {
                document_id: doc.id.clone(),
                event_id: doc.head_event().into(),
                revision: doc.head_revision().into(),
            }));
        }
        metadata.entry(doc.head_revision().into()).or_insert(info);
    }
    let mut historical = 0usize;
    let mut omitted = 0usize;
    for doc in docs {
        for event in doc.summaries() {
            if event.event_id == doc.head_event() || !required.contains(&event.revision) {
                continue;
            }
            if !metadata.contains_key(&event.revision) {
                let m = doc
                    .revision_molecule_ref(&event.event_id)
                    .map_err(|e| e.to_string())?;
                metadata.insert(event.revision.clone(), Metadata::from(m));
            }
            let info = metadata
                .get(&event.revision)
                .expect("retained metadata populated");
            if info.alphabet != Alphabet::Dna {
                continue;
            }
            if historical >= max_historical {
                omitted = omitted
                    .checked_add(1)
                    .ok_or("Source catalog count overflow.")?;
                continue;
            }
            choices.push(info.choice(EventBinding {
                document_id: doc.id.clone(),
                event_id: event.event_id,
                revision: event.revision,
            }));
            historical += 1;
        }
    }
    Ok(Catalog { choices, omitted })
}
/// Exact membership only. A same-content event in another document, another
/// event of this document or its current head cannot satisfy this binding.
pub fn resolve<'a>(
    docs: &'a [DocumentHistory],
    binding: &EventBinding,
) -> Result<&'a Molecule, String> {
    bounds(docs)?;
    if !digest(&binding.revision) || binding.event_id.is_empty() || binding.event_id.len() > 128 {
        return Err("Source event binding is malformed.".into());
    }
    let doc = docs
        .iter()
        .find(|d| d.id == binding.document_id)
        .ok_or("The exact source document is unavailable.")?;
    let molecule = doc
        .revision_molecule_ref(&binding.event_id)
        .map_err(|_| "The exact source event is not retained by this document.")?;
    if molecule.revision() != binding.revision {
        return Err("The exact source event revision does not match.".into());
    }
    Ok(molecule)
}

/// A target is an exact independently selected retained event, not the viewer head.
pub fn resolve_target<'a>(
    docs: &'a [DocumentHistory],
    binding: &EventBinding,
    expected_revision: &str,
) -> Result<&'a Molecule, String> {
    if binding.revision != expected_revision {
        return Err("The pinned independent target revision differs from this plan.".into());
    }
    let molecule = resolve(docs, binding)?;
    if molecule.sequence.alphabet() != Alphabet::Dna {
        return Err("Choose an independent DNA target event.".into());
    }
    Ok(molecule)
}
/// Resolve the same exact event from a native immutable archive before recording outputs.
pub fn workspace_target<'a>(
    workspace: &'a collagene_core::workspace::Workspace,
    binding: &EventBinding,
    expected_revision: &str,
) -> Result<&'a Molecule, String> {
    if binding.revision != expected_revision
        || !workspace
            .documents
            .iter()
            .any(|d| d.id == binding.document_id)
        || !workspace.history.iter().any(|e| {
            e.id == binding.event_id
                && e.document_id == binding.document_id
                && e.revision == binding.revision
        })
    {
        return Err(
            "The exact independent target event is unavailable or belongs to another plan.".into(),
        );
    }
    let molecule = workspace
        .revisions
        .get(&binding.revision)
        .ok_or("The independent target revision is unavailable.")?;
    if molecule.revision() != expected_revision || molecule.sequence.alphabet() != Alphabet::Dna {
        return Err("The independent target snapshot differs from this plan.".into());
    }
    Ok(molecule)
}
/// Capture an explicit target once when an unbound plan is opened. Later calls resolve it.
pub fn capture_target(
    workspace: &collagene_core::workspace::Workspace,
    preferred: Option<&EventBinding>,
    expected_revision: &str,
) -> Result<EventBinding, String> {
    let selected = match preferred {
        Some(binding) => binding.clone(),
        None => {
            let d = workspace
                .documents
                .iter()
                .find(|d| d.id == workspace.active_document)
                .ok_or("Choose an independent target document.")?;
            EventBinding {
                document_id: d.id.clone(),
                event_id: d
                    .head_event
                    .clone()
                    .ok_or("The target has no exact retained event.")?,
                revision: d.head.clone(),
            }
        }
    };
    workspace_target(workspace, &selected, expected_revision)?;
    Ok(selected)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn molecule(s: &str) -> Molecule {
        Molecule::dna("Same sequence name", s, Topology::Linear).unwrap()
    }
    fn binding(d: &DocumentHistory) -> EventBinding {
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event().into(),
            revision: d.head_revision().into(),
        }
    }
    #[test]
    fn current_a_and_retained_b_with_same_dna_remain_two_explicit_event_choices() {
        let a = DocumentHistory::new("doc-A", molecule("ACGTAC")).unwrap();
        let mut b = DocumentHistory::new("doc-B", molecule("ACGTAC")).unwrap();
        let old_b = binding(&b);
        b.commit(molecule("ACGTACA"), "Edited B").unwrap();
        let required = BTreeSet::from([a.head_revision().to_string()]);
        let docs = vec![a, b];
        let c = catalog(&docs, &required, 4096).unwrap();
        assert_eq!(c.choices.len(), 3);
        assert_eq!(c.omitted, 0);
        let matching = c
            .choices
            .iter()
            .filter(|s| required.contains(&s.binding.revision))
            .collect::<Vec<_>>();
        assert_eq!(matching.len(), 2);
        assert_ne!(matching[0].binding.event_id, matching[1].binding.event_id);
        assert_eq!(matching[1].binding, old_b);
        assert_eq!(resolve(&docs, &old_b).unwrap().sequence.as_str(), "ACGTAC");
        assert!(
            std::ptr::eq(
                resolve(&docs, &old_b).unwrap(),
                docs[1].revision_molecule_ref(&old_b.event_id).unwrap()
            ),
            "resolution borrows the immutable snapshot, not a DNA clone"
        );
    }
    #[test]
    fn repeated_content_in_one_document_keeps_distinct_historical_events() {
        let mut d = DocumentHistory::new("doc-A", molecule("ACGTAC")).unwrap();
        let first = binding(&d);
        d.commit(molecule("ACGTACA"), "Temporary edit").unwrap();
        d.commit(molecule("ACGTAC"), "Return to earlier content")
            .unwrap();
        let final_event = binding(&d);
        assert_eq!(first.revision, final_event.revision);
        assert_ne!(first.event_id, final_event.event_id);
        let docs = vec![d];
        let required = BTreeSet::from([first.revision.clone()]);
        let c = catalog(&docs, &required, 4096).unwrap();
        assert_eq!(c.choices.len(), 2);
        assert_eq!(c.choices[0].binding, final_event);
        assert_eq!(c.choices[1].binding, first);
        resolve(&docs, &c.choices[0].binding).unwrap();
        resolve(&docs, &c.choices[1].binding).unwrap();
    }
    #[test]
    fn forged_document_foreign_event_and_revision_mismatch_never_fall_back_to_content() {
        let a = DocumentHistory::new("doc-A", molecule("ACGTAC")).unwrap();
        let b = DocumentHistory::new("doc-B", molecule("ACGTAC")).unwrap();
        let exact = binding(&a);
        let docs = vec![a, b];
        let mut forged = exact.clone();
        forged.document_id = "missing".into();
        assert!(resolve(&docs, &forged).is_err());
        forged = exact.clone();
        forged.document_id = "doc-B".into();
        assert!(
            resolve(&docs, &forged).is_err(),
            "same content in B does not own A's event"
        );
        forged = exact.clone();
        forged.event_id = "event-unknown".into();
        assert!(resolve(&docs, &forged).is_err());
        forged = exact;
        forged.revision = "0".repeat(64);
        assert!(resolve(&docs, &forged).is_err());
    }
    #[test]
    fn historical_cap_reports_all_omissions_without_hiding_current_heads_or_nonrequired_events() {
        let mut a = DocumentHistory::new("doc-A", molecule("ACGTAC")).unwrap();
        let required = BTreeSet::from([a.head_revision().to_string()]);
        a.commit(molecule("ACGTACA"), "Other revision").unwrap();
        a.commit(molecule("ACGTAC"), "Same original revision")
            .unwrap();
        a.commit(molecule("ACGTACAA"), "Current revision").unwrap();
        let b = DocumentHistory::new("doc-B", molecule("GGCC")).unwrap();
        let docs = vec![a, b];
        let none = catalog(&docs, &required, 0).unwrap();
        assert_eq!(none.choices.len(), 2);
        assert_eq!(none.omitted, 2);
        let one = catalog(&docs, &required, 1).unwrap();
        assert_eq!(one.choices.len(), 3);
        assert_eq!(one.omitted, 1);
        let all = catalog(&docs, &required, 2).unwrap();
        assert_eq!(all.choices.len(), 4);
        assert_eq!(all.omitted, 0);
        assert_eq!(
            catalog(&docs, &BTreeSet::new(), 0).unwrap().choices.len(),
            2
        );
        assert!(catalog(&docs, &required, 4097).is_err());
        let too_many = vec![docs[0].clone(); 1001];
        assert!(catalog(&too_many, &required, 0).is_err());
        let duplicate = vec![docs[0].clone(), docs[0].clone()];
        assert!(resolve(&duplicate, &binding(&docs[0])).is_err());
    }
    #[test]
    fn retained_target_and_source_do_not_follow_later_heads_or_equal_content_aliases() {
        use collagene_core::workspace::Workspace;
        let mut source = DocumentHistory::new("source", molecule("AGTCCGAT")).unwrap();
        let source_event = binding(&source);
        let mut target = DocumentHistory::new("target", molecule("AGTCCGAT")).unwrap();
        let target_event = binding(&target);
        source
            .commit(molecule("TTTTAAAA"), "Later source edit")
            .unwrap();
        target
            .commit(molecule("CCCCGGGG"), "Later target edit")
            .unwrap();
        let docs = vec![source, target];
        let before = docs.iter().map(binding).collect::<Vec<_>>();
        assert_eq!(
            resolve(&docs, &source_event).unwrap().sequence.as_str(),
            "AGTCCGAT"
        );
        assert_eq!(
            resolve_target(&docs, &target_event, &target_event.revision)
                .unwrap()
                .sequence
                .as_str(),
            "AGTCCGAT"
        );
        let mut workspace =
            Workspace::from_documents("Unrelated viewer", vec![molecule("AAAATTTT")], 0).unwrap();
        for doc in &docs {
            doc.contribute(&mut workspace).unwrap();
        }
        workspace.validate().unwrap();
        assert_eq!(
            capture_target(&workspace, Some(&target_event), &target_event.revision).unwrap(),
            target_event
        );
        assert!(capture_target(&workspace, None, &target_event.revision).is_err());
        assert_eq!(
            workspace_target(&workspace, &target_event, &target_event.revision)
                .unwrap()
                .sequence
                .as_str(),
            "AGTCCGAT"
        );
        let mut wrong = target_event.clone();
        wrong.document_id = source_event.document_id.clone();
        assert!(resolve_target(&docs, &wrong, &target_event.revision).is_err());
        assert!(workspace_target(&workspace, &wrong, &target_event.revision).is_err());
        assert!(workspace_target(&workspace, &target_event, &"0".repeat(64)).is_err());
        assert_eq!(docs.iter().map(binding).collect::<Vec<_>>(), before);
    }
}
