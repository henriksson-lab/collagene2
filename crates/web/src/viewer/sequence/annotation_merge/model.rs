//! UI event/epoch/proposal guards; native merge owns all relationship semantics.
use collagene_core::{
    material::lineage::EventBinding,
    sequence::annotations::merge::{Decision, Selection},
};
use std::collections::BTreeMap;

#[derive(Clone)]
pub(super) struct Context {
    pub event: EventBinding,
    pub epoch: u64,
    pub input: (u32, Vec<String>),
}
impl Context {
    pub fn validate_current(
        &self,
        event: &EventBinding,
        epoch: u64,
        input: &(u32, Vec<String>),
        locked: bool,
    ) -> Result<(), String> {
        if locked || &self.event != event || self.epoch != epoch || &self.input != input {
            return Err("The source event, proposals or workspace changed. Preview again.".into());
        }
        Ok(())
    }
}
pub(super) fn selections(
    original: &[Selection],
    choices: &BTreeMap<String, Decision>,
) -> Result<Vec<Selection>, String> {
    if original.len() != choices.len() {
        return Err("Choose one explicit decision for every current proposal.".into());
    }
    original
        .iter()
        .map(|selection| {
            let decision = choices
                .get(&selection.proposal_id)
                .cloned()
                .ok_or("A current proposal has no explicit decision.")?;
            Ok(Selection {
                proposal_id: selection.proposal_id.clone(),
                decision,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::documents::workspace::Workspace;
    use collagene_model::{Molecule, Topology};

    #[test]
    fn identical_calls_do_not_allow_another_document_event_or_workspace_context() {
        let molecule = Molecule::dna("Identical calls", "ACGT", Topology::Linear).unwrap();
        let workspace =
            Workspace::from_documents("Merge", vec![molecule.clone(), molecule], 0).unwrap();
        let event = |index: usize| {
            let document = &workspace.documents[index];
            EventBinding {
                document_id: document.id.clone(),
                event_id: document.head_event.clone().unwrap(),
                revision: document.head.clone(),
            }
        };
        assert_eq!(event(0).revision, event(1).revision);
        let context = Context {
            event: event(0),
            epoch: 5,
            input: (2, vec!["one".into()]),
        };
        context
            .validate_current(&event(0), 5, &context.input, false)
            .unwrap();
        assert!(
            context
                .validate_current(&event(1), 5, &context.input, false)
                .is_err()
        );
        assert!(
            context
                .validate_current(&event(0), 6, &context.input, false)
                .is_err()
        );
        assert!(
            context
                .validate_current(&event(0), 5, &(3, vec!["one".into()]), false)
                .is_err()
        );
        assert!(
            context
                .validate_current(&event(0), 5, &(2, vec!["two".into()]), false)
                .is_err()
        );
        assert!(
            context
                .validate_current(&event(0), 5, &context.input, true)
                .is_err()
        );
    }
    #[test]
    fn missing_extra_and_changed_proposal_choices_never_default_to_append() {
        let original = vec![
            Selection {
                proposal_id: "one".into(),
                decision: Decision::Skip,
            },
            Selection {
                proposal_id: "two".into(),
                decision: Decision::Skip,
            },
        ];
        let mut choices = BTreeMap::from([("one".into(), Decision::Append)]);
        assert!(selections(&original, &choices).is_err());
        choices.insert(
            "two".into(),
            Decision::AddProvenance {
                target_index: 3,
                annotation_revision: "a".repeat(64),
            },
        );
        let selected = selections(&original, &choices).unwrap();
        assert_eq!(selected[0].decision, Decision::Append);
        assert!(matches!(
            selected[1].decision,
            Decision::AddProvenance {
                target_index: 3,
                ..
            }
        ));
        choices.insert("other".into(), Decision::Skip);
        assert!(selections(&original, &choices).is_err());
        choices.remove("one");
        assert!(selections(&original, &choices).is_err());
    }
}
