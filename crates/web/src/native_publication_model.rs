//! Pure saved-source publication controls; no transport or biological authority.
use collagene_desktop_sync_protocol::{self as wire, library as base};
use std::collections::{BTreeMap, BTreeSet};
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SavedPin {
    pub document: base::DocumentId,
    pub revision: base::RevisionId,
    pub dirty: bool,
    pub name: String,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Context {
    pub source: Option<SavedPin>,
    pub remote: wire::Remote,
    pub destination: Option<base::FolderId>,
}
impl Context {
    pub fn clean_source(&self) -> Result<&SavedPin, String> {
        self.source.as_ref().filter(|s|!s.dirty).ok_or_else(||"Save the current plasmid in the local library first, then review its saved revision.".into())
    }
}
pub(crate) fn configuration_matches(
    actual: &wire::PublicationRemote,
    expected: &wire::Remote,
) -> bool {
    actual.id == expected.id
        && actual.name == expected.name
        && actual.endpoint == expected.endpoint
        && actual.lab == expected.expected_lab
}
pub(crate) fn deleted_review_matches(
    item: &wire::DeletedFolderPublication,
    description: &wire::PublicationDescription,
) -> bool {
    matches!((&description.target,&description.action),
        (wire::PublicationTarget::DeletedFolder(local),wire::MetadataAction::DeleteEmptyFolder{id,expected_generation})
        if local==&item.local_folder && id==&item.remote_folder.id && *expected_generation==item.remote_folder.generation)
}
pub(crate) type ChoiceKey = (String, String);
pub(crate) type Chosen = BTreeMap<ChoiceKey, BTreeSet<String>>;
pub(crate) fn choices(
    requirements: &wire::PushChoices,
    chosen: &Chosen,
) -> Result<base::Rows<wire::PushChoice, 128>, String> {
    let mut output = Vec::new();
    for requirement in requirements.requirements.iter() {
        let key = (
            requirement.original_hash.as_str().into(),
            requirement.component_iri.as_str().into(),
        );
        let mut selected = chosen.get(&key).cloned().unwrap_or_default();
        selected.extend(
            requirement
                .required_sequences
                .iter()
                .map(|s| s.as_str().to_owned()),
        );
        if selected.is_empty() {
            return Err("Choose at least one sequence for each reachable component with alternatives. No first sequence is selected automatically.".into());
        }
        if selected.iter().any(|s| {
            !requirement
                .available_sequences
                .iter()
                .any(|v| v.as_str() == s)
        }) {
            return Err("A sequence choice is outside the inspected saved-source alternatives. Inspect again.".into());
        }
        output.push(wire::PushChoice {
            original_hash: requirement.original_hash.clone(),
            component_iri: requirement.component_iri.clone(),
            sequences: base::Rows::new(
                selected
                    .iter()
                    .map(base::Text::new)
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|_| "Sequence identifiers exceed their native bound.")?,
            )
            .map_err(|_| "Too many selected alternatives.")?,
        });
    }
    if chosen.keys().any(|key| {
        !requirements
            .requirements
            .iter()
            .any(|r| r.original_hash.as_str() == key.0 && r.component_iri.as_str() == key.1)
    }) {
        return Err("Sequence choices belong to a different saved-source inspection.".into());
    }
    base::Rows::new(output).map_err(|_| "Too many reachable components.".into())
}
#[cfg(test)]
mod tests {
    use super::*;
    fn requirements() -> wire::PushChoices {
        wire::PushChoices {
            local_document: base::DocumentId::new(format!("doc-{:032x}", 1)).unwrap(),
            local_revision: base::RevisionId::new("11".repeat(32)).unwrap(),
            remote: wire::RemoteId::new(format!("rem-{:032x}", 2)).unwrap(),
            requirements: base::Rows::new(vec![wire::PushChoiceRequirement {
                original_hash: base::Digest::new("22".repeat(32)).unwrap(),
                component_iri: base::Text::new("https://example.test/c").unwrap(),
                available_sequences: base::Rows::new(vec![
                    base::Text::new("https://example.test/a").unwrap(),
                    base::Text::new("https://example.test/b").unwrap(),
                ])
                .unwrap(),
                required_sequences: base::Rows::new(vec![]).unwrap(),
            }])
            .unwrap(),
        }
    }
    #[test]
    fn alternatives_require_explicit_selection_and_preserve_required_native_sources() {
        let mut r = requirements();
        assert!(choices(&r, &Chosen::new()).is_err());
        let key = ("22".repeat(32), "https://example.test/c".into());
        let mut chosen = Chosen::from([(
            key.clone(),
            BTreeSet::from(["https://example.test/b".into()]),
        )]);
        let out = choices(&r, &chosen).unwrap();
        assert_eq!(out[0].sequences[0].as_str(), "https://example.test/b");
        let mut row = r.requirements[0].clone();
        row.required_sequences =
            base::Rows::new(vec![base::Text::new("https://example.test/a").unwrap()]).unwrap();
        r.requirements = base::Rows::new(vec![row]).unwrap();
        assert_eq!(choices(&r, &chosen).unwrap()[0].sequences.len(), 2);
        chosen
            .get_mut(&key)
            .unwrap()
            .insert("https://private.invalid/unavailable".into());
        assert!(choices(&r, &chosen).is_err());
    }
    #[test]
    fn named_deleted_choice_requires_the_exact_remote_folder_and_generation() {
        let local = base::FolderId::new(format!("fld-{:032x}", 1)).unwrap();
        let remote = base::FolderId::new(format!("fld-{:032x}", 2)).unwrap();
        let lab = base::LabId::new(format!("lab-{:032x}", 3)).unwrap();
        let item = wire::DeletedFolderPublication {
            local_folder: local.clone(),
            remote_folder: base::Folder {
                id: remote.clone(),
                parent: None,
                name: base::Text::new("Known previously synchronized folder").unwrap(),
                generation: 7,
            },
        };
        let mut description = wire::PublicationDescription {
            local_lab: lab.clone(),
            target: wire::PublicationTarget::DeletedFolder(local),
            remote: wire::PublicationRemote {
                id: wire::RemoteId::new(format!("rem-{:032x}", 4)).unwrap(),
                name: base::Text::new("Lab").unwrap(),
                endpoint: base::Text::new("https://lab.example").unwrap(),
                lab,
            },
            expected_checkpoint: None,
            local_generation: None,
            local_name: None,
            local_parent: None,
            remote_parent: None,
            parent_binding: None,
            action: wire::MetadataAction::DeleteEmptyFolder {
                id: remote.clone(),
                expected_generation: 7,
            },
        };
        assert!(deleted_review_matches(&item, &description));
        description.action = wire::MetadataAction::DeleteEmptyFolder {
            id: remote,
            expected_generation: 8,
        };
        assert!(!deleted_review_matches(&item, &description));
        description.action = wire::MetadataAction::DeleteEmptyFolder {
            id: base::FolderId::new(format!("fld-{:032x}", 9)).unwrap(),
            expected_generation: 7,
        };
        assert!(!deleted_review_matches(&item, &description));
    }
    #[test]
    fn source_and_full_remote_configuration_bind_the_review() {
        let r = requirements();
        let mut context = Context {
            source: Some(SavedPin {
                document: r.local_document,
                revision: r.local_revision,
                dirty: false,
                name: "Saved plasmid".into(),
            }),
            remote: wire::Remote {
                id: r.remote,
                name: base::Text::new("Lab").unwrap(),
                endpoint: base::Text::new("https://lab.example").unwrap(),
                expected_lab: base::LabId::new(format!("lab-{:032x}", 3)).unwrap(),
            },
            destination: None,
        };
        let original = context.clone();
        let mut actual = wire::PublicationRemote {
            id: context.remote.id.clone(),
            name: context.remote.name.clone(),
            endpoint: context.remote.endpoint.clone(),
            lab: context.remote.expected_lab.clone(),
        };
        assert!(configuration_matches(&actual, &context.remote));
        actual.endpoint = base::Text::new("https://different.example").unwrap();
        assert!(!configuration_matches(&actual, &context.remote));
        assert!(context.clean_source().is_ok());
        context.source.as_mut().unwrap().dirty = true;
        assert!(context.clean_source().is_err());
        assert_ne!(context, original);
        context = original.clone();
        context.remote.expected_lab = base::LabId::new(format!("lab-{:032x}", 4)).unwrap();
        assert_ne!(context, original);
        context = original.clone();
        context.source.as_mut().unwrap().revision = base::RevisionId::new("33".repeat(32)).unwrap();
        assert_ne!(context, original);
    }
}
