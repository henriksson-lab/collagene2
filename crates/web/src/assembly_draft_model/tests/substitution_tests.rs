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
    let library = crate::assembly_alternatives::enable(&draft, "patch-library").unwrap();
    assert!(decode(&encode_raw(&library).unwrap()).unwrap() == library);
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
