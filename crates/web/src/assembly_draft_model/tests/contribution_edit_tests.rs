use super::*;
use collagene_core::{assembly_edit_intent::Edit, assembly_substitution_intent::Substitution};
fn edited() -> Draft {
    let mut d = Draft::default();
    let roles = d.methods.role_ids(d.slots.len()).unwrap();
    d.edits.push(Edit {
        role: roles[0].role.clone(),
        start: 2,
        end: 4,
        expected: "CC".into(),
        replacement: "N".into(),
    });
    d.methods.source_roles = Some(roles);
    d.sync_rows();
    d
}
#[test]
fn persisted_general_edits_require_stable_roles_and_bounded_literals() {
    let d = edited();
    assert_eq!(d.version, CURRENT_VERSION);
    assert!(d.validate().is_ok());
    assert!(decode(&encode_raw(&d).unwrap()).unwrap() == d);
    let mut wrong = d.clone();
    wrong.version = CURRENT_VERSION - 1;
    assert!(wrong.validate().is_err());
    wrong = d.clone();
    wrong.methods.source_roles = None;
    assert!(wrong.validate().is_err());
    wrong = d.clone();
    wrong.edits[0].expected = "C".into();
    assert!(wrong.validate().is_err());
    wrong = d.clone();
    wrong.edits[0].replacement = "Z".into();
    assert!(wrong.validate().is_err());
    wrong = d.clone();
    wrong.edits[0].replacement = "A".repeat(100_001);
    assert!(wrong.validate().is_err());
    let mut oversized = serde_json::to_value(&d).unwrap();
    oversized["edits"] = serde_json::to_value(vec![d.edits[0].clone(); 257]).unwrap();
    assert!(decode(&serde_json::to_vec(&oversized).unwrap()).is_err());
}
#[test]
fn general_edits_reject_ambiguous_order_and_cross_kind_overlap() {
    let d = edited();
    let role = d.edits[0].role.clone();
    let mut next = d.clone();
    next.edits.push(Edit {
        role: role.clone(),
        start: 4,
        end: 5,
        expected: "G".into(),
        replacement: "T".into(),
    });
    assert!(next.validate().is_ok()); // Adjacent nonempty intervals have an unambiguous order.
    next.edits[1] = Edit {
        role: role.clone(),
        start: 4,
        end: 4,
        expected: "".into(),
        replacement: "T".into(),
    };
    assert!(next.validate().is_err());
    next = d.clone();
    next.substitutions.push(Substitution {
        role: role.clone(),
        offset: 3,
        expected: "C".into(),
        replacement: "G".into(),
    });
    assert!(next.validate().is_err());
    next = d.clone();
    next.substitutions.push(Substitution {
        role: role.clone(),
        offset: 4,
        expected: "G".into(),
        replacement: "T".into(),
    });
    assert!(next.validate().is_ok());
    next.edits[0] = Edit {
        role,
        start: 4,
        end: 4,
        expected: "".into(),
        replacement: "C".into(),
    };
    assert!(next.validate().is_err());
}
