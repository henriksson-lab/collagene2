use super::*;
use crate::cloning::restriction::catalog::select;

fn options(pairs: &[(&str, &str)]) -> ReactionConditions {
    ReactionConditions {
        schema_version: 1,
        formulations: pairs
            .iter()
            .map(|(e, f)| ChosenFormulation {
                enzyme_id: (*e).into(),
                formulation_id: Some((*f).into()),
            })
            .collect(),
        buffer: "rCutSmart".into(),
        temperature_milli_c: 37000,
        incubation_minutes: 60,
        minimum_activity_percent: 75,
        activity_context: Some("NEB unit-assay relative activity".into()),
        temperature_context: Some("vendor recommended incubation".into()),
    }
}
fn selected(c: &catalog::Catalog, o: &ReactionConditions) -> CatalogSelection {
    select(
        c,
        &o.formulations
            .iter()
            .map(|f| f.enzyme_id.clone())
            .collect::<Vec<_>>(),
    )
    .unwrap()
}
fn modified(mut change: impl FnMut(&mut catalog::Edition)) -> catalog::Catalog {
    let c = catalog::Catalog::builtin().unwrap();
    let mut edition = c.edition().clone();
    change(&mut edition);
    catalog::Catalog::from_edition(catalog::seal(edition).unwrap()).unwrap()
}

#[test]
fn ecori_and_hf_are_distinct_assay_facts_and_star_absence_is_never_inferred() {
    let c = catalog::Catalog::builtin().unwrap();
    let standard = options(&[("ecori", "neb-r0101")]);
    let hf = options(&[("ecori", "neb-r3101")]);
    let selection = selected(&c, &standard);
    let a = assess(&c, &selection, &standard).unwrap();
    let b = assess(&c, &selection, &hf).unwrap();
    assert_eq!(a.documentary_status, EvidenceStatus::BelowThreshold); // literal 50%
    assert_eq!(b.documentary_status, EvidenceStatus::Known); // literal 100%
    let f = &a.formulations[0];
    let row = &f.formulation.as_ref().unwrap().evidence[f.activity.evidence_indices[0]];
    assert_eq!(row.value, catalog::FactValue::Percent(50));
    assert!(
        row.qualifiers
            .iter()
            .any(|q| q.contains("possible star activity"))
    );
    assert_eq!(f.star_risk, Unassessed::Unknown);
    assert_eq!(b.formulations[0].star_risk, Unassessed::Unknown);
    assert_eq!(b.formulations[0].incubation_duration, Unassessed::Unknown);
    assert_ne!(a.options_revision, b.options_revision);
    assert_ne!(a.batch_revision, b.batch_revision);
    assert!(!a.sources.is_empty());
}

#[test]
fn joint_digest_distinguishes_recorded_failure_missing_evidence_and_common_candidates() {
    let c = catalog::Catalog::builtin().unwrap();
    let mut o = options(&[("ecori", "neb-r0101"), ("smai", "neb-r0141")]);
    let selection = selected(&c, &o);
    let r = assess(&c, &selection, &o).unwrap();
    assert_eq!(r.documentary_status, EvidenceStatus::BelowThreshold);
    assert!(r.common_conditions.is_empty()); // EcoRI 50% rCutSmart; SmaI <10% elsewhere
    assert!(r.suggestions_complete);
    o.buffer = "unrecorded buffer".into();
    let unknown = assess(&c, &selection, &o).unwrap();
    assert_eq!(unknown.documentary_status, EvidenceStatus::Unknown);
    let hf = options(&[("ecori", "neb-r3101"), ("smai", "neb-r0141")]);
    let r = assess(&c, &selection, &hf).unwrap();
    assert_eq!(r.documentary_status, EvidenceStatus::Known);
    assert_eq!(r.common_conditions.len(), 1);
    assert_eq!(r.common_conditions[0].buffer, "rCutSmart");
    assert_eq!(r.common_conditions[0].temperature_milli_c, 37000);
    assert_eq!(
        r.common_conditions[0].activity_context,
        "NEB unit-assay relative activity"
    );
}

#[test]
fn changing_temperature_context_or_time_cannot_create_activity_evidence() {
    let c = catalog::Catalog::builtin().unwrap();
    let mut o = options(&[("ecori", "neb-r3101")]);
    let s = selected(&c, &o);
    o.temperature_milli_c = 25000;
    let r = assess(&c, &s, &o).unwrap();
    assert_eq!(r.formulations[0].activity.status, EvidenceStatus::Unknown);
    assert_eq!(
        r.formulations[0].temperature.status,
        EvidenceStatus::TemperatureMismatch
    );
    assert!(
        r.common_conditions
            .iter()
            .all(|v| v.temperature_milli_c == 37000)
    );
    o.temperature_milli_c = 37000;
    o.activity_context = Some("unrecorded assay context".into());
    assert_eq!(
        assess(&c, &s, &o).unwrap().formulations[0].activity.status,
        EvidenceStatus::Unknown
    );
    o.activity_context = Some("NEB unit-assay relative activity".into());
    o.incubation_minutes = 1;
    let short = assess(&c, &s, &o).unwrap();
    o.incubation_minutes = 100_800;
    let long = assess(&c, &s, &o).unwrap();
    assert_eq!(short.formulations, long.formulations);
    assert_ne!(short.options_revision, long.options_revision);
    assert_eq!(
        long.formulations[0].incubation_duration,
        Unassessed::Unknown
    );
}

#[test]
fn missing_formulation_remains_unknown_without_affecting_other_selected_evidence() {
    let c = catalog::Catalog::builtin().unwrap();
    let mut o = options(&[("ecori", "neb-r3101"), ("smai", "neb-r0141")]);
    o.formulations[1].formulation_id = None;
    let s = selected(&c, &o);
    let r = assess(&c, &s, &o).unwrap();
    assert_eq!(r.formulations[0].activity.status, EvidenceStatus::Known);
    assert_eq!(r.formulations[1].activity.status, EvidenceStatus::Unknown);
    assert_eq!(
        r.formulations[1].temperature.status,
        EvidenceStatus::Unknown
    );
    assert!(r.formulations[1].formulation.is_none());
    assert_eq!(r.documentary_status, EvidenceStatus::Unknown);
    assert!(r.common_conditions.is_empty());
    r.validate_for(&c, &s, &o).unwrap();
}

#[test]
fn multiple_contexts_require_choice_and_conflicting_facts_remain_conflicts() {
    let c = modified(|e| {
        let f = e
            .formulations
            .iter_mut()
            .find(|f| f.id == "neb-r3101")
            .unwrap();
        let mut row=f.evidence.iter().find(|e|matches!(&e.condition,catalog::Condition::BufferActivity {buffer,..} if buffer=="rCutSmart")).unwrap().clone();
        row.context = "another unit assay".into();
        row.value = catalog::FactValue::Percent(25);
        f.evidence.push(row);
    });
    let mut o = options(&[("ecori", "neb-r3101")]);
    let s = selected(&c, &o);
    o.activity_context = None;
    assert_eq!(
        assess(&c, &s, &o).unwrap().formulations[0].activity.status,
        EvidenceStatus::ContextRequired
    );
    o.activity_context = Some("another unit assay".into());
    assert_eq!(
        assess(&c, &s, &o).unwrap().documentary_status,
        EvidenceStatus::BelowThreshold
    );
    let c = modified(|e| {
        let f = e
            .formulations
            .iter_mut()
            .find(|f| f.id == "neb-r3101")
            .unwrap();
        let mut row=f.evidence.iter().find(|e|matches!(&e.condition,catalog::Condition::BufferActivity {buffer,..} if buffer=="rCutSmart")).unwrap().clone();
        row.value = catalog::FactValue::Percent(25);
        f.evidence.push(row);
    });
    let o = options(&[("ecori", "neb-r3101")]);
    let s = selected(&c, &o);
    assert_eq!(
        assess(&c, &s, &o).unwrap().documentary_status,
        EvidenceStatus::Conflict
    );
}

#[test]
fn less_than_values_do_not_imply_a_lower_bound() {
    let c = catalog::Catalog::builtin().unwrap();
    let mut o = options(&[("ecori", "neb-r3101")]);
    let s = selected(&c, &o);
    o.buffer = "NEBuffer r3.1".into();
    o.minimum_activity_percent = 10;
    assert_eq!(
        assess(&c, &s, &o).unwrap().documentary_status,
        EvidenceStatus::BelowThreshold
    );
    o.minimum_activity_percent = 5;
    assert_eq!(
        assess(&c, &s, &o).unwrap().documentary_status,
        EvidenceStatus::Unknown
    );
}

#[test]
fn exact_membership_bounds_and_fresh_validation_reject_tampered_reports() {
    let c = catalog::Catalog::builtin().unwrap();
    let o = options(&[("ecori", "neb-r3101")]);
    let s = selected(&c, &o);
    let mut r = assess(&c, &s, &o).unwrap();
    r.validate_for(&c, &s, &o).unwrap();
    assert_eq!(
        revision(&c, &s, &o).unwrap(),
        documentary_revision(&s, &o).unwrap()
    );
    r.validate_binding(&revision(&c, &s, &o).unwrap()).unwrap();
    let mut other = o.clone();
    other.incubation_minutes += 1;
    assert!(
        r.validate_binding(&documentary_revision(&s, &other).unwrap())
            .is_err()
    );
    r.formulations[0]
        .formulation
        .as_mut()
        .unwrap()
        .formulation_scope = "fabricated scope".into();
    r.batch_revision = r.revision().unwrap();
    r.validate().unwrap(); // internally consistent history is not fresh authority
    assert!(r.validate_for(&c, &s, &o).is_err());
    let mut wrong = o.clone();
    wrong.formulations[0].formulation_id = Some("EcoRI-HF".into()); // name is not exact ID
    assert!(assess(&c, &s, &wrong).is_err());
    wrong.formulations[0].formulation_id = Some("neb-r0141".into()); // wrong enzyme
    assert!(assess(&c, &s, &wrong).is_err());
    wrong = o.clone();
    wrong.formulations = vec![wrong.formulations[0].clone(); 65];
    assert!(wrong.validate(&s).is_err());
    wrong = o.clone();
    wrong.buffer = "x".repeat(MAX_OPTIONS_BYTES + 1);
    assert!(matches!(wrong.validate(&s), Err(Error::TooLarge)));
    let c = modified(|e| {
        let f = e
            .formulations
            .iter_mut()
            .find(|f| f.id == "neb-r3101")
            .unwrap();
        let mut row = f.evidence[0].clone();
        row.qualifiers = vec!["x".repeat(1024); 16];
        f.evidence = vec![row; 64];
    });
    let s = selected(&c, &o);
    assert!(matches!(assess(&c, &s, &o), Err(Error::TooLarge)));
}
