//! Source-authored tests, deliberately not compiled/executed before registration.
use collagene_assembly_standards::{
    Error,
    assessment::{BudgetHit, Certainty, FindingKind, Limits, Outcome, Part, Severity},
    grammar::{BodyRule, GrammarIssue},
    profile::{CatalogPin, Definition, Profile, ProfileKind, RFC10_ID},
};
use collagene_enzyme_catalog::Catalog;
use collagene_model::Sequence;

fn profile() -> (Catalog, Profile) {
    let catalog = Catalog::builtin().unwrap();
    let profile = Profile::builtin(RFC10_ID, &catalog).unwrap();
    (catalog, profile)
}
fn part(id: &str, slot: &str, body: &str) -> Part {
    // Independent literal RFC10 flank fixtures, not read from the profile under
    // test. Coding-prefix overlap uses the A of the body's mandatory ATG.
    let left = if slot == "coding" {
        "GAATTCGCGGCCGCTTCTAG"
    } else {
        "GAATTCGCGGCCGCTTCTAGAG"
    };
    let right = "TACTAGTAGCGGCCGCTGCAG";
    Part {
        id: id.into(),
        slot: slot.into(),
        sequence: Sequence::dna(&format!("{left}{body}{right}")).unwrap(),
        body_start: left.len(),
        body_end: left.len() + body.len(),
    }
}

#[test]
fn source_profile_distinguishes_coding_prefix_and_mandatory_noti_from_advisory_sites() {
    let (_, p) = profile();
    assert_eq!(p.kind(), ProfileKind::Builtin);
    assert_eq!(p.definition().grammar.slots.len(), 2);
    assert_eq!(
        p.definition()
            .grammar
            .slot("noncoding")
            .unwrap()
            .left_addition
            .len(),
        22
    );
    assert_eq!(
        p.definition()
            .grammar
            .slot("coding")
            .unwrap()
            .left_addition
            .len(),
        20
    );
    assert_eq!(
        p.definition()
            .grammar
            .slot("coding")
            .unwrap()
            .right_addition
            .len(),
        21
    );
    let plain = p
        .assess(
            &[part("plain", "noncoding", "AACCGGTT")],
            &Limits::default(),
        )
        .unwrap();
    assert_eq!(plain.outcome, Outcome::Conforms);
    assert!(plain.complete && plain.findings.is_empty());
    let noti = p
        .assess(
            &[part("noti", "noncoding", "AAGCGGCCGCTT")],
            &Limits::default(),
        )
        .unwrap();
    assert_eq!(noti.outcome, Outcome::Violations);
    assert!(noti.findings.iter().any(|f| f.severity == Severity::Violation && matches!(&f.kind, FindingKind::AvoidedSite { rule_id, start: 24, end: 32, certainty: Certainty::Certain, .. } if rule_id == "noti")));
    let optional = p
        .assess(
            &[part("optional", "noncoding", "AACAGCTGTT")],
            &Limits::default(),
        )
        .unwrap();
    assert_eq!(optional.outcome, Outcome::Conforms);
    assert!(optional.findings.iter().any(|f| f.severity == Severity::Advisory && matches!(&f.kind, FindingKind::AvoidedSite { rule_id, .. } if rule_id == "pvuii")));
}

#[test]
fn coding_start_is_required_but_preferred_stop_is_advisory_not_a_functional_certificate() {
    let (_, p) = profile();
    let valid = p
        .assess(
            &[part("coding", "coding", "ATGAAATAATAA")],
            &Limits::default(),
        )
        .unwrap();
    assert_eq!(valid.outcome, Outcome::Conforms);
    assert_eq!(valid.advisories_observed, 0);
    let no_start = p
        .assess(
            &[part("coding", "coding", "GTGAAATAATAA")],
            &Limits::default(),
        )
        .unwrap();
    assert_eq!(no_start.outcome, Outcome::Violations);
    assert!(no_start.findings.iter().any(|f| matches!(&f.kind, FindingKind::BodyRule { rule: BodyRule::RequiredPrefix(s), .. } if s == "ATG")));
    let alternate_stop = p
        .assess(&[part("coding", "coding", "ATGAAATGA")], &Limits::default())
        .unwrap();
    assert_eq!(alternate_stop.outcome, Outcome::Conforms);
    assert!(alternate_stop.advisories_observed > 0);
    let unknown_start = p
        .assess(
            &[part("coding", "coding", "NTGAAATAATAA")],
            &Limits::default(),
        )
        .unwrap();
    assert_eq!(unknown_start.outcome, Outcome::Indeterminate);
}

#[test]
fn forbidden_sites_are_scanned_on_both_orientations_and_possible_is_not_certain() {
    let (_, p) = profile();
    let reverse_sapi = p
        .assess(
            &[part("reverse", "noncoding", "AAGAAGAGCTT")],
            &Limits::default(),
        )
        .unwrap();
    assert!(reverse_sapi.findings.iter().any(|f| matches!(&f.kind, FindingKind::AvoidedSite { rule_id, reverse: true, .. } if rule_id == "sapi")));
    let uncertain = p
        .assess(
            &[part("ambiguous", "noncoding", "AAGAATTNTT")],
            &Limits::default(),
        )
        .unwrap();
    assert_eq!(uncertain.outcome, Outcome::Indeterminate);
    assert!(uncertain.findings.iter().any(|f| f.severity == Severity::Indeterminate && matches!(&f.kind, FindingKind::AvoidedSite { rule_id, certainty: Certainty::Possible, .. } if rule_id == "ecori")));
    // Expected format sites in the declared additions are not illegal internal
    // sites. A malformed suffix is independently rejected as a format error.
    let mut malformed = part("flank", "noncoding", "AACCGGTT");
    malformed.sequence =
        Sequence::dna(&malformed.sequence.as_str().replace("CTGCAG", "CTGCAA")).unwrap();
    assert_eq!(
        p.assess(&[malformed], &Limits::default()).unwrap().outcome,
        Outcome::Violations
    );
}

#[test]
fn rfc10_has_no_invented_functional_order_but_custom_grammars_enforce_declared_cardinality() {
    let (catalog, p) = profile();
    for prefix in [
        vec![],
        vec!["coding".into()],
        vec!["noncoding".into(), "coding".into()],
    ] {
        assert_eq!(
            p.legal_successors(&prefix).unwrap().allowed_slots,
            ["coding", "noncoding"]
        );
    }
    assert!(!p.legal_successors(&[]).unwrap().may_finish);
    assert!(p.legal_successors(&["coding".into()]).unwrap().may_finish);
    let mut definition = p.definition().clone();
    definition.id = "declared-custom-two-slot-order".into();
    definition.grammar.starts = vec!["noncoding".into()];
    definition.grammar.ends = vec!["coding".into()];
    definition
        .grammar
        .transitions
        .retain(|e| e.from == "noncoding" && e.to == "coding");
    definition.grammar.slot("coding").unwrap();
    for slot in &mut definition.grammar.slots {
        slot.minimum = 1;
        slot.maximum = Some(1);
    }
    let custom = Profile::new(definition, CatalogPin::current(&catalog), &catalog).unwrap();
    assert_eq!(custom.kind(), ProfileKind::Custom);
    assert_eq!(
        custom.legal_successors(&[]).unwrap().allowed_slots,
        ["noncoding"]
    );
    assert_eq!(
        custom
            .legal_successors(&["noncoding".into()])
            .unwrap()
            .allowed_slots,
        ["coding"]
    );
    assert!(
        custom
            .legal_successors(&["noncoding".into(), "coding".into()])
            .unwrap()
            .may_finish
    );
    assert!(
        custom
            .legal_successors(&["coding".into()])
            .unwrap()
            .prefix_issues
            .iter()
            .any(|i| matches!(i, GrammarIssue::InvalidStart { .. }))
    );
}

#[test]
fn built_in_identity_requires_current_exact_facts_and_catalog_pin() {
    let (catalog, p) = profile();
    let mut definition = p.definition().clone();
    definition.grammar.slots[0].left_addition.push('A');
    assert!(Profile::new(definition.clone(), CatalogPin::current(&catalog), &catalog).is_err());
    definition.id = "user-declared-change".into();
    assert_eq!(
        Profile::new(definition, CatalogPin::current(&catalog), &catalog)
            .unwrap()
            .kind(),
        ProfileKind::Custom
    );
    let mut pin = CatalogPin::current(&catalog);
    pin.content_hash = "0".repeat(64);
    assert_eq!(
        Profile::new(p.definition().clone(), pin, &catalog).err(),
        Some(Error::CatalogMismatch)
    );
    let mut geometry = p.definition().clone();
    geometry.enzymes[0].geometry.top_cut += 1;
    assert_eq!(
        Profile::new(geometry, CatalogPin::current(&catalog), &catalog).err(),
        Some(Error::CatalogMismatch)
    );
    assert_eq!(
        Profile::builtin("unverified-moclo", &catalog).err(),
        Some(Error::UnsupportedProfile)
    );
    let mut extra = serde_json::to_value(p.definition()).unwrap();
    extra["verified"] = true.into();
    assert!(
        Profile::from_json(
            &serde_json::to_vec(&extra).unwrap(),
            CatalogPin::current(&catalog),
            &catalog
        )
        .is_err()
    );
}

#[test]
fn interrupted_or_truncated_checks_never_establish_conformance() {
    let (_, p) = profile();
    let parts = [part("long", "noncoding", &"A".repeat(2000))];
    let mut limits = Limits::default();
    limits.max_comparisons = 1;
    let limited = p.assess(&parts, &limits).unwrap();
    assert!(!limited.complete);
    assert_eq!(limited.outcome, Outcome::Indeterminate);
    assert!(limited.budgets_hit.contains(&BudgetHit::Comparisons));
    let many = [part("many", "noncoding", &"GAATTC".repeat(20))];
    let mut limits = Limits::default();
    limits.max_findings = 1;
    let limited = p.assess(&many, &limits).unwrap();
    assert_eq!(limited.outcome, Outcome::Violations);
    assert!(!limited.complete && limited.budgets_hit.contains(&BudgetHit::Findings));
    let mut limits = Limits::default();
    limits.max_bases = 1;
    assert!(matches!(p.assess(&parts, &limits), Err(Error::Limit(_))));
    let parts = [part("report-limit", "noncoding", &"GAATTC".repeat(200))];
    let mut limits = Limits::default();
    limits.max_report_bytes = 4096;
    let clipped = p.assess(&parts, &limits).unwrap();
    assert_eq!(clipped.outcome, Outcome::Violations);
    assert!(!clipped.complete && clipped.budgets_hit.contains(&BudgetHit::ReportBytes));
    assert!(serde_json::to_vec(&clipped).unwrap().len() <= 4096);
}

#[test]
fn replay_and_binary_roundtrip_bind_scope_source_order_body_and_options() {
    let (_, p) = profile();
    let parts = [
        part("one", "coding", "ATGAAATAATAA"),
        part("two", "noncoding", "AACCGGTT"),
    ];
    let limits = Limits::default();
    let report = p.assess(&parts, &limits).unwrap();
    report.validate_for(&p, &parts, &limits).unwrap();
    let encoded = postcard::to_allocvec(&report).unwrap();
    let decoded: collagene_assembly_standards::assessment::Assessment =
        postcard::from_bytes(&encoded).unwrap();
    assert_eq!(decoded, report);
    let definition = postcard::to_allocvec(p.definition()).unwrap();
    let decoded: Definition = postcard::from_bytes(&definition).unwrap();
    assert_eq!(&decoded, p.definition());
    let mut changed = parts.clone();
    changed.swap(0, 1);
    assert!(report.validate_for(&p, &changed, &limits).is_err());
    changed[0].body_start += 1;
    assert!(report.validate_for(&p, &changed, &limits).is_err());
    let mut changed_limits = limits;
    changed_limits.max_comparisons -= 1;
    assert!(report.validate_for(&p, &parts, &changed_limits).is_err());
}
