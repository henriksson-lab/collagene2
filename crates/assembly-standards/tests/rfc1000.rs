//! Independent literal fixtures for the scoped iGEM Level-0 factual profile.
use collagene_assembly_standards::{
    Error,
    assessment::{Certainty, FindingKind, Limits, Outcome, Part, Severity},
    grammar::{BodyRule, GrammarIssue},
    profile::{
        BUILTIN_IDS, CatalogPin, Profile, ProfileKind, RFC10_ID, RFC1000_L0_ID, Requirement,
    },
};
use collagene_enzyme_catalog::Catalog;
use collagene_model::Sequence;

fn installed() -> (Catalog, Profile) {
    let catalog = Catalog::builtin().unwrap();
    let profile = Profile::builtin(RFC1000_L0_ID, &catalog).unwrap();
    (catalog, profile)
}
fn part(id: &str, slot: &str, body: &str) -> Part {
    Part {
        id: id.into(),
        slot: slot.into(),
        sequence: Sequence::dna(body).unwrap(),
        body_start: 0,
        body_end: body.len(),
    }
}
fn roles() -> Vec<Part> {
    vec![
        part("p", "promoter", "AACCGGTT"),
        part("u", "five-prime-utr", "AGGAGAAA"),
        part("c", "cds", "ATGAAATAA"),
        part("t", "terminator", "TTACCGGT"),
    ]
}

#[test]
fn independent_primary_table_has_reference_oriented_four_base_fusions_and_no_invented_flanks() {
    let (_, p) = installed();
    assert_eq!(p.kind(), ProfileKind::Builtin);
    assert!(BUILTIN_IDS.contains(&RFC1000_L0_ID));
    let literal = [
        ("promoter", "GGAG", "TACT"),
        ("five-prime-utr", "TACT", "AATG"),
        ("cds", "AATG", "GCTT"),
        ("terminator", "GCTT", "CGCT"),
    ];
    assert_eq!(p.definition().grammar.slots.len(), literal.len());
    for (role, left, right) in literal {
        let slot = p.definition().grammar.slot(role).unwrap();
        assert_eq!((&*slot.left_overhang, &*slot.right_overhang), (left, right));
        assert_eq!((&*slot.left_enzyme, &*slot.right_enzyme), ("bsai", "bsai"));
        assert_eq!((&*slot.left_addition, &*slot.right_addition), ("", ""));
        assert_eq!(slot.minimum, 0);
        assert_eq!(slot.maximum, Some(1));
    }
    assert_eq!(
        p.definition().grammar.slot("cds").unwrap().right_overhang,
        "GCTT"
    );
    // A reverse-complemented right physical strand label would be AAGC;
    // that is deliberately not the source table's reference-oriented label.
    assert_ne!(
        p.definition().grammar.slot("cds").unwrap().right_overhang,
        "AAGC"
    );
}

#[test]
fn declared_body_retains_single_start_codon_overlapping_the_fusion_label() {
    let (_, p) = installed();
    let cds = p.definition().grammar.slot("cds").unwrap();
    assert_eq!(&cds.left_overhang[1..], "ATG");
    assert_eq!(cds.body_rules, vec![BodyRule::RequiredPrefix("ATG".into())]);
    let input = part("coding", "cds", "ATGAAATAA");
    assert_eq!(input.sequence.len(), 9);
    let assessment = p.assess(&[input.clone()], &Limits::default()).unwrap();
    assert_eq!(assessment.outcome, Outcome::Conforms);
    assessment
        .validate_for(&p, &[input], &Limits::default())
        .unwrap();
    let duplicated_fusion = p
        .assess(&[part("wrong", "cds", "AATGATGAAATAA")], &Limits::default())
        .unwrap();
    assert_eq!(duplicated_fusion.outcome, Outcome::Violations);
    assert!(duplicated_fusion.findings.iter().any(|f| matches!(&f.kind, FindingKind::BodyRule {rule:BodyRule::RequiredPrefix(s),..} if s=="ATG")));
    let ambiguous = p
        .assess(&[part("unknown", "cds", "NTGAAATAA")], &Limits::default())
        .unwrap();
    assert_eq!(ambiguous.outcome, Outcome::Indeterminate);
}

#[test]
fn exact_source_site_rules_are_required_in_both_orientations_and_extra_sites_are_not_invented() {
    let (_, p) = installed();
    assert_eq!(p.definition().avoided_sites.len(), 2);
    for (id, motif) in [("bsai", "GGTCTC"), ("sapi", "GCTCTTC")] {
        let rule = p
            .definition()
            .avoided_sites
            .iter()
            .find(|s| s.id == id)
            .unwrap();
        assert_eq!(rule.recognition, motif);
        assert_eq!(rule.requirement, Requirement::Required);
    }
    for (id, literal, reverse) in [
        ("bsai", "AAGGTCTCTT", false),
        ("bsai", "AAGAGACCTT", true),
        ("sapi", "AAGCTCTTCTT", false),
        ("sapi", "AAGAAGAGCTT", true),
    ] {
        let report = p
            .assess(&[part("site", "promoter", literal)], &Limits::default())
            .unwrap();
        assert_eq!(report.outcome, Outcome::Violations);
        assert!(report.findings.iter().any(|f| f.severity==Severity::Violation && matches!(&f.kind, FindingKind::AvoidedSite {rule_id,reverse:orientation,certainty:Certainty::Certain,..} if rule_id==id && *orientation==reverse)));
    }
    let extra = p
        .assess(
            &[part("extra", "promoter", "AAGAATTCCGTCTCTT")],
            &Limits::default(),
        )
        .unwrap();
    assert_eq!(extra.outcome, Outcome::Conforms);
    assert!(extra.findings.is_empty()); // EcoRI and BsmBI are not this profile's rules.
    let uncertain = p
        .assess(
            &[part("maybe", "promoter", "AAGGTCTNTT")],
            &Limits::default(),
        )
        .unwrap();
    assert_eq!(uncertain.outcome, Outcome::Indeterminate);
}

#[test]
fn ordered_role_chain_and_partial_subpaths_have_explicitly_local_grammar() {
    let (_, p) = installed();
    let inputs = roles();
    assert_eq!(
        p.assess(&inputs, &Limits::default()).unwrap().outcome,
        Outcome::Conforms
    );
    assert_eq!(
        p.legal_successors(&["promoter".into()])
            .unwrap()
            .allowed_slots,
        ["five-prime-utr"]
    );
    assert_eq!(
        p.legal_successors(&["promoter".into(), "five-prime-utr".into()])
            .unwrap()
            .allowed_slots,
        ["cds"]
    );
    assert!(
        p.legal_successors(&["terminator".into()])
            .unwrap()
            .may_finish
    );
    assert!(
        p.legal_successors(&["terminator".into()])
            .unwrap()
            .allowed_slots
            .is_empty()
    );
    assert!(p.legal_successors(&["cds".into()]).unwrap().may_finish);
    assert!(!p.legal_successors(&[]).unwrap().may_finish);
    assert_eq!(
        p.assess(&inputs[1..3], &Limits::default()).unwrap().outcome,
        Outcome::Conforms
    );
    let mut wrong = inputs.clone();
    wrong.swap(1, 2);
    assert_eq!(
        p.assess(&wrong, &Limits::default()).unwrap().outcome,
        Outcome::Violations
    );
    assert!(
        p.legal_successors(&["cds".into(), "promoter".into()])
            .unwrap()
            .prefix_issues
            .iter()
            .any(|issue| matches!(issue, GrammarIssue::InvalidTransition { .. }))
    );
    assert_eq!(
        p.assess(&[part("tag", "n-tag", "ATGAAATAA")], &Limits::default())
            .unwrap()
            .outcome,
        Outcome::Violations
    );
}

#[test]
fn documented_body_inputs_refuse_unclaimed_vector_or_fusion_additions() {
    let (_, p) = installed();
    let mut formatted = part("with-flanks", "promoter", "GGTCTCAGGAGAACCGGTTTACTAGAGACC");
    formatted.body_start = 11;
    formatted.body_end = 19;
    let report = p.assess(&[formatted], &Limits::default()).unwrap();
    assert_eq!(report.outcome, Outcome::Violations);
}

#[test]
fn reserved_profile_ids_and_assessments_cannot_be_cross_rebound() {
    let (catalog, p) = installed();
    let rfc10 = Profile::builtin(RFC10_ID, &catalog).unwrap();
    assert_ne!(p.identity(), rfc10.identity());
    for id in [RFC10_ID, RFC1000_L0_ID] {
        let mut definition = if id == RFC10_ID {
            p.definition().clone()
        } else {
            rfc10.definition().clone()
        };
        definition.id = id.into();
        assert!(Profile::new(definition, CatalogPin::current(&catalog), &catalog).is_err());
    }
    let mut definition = p.definition().clone();
    definition.grammar.slots[2].left_overhang = "CATG".into();
    assert!(Profile::new(definition.clone(), CatalogPin::current(&catalog), &catalog).is_err());
    definition.id = "user-defined-different-fusion".into();
    let custom = Profile::new(definition, CatalogPin::current(&catalog), &catalog).unwrap();
    assert_eq!(custom.kind(), ProfileKind::Custom);
    assert_ne!(p.identity(), custom.identity());
    let inputs = roles();
    let limits = Limits::default();
    let report = p.assess(&inputs, &limits).unwrap();
    assert!(report.validate_for(&rfc10, &inputs, &limits).is_err());
    assert!(report.validate_for(&custom, &inputs, &limits).is_err());
    let mut changed = inputs;
    changed[2].sequence = Sequence::dna("ATGCCCTAA").unwrap();
    assert!(report.validate_for(&p, &changed, &limits).is_err());
    assert_eq!(
        Profile::builtin("igem-rfc1000-level1-unreviewed", &catalog).err(),
        Some(Error::UnsupportedProfile)
    );
}
