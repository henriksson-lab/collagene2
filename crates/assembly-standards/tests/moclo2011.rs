//! Independent literal transcription of Weber2011 Figure 1B/C, not another toolkit.
use collagene_assembly_standards::{
    Error,
    assessment::{Certainty, FindingKind, Limits, Outcome, Part, Severity},
    grammar::GrammarIssue,
    profile::{
        BUILTIN_IDS, CatalogPin, MOCLO_2011_L0_ID, Profile, ProfileKind, RFC10_ID, RFC1000_L0_ID,
        Requirement,
    },
};
use collagene_enzyme_catalog::Catalog;
use collagene_model::Sequence;

fn installed() -> (Catalog, Profile) {
    let catalog = Catalog::builtin().unwrap();
    let profile = Profile::builtin(MOCLO_2011_L0_ID, &catalog).unwrap();
    (catalog, profile)
}
fn part(slot: &str, body: &str) -> Part {
    Part {
        id: slot.into(),
        slot: slot.into(),
        sequence: Sequence::dna(body).unwrap(),
        body_start: 0,
        body_end: body.len(),
    }
}
fn cytosolic() -> Vec<Part> {
    vec![
        part("promoter", "AACCGGTT"),
        part("five-prime-utr", "AAGGAAAA"),
        part("cytosolic-cds", "ATGAAATAA"),
        part("terminator", "TTACCGGT"),
    ]
}
fn signal_path() -> Vec<Part> {
    vec![
        part("promoter", "AACCGGTT"),
        part("five-prime-utr", "AAGGAAAA"),
        part("signal-peptide", "ATGAAGGGT"),
        part("cds-after-signal", "AAACCCTAA"),
        part("terminator", "TTACCGGT"),
    ]
}

#[test]
fn primary_figure_has_six_reference_oriented_pairs_and_explicit_body_only_scope() {
    let (catalog, profile) = installed();
    assert_eq!(profile.kind(), ProfileKind::Builtin);
    assert!(BUILTIN_IDS.contains(&MOCLO_2011_L0_ID));
    assert_eq!(profile.catalog_pin(), &CatalogPin::current(&catalog));
    let literal = [
        ("promoter", "GGAG", "TACT", 1),
        ("five-prime-utr", "TACT", "AATG", 1),
        ("signal-peptide", "AATG", "AGGT", 0),
        ("cds-after-signal", "AGGT", "GCTT", 0),
        ("cytosolic-cds", "AATG", "GCTT", 0),
        ("terminator", "GCTT", "CGCT", 1),
    ];
    let definition = profile.definition();
    assert_eq!(definition.grammar.slots.len(), literal.len());
    for (role, left, right, minimum) in literal {
        let slot = definition.grammar.slot(role).unwrap();
        assert_eq!(
            (slot.left_overhang.as_str(), slot.right_overhang.as_str()),
            (left, right)
        );
        assert_eq!(
            (slot.left_enzyme.as_str(), slot.right_enzyme.as_str()),
            ("bsai", "bsai")
        );
        assert_eq!(
            (slot.left_addition.as_str(), slot.right_addition.as_str()),
            ("", "")
        );
        assert_eq!(slot.minimum, minimum);
        assert_eq!(slot.maximum, Some(1));
        assert!(
            slot.body_rules.is_empty(),
            "fusion syntax is not translation/frame evidence"
        );
    }
    assert_eq!(
        &definition
            .grammar
            .slot("signal-peptide")
            .unwrap()
            .left_overhang[1..],
        "ATG"
    );
    assert_eq!(
        &definition
            .grammar
            .slot("signal-peptide")
            .unwrap()
            .right_overhang[1..],
        "GGT"
    );
    assert_ne!(
        definition
            .grammar
            .slot("cytosolic-cds")
            .unwrap()
            .right_overhang,
        "AAGC"
    );
    assert_eq!(
        definition.sources[0].upstream_sha256.as_deref(),
        Some("cbc7be2584380029a96188834bff83a13eec15c0d35e83dc830101ad25ade0ad")
    );
}

#[test]
fn both_literal_paths_are_complete_and_successor_hints_keep_the_real_branch() {
    let (_, profile) = installed();
    let limits = Limits::default();
    for input in [cytosolic(), signal_path()] {
        let before = input.clone();
        let assessment = profile.assess(&input, &limits).unwrap();
        assert_eq!(assessment.outcome, Outcome::Conforms);
        assert!(assessment.complete);
        assessment.validate_for(&profile, &input, &limits).unwrap();
        assert_eq!(
            input, before,
            "ATG/GGT overlap is never prepended or edited"
        );
        let bytes = postcard::to_allocvec(&assessment).unwrap();
        let decoded: collagene_assembly_standards::assessment::Assessment =
            postcard::from_bytes(&bytes).unwrap();
        decoded.validate_for(&profile, &input, &limits).unwrap();
    }
    let empty = profile.legal_successors(&[]).unwrap();
    assert_eq!(empty.allowed_slots, vec!["promoter"]);
    assert!(!empty.may_finish);
    let branch = profile
        .legal_successors(&["promoter".into(), "five-prime-utr".into()])
        .unwrap();
    assert_eq!(
        branch.allowed_slots,
        vec!["cytosolic-cds", "signal-peptide"]
    );
    assert!(!branch.may_finish);
    let secreted = profile
        .legal_successors(&[
            "promoter".into(),
            "five-prime-utr".into(),
            "signal-peptide".into(),
        ])
        .unwrap();
    assert_eq!(secreted.allowed_slots, vec!["cds-after-signal"]);
    assert!(!secreted.may_finish);
    assert!(
        profile
            .legal_successors(
                &cytosolic()
                    .iter()
                    .map(|p| p.slot.clone())
                    .collect::<Vec<_>>()
            )
            .unwrap()
            .may_finish
    );
}

#[test]
fn omitted_swapped_and_cross_branch_roles_cannot_masquerade_as_complete_paths() {
    let (_, profile) = installed();
    let mut swapped = cytosolic();
    swapped.swap(1, 2);
    let mut omitted = signal_path();
    omitted.remove(3);
    let mut mixed = signal_path();
    mixed[3] = part("cytosolic-cds", "ATGAAATAA");
    for input in [
        swapped,
        omitted,
        mixed,
        cytosolic()[1..].to_vec(),
        vec![part("cytosolic-cds", "ATGAAATAA")],
    ] {
        let result = profile.assess(&input, &Limits::default()).unwrap();
        assert_eq!(result.outcome, Outcome::Violations);
        assert!(
            result
                .findings
                .iter()
                .any(|finding| matches!(finding.kind, FindingKind::Grammar(_)))
        );
    }
    let issue = profile
        .definition()
        .grammar
        .violations(&[
            "promoter".into(),
            "five-prime-utr".into(),
            "signal-peptide".into(),
            "terminator".into(),
        ])
        .unwrap();
    assert!(
        issue
            .iter()
            .any(|i| matches!(i, GrammarIssue::InvalidTransition { .. }))
    );
}

#[test]
fn literal_internal_sites_are_scanned_both_strands_with_honest_hierarchy_severity() {
    let (_, profile) = installed();
    for (id, forward, reverse, severity) in [
        ("bsai", "GGTCTC", "GAGACC", Severity::Violation),
        ("bpii", "GAAGAC", "GTCTTC", Severity::Violation),
        ("esp3i", "CGTCTC", "GAGACG", Severity::Advisory),
    ] {
        let rule = profile
            .definition()
            .avoided_sites
            .iter()
            .find(|rule| rule.id == id)
            .unwrap();
        assert_eq!(rule.recognition, forward);
        assert_eq!(
            rule.requirement,
            if severity == Severity::Violation {
                Requirement::Required
            } else {
                Requirement::Advisory
            }
        );
        for (motif, reversed) in [(forward, false), (reverse, true)] {
            let mut input = cytosolic();
            input[0] = part("promoter", &format!("AA{motif}TT"));
            let result = profile.assess(&input, &Limits::default()).unwrap();
            assert_eq!(
                result.outcome,
                if severity == Severity::Violation {
                    Outcome::Violations
                } else {
                    Outcome::Conforms
                }
            );
            assert!(result.findings.iter().any(|finding| {
                finding.part_index == Some(0) && finding.severity == severity
                    && matches!(&finding.kind, FindingKind::AvoidedSite {rule_id,start:2,end:8,reverse,certainty:Certainty::Certain} if rule_id==id && *reverse==reversed)
            }));
        }
    }
}

#[test]
fn ambiguous_sites_partial_work_and_unclaimed_format_flanks_never_gain_conformance() {
    let (_, profile) = installed();
    let mut input = cytosolic();
    input[0] = part("promoter", "AAGGTCTNTT");
    let uncertain = profile.assess(&input, &Limits::default()).unwrap();
    assert_eq!(uncertain.outcome, Outcome::Indeterminate);
    assert!(uncertain.findings.iter().any(|finding| matches!(&finding.kind, FindingKind::AvoidedSite {rule_id,certainty:Certainty::Possible,..} if rule_id=="bsai")));
    let mut limits = Limits::default();
    limits.max_comparisons = 1;
    let limited = profile.assess(&cytosolic(), &limits).unwrap();
    assert!(!limited.complete);
    assert_ne!(limited.outcome, Outcome::Conforms);
    let mut input = cytosolic();
    input[0] = part("promoter", "GGTCTCAACCGGTTGAGACC");
    input[0].body_start = 6;
    input[0].body_end = 14;
    let flanked = profile.assess(&input, &Limits::default()).unwrap();
    assert_eq!(flanked.outcome, Outcome::Violations);
    assert!(
        flanked
            .findings
            .iter()
            .any(|finding| matches!(finding.kind, FindingKind::Addition { .. }))
    );
}

#[test]
fn reserved_ids_catalog_pins_and_retained_assessments_cannot_be_cross_rebound() {
    let (catalog, profile) = installed();
    let input = cytosolic();
    let limits = Limits::default();
    let result = profile.assess(&input, &limits).unwrap();
    for id in [RFC10_ID, RFC1000_L0_ID] {
        let other = Profile::builtin(id, &catalog).unwrap();
        assert_ne!(other.identity(), profile.identity());
        assert!(result.validate_for(&other, &input, &limits).is_err());
        let mut disguised = profile.definition().clone();
        disguised.id = id.into();
        assert!(Profile::new(disguised, CatalogPin::current(&catalog), &catalog).is_err());
        let mut inverse = other.definition().clone();
        inverse.id = MOCLO_2011_L0_ID.into();
        assert!(Profile::new(inverse, CatalogPin::current(&catalog), &catalog).is_err());
    }
    let mut changed = profile.definition().clone();
    changed.grammar.slots[2].right_overhang = "CGAC".into();
    assert!(Profile::new(changed.clone(), CatalogPin::current(&catalog), &catalog).is_err());
    changed.id = "user-authored-different-junction".into();
    let custom = Profile::new(changed, CatalogPin::current(&catalog), &catalog).unwrap();
    assert_eq!(custom.kind(), ProfileKind::Custom);
    assert!(result.validate_for(&custom, &input, &limits).is_err());
    let mut pin = CatalogPin::current(&catalog);
    pin.content_hash = "0".repeat(64);
    assert!(matches!(
        Profile::new(profile.definition().clone(), pin, &catalog),
        Err(Error::CatalogMismatch)
    ));
    assert_eq!(
        Profile::builtin("moclo-2011-level2-unreviewed", &catalog).err(),
        Some(Error::UnsupportedProfile)
    );
}
