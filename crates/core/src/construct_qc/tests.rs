use super::*;
use collagene_model::Topology;

fn feature(kind: &str, location: Location) -> Feature {
    Feature {
        kind: kind.into(),
        location,
        qualifiers: vec![],
        needs_review: false,
    }
}
fn fixture(dna: &str) -> (Molecule, Options) {
    let mut m = Molecule::dna("literal construct", dna, Topology::Linear).unwrap();
    m.features
        .push(feature("CDS", Location::range(0, dna.len())));
    let mut o = Options::default();
    declare(&m, &mut o, 0);
    (m, o)
}
fn declare(m: &Molecule, o: &mut Options, i: usize) {
    o.annotations.push(AnnotationAssumption {feature_index:i,source_revision:m.revision(),
        annotation_revision:annotation_revision(&m.features[i]).unwrap(),rationale:"Explicit CDS/component role for a literal structural fixture; no biological validation claimed".into()});
}
fn refresh(m: &Molecule, o: &mut Options) {
    for a in &mut o.annotations {
        a.source_revision = m.revision();
        a.annotation_revision = annotation_revision(&m.features[a.feature_index]).unwrap();
    }
}
fn verdict(r: &Report, rule: Rule, i: usize) -> Verdict {
    r.findings
        .iter()
        .find(|f| f.rule == rule && f.feature_indices == [i])
        .unwrap()
        .verdict
        .clone()
}

#[test]
fn complete_literal_cds_has_independent_structural_facts_not_an_expression_verdict() {
    let (mut m, mut o) = fixture("ATGAAATAA");
    m.features[0]
        .qualifiers
        .push(("translation".into(), Some("MK".into())));
    refresh(&m, &mut o);
    let r = assess(&m, &o).unwrap();
    for rule in [
        Rule::CdsFrame,
        Rule::CdsStart,
        Rule::CdsTerminalStop,
        Rule::CdsInternalStop,
        Rule::StoredTranslation,
    ] {
        assert_eq!(verdict(&r, rule, 0), Verdict::Pass);
    }
    assert!(
        r.findings
            .iter()
            .filter(|f| f.verdict == Verdict::Pass)
            .all(|f| f.under_declared_annotation_assumptions)
    );
    assert_eq!(r.statement, STATEMENT);
    assert_eq!(r.codes[0].descriptor.id, 1);
    assert_eq!(r.codes[0].ordinary_amino_acids.len(), 64);
    r.validate_for(&m, &o).unwrap();
}

#[test]
fn origin_spanning_cds_reuses_existing_biological_join_order() {
    let mut m = Molecule::dna("joined literal", "AAATAACCCATG", Topology::Circular).unwrap();
    m.features.push(feature(
        "CDS",
        Location::Join(vec![Location::range(9, 12), Location::range(0, 6)]),
    ));
    let mut o = Options::default();
    declare(&m, &mut o, 0);
    let r = assess(&m, &o).unwrap();
    for rule in [
        Rule::CdsFrame,
        Rule::CdsStart,
        Rule::CdsTerminalStop,
        Rule::CdsInternalStop,
    ] {
        assert_eq!(verdict(&r, rule, 0), Verdict::Pass);
    }
    assert_eq!(r.source_topology, Topology::Circular);
    assert_eq!(r.findings[0].coordinates, [m.features[0].location.clone()]);
}

#[test]
fn literal_frame_start_terminal_and_internal_stop_failures_are_distinct() {
    for (dna, rule, expected) in [
        ("ATGAAAATAA", Rule::CdsFrame, Verdict::Fail),
        ("AAAAAATAA", Rule::CdsStart, Verdict::Fail),
        ("ATGAAA", Rule::CdsTerminalStop, Verdict::Fail),
        ("ATGTAAGGGTAA", Rule::CdsInternalStop, Verdict::Fail),
        ("ATGAAATAN", Rule::CdsTerminalStop, Verdict::Indeterminate),
        (
            "ATGTANAAATAA",
            Rule::CdsInternalStop,
            Verdict::Indeterminate,
        ),
        ("NTGAAATAA", Rule::CdsStart, Verdict::Indeterminate),
    ] {
        let (m, o) = fixture(dna);
        let r = assess(&m, &o).unwrap();
        assert_eq!(verdict(&r, rule, 0), expected, "{dna}");
        assert!(
            r.findings
                .iter()
                .filter(|f| f.feature_indices == [0])
                .all(|f| f.coordinates == [Location::range(0, dna.len())])
        );
    }
}

#[test]
fn selected_and_annotated_genetic_codes_and_context_dependent_stops_are_pinned() {
    let (mut m, mut o) = fixture("ATGAGATAA");
    assert_eq!(
        verdict(&assess(&m, &o).unwrap(), Rule::CdsInternalStop, 0),
        Verdict::Pass
    );
    m.features[0]
        .qualifiers
        .push(("transl_table".into(), Some(" 2 ".into())));
    refresh(&m, &mut o);
    let r = assess(&m, &o).unwrap();
    assert_eq!(verdict(&r, Rule::CdsInternalStop, 0), Verdict::Fail); // AGA is stop in table 2
    assert!(r.codes.iter().any(|c| c.descriptor.id == 2));
    let (m, mut o) = fixture("ATGAAATGA");
    o.genetic_code = 28;
    assert_eq!(
        verdict(&assess(&m, &o).unwrap(), Rule::CdsTerminalStop, 0),
        Verdict::Indeterminate
    );
}

#[test]
fn missing_stale_unreviewed_partial_or_exception_annotations_never_become_passes() {
    let (m, mut o) = fixture("ATGAAATAA");
    o.annotations.clear();
    let r = assess(&m, &o).unwrap();
    assert_eq!(verdict(&r, Rule::CdsStart, 0), Verdict::Indeterminate);
    for (key, value) in [
        ("partial", None),
        ("transl_except", Some("(pos:4..6,aa:Sec)")),
    ] {
        let (mut m, mut o) = fixture("ATGAAATAA");
        m.features[0]
            .qualifiers
            .push((key.into(), value.map(str::to_string)));
        refresh(&m, &mut o);
        let r = assess(&m, &o).unwrap();
        assert_eq!(verdict(&r, Rule::CdsStart, 0), Verdict::Indeterminate);
    }
    let (mut m, mut o) = fixture("ATGAAATAA");
    m.features[0].needs_review = true;
    refresh(&m, &mut o);
    assert_eq!(
        verdict(&assess(&m, &o).unwrap(), Rule::CdsFrame, 0),
        Verdict::Indeterminate
    );
    let (mut m, mut o) = fixture("ATGAAATAA");
    m.features.push(feature("CDS", Location::range(0, 9)));
    refresh(&m, &mut o);
    declare(&m, &mut o, 1);
    o.annotations[0].source_revision = "0".repeat(64);
    let r = assess(&m, &o).unwrap();
    assert_eq!(verdict(&r, Rule::CdsStart, 0), Verdict::Indeterminate);
    assert_eq!(verdict(&r, Rule::CdsStart, 1), Verdict::Pass); // unrelated check still runs
    let mut no_annotation = m;
    no_annotation.features.clear();
    let r = assess(&no_annotation, &Options::default()).unwrap();
    assert!(
        r.findings
            .iter()
            .any(|f| f.rule == Rule::AnnotationBasis && f.verdict == Verdict::Indeterminate)
    );
    assert!(!r.findings.iter().any(|f| f.verdict == Verdict::Pass));
}

#[test]
fn stored_translation_disagreement_is_a_documentary_finding() {
    let (mut m, mut o) = fixture("ATGAAATAA");
    m.features[0]
        .qualifiers
        .push(("translation".into(), Some("MR".into())));
    refresh(&m, &mut o);
    let r = assess(&m, &o).unwrap();
    assert_eq!(
        verdict(&r, Rule::StoredTranslation, 0),
        Verdict::Informational
    );
    assert!(
        r.findings
            .iter()
            .any(|f| f.rule == Rule::StoredTranslation && f.reason.contains("differs literally"))
    );
}

#[test]
fn fusion_and_junction_checks_report_declared_register_without_edit_causality() {
    let (mut m, mut o) = fixture("ATGAAACCCGGGTAA");
    m.features
        .push(feature("misc_feature", Location::range(3, 9)));
    refresh(&m, &mut o);
    declare(&m, &mut o, 1);
    o.fusions.push(FusionClaim {
        cds_feature: 0,
        part_feature: 1,
        expected_offset_mod3: 0,
        rationale: "Declared full-codon component boundaries".into(),
    });
    o.junctions.push(JunctionClaim {
        cds_feature: 0,
        boundary: 4,
        expected_offset_mod3: 0,
        rationale: "Declared current boundary".into(),
        history_document_id: Some("document".into()),
        history_event_id: Some("event".into()),
    });
    let r = assess(&m, &o).unwrap();
    assert!(
        r.findings
            .iter()
            .any(|f| f.rule == Rule::FusionRegister && f.verdict == Verdict::Pass)
    );
    assert_eq!(
        verdict(&r, Rule::DeclaredJunctionRegister, 0),
        Verdict::Fail
    );
    assert!(r.findings.iter().any(|f|f.rule==Rule::DeclaredJunctionRegister && f.reason.contains("does not show")));
    m.features[1].location = Location::range(4, 9);
    refresh(&m, &mut o);
    assert!(
        assess(&m, &o)
            .unwrap()
            .findings
            .iter()
            .any(|f| f.rule == Rule::FusionRegister && f.verdict == Verdict::Fail)
    );
    let mut reverse =
        Molecule::dna("reverse literal", "TTACCCGGGTTTCAT", Topology::Linear).unwrap();
    reverse.features = vec![
        feature(
            "CDS",
            Location::Complement(Box::new(Location::range(0, 15))),
        ),
        feature(
            "misc_feature",
            Location::Complement(Box::new(Location::range(6, 12))),
        ),
    ];
    let mut ro = Options::default();
    declare(&reverse, &mut ro, 0);
    declare(&reverse, &mut ro, 1);
    ro.fusions = o.fusions.clone();
    assert!(
        assess(&reverse, &ro)
            .unwrap()
            .findings
            .iter()
            .any(|f| f.rule == Rule::FusionRegister && f.verdict == Verdict::Pass)
    );
}

#[test]
fn pinned_reference_identity_is_separate_from_role_and_provenance_truth() {
    let (m, mut o) = fixture("ATGAAATAA");
    let sequence = Sequence::dna("ATGAAATAA").unwrap();
    o.references.push(ReferenceElement {
        id: "literal-reference".into(),
        sequence_revision: reference_revision(&sequence).unwrap(),
        sequence,
        provenance: Provenance {
            source_url: "https://example.test/reference".into(),
            accession: "literal-fixture-1".into(),
            retrieved_on: "2026-11-28".into(),
            license: "MIT".into(),
        },
    });
    o.reference_claims.push(ReferenceClaim {
        feature_index: 0,
        reference_id: "literal-reference".into(),
    });
    assert_eq!(
        verdict(&assess(&m, &o).unwrap(), Rule::ReferenceIntegrity, 0),
        Verdict::Pass
    );
    o.references[0].sequence = Sequence::dna("ATGCCCTAA").unwrap();
    assert_eq!(
        verdict(&assess(&m, &o).unwrap(), Rule::ReferenceIntegrity, 0),
        Verdict::Indeterminate
    );
    o.references[0].sequence_revision = reference_revision(&o.references[0].sequence).unwrap();
    assert_eq!(
        verdict(&assess(&m, &o).unwrap(), Rule::ReferenceIntegrity, 0),
        Verdict::Fail
    );
    o.references.clear();
    let r = assess(&m, &o).unwrap();
    assert_eq!(
        verdict(&r, Rule::ReferenceIntegrity, 0),
        Verdict::Indeterminate
    );
    assert_eq!(verdict(&r, Rule::CdsStart, 0), Verdict::Pass);
}

#[test]
fn bindings_tampering_and_repeated_location_amplification_are_bounded() {
    let (m, o) = fixture("ATGAAATAA");
    let r = assess(&m, &o).unwrap();
    let mut other = m.clone();
    other.name = "changed annotation/source envelope".into();
    assert!(r.validate_binding(&revision(&other, &o).unwrap()).is_err());
    let mut tampered = r;
    tampered.findings[0].reason = "fabricated documentary content".into();
    assert!(tampered.validate().is_err());
    tampered.report_revision = tampered.content_revision().unwrap();
    tampered.validate().unwrap();
    assert!(tampered.validate_for(&m, &o).is_err());
    let mut excessive = o.clone();
    excessive.annotations = vec![o.annotations[0].clone(); 257];
    assert!(matches!(revision(&m, &excessive), Err(Error::TooLarge)));
    let mut repeated =
        Molecule::dna("repeated CDS path", &"ATG".repeat(3000), Topology::Linear).unwrap();
    repeated.features.push(feature(
        "CDS",
        Location::Join(vec![Location::range(0, 9000); 200]),
    ));
    assert!(matches!(
        assess(&repeated, &Options::default()),
        Err(Error::TooLarge)
    ));
}
