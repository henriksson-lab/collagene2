use super::*;
use element_chain::{ChainElement, ElementChain, ElementRole};
use required_parts::{RequiredOrientation, RequiredPart};

fn feature(kind: &str, start: usize, end: usize, reverse: bool) -> Feature {
    let location = Location::range(start, end);
    Feature {
        kind: kind.into(),
        location: if reverse {
            Location::Complement(Box::new(location))
        } else {
            location
        },
        qualifiers: vec![],
        needs_review: false,
    }
}
fn reference(options: &mut Options, id: &str, dna: &str) {
    let sequence = Sequence::dna(dna).unwrap();
    options.references.push(ReferenceElement {
        id: id.into(),
        sequence_revision: reference_revision(&sequence).unwrap(),
        sequence,
        provenance: Provenance {
            source_url: "https://example.test/literal-fixtures".into(),
            accession: id.into(),
            retrieved_on: "2026-11-29".into(),
            license: "MIT".into(),
        },
    });
}
fn assume_all(m: &Molecule, options: &mut Options) {
    options.annotations = m
        .features
        .iter()
        .enumerate()
        .map(|(i, f)| AnnotationAssumption {
            feature_index: i,
            source_revision: m.revision(),
            annotation_revision: annotation_revision(f).unwrap(),
            rationale:
                "Literal fixture roles are declared assumptions, not biological identity evidence."
                    .into(),
        })
        .collect();
}
fn chain_fixture(reverse: bool) -> (Molecule, Options) {
    let dna = if reverse {
        "CCTTATTTCATGG"
    } else {
        "CCATGAAATAAGG"
    };
    let mut m = Molecule::dna("declared linear chain", dna, Topology::Linear).unwrap();
    m.features = if reverse {
        vec![
            feature("CDS", 2, 11, true),
            feature("misc_feature", 11, 13, true),
            feature("misc_feature", 0, 2, true),
        ]
    } else {
        vec![
            feature("CDS", 2, 11, false),
            feature("misc_feature", 0, 2, false),
            feature("misc_feature", 11, 13, false),
        ]
    };
    let mut o = Options::default();
    assume_all(&m, &mut o);
    reference(&mut o, "upstream", "CC");
    reference(&mut o, "downstream", "GG");
    o.element_chains.push(ElementChain {
        id: "literal-chain".into(),
        cds_feature: 0,
        upstream: Some(ChainElement {
            feature_index: 1,
            reference_id: Some("upstream".into()),
            role: ElementRole::Promoter,
        }),
        downstream: Some(ChainElement {
            feature_index: 2,
            reference_id: Some("downstream".into()),
            role: ElementRole::Terminator,
        }),
        max_upstream_distance: Some(0),
        max_downstream_distance: Some(0),
    });
    (m, o)
}
fn requirement(id: &str, reference_id: &str) -> RequiredPart {
    RequiredPart {
        id: id.into(),
        role: "declared marker role (not functionally verified)".into(),
        reference_id: reference_id.into(),
        minimum_copies: 1,
        maximum_copies: Some(1),
        orientation: RequiredOrientation::Either,
    }
}

#[test]
fn exact_flanking_chains_measure_zero_distance_on_both_strands() {
    for reverse in [false, true] {
        let (m, o) = chain_fixture(reverse);
        let r = assess(&m, &o).unwrap();
        let chain = &r.element_chains[0];
        assert_eq!(chain.verdict, Verdict::Pass);
        assert_eq!(chain.upstream.distance_bases, Some(0));
        assert_eq!(chain.downstream.distance_bases, Some(0));
        assert_eq!(chain.upstream.reverse, Some(reverse));
        assert!(
            r.findings
                .iter()
                .any(|f| f.rule == Rule::ElementChain && f.coordinates.len() == 3)
        );
        r.validate_for(&m, &o).unwrap();
    }
}

#[test]
fn incomplete_chain_premises_remain_indeterminate_and_unrelated_cds_checks_run() {
    let (m, mut o) = chain_fixture(false);
    o.element_chains[0].downstream = None;
    let r = assess(&m, &o).unwrap();
    assert_eq!(r.element_chains[0].verdict, Verdict::Indeterminate);
    assert_eq!(r.element_chains[0].upstream.verdict, Verdict::Pass);
    assert_eq!(r.element_chains[0].downstream.distance_bases, None);
    assert!(
        r.findings
            .iter()
            .any(|f| f.rule == Rule::CdsStart && f.verdict == Verdict::Pass)
    );
    let (mut m, mut o) = chain_fixture(false);
    m.features[1].needs_review = true;
    assume_all(&m, &mut o);
    assert_eq!(
        assess(&m, &o).unwrap().element_chains[0].upstream.verdict,
        Verdict::Indeterminate
    );
    m.features[1].needs_review = false;
    m.features[0].location = Location::Join(vec![Location::range(2, 5), Location::range(5, 11)]);
    assume_all(&m, &mut o);
    assert_eq!(
        assess(&m, &o).unwrap().element_chains[0].upstream.verdict,
        Verdict::Indeterminate
    );
    let (m, mut o) = chain_fixture(false);
    o.references[0].sequence_revision = "0".repeat(64);
    assert_eq!(
        assess(&m, &o).unwrap().element_chains[0].upstream.verdict,
        Verdict::Indeterminate
    );
}

#[test]
fn orientation_and_explicit_spacing_failures_are_structural() {
    let (mut m, mut o) = chain_fixture(false);
    m.features[1].location = Location::Complement(Box::new(Location::range(0, 2)));
    o.references[0].sequence = Sequence::dna("GG").unwrap();
    o.references[0].sequence_revision = reference_revision(&o.references[0].sequence).unwrap();
    assume_all(&m, &mut o);
    assert_eq!(
        assess(&m, &o).unwrap().element_chains[0].upstream.verdict,
        Verdict::Fail
    );
    let (mut m, mut o) = chain_fixture(false);
    m.features[1].location = Location::range(0, 1);
    o.references[0].sequence = Sequence::dna("C").unwrap();
    o.references[0].sequence_revision = reference_revision(&o.references[0].sequence).unwrap();
    assume_all(&m, &mut o);
    let r = assess(&m, &o).unwrap();
    assert_eq!(r.element_chains[0].upstream.distance_bases, Some(1));
    assert_eq!(r.element_chains[0].verdict, Verdict::Fail);
}

#[test]
fn circular_chain_requires_a_nonreusing_single_turn() {
    let (mut m, mut o) = chain_fixture(false);
    m.topology = Topology::Circular;
    assume_all(&m, &mut o);
    assert_eq!(
        assess(&m, &o).unwrap().element_chains[0].verdict,
        Verdict::Pass
    );
    o.element_chains[0].upstream.as_mut().unwrap().feature_index = 2;
    o.element_chains[0].upstream.as_mut().unwrap().reference_id = Some("downstream".into());
    o.element_chains[0]
        .downstream
        .as_mut()
        .unwrap()
        .feature_index = 1;
    o.element_chains[0]
        .downstream
        .as_mut()
        .unwrap()
        .reference_id = Some("upstream".into());
    o.element_chains[0].max_upstream_distance = None;
    o.element_chains[0].max_downstream_distance = None;
    let r = assess(&m, &o).unwrap();
    assert_eq!(r.element_chains[0].upstream.verdict, Verdict::Pass);
    assert_eq!(r.element_chains[0].downstream.verdict, Verdict::Pass);
    assert_eq!(r.element_chains[0].verdict, Verdict::Fail);
    assert!(r.element_chains[0].reason.contains("one source turn"));
}

#[test]
fn required_intervals_count_overlapping_hits_and_only_complete_searches_prove_absence() {
    let m = Molecule::dna("literal overlapping intervals", "AAAAA", Topology::Linear).unwrap();
    let mut o = Options::default();
    reference(&mut o, "motif", "AAA");
    o.required_parts.push(requirement("marker", "motif"));
    let r = assess(&m, &o).unwrap();
    assert_eq!(r.required_parts[0].certain_copies, 3);
    assert_eq!(r.required_parts[0].possible_copies, 3);
    assert_eq!(r.required_parts[0].verdict, Verdict::Fail);
    assert!(r.required_part_searches[0].complete);
    assert!(r.findings.iter().any(|f| f.rule == Rule::RequiredPart
        && f.feature_indices.is_empty()
        && !f.coordinates.is_empty()));
    let absent = Molecule::dna("literal absence", "CCCCC", Topology::Linear).unwrap();
    let r = assess(&absent, &o).unwrap();
    assert_eq!(r.required_parts[0].certain_copies, 0);
    assert_eq!(r.required_parts[0].verdict, Verdict::Fail);
    o.required_part_limits.max_work = 1;
    let r = assess(&absent, &o).unwrap();
    assert!(!r.required_part_searches[0].complete);
    assert_eq!(r.required_parts[0].verdict, Verdict::Indeterminate);
}

#[test]
fn ambiguous_source_and_palindromic_orientation_are_never_guessed() {
    let m = Molecule::dna("possible motif", "ATNAAA", Topology::Linear).unwrap();
    let mut o = Options::default();
    reference(&mut o, "motif", "ATG");
    o.required_parts.push(requirement("marker", "motif"));
    let r = assess(&m, &o).unwrap();
    assert!(r.required_part_searches[0].complete);
    assert_eq!(
        (
            r.required_parts[0].certain_copies,
            r.required_parts[0].possible_copies
        ),
        (0, 1)
    );
    assert_eq!(r.required_parts[0].verdict, Verdict::Indeterminate);
    let m = Molecule::dna("one palindrome", "GAATTCAAAA", Topology::Linear).unwrap();
    o.references.clear();
    reference(&mut o, "motif", "GAATTC");
    let r = assess(&m, &o).unwrap();
    assert_eq!(r.required_parts[0].certain_copies, 1);
    assert_eq!(r.required_parts[0].verdict, Verdict::Pass);
    assert_eq!(r.required_part_searches[0].matches[0].reverse, None);
    o.required_parts[0].orientation = RequiredOrientation::Reverse;
    assert_eq!(
        assess(&m, &o).unwrap().required_parts[0].verdict,
        Verdict::Indeterminate
    );
}

#[test]
fn origin_spanning_reference_matches_retain_original_selectable_spans() {
    let m = Molecule::dna("origin crossing literal", "TGCCCA", Topology::Circular).unwrap();
    let mut o = Options::default();
    reference(&mut o, "motif", "ATG");
    o.required_parts.push(requirement("marker", "motif"));
    o.required_parts[0].orientation = RequiredOrientation::Forward;
    let r = assess(&m, &o).unwrap();
    assert_eq!(r.required_parts[0].verdict, Verdict::Pass);
    let hit = r.required_part_searches[0]
        .matches
        .iter()
        .find(|h| h.reverse == Some(false))
        .unwrap();
    assert_eq!(
        hit.coordinates,
        [Location::range(5, 6), Location::range(0, 2)]
    );
}

#[test]
fn shared_query_work_and_match_caps_remain_explicit_and_bound() {
    let m = Molecule::dna("two literal copies", "ATGAAAATG", Topology::Linear).unwrap();
    let mut o = Options::default();
    reference(&mut o, "motif", "ATG");
    o.required_parts.push(requirement("role-one", "motif"));
    o.required_parts.push(requirement("role-two", "motif"));
    let r = assess(&m, &o).unwrap();
    assert_eq!(r.required_parts.len(), 2);
    assert_eq!(r.required_part_searches.len(), 1);
    assert_eq!(
        r.required_parts[0].search_index,
        r.required_parts[1].search_index
    );
    assert_eq!(
        r.required_part_search_work,
        r.required_part_searches[0].work
    );
    o.required_part_limits.max_matches = 1;
    o.required_part_limits.max_matches_per_reference = 1;
    reference(&mut o, "other", "AAA");
    o.required_parts.push(requirement("role-three", "other"));
    let r = assess(&m, &o).unwrap();
    assert_eq!(r.required_part_searches[0].matches.len(), 1);
    assert!(!r.required_part_searches[0].complete);
    assert!(r.required_part_searches[1].matches.is_empty());
    assert!(!r.required_part_searches[1].complete);
    assert!(
        r.required_parts
            .iter()
            .all(|p| p.verdict == Verdict::Indeterminate)
    );
    r.validate_for(&m, &o).unwrap();
    let mut tampered = r;
    tampered.required_part_searches[0].complete = true;
    assert!(tampered.validate().is_err());
    let mut old = o;
    old.schema_version = 1;
    assert!(revision(&m, &old).is_err());
}

#[test]
fn missing_references_declaration_limits_and_changed_bindings_do_not_grant_a_pass() {
    let m = Molecule::dna(
        "unresolved declared requirement",
        "ATGAAA",
        Topology::Linear,
    )
    .unwrap();
    let mut o = Options::default();
    o.required_parts.push(requirement("marker", "missing"));
    let r = assess(&m, &o).unwrap();
    assert_eq!(r.required_parts[0].verdict, Verdict::Indeterminate);
    assert_eq!(r.required_part_searches[0].reference_revision, None);
    let mut changed = o.clone();
    changed.required_parts[0].minimum_copies = 0;
    assert!(
        r.validate_binding(&revision(&m, &changed).unwrap())
            .is_err()
    );
    reference(&mut o, "missing", "ATN");
    assert_eq!(
        assess(&m, &o).unwrap().required_parts[0].verdict,
        Verdict::Indeterminate
    );
    o.required_parts = vec![o.required_parts[0].clone(); 65];
    assert!(matches!(revision(&m, &o), Err(Error::TooLarge)));
    let (m, mut o) = chain_fixture(false);
    o.element_chains = vec![o.element_chains[0].clone(); 65];
    assert!(matches!(revision(&m, &o), Err(Error::TooLarge)));
    let (m, o) = chain_fixture(false);
    let mut r = assess(&m, &o).unwrap();
    r.element_chains[0].upstream.distance_bases = Some(1);
    assert!(r.validate().is_err());
    r.report_revision = r.content_revision().unwrap();
    r.validate().unwrap();
    assert!(r.validate_for(&m, &o).is_err());
}

#[test]
fn repeated_chain_identity_extraction_is_admitted_as_one_aggregate_budget() {
    let mut m = Molecule::dna(
        "bounded repeated claims",
        &"ACGT".repeat(5000),
        Topology::Linear,
    )
    .unwrap();
    m.features = vec![
        feature("CDS", 0, 9, false),
        feature("misc_feature", 0, 20_000, false),
    ];
    let mut o = Options::default();
    o.element_chains = (0..64)
        .map(|i| ElementChain {
            id: format!("chain-{i}"),
            cds_feature: 0,
            upstream: Some(ChainElement {
                feature_index: 1,
                reference_id: Some("missing".into()),
                role: ElementRole::Promoter,
            }),
            downstream: None,
            max_upstream_distance: None,
            max_downstream_distance: None,
        })
        .collect();
    assert!(matches!(revision(&m, &o), Err(Error::TooLarge)));
}
