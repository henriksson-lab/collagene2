use super::*;
use collagene_model::{Feature, Location, Sequence};

const GUIDE: &str = "AAACCCCAAAAAAAAAAAAA";
fn molecule() -> Molecule {
    Molecule::dna(
        "literal editor fixture",
        &format!("{GUIDE}AGG"),
        Topology::Linear,
    )
    .unwrap()
}
fn make_request(molecule: &Molecule, reverse: bool, guide: &str) -> Request {
    let guide_options = GuideOptions::default();
    let guide = super::super::design_guides(molecule, &guide_options)
        .unwrap()
        .candidates
        .into_iter()
        .find(|candidate| candidate.reverse == reverse && candidate.guide_sequence == guide)
        .unwrap();
    Request {
        guide_options,
        guide,
        profile: be1_2016().unwrap(),
        intended_positions: vec![4],
        max_outcomes: 4096,
        cds: None,
    }
}
fn cds(molecule: &mut Molecule, location: Location) -> CdsRequest {
    molecule.features.push(Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![],
        needs_review: false,
    });
    CdsRequest {annotation: crate::construct_qc::AnnotationAssumption {
        feature_index: 0, source_revision: molecule.revision(),
        annotation_revision: crate::construct_qc::annotation_revision(&molecule.features[0]).unwrap(),
        rationale: "Explicit literal fixture CDS role/coordinate assumption; no annotation truth is established".into(),
    }, genetic_code: 1}
}
#[test]
fn literal_window_subsets_retain_intended_and_bystander_changes_without_probabilities() {
    let molecule = molecule();
    let request = make_request(&molecule, false, GUIDE);
    let report = assess(&molecule, &request).unwrap();
    assert_eq!(
        report
            .sites
            .iter()
            .map(|site| site.guide_position)
            .collect::<Vec<_>>(),
        [4, 5, 6, 7]
    );
    assert_eq!(
        report
            .sites
            .iter()
            .map(|site| site.reference_position)
            .collect::<Vec<_>>(),
        [3, 4, 5, 6]
    );
    assert_eq!(report.finite_state_count, 16);
    assert_eq!(report.states.len(), 16);
    assert!(report.finite_enumeration_complete);
    assert!(!report.biological_outcomes_exhaustive);
    assert!(report.states[0].edited_positions.is_empty());
    assert_eq!(report.states[0].reference_guide, GUIDE);
    assert_eq!(report.states[0].satisfies_intent, Some(false));
    assert_eq!(report.states[1].edited_positions, [4]);
    assert_eq!(report.states[1].satisfies_intent, Some(true));
    assert!(report.states[1].bystander_positions.is_empty());
    assert_eq!(report.states[2].edited_positions, [5]);
    assert_eq!(report.states[2].satisfies_intent, Some(false));
    assert_eq!(report.states[15].bystander_positions, [5, 6, 7]);
    assert_eq!(report.states[15].reference_guide, "AAATTTTAAAAAAAAAAAAA");
    report.validate_for(&molecule, &request).unwrap();
}
#[test]
fn reverse_cytosines_map_to_source_forward_guanine_to_adenine() {
    let molecule = Molecule::dna(
        "reverse fixture",
        &format!("CCT{}GGGGTTT", "T".repeat(13)),
        Topology::Linear,
    )
    .unwrap();
    let request = make_request(&molecule, true, GUIDE);
    let report = assess(&molecule, &request).unwrap();
    assert_eq!(
        report
            .sites
            .iter()
            .map(|site| site.reference_position)
            .collect::<Vec<_>>(),
        [19, 18, 17, 16]
    );
    assert!(
        report
            .sites
            .iter()
            .all(|site| site.reference_from == 'G' && site.reference_to == 'A')
    );
    assert_eq!(report.states[1].reference_guide, "AAATCCCAAAAAAAAAAAAA");
    assert!(
        report
            .sites
            .iter()
            .all(|site| molecule.sequence.as_str().as_bytes()[site.reference_position] == b'G')
    );
}
#[test]
fn circular_origin_positions_come_from_the_canonical_guide_owner() {
    let source = molecule();
    let bases = source.sequence.as_str();
    let circular = Molecule::dna(
        "rotated literal",
        &format!("{}{}", &bases[7..], &bases[..7]),
        Topology::Circular,
    )
    .unwrap();
    let request = make_request(&circular, false, GUIDE);
    assert_eq!(request.guide.guide_positions[0], 16);
    let report = assess(&circular, &request).unwrap();
    assert_eq!(
        report
            .sites
            .iter()
            .map(|site| site.reference_position)
            .collect::<Vec<_>>(),
        [19, 20, 21, 22]
    );
    assert_eq!(report.states[15].reference_guide, "AAATTTTAAAAAAAAAAAAA");
}
#[test]
fn capped_enumeration_is_not_a_complete_finite_or_biological_outcome_set() {
    let molecule = molecule();
    let mut request = make_request(&molecule, false, GUIDE);
    request.max_outcomes = 1;
    let report = assess(&molecule, &request).unwrap();
    assert_eq!(report.finite_state_count, 16);
    assert_eq!(report.states.len(), 1);
    assert!(!report.finite_enumeration_complete);
    assert!(!report.biological_outcomes_exhaustive);
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.rule == Rule::FiniteEnumeration
                && finding.verdict == Verdict::Indeterminate)
    );
}
#[test]
fn an_unreachable_intended_edit_fails_only_the_declared_finite_model() {
    let molecule = molecule();
    let mut request = make_request(&molecule, false, GUIDE);
    request.intended_positions = vec![3]; // A, outside the C4–8 window.
    let report = assess(&molecule, &request).unwrap();
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.rule == Rule::IntendedEdits && finding.verdict == Verdict::Fail)
    );
    assert!(
        report
            .states
            .iter()
            .all(|state| state.satisfies_intent == Some(false))
    );
    request.intended_positions.clear();
    let report = assess(&molecule, &request).unwrap();
    assert!(
        report
            .states
            .iter()
            .all(|state| state.satisfies_intent.is_none() && state.bystander_positions.is_empty())
    );
}
#[test]
fn source_guide_and_profile_tampering_cannot_reuse_evidence() {
    let molecule = molecule();
    let mut request = make_request(&molecule, false, GUIDE);
    let report = assess(&molecule, &request).unwrap();
    request.guide.guide_positions[0] = 1;
    assert!(assess(&molecule, &request).is_err());
    let mut changed = molecule.clone();
    changed.sequence = Sequence::dna(&format!("AAATCCCAAAAAAAAAAAAAAGG")).unwrap();
    assert!(report.validate_for(&changed, &report.request).is_err());
    let mut profile = request.profile.clone();
    profile.window_first = 5;
    assert!(profile.validate().is_err());
    profile.revision = profile.compute_revision().unwrap();
    assert!(
        profile.validate().is_err(),
        "altered builtin is not admitted as a current profile"
    );
    let mut tampered = report;
    tampered.biological_outcomes_exhaustive = true;
    assert!(
        tampered
            .validate_binding(&tampered.context_revision)
            .is_err()
    );
}
#[test]
fn assumed_cds_ordinary_codon_change_uses_the_canonical_genetic_code() {
    let guide = format!("ATGCAA{}", "A".repeat(14));
    let mut molecule = Molecule::dna(
        "stop-gain literal",
        &format!("{guide}AGG"),
        Topology::Linear,
    )
    .unwrap();
    let selected_cds = cds(&mut molecule, Location::range(0, 21));
    let mut request = make_request(&molecule, false, &guide);
    request.cds = Some(selected_cds);
    let report = assess(&molecule, &request).unwrap();
    assert!(report.cds.as_ref().unwrap().available);
    assert!(report.states[0].codon_consequences.is_empty());
    let consequence = &report.states[1].codon_consequences[0];
    assert_eq!(consequence.codon_number, 2);
    assert_eq!(consequence.reference_positions, [3, 4, 5]);
    assert_eq!(consequence.before, "CAA");
    assert_eq!(consequence.after, "TAA");
    assert_eq!(consequence.before_amino_acid, 'Q');
    assert_eq!(consequence.after_amino_acid, '*');
    assert_eq!(consequence.kind, ConsequenceKind::StopGained);
}
#[test]
fn reverse_cds_and_reverse_guide_preserve_the_actual_codon_coordinate_order() {
    let guide = format!("ATGCAA{}", "A".repeat(14));
    // Reverse complement of the independently specified 23-base source above.
    let mut molecule = Molecule::dna(
        "reverse coding literal",
        &format!("CCT{}TTGCAT", "T".repeat(14)),
        Topology::Linear,
    )
    .unwrap();
    let selected_cds = cds(
        &mut molecule,
        Location::Complement(Box::new(Location::range(2, 23))),
    );
    let mut request = make_request(&molecule, true, &guide);
    request.cds = Some(selected_cds);
    let report = assess(&molecule, &request).unwrap();
    let consequence = &report.states[1].codon_consequences[0];
    assert_eq!(consequence.reference_positions, [19, 18, 17]);
    assert_eq!(consequence.before, "CAA");
    assert_eq!(consequence.after, "TAA");
    assert_eq!(consequence.kind, ConsequenceKind::StopGained);
}
#[test]
fn stale_review_required_joined_or_wrong_code_cds_never_claims_consequences() {
    let mut molecule = molecule();
    let mut selected = cds(&mut molecule, Location::range(0, 21));
    let mut request = make_request(&molecule, false, GUIDE);
    selected.annotation.source_revision = "0".repeat(64);
    request.cds = Some(selected);
    let report = assess(&molecule, &request).unwrap();
    assert!(!report.cds.as_ref().unwrap().available);
    assert!(
        report
            .states
            .iter()
            .all(|state| state.codon_consequences.is_empty())
    );
    let mut selected = request.cds.as_ref().unwrap().clone();
    selected.annotation.source_revision = molecule.revision();
    selected.annotation.annotation_revision =
        crate::construct_qc::annotation_revision(&molecule.features[0]).unwrap();
    selected.genetic_code = 11;
    request.cds = Some(selected);
    assert!(
        !assess(&molecule, &request)
            .unwrap()
            .cds
            .as_ref()
            .unwrap()
            .available
    );
    molecule.features[0].needs_review = true;
    let mut request = make_request(&molecule, false, GUIDE);
    request.cds = Some(CdsRequest {
        annotation: crate::construct_qc::AnnotationAssumption {
            feature_index: 0,
            source_revision: molecule.revision(),
            annotation_revision: crate::construct_qc::annotation_revision(&molecule.features[0])
                .unwrap(),
            rationale: "This declaration does not override the feature's review-required state"
                .into(),
        },
        genetic_code: 1,
    });
    assert!(
        !assess(&molecule, &request)
            .unwrap()
            .cds
            .as_ref()
            .unwrap()
            .available
    );
    molecule.features[0].needs_review = false;
    molecule.features[0].location =
        Location::Join(vec![Location::range(0, 3), Location::range(3, 21)]);
    let mut request = make_request(&molecule, false, GUIDE);
    request.cds = Some(CdsRequest {annotation: crate::construct_qc::AnnotationAssumption {
        feature_index: 0, source_revision: molecule.revision(),
        annotation_revision: crate::construct_qc::annotation_revision(&molecule.features[0]).unwrap(),
        rationale: "The joined path is explicit but not supported by this consequence-coordinate profile".into(),
    }, genetic_code: 1});
    assert!(
        !assess(&molecule, &request)
            .unwrap()
            .cds
            .as_ref()
            .unwrap()
            .available
    );
}
#[test]
fn current_request_requires_the_nullable_cds_field_and_explicit_limits() {
    let molecule = molecule();
    let mut request = make_request(&molecule, false, GUIDE);
    let mut encoded = serde_json::to_value(&request).unwrap();
    encoded.as_object_mut().unwrap().remove("cds");
    assert!(serde_json::from_value::<Request>(encoded).is_err());
    request.max_outcomes = 0;
    assert!(assess(&molecule, &request).is_err());
    request.max_outcomes = 4097;
    assert!(assess(&molecule, &request).is_err());
}

#[test]
fn explicitly_declared_adenine_transition_uses_the_same_strand_and_subset_owner() {
    let molecule = molecule();
    let mut request = make_request(&molecule, false, GUIDE);
    request.profile.id = "declared-adenine-fixture".into();
    request.profile.custom = true;
    request.profile.conversion = Conversion::AdenineTransition;
    request.profile.window_first = 8;
    request.profile.window_last = 9;
    request.profile.rationale = "Explicit synthetic finite A→G window fixture; the attached BE1 citation does not establish this custom editor's biological activity or geometry.".into();
    request.profile.revision = request.profile.compute_revision().unwrap();
    request.intended_positions = vec![8];
    let report = assess(&molecule, &request).unwrap();
    assert_eq!(report.finite_state_count, 4);
    assert_eq!(
        report
            .sites
            .iter()
            .map(|site| site.guide_position)
            .collect::<Vec<_>>(),
        [8, 9]
    );
    assert!(
        report
            .sites
            .iter()
            .all(|site| site.reference_from == 'A' && site.reference_to == 'G')
    );
    assert_eq!(report.states[3].bystander_positions, [9]);
}

#[test]
fn report_bytes_bound_stops_large_finite_sets_without_claiming_completion() {
    let guide = format!("AAA{}AAAAA", "C".repeat(12));
    let mut molecule = Molecule::dna(
        "bounded consequence set",
        &format!("{guide}AGG"),
        Topology::Linear,
    )
    .unwrap();
    let selected_cds = cds(&mut molecule, Location::range(0, 21));
    let mut request = make_request(&molecule, false, &guide);
    request.cds = Some(selected_cds);
    request.profile.id = "declared-twelve-position-fixture".into();
    request.profile.custom = true;
    request.profile.window_last = 15;
    request.profile.rationale = "Explicit synthetic combinatorial bound fixture; not a published 12-base editor activity claim.".into();
    request.profile.revision = request.profile.compute_revision().unwrap();
    let report = assess(&molecule, &request).unwrap();
    assert_eq!(report.finite_state_count, 4096);
    assert!(!report.finite_enumeration_complete);
    assert!(!report.states.is_empty());
    assert!(report.states.len() < 4096);
    assert!(serde_json::to_vec(&report).unwrap().len() <= MAX_REPORT_BYTES);
}

#[test]
fn ambiguous_guides_are_refused_by_the_existing_selection_owner() {
    let mut molecule = molecule();
    molecule.sequence = Sequence::dna("AAANCCCAAAAAAAAAAAAAAGG").unwrap();
    let options = GuideOptions::default();
    let guide = super::super::design_guides(&molecule, &options)
        .unwrap()
        .candidates
        .into_iter()
        .find(|candidate| !candidate.reverse && candidate.guide_positions[0] == 0)
        .unwrap();
    assert!(!guide.concrete);
    let request = Request {
        guide_options: options,
        guide,
        profile: be1_2016().unwrap(),
        intended_positions: vec![4],
        max_outcomes: 32,
        cds: None,
    };
    assert!(assess(&molecule, &request).is_err());
}
