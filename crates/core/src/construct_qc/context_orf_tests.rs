use super::*;
use initiation_context::{ContextClaim, kozak_1986_pattern};
use unannotated_orfs::{AnnotationRelation, ScanOptions};

fn source(dna: &str, topology: Topology, location: Option<Location>) -> Molecule {
    let mut molecule = Molecule::dna("literal QC context fixture", dna, topology).unwrap();
    if let Some(location) = location {
        molecule.features.push(Feature {
            kind: "CDS".into(),
            location,
            qualifiers: vec![],
            needs_review: false,
        });
    }
    molecule
}
fn options(molecule: &Molecule) -> Options {
    let mut options = Options::default();
    for (index, feature) in molecule.features.iter().enumerate() {
        options.annotations.push(AnnotationAssumption {
            feature_index: index,
            source_revision: molecule.revision(),
            annotation_revision: annotation_revision(feature).unwrap(),
            rationale:
                "Declared fixture CDS role and coordinates, not verified biological annotation"
                    .into(),
        });
    }
    options
}
fn context_options(molecule: &Molecule) -> Options {
    let mut options = options(molecule);
    options.initiation_contexts.push(ContextClaim {
        cds_feature: 0,
        profile: kozak_1986_pattern().unwrap(),
        rationale: "Declare mammalian AUG context and current sequence as the contiguous transcript-context model for this literal comparison".into(),
    });
    options
}
fn scan_options(molecule: &Molecule) -> Options {
    let mut options = options(molecule);
    options.orf_scan = Some(ScanOptions {
        start_policy: translation::StartPolicy::AtgOnly,
        incomplete: translation::IncompletePolicy::CompleteOnly,
        minimum_amino_acids: 2,
        max_candidates: 128,
        both_strands: false,
    });
    options
}

#[test]
fn published_context_literal_conformity_is_pinned_and_not_an_expression_score() {
    let molecule = source(
        "ACCATGGAATAA",
        Topology::Linear,
        Some(Location::range(3, 12)),
    );
    let options = context_options(&molecule);
    let report = assess(&molecule, &options).unwrap();
    let evidence = &report.initiation_contexts[0];
    assert_eq!(evidence.observed.as_deref(), Some("ACCATGG"));
    assert_eq!(evidence.coordinates, [Location::range(0, 7)]);
    assert_eq!(evidence.verdict, Verdict::Pass);
    assert_eq!(
        evidence.profile_revision,
        options.initiation_contexts[0].profile.revision
    );
    assert_eq!(
        evidence.provenance.source,
        "https://doi.org/10.1016/0092-8674(86)90762-2"
    );
    assert!(evidence.reason.contains("not an initiation-strength"));
    assert!(
        report
            .findings
            .iter()
            .find(|finding| finding.rule == Rule::InitiationContext)
            .unwrap()
            .citations
            .iter()
            .any(|citation| citation == &evidence.provenance)
    );
    report.validate_for(&molecule, &options).unwrap();
}

#[test]
fn reverse_and_origin_flanks_use_biological_orientation() {
    let reverse = source(
        "TTATTCCATGGT",
        Topology::Linear,
        Some(Location::Complement(Box::new(Location::range(0, 9)))),
    );
    let report = assess(&reverse, &context_options(&reverse)).unwrap();
    assert_eq!(
        report.initiation_contexts[0].observed.as_deref(),
        Some("ACCATGG")
    );
    assert_eq!(report.initiation_contexts[0].verdict, Verdict::Pass);
    let circular = source(
        "ATGGAATAAACC",
        Topology::Circular,
        Some(Location::range(0, 9)),
    );
    let report = assess(&circular, &context_options(&circular)).unwrap();
    assert_eq!(
        report.initiation_contexts[0].observed.as_deref(),
        Some("ACCATGG")
    );
    assert_eq!(
        report.initiation_contexts[0].coordinates,
        [Location::Join(vec![
            Location::range(9, 12),
            Location::range(0, 4)
        ])]
    );
    assert_eq!(report.initiation_contexts[0].verdict, Verdict::Pass);
}

#[test]
fn iupac_context_certainty_is_subset_or_disjoint_not_arbitrary_resolution() {
    for (dna, expected) in [
        ("NCCATGGAATAA", Verdict::Indeterminate),
        ("GCCATGGAATAA", Verdict::Fail),
    ] {
        let molecule = source(dna, Topology::Linear, Some(Location::range(3, 12)));
        let report = assess(&molecule, &context_options(&molecule)).unwrap();
        assert_eq!(report.initiation_contexts[0].verdict, expected);
    }
    let molecule = source(
        "GCCATGGAATAA",
        Topology::Linear,
        Some(Location::range(3, 12)),
    );
    let mut options = context_options(&molecule);
    let profile = &mut options.initiation_contexts[0].profile;
    profile.id = "declared-purine-fixture".into();
    profile.pattern = Sequence::dna("RCCATGG").unwrap();
    profile.revision = profile.compute_revision().unwrap();
    assert_eq!(
        assess(&molecule, &options).unwrap().initiation_contexts[0].verdict,
        Verdict::Pass
    );
}

#[test]
fn rna_pattern_comparison_preserves_observed_rna() {
    let mut molecule = source(
        "ACCATGGAATAA",
        Topology::Linear,
        Some(Location::range(3, 12)),
    );
    molecule.sequence = Sequence::new(Alphabet::Rna, "ACCAUGGAAUAA").unwrap();
    let report = assess(&molecule, &context_options(&molecule)).unwrap();
    assert_eq!(
        report.initiation_contexts[0].observed.as_deref(),
        Some("ACCAUGG")
    );
    assert_eq!(report.initiation_contexts[0].verdict, Verdict::Pass);
}

#[test]
fn missing_flanks_partial_paths_unreviewed_annotations_and_alternative_starts_are_unresolved() {
    let molecule = source("ATGGAATAA", Topology::Linear, Some(Location::range(0, 9)));
    let evidence = assess(&molecule, &context_options(&molecule))
        .unwrap()
        .initiation_contexts
        .remove(0);
    assert_eq!(evidence.verdict, Verdict::Indeterminate);
    assert!(evidence.observed.is_none());
    let mut alternative = source(
        "ACCTTGGAATAA",
        Topology::Linear,
        Some(Location::range(3, 12)),
    );
    alternative.features[0]
        .qualifiers
        .push(("transl_table".into(), Some("11".into())));
    let mut options = context_options(&alternative);
    options.initiation_contexts[0].profile.genetic_code = 11;
    options.initiation_contexts[0].profile.revision = options.initiation_contexts[0]
        .profile
        .compute_revision()
        .unwrap();
    let evidence = assess(&alternative, &options)
        .unwrap()
        .initiation_contexts
        .remove(0);
    assert_eq!(evidence.verdict, Verdict::Indeterminate);
    assert_eq!(evidence.observed.as_deref(), Some("ACCTTGG"));
    let mut partial = source(
        "ACCATGGAATAA",
        Topology::Linear,
        Some(Location::range(3, 12)),
    );
    partial.features[0]
        .qualifiers
        .push(("partial".into(), None));
    assert_eq!(
        assess(&partial, &context_options(&partial))
            .unwrap()
            .initiation_contexts[0]
            .verdict,
        Verdict::Indeterminate
    );
    partial.features[0].needs_review = true;
    assert_eq!(
        assess(&partial, &context_options(&partial))
            .unwrap()
            .initiation_contexts[0]
            .verdict,
        Verdict::Indeterminate
    );
}

#[test]
fn context_profile_pins_and_fresh_source_replay_reject_tampering() {
    let molecule = source(
        "ACCATGGAATAA",
        Topology::Linear,
        Some(Location::range(3, 12)),
    );
    let mut options = context_options(&molecule);
    let report = assess(&molecule, &options).unwrap();
    options.initiation_contexts[0]
        .profile
        .host_context
        .push_str(" changed");
    assert!(assess(&molecule, &options).is_err());
    let profile = &mut options.initiation_contexts[0].profile;
    profile.revision = profile.compute_revision().unwrap();
    assert!(report.validate_for(&molecule, &options).is_err());
    let mut changed = molecule.clone();
    changed.sequence = Sequence::dna("GCCATGGAATAA").unwrap();
    assert!(
        report
            .validate_for(&changed, &context_options(&changed))
            .is_err()
    );
    let mut tampered = report;
    tampered.initiation_contexts[0].observed = Some("GCCATGG".into());
    assert!(tampered.validate().is_err());
}

#[test]
fn unannotated_literal_orf_is_informational_and_retains_source_coordinates() {
    let molecule = source("ATGAAATAA", Topology::Linear, None);
    let options = scan_options(&molecule);
    let report = assess(&molecule, &options).unwrap();
    let scan = report.orf_scan.as_ref().unwrap();
    assert!(scan.enumeration_complete);
    assert_eq!(scan.candidates.len(), 1);
    let candidate = &scan.candidates[0];
    assert_eq!(candidate.location, Location::range(0, 9));
    assert_eq!(candidate.amino_acids, 2);
    assert_eq!(
        candidate.annotation_relation,
        AnnotationRelation::NoDeclaredCdsCoverage
    );
    assert!(
        report
            .findings
            .iter()
            .filter(|finding| finding.rule == Rule::UnannotatedOrf)
            .all(|finding| finding.verdict == Verdict::Informational)
    );
    report.validate_for(&molecule, &options).unwrap();
}

#[test]
fn nested_starts_are_retained_and_same_register_coverage_is_not_simple_overlap() {
    let molecule = source(
        "ATGATGAAATAA",
        Topology::Linear,
        Some(Location::range(0, 12)),
    );
    let report = assess(&molecule, &scan_options(&molecule)).unwrap();
    let candidates = &report.orf_scan.as_ref().unwrap().candidates;
    assert_eq!(candidates.len(), 2);
    assert_eq!(candidates[0].location, Location::range(0, 12));
    assert_eq!(candidates[1].location, Location::range(3, 12));
    assert!(
        candidates
            .iter()
            .all(|candidate| candidate.annotation_relation
                == AnnotationRelation::CoveredByDeclaredCds)
    );
    // The declared CDS starts one base earlier, so mere interval containment
    // cannot cover the definite ATG candidate's coding register.
    let offset = source(
        "CATGAAATAACC",
        Topology::Linear,
        Some(Location::range(0, 12)),
    );
    let report = assess(&offset, &scan_options(&offset)).unwrap();
    assert_eq!(
        report.orf_scan.as_ref().unwrap().candidates[0].annotation_relation,
        AnnotationRelation::NoDeclaredCdsCoverage
    );
}

#[test]
fn unreviewed_cds_coverage_is_unresolved_instead_of_silently_unannotated() {
    let mut molecule = source("ATGAAATAA", Topology::Linear, Some(Location::range(0, 9)));
    molecule.features[0].needs_review = true;
    let report = assess(&molecule, &scan_options(&molecule)).unwrap();
    assert_eq!(
        report.orf_scan.as_ref().unwrap().candidates[0].annotation_relation,
        AnnotationRelation::Unresolved
    );
}

#[test]
fn reverse_and_circular_orf_candidates_keep_canonical_locations() {
    let reverse = source("TTATTTCAT", Topology::Linear, None);
    let mut options = scan_options(&reverse);
    options.orf_scan.as_mut().unwrap().both_strands = true;
    let report = assess(&reverse, &options).unwrap();
    let candidate = report
        .orf_scan
        .as_ref()
        .unwrap()
        .candidates
        .iter()
        .find(|candidate| candidate.reverse)
        .unwrap();
    assert_eq!(
        candidate.location,
        Location::Complement(Box::new(Location::range(0, 9)))
    );
    assert_eq!(candidate.amino_acids, 2);
    let circular = source(
        "AAATAACCCATG",
        Topology::Circular,
        Some(Location::Join(vec![
            Location::range(9, 12),
            Location::range(0, 6),
        ])),
    );
    let report = assess(&circular, &scan_options(&circular)).unwrap();
    let candidate = &report.orf_scan.as_ref().unwrap().candidates[0];
    assert_eq!(candidate.location, circular.features[0].location);
    assert!(candidate.crosses_origin);
    assert_eq!(
        candidate.annotation_relation,
        AnnotationRelation::CoveredByDeclaredCds
    );
}

#[test]
fn capped_or_ambiguous_scans_never_establish_absence_and_keep_available_candidates() {
    let molecule = source("ATGATGAAATAA", Topology::Linear, None);
    let mut options = scan_options(&molecule);
    options.orf_scan.as_mut().unwrap().max_candidates = 1;
    let report = assess(&molecule, &options).unwrap();
    assert!(!report.orf_scan.as_ref().unwrap().enumeration_complete);
    assert_eq!(report.orf_scan.as_ref().unwrap().candidates.len(), 1);
    assert_eq!(
        report
            .findings
            .iter()
            .find(|finding| finding.rule == Rule::UnannotatedOrf)
            .unwrap()
            .verdict,
        Verdict::Indeterminate
    );
    let ambiguous = source("NTGAAATAN", Topology::Linear, None);
    let report = assess(&ambiguous, &scan_options(&ambiguous)).unwrap();
    let scan = report.orf_scan.as_ref().unwrap();
    assert!(scan.ambiguity_present);
    assert!(scan.candidates.is_empty());
    assert_eq!(
        report
            .findings
            .iter()
            .find(|finding| finding.rule == Rule::UnannotatedOrf)
            .unwrap()
            .verdict,
        Verdict::Indeterminate
    );
}

#[test]
fn new_options_are_required_and_bounded_and_previous_schema_is_refused() {
    let molecule = source("ATGAAATAA", Topology::Linear, None);
    let options = scan_options(&molecule);
    let mut json = serde_json::to_value(&options).unwrap();
    json.as_object_mut().unwrap().remove("orf_scan");
    assert!(serde_json::from_value::<Options>(json).is_err());
    let mut old = options.clone();
    old.schema_version = 2;
    assert!(assess(&molecule, &old).is_err());
    let mut too_many = options.clone();
    too_many.orf_scan.as_mut().unwrap().max_candidates = 129;
    assert!(assess(&molecule, &too_many).is_err());
    let mut no_threshold = options;
    no_threshold.orf_scan.as_mut().unwrap().minimum_amino_acids = 0;
    assert!(assess(&molecule, &no_threshold).is_err());
}

#[test]
fn uncertain_stop_and_explicit_partial_policies_remain_visible() {
    let ambiguous = source("ATGAAATAN", Topology::Linear, None);
    let mut options = scan_options(&ambiguous);
    options.orf_scan.as_mut().unwrap().incomplete =
        translation::IncompletePolicy::IncludeAllPartials;
    let report = assess(&ambiguous, &options).unwrap();
    let scan = report.orf_scan.as_ref().unwrap();
    assert!(scan.ambiguity_present);
    assert_eq!(
        scan.candidates[0].boundary,
        unannotated_orfs::Boundary::AmbiguousStop
    );
    assert!(!scan.candidates[0].has_stop_codon);
    assert!(scan.candidates[0].partial_3prime);
    let open = source("ATGAAA", Topology::Circular, None);
    let mut options = scan_options(&open);
    options.orf_scan.as_mut().unwrap().incomplete =
        translation::IncompletePolicy::IncludeAllPartials;
    let report = assess(&open, &options).unwrap();
    assert_eq!(
        report.orf_scan.as_ref().unwrap().candidates[0].boundary,
        unannotated_orfs::Boundary::OneTurnLimit
    );
}

#[test]
fn context_mismatched_code_and_compound_start_paths_do_not_claim_conformity() {
    let mut molecule = source(
        "ACCATGGAATAA",
        Topology::Linear,
        Some(Location::range(3, 12)),
    );
    molecule.features[0]
        .qualifiers
        .push(("transl_table".into(), Some("11".into())));
    let report = assess(&molecule, &context_options(&molecule)).unwrap();
    assert_eq!(
        report.initiation_contexts[0].verdict,
        Verdict::Indeterminate
    );
    let molecule = source(
        "GAATAAATGACC",
        Topology::Circular,
        Some(Location::Join(vec![
            Location::range(6, 9),
            Location::range(0, 6),
        ])),
    );
    let report = assess(&molecule, &context_options(&molecule)).unwrap();
    assert_eq!(
        report.initiation_contexts[0].verdict,
        Verdict::Indeterminate
    );
}

#[test]
fn orf_reports_bind_the_declared_threshold_code_and_current_annotations() {
    let molecule = source("ATGAAATAA", Topology::Linear, None);
    let mut options = scan_options(&molecule);
    let report = assess(&molecule, &options).unwrap();
    options.orf_scan.as_mut().unwrap().minimum_amino_acids = 3;
    assert!(report.validate_for(&molecule, &options).is_err());
    options.orf_scan.as_mut().unwrap().minimum_amino_acids = 2;
    options.genetic_code = 11;
    assert!(report.validate_for(&molecule, &options).is_err());
    let mut changed = molecule.clone();
    changed.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 9),
        qualifiers: vec![],
        needs_review: false,
    });
    assert!(
        report
            .validate_for(&changed, &scan_options(&changed))
            .is_err()
    );
    let mut tampered = report;
    tampered.orf_scan.as_mut().unwrap().candidates[0].amino_acids = 3;
    assert!(tampered.validate().is_err());
}
