//! Catalog-pinned gel lanes replay the complete recipe; saved reports are not bands.
use collagene_core::{
    analysis,
    catalog_digest::{self, InspectionOptions},
    digest_setups::{CleavageMode, Recipe},
    gel::{self, Calibration, CalibrationPoint, LaneInput, MobilityState, Preparation},
    run_lineage::EventBinding,
    workspace::Workspace,
};
use collagene_model::{Molecule, Topology};

fn wide() -> Calibration {
    Calibration {
        label: "Independent illustrative span".into(),
        points: vec![
            CalibrationPoint {
                size_bp: 1.,
                migration: 0.9,
            },
            CalibrationPoint {
                size_bp: 10000.,
                migration: 0.1,
            },
        ],
        resolution: 0.01,
    }
}
fn recipe(molecule: &Molecule, ids: &[&str], options: InspectionOptions) -> Recipe {
    let frame = Workspace::from_documents("gel source", vec![molecule.clone()], 0).unwrap();
    let doc = &frame.documents[0];
    Recipe {
        source: EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        },
        selection: catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &ids.iter().map(|id| id.to_string()).collect::<Vec<_>>(),
        )
        .unwrap(),
        options,
        cleavage_mode: CleavageMode::AssumedCompleteMethylationUnassessed,
    }
}
fn input(molecule: Molecule, recipe: Recipe) -> LaneInput {
    LaneInput {
        label: "Pinned diagnostic digest".into(),
        molecule,
        preparation: Preparation::CatalogDigest { recipe },
        mass_ng: Some(100.),
    }
}

#[test]
fn all_starter_catalog_geometries_preserve_existing_gel_fragment_pipeline() {
    for enzyme in analysis::starter_enzymes() {
        let molecule = Molecule::dna(
            "starter comparison",
            &format!("AAA{}AAAAAAAAAAAAAAAA", enzyme.recognition),
            Topology::Linear,
        )
        .unwrap();
        let recipe = recipe(
            &molecule,
            &[&enzyme.name.to_ascii_lowercase()],
            InspectionOptions::default(),
        );
        let pinned = input(molecule.clone(), recipe);
        let ordinary = LaneInput {
            label: pinned.label.clone(),
            molecule,
            preparation: Preparation::Digest {
                enzymes: vec![enzyme.name.clone()],
            },
            mass_ng: pinned.mass_ng,
        };
        let report = gel::predict_gel(&[ordinary, pinned], &wide()).unwrap();
        let (old, new) = (&report.lanes[0], &report.lanes[1]);
        assert!(old.error.is_none(), "{}: {:?}", enzyme.name, old.error);
        assert!(new.error.is_none(), "{}: {:?}", enzyme.name, new.error);
        assert_eq!(
            serde_json::to_value(&old.fragments).unwrap(),
            serde_json::to_value(&new.fragments).unwrap(),
            "{}",
            enzyme.name
        );
        assert_eq!(
            serde_json::to_value(&old.bands).unwrap(),
            serde_json::to_value(&new.bands).unwrap(),
            "{}",
            enzyme.name
        );
        assert!(new.complete_under_model);
        assert!(
            new.warnings
                .iter()
                .any(|warning| warning.contains("methylation remains unassessed"))
        );
    }
}

#[test]
fn nonstarter_catalog_identity_works_without_name_lookup_fallback() {
    let molecule = Molecule::dna("AluI source", "TTTAGCTAAA", Topology::Linear).unwrap();
    let pinned = input(
        molecule.clone(),
        recipe(&molecule, &["alui"], InspectionOptions::default()),
    );
    let ordinary = LaneInput {
        label: "Names-only legacy path".into(),
        molecule,
        preparation: Preparation::Digest {
            enzymes: vec!["AluI".into()],
        },
        mass_ng: None,
    };
    let report = gel::predict_gel(&[pinned, ordinary], &wide()).unwrap();
    assert!(report.lanes[0].error.is_none());
    assert_eq!(
        report.lanes[0]
            .fragments
            .iter()
            .map(|f| (f.top_nt, f.bottom_nt))
            .collect::<Vec<_>>(),
        vec![(5, 5), (5, 5)]
    );
    assert!(report.lanes[0].preparation.contains("AluI"));
    assert!(report.lanes[1].error.is_some());
    assert!(report.lanes[1].fragments.is_empty());
}

#[test]
fn staggered_lengths_use_both_physical_strands_and_conserve_mass() {
    let molecule = Molecule::dna("HinfI stagger", "AAAGACTCAAA", Topology::Linear).unwrap();
    let lane = input(
        molecule.clone(),
        recipe(&molecule, &["hinfi"], InspectionOptions::default()),
    );
    let report = gel::predict_gel(&[lane], &wide()).unwrap();
    let lane = &report.lanes[0];
    assert!(lane.error.is_none());
    assert_eq!(
        lane.fragments
            .iter()
            .map(|f| (f.top_nt, f.bottom_nt))
            .collect::<Vec<_>>(),
        vec![(4, 7), (7, 4)]
    );
    assert!(
        lane.fragments
            .iter()
            .all(|f| f.size_bp == 5.5 && f.state == MobilityState::LinearDuplex)
    );
    assert_eq!(
        lane.fragments
            .iter()
            .map(|f| f.amount_ng.unwrap())
            .sum::<f64>(),
        100.
    );
    assert!(
        lane.warnings
            .iter()
            .any(|warning| warning.contains("mean strand length"))
    );
    assert_eq!(lane.bands.len(), 1);
}

#[test]
fn uncut_circle_remains_unpositioned_and_one_cut_circle_is_linear() {
    let uncut = Molecule::dna("uncut circle", "AAAACCCC", Topology::Circular).unwrap();
    let cut = Molecule::dna("one cut circle", "TTTGAATTCAAA", Topology::Circular).unwrap();
    let report = gel::predict_gel(
        &[
            input(
                uncut.clone(),
                recipe(&uncut, &["ecori"], InspectionOptions::default()),
            ),
            input(
                cut.clone(),
                recipe(&cut, &["ecori"], InspectionOptions::default()),
            ),
        ],
        &wide(),
    )
    .unwrap();
    let circle = &report.lanes[0];
    assert!(circle.error.is_none());
    assert_eq!(circle.fragments.len(), 1);
    assert_eq!(
        circle.fragments[0].state,
        MobilityState::CircularConformationUnknown
    );
    assert!(circle.fragments[0].migration.is_none());
    assert!(circle.bands.is_empty());
    assert_eq!(circle.unpositioned_fragments, vec![0]);
    assert!(
        circle
            .warnings
            .iter()
            .any(|warning| warning.contains("Uncut circular conformation is unknown"))
    );
    let linear = &report.lanes[1];
    assert!(linear.error.is_none());
    assert_eq!(
        (linear.fragments[0].top_nt, linear.fragments[0].bottom_nt),
        (12, 12)
    );
    assert_eq!(linear.fragments[0].state, MobilityState::LinearDuplex);
    assert!(linear.fragments[0].migration.is_some());
}

#[test]
fn partial_ambiguous_outside_and_unmodeled_physical_products_are_withheld() {
    let exact = Molecule::dna("partial source", "TTTGAATTCAAA", Topology::Linear).unwrap();
    let ambiguous = Molecule::dna("possible EcoRI", "TTTGAANTCAAA", Topology::Linear).unwrap();
    let outside = Molecule::dna("outside Type IIS", "AAAGGTCTC", Topology::Linear).unwrap();
    let mut physical = exact.clone();
    physical
        .metadata
        .insert("duplex_model".into(), "untrusted physical evidence".into());
    let mut comparison_limits = InspectionOptions::default();
    comparison_limits.max_comparisons = 1;
    let mut fragment_limits = InspectionOptions::default();
    fragment_limits.max_fragments = 1;
    let report = gel::predict_gel(
        &[
            input(exact.clone(), recipe(&exact, &["ecori"], comparison_limits)),
            input(exact.clone(), recipe(&exact, &["ecori"], fragment_limits)),
            input(
                ambiguous.clone(),
                recipe(&ambiguous, &["ecori"], InspectionOptions::default()),
            ),
            input(
                outside.clone(),
                recipe(&outside, &["bsai"], InspectionOptions::default()),
            ),
            input(
                physical.clone(),
                recipe(&physical, &["ecori"], InspectionOptions::default()),
            ),
        ],
        &wide(),
    )
    .unwrap();
    for lane in &report.lanes {
        assert!(lane.error.as_ref().unwrap().contains("withheld"));
        assert!(!lane.complete_under_model);
        assert!(lane.fragments.is_empty());
        assert!(lane.bands.is_empty());
    }
}

#[test]
fn current_catalog_and_actual_source_revision_must_match_the_pinned_recipe() {
    let molecule = Molecule::dna("source pins", "TTTGAATTCAAA", Topology::Linear).unwrap();
    let valid = recipe(&molecule, &["ecori"], InspectionOptions::default());
    let mut drift = valid.clone();
    drift.selection.catalog_edition = "unavailable future catalog".into();
    let mut geometry = valid.clone();
    geometry.selection.enzymes[0].enzyme.top_cut = 2;
    let mut source = valid.clone();
    source.source.revision = "0".repeat(64);
    let report = gel::predict_gel(
        &[
            input(molecule.clone(), drift),
            input(molecule.clone(), geometry),
            input(molecule.clone(), source),
        ],
        &wide(),
    )
    .unwrap();
    assert!(report.lanes.iter().all(|lane| lane.error.is_some()
        && lane.fragments.is_empty()
        && !lane.complete_under_model));
    let mut event = valid.clone();
    event.source.event_id = "not-an-exact-event".into();
    assert!(gel::predict_gel(&[input(molecule.clone(), event)], &wide()).is_err());
    let mut dimension = valid.clone();
    dimension.options.max_sites = 0;
    assert!(gel::predict_gel(&[input(molecule, dimension)], &wide()).is_err());
    // LaneInput has no workspace ancestry. The host must resolve the exact event;
    // matching DNA alone cannot authenticate another document/event identity.
}

#[test]
fn gel_fragment_limit_is_an_explicit_lane_error_without_rewriting_recipe_limits() {
    let molecule = Molecule::dna("many sites", &"GAATTC".repeat(1025), Topology::Circular).unwrap();
    let recipe = recipe(&molecule, &["ecori"], InspectionOptions::default());
    assert!(recipe.options.max_fragments > 1024);
    let lane = input(molecule, recipe);
    let original = serde_json::to_value(&lane).unwrap();
    let report = gel::predict_gel(std::slice::from_ref(&lane), &wide()).unwrap();
    let result = &report.lanes[0];
    assert!(result.error.as_ref().unwrap().contains("1024"));
    assert!(result.error.as_ref().unwrap().contains("not changed"));
    assert!(result.fragments.is_empty() && result.bands.is_empty());
    assert!(!result.complete_under_model);
    assert_eq!(serde_json::to_value(lane).unwrap(), original);
}

#[test]
fn oversized_borrowed_recipe_or_source_rejects_before_prediction() {
    let molecule = Molecule::dna("bounded source", "TTTGAATTCAAA", Topology::Linear).unwrap();
    let mut malformed = recipe(&molecule, &["ecori"], InspectionOptions::default());
    malformed.selection.enzymes = vec![malformed.selection.enzymes[0].clone(); 65];
    assert!(gel::predict_gel(&[input(molecule.clone(), malformed)], &wide()).is_err());
    let mut oversized = molecule.clone();
    oversized
        .metadata
        .insert("oversized".into(), "X".repeat(1_048_577));
    let original_recipe = recipe(&molecule, &["ecori"], InspectionOptions::default());
    assert!(gel::predict_gel(&[input(oversized, original_recipe)], &wide()).is_err());
}

#[test]
fn legacy_large_annotation_lane_keeps_native_capacity_alongside_pinned_catalog_lane() {
    let mut legacy = Molecule::dna("legacy annotation capacity", "ACGT", Topology::Linear).unwrap();
    legacy.features = vec![
        collagene_model::Feature {
            kind: "misc_feature".into(),
            location: collagene_model::Location::range(0, 4),
            qualifiers: Vec::new(),
            needs_review: false,
        };
        10_000
    ];
    let pinned = Molecule::dna("pinned AluI source", "TTTAGCTAAA", Topology::Linear).unwrap();
    let recipe = recipe(&pinned, &["alui"], InspectionOptions::default());
    let report = gel::predict_gel(
        &[
            LaneInput {
                label: "Legacy uncut lane".into(),
                molecule: legacy.clone(),
                preparation: Preparation::Uncut,
                mass_ng: None,
            },
            input(pinned, recipe),
        ],
        &wide(),
    )
    .unwrap();
    assert!(report.lanes.iter().all(|lane| lane.error.is_none()));
    assert_eq!(
        report.lanes[0]
            .fragments
            .iter()
            .map(|fragment| (fragment.top_nt, fragment.bottom_nt))
            .collect::<Vec<_>>(),
        vec![(4, 4)]
    );
    assert_eq!(report.lanes[1].fragments.len(), 2);
    assert_eq!(legacy.features.len(), 10_000);
    let diagnostic = gel::diagnostic_digests(
        &[legacy.clone(), legacy],
        &gel::DiagnosticOptions {
            enzymes: vec!["EcoRI".into()],
            include_double_digests: false,
            calibration: wide(),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(diagnostic.evaluated_digests, 1);
    assert!(
        diagnostic.candidates[0]
            .profiles
            .iter()
            .all(|lane| lane.error.is_none())
    );
}
