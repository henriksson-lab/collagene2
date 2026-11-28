use collagene_core::{
    cloning::{
        assembly::recipe::{FragmentSpec, Preparation},
        golden_gate::{GoldenGateEnzyme, GoldenGateOptions, prepared::*, primers},
        ligation::recipe as ligation_recipe,
    },
    planning::plan::{Operation, PlanInput},
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
const A: &str = "ATGACCCCGGGGTTTCCAGA";
const B: &str = "ATGAGTTAGCGACCTTGGAA";
fn input(id: &str, dna: &str) -> PlanInput {
    let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fixture(a: &str, b: &str, right_a: &str, right_b: &str) -> (Vec<PlanInput>, Molecule, Request) {
    let inputs = vec![input("a", a), input("b", b)];
    let target = Molecule::dna(
        "independent intended circle",
        &format!("{a}{b}"),
        Topology::Circular,
    )
    .unwrap();
    let request = primers::Request {
        name: "Two actual BsaI preparations".into(),
        fragments: vec![
            primers::Fragment {
                role: "first".into(),
                input_id: "a".into(),
                start: 0,
                end: a.len(),
                reverse: false,
                forward_length: 8,
                reverse_length: 8,
                right_overhang: right_a.into(),
            },
            primers::Fragment {
                role: "second".into(),
                input_id: "b".into(),
                start: 0,
                end: b.len(),
                reverse: false,
                forward_length: 8,
                reverse_length: 8,
                right_overhang: right_b.into(),
            },
        ],
        enzyme: GoldenGateEnzyme::BsaI,
        clamp: "ACGT".into(),
        spacer: "A".into(),
    };
    let designed = primers::design(&inputs, &target, &request).unwrap();
    let recipe = ligation_recipe::Recipe {
        standard: None,
        name: "Ordered substrate declarations".into(),
        fragments: designed
            .fragments
            .into_iter()
            .map(|source| ligation_recipe::Fragment {
                source,
                treatment: None,
            })
            .collect(),
        circular: true,
        allow_nicked: false,
    };
    (
        inputs,
        target,
        Request {
            recipe,
            options: GoldenGateOptions {
                both_orientations: false,
                max_fragments: 2,
                ..Default::default()
            },
        },
    )
}
fn selection(report: &Report, index: usize) -> Selection {
    let c = &report.pool.as_ref().unwrap().candidates[index];
    Selection {
        context_revision: report.context_revision.clone(),
        candidate_index: index,
        directed_sequence_key: c.directed_sequence_key.clone(),
        physical_sequence_key: c.physical_sequence_key.clone(),
    }
}
#[test]
fn physical_pool_keeps_self_subset_and_full_alternatives_and_selection_preserves_all_preparations()
{
    let (inputs, target, r) = fixture(A, B, "ATGA", "ATGA");
    let before = serde_json::to_value(&inputs).unwrap();
    assert_eq!(
        target.sequence.as_str(),
        "ATGACCCCGGGGTTTCCAGAATGAGTTAGCGACCTTGGAA"
    );
    let report = inspect(&inputs, &target, &r).unwrap();
    report
        .validate_binding(&inspect_revision(&inputs, &target, &r).unwrap())
        .unwrap();
    let pool = report.pool.as_ref().unwrap();
    assert!(pool.enumeration_complete);
    assert!(pool.stable_physical_product_count >= 3);
    assert!(report.assessments.iter().any(|a| a.source_roles == vec![0]));
    assert!(report.assessments.iter().any(|a| a.source_roles == vec![1]));
    let full = report
        .assessments
        .iter()
        .find(|a| a.target_matches && a.source_roles.len() == 2)
        .unwrap()
        .candidate_index;
    let selected = select(&report, &selection(&report, full)).unwrap();
    selected
        .validate_binding(&selection_revision(&report, &selection(&report, full)).unwrap())
        .unwrap();
    let old = &report.preview.compiled.plan;
    let n = old.steps.len() - 1;
    assert_eq!(
        serde_json::to_value(&selected.plan.inputs).unwrap(),
        serde_json::to_value(&old.inputs).unwrap()
    );
    assert_eq!(
        serde_json::to_value(&selected.plan.steps[..n]).unwrap(),
        serde_json::to_value(&old.steps[..n]).unwrap()
    );
    assert!(matches!(
        selected.plan.steps.last(),
        Some(Operation::GoldenGate { .. })
    ));
    assert_eq!(
        verification::engine::verify(&selected.plan, &target).status,
        Status::ConsistentUnderAssumptions
    );
    let self_a = report
        .assessments
        .iter()
        .find(|a| a.source_roles == vec![0])
        .unwrap()
        .candidate_index;
    let selected_self = select(&report, &selection(&report, self_a)).unwrap();
    assert_eq!(
        verification::engine::verify(&selected_self.plan, &target).status,
        Status::Invalid
    );
    assert_eq!(selected_self.plan.inputs.len(), 2);
    assert_eq!(serde_json::to_value(&inputs).unwrap(), before);
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("Discarded digest fragments"))
    );
}
#[test]
fn missing_middle_preparation_and_metadata_only_end_claims_withhold_the_entire_pool() {
    let (inputs, target, mut r) = fixture(A, B, "ATGA", "ATGA");
    if let Preparation::PcrCatalogDigest { forward, .. } =
        &mut r.recipe.fragments[1].source.preparation
    {
        *forward = Primer::new("missing actual forward", "AAAAAAAA", "").unwrap();
    }
    let report = inspect(&inputs, &target, &r).unwrap();
    assert_eq!(report.preview.prepared.len(), 2);
    assert!(report.preview.prepared[0].is_some());
    assert!(report.preview.prepared[1].is_none());
    assert!(!report.preview.preparation_complete);
    assert_eq!(report.missing_physical, vec![1]);
    assert!(report.pool.is_none());
    assert!(report.assessments.is_empty());
    assert_eq!(report.preview.verification.reconstructed_count(), 2);
    let mut claimed = inputs.clone();
    claimed[0].molecule.metadata.insert(
        "duplex_model".into(),
        "Imported end summary cannot create a live duplex".into(),
    );
    claimed[0].revision = claimed[0].molecule.revision();
    let rows = claimed
        .iter()
        .enumerate()
        .map(|(i, s)| ligation_recipe::Fragment {
            source: FragmentSpec {
                generation: None,
                role: format!("raw-{i}"),
                input_id: s.id.clone(),
                preparation: Preparation::WholeLinear,
                reverse: false,
            },
            treatment: None,
        })
        .collect();
    let request = Request {
        recipe: ligation_recipe::Recipe {
            standard: None,
            name: "Only selected raw references".into(),
            fragments: rows,
            circular: true,
            allow_nicked: false,
        },
        options: r.options,
    };
    let report = inspect(&claimed, &target, &request).unwrap();
    assert!(report.preview.preparation_complete);
    assert_eq!(report.preview.prepared.len(), 2);
    assert!(report.missing_physical.contains(&0));
    assert!(report.pool.is_none());
}
#[test]
fn bounded_search_retains_candidates_without_claiming_complete_or_exhaustive_mixture_coverage() {
    let (inputs, target, mut r) = fixture(A, B, "ATGA", "ATGA");
    r.options.max_states = 1;
    let report = inspect(&inputs, &target, &r).unwrap();
    let pool = report.pool.as_ref().unwrap();
    assert!(!pool.enumeration_complete);
    assert!(!pool.exhaustive_mixture_search);
    assert!(!pool.truncation_reasons.is_empty());
    assert_eq!(report.assessments.len(), pool.candidates.len());
    if let Some(a) = report.assessments.first() {
        let selected = select(&report, &selection(&report, a.candidate_index)).unwrap();
        assert!(selected.warnings.iter().any(|w| w.contains("truncated")));
    }
    // An explicitly prepared uncut circle is carryover of that substrate, not an
    // invented newly assembled candidate or omitted original reaction material.
    let mut original = input("uncut", A);
    original.molecule.topology = Topology::Circular;
    original.revision = original.molecule.revision();
    let target = Molecule::dna("independent original circle", A, Topology::Circular).unwrap();
    let r = Request {
        recipe: ligation_recipe::Recipe {
            standard: None,
            name: "Explicit intact prepared circle".into(),
            fragments: vec![ligation_recipe::Fragment {
                source: FragmentSpec {
                    generation: None,
                    role: "carryover".into(),
                    input_id: "uncut".into(),
                    preparation: Preparation::Digest {
                        enzymes: vec!["EcoRV".into()],
                        fragment_index: 0,
                        repair: None,
                    },
                    reverse: false,
                },
                treatment: None,
            }],
            circular: true,
            allow_nicked: false,
        },
        options: GoldenGateOptions {
            both_orientations: false,
            ..Default::default()
        },
    };
    let report = inspect(&[original], &target, &r).unwrap();
    let pool = report.pool.as_ref().unwrap();
    assert_eq!(pool.carryover_sources, vec![0]);
    assert!(pool.candidates.is_empty());
}
#[test]
fn newly_recreated_recognition_site_stays_visible_but_cannot_be_selected_as_site_free() {
    // Last A base G + first five B bases GTCTC form GGTCTC only after
    // joining: the four-base GTCT fusion alone does not contain that full site.
    let (inputs, target, r) = fixture(
        "ATGACCCCGGGGTTTCCAGG",
        "GTCTCGTAGCGACCTTGGAA",
        "GTCT",
        "ATGA",
    );
    let report = inspect(&inputs, &target, &r).unwrap();
    let pool = report.pool.as_ref().unwrap();
    let recut = pool
        .candidates
        .iter()
        .position(|c| !c.restriction_site_free)
        .expect("literal boundary creates BsaI recognition");
    assert!(!pool.candidates[recut].recut_sites.is_empty());
    assert!(select(&report, &selection(&report, recut)).is_err());
}
#[test]
fn raw_orientation_scope_ignores_ordered_reverse_flags_and_stale_or_forged_bindings_fail() {
    let (inputs, target, r) = fixture(A, B, "ATGA", "ATGA");
    let first = inspect(&inputs, &target, &r).unwrap();
    let mut reversed = r.clone();
    reversed.recipe.fragments[1].source.reverse = true;
    let reverse = inspect(&inputs, &target, &reversed).unwrap();
    let keys = |r: &Report| {
        r.pool
            .as_ref()
            .unwrap()
            .candidates
            .iter()
            .map(|c| (&c.directed_sequence_key, &c.physical_sequence_key))
            .map(|(a, b)| (a.clone(), b.clone()))
            .collect::<Vec<_>>()
    };
    assert_eq!(keys(&first), keys(&reverse));
    assert_ne!(first.context_revision, reverse.context_revision);
    assert!(
        reverse
            .pool
            .as_ref()
            .unwrap()
            .candidates
            .iter()
            .all(|c| c.path.iter().all(|p| !p.reverse))
    );
    let index = first
        .assessments
        .iter()
        .find(|a| a.target_matches)
        .unwrap()
        .candidate_index;
    let mut stale = selection(&first, index);
    stale.context_revision = "0".repeat(64);
    assert!(select(&first, &stale).is_err());
    let mut key = selection(&first, index);
    key.physical_sequence_key = "0".repeat(64);
    assert!(select(&first, &key).is_err());
    let mut forged = first.clone();
    forged.assessments[index].target_matches = !forged.assessments[index].target_matches;
    assert!(forged.validate_binding(&first.context_revision).is_err());
    let mut path = first.clone();
    path.assessments[index].source_roles = vec![31];
    assert!(path.validate_binding(&first.context_revision).is_err());
    let mut candidate = first.clone();
    candidate.pool.as_mut().unwrap().candidates[index].directed_sequence_key = "0".repeat(64);
    assert!(candidate.validate_binding(&first.context_revision).is_err());
    let mut stale_input = inputs;
    stale_input[0].revision = "0".repeat(64);
    assert!(inspect_revision(&stale_input, &target, &r).is_err());
}
