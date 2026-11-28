use collagene_core::{
    analysis::Primer,
    assembly_intent::{self, Contribution, Selector},
    assembly_junction_edit::Edit,
    assembly_recipe::{self, FragmentSpec, Preparation},
    catalog_digest,
    primers::DesignOptions,
    selective_junction_edit::*,
    verification::{self, PlanInput, Status},
};
use collagene_model::{Molecule, Topology};
const A: &str = "ATGACCCCGGGGTTTCCAGA";
const B: &str = "ATCCGTTAGCGACCTTGGAA";
const C: &str = "GATATCGTCGATGCTAGCGGTACCTTAGC";
fn input(id: &str, dna: &str) -> PlanInput {
    let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fixture() -> (Vec<PlanInput>, Request) {
    let inputs = vec![input("a", A), input("b", B), input("c", C)];
    let base = assembly_intent::Request {
        material_bindings: vec![],
        name: "Selective mixed preparation".into(),
        fragments: vec![
            Contribution {
                role: "first".into(),
                input_id: "a".into(),
                selector: Selector::ReferenceRange {
                    start: 0,
                    end: 20,
                    reverse: false,
                },
            },
            Contribution {
                role: "manual".into(),
                input_id: "b".into(),
                selector: Selector::ReferenceRange {
                    start: 0,
                    end: 20,
                    reverse: false,
                },
            },
            Contribution {
                role: "cut".into(),
                input_id: "c".into(),
                selector: Selector::ReferenceRange {
                    start: 3,
                    end: 29,
                    reverse: false,
                },
            },
        ],
        circular: false,
    };
    let manual = FragmentSpec {
        generation: None,
        role: "manual".into(),
        input_id: "b".into(),
        preparation: Preparation::ExactPcr {
            forward: Primer::new("operator named B forward", "ATCCGTTA", "").unwrap(),
            reverse: Primer::new("operator named B reverse", "TTCCAAGG", "ACCGCTAGCATCGACGAT")
                .unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        },
        reverse: false,
    };
    let cut = FragmentSpec {
        generation: None,
        role: "cut".into(),
        input_id: "c".into(),
        preparation: Preparation::CatalogDigest {
            selection: catalog_digest::select(
                &catalog_digest::builtin_catalog().unwrap(),
                &["ecorv".into()],
            )
            .unwrap(),
            fragment_index: 1,
            repair: None,
        },
        reverse: false,
    };
    let r = Request {
        base,
        edits: vec![Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "GC".into(),
        }],
        rows: vec![
            Row::DesignPcr {
                role: "first".into(),
                input_id: "a".into(),
                reverse: false,
            },
            Row::Retain(manual),
            Row::Retain(cut),
        ],
        overlaps: vec![15, 18],
        primer_options: DesignOptions {
            min_length: 8,
            max_length: 8,
            min_gc_fraction: 0.,
            max_gc_fraction: 1.,
            max_tm_difference_c: 100.,
            max_homopolymer: 200,
            max_three_prime_complementarity: 200,
            ..Default::default()
        },
        profile_options: Default::default(),
    };
    (inputs, r)
}
fn fresh(inputs: &[PlanInput], report: &Report) -> verification::VerificationReport {
    let plan = assembly_recipe::compile(inputs, &report.intent.target, &report.recipe)
        .unwrap()
        .plan;
    verification::verify(&plan, &report.intent.target)
}
#[test]
fn one_selected_row_preserves_manual_and_catalog_bytes_and_matches_literal_independent_target() {
    let (inputs, r) = fixture();
    let original = serde_json::to_value(&inputs).unwrap();
    let retained = r
        .rows
        .iter()
        .skip(1)
        .map(|r| {
            if let Row::Retain(s) = r {
                serde_json::to_vec(s).unwrap()
            } else {
                panic!()
            }
        })
        .collect::<Vec<_>>();
    let report = design(&inputs, &r).unwrap();
    report
        .validate_binding(&design_revision(&inputs, &r).unwrap())
        .unwrap();
    assert_eq!(
        report.intent.target.sequence.as_str(),
        "ATGACCCCGGGGTTTCCAGAGCATCCGTTAGCGACCTTGGAAATCGTCGATGCTAGCGGTACCTTAGC"
    );
    assert_eq!(report.retained_rows, vec![1, 2]);
    assert!(report.designs[0].is_some());
    assert!(report.designs[1].is_none());
    assert!(report.designs[2].is_none());
    assert_eq!(
        serde_json::to_vec(&report.recipe.fragments[1]).unwrap(),
        retained[0]
    );
    assert_eq!(
        serde_json::to_vec(&report.recipe.fragments[2]).unwrap(),
        retained[1]
    );
    assert_eq!(
        report.designs[0].as_ref().unwrap().pairs[0]
            .reverse
            .primer
            .tail,
        "AGGTCGCTAACGGATGC"
    );
    assert_eq!(serde_json::to_value(&inputs).unwrap(), original);
    let checked = fresh(&inputs, &report);
    assert_eq!(
        checked.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        checked.findings
    );
    assert_eq!(
        checked.products.last().unwrap().sequence,
        report.intent.target.sequence
    );
}
#[test]
fn incompatible_manual_row_is_kept_literal_and_only_fresh_full_reconstruction_reports_mismatch() {
    let (inputs, mut r) = fixture();
    if let Row::Retain(s) = &mut r.rows[1] {
        if let Preparation::ExactPcr { reverse, .. } = &mut s.preparation {
            reverse.tail = "A".repeat(18);
        }
    }
    let before = serde_json::to_vec(&r.rows[1]).unwrap();
    let report = design(&inputs, &r).unwrap();
    assert_eq!(serde_json::to_vec(&report.request.rows[1]).unwrap(), before);
    assert!(report.designs[1].is_none());
    assert_ne!(
        fresh(&inputs, &report).status,
        Status::ConsistentUnderAssumptions
    );
    assert!(report.warnings.iter().any(|w| w.contains("incompatible")));
    let mut wrong = report.intent.target.clone();
    wrong.sequence = collagene_model::Sequence::dna(&"A".repeat(68)).unwrap();
    let plan = assembly_recipe::compile(&inputs, &wrong, &report.recipe)
        .unwrap()
        .plan;
    assert_ne!(
        verification::verify(&plan, &wrong).status,
        Status::ConsistentUnderAssumptions
    );
}
#[test]
fn reverse_and_circular_variable_homologies_change_only_explicitly_selected_originals() {
    let (inputs, mut r) = fixture();
    let Selector::ReferenceRange {
        reverse: selected_reverse,
        ..
    } = &mut r.base.fragments[0].selector
    else {
        panic!("reference selector")
    };
    *selected_reverse = true;
    r.rows[0] = Row::DesignPcr {
        role: "first".into(),
        input_id: "a".into(),
        reverse: true,
    };
    let report = design(&inputs, &r).unwrap();
    assert_eq!(
        report.intent.target.sequence.as_str(),
        "TCTGGAAACCCCGGGGTCATGCATCCGTTAGCGACCTTGGAAATCGTCGATGCTAGCGGTACCTTAGC"
    );
    assert_eq!(
        report.designs[0].as_ref().unwrap().pairs[0]
            .forward
            .primer
            .tail,
        "AGGTCGCTAACGGATGC"
    );
    assert!(
        report.designs[0].as_ref().unwrap().pairs[0]
            .reverse
            .primer
            .tail
            .is_empty()
    );
    assert_eq!(
        fresh(&inputs, &report).status,
        Status::ConsistentUnderAssumptions
    );
    r.base.circular = true;
    r.overlaps.push(8);
    r.rows[2] = Row::DesignPcr {
        role: "cut".into(),
        input_id: "c".into(),
        reverse: false,
    };
    let report = design(&inputs, &r).unwrap();
    assert_eq!(report.recipe.overlaps, vec![15, 18, 8]);
    assert_eq!(report.retained_rows, vec![1]);
    assert_eq!(
        report.designs[2].as_ref().unwrap().pairs[0]
            .reverse
            .primer
            .tail,
        "TTTCCAGA"
    );
    assert_eq!(
        fresh(&inputs, &report).status,
        Status::ConsistentUnderAssumptions
    );
}
#[test]
fn retained_manual_domains_do_not_inherit_selected_primer_length_constraints() {
    let (inputs, mut r) = fixture();
    if let Row::Retain(s) = &mut r.rows[1] {
        s.preparation = Preparation::ExactPcr {
            forward: Primer::new("manual four-base F", "ATCC", "").unwrap(),
            reverse: Primer::new("manual four-base R", "TTCC", "ACCGCTAGCATCGACGAT").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        };
    }
    let report = design(&inputs, &r).unwrap();
    let Preparation::ExactPcr {
        forward, reverse, ..
    } = &report.recipe.fragments[1].preparation
    else {
        panic!()
    };
    assert_eq!(forward.annealing.len(), 4);
    assert_eq!(reverse.annealing.len(), 4);
    assert!(report.designs[1].is_none());
    assert_eq!(
        fresh(&inputs, &report).status,
        Status::ConsistentUnderAssumptions
    );
    // All-retained proposals keep declarations even if primer-design parameters
    // would be unsuitable for new designs; they never fabricate positive evidence.
    let mut retained = r.clone();
    retained.rows = report
        .recipe
        .fragments
        .iter()
        .cloned()
        .map(Row::Retain)
        .collect();
    retained.primer_options.min_length = 1;
    retained.primer_options.max_length = 1;
    let proposal = design(&inputs, &retained).unwrap();
    assert!(proposal.designs.iter().all(Option::is_none));
    assert_eq!(proposal.retained_rows, vec![0, 1, 2]);
}
#[test]
fn role_order_stale_sources_unbounded_work_and_receipt_tampering_fail_before_application() {
    let (inputs, r) = fixture();
    let borrowed = inputs
        .iter()
        .map(|i| (i.id.as_str(), &i.molecule))
        .collect::<Vec<_>>();
    let mut swapped = r.clone();
    swapped.rows.swap(0, 1);
    assert!(preflight(&borrowed, &swapped).is_err());
    let mut wrong = r.clone();
    if let Row::DesignPcr { role, .. } = &mut wrong.rows[0] {
        *role = "different".into();
    }
    assert!(preflight(&borrowed, &wrong).is_err());
    let mut stale = inputs.clone();
    stale[0].revision = "0".repeat(64);
    assert!(design_revision(&stale, &r).is_err());
    let mut wide = r.clone();
    wide.primer_options.max_length = 71;
    assert!(preflight(&borrowed, &wide).is_err());
    let mut overflow = r.clone();
    overflow.overlaps[0] = 10_000;
    assert!(preflight(&borrowed, &overflow).is_err());
    let report = design(&inputs, &r).unwrap();
    let mut tampered = report.clone();
    if let Preparation::ExactPcr { reverse, .. } = &mut tampered.recipe.fragments[1].preparation {
        reverse.name = "silently renamed".into();
    }
    assert!(tampered.validate_binding(&report.context_revision).is_err());
    let mut indices = report.clone();
    indices.retained_rows = vec![2, 1];
    assert!(indices.validate_binding(&report.context_revision).is_err());
    let mut metric = report.clone();
    metric.designs[0].as_mut().unwrap().pairs[0]
        .reverse
        .analysis
        .tail = "AC".into();
    assert!(metric.validate_binding(&report.context_revision).is_err());
}
