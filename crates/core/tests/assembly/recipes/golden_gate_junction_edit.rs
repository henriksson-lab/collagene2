use collagene_core::{
    chemistry::EndRepairMethod,
    cloning::{
        assembly::{
            intent::{self as assembly_intent, Contribution, Selector},
            junction_edit::Edit,
            recipe::{FragmentSpec, Preparation},
        },
        golden_gate::{GoldenGateEnzyme, junction_edit::*},
        ligation::recipe as ligation_recipe,
        restriction::catalog as catalog_digest,
    },
    planning::plan::PlanInput,
    sequence::analysis::Primer,
    verification::{self, evidence::Status},
};
use collagene_model::{Molecule, Topology};
const A: &str = "ATGACCCCGGGGTTTCCAGA";
const B: &str = "ATCCGTTAGCGACCTTGGAA";
fn input(id: &str, dna: &str) -> PlanInput {
    let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
    PlanInput {
        id: id.into(),
        revision: molecule.revision(),
        molecule,
    }
}
fn fixture(enzyme: GoldenGateEnzyme) -> (Vec<PlanInput>, Request) {
    let inputs = vec![input("a", A), input("b", B)];
    let p = enzyme.profile();
    let k = (p.bottom_cut - p.top_cut) as usize;
    let spacer = "A".repeat(p.top_cut as usize - p.recognition.len());
    let t = format!("ACGT{}{}", p.recognition, spacer);
    let closing = if k == 3 { "CAT" } else { "TCAT" };
    let retained = ligation_recipe::Fragment {
        source: FragmentSpec {
            generation: None,
            role: "retained".into(),
            input_id: "b".into(),
            reverse: false,
            preparation: Preparation::PcrCatalogDigest {
                conditions: None,
                forward: Primer::new("manual B forward", "ATCCGTTA", &t).unwrap(),
                reverse: Primer::new("manual B reverse", "TTCCAAGG", &format!("{t}{closing}"))
                    .unwrap(),
                selection: catalog_digest::select(
                    &catalog_digest::builtin_catalog().unwrap(),
                    &[p.name.to_ascii_lowercase()],
                )
                .unwrap(),
                fragment_index: 1,
                repair: None,
                template_strand: collagene_core::chemistry::Strand::Top,
            },
        },
        treatment: None,
    };
    let r = Request {
        base: assembly_intent::Request {
            material_bindings: vec![],
            name: "Independent inserted circle".into(),
            fragments: vec![
                Contribution {
                    role: "selected".into(),
                    input_id: "a".into(),
                    selector: Selector::ReferenceRange {
                        start: 0,
                        end: 20,
                        reverse: false,
                    },
                },
                Contribution {
                    role: "retained".into(),
                    input_id: "b".into(),
                    selector: Selector::ReferenceRange {
                        start: 0,
                        end: 20,
                        reverse: false,
                    },
                },
            ],
            circular: true,
        },
        edits: vec![Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "GC".into(),
        }],
        rows: vec![
            Row::DesignPcr {
                role: "selected".into(),
                input_id: "a".into(),
                reverse: false,
                forward_length: 8,
                reverse_length: 8,
            },
            Row::Retain(retained),
        ],
        enzyme,
        clamp: "ACGT".into(),
        spacer,
        terminal_right_overhang: None,
        allow_nicked: false,
    };
    (inputs, r)
}
fn fresh(inputs: &[PlanInput], report: &Report) -> verification::evidence::VerificationReport {
    let plan = ligation_recipe::compile(
        inputs,
        &report.intent.target,
        &report.recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap()
    .plan;
    verification::engine::verify(&plan, &report.intent.target)
}
#[test]
fn literal_inserted_circle_uses_real_bsai_or_sapi_ends_and_keeps_retained_oligos_exact() {
    for enzyme in [GoldenGateEnzyme::BsaI, GoldenGateEnzyme::SapI] {
        let (inputs, r) = fixture(enzyme);
        let original = serde_json::to_value(&inputs).unwrap();
        let retained = if let Row::Retain(f) = &r.rows[1] {
            serde_json::to_vec(f).unwrap()
        } else {
            panic!()
        };
        let p = design(&inputs, &r).unwrap();
        p.validate_binding(&design_revision(&inputs, &r).unwrap())
            .unwrap();
        assert_eq!(
            p.intent.target.sequence.as_str(),
            "ATGACCCCGGGGTTTCCAGAGCATCCGTTAGCGACCTTGGAA"
        );
        assert_eq!(p.retained_rows, vec![1]);
        assert!(p.designs[1].is_none());
        assert_eq!(
            serde_json::to_vec(&p.recipe.fragments[1]).unwrap(),
            retained
        );
        let d = p.designs[0].as_ref().unwrap();
        assert_eq!(d.oriented_top, "ATGACCCCGGGGTTTCCAGAGC");
        assert_eq!(d.reverse_tail.inserted, "GC");
        assert_eq!(d.reverse.tail, d.reverse_tail.sequence());
        assert_eq!(d.fragment_index, 1);
        if enzyme == GoldenGateEnzyme::BsaI {
            assert_eq!(d.reverse_tail.fusion, "GGAT");
            assert_eq!(d.left.overhang, "ATGA");
            assert_eq!(d.right.overhang, "GGAT");
            assert_eq!(d.reverse.tail, "ACGTGGTCTCAGGATGC");
        } else {
            assert_eq!(d.reverse_tail.fusion, "GAT");
            assert_eq!(d.left.overhang, "ATG");
            assert_eq!(d.right.overhang, "GAT");
            assert_eq!(d.reverse.tail, "ACGTGCTCTTCAGATGC");
        }
        assert_eq!(serde_json::to_value(&inputs).unwrap(), original);
        let report = fresh(&inputs, &p);
        assert_eq!(
            report.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            report.findings
        );
        assert_eq!(
            report.products.last().unwrap().sequence,
            p.intent.target.sequence
        );
    }
}
#[test]
fn reverse_row_places_nonpalindromic_insert_on_original_forward_primer_and_actual_bottom_projection()
 {
    let (inputs, mut r) = fixture(GoldenGateEnzyme::BsaI);
    let Selector::ReferenceRange {
        reverse: selected_reverse,
        ..
    } = &mut r.base.fragments[0].selector
    else {
        panic!("reference selector")
    };
    *selected_reverse = true;
    r.rows[0] = Row::DesignPcr {
        role: "selected".into(),
        input_id: "a".into(),
        reverse: true,
        forward_length: 8,
        reverse_length: 8,
    };
    r.edits[0].inserted = "AG".into();
    if let Row::Retain(row) = &mut r.rows[1] {
        if let Preparation::PcrCatalogDigest { reverse, .. } = &mut row.source.preparation {
            reverse.tail = "ACGTGGTCTCACAGA".into();
        }
    }
    let p = design(&inputs, &r).unwrap();
    assert_eq!(
        p.intent.target.sequence.as_str(),
        "TCTGGAAACCCCGGGGTCATAGATCCGTTAGCGACCTTGGAA"
    );
    let d = p.designs[0].as_ref().unwrap();
    assert_eq!(d.oriented_top, "TCTGGAAACCCCGGGGTCATAG");
    assert_eq!(d.forward_tail.fusion, "GGAT");
    assert_eq!(d.forward_tail.inserted, "CT");
    assert_eq!(d.forward.tail, "ACGTGGTCTCAGGATCT");
    assert_eq!(d.reverse_tail.fusion, "TCTG");
    assert!(d.reverse_tail.inserted.is_empty());
    assert_eq!(
        fresh(&inputs, &p).status,
        Status::ConsistentUnderAssumptions
    );
}
#[test]
fn trimmed_linear_intent_and_explicit_last_pcr_terminal_are_separate_from_retained_choices() {
    let (inputs, mut r) = fixture(GoldenGateEnzyme::BsaI);
    r.base.circular = false;
    r.rows[1] = Row::DesignPcr {
        role: "retained".into(),
        input_id: "b".into(),
        reverse: false,
        forward_length: 8,
        reverse_length: 8,
    };
    r.terminal_right_overhang = Some("ACGT".into());
    r.edits[0].remove_left = 2;
    r.edits[0].remove_right = 2;
    let p = design(&inputs, &r).unwrap();
    assert_eq!(
        p.intent.target.sequence.as_str(),
        "ATGACCCCGGGGTTTCCAGCCCGTTAGCGACCTTGGAA"
    );
    assert!(p.retained_rows.is_empty());
    assert_eq!(p.designs[0].as_ref().unwrap().right.overhang, "ACGG");
    assert_eq!(p.designs[1].as_ref().unwrap().right.overhang, "ACGT");
    assert_eq!(
        fresh(&inputs, &p).status,
        Status::ConsistentUnderAssumptions
    );
    let mut missing = r.clone();
    missing.terminal_right_overhang = None;
    assert!(
        preflight(
            &inputs
                .iter()
                .map(|i| (i.id.as_str(), &i.molecule))
                .collect::<Vec<_>>(),
            &missing
        )
        .is_err()
    );
    let (_, mut retained) = fixture(GoldenGateEnzyme::BsaI);
    retained.base.circular = false;
    assert!(design(&inputs, &retained).is_ok());
    retained.terminal_right_overhang = Some("ACGT".into());
    assert!(design(&inputs, &retained).is_err());
}
#[test]
fn retained_repairs_and_names_stay_unchanged_while_fresh_verification_can_reject_the_recipe() {
    let (inputs, mut r) = fixture(GoldenGateEnzyme::BsaI);
    if let Row::Retain(row) = &mut r.rows[1] {
        if let Preparation::PcrCatalogDigest { repair, .. } = &mut row.source.preparation {
            *repair = Some(EndRepairMethod::Blunt);
        }
    }
    let before = if let Row::Retain(row) = &r.rows[1] {
        serde_json::to_vec(row).unwrap()
    } else {
        panic!()
    };
    let p = design(&inputs, &r).unwrap();
    assert_eq!(serde_json::to_vec(&p.recipe.fragments[1]).unwrap(), before);
    assert_ne!(
        fresh(&inputs, &p).status,
        Status::ConsistentUnderAssumptions
    );
    let wrong = Molecule::dna(
        "wrong independent unedited target",
        &format!("{A}{B}"),
        Topology::Circular,
    )
    .unwrap();
    let plan = ligation_recipe::compile(
        &inputs,
        &wrong,
        &p.recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap()
    .plan;
    assert_ne!(
        verification::engine::verify(&plan, &wrong).status,
        Status::ConsistentUnderAssumptions
    );
}
#[test]
fn full_ordered_budget_allows_longer_inserts_but_rejects201nt_oligos_stale_pins_and_internal_sites()
{
    let (inputs, mut r) = fixture(GoldenGateEnzyme::BsaI);
    r.edits[0].inserted = "A".repeat(110);
    let p = design(&inputs, &r).unwrap();
    assert_eq!(
        p.designs[0]
            .as_ref()
            .unwrap()
            .reverse
            .ordered_sequence()
            .len(),
        133
    );
    assert_eq!(p.intent.target.sequence.len(), 150);
    assert_eq!(
        fresh(&inputs, &p).status,
        Status::ConsistentUnderAssumptions
    );
    r.edits[0].inserted = "A".repeat(178);
    assert!(
        preflight(
            &inputs
                .iter()
                .map(|i| (i.id.as_str(), &i.molecule))
                .collect::<Vec<_>>(),
            &r
        )
        .is_err()
    );
    let (inputs, r) = fixture(GoldenGateEnzyme::BsaI);
    let mut stale = inputs.clone();
    stale[0].revision = "0".repeat(64);
    assert!(design_revision(&stale, &r).is_err());
    let mut pin = r.clone();
    if let Row::Retain(row) = &mut pin.rows[1] {
        if let Preparation::PcrCatalogDigest { selection, .. } = &mut row.source.preparation {
            selection.catalog_content_hash = "0".repeat(64);
        }
    }
    assert!(design(&inputs, &pin).is_err());
    let mut internal_inputs = inputs.clone();
    internal_inputs[0] = input("a", "ATGACCCCGGTCTCTCCAGA");
    assert!(design(&internal_inputs, &r).is_err());
    let mut swapped = r.clone();
    swapped.rows.swap(0, 1);
    assert!(design(&inputs, &swapped).is_err());
    let p = design(&inputs, &r).unwrap();
    let mut forged = p.clone();
    forged.designs[0].as_mut().unwrap().reverse_tail.inserted = "AT".into();
    assert!(forged.validate_binding(&p.context_revision).is_err());
    let mut renamed = p.clone();
    if let Preparation::PcrCatalogDigest { forward, .. } =
        &mut renamed.recipe.fragments[1].source.preparation
    {
        forward.name = "changed retained name".into();
    }
    assert!(renamed.validate_binding(&p.context_revision).is_err());
    // Existing old DTO/public design still uses empty literal inserts and its
    // original four-part segments; the helper refactor does not change that path.
    let target = Molecule::dna("old independent target", A, Topology::Linear).unwrap();
    let old = collagene_core::cloning::golden_gate::primers::Request {
        name: "Old empty-insert regression".into(),
        fragments: vec![collagene_core::cloning::golden_gate::primers::Fragment {
            role: "selected".into(),
            input_id: "a".into(),
            start: 0,
            end: 20,
            reverse: false,
            forward_length: 8,
            reverse_length: 8,
            right_overhang: "ATCC".into(),
        }],
        enzyme: GoldenGateEnzyme::BsaI,
        clamp: "ACGT".into(),
        spacer: "A".into(),
    };
    let report =
        collagene_core::cloning::golden_gate::primers::design(&inputs[..1], &target, &old).unwrap();
    assert_eq!(report.rows[0].oriented_top, A);
    assert_eq!(report.rows[0].reverse_tail.fusion, "GGAT");
    assert_eq!(report.rows[0].reverse.tail, "ACGTGGTCTCAGGAT");
    assert!(matches!(
        report.fragments[0].preparation,
        Preparation::PcrCatalogDigest { repair: None, .. }
    ));
}
