use collagene_core::{
    chemistry::{self, FivePrime, Terminus, ThreePrime, tracked_duplex as modified_duplex},
    cloning::{assembly::recipe as a, overlap::profiles::ProfileOptions},
    formats::sequence as formats,
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    verification::{self as v, evidence::Status},
};
use collagene_model::{Molecule, Topology};
fn molecule(name: &str, dna: &str) -> Molecule {
    Molecule::dna(name, dna, Topology::Linear).unwrap()
}
fn input(id: &str, m: Molecule) -> PlanInput {
    PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn terminus(five_prime: FivePrime, three_prime: ThreePrime) -> Terminus {
    Terminus {
        five_prime,
        three_prime,
    }
}
fn declaration(id: &str, t: Terminus) -> Operation {
    Operation::DeclareSourceTermini {
        template: Source::Input(id.into()),
        left: t,
        right: t,
    }
}
fn plan(a: Terminus, b: Terminus, allow_nicked: bool) -> (DesignPlan, Molecule) {
    let target = molecule("independent target", "AAAACCCC");
    (
        DesignPlan {
            name: "explicit original source hypotheses".into(),
            target_revision: target.revision(),
            inputs: vec![
                input("a", molecule("first original", "AAAA")),
                input("b", molecule("second original", "CCCC")),
            ],
            steps: vec![
                declaration("a", a),
                declaration("b", b),
                Operation::Ligate {
                    fragments: vec![
                        OrientedSource {
                            source: Source::Step(0),
                            reverse: false,
                        },
                        OrientedSource {
                            source: Source::Step(1),
                            reverse: false,
                        },
                    ],
                    circular: false,
                    allow_nicked,
                    options: Default::default(),
                },
            ],
        },
        target,
    )
}
#[test]
fn declared_phosphate_oh_seals_while_oh_nicks_unknown_and_blocked_remain_explicit() {
    let known = terminus(FivePrime::Phosphate, ThreePrime::Hydroxyl);
    let hydroxyl = terminus(FivePrime::Hydroxyl, ThreePrime::Hydroxyl);
    let (p, t) = plan(known, known, false);
    let r = v::engine::verify(&p, &t);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    let physical: chemistry::DuplexSummary = serde_json::from_str(
        &r.product(2)
            .expect("selected step must reconstruct")
            .metadata["duplex_model"],
    )
    .unwrap();
    assert!(physical.nicks.is_empty());
    assert!(
        r.findings
            .iter()
            .any(|s| s.message.contains("operator hypothesis"))
    );
    let (p, t) = plan(known, hydroxyl, true);
    let r = v::engine::verify(&p, &t);
    assert_eq!(r.status, Status::ConsistentUnderAssumptions);
    let physical: chemistry::DuplexSummary = serde_json::from_str(
        &r.product(2)
            .expect("selected step must reconstruct")
            .metadata["duplex_model"],
    )
    .unwrap();
    assert_eq!(physical.nicks.len(), 1);
    let (p, t) = plan(hydroxyl, hydroxyl, true);
    assert_eq!(v::engine::verify(&p, &t).status, Status::Invalid);
    let (p, t) = plan(
        terminus(FivePrime::Phosphate, ThreePrime::Unknown),
        known,
        true,
    );
    assert_eq!(v::engine::verify(&p, &t).status, Status::Incomplete);
    let blocked = terminus(FivePrime::Phosphate, ThreePrime::Blocked);
    let (p, t) = plan(blocked, blocked, true);
    assert_eq!(v::engine::verify(&p, &t).status, Status::Invalid);
}
#[test]
fn dna_annotations_citations_and_unknown_modification_state_are_preserved() {
    let raw = "LOCUS       reference                 12 bp    DNA     linear   SYN 01-JAN-2026\nDEFINITION  Literal original source.\nACCESSION   AB123456\nVERSION     AB123456.2\nREFERENCE   1  (bases 3 to 8; 10 to 12)\n  AUTHORS   Example,A.\n  TITLE     Literal citation\n  JOURNAL   Unpublished\nFEATURES             Location/Qualifiers\n     misc_feature    3..8\n                     /label=\"literal feature\"\nORIGIN\n        1 aaaaccccgg gg\n//\n";
    let m = formats::import_genbank(raw).unwrap().molecules.remove(0);
    let target = molecule("independent unchanged sequence", m.sequence.as_str());
    let known = Terminus::CLEAVED;
    let plan = DesignPlan {
        name: "hypothesis only".into(),
        target_revision: target.revision(),
        inputs: vec![input("source", m.clone())],
        steps: vec![declaration("source", known)],
    };
    let report = v::engine::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    let p = report.product(0).expect("selected step must reconstruct");
    assert_eq!(p.sequence, m.sequence);
    assert_eq!(p.features, m.features);
    assert_eq!(p.original_genbank, m.original_genbank);
    assert!(
        formats::export_genbank(p)
            .unwrap()
            .contains("REFERENCE   1")
    );
    let mods: modified_duplex::ModificationSummary =
        serde_json::from_str(&p.metadata["tracked_modifications"]).unwrap();
    assert!(mods.unknown_eligible > 0);
    assert_eq!(mods.known_unmodified, 0);
    assert_eq!(mods.known_modified, 0);
    assert!(mods.generation.is_empty());
    assert!(v::sources::declared_source_termini_preflight(p).is_err());
}
#[test]
fn descendants_circles_and_restrictive_markers_cannot_be_promoted_by_declaration() {
    let known = Terminus::CLEAVED;
    let (mut p, t) = plan(known, known, false);
    p.steps = vec![
        Operation::TreatEnds {
            template: Source::Input("a".into()),
            treatment: chemistry::EndTreatment::Phosphorylate,
        },
        Operation::DeclareSourceTermini {
            template: Source::Step(0),
            left: known,
            right: known,
        },
    ];
    let r = v::engine::verify(&p, &t);
    assert_eq!(r.status, Status::Incomplete);
    assert_eq!(r.reconstructed_count(), 1);
    for (key, value) in [
        ("duplex_model", "{}"),
        ("source:legacy:chemistry_unmapped", "false"),
        ("variant:reference_only", "false"),
        ("recombination:reference_only", "false"),
        ("topo:physical_graph_required", "false"),
        ("tracked_modifications", "{}"),
        ("source:restriction_catalog:selection", "{}"),
        ("genbank.molecule_type", "ss-DNA"),
    ] {
        let mut m = molecule("untrusted imported source", "AAAA");
        m.metadata.insert(key.into(), value.into());
        assert!(
            v::sources::declared_source_termini_preflight(&m).is_err(),
            "{key}"
        );
    }
    let mut circle = molecule("circle", "AAAA");
    circle.topology = Topology::Circular;
    assert!(v::sources::declared_source_termini_preflight(&circle).is_err());
    let mut raw = molecule("raw marker", "AAAA");
    raw.original_genbank=Some("LOCUS       raw 4 bp DNA linear\nCOMMENT     Collagene physical product v1: unverified sequence-only evidence\nORIGIN\n        1 aaaa\n//\n".into());
    assert!(v::sources::declared_source_termini_preflight(&raw).is_err());
}
#[test]
fn preparation_and_operation_are_append_only_and_whole_linear_is_not_upgraded() {
    let inputs = vec![
        input("a", molecule("source a", "AAAT")),
        input("b", molecule("source b", "TGATTT")),
    ];
    let target = molecule("independent overlap target", "AAATGATTT");
    let known = Terminus::CLEAVED;
    let recipe = a::Recipe {
        name: "declaration reuse".into(),
        fragments: vec![
            a::FragmentSpec {
                generation: None,
                role: "left".into(),
                input_id: "a".into(),
                preparation: a::Preparation::DeclaredDuplex {
                    left: known,
                    right: known,
                },
                reverse: false,
            },
            a::FragmentSpec {
                generation: None,
                role: "right".into(),
                input_id: "b".into(),
                preparation: a::Preparation::WholeLinear,
                reverse: false,
            },
        ],
        overlaps: vec![1],
        circular: false,
        options: ProfileOptions::default(),
    };
    let preview = a::preview(&inputs, &target, &recipe).unwrap();
    assert_eq!(preview.compiled.rows[0].preparation_step, Some(0));
    assert!(preview.compiled.rows[1].preparation_step.is_none());
    assert!(
        matches!(preview.compiled.plan.steps[0],Operation::DeclareSourceTermini{template:Source::Input(ref id),..}if id=="a")
    );
    assert_eq!(
        postcard::to_allocvec(&preview.compiled.plan.steps[0]).unwrap()[0],
        16
    );
    assert_eq!(
        postcard::to_allocvec(&recipe.fragments[0].preparation).unwrap()[0],
        3
    );
    assert!(
        !preview.prepared[1]
            .as_ref()
            .expect("prepared row must reconstruct")
            .metadata
            .contains_key("source:declared_termini")
    );
}
