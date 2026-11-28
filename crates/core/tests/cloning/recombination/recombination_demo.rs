//! The user-facing import fixture retains two independent expected components.
use collagene_core::formats::sequence as formats;
use collagene_model::Topology;
use collagene_recombination::{
    Input, Limits, Reaction, TargetPolicy, check_batch_revision, check_target,
    manufacturer_context_profile, recombine,
};

#[test]
fn bp_demo_imports_and_matches_both_independently_specified_components() {
    let imported = formats::import_genbank(include_str!(
        "../../../../../fixtures/recombination-bp-demo.gb"
    ))
    .unwrap();
    assert_eq!(imported.molecules.len(), 4);
    let molecules = imported.molecules;
    assert_eq!(
        molecules
            .iter()
            .map(|m| m.sequence.len())
            .collect::<Vec<_>>(),
        [59, 4761, 2553, 2267]
    );
    assert_eq!(molecules[2].topology, Topology::Circular);
    assert_eq!(molecules[3].topology, Topology::Linear);
    let insert = molecules[0]
        .features
        .iter()
        .find(|f| f.label() == "synthetic insert")
        .unwrap();
    assert_eq!(
        insert.location.extract(&molecules[0].sequence).unwrap(),
        "ATGCCGTAA"
    );
    let target_insert = molecules[2]
        .features
        .iter()
        .find(|f| f.label() == "synthetic insert")
        .unwrap();
    assert_eq!(
        target_insert
            .location
            .extract(&molecules[2].sequence)
            .unwrap(),
        "ATGCCGTAA"
    );
    let inputs = vec![
        Input {
            role: "insert".into(),
            molecule: molecules[0].clone(),
        },
        Input {
            role: "donor".into(),
            molecule: molecules[1].clone(),
        },
    ];
    let profile = manufacturer_context_profile().unwrap();
    let limits = Limits::default();
    let report = recombine(&inputs, &profile, Reaction::BP, &limits).unwrap();
    assert!(report.complete);
    assert_eq!(report.candidates.len(), 1);
    assert_eq!(report.candidates[0].products.len(), 2);
    for target in &molecules[2..] {
        let check = check_target(
            &inputs,
            &profile,
            Reaction::BP,
            &limits,
            target,
            TargetPolicy {
                allow_reverse_complement: false,
            },
        )
        .unwrap();
        assert!(check.search_complete);
        assert_eq!(
            check.batch_revision,
            check_batch_revision(
                &inputs,
                &profile,
                Reaction::BP,
                &limits,
                target,
                TargetPolicy {
                    allow_reverse_complement: false
                },
            )
            .unwrap()
        );
        assert_eq!(check.matches.len(), 1);
        assert!(!check.no_match_proven_under_model);
    }
}
