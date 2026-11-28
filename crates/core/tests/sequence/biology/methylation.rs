//! Hand-specified source-context and provenance fixtures; no scientific executable.
use collagene_core::{methylation::*, *};
use collagene_model::{Error, Molecule, Topology};
fn mol(s: &str) -> Molecule {
    Molecule::dna("source", s, Topology::Linear).unwrap()
}
fn opts(name: &str) -> MethylationOptions {
    MethylationOptions {
        enzymes: vec![name.into()],
        dam: MethylationState::FullyModified,
        dcm: MethylationState::FullyModified,
        cpg: MethylationState::FullyModified,
        ..Default::default()
    }
}
#[test]
fn all_starter_identities_have_explicit_versioned_rules_without_formulation_aliases() {
    let ps = profiles();
    assert_eq!(ps.len(), analysis::starter_enzymes().len());
    let cases = [
        ("BamHI", [Sensitivity::NotSensitive; 3]),
        ("HindIII", [Sensitivity::NotSensitive; 3]),
        ("PstI", [Sensitivity::NotSensitive; 3]),
        ("BbsI", [Sensitivity::NotSensitive; 3]),
        ("SapI", [Sensitivity::NotSensitive; 3]),
        (
            "NotI",
            [
                Sensitivity::NotSensitive,
                Sensitivity::NotSensitive,
                Sensitivity::Blocked,
            ],
        ),
        (
            "SmaI",
            [
                Sensitivity::NotSensitive,
                Sensitivity::NotSensitive,
                Sensitivity::Blocked,
            ],
        ),
        (
            "EcoRI",
            [
                Sensitivity::NotSensitive,
                Sensitivity::NotSensitive,
                Sensitivity::BlockedBySomeOverlaps,
            ],
        ),
        (
            "EcoRV",
            [
                Sensitivity::NotSensitive,
                Sensitivity::NotSensitive,
                Sensitivity::ImpairedBySomeOverlaps,
            ],
        ),
        (
            "XbaI",
            [
                Sensitivity::BlockedByOverlapping,
                Sensitivity::NotSensitive,
                Sensitivity::NotSensitive,
            ],
        ),
        (
            "BsaI",
            [
                Sensitivity::NotSensitive,
                Sensitivity::Unsupported,
                Sensitivity::BlockedBySomeOverlaps,
            ],
        ),
        ("BsmBI", [Sensitivity::Unsupported; 3]),
    ];
    for (name, expected) in cases {
        let p = ps.iter().find(|p| p.enzyme == name).unwrap();
        assert_eq!(p.sensitivities, expected);
        assert!(!p.sources.is_empty());
    }
    assert!(
        ps.iter()
            .find(|p| p.enzyme == "BsmBI")
            .unwrap()
            .supplier_identity
            .contains("Generic")
    );
    assert!(
        ps.iter()
            .find(|p| p.enzyme == "BsaI")
            .unwrap()
            .supplier_identity
            .contains("R0535")
    );
    assert_eq!(profile_revision().len(), 64);
}
#[test]
fn xbai_context_is_site_specific_and_not_a_global_enzyme_block() {
    let source = mol("AATCTAGACCAAGATCTAGACCAATCTAGATCAA");
    let r = assess(&source, &opts("XbaI")).unwrap();
    assert_eq!(r.sites.len(), 3);
    assert_eq!(
        r.sites.iter().map(|s| s.outcome).collect::<Vec<_>>(),
        vec![
            Outcome::PermittedUnderProfile,
            Outcome::Blocked,
            Outcome::Blocked
        ]
    );
    let o = &r.sites[1].effects[0].overlaps[0];
    assert_eq!(o.motif, "GATC");
    assert_eq!(o.positions, vec![12, 13, 14, 15]);
    assert_eq!(
        o.modified_bases,
        vec![
            ModifiedBase {
                position: 13,
                reverse: false,
                base: 'A'
            },
            ModifiedBase {
                position: 14,
                reverse: true,
                base: 'A'
            }
        ]
    );
    r.validate_bindings_for(&source, &opts("XbaI")).unwrap();
    replay(&source, &opts("XbaI"), &r).unwrap();
}
#[test]
fn insensitive_sites_do_not_become_blocked_merely_because_a_motif_overlaps() {
    for name in ["BamHI", "HindIII", "PstI", "BbsI", "SapI"] {
        let p = profiles().into_iter().find(|p| p.enzyme == name).unwrap();
        let source = mol(&format!("AAGATCCWGGCG{}AAGATCCWGGCG", p.recognition));
        let r = assess(&source, &opts(name)).unwrap();
        assert!(
            r.sites
                .iter()
                .filter(|s| s.certain)
                .all(|s| s.outcome == Outcome::PermittedUnderProfile)
        );
    }
    let source = mol("AAGGATCCAA");
    let r = assess(&source, &opts("BamHI")).unwrap();
    assert_eq!(r.sites[0].effects[0].overlaps.len(), 1);
    assert_eq!(r.sites[0].outcome, Outcome::PermittedUnderProfile);
}
#[test]
fn incomplete_context_and_unsupported_hemi_mixed_states_never_silently_permit() {
    let source = mol("AATCTAGANC");
    let r = assess(&source, &opts("XbaI")).unwrap();
    assert_eq!(r.sites[0].outcome, Outcome::Unknown);
    assert!(!r.sites[0].effects[0].overlaps[0].certain);
    let source = mol("AATCTAGATCAA");
    for s in [MethylationState::Unknown, MethylationState::HemiOrMixed] {
        let mut o = opts("XbaI");
        o.dam = s;
        assert_eq!(
            assess(&source, &o).unwrap().sites[0].outcome,
            Outcome::Unknown
        );
    }
    let mut o = opts("XbaI");
    o.dam = MethylationState::Unmodified;
    assert_eq!(
        assess(&source, &o).unwrap().sites[0].outcome,
        Outcome::PermittedUnderProfile
    );
}
#[test]
fn some_combination_profiles_retain_unknown_and_generic_bsai_conflict_is_disclosed() {
    for (name, sequence, category) in [
        ("EcoRI", "CGAATTCG", 2),
        ("EcoRV", "CGATATCG", 2),
        ("BsaI", "GGTCTCCAGGAAAAAA", 1),
        ("BsaI", "GGTCTCGAAAAAA", 2),
    ] {
        let r = assess(&mol(sequence), &opts(name)).unwrap();
        assert_eq!(r.sites.len(), 1);
        assert!(!r.sites[0].effects[category].overlaps.is_empty());
        assert_eq!(r.sites[0].effects[category].outcome, Outcome::Unknown);
    }
    let r = assess(&mol("ACGTCTCAAAAAA"), &opts("BsmBI")).unwrap();
    assert_eq!(r.sites[0].outcome, Outcome::Unknown);
}
#[test]
fn recognition_certainty_and_available_cuts_are_separate_from_methylation() {
    let r = assess(&mol("AAGGTCTC"), &opts("BsaI")).unwrap();
    assert!(r.sites[0].certain);
    assert_eq!(r.sites[0].top_cut, None);
    assert_eq!(r.sites[0].bottom_cut, None);
    let mut o = opts("BamHI");
    o.dam = MethylationState::Unmodified;
    let r = assess(&mol("AAGGNTCCAA"), &o).unwrap();
    assert!(!r.sites[0].certain);
    assert_eq!(r.sites[0].outcome, Outcome::Unknown);
}
#[test]
fn bounded_work_and_source_preflight_do_not_assert_absent_sites() {
    let source = mol("TCTAGATC".repeat(40).as_str());
    let mut o = opts("XbaI");
    o.max_comparisons = 1;
    let r = assess(&source, &o).unwrap();
    assert!(!r.complete && !r.limits.is_empty());
    assert!(r.sites.is_empty());
    let mut o = opts("XbaI");
    o.max_sites = 1;
    let r = assess(&source, &o).unwrap();
    assert!(!r.complete);
    assert_eq!(r.sites.len(), 1);
    let mut too_big = mol("AAAA");
    too_big.original_genbank = Some("x".repeat(8 * 1024 * 1024));
    assert!(matches!(batch_revision(&too_big, &o), Err(Error::TooLarge)));
    let mut o = opts("XbaI");
    o.enzymes.push("xbai".into());
    assert!(source_preflight(&source, &o).is_err());
}
#[test]
fn stale_source_options_profiles_and_forged_evidence_are_rejected_without_imported_authority() {
    let source = mol("AATCTAGATCAA");
    let o = opts("XbaI");
    let r = assess(&source, &o).unwrap();
    let mut stale = source.clone();
    stale.name.push('!');
    assert!(r.validate_bindings_for(&stale, &o).is_err());
    let mut wrong = o.clone();
    wrong.dam = MethylationState::Unmodified;
    assert!(r.validate_bindings_for(&source, &wrong).is_err());
    let mut bad = r.clone();
    bad.profiles[0].sources[0] = "https://example.com/forged".into();
    bad.validate().unwrap();
    assert!(bad.validate_bindings_for(&source, &o).is_err());
    let mut forged = r.clone();
    forged.sites[0].top_cut = Some(0);
    forged.validate_bindings_for(&source, &o).unwrap();
    assert!(replay(&source, &o, &forged).is_err());
    let mut old = r.clone();
    old.engine_version = "historical-methylation".into();
    old.validate().unwrap();
    assert!(old.validate_bindings_for(&source, &o).is_err());
    let mut value = serde_json::to_value(o).unwrap();
    value["infer_from_host"] = serde_json::json!(true);
    assert!(serde_json::from_value::<MethylationOptions>(value).is_err());
}
