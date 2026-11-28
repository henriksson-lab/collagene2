use collagene_core::sequence::codon_design::*;
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
fn protein(s: &str) -> Molecule {
    Molecule {
        name: "required peptide".into(),
        sequence: Sequence::new(Alphabet::Protein, s).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn dna(s: &str) -> Molecule {
    Molecule::dna("candidate", s, Topology::Linear).unwrap()
}
fn chosen(r: &DesignReport) -> &str {
    r.candidate.as_ref().unwrap().sequence.as_str()
}
#[test]
fn reverse_translation_is_exact_deterministic_and_retains_declared_provenance() {
    let p = protein("MK*");
    let o = DesignOptions::default();
    let r = design(&p, &o).unwrap();
    assert_eq!(chosen(&r), "ATGAAATAA");
    assert_eq!(r.target_protein, "MK*");
    assert!(r.search_complete && r.optimality_proven);
    assert_eq!(r.score, Some(3));
    r.validate_for_source(&p, &o).unwrap();
    let c = r.candidate.as_ref().unwrap();
    assert_eq!(c.metadata["codon_design:required_protein"], "MK*");
    assert_eq!(c.metadata["codon_design:source_revision"], p.revision());
    assert_eq!(
        serde_json::from_str::<DesignOptions>(&c.metadata["codon_design:options"]).unwrap(),
        o
    );
    assert_eq!(design(&p, &o).unwrap(), r);
}
#[test]
fn mitochondrial_and_bacterial_decoding_is_ordinary_not_initiation_reencoding() {
    let o = DesignOptions {
        genetic_code: 2,
        ..Default::default()
    };
    let r = design(&dna("ATAAGA"), &o).unwrap();
    assert_eq!(r.target_protein, "M*");
    assert_eq!(chosen(&r), "ATAAGA");
    let o = DesignOptions {
        genetic_code: 11,
        ..Default::default()
    };
    let r = design(&dna("GTGAAATAG"), &o).unwrap();
    assert_eq!(r.target_protein, "VK*");
    assert!(chosen(&r).starts_with("GTA"));
    assert!(
        design(&dna("ATAAGA"), &DesignOptions::default())
            .unwrap()
            .target_protein
            == "IR"
    );
}
#[test]
fn weights_are_explicit_scores_and_zero_disallows_without_fitness_claims() {
    let p = protein("KK");
    let o = DesignOptions {
        preferences: collagene_core::sequence::codon_design::CodonPreferences::Explicit {
            weights: vec![
                CodonWeight {
                    codon: "AAA".into(),
                    weight: 0,
                },
                CodonWeight {
                    codon: "AAG".into(),
                    weight: 17,
                },
            ],
        },
        ..Default::default()
    };
    let r = design(&p, &o).unwrap();
    assert_eq!(chosen(&r), "AAGAAG");
    assert_eq!(r.score, Some(34));
    let c = check(&p, &dna("AAAAAG"), &o).unwrap();
    assert_eq!(c.satisfies, Some(false));
    assert!(c.violations.iter().any(|s| s.contains("zero-weight")));
    let o = DesignOptions {
        preferences: collagene_core::sequence::codon_design::CodonPreferences::Explicit {
            weights: vec![CodonWeight {
                codon: "ATG".into(),
                weight: 0,
            }],
        },
        ..Default::default()
    };
    let r = design(&protein("M"), &o).unwrap();
    assert_eq!(r.status, DesignStatus::Unsatisfiable);
    assert!(r.search_complete);
}
#[test]
fn motifs_check_both_strands_and_iupac_and_context_only_matches() {
    let p = protein("M");
    let o = DesignOptions {
        forbidden_motifs: vec!["CAT".into()],
        ..Default::default()
    };
    assert_eq!(design(&p, &o).unwrap().status, DesignStatus::Unsatisfiable);
    assert_eq!(check(&p, &dna("ATG"), &o).unwrap().satisfies, Some(false));
    let o = DesignOptions {
        forbidden_motifs: vec!["NAT".into()],
        left_context: "A".into(),
        ..Default::default()
    };
    assert_eq!(design(&p, &o).unwrap().status, DesignStatus::Unsatisfiable);
    let o = DesignOptions {
        forbidden_motifs: vec!["ACGT".into()],
        right_context: "ACGT".into(),
        ..Default::default()
    };
    assert_eq!(design(&p, &o).unwrap().status, DesignStatus::Unsatisfiable);
}
#[test]
fn homopolymers_cross_left_and_right_junctions_without_including_context_in_product() {
    let p = protein("K");
    let o = DesignOptions {
        left_context: "AA".into(),
        max_homopolymer: 4,
        ..Default::default()
    };
    let r = design(&p, &o).unwrap();
    assert_eq!(chosen(&r), "AAG");
    assert_eq!(r.candidate_check.as_ref().unwrap().coding_length, 3);
    let o = DesignOptions {
        right_context: "AAA".into(),
        max_homopolymer: 3,
        ..Default::default()
    };
    assert_eq!(chosen(&design(&p, &o).unwrap()), "AAG");
    let o = DesignOptions {
        left_context: "AAAA".into(),
        max_homopolymer: 4,
        ..Default::default()
    };
    assert_eq!(design(&p, &o).unwrap().status, DesignStatus::Unsatisfiable);
}
#[test]
fn window_gc_and_global_gc_are_exact_integer_counts_and_full_windows_only() {
    let p = protein("M");
    let o = DesignOptions {
        window_gc: Some(WindowGc {
            window: 2,
            bounds: GcBounds {
                minimum_bp: 1,
                maximum_bp: 10_000,
            },
        }),
        ..Default::default()
    };
    assert_eq!(design(&p, &o).unwrap().status, DesignStatus::Unsatisfiable);
    let o = DesignOptions {
        window_gc: Some(WindowGc {
            window: 20,
            bounds: GcBounds {
                minimum_bp: 10_000,
                maximum_bp: 10_000,
            },
        }),
        ..Default::default()
    };
    assert_eq!(check(&p, &dna("ATG"), &o).unwrap().satisfies, Some(true));
    let o = DesignOptions {
        gc: GcBounds {
            minimum_bp: 3334,
            maximum_bp: 10_000,
        },
        ..Default::default()
    };
    assert_eq!(design(&p, &o).unwrap().status, DesignStatus::Unsatisfiable);
    let o = DesignOptions {
        gc: GcBounds {
            minimum_bp: 3333,
            maximum_bp: 3334,
        },
        ..Default::default()
    };
    assert_eq!(chosen(&design(&p, &o).unwrap()), "ATG");
}
#[test]
fn beam_pruning_cannot_become_a_false_unsatisfiable_or_optimality_claim() {
    let p = protein("KK");
    let mut o = DesignOptions {
        beam_width: 1,
        max_homopolymer: 4,
        preferences: collagene_core::sequence::codon_design::CodonPreferences::Explicit {
            weights: vec![CodonWeight {
                codon: "AAA".into(),
                weight: 10,
            }],
        },
        ..Default::default()
    };
    let r = design(&p, &o).unwrap();
    assert_eq!(r.status, DesignStatus::SearchIncomplete);
    assert!(!r.search_complete && !r.optimality_proven);
    assert!(r.pruned_states > 0);
    o.beam_width = 16;
    let r = design(&p, &o).unwrap();
    assert_eq!(chosen(&r), "AAGAAA");
    assert!(r.search_complete && r.optimality_proven);
    let r = design(
        &protein("LLLL"),
        &DesignOptions {
            beam_width: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(r.status, DesignStatus::CandidateFound);
    assert!(!r.search_complete && !r.optimality_proven);
    assert!(r.candidate_check.unwrap().satisfies == Some(true));
}
#[test]
fn search_and_checker_work_exhaustion_are_explicit_and_never_release_unchecked_dna() {
    let p = protein("MK");
    let o = DesignOptions {
        max_work: 1,
        ..Default::default()
    };
    let r = design(&p, &o).unwrap();
    assert_eq!(r.status, DesignStatus::SearchIncomplete);
    assert!(r.candidate.is_none());
    assert!(r.work_used <= o.max_work);
    let r = check(&p, &dna("ATGAAA"), &o).unwrap();
    assert!(!r.complete);
    assert_eq!(r.satisfies, None);
    assert!(r.work_used <= o.max_work);
}
#[test]
fn wrong_translation_source_options_and_tampered_reports_are_rejected() {
    let p = protein("MK");
    let o = DesignOptions::default();
    let r = check(&p, &dna("ATGAAG"), &o).unwrap();
    r.validate_for(&p, &dna("ATGAAG"), &o).unwrap();
    assert_eq!(
        check(&p, &dna("ATGTAG"), &o).unwrap().satisfies,
        Some(false)
    );
    let mut forged = r.clone();
    forged.score += 1;
    assert!(forged.validate_for(&p, &dna("ATGAAG"), &o).is_err());
    assert!(r.validate_for(&p, &dna("ATGAAA"), &o).is_err());
    let mut r = design(&p, &o).unwrap();
    r.target_protein = "MW".into();
    assert!(r.validate_for_source(&p, &o).is_err());
}
#[test]
fn strict_json_and_option_bounds_do_not_silently_normalize_or_ignore_declarations() {
    let mut value = serde_json::to_value(DesignOptions::default()).unwrap();
    value["unknown_chemistry"] = true.into();
    assert!(serde_json::from_value::<DesignOptions>(value).is_err());
    let mut value = serde_json::to_value(DesignOptions::default()).unwrap();
    value["gc"]["extra"] = 1.into();
    assert!(serde_json::from_value::<DesignOptions>(value).is_err());
    for o in [
        DesignOptions {
            left_context: "atg".into(),
            ..Default::default()
        },
        DesignOptions {
            forbidden_motifs: vec!["A T".into()],
            ..Default::default()
        },
        DesignOptions {
            beam_width: 0,
            ..Default::default()
        },
        DesignOptions {
            max_work: 100_000_001,
            ..Default::default()
        },
        DesignOptions {
            genetic_code: 99,
            ..Default::default()
        },
        DesignOptions {
            preferences: collagene_core::sequence::codon_design::CodonPreferences::Explicit {
                weights: vec![
                    CodonWeight {
                        codon: "ATG".into(),
                        weight: 1,
                    },
                    CodonWeight {
                        codon: "ATG".into(),
                        weight: 2,
                    },
                ],
            },
            ..Default::default()
        },
    ] {
        assert!(o.validate().is_err());
    }
}
#[test]
fn unsupported_biology_and_retained_translation_exceptions_do_not_become_plain_coding_targets() {
    for p in ["*", "MX", "M*W", "MU", "MB", "MO", "M**"] {
        assert!(design(&protein(p), &DesignOptions::default()).is_err());
    }
    assert!(design(&dna("ATGNAA"), &DesignOptions::default()).is_err());
    assert!(design(&dna("ATGA"), &DesignOptions::default()).is_err());
    let mut circular = dna("ATG");
    circular.topology = Topology::Circular;
    assert!(design(&circular, &DesignOptions::default()).is_err());
    let mut m = protein("MK");
    m.metadata
        .insert("protein:modified".into(), "unknown".into());
    assert!(design(&m, &DesignOptions::default()).is_err());
    let mut m = dna("ATGAAA");
    m.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 6),
        qualifiers: vec![("transl_except".into(), Some("(pos:4..6,aa:Sec)".into()))],
        needs_review: false,
    });
    assert!(design(&m, &DesignOptions::default()).is_err());
    m.features[0].qualifiers = vec![("transl_table".into(), Some("2".into()))];
    assert!(design(&m, &DesignOptions::default()).is_err());
}
#[test]
fn sequence_design_does_not_copy_source_annotation_coordinates_or_edit_original() {
    let mut m = protein("MK");
    m.features.push(Feature {
        kind: "domain".into(),
        location: Location::range(0, 2),
        qualifiers: vec![("note".into(), Some("peptide annotation".into()))],
        needs_review: false,
    });
    let original = m.clone();
    let r = design(&m, &DesignOptions::default()).unwrap();
    assert_eq!(m, original);
    assert!(r.candidate.as_ref().unwrap().features.is_empty());
    assert!(r.warnings.iter().any(|w| w.contains("not copied")));
}
#[test]
fn oversized_metadata_evidence_and_deep_locations_reject_before_revision_hashing() {
    let mut m = protein("MK");
    m.original_genbank = Some("a".repeat(4 * 1024 * 1024));
    assert!(matches!(
        design(&m, &DesignOptions::default()),
        Err(collagene_model::Error::TooLarge)
    ));
    let mut m = protein("MK");
    let mut location = Location::range(0, 1);
    for _ in 0..34 {
        location = Location::Complement(Box::new(location));
    }
    m.features.push(Feature {
        kind: "domain".into(),
        location,
        qualifiers: vec![],
        needs_review: false,
    });
    assert!(matches!(
        design(&m, &DesignOptions::default()),
        Err(collagene_model::Error::TooLarge)
    ));
    let m = protein(&"M".repeat(MAX_RESIDUES + 1));
    assert!(matches!(
        design(&m, &DesignOptions::default()),
        Err(collagene_model::Error::TooLarge)
    ));
}

#[test]
fn violation_detail_truncation_and_binding_preflight_are_explicit() {
    let p = protein(&"K".repeat(500));
    let o = DesignOptions {
        preferences: collagene_core::sequence::codon_design::CodonPreferences::Explicit {
            weights: vec![CodonWeight {
                codon: "AAA".into(),
                weight: 0,
            }],
        },
        ..Default::default()
    };
    let c = check(&p, &dna(&"AAA".repeat(500)), &o).unwrap();
    assert!(c.complete);
    assert_eq!(c.satisfies, Some(false));
    assert_eq!(c.violations.len(), 256);
    assert!(c.warnings.iter().any(|s| s.contains("256")));
    let p = protein("MK");
    let o = DesignOptions::default();
    let r = design(&p, &o).unwrap();
    let mut oversized = p.clone();
    oversized.original_genbank = Some("x".repeat(4 * 1024 * 1024));
    assert!(matches!(
        r.validate_bindings_for_source(&oversized, &o),
        Err(collagene_model::Error::TooLarge)
    ));
    assert!(matches!(
        source_revision(&oversized, &o),
        Err(collagene_model::Error::TooLarge)
    ));
}

#[test]
fn expanded_tables_are_usable_but_dual_stop_codons_never_become_safe_design_choices() {
    let fixed = DesignOptions {
        genetic_code: 4,
        ..Default::default()
    };
    let r = design(&dna("TGA"), &fixed).unwrap();
    assert_eq!(r.target_protein, "W");
    assert_eq!(
        check(&protein("W"), &dna("TGA"), &fixed).unwrap().satisfies,
        Some(true)
    );
    for (id, residue, uncertain) in [(27, "W", "TGA"), (28, "Q", "TAA"), (31, "E", "TAG")] {
        let o = DesignOptions {
            genetic_code: id,
            preferences: collagene_core::sequence::codon_design::CodonPreferences::Explicit {
                weights: vec![CodonWeight {
                    codon: uncertain.into(),
                    weight: 100_000,
                }],
            },
            ..Default::default()
        };
        let r = design(&protein(residue), &o).unwrap();
        assert_ne!(chosen(&r), uncertain);
        assert_eq!(
            check(&protein(residue), &dna(uncertain), &o)
                .unwrap()
                .satisfies,
            Some(false)
        );
        assert!(design(&dna(uncertain), &o).is_err());
        r.validate_for_source(&protein(residue), &o).unwrap();
    }
    let no_definite_stop = DesignOptions {
        genetic_code: 28,
        ..Default::default()
    };
    let report = design(&protein("M*"), &no_definite_stop).unwrap();
    assert_eq!(report.status, DesignStatus::Unsatisfiable);
    assert!(report.candidate.is_none());
}
