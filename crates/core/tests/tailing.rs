use collagene_core::{
    Sequence, Topology,
    chemistry::{self, Duplex, FivePrime, Nick, Polarity, Strand, Terminus, ThreePrime},
    tailing::*,
};
fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => panic!("literal DNA oracle"),
        })
        .collect()
}
fn blunt(sequence: &str, chemistry: Terminus) -> Duplex {
    Duplex {
        top: Sequence::dna(sequence).unwrap(),
        bottom: Sequence::dna(&rc(sequence)).unwrap(),
        bottom_start: 0,
        topology: Topology::Linear,
        left: chemistry,
        right: chemistry,
        nicks: vec![],
    }
}
fn phosphorylated(sequence: &str) -> Duplex {
    blunt(sequence, Terminus::CLEAVED)
}
fn check_coverage(report: &TailingReport) {
    let mut top = vec![0; report.duplex.top.len()];
    let mut bottom = vec![0; report.duplex.bottom.len()];
    for mapping in &report.mappings {
        assert!(!mapping.synthesized && !mapping.reverse);
        assert_eq!(mapping.source_strand, mapping.product_strand);
        let from = match mapping.source_strand {
            Strand::Top => report.source.top.as_str(),
            Strand::Bottom => report.source.bottom.as_str(),
        };
        let to = match mapping.product_strand {
            Strand::Top => report.duplex.top.as_str(),
            Strand::Bottom => report.duplex.bottom.as_str(),
        };
        assert_eq!(
            &from[mapping.source_start..mapping.source_end],
            &to[mapping.product_start..mapping.product_end]
        );
        let coverage = match mapping.product_strand {
            Strand::Top => &mut top,
            Strand::Bottom => &mut bottom,
        };
        for item in &mut coverage[mapping.product_start..mapping.product_end] {
            *item += 1;
        }
    }
    for addition in &report.additions {
        let to = match addition.strand {
            Strand::Top => report.duplex.top.as_str(),
            Strand::Bottom => report.duplex.bottom.as_str(),
        };
        assert_eq!(
            addition.sequence,
            &to[addition.product_start..addition.product_end]
        );
        assert_eq!(addition.product_end - addition.product_start, 1);
        let coverage = match addition.strand {
            Strand::Top => &mut top,
            Strand::Bottom => &mut bottom,
        };
        for item in &mut coverage[addition.product_start..addition.product_end] {
            *item += 1;
        }
    }
    assert!(top.iter().chain(&bottom).all(|n| *n == 1));
}
#[test]
fn literal_a_tailing_extends_actual_three_prime_suffixes_without_template_ancestry() {
    let source = phosphorylated("ACGC");
    let report = tail_blunt_duplex(&source, TailingBase::A).unwrap();
    assert_eq!(report.duplex.top.as_str(), "ACGCA");
    assert_eq!(report.duplex.bottom.as_str(), "GCGTA");
    assert_eq!(report.duplex.bottom_start, -1);
    assert_eq!(report.source, source);
    assert_eq!(report.duplex.left, source.left);
    assert_eq!(report.duplex.right, source.right);
    let (left, right) = report.duplex.ends().unwrap();
    assert_eq!(left.polarity, Polarity::ThreePrime);
    assert_eq!(right.polarity, Polarity::ThreePrime);
    assert_eq!(
        (left.overhang.as_str(), right.overhang.as_str()),
        ("A", "A")
    );
    assert_eq!(report.mappings.len(), 2);
    assert_eq!(report.additions.len(), 2);
    check_coverage(&report);
    report.validate().unwrap();
}
#[test]
fn all_1360_short_source_and_tail_pairs_follow_independent_coordinates() {
    let alphabet = b"ACGT";
    for length in 1..=4 {
        for code in 0..4usize.pow(length as u32) {
            let mut code = code;
            let text = (0..length)
                .map(|_| {
                    let b = alphabet[code % 4];
                    code /= 4;
                    b as char
                })
                .collect::<String>();
            let source = phosphorylated(&text);
            for base in [
                TailingBase::A,
                TailingBase::C,
                TailingBase::G,
                TailingBase::T,
            ] {
                let report = tail_blunt_duplex(&source, base).unwrap();
                let symbol = base.as_str();
                assert_eq!(report.duplex.top.as_str(), format!("{text}{symbol}"));
                assert_eq!(
                    report.duplex.bottom.as_str(),
                    format!("{}{symbol}", rc(&text))
                );
                for i in 0..length {
                    let bottom_index = length - i - 1;
                    assert_eq!(
                        report.duplex.top.as_str().as_bytes()[i],
                        rc(&report.duplex.bottom.as_str()[bottom_index..bottom_index + 1])
                            .as_bytes()[0]
                    );
                }
                let reverse = tail_blunt_duplex(&source.reverse().unwrap(), base).unwrap();
                assert_eq!(report.duplex.reverse().unwrap(), reverse.duplex);
                check_coverage(&report);
            }
        }
    }
}
#[test]
fn nicked_sources_keep_exact_boundaries_and_five_prime_states() {
    let mut source = blunt("ACGC", Terminus::PCR);
    source.left.five_prime = FivePrime::Unknown;
    source.nicks = vec![
        Nick {
            strand: Strand::Top,
            position: 1,
        },
        Nick {
            strand: Strand::Bottom,
            position: 3,
        },
    ];
    let report = tail_blunt_duplex(&source, TailingBase::G).unwrap();
    assert_eq!(report.duplex.nicks, source.nicks);
    assert_eq!(report.duplex.left, source.left);
    assert_eq!(report.duplex.right, source.right);
    assert_eq!(
        report.duplex.reverse().unwrap(),
        tail_blunt_duplex(&source.reverse().unwrap(), TailingBase::G)
            .unwrap()
            .duplex
    );
    assert!(!report.duplex.summary().unwrap().covalently_closed_in_model);
}
#[test]
fn single_base_pairs_are_valid_and_nucleotide_complements_are_explicit() {
    for base in [
        TailingBase::A,
        TailingBase::C,
        TailingBase::G,
        TailingBase::T,
    ] {
        let report = tail_blunt_duplex(&phosphorylated("A"), base).unwrap();
        assert_eq!(report.duplex.top.len(), 2);
        assert_eq!(report.duplex.bottom.len(), 2);
        assert_eq!(base.complement().complement(), base);
        assert_eq!(rc(base.as_str()), base.complement().as_str());
        report.validate().unwrap();
    }
}
#[test]
fn both_ta_junctions_add_two_reference_bases_and_allow_either_insert_orientation() {
    let vector = tail_blunt_duplex(&phosphorylated("ACGC"), TailingBase::T)
        .unwrap()
        .duplex;
    let insert = tail_blunt_duplex(&phosphorylated("GGTA"), TailingBase::A)
        .unwrap()
        .duplex;
    let product = chemistry::ligate(&[vector.clone(), insert.clone()], true, false).unwrap();
    assert_eq!(product.top.as_str(), "ACGCTGGTAA");
    assert_eq!(product.bottom.as_str(), "TTACCAGCGT");
    assert_eq!(product.top.len(), 4 + 4 + 2);
    assert!(product.nicks.is_empty());
    assert_eq!(product.topology, Topology::Circular);
    let reversed = chemistry::ligate(&[vector, insert.reverse().unwrap()], true, false).unwrap();
    assert_eq!(reversed.top.as_str(), "ACGCTTACCA");
    assert_ne!(product.top, reversed.top);
}
#[test]
fn gc_junctions_and_dephosphorylated_vector_keep_physical_nicks() {
    let vector = tail_blunt_duplex(&blunt("ACGC", Terminus::PCR), TailingBase::C)
        .unwrap()
        .duplex;
    let insert = tail_blunt_duplex(&phosphorylated("GGTA"), TailingBase::G)
        .unwrap()
        .duplex;
    assert!(chemistry::ligate(&[vector.clone(), insert.clone()], true, false).is_err());
    let nicked = chemistry::ligate(&[vector.clone(), insert.clone()], true, true).unwrap();
    assert_eq!(nicked.top.as_str(), "ACGCCGGTAG");
    assert_eq!(
        nicked.nicks,
        vec![
            Nick {
                strand: Strand::Bottom,
                position: 4
            },
            Nick {
                strand: Strand::Top,
                position: 0
            }
        ]
    );
    assert!(!nicked.summary().unwrap().covalently_closed_in_model);
    let sealed = chemistry::ligate(
        &[
            vector
                .treat(chemistry::EndTreatment::Phosphorylate)
                .unwrap(),
            insert,
        ],
        true,
        false,
    )
    .unwrap();
    assert_eq!(sealed.top, nicked.top);
    assert!(sealed.nicks.is_empty());
}
#[test]
fn no_bonds_without_explicit_phosphates_and_noncomplementary_tails_reject() {
    let a = tail_blunt_duplex(&blunt("ACGC", Terminus::PCR), TailingBase::A)
        .unwrap()
        .duplex;
    let t = tail_blunt_duplex(&blunt("GGTA", Terminus::PCR), TailingBase::T)
        .unwrap()
        .duplex;
    assert!(chemistry::ligate(&[a.clone(), t], true, true).is_err());
    let a_p = tail_blunt_duplex(&phosphorylated("GGTA"), TailingBase::A)
        .unwrap()
        .duplex;
    assert!(chemistry::ligate(&[a, a_p], true, true).is_err());
    let mut unknown = phosphorylated("ACGC");
    unknown.left.five_prime = FivePrime::Unknown;
    let u = tail_blunt_duplex(&unknown, TailingBase::T).unwrap().duplex;
    let a = tail_blunt_duplex(&phosphorylated("GGTA"), TailingBase::A)
        .unwrap()
        .duplex;
    assert!(chemistry::ligate(&[u, a], true, true).is_err());
}
#[test]
fn staggered_circular_unpaired_and_inaccessible_ends_never_silently_repair() {
    let mut circular = phosphorylated("ACGC");
    circular.topology = Topology::Circular;
    assert!(tail_blunt_duplex(&circular, TailingBase::A).is_err());
    for which in [0, 1] {
        for state in [ThreePrime::Blocked, ThreePrime::Unknown] {
            let mut d = phosphorylated("ACGC");
            if which == 0 {
                d.left.three_prime = state
            } else {
                d.right.three_prime = state
            }
            assert!(tail_blunt_duplex(&d, TailingBase::A).is_err());
        }
    }
    let staggered = Duplex {
        top: Sequence::dna("AATGCC").unwrap(),
        bottom: Sequence::dna("ACGGCA").unwrap(),
        bottom_start: 2,
        topology: Topology::Linear,
        left: Terminus::CLEAVED,
        right: Terminus::CLEAVED,
        nicks: vec![],
    };
    staggered.validate().unwrap();
    assert!(tail_blunt_duplex(&staggered, TailingBase::A).is_err());
    let mut discordant = phosphorylated("ACGC");
    discordant.bottom = Sequence::dna("AAAA").unwrap();
    assert!(tail_blunt_duplex(&discordant, TailingBase::A).is_err());
    let mut unpaired = phosphorylated("ACGC");
    unpaired.bottom_start = 10;
    assert!(tail_blunt_duplex(&unpaired, TailingBase::A).is_err());
    let mut invalid_nick = phosphorylated("ACGC");
    invalid_nick.nicks.push(Nick {
        strand: Strand::Top,
        position: 0,
    });
    assert!(tail_blunt_duplex(&invalid_nick, TailingBase::A).is_err());
}
#[test]
fn retained_bases_must_be_concrete_dna() {
    let ambiguous = Duplex {
        top: Sequence::dna("AN").unwrap(),
        bottom: Sequence::dna("NT").unwrap(),
        bottom_start: 0,
        topology: Topology::Linear,
        left: Terminus::CLEAVED,
        right: Terminus::CLEAVED,
        nicks: vec![],
    };
    assert!(tail_blunt_duplex(&ambiguous, TailingBase::A).is_err());
    let mut rna = phosphorylated("ACGC");
    rna.top = Sequence::new(collagene_core::Alphabet::Rna, "ACGC").unwrap();
    assert!(tail_blunt_duplex(&rna, TailingBase::A).is_err());
}
#[test]
fn source_binding_and_full_replay_reject_tampered_derived_claims() {
    let report = tail_blunt_duplex(&phosphorylated("ACGC"), TailingBase::A).unwrap();
    let json = serde_json::to_vec(&report).unwrap();
    let loaded: TailingReport = serde_json::from_slice(&json).unwrap();
    assert_eq!(loaded, report);
    loaded.validate().unwrap();
    for field in 0..8 {
        let mut bad = report.clone();
        match field {
            0 => bad.base = TailingBase::T,
            1 => bad.source_duplex_blake3 = "0".repeat(64),
            2 => bad.engine_version = "unknown".into(),
            3 => bad.duplex.top = Sequence::dna("ACGCT").unwrap(),
            4 => bad.mappings[0].synthesized = true,
            5 => bad.additions[0].product_start = 0,
            6 => bad.warnings.clear(),
            _ => bad.source.left.five_prime = FivePrime::Hydroxyl,
        }
        assert!(bad.validate().is_err(), "field{field}");
    }
    let changed = tail_blunt_duplex(&blunt("ACGC", Terminus::PCR), TailingBase::A).unwrap();
    assert_ne!(changed.source_duplex_blake3, report.source_duplex_blake3);
    let mut value = serde_json::to_value(&report).unwrap();
    value["additions"][0]["source_start"] = serde_json::json!(0);
    assert!(serde_json::from_value::<TailingReport>(value).is_err());
}
#[test]
fn oversized_report_is_rejected_before_source_product_cloning() {
    assert_eq!(MAX_TAILING_REPORT_BYTES, 32 * 1024 * 1024);
    let source = blunt(&"A".repeat(8_400_000), Terminus::CLEAVED);
    assert!(matches!(
        tail_blunt_duplex(&source, TailingBase::A),
        Err(collagene_core::Error::TooLarge)
    ));
}
#[test]
fn strict_json_loader_does_not_discard_unknown_physical_or_ancestry_evidence() {
    let report = tail_blunt_duplex(&phosphorylated("ACGC"), TailingBase::A).unwrap();
    let bytes = serde_json::to_vec(&report).unwrap();
    assert_eq!(load_tailing_report(&bytes).unwrap(), report);
    for location in 0..3 {
        let mut raw = serde_json::to_value(&report).unwrap();
        match location {
            0 => raw["source"]["protein_adduct"] = serde_json::json!(true),
            1 => raw["duplex"]["ligase_verified"] = serde_json::json!(true),
            _ => raw["mappings"][0]["unrecorded_source"] = serde_json::json!("other"),
        };
        assert!(load_tailing_report(&serde_json::to_vec(&raw).unwrap()).is_err());
    }
    let mut raw = serde_json::to_value(&report).unwrap();
    raw["source"]["top"]["bases"] = serde_json::json!("acgc");
    assert!(load_tailing_report(&serde_json::to_vec(&raw).unwrap()).is_err());
}
