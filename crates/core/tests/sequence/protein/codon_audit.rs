//! Independent finite-space oracle: literal NCBI assignments, exhaustive products,
//! integer inequalities and direct string scans; no optimizer helpers are used.
use collagene_core::sequence::codon_design::*;
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
fn protein(s: &str) -> Molecule {
    Molecule {
        name: "independent target peptide".into(),
        sequence: Sequence::new(Alphabet::Protein, s).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn codons(a: u8) -> &'static [&'static str] {
    match a {
        b'A' => &["GCA", "GCC", "GCG", "GCT"],
        b'C' => &["TGC", "TGT"],
        b'D' => &["GAC", "GAT"],
        b'E' => &["GAA", "GAG"],
        b'F' => &["TTC", "TTT"],
        b'G' => &["GGA", "GGC", "GGG", "GGT"],
        b'H' => &["CAC", "CAT"],
        b'I' => &["ATA", "ATC", "ATT"],
        b'K' => &["AAA", "AAG"],
        b'L' => &["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"],
        b'M' => &["ATG"],
        b'N' => &["AAC", "AAT"],
        b'P' => &["CCA", "CCC", "CCG", "CCT"],
        b'Q' => &["CAA", "CAG"],
        b'R' => &["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"],
        b'S' => &["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"],
        b'T' => &["ACA", "ACC", "ACG", "ACT"],
        b'V' => &["GTA", "GTC", "GTG", "GTT"],
        b'W' => &["TGG"],
        b'Y' => &["TAC", "TAT"],
        b'*' => &["TAA", "TAG", "TGA"],
        _ => panic!("oracle unknown amino acid"),
    }
}
fn enumerate(peptide: &str) -> Vec<String> {
    let mut rows = vec![String::new()];
    for aa in peptide.bytes() {
        rows = rows
            .into_iter()
            .flat_map(|prefix| codons(aa).iter().map(move |c| format!("{prefix}{c}")))
            .collect();
    }
    rows
}
fn mask(b: u8) -> u8 {
    match b {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}
fn reverse(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => panic!("not concrete DNA"),
        })
        .collect()
}
fn gc(s: &str, bounds: &GcBounds) -> bool {
    let count = s.bytes().filter(|b| matches!(b, b'G' | b'C')).count();
    count * 10000 >= s.len() * usize::from(bounds.minimum_bp)
        && count * 10000 <= s.len() * usize::from(bounds.maximum_bp)
}
fn allowed(s: &str, o: &DesignOptions) -> bool {
    if !gc(s, &o.gc) {
        return false;
    }
    if s.as_bytes().chunks_exact(3).any(|c| {
        declared_weights(o)
            .iter()
            .any(|w| w.codon.as_bytes() == c && w.weight == 0)
    }) {
        return false;
    }
    let whole = format!("{}{s}{}", o.left_context, o.right_context);
    let mut previous = 0;
    let mut run = 0;
    for b in whole.bytes() {
        if previous == b {
            run += 1;
        } else {
            previous = b;
            run = 1;
        }
        if run > o.max_homopolymer {
            return false;
        }
    }
    for orientation in [&whole, &reverse(&whole)] {
        for motif in &o.forbidden_motifs {
            if orientation.as_bytes().windows(motif.len()).any(|window| {
                window
                    .iter()
                    .zip(motif.bytes())
                    .all(|(b, m)| mask(*b) & mask(m) != 0)
            }) {
                return false;
            }
        }
    }
    if let Some(window) = &o.window_gc {
        for w in whole.as_bytes().windows(window.window) {
            if !gc(std::str::from_utf8(w).unwrap(), &window.bounds) {
                return false;
            }
        }
    }
    true
}
fn score(s: &str, o: &DesignOptions) -> u64 {
    s.as_bytes()
        .chunks_exact(3)
        .map(|c| {
            u64::from(
                declared_weights(o)
                    .iter()
                    .find(|w| w.codon.as_bytes() == c)
                    .map_or(1, |w| w.weight),
            )
        })
        .sum()
}
fn oracle(peptide: &str, o: &DesignOptions) -> Option<String> {
    enumerate(peptide)
        .into_iter()
        .filter(|s| allowed(s, o))
        .max_by(|a, b| score(a, o).cmp(&score(b, o)).then_with(|| b.cmp(a)))
}
fn compare(peptide: &str, o: &DesignOptions) {
    let source = protein(peptide);
    let expected = oracle(peptide, o);
    let report = design(&source, o).unwrap();
    report.validate_for_source(&source, o).unwrap();
    assert!(
        report.search_complete,
        "exhaustive finite case was truncated: {peptide} {report:?}"
    );
    match expected {
        Some(expected) => {
            assert_eq!(report.status, DesignStatus::CandidateFound);
            assert!(report.optimality_proven);
            let candidate = report.candidate.as_ref().unwrap();
            assert_eq!(candidate.sequence.as_str(), expected, "{peptide}");
            assert_eq!(report.score, Some(score(&expected, o)));
            let c = check(&source, candidate, o).unwrap();
            assert_eq!(c.satisfies, Some(true));
            assert!(c.complete);
        }
        None => {
            assert_eq!(report.status, DesignStatus::Unsatisfiable);
            assert!(report.candidate.is_none());
        }
    }
}
#[test]
fn every_canonical_pair_matches_an_independent_exhaustive_codon_oracle() {
    let alphabet = "ACDEFGHIKLMNPQRSTVWY";
    let mut o = DesignOptions::default();
    o.beam_width = 256;
    for a in alphabet.chars() {
        for b in alphabet.chars() {
            compare(&format!("{a}{b}"), &o);
        }
    }
}
#[test]
fn constrained_weighted_stop_terminated_design_matches_all_216_literal_options() {
    let mut o = DesignOptions::default();
    o.beam_width = 256;
    o.max_homopolymer = 4;
    o.preferences = collagene_core::sequence::codon_design::CodonPreferences::Explicit {
        weights: vec![
            CodonWeight {
                codon: "AAA".into(),
                weight: 40,
            },
            CodonWeight {
                codon: "CTG".into(),
                weight: 17,
            },
            CodonWeight {
                codon: "TCG".into(),
                weight: 11,
            },
            CodonWeight {
                codon: "TGA".into(),
                weight: 5,
            },
        ],
    };
    o.left_context = "GTC".into();
    o.right_context = "AAC".into();
    for (lo, hi) in [(0, 10000), (2500, 5000), (5000, 7500), (8333, 10000)] {
        for motif in ["ACGA", "AAR", "CTGA"] {
            o.gc = GcBounds {
                minimum_bp: lo,
                maximum_bp: hi,
            };
            o.forbidden_motifs = vec![motif.into()];
            compare("KLS*", &o);
        }
    }
}
#[test]
fn fixed_junction_contexts_reverse_motifs_and_integer_gc_boundaries_are_not_lost() {
    let mut o = DesignOptions::default();
    o.beam_width = 256;
    o.max_homopolymer = 3;
    o.right_context = "AAA".into();
    compare("K", &o);
    assert_eq!(oracle("K", &o).as_deref(), Some("AAG"));
    o.right_context.clear();
    o.left_context = "G".into();
    o.forbidden_motifs = vec!["GTT".into()];
    compare("F", &o);
    assert!(oracle("F", &o).is_none());
    o.left_context.clear();
    o.forbidden_motifs = vec!["ACGA".into()];
    let bad = Molecule::dna("literal reverse-site candidate", "TCGTTT", Topology::Linear).unwrap();
    let c = check(&protein("SF"), &bad, &o).unwrap();
    assert!(c.complete);
    assert_eq!(c.satisfies, Some(false));
    o.forbidden_motifs.clear();
    o.max_homopolymer = 12;
    o.gc = GcBounds {
        minimum_bp: 3333,
        maximum_bp: 3333,
    };
    compare("K", &o);
    assert!(oracle("K", &o).is_none());
    o.gc = GcBounds {
        minimum_bp: 3333,
        maximum_bp: 3334,
    };
    compare("K", &o);
    assert_eq!(oracle("K", &o).as_deref(), Some("AAG"));
}
#[test]
fn every_context_window_uses_the_full_construct_and_no_partial_short_window_is_invented() {
    let mut o = DesignOptions::default();
    o.beam_width = 256;
    o.left_context = "CG".into();
    o.right_context = "GC".into();
    for window in [2, 3, 5, 8] {
        o.window_gc = Some(WindowGc {
            window,
            bounds: GcBounds {
                minimum_bp: 2000,
                maximum_bp: 8000,
            },
        });
        compare("KF", &o);
    }
}
#[test]
fn independent_candidate_check_does_not_accept_target_or_options_substitution() {
    let o = DesignOptions::default();
    let source = protein("K");
    let candidate = Molecule::dna("untrusted proposed DNA", "ATG", Topology::Linear).unwrap();
    let r = check(&source, &candidate, &o).unwrap();
    assert!(r.complete);
    assert_eq!(r.satisfies, Some(false));
    let correct = Molecule::dna("lysine", "AAA", Topology::Linear).unwrap();
    let good = check(&source, &correct, &o).unwrap();
    good.validate_for(&source, &correct, &o).unwrap();
    assert!(good.validate_for(&protein("M"), &correct, &o).is_err());
    assert!(good.validate_for(&source, &candidate, &o).is_err());
    let mut changed = o.clone();
    changed.forbidden_motifs = vec!["AAA".into()];
    assert!(good.validate_for(&source, &correct, &changed).is_err());
}
#[test]
fn alternate_code_assignments_and_bacterial_internal_decoding_have_literal_targets() {
    let mut o = DesignOptions::default();
    o.genetic_code = 2;
    o.preferences = collagene_core::sequence::codon_design::CodonPreferences::Explicit {
        weights: vec![
            CodonWeight {
                codon: "ATA".into(),
                weight: 10,
            },
            CodonWeight {
                codon: "TGA".into(),
                weight: 10,
            },
        ],
    };
    let r = design(&protein("MW"), &o).unwrap();
    assert_eq!(r.candidate.unwrap().sequence.as_str(), "ATATGA");
    let dna = Molecule::dna("mitochondrial coding sequence", "ATATGA", Topology::Linear).unwrap();
    let r = design(&dna, &o).unwrap();
    assert_eq!(r.target_protein, "MW");
    o.genetic_code = 11;
    o.preferences = collagene_core::sequence::codon_design::CodonPreferences::Uniform;
    let dna = Molecule::dna(
        "ordinary codons, no inferred initiation",
        "GTGAAA",
        Topology::Linear,
    )
    .unwrap();
    let r = design(&dna, &o).unwrap();
    assert_eq!(r.target_protein, "VK");
    assert_eq!(r.candidate.unwrap().sequence.as_str(), "GTAAAA");
}
#[test]
fn default_budget_designs_a_full_kiloresidue_gene_with_bound_source_and_check() {
    for peptide in ["L".repeat(1000), "ACDEFGHIKLMNPQRSTVWY".repeat(50)] {
        let source = protein(&peptide);
        let options = DesignOptions::default();
        let report = design(&source, &options).unwrap();
        assert_eq!(
            report.status,
            DesignStatus::CandidateFound,
            "full-gene design failed: {:?}",
            report.limits
        );
        let candidate = report.candidate.as_ref().unwrap();
        assert_eq!(candidate.sequence.len(), 3000);
        assert_eq!(report.target_protein, peptide);
        assert!(!report.optimality_proven);
        assert!(report.work_used <= options.max_work);
        let checked = check(&source, candidate, &options).unwrap();
        assert_eq!(checked.satisfies, Some(true));
        assert!(checked.complete);
    }
}

fn declared_weights(options: &DesignOptions) -> &[CodonWeight] {
    use collagene_core::sequence::codon_design::CodonPreferences;
    match &options.preferences {
        CodonPreferences::Uniform => &[],
        CodonPreferences::Explicit { weights } => weights,
        CodonPreferences::ReferenceHost { .. } => {
            panic!("this independent audit fixture declares manual weights")
        }
    }
}
