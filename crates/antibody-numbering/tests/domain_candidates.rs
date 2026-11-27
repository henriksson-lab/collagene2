use collagene_antibody_numbering::{
    Error, Insertion,
    domain_candidates::*,
    profile,
    references::human_references,
    v_mapping::{self, Locus, ProjectionIssue},
};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
// OGRDB IGHV3-23*01 V1..104 and IGHJ4*02 FR4. Joins/loops below are authored
// calculation fixtures, not naturally observed rearrangements. CC0 attribution:
// data/ogrdb-human-vj/ATTRIBUTION.md.
const V: &str = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC";
const J: &str = "WGQGTLVTVSS";
const LOOP: &str = "ARDYGMDV";
fn source(s: &str) -> Molecule {
    Molecule {
        name: "Authored V/J pairing fixture".into(),
        sequence: Sequence::new(Alphabet::Protein, s).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn options() -> Options {
    Options {
        v_scoring: v_mapping::Scoring {
            observed_reward: 5,
            frequency_bonus: 8,
            unseen_penalty: -6,
            gap_open: -12,
            gap_extend: -2,
            minimum_score: 300,
        },
        min_cdr3_residues: 2,
        max_cdr3_residues: 60,
        max_fr4_mismatches: 0,
    }
}
fn limits() -> Limits {
    Limits {
        v_mapping: v_mapping::Limits {
            max_cells: profile::MAX_CELLS,
            max_endpoints_per_locus: profile::MAX_ENDPOINTS,
        },
        max_scan_windows: MAX_SCAN_WINDOWS,
        max_candidates: MAX_CANDIDATES,
    }
}
fn sequence() -> String {
    format!("{V}{LOOP}{J}")
}
#[test]
fn heavy_pairing_preserves_all_identical_j_alleles_and_exact_conditional_loop_labels() {
    let m = source(&sequence());
    let r = find_candidates(&m, &options(), &limits(), &mut || false).unwrap();
    assert_eq!(r.qualified_pairs_in_retained_v, 1);
    assert_eq!(r.candidates.len(), 1);
    assert!(!r.v_retention_truncated && !r.pair_retention_truncated);
    let c = &r.candidates[0];
    assert_eq!(c.locus, Locus::Heavy);
    assert_eq!(
        (c.v_endpoint, c.cdr3.start, c.cdr3.end, c.source_end),
        (V.len(), V.len(), V.len() + LOOP.len(), sequence().len())
    );
    assert_eq!(
        c.cdr3
            .residues
            .iter()
            .map(|r| r.residue)
            .collect::<String>(),
        LOOP
    );
    assert_eq!(
        c.cdr3
            .residues
            .iter()
            .map(|r| r.position.number)
            .collect::<Vec<_>>(),
        [105, 106, 107, 108, 114, 115, 116, 117]
    );
    assert!(
        c.cdr3
            .residues
            .iter()
            .all(|r| r.position.insertion == Insertion::None)
    );
    assert_eq!(c.fr4.iter().map(|r| r.residue).collect::<String>(), J);
    assert_eq!(c.fr4.first().unwrap().position.number, 118);
    assert_eq!(c.fr4.last().unwrap().position.number, 128);
    assert_eq!(c.j_matches.len(), 1);
    assert_eq!(c.j_matches[0].allele_refs.len(), 3);
    for gene in ["IGHJ1*01", "IGHJ4*02", "IGHJ5*02"] {
        assert!(c.j_matches[0].allele_refs.iter().any(|s| s.contains(gene)));
    }
    assert!(c.j_matches[0].mismatches.is_empty());
    r.validate_for(&m, &options(), &limits()).unwrap();
}
#[test]
fn light_loci_and_tandem_constructs_keep_absolute_source_ownership() {
    let corpus = human_references().unwrap();
    for (gene, tail, locus) in [
        ("IGKV1-39*01", "FGQGTKVEIK", Locus::Kappa),
        ("IGLV2-14*01", "FGGGTKLTVL", Locus::Lambda),
    ] {
        let prefix = corpus
            .v_entries()
            .iter()
            .find(|v| v.label == gene)
            .unwrap()
            .aligned_protein
            .replace('-', "");
        let m = source(&format!("{prefix}QQSYSTP{tail}"));
        let r = find_candidates(&m, &options(), &limits(), &mut || false).unwrap();
        let c = r
            .candidates
            .iter()
            .find(|c| c.locus == locus)
            .expect("literal light reference join should remain a candidate");
        assert_eq!(c.cdr3.start, prefix.len());
        assert_eq!(c.cdr3.end, prefix.len() + 7);
        assert_eq!(c.fr4.last().unwrap().position.number, 127);
    }
    let linker = "GGGGSGGGGS";
    let block = sequence();
    let m = source(&format!("{block}{linker}{block}"));
    let r = find_candidates(&m, &options(), &limits(), &mut || false).unwrap();
    let expected = [V.len(), block.len() + linker.len() + V.len()];
    for start in expected {
        let c = r
            .candidates
            .iter()
            .find(|c| c.locus == Locus::Heavy && c.cdr3.start == start)
            .unwrap();
        assert_eq!(c.cdr3.end, start + LOOP.len());
        assert!(
            c.fr4
                .iter()
                .all(|r| r.source_index >= c.cdr3.end && r.source_index < c.source_end)
        );
    }
    assert!(r.candidates.windows(2).all(|w| w[0].id != w[1].id));
}
#[test]
fn competing_j_sites_and_substitution_matches_are_retained_without_silent_best_selection() {
    let m = source(&format!("{V}{LOOP}{J}AA{J}"));
    let r = find_candidates(&m, &options(), &limits(), &mut || false).unwrap();
    assert_eq!(r.qualified_pairs_in_retained_v, 2);
    assert_eq!(r.candidates.len(), 2);
    assert_eq!(
        r.candidates.iter().map(|c| c.cdr3.end).collect::<Vec<_>>(),
        [V.len() + LOOP.len(), V.len() + LOOP.len() + J.len() + 2]
    );
    let mut varied = options();
    varied.max_fr4_mismatches = 1;
    let r = find_candidates(&source(&sequence()), &varied, &limits(), &mut || false).unwrap();
    let c = &r.candidates[0];
    assert_eq!(c.j_matches.len(), 4);
    assert!(
        c.j_matches
            .iter()
            .any(|j| j.reference_protein == J && j.mismatches.is_empty())
    );
    let alternate = c
        .j_matches
        .iter()
        .find(|j| j.reference_protein == "WGQGTMVTVSS")
        .unwrap();
    assert_eq!(alternate.mismatches.len(), 1);
    assert_eq!(
        (
            alternate.mismatches[0].source_index,
            alternate.mismatches[0].observed,
            alternate.mismatches[0].reference
        ),
        (V.len() + LOOP.len() + 5, 'L', 'M')
    );
}
#[test]
fn counts_remain_scoped_to_retained_v_ends_and_pair_caps_never_claim_absence() {
    let m = source(&format!("{V}{LOOP}{J}AA{J}"));
    let mut capped = limits();
    capped.max_candidates = 1;
    let r = find_candidates(&m, &options(), &capped, &mut || false).unwrap();
    assert_eq!(
        (r.qualified_pairs_in_retained_v, r.candidates.len()),
        (2, 1)
    );
    assert!(r.pair_retention_truncated);
    let scanned = r.windows_examined;
    capped.max_candidates = 0;
    let r = find_candidates(&m, &options(), &capped, &mut || false).unwrap();
    assert_eq!(
        (r.qualified_pairs_in_retained_v, r.candidates.len()),
        (2, 0)
    );
    assert!(r.pair_retention_truncated);
    assert_eq!(r.windows_examined, scanned);
    capped.v_mapping.max_endpoints_per_locus = 0;
    let r = find_candidates(&m, &options(), &capped, &mut || false).unwrap();
    assert!(r.v_retention_truncated);
    assert_eq!(r.qualified_pairs_in_retained_v, 0);
    assert_eq!(r.windows_examined, 0);
}
#[test]
fn stops_partial_tails_conserved_motif_changes_and_loop_scope_do_not_create_pairs() {
    for s in [
        format!("{V}{LOOP}*{J}"),
        format!("{V}{LOOP}{}", &J[..5]),
        format!("{V}{LOOP}WAQGTLVTVSS"),
        format!("{V}A{J}"),
    ] {
        let mut o = options();
        o.max_fr4_mismatches = 3;
        let r = find_candidates(&source(&s), &o, &limits(), &mut || false).unwrap();
        assert!(r.candidates.is_empty(), "{s}");
        assert!(
            r.v_endpoints_without_j
                .iter()
                .any(|v| v.locus == Locus::Heavy)
        );
    }
    let mut o = options();
    o.max_cdr3_residues = 4;
    assert!(
        find_candidates(&source(&sequence()), &o, &limits(), &mut || false)
            .unwrap()
            .candidates
            .is_empty()
    );
    let r = find_candidates(
        &source(&format!("X*{}", sequence())),
        &options(),
        &limits(),
        &mut || false,
    )
    .unwrap();
    assert!(r.candidates.is_empty());
    assert!(!r.blocked_v_endpoints.is_empty());
    assert!(
        r.blocked_v_endpoints
            .iter()
            .all(|v| v.issue == ProjectionIssue::AmbiguousSource)
    );
}
#[test]
fn resource_bounds_cancellation_and_source_report_reconstruction_are_current() {
    let m = source(&sequence());
    let r = find_candidates(&m, &options(), &limits(), &mut || false).unwrap();
    let round: Report = serde_json::from_slice(&serde_json::to_vec(&r).unwrap()).unwrap();
    round.validate_for(&m, &options(), &limits()).unwrap();
    for mutation in 0..4 {
        let mut forged = r.clone();
        let c = &mut forged.candidates[0];
        match mutation {
            0 => c.cdr3.residues[0].position.number += 1,
            1 => c.fr4[0].source_index += 1,
            2 => c.j_matches[0].allele_refs[0] = "forged".into(),
            _ => c.id = "forged".into(),
        }
        assert_eq!(
            forged.validate_for(&m, &options(), &limits()),
            Err(Error::Binding)
        );
    }
    let mut renamed = m.clone();
    renamed.name = "Different document".into();
    assert_eq!(
        r.validate_for(&renamed, &options(), &limits()),
        Err(Error::Binding)
    );
    let mut o = options();
    o.max_fr4_mismatches = 4;
    assert!(find_candidates(&m, &o, &limits(), &mut || false).is_err());
    o = options();
    o.min_cdr3_residues = 1;
    assert!(find_candidates(&m, &o, &limits(), &mut || false).is_err());
    let mut l = limits();
    l.max_scan_windows = 1;
    assert_eq!(
        find_candidates(&m, &options(), &l, &mut || false).unwrap_err(),
        Error::Limit("V/J scan windows")
    );
    l = limits();
    l.max_candidates = MAX_CANDIDATES + 1;
    assert!(find_candidates(&m, &options(), &l, &mut || false).is_err());
    let mut calls = 0;
    assert_eq!(
        find_candidates(&m, &options(), &limits(), &mut || {
            calls += 1;
            calls == 10
        })
        .unwrap_err(),
        Error::Cancelled
    );
    assert!(calls >= 10);
    let dna = Molecule::dna("DNA", "ACGT", Topology::Linear).unwrap();
    assert!(find_candidates(&dna, &options(), &limits(), &mut || false).is_err());
}

#[test]
fn observed_antibody_with_tied_v_positions_keeps_conditional_cdr3_and_explicit_ambiguity() {
    // Independent archival 1BJ1 entity2 sequence; this checks conditional anchor
    // reconstruction, not complete IMGT assignment from depositor residue IDs.
    let panel: serde_json::Value =
        serde_json::from_str(include_str!("../data/rcsb-assessment/panel.json")).unwrap();
    let row = panel["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "1BJ1:2")
        .unwrap();
    let m = source(row["sequence"].as_str().unwrap());
    let r = find_candidates(&m, &options(), &limits(), &mut || false).unwrap();
    let c = r
        .candidates
        .iter()
        .find(|c| c.locus == Locus::Heavy && c.cdr3.start == 96)
        .unwrap();
    assert_eq!(
        c.v_projection_issue,
        Some(ProjectionIssue::TiedOptimalPaths)
    );
    assert_eq!(
        c.cdr3
            .residues
            .iter()
            .map(|r| r.residue)
            .collect::<String>(),
        "AKYPHYYGSSHWYFDV"
    );
    assert_eq!(c.fr4.iter().map(|r| r.residue).collect::<String>(), J);
    let v = r.v_evidence.mappings[0]
        .projections
        .iter()
        .find(|p| p.end == 96)
        .unwrap();
    assert!(v.projection.is_none());
    assert!(v.terminal_anchor.is_some());
    r.validate_for(&m, &options(), &limits()).unwrap();
}
