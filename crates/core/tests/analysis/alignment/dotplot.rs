use collagene_core::alignment::dotplot::*;
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
fn dna(s: &str, circular: bool) -> Molecule {
    Molecule::dna(
        "independent dotplot input",
        s,
        if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )
    .unwrap()
}
fn options(w: usize) -> DotPlotOptions {
    DotPlotOptions {
        window: w,
        min_matches: w,
        max_work: 10_000_000,
        ..Default::default()
    }
}
fn key(r: &DotPlotReport) -> Vec<(usize, usize, bool, usize, usize)> {
    r.hits
        .iter()
        .map(|h| {
            (
                h.reference_start,
                h.query_start,
                h.reverse,
                h.compatible_matches,
                h.certain_matches,
            )
        })
        .collect()
}
// Small independent scalar oracle: direct biological coordinates, literal masks,
// no indexing, core comparison or span helpers.
fn nucleotide(b: u8) -> u8 {
    match b {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
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
        _ => panic!("not nucleotide"),
    }
}
fn oracle(
    a: &str,
    b: &str,
    o: &DotPlotOptions,
    ac: bool,
    bc: bool,
) -> Vec<(usize, usize, bool, usize, usize)> {
    let mut out = vec![];
    let nx = if ac { a.len() } else { a.len() - o.window + 1 };
    let ny = if bc { b.len() } else { b.len() - o.window + 1 };
    for reverse in [false, true]
        .into_iter()
        .take(if o.both_strands { 2 } else { 1 })
    {
        for x in (0..nx).step_by(o.stride) {
            for y in (0..ny).step_by(o.stride) {
                let mut compatible = 0;
                let mut certain = 0;
                for i in 0..o.window {
                    let u = nucleotide(a.as_bytes()[(x + i) % a.len()]);
                    let mut v = nucleotide(
                        b.as_bytes()[if reverse {
                            b.len() - 1 - (y + i) % b.len()
                        } else {
                            (y + i) % b.len()
                        }],
                    );
                    if reverse {
                        v = ((v & 1) << 3) | ((v & 2) << 1) | ((v & 4) >> 1) | ((v & 8) >> 3);
                    }
                    if u & v != 0 {
                        compatible += 1;
                        certain += usize::from(u == v && u.count_ones() == 1);
                    }
                }
                if compatible >= o.min_matches {
                    out.push((
                        x,
                        if reverse { b.len() - 1 - y } else { y },
                        reverse,
                        compatible,
                        certain,
                    ));
                }
            }
        }
    }
    out.sort();
    out
}
#[test]
fn exhaustive_small_concrete_and_ambiguous_windows_match_independent_full_cell_oracle() {
    for a in ["ACG", "GAT", "AAA", "RYN", "WSK"] {
        for b in ["ACG", "CAT", "AAA", "NNR"] {
            for circular in [false, true] {
                for w in 1..=3 {
                    for min_matches in 1..=w {
                        let mut o = options(w);
                        o.min_matches = min_matches;
                        let r = analyze(&dna(a, circular), &dna(b, circular), &o).unwrap();
                        r.validate().unwrap();
                        assert!(r.complete);
                        assert_eq!(r.examined_windows, r.total_windows);
                        assert_eq!(
                            key(&r),
                            oracle(a, b, &o, circular, circular),
                            "{a} {b} {w} {min_matches} {circular}"
                        );
                    }
                }
            }
        }
    }
}
#[test]
fn identical_palindromic_queries_preserve_both_orientation_records() {
    let o = options(2);
    let r = analyze(&dna("AT", false), &dna("AT", false), &o).unwrap();
    assert_eq!(r.hits.len(), 2);
    assert_eq!(r.hits[0].query_start, 0);
    assert!(!r.hits[0].reverse);
    assert_eq!(r.hits[1].query_start, 1);
    assert!(r.hits[1].reverse);
}
#[test]
fn reverse_circular_hits_retain_biological_span_order_without_extra_laps() {
    let o = options(4);
    let a = dna("ACGT", true);
    let r = analyze(&a, &a, &o).unwrap();
    r.validate_bindings_for(&a, &a, &o).unwrap();
    let h = r
        .hits
        .iter()
        .find(|h| h.reference_start == 0 && h.reverse && h.query_start == 3)
        .unwrap();
    assert_eq!(h.query_spans, vec![Span { start: 0, end: 4 }]);
    for h in &r.hits {
        assert_eq!(
            h.reference_spans
                .iter()
                .map(|s| s.end - s.start)
                .sum::<usize>(),
            4
        );
        assert_eq!(
            h.query_spans.iter().map(|s| s.end - s.start).sum::<usize>(),
            4
        );
    }
    let mut long = o.clone();
    long.window = 5;
    long.min_matches = 5;
    assert!(analyze(&a, &a, &long).is_err());
}
#[test]
fn stride_reports_grid_completion_separately_from_full_domain_completion() {
    let mut o = options(2);
    o.stride = 2;
    let a = dna("ACGTAC", true);
    let r = analyze(&a, &a, &o).unwrap();
    assert!(r.grid_complete && r.sampled && !r.complete);
    assert_eq!(r.total_windows, 18);
    assert_eq!(key(&r), oracle("ACGTAC", "ACGTAC", &o, true, true));
    r.validate().unwrap();
}
#[test]
fn work_hit_and_index_caps_never_claim_absence_or_exhaustiveness() {
    let a = dna("AAAAAAAA", false);
    for kind in 0..3 {
        let mut o = options(2);
        match kind {
            0 => o.max_work = 1,
            1 => o.max_hits = 1,
            _ => o.max_index_entries = 1,
        };
        let r = analyze(&a, &a, &o).unwrap();
        assert!(!r.complete && !r.grid_complete);
        assert!(!r.limits.is_empty());
        assert!(r.work_used <= o.max_work);
        r.validate().unwrap();
    }
    let mut o = options(2);
    o.min_matches = 1;
    o.max_work = 5;
    let r = analyze(&dna("NNNNNNNN", false), &a, &o).unwrap();
    assert!(!r.complete);
    r.validate().unwrap();
}
#[test]
fn protein_masks_keep_unknown_residues_possible_not_certain_and_stop_distinct() {
    fn protein(s: &str) -> Molecule {
        Molecule {
            name: "protein".into(),
            sequence: Sequence::new(Alphabet::Protein, s).unwrap(),
            topology: Topology::Linear,
            features: vec![],
            metadata: Default::default(),
            original_genbank: None,
        }
    }
    let mut o = options(1);
    o.both_strands = false;
    let a = protein("BXUO*");
    let b = protein("DNXUO*");
    let r = analyze(&a, &b, &o).unwrap();
    r.validate().unwrap();
    assert!(
        r.hits
            .iter()
            .any(|h| h.reference_start == 0 && h.query_start == 0 && !h.certain)
    );
    assert!(
        r.hits
            .iter()
            .any(|h| h.reference_start == 1 && h.query_start == 3 && !h.certain)
    );
    assert!(
        r.hits
            .iter()
            .any(|h| h.reference_start == 1 && h.query_start == 4 && !h.certain)
    );
    assert!(
        !r.hits
            .iter()
            .any(|h| h.reference_start == 1 && h.query_start == 5)
    );
    assert!(
        r.hits
            .iter()
            .any(|h| h.reference_start == 4 && h.query_start == 5 && h.certain)
    );
    o.both_strands = true;
    assert!(analyze(&a, &b, &o).is_err());
}
#[test]
fn rna_reverse_complement_and_mixed_alphabet_diagnostics_are_explicit() {
    let mut a = dna("ACG", false);
    a.sequence = Sequence::new(Alphabet::Rna, "ACG").unwrap();
    let mut b = a.clone();
    b.sequence = Sequence::new(Alphabet::Rna, "CGU").unwrap();
    let r = analyze(&a, &b, &options(3)).unwrap();
    assert!(r.hits.iter().any(|h| h.reverse && h.certain));
    assert!(analyze(&a, &dna("ACG", false), &options(3)).is_err());
}
#[test]
fn source_options_and_coordinate_report_tampering_is_rejected() {
    let a = dna("ACGT", true);
    let o = options(2);
    let r = analyze(&a, &a, &o).unwrap();
    let mut forged = r.clone();
    forged.hits[0].query_start = 99;
    assert!(forged.validate().is_err());
    let mut forged = r.clone();
    forged.examined_windows += 1;
    assert!(forged.validate().is_err());
    let mut changed = o.clone();
    changed.stride = 2;
    assert!(r.validate_bindings_for(&a, &a, &changed).is_err());
    let mut renamed = a.clone();
    renamed.name = "different exact document revision".into();
    assert!(r.validate_bindings_for(&renamed, &a, &o).is_err());
    let mut ignored = serde_json::to_value(&o).unwrap();
    ignored["invented_identity_mode"] = serde_json::json!(true);
    assert!(serde_json::from_value::<DotPlotOptions>(ignored).is_err());
}
#[test]
fn evidence_size_and_location_depth_fail_before_hashing() {
    let mut a = dna("ACGT", false);
    a.original_genbank = Some("x".repeat(16 * 1024 * 1024));
    assert!(matches!(
        source_revisions(&a, &a, &options(2)),
        Err(collagene_model::Error::TooLarge)
    ));
    let mut a = dna("ACGT", false);
    let mut l = Location::range(0, 1);
    for _ in 0..34 {
        l = Location::Complement(Box::new(l));
    }
    a.features.push(Feature {
        kind: "region".into(),
        location: l,
        qualifiers: vec![],
        needs_review: false,
    });
    assert!(matches!(
        analyze(&a, &a, &options(2)),
        Err(collagene_model::Error::TooLarge)
    ));
}
#[test]
fn svg_maps_reverse_and_origin_crossings_as_bounded_vector_geometry() {
    let a = dna("ACGT", true);
    let r = analyze(&a, &a, &options(3)).unwrap();
    let svg = export_svg(&r).unwrap();
    assert!(svg.starts_with("<svg "));
    assert!(svg.contains("reference "));
    assert!(svg.contains("#ad426c"));
    assert!(svg.contains("#267b61"));
    assert!(svg.contains("complete true"));
    assert!(!svg.contains("NaN"));
}
#[test]
fn indexed_plot_is_practical_for_large_sampled_documents_without_discarding_repeats() {
    let sequence = "ACGTGATCAGTCGACATGCA".repeat(1000);
    let a = dna(&sequence, false);
    let mut o = options(15);
    o.stride = 100;
    o.max_hits = 100_000;
    let r = analyze(&a, &a, &o).unwrap();
    assert!(r.grid_complete);
    assert!(!r.complete && r.sampled);
    assert!(r.hits.len() >= 1000);
    r.validate().unwrap();
}
#[test]
fn dense_output_budget_is_explicit_and_the_retained_report_remains_valid() {
    let a = dna(&"A".repeat(500), false);
    let o = options(1);
    let r = analyze(&a, &a, &o).unwrap();
    assert!(!r.grid_complete && !r.complete);
    assert!(r.hits.len() < o.max_hits);
    assert!(r.limits.iter().any(|v| v.contains("Output JSON")));
    assert!(r.examined_windows >= r.hits.len() as u64);
    r.validate().unwrap();
    assert!(serde_json::to_vec(&r).unwrap().len() <= 16 * 1024 * 1024);
}
#[test]
fn one_megabase_sampled_plot_has_full_declared_grid_with_bounded_work() {
    let mut state = 17u32;
    let mut text = String::with_capacity(1_000_000);
    for _ in 0..1_000_000 {
        state = state.wrapping_mul(1664525).wrapping_add(1013904223);
        text.push(b"ACGT"[(state >> 30) as usize] as char);
    }
    let a = dna(&text, false);
    let o = DotPlotOptions {
        stride: 1000,
        ..Default::default()
    };
    let start = std::time::Instant::now();
    let r = analyze(&a, &a, &o).unwrap();
    eprintln!(
        "1Mb sampled plot: {:?}, work {}, hits {}, grid {}",
        start.elapsed(),
        r.work_used,
        r.hits.len(),
        r.total_windows
    );
    assert!(r.grid_complete && r.sampled && !r.complete);
    assert_eq!(r.total_windows, 2_000_000);
    assert!(r.work_used <= o.max_work);
    for x in (0..1_000_000 - 14).step_by(1000) {
        assert!(
            r.hits
                .iter()
                .any(|h| h.reference_start == x && h.query_start == x && !h.reverse)
        );
    }
    r.validate_bindings_for(&a, &a, &o).unwrap();
}
