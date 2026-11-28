//! Independent selected-record gap-placement and original-coordinate fixtures.
use collagene_core::sequence::crispr::{
    self,
    alignment::{AlignmentColumn, AlignmentLimitReached},
    *,
};
use collagene_model::{Molecule, Topology};
use std::collections::BTreeSet;
const QUERY: &str = "ACGTCAGTACGA";
fn profile(side: PamSide) -> NucleaseProfile {
    NucleaseProfile {
        id: "literal-gap-fixture".into(),
        version: "1".into(),
        guide_length: 12,
        pam: "GGA".into(),
        pam_side: side,
        nominal_top_offset: 5,
        nominal_bottom_offset: 8,
        custom: true,
        sources: vec![],
    }
}
fn options(side: PamSide) -> OffTargetOptions {
    OffTargetOptions {
        profile: profile(side),
        max_mismatches: 0,
        bulges: BulgePolicy {
            max_dna_bulge_bases: 1,
            max_rna_bulge_bases: 1,
            max_total_bulge_bases: 1,
            ..Default::default()
        },
        ..Default::default()
    }
}
fn reference(text: &str, circular: bool) -> ReferenceInput {
    ReferenceInput {
        id: "literal".into(),
        molecule: Molecule::dna(
            "original",
            text,
            if circular {
                Topology::Circular
            } else {
                Topology::Linear
            },
        )
        .unwrap(),
        build: None,
        accession: None,
    }
}
fn query() -> GuideQuery {
    GuideQuery {
        id: "query".into(),
        sequence: QUERY.into(),
    }
}
fn complement(b: u8) -> char {
    match b {
        b'A' => 'T',
        b'C' => 'G',
        b'G' => 'C',
        b'T' => 'A',
        _ => panic!(),
    }
}
fn rc(s: &str) -> String {
    s.bytes().rev().map(complement).collect()
}
fn operations(h: &OffTargetHit) -> BTreeSet<String> {
    h.alignment
        .alignments
        .iter()
        .map(|a| {
            a.columns
                .iter()
                .map(|c| match c {
                    AlignmentColumn::Paired { .. } => 'P',
                    AlignmentColumn::DnaBulge { .. } => 'D',
                    AlignmentColumn::RnaBulge { .. } => 'R',
                })
                .collect()
        })
        .collect()
}
type Key = (bool, Vec<usize>, Vec<usize>, BTreeSet<String>);
/// For total<=1, independently delete each single base, or compare the equal-width
/// words. This does not use the kernel's recursive search, masks or geometry helpers.
fn oracle(source: &ReferenceInput, side: PamSide) -> BTreeSet<Key> {
    let n = source.molecule.sequence.len();
    let circle = source.molecule.topology == Topology::Circular;
    let mut result = BTreeSet::new();
    for reverse in [false, true] {
        let oriented = if reverse {
            rc(source.molecule.sequence.as_str())
        } else {
            source.molecule.sequence.as_str().into()
        };
        for width in 11..=13 {
            let extent = width + 3;
            if extent > n {
                continue;
            }
            for start in 0..if circle { n } else { n - extent + 1 } {
                let (go, po) = if side == PamSide::ThreePrime {
                    (0, width)
                } else {
                    (3, 0)
                };
                let at = |j: usize| oriented.as_bytes()[(start + j) % n];
                if (0..3).map(|j| at(po + j)).collect::<Vec<_>>() != b"GGA" {
                    continue;
                }
                let target = (0..width).map(|j| at(go + j)).collect::<Vec<_>>();
                let mut paths = BTreeSet::new();
                match width {
                    12 if target == QUERY.as_bytes() => {
                        paths.insert("P".repeat(12));
                    }
                    13 => {
                        for extra in 0..13 {
                            let reduced = target
                                .iter()
                                .enumerate()
                                .filter_map(|(j, b)| (j != extra).then_some(*b))
                                .collect::<Vec<_>>();
                            if reduced == QUERY.as_bytes() {
                                paths.insert(format!(
                                    "{}D{}",
                                    "P".repeat(extra),
                                    "P".repeat(12 - extra)
                                ));
                            }
                        }
                    }
                    11 => {
                        for extra in 0..12 {
                            let reduced = QUERY
                                .bytes()
                                .enumerate()
                                .filter_map(|(j, b)| (j != extra).then_some(b))
                                .collect::<Vec<_>>();
                            if reduced == target {
                                paths.insert(format!(
                                    "{}R{}",
                                    "P".repeat(extra),
                                    "P".repeat(11 - extra)
                                ));
                            }
                        }
                    }
                    _ => {}
                }
                if paths.is_empty() {
                    continue;
                }
                let coord = |j: usize| {
                    if reverse {
                        n - 1 - (start + j) % n
                    } else {
                        (start + j) % n
                    }
                };
                result.insert((
                    reverse,
                    (0..width).map(|j| coord(go + j)).collect(),
                    (0..3).map(|j| coord(po + j)).collect(),
                    paths,
                ));
            }
        }
    }
    result
}
#[test]
fn selected_search_matches_independent_one_gap_oracle_both_pam_sides_strands_and_circular_origins()
{
    for side in [PamSide::ThreePrime, PamSide::FivePrime] {
        for target in ["ACGTTCAGTACGA", "ACGTAGTACGA", QUERY] {
            let text = if side == PamSide::ThreePrime {
                format!("CC{target}GGATT")
            } else {
                format!("CCGGA{target}TT")
            };
            for circular in [false, true] {
                for rotate in 0..if circular { text.len() } else { 1 } {
                    let rotated = format!("{}{}", &text[rotate..], &text[..rotate]);
                    let source = reference(&rotated, circular);
                    let expected = oracle(&source, side);
                    let report =
                        search_off_targets(&query(), &[source.clone()], &options(side)).unwrap();
                    assert!(report.complete, "{:?}", report.limits);
                    assert!(!report.index_used);
                    let actual = report
                        .hits
                        .iter()
                        .map(|h| {
                            (
                                h.reverse,
                                h.guide_positions.clone(),
                                h.pam_positions.clone(),
                                operations(h),
                            )
                        })
                        .collect::<BTreeSet<_>>();
                    assert_eq!(
                        actual, expected,
                        "{side:?} circle{circular} rotate{rotate} target{target}"
                    );
                    assert!(
                        report
                            .hits
                            .iter()
                            .all(|h| h.alignment.complete && h.definite)
                    );
                    replay_off_targets(&query(), &[source], &options(side), &report).unwrap();
                }
            }
        }
    }
}
#[test]
fn repeat_ties_caps_shared_columns_and_work_never_establish_absence() {
    let q = GuideQuery {
        id: "repeat".into(),
        sequence: "A".repeat(12),
    };
    let refs = vec![reference(&format!("{}GGA", "A".repeat(13)), false)];
    let mut o = options(PamSide::ThreePrime);
    let full = search_off_targets(&q, &refs, &o).unwrap();
    let h = full
        .hits
        .iter()
        .find(|h| !h.reverse && h.reference_guide.len() == 13)
        .unwrap();
    assert_eq!(h.alignment.alignments.len(), 13);
    assert!(h.alignment.complete);
    assert_eq!((h.nominal_top_cut, h.nominal_bottom_cut), (None, None));
    o.bulges.max_alignments_per_hit = 2;
    let capped = search_off_targets(&q, &refs, &o).unwrap();
    assert!(!capped.complete && !capped.exhaustive_for_selected_records);
    let h = capped
        .hits
        .iter()
        .find(|h| !h.reverse && h.reference_guide.len() == 13)
        .unwrap();
    assert_eq!(h.alignment.alignments.len(), 2);
    assert!(!h.alignment.complete);
    assert_eq!(
        h.alignment.limit_reached,
        Some(AlignmentLimitReached::Alignments)
    );
    assert_eq!((h.nominal_top_cut, h.nominal_bottom_cut), (None, None));
    let mut columns = options(PamSide::ThreePrime);
    columns.bulges.max_retained_alignment_columns = 1;
    let r = search_off_targets(&q, &refs, &columns).unwrap();
    assert!(!r.complete && r.hits.is_empty());
    assert!(r.limits.iter().any(|s| s.contains("column budget")));
    let mut work = options(PamSide::ThreePrime);
    work.max_comparisons = 1;
    let r = search_off_targets(&q, &refs, &work).unwrap();
    assert!(!r.complete && r.hits.is_empty());
    assert_eq!(r.comparisons, 1);
}
#[test]
fn intended_locus_requires_exact_ungapped_alignment_but_keeps_other_paths() {
    let q = GuideQuery {
        id: "repeat".into(),
        sequence: "A".repeat(12),
    };
    let source = reference(&format!("{}GGA", "A".repeat(12)), false);
    let mut o = options(PamSide::ThreePrime);
    o.bulges.max_total_bulge_bases = 2;
    o.intended = Some(IntendedTarget {
        reference_id: source.id.clone(),
        reference_revision: source.molecule.revision(),
        guide_positions: (0..12).collect(),
        pam_positions: vec![12, 13, 14],
        reverse: false,
    });
    let report = search_off_targets(&q, &[source.clone()], &o).unwrap();
    let h = report.hits.iter().find(|h| h.is_intended).unwrap();
    assert!(
        h.alignment
            .alignments
            .iter()
            .any(|a| a.dna_bulge_bases == 0 && a.rna_bulge_bases == 0)
    );
    assert!(
        h.alignment
            .alignments
            .iter()
            .any(|a| a.dna_bulge_bases == 1 && a.rna_bulge_bases == 1)
    );
    assert_eq!((h.nominal_top_cut, h.nominal_bottom_cut), (None, None));
    let mut wrong = query();
    wrong.sequence = "CAAAAAAAAAAA".into();
    o.max_mismatches = 1;
    let report = search_off_targets(&wrong, &[source], &o).unwrap();
    assert!(report.hits.iter().all(|h| !h.is_intended));
}
#[test]
fn declared_screening_pam_is_distinct_from_design_admission_and_is_bound() {
    let q = query();
    let refs = vec![reference(&format!("{QUERY}AGA"), false)];
    let mut o = options(PamSide::ThreePrime);
    o.bulges = BulgePolicy::default();
    assert!(search_off_targets(&q, &refs, &o).unwrap().hits.is_empty());
    o.pam_policy = PamPolicy::DeclaredScreening {
        pam: "RGA".into(),
        rationale: "Literal alternate screening scope".into(),
        sources: vec![],
    };
    let r = search_off_targets(&q, &refs, &o).unwrap();
    assert!(r.complete);
    let h = r.hits.iter().find(|h| !h.reverse).unwrap();
    assert_eq!(h.design_pam, PamEvidence::Excluded);
    assert!(h.pam_certain && h.definite && !h.is_intended);
    let mut changed = o.clone();
    changed.pam_policy = PamPolicy::DesignProfile;
    assert!(r.validate_bindings_for(&q, &refs, &changed).is_err());
    let mut altered = r.clone();
    altered.hits[0].design_pam = PamEvidence::Certain;
    assert!(altered.validate().is_err());
    let mut invalid = o;
    invalid.pam_policy = PamPolicy::DeclaredScreening {
        pam: "GG".into(),
        rationale: "invalid length".into(),
        sources: vec![],
    };
    assert!(search_off_targets(&q, &refs, &invalid).is_err());
}
#[test]
fn current_model_path_integrity_postcard_and_bounds_reject_fabricated_evidence() {
    let refs = vec![reference(&format!("ACGTTCAGTACGAGGA"), false)];
    let o = options(PamSide::ThreePrime);
    let r = search_off_targets(&query(), &refs, &o).unwrap();
    assert_eq!(r.alignment_model, crispr::alignment::MODEL_ID);
    assert_eq!(r.engine_version, "crispr-2");
    let encoded = postcard::to_allocvec(&r).unwrap();
    let round: OffTargetReport = postcard::from_bytes(&encoded).unwrap();
    assert_eq!(round, r);
    let mut old = r.clone();
    old.engine_version = "crispr-1".into();
    assert!(old.validate().is_err());
    let mut wrong = r.clone();
    wrong.hits[0].alignment.alignments[0].columns[0] = AlignmentColumn::RnaBulge { guide: 11 };
    assert!(wrong.validate().is_err());
    let mut wrong = r.clone();
    wrong.hits[0].alignment.alignments[0].dna_bulge_bases = 0;
    assert!(wrong.validate().is_err());
    let mut wrong = r.clone();
    wrong.hits[0].guide_positions[0] = usize::MAX;
    assert!(wrong.validate().is_err());
    let mut invalid = o.clone();
    invalid.bulges.max_total_bulge_bases = 5;
    assert!(invalid.validate().is_err());
    let mut json = serde_json::to_value(o).unwrap();
    json.as_object_mut().unwrap().remove("bulges");
    assert!(serde_json::from_value::<OffTargetOptions>(json).is_err());
}
#[test]
fn source_bound_added_profiles_keep_literal_lengths_pams_offsets_and_unsupported_names_out() {
    let profiles = catalogue();
    assert_eq!(profiles.len(), 5);
    for (id, length, pam, cut) in [
        ("sacas9-21", 21, "NNGRRT", 18),
        ("sacas9-22", 22, "NNGRRT", 19),
        ("spcas9-vrer", 20, "NGCG", 17),
    ] {
        let p = profiles.iter().find(|p| p.id == id).unwrap();
        assert_eq!(
            (
                p.guide_length,
                p.pam.as_str(),
                p.nominal_top_offset,
                p.nominal_bottom_offset
            ),
            (length, pam, cut, cut)
        );
        assert!(!p.custom && !p.sources.is_empty());
        let observed = if id.starts_with("sacas9") {
            "AAGAGT"
        } else {
            "AGCG"
        };
        let guide = "A".repeat(length);
        let text = format!("CC{guide}{observed}TT");
        let mut o = GuideOptions::default();
        o.profile = p.clone();
        let f = design_guides(&reference(&text, false).molecule, &o).unwrap();
        let h = f
            .candidates
            .iter()
            .find(|h| !h.reverse && h.guide_positions[0] == 2)
            .unwrap();
        assert_eq!(h.nominal_top_cut, Some(2 + cut as usize));
        let reversed = design_guides(&reference(&rc(&text), false).molecule, &o).unwrap();
        let h = reversed
            .candidates
            .iter()
            .find(|h| h.reverse && h.guide_sequence == guide)
            .unwrap();
        assert_eq!(h.nominal_top_cut, Some(text.len() - 2 - cut as usize));
        let mut altered = o.clone();
        altered.profile.pam = "N".repeat(pam.len());
        assert!(design_guides(&reference(&text, false).molecule, &altered).is_err());
    }
    assert!(!profiles.iter().any(|p| p.id.contains("vqr")));
}
