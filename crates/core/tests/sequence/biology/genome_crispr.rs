//! Independent literal/oriented-window oracle; no Molecule is built for large contigs.
use collagene_core::{
    sequence::crispr::{self, GuideOptions, GuideQuery},
    tracks::{crispr::*, genome::*},
};
use collagene_model::Topology;
const G: &str = "ACGTCAGTACGATCGATCGA";
fn build(
    records: &[(&str, String, Topology)],
    expected: bool,
) -> (GenomeManifest, Vec<EncodedBlock>) {
    let info = BuildInfo {
        expected_components: if expected {
            records
                .iter()
                .map(|(id, dna, t)| ExpectedComponent {
                    id: id.to_string(),
                    source_id: id.to_string(),
                    accession: None,
                    length: dna.len() as u64,
                    sequence_hash: Some(blake3::hash(dna.as_bytes()).to_hex().to_string()),
                    topology: *t,
                })
                .collect()
        } else {
            vec![]
        },
        ..Default::default()
    };
    let mut b = FastaBuilder::begin(info, IngestOptions::default()).unwrap();
    let mut out = vec![];
    for (id, dna, _) in records {
        out.extend(b.push(format!(">{id}\n").as_bytes()).unwrap());
        for s in dna.as_bytes().chunks(8191) {
            out.extend(b.push(s).unwrap())
        }
        out.extend(b.push(b"\n").unwrap())
    }
    let (last, f) = b.finish().unwrap();
    out.extend(last);
    (f.manifest, out)
}
fn get(req: &[BlockRequest], all: &[EncodedBlock]) -> Vec<EncodedBlock> {
    req.iter()
        .map(|r| {
            all.iter()
                .find(|b| b.contig_id == r.contig_id && b.descriptor.index == r.index)
                .unwrap()
                .clone()
        })
        .collect()
}
fn scan(m: &GenomeManifest, all: &[EncodedBlock], o: &GenomeSearchOptions) -> GenomeSearchSession {
    let mut s = GenomeSearchSession::begin(
        m,
        &GuideQuery {
            id: "guide".into(),
            sequence: G.into(),
        },
        o,
    )
    .unwrap();
    while !s.progress().finished {
        let b = get(&s.progress().next_blocks, all);
        s.advance(&b).unwrap();
    }
    s
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
        _ => panic!("fixture DNA"),
    }
}
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        _ => panic!("fixture complement"),
    }
}
fn rc(s: &str) -> String {
    s.bytes().rev().map(|b| complement(b) as char).collect()
}
type Tuple = (String, Vec<u64>, Vec<u64>, bool, usize, usize, bool);
fn oracle(id: &str, dna: &str, circle: bool, max: usize) -> Vec<Tuple> {
    let n = dna.len();
    let mut out = vec![];
    if n < 23 {
        return out;
    }
    for f in 0..if circle { n } else { n - 22 } {
        for reverse in [false, true] {
            let positions = (0..23)
                .map(|j| {
                    if reverse {
                        (f + 22 - j) % n
                    } else {
                        (f + j) % n
                    }
                })
                .collect::<Vec<_>>();
            let bases = positions
                .iter()
                .map(|p| {
                    if reverse {
                        complement(dna.as_bytes()[*p])
                    } else {
                        dna.as_bytes()[*p]
                    }
                })
                .collect::<Vec<_>>();
            if mask(bases[21]) & 4 == 0 || mask(bases[22]) & 4 == 0 {
                continue;
            }
            let mut lo = 0;
            let mut hi = 0;
            for (a, b) in bases[..20].iter().zip(G.bytes()) {
                if mask(*a) & mask(b) == 0 {
                    lo += 1
                }
                if mask(*a) != mask(b) {
                    hi += 1
                }
            }
            if lo <= max {
                out.push((
                    id.into(),
                    positions[..20].iter().map(|p| *p as u64).collect(),
                    positions[20..].iter().map(|p| *p as u64).collect(),
                    reverse,
                    lo,
                    hi,
                    mask(bases[21]) == 4 && mask(bases[22]) == 4,
                ));
            }
        }
    }
    out.sort();
    out
}
fn tuples(s: &GenomeSearchSession) -> Vec<Tuple> {
    let mut out = vec![];
    let mut at = 0;
    while at < s.progress().hits_retained {
        let page = s.page(at, 256).unwrap();
        at += page.len();
        out.extend(page.into_iter().map(|h| {
            (
                h.contig_id,
                h.guide_positions,
                h.pam_positions,
                h.reverse,
                h.min_mismatches,
                h.max_mismatches,
                h.pam_certain,
            )
        }))
    }
    out.sort();
    out
}
#[test]
fn streamed_both_strands_block_seams_and_small_final_halo_match_independent_oracle() {
    for tail in [1, 7, 23, 79] {
        let mut dna = vec![b'A'; BLOCK_BASES * 2 + tail];
        let f = BLOCK_BASES - 9;
        dna[f..f + 23].copy_from_slice(format!("{G}TGG").as_bytes());
        let r = BLOCK_BASES * 2 - 22;
        let reverse = rc(&format!("{G}AGG"));
        dna[r..r + 23].copy_from_slice(reverse.as_bytes());
        let text = String::from_utf8(dna).unwrap();
        let (m, blocks) = build(&[("chr", text.clone(), Topology::Linear)], true);
        let s = scan(
            &m,
            &blocks,
            &GenomeSearchOptions {
                max_mismatches: 0,
                ..Default::default()
            },
        );
        assert!(s.progress().complete);
        assert_eq!(tuples(&s), oracle("chr", &text, false, 0));
    }
}
#[test]
fn circular_origin_both_strands_and_global_cut_boundaries() {
    for reverse in [false, true] {
        let mut dna = vec![b'A'; BLOCK_BASES + 1];
        let f = dna.len() - 8;
        let text = if reverse {
            rc(&format!("{G}TGG"))
        } else {
            format!("{G}TGG")
        };
        let n = dna.len();
        for (i, b) in text.bytes().enumerate() {
            dna[(f + i) % n] = b
        }
        let text = String::from_utf8(dna).unwrap();
        let (m, blocks) = build(&[("chr", text.clone(), Topology::Circular)], true);
        let s = scan(
            &m,
            &blocks,
            &GenomeSearchOptions {
                max_mismatches: 0,
                ..Default::default()
            },
        );
        assert!(s.progress().complete);
        assert_eq!(tuples(&s), oracle("chr", &text, true, 0));
        let hit = s
            .page(0, 256)
            .unwrap()
            .into_iter()
            .find(|h| h.reverse == reverse && h.reference_guide == G)
            .unwrap();
        let expected = (f + if reverse { 6 } else { 17 }) % n;
        assert_eq!(hit.nominal_top_cut, Some(expected as u64));
    }
}
#[test]
fn ambiguous_owner_fallback_preserves_min_max_hamming_evidence() {
    let mut dna = vec![b'A'; BLOCK_BASES + 80];
    let mut motif = format!("{G}TGG").into_bytes();
    motif[7] = b'N';
    dna[BLOCK_BASES - 5..BLOCK_BASES + 18].copy_from_slice(&motif);
    let text = String::from_utf8(dna).unwrap();
    let (m, b) = build(&[("chr", text.clone(), Topology::Linear)], true);
    let s = scan(
        &m,
        &b,
        &GenomeSearchOptions {
            max_mismatches: 0,
            ..Default::default()
        },
    );
    assert!(s.progress().complete);
    assert_eq!(tuples(&s), oracle("chr", &text, false, 0));
    assert!(
        s.page(0, 256)
            .unwrap()
            .iter()
            .any(|h| h.min_mismatches == 0 && h.max_mismatches == 1 && !h.definite)
    );
}
#[test]
fn interval_footprint_ownership_stable_ids_and_original_positions() {
    let mut dna = vec![b'A'; BLOCK_BASES + 70];
    let f = BLOCK_BASES - 10;
    dna[f..f + 23].copy_from_slice(format!("{G}TGG").as_bytes());
    let text = String::from_utf8(dna).unwrap();
    let (m, b) = build(&[("chr", text, Topology::Linear)], true);
    let options = |start, end| GenomeGuideOptions {
        source_role: "source-role".into(),
        interval: GenomicInterval {
            contig_id: "chr".into(),
            start,
            end,
        },
        guide: GuideOptions::default(),
    };
    let a = options(f as u64, f as u64 + 1);
    let first =
        design_interval(&m, &a, &get(&interval_required_blocks(&m, &a).unwrap(), &b)).unwrap();
    assert!(first.complete);
    let h = first.candidates.iter().find(|h| !h.reverse).unwrap();
    assert_eq!(
        h.guide_positions,
        (f as u64..f as u64 + 20).collect::<Vec<_>>()
    );
    assert_eq!(h.nominal_top_cut, Some(f as u64 + 17));
    let a = options(f as u64 - 8, f as u64 + 5);
    let second =
        design_interval(&m, &a, &get(&interval_required_blocks(&m, &a).unwrap(), &b)).unwrap();
    assert!(second.candidates.iter().any(|c| c.id == h.id));
    let a = options(f as u64 + 1, f as u64 + 2);
    let excluded =
        design_interval(&m, &a, &get(&interval_required_blocks(&m, &a).unwrap(), &b)).unwrap();
    assert!(!excluded.candidates.iter().any(|c| c.id == h.id));
}
#[test]
fn complete_manifest_is_distinct_from_declared_expected_build_and_no_token_can_resume() {
    let (m, b) = build(&[("chr", format!("{G}TGG"), Topology::Linear)], false);
    let mut s = GenomeSearchSession::begin(
        &m,
        &GuideQuery {
            id: "q".into(),
            sequence: G.into(),
        },
        &GenomeSearchOptions::default(),
    )
    .unwrap();
    let progress = serde_json::to_value(s.progress()).unwrap();
    assert!(!progress["finished"].as_bool().unwrap());
    s.advance(&get(&s.progress().next_blocks, &b)).unwrap();
    assert!(s.progress().exhaustive_for_manifest);
    assert!(s.progress().source_integrity_verified);
    assert!(!s.progress().complete);
    assert!(!s.progress().expected_components_complete);
    assert!(s.advance(&b).is_err());
}
#[test]
fn tampered_missing_and_out_of_order_blocks_poison_private_search() {
    let (m, b) = build(
        &[("chr", "A".repeat(BLOCK_BASES + 100), Topology::Linear)],
        true,
    );
    for kind in 0..3 {
        let mut s = GenomeSearchSession::begin(
            &m,
            &GuideQuery {
                id: "q".into(),
                sequence: G.into(),
            },
            &GenomeSearchOptions::default(),
        )
        .unwrap();
        let mut blocks = get(&s.progress().next_blocks, &b);
        match kind {
            0 => blocks[0].encoded[0] ^= 16,
            1 => {
                blocks.pop();
            }
            _ => blocks.reverse(),
        }
        assert!(s.advance(&blocks).is_err());
        assert!(
            s.progress().finished
                && !s.progress().complete
                && !s.progress().source_integrity_verified
        );
        assert!(s.progress().next_blocks.is_empty());
        assert!(s.advance(&b).is_err());
    }
}
#[test]
fn canonical_hash_forgery_is_rejected_even_when_each_encoded_block_hash_matches() {
    let (mut m, b) = build(&[("chr", format!("{G}TGG"), Topology::Linear)], true);
    m.contigs[0].sequence_hash = "1".repeat(64);
    m.build.expected_components[0].sequence_hash = Some("1".repeat(64));
    let roles = vec![(
        &m.contigs[0].id,
        &m.contigs[0].accession,
        m.contigs[0].length,
        m.contigs[0].topology,
        &m.contigs[0].sequence_hash,
    )];
    m.content_revision =
        blake3::hash(&postcard::to_allocvec(&("cgref-canonical-content-v1", roles)).unwrap())
            .to_hex()
            .to_string();
    m.validate().unwrap();
    let mut s = GenomeSearchSession::begin(
        &m,
        &GuideQuery {
            id: "q".into(),
            sequence: G.into(),
        },
        &GenomeSearchOptions::default(),
    )
    .unwrap();
    assert!(s.advance(&b).is_err());
    assert!(!s.progress().complete && !s.progress().source_integrity_verified);
}
#[test]
fn intended_tags_bind_role_and_locus_and_never_drop_identical_dna() {
    let text = format!("{G}TGGAAA{G}AGG");
    let (m, b) = build(
        &[
            ("one", text.clone(), Topology::Linear),
            ("two", text, Topology::Linear),
        ],
        true,
    );
    let mut o = GenomeSearchOptions {
        reference_role: "dataset-A".into(),
        max_mismatches: 0,
        ..Default::default()
    };
    let s = scan(&m, &b, &o);
    let all = s.page(0, 256).unwrap();
    assert_eq!(all.len(), 4);
    let h = &all[0];
    o.intended = Some(GenomeIntended {
        reference_role: o.reference_role.clone(),
        manifest_hash: m.manifest_hash().unwrap(),
        content_revision: m.content_revision.clone(),
        contig_id: h.contig_id.clone(),
        guide_positions: h.guide_positions.clone(),
        pam_positions: h.pam_positions.clone(),
        reverse: h.reverse,
    });
    let s = scan(&m, &b, &o);
    assert_eq!(s.progress().hits_retained, 4);
    assert_eq!(
        s.page(0, 256)
            .unwrap()
            .iter()
            .filter(|h| h.is_intended)
            .count(),
        1
    );
    o.reference_role = "dataset-B".into();
    assert!(
        GenomeSearchSession::begin(
            &m,
            &GuideQuery {
                id: "q".into(),
                sequence: G.into()
            },
            &o
        )
        .is_err()
    );
}
#[test]
fn cancellation_and_resource_caps_never_claim_absence() {
    let (m, b) = build(
        &[("chr", format!("{G}TGGAAA{G}TGG"), Topology::Linear)],
        true,
    );
    for o in [
        GenomeSearchOptions {
            max_total_comparisons: 1,
            ..Default::default()
        },
        GenomeSearchOptions {
            max_hits: 1,
            ..Default::default()
        },
        GenomeSearchOptions {
            max_index_entries_per_block: 1,
            ..Default::default()
        },
    ] {
        let s = scan(&m, &b, &o);
        assert!(s.progress().finished && !s.progress().complete && !s.progress().limits.is_empty());
        assert!(!s.progress().exhaustive_for_manifest);
        assert!(!s.progress().source_integrity_verified);
        s.progress().validate().unwrap();
    }
    let mut s = GenomeSearchSession::begin(
        &m,
        &GuideQuery {
            id: "q".into(),
            sequence: G.into(),
        },
        &GenomeSearchOptions::default(),
    )
    .unwrap();
    assert!(s.cancel().cancelled);
    assert!(!s.progress().complete);
    assert!(!s.progress().exhaustive_for_manifest);
    assert!(!s.progress().source_integrity_verified);
    assert!(s.advance(&b).is_err());
}
#[test]
fn interval_read_circular_unwrapped_and_nominal_cuts_outside_local_window() {
    let text = format!("AAA{G}TGGCCC");
    let (m, b) = build(&[("chr", text.clone(), Topology::Circular)], true);
    let o = GenomeIntervalReadOptions {
        source_role: "persistent-role".into(),
        interval: GenomicInterval {
            contig_id: "chr".into(),
            start: text.len() as u64 - 3,
            end: text.len() as u64 + 5,
        },
    };
    let read = read_interval(
        &m,
        &o,
        &get(&read_interval_required_blocks(&m, &o).unwrap(), &b),
    )
    .unwrap();
    assert_eq!(
        read.sequence,
        format!("{}{}", &text[text.len() - 3..], &text[..5])
    );
    let mut o = GenomeGuideOptions {
        source_role: "persistent-role".into(),
        interval: GenomicInterval {
            contig_id: "chr".into(),
            start: 3,
            end: 4,
        },
        guide: GuideOptions::default(),
    };
    o.guide.profile.nominal_top_offset = 40;
    o.guide.profile.nominal_bottom_offset = 40;
    o.guide.profile.id = "custom-fixture".into();
    o.guide.profile.custom = true;
    let report =
        design_interval(&m, &o, &get(&interval_required_blocks(&m, &o).unwrap(), &b)).unwrap();
    let h = report
        .candidates
        .iter()
        .find(|h| !h.reverse && h.guide_sequence == G)
        .unwrap();
    assert_eq!(h.nominal_top_cut, Some(43 % text.len() as u64));
}
#[test]
fn stream_reference_larger_than_document_and_browser_file_caps() {
    let prefix = 33 * 1024 * 1024 + 17;
    let total = prefix + 23;
    let info = BuildInfo {
        expected_components: vec![ExpectedComponent {
            id: "large".into(),
            source_id: "large".into(),
            accession: None,
            length: total as u64,
            sequence_hash: None,
            topology: Topology::Linear,
        }],
        ..Default::default()
    };
    let mut builder = FastaBuilder::begin(info, IngestOptions::default()).unwrap();
    let mut all = builder.push(b">large\n").unwrap();
    let chunk = vec![b'A'; 1024 * 1024];
    let mut left = prefix;
    while left > 0 {
        let n = left.min(chunk.len());
        all.extend(builder.push(&chunk[..n]).unwrap());
        left -= n
    }
    all.extend(builder.push(format!("{G}TGG").as_bytes()).unwrap());
    let (tail, finished) = builder.finish().unwrap();
    all.extend(tail);
    assert!(finished.manifest.source_bytes > 32 * 1024 * 1024);
    assert!(finished.manifest.contigs[0].length > 10_000_000);
    assert!(all.iter().all(|b| b.encoded.len() <= BLOCK_BASES / 2));
    let o = GenomeGuideOptions {
        source_role: "large-role".into(),
        interval: GenomicInterval {
            contig_id: "large".into(),
            start: prefix as u64,
            end: prefix as u64 + 1,
        },
        guide: GuideOptions::default(),
    };
    let report = design_interval(
        &finished.manifest,
        &o,
        &get(
            &interval_required_blocks(&finished.manifest, &o).unwrap(),
            &all,
        ),
    )
    .unwrap();
    assert!(
        report
            .candidates
            .iter()
            .any(|h| !h.reverse && h.guide_positions[0] == prefix as u64 && h.guide_sequence == G)
    );
    let s = scan(
        &finished.manifest,
        &all,
        &GenomeSearchOptions {
            max_mismatches: 0,
            ..Default::default()
        },
    );
    assert!(s.progress().complete);
    assert_eq!(s.progress().bases_verified, total as u64);
    assert_eq!(s.progress().blocks_processed, all.len() as u64);
    assert_eq!(s.progress().hits_retained, 1);
    assert_eq!(s.page(0, 1).unwrap()[0].guide_positions[0], prefix as u64);
}
#[test]
fn short_circular_owner_reuses_one_verified_block_without_duplicate_loci() {
    let text = format!("{G}TGG");
    for shift in 0..23 {
        let rotated = format!("{}{}", &text[shift..], &text[..shift]);
        let (m, b) = build(&[("short", rotated.clone(), Topology::Circular)], true);
        let s = scan(
            &m,
            &b,
            &GenomeSearchOptions {
                max_mismatches: 0,
                ..Default::default()
            },
        );
        assert!(s.progress().complete);
        assert_eq!(tuples(&s), oracle("short", &rotated, true, 0));
        let options = GenomeGuideOptions {
            source_role: "short-role".into(),
            interval: GenomicInterval {
                contig_id: "short".into(),
                start: 0,
                end: 23,
            },
            guide: GuideOptions::default(),
        };
        let report = design_interval(
            &m,
            &options,
            &get(&interval_required_blocks(&m, &options).unwrap(), &b),
        )
        .unwrap();
        assert!(report.complete);
        assert_eq!(
            report
                .candidates
                .iter()
                .filter(|h| h.guide_sequence == G && !h.reverse)
                .count(),
            1
        );
        report.validate_bindings_for(&m, &options).unwrap();
    }
}
#[test]
fn cas12a_five_prime_pam_crosses_block_boundary_in_both_orientations() {
    let guide = format!("{G}TCA");
    let profile = crispr::cas12a();
    for reverse in [false, true] {
        let f = BLOCK_BASES - 8;
        let mut dna = vec![b'C'; BLOCK_BASES + 70];
        let motif = format!("TTTA{guide}");
        let motif = if reverse { rc(&motif) } else { motif };
        dna[f..f + 27].copy_from_slice(motif.as_bytes());
        let (m, b) = build(
            &[("cas12", String::from_utf8(dna).unwrap(), Topology::Linear)],
            true,
        );
        let o = GenomeGuideOptions {
            source_role: "source".into(),
            interval: GenomicInterval {
                contig_id: "cas12".into(),
                start: f as u64,
                end: f as u64 + 1,
            },
            guide: GuideOptions {
                profile: profile.clone(),
                ..Default::default()
            },
        };
        let report =
            design_interval(&m, &o, &get(&interval_required_blocks(&m, &o).unwrap(), &b)).unwrap();
        let h = report
            .candidates
            .iter()
            .find(|h| h.reverse == reverse && h.guide_sequence == guide)
            .unwrap();
        assert_eq!(h.pam_sequence, "TTTA");
        assert_eq!(
            h.nominal_top_cut,
            Some(if reverse { f as u64 } else { f as u64 + 22 })
        );
        assert_eq!(
            h.nominal_bottom_cut,
            Some(if reverse { f as u64 + 5 } else { f as u64 + 27 })
        );
        let query = GuideQuery {
            id: "cas12query".into(),
            sequence: guide.clone(),
        };
        let options = GenomeSearchOptions {
            profile: profile.clone(),
            max_mismatches: 0,
            ..Default::default()
        };
        let mut session = GenomeSearchSession::begin(&m, &query, &options).unwrap();
        while !session.progress().finished {
            session
                .advance(&get(&session.progress().next_blocks, &b))
                .unwrap();
        }
        assert!(session.progress().complete);
        let hits = session.page(0, 256).unwrap();
        let found = hits
            .iter()
            .find(|x| x.reverse == reverse && x.reference_guide == guide)
            .unwrap();
        assert_eq!(found.guide_positions, h.guide_positions);
        assert_eq!(found.pam_positions, h.pam_positions);
        assert_eq!(found.nominal_top_cut, h.nominal_top_cut);
        assert_eq!(found.nominal_bottom_cut, h.nominal_bottom_cut);
    }
}
#[test]
fn report_binding_checks_do_not_accept_another_role_or_manifest() {
    let (m, b) = build(&[("chr", format!("{G}TGG"), Topology::Linear)], true);
    let o = GenomeGuideOptions {
        source_role: "a".into(),
        interval: GenomicInterval {
            contig_id: "chr".into(),
            start: 0,
            end: 1,
        },
        guide: GuideOptions::default(),
    };
    let report =
        design_interval(&m, &o, &get(&interval_required_blocks(&m, &o).unwrap(), &b)).unwrap();
    report.validate_bindings_for(&m, &o).unwrap();
    let mut other = o.clone();
    other.source_role = "same-DNA-other-role".into();
    assert!(report.validate_bindings_for(&m, &other).is_err());
    let mut stale = m.clone();
    stale.build.name = "Another declared build".into();
    assert_eq!(stale.content_revision, m.content_revision);
    assert_ne!(stale.manifest_hash().unwrap(), m.manifest_hash().unwrap());
    assert!(report.validate_bindings_for(&stale, &o).is_err());
    let query = GuideQuery {
        id: "q".into(),
        sequence: G.into(),
    };
    let options = GenomeSearchOptions::default();
    let s = scan(&m, &b, &options);
    let mut report = s.progress();
    report.query.id = "q".into();
    report.validate_bindings_for(&m, &query, &options).unwrap();
    report.bases_verified -= 1;
    assert!(report.validate_bindings_for(&m, &query, &options).is_err());
    report.bases_verified += 1;
    report.finished = false;
    assert!(report.validate().is_err());
}
#[test]
fn an_explicit_missing_expected_chromosome_prevents_build_completeness() {
    let (mut m, b) = build(&[("present", format!("{G}TGG"), Topology::Linear)], true);
    m.build.expected_components.push(ExpectedComponent {
        id: "missing-role".into(),
        source_id: "missing-accession".into(),
        accession: Some("fixture-missing.1".into()),
        length: 100,
        sequence_hash: None,
        topology: Topology::Linear,
    });
    m.missing_expected.push("missing-role".into());
    m.expected_complete = false;
    m.validate().unwrap();
    let s = scan(&m, &b, &GenomeSearchOptions::default());
    let p = s.progress();
    assert!(p.exhaustive_for_manifest && p.source_integrity_verified);
    assert!(!p.complete && !p.expected_components_complete);
    assert_eq!(p.hits_retained, 1);
}
