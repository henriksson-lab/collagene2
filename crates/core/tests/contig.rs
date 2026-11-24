//! Synthetic read fixtures with independently specified sequence/call evidence.
use collagene_core::{
    contig::*,
    oligo_formats::{self, QualityEncoding},
    traces, *,
};
fn options() -> ContigOptions {
    ContigOptions {
        min_overlap: 12,
        trim: ContigTrimOptions {
            min_retained_bases: 1,
            ..Default::default()
        },
        ..Default::default()
    }
}
fn genome(n: usize) -> String {
    let mut x = 173u32;
    (0..n)
        .map(|_| {
            x = x.wrapping_mul(1664525).wrapping_add(1013904223);
            b"ACGT"[((x >> 24) & 3) as usize] as char
        })
        .collect()
}
fn plain(id: &str, s: &str) -> ReadInput {
    ReadInput {
        id: id.into(),
        molecule: Molecule::dna(id, s, Topology::Linear).unwrap(),
    }
}
fn fastq(id: &str, s: &str, scores: &[u8]) -> ReadInput {
    let q = scores.iter().map(|x| (x + 33) as char).collect::<String>();
    let text = format!("@{id}\n{s}\n+\n{q}\n");
    ReadInput {
        id: id.into(),
        molecule: oligo_formats::import_fastq(text.as_bytes(), QualityEncoding::Phred33)
            .unwrap()
            .molecules
            .remove(0),
    }
}
fn best(r: &ContigReport) -> &Contig {
    r.contigs.iter().max_by_key(|c| c.placements.len()).unwrap()
}
fn matches_gauge(s: &str, expected: &str) -> bool {
    s == expected
        || Sequence::dna(s)
            .unwrap()
            .reverse_complement()
            .unwrap()
            .as_str()
            == expected
}
#[test]
fn three_reads_with_transitive_overlap_form_a_source_supported_contig() {
    let g = genome(150);
    let reads = vec![
        plain("a", &g[..90]),
        plain("b", &g[40..120]),
        plain("c", &g[70..]),
    ];
    let r = assemble(&reads, &options()).unwrap();
    r.validate_bindings_for(&reads, &options()).unwrap();
    assert!(r.complete, "{:?}", r.limits);
    let c = best(&r);
    assert_eq!(c.placements.len(), 3, "{r:#?}");
    assert!(matches_gauge(&c.sequence, &g));
    assert!(c.complete && c.layout_unique);
    assert!(r.overlaps.iter().any(|e| e.accepted && !e.selected));
    assert_eq!(
        consensus_documents(&reads, &options(), &r).unwrap().len(),
        r.contigs.len()
    );
}
#[test]
fn contained_reverse_read_keeps_every_current_and_original_call() {
    let g = genome(150);
    let rc = Sequence::dna(&g[50..90])
        .unwrap()
        .reverse_complement()
        .unwrap();
    let reads = vec![
        plain("a", &g[..90]),
        plain("b", &g[40..120]),
        plain("c", &g[70..]),
        plain("contained", rc.as_str()),
    ];
    let r = assemble(&reads, &options()).unwrap();
    let c = best(&r);
    assert_eq!(c.placements.len(), 4, "{r:#?}");
    assert!(matches_gauge(&c.sequence, &g));
    let p = c.placements.iter().find(|p| p.read_index == 3).unwrap();
    let root = c.placements.iter().find(|p| p.read_index == 0).unwrap();
    assert_ne!(p.reverse, root.reverse);
    for call in &p.calls {
        assert_eq!(call.call_index, call.stored_call_index);
        assert!(
            c.columns[call.column]
                .observations
                .iter()
                .any(|o| o.read_index == 3 && o.call_index == call.call_index)
        );
    }
}
#[test]
fn low_quality_substitution_does_not_overrule_an_observed_eligible_call() {
    let g = genome(60);
    let mut changed = g.as_bytes().to_vec();
    changed[25] = if changed[25] == b'A' { b'C' } else { b'A' };
    let mut scores = vec![40; 60];
    scores[25] = 0;
    let reads = vec![
        fastq("bad", std::str::from_utf8(&changed).unwrap(), &scores),
        fastq("good", &g, &vec![40; 60]),
    ];
    let r = assemble(&reads, &options()).unwrap();
    let c = best(&r);
    assert_eq!(c.placements.len(), 2, "{r:#?}");
    assert!(matches_gauge(&c.sequence, &g));
    let at = c
        .columns
        .iter()
        .find(|c| {
            c.observations
                .iter()
                .any(|o| o.read_index == 0 && o.call_index == 25)
        })
        .unwrap();
    assert_eq!(at.state, ConsensusState::SupportedBase);
    assert!(
        !at.observations
            .iter()
            .find(|o| o.read_index == 0)
            .unwrap()
            .eligible
    );
}
#[test]
fn conflicting_eligible_substitution_stays_n_instead_of_majority_or_reference_fill() {
    let g = genome(60);
    let mut changed = g.as_bytes().to_vec();
    changed[25] = if changed[25] == b'A' { b'C' } else { b'A' };
    let reads = vec![
        plain("a", &g),
        plain("b", std::str::from_utf8(&changed).unwrap()),
        plain("c", &g),
    ];
    let r = assemble(&reads, &options()).unwrap();
    let c = best(&r);
    assert!(
        c.columns
            .iter()
            .any(|c| c.state == ConsensusState::ConflictingCalls)
    );
    assert!(c.sequence.contains('N'));
    assert!(!c.complete);
}
#[test]
fn low_quality_insertion_is_removed_only_by_eligible_flanking_gap_evidence() {
    let g = genome(60);
    let position = 30;
    let inserted = *b"ACGT"
        .iter()
        .find(|b| **b != g.as_bytes()[position - 1] && **b != g.as_bytes()[position])
        .unwrap();
    let mut bad = g.clone();
    bad.insert(position, inserted as char);
    let mut scores = vec![40; 61];
    scores[position] = 0;
    let reads = vec![
        fastq("good", &g, &vec![40; 60]),
        fastq("bad-insert", &bad, &scores),
    ];
    let r = assemble(&reads, &options()).unwrap();
    let c = best(&r);
    assert_eq!(c.placements.len(), 2, "{r:#?}");
    assert!(matches_gauge(&c.sequence, &g), "{r:#?}");
    let gap = c
        .columns
        .iter()
        .find(|c| c.state == ConsensusState::SupportedGap)
        .unwrap();
    assert_eq!(gap.base, None);
    assert_eq!(gap.consensus_position, None);
    assert!(!gap.eligible_gap_reads.is_empty());
    assert!(gap.observations.iter().all(|o| !o.eligible));
}
#[test]
fn good_insertion_versus_good_gap_remains_explicitly_conflicting() {
    let g = genome(60);
    let position = 30;
    let inserted = *b"ACGT"
        .iter()
        .find(|b| **b != g.as_bytes()[position - 1] && **b != g.as_bytes()[position])
        .unwrap();
    let mut bad = g.clone();
    bad.insert(position, inserted as char);
    let reads = vec![plain("a", &g), plain("b", &bad)];
    let r = assemble(&reads, &options()).unwrap();
    let c = best(&r);
    assert_eq!(c.placements.len(), 2, "{r:#?}");
    assert!(
        c.columns
            .iter()
            .any(|c| c.state == ConsensusState::GapConflict)
    );
    assert!(c.sequence.contains('N'));
    assert!(!c.complete);
}
#[test]
fn periodic_duplicates_do_not_prove_copy_number_or_multiply_support() {
    let text = "ACGT".repeat(12);
    let reads = vec![plain("copy-a", &text), plain("copy-b", &text)];
    let mut o = options();
    o.min_support = 2;
    let r = assemble(&reads, &o).unwrap();
    assert!(r.reads.iter().all(|r| r.duplicate_group.is_some()));
    for c in &r.contigs {
        assert!(!c.complete && !c.layout_unique);
        for column in &c.columns {
            assert!(column.counts.iter().all(|n| *n <= 1));
        }
    }
}
#[test]
fn read_and_pair_budgets_never_turn_partial_graphs_into_unique_layouts() {
    let g = genome(100);
    let reads = vec![plain("a", &g[..60]), plain("b", &g[40..])];
    let mut o = options();
    o.max_total_cells = 1;
    let r = assemble(&reads, &o).unwrap();
    assert!(!r.complete && !r.limits.is_empty());
    assert!(r.contigs.iter().all(|c| !c.complete && !c.layout_unique));
    o.max_total_cells = 128_000_000;
    o.max_edges = 1;
    let r = assemble(&reads, &o).unwrap();
    assert!(!r.complete && !r.limits.is_empty());
    assert!(r.contigs.iter().all(|c| !c.layout_unique));
}
#[test]
fn source_options_or_generated_metadata_tampering_cannot_be_replayed() {
    let g = genome(100);
    let reads = vec![plain("a", &g[..60]), plain("b", &g[40..])];
    let o = options();
    let r = assemble(&reads, &o).unwrap();
    let mut bad = r.clone();
    bad.contigs[0]
        .molecule
        .as_mut()
        .unwrap()
        .metadata
        .insert("contig:options".into(), "{}".into());
    assert!(bad.validate().is_err());
    let mut altered = reads.clone();
    altered[0].molecule.name.push_str(" changed");
    assert!(r.validate_bindings_for(&altered, &o).is_err());
    let mut bad = r.clone();
    bad.contigs[0].sequence.replace_range(..1, "N");
    assert!(consensus_documents(&reads, &o, &bad).is_err());
    // A structurally plausible imported score is not source authority. The
    // lightweight binding check deliberately avoids replay; fresh science must
    // reject even a mutation that leaves the sequence and source hashes intact.
    let mut forged = r.clone();
    forged.contigs[0].columns[0].observations[0].raw_score = Some(40);
    forged.validate_bindings_for(&reads, &o).unwrap();
    assert!(consensus_documents(&reads, &o, &forged).is_err());
    let mut old = r.clone();
    old.engine_version = "historical-contig-engine".into();
    for c in &mut old.contigs {
        if let Some(m) = &mut c.molecule {
            m.metadata
                .insert("contig:engine".into(), old.engine_version.clone());
        }
    }
    old.validate().unwrap();
    assert!(old.validate_bindings_for(&reads, &o).is_err());
}
#[test]
fn missing_quality_requires_declared_permission_and_stale_fastq_never_falls_back() {
    let g = genome(30);
    let reads = vec![plain("plain", &g)];
    let mut o = options();
    o.trim.allow_unscored = false;
    let r = assemble(&reads, &o).unwrap();
    assert!(!r.complete && r.reads[0].error.is_some());
    let mut read = fastq("stale", &g, &vec![40; 30]);
    read.molecule.sequence = Sequence::dna(&"A".repeat(30)).unwrap();
    let r = assemble(&[read], &options()).unwrap();
    assert!(!r.complete && r.reads[0].error.is_some());
    assert!(r.contigs.is_empty());
}
fn scf(calls: &str, scores: &[u8]) -> Vec<u8> {
    let n = calls.len();
    let samples = n * 3 + 2;
    let offset = 128 + samples * 8;
    let mut bytes = vec![0; 128];
    bytes[..4].copy_from_slice(b".scf");
    for (at, v) in [
        (4, samples),
        (8, 128),
        (12, n),
        (24, offset),
        (28, 0),
        (32, offset + n * 12),
        (40, 2),
        (44, 2),
    ] {
        bytes[at..at + 4].copy_from_slice(&(v as u32).to_be_bytes());
    }
    bytes[36..40].copy_from_slice(b"2.02");
    bytes.resize(offset, 0);
    for (i, b) in calls.bytes().enumerate() {
        bytes.extend_from_slice(&((i * 3 + 1) as u32).to_be_bytes());
        for x in b"ACGT" {
            bytes.push(if *x == b { scores[i] } else { 0 });
        }
        bytes.push(b);
        bytes.extend_from_slice(&[0; 3]);
    }
    bytes
}
#[test]
fn reverse_trace_views_keep_original_call_samples_and_cannot_fake_opposite_strands() {
    let g = genome(40);
    let raw = scf(&g, &vec![40; 40]);
    let trace = traces::import_trace(&raw).unwrap();
    let a = ReadInput {
        id: "stored".into(),
        molecule: trace.attached_molecule().unwrap(),
    };
    let b = ReadInput {
        id: "display-reversed".into(),
        molecule: trace
            .reverse_complement()
            .unwrap()
            .attached_molecule()
            .unwrap(),
    };
    let mut o = options();
    o.require_both_strands = true;
    let r = assemble(&[a, b], &o).unwrap();
    let c = best(&r);
    assert_eq!(c.placements.len(), 2, "{r:#?}");
    assert!(c.sequence.bytes().all(|b| b == b'N'));
    for p in &c.placements {
        for call in &p.calls {
            let obs = c.columns[call.column]
                .observations
                .iter()
                .find(|x| x.read_index == p.read_index)
                .unwrap();
            assert_eq!(
                obs.stored_sample_index,
                Some((obs.stored_call_index * 3 + 1) as u32)
            );
            assert_eq!(
                p.reverse ^ r.reads[p.read_index].source_reverse,
                c.placements[0].reverse ^ r.reads[c.placements[0].read_index].source_reverse
            );
        }
    }
}
#[test]
fn quality_windows_trim_source_indices_without_deleting_evidence_or_realigning_scores() {
    let g = genome(60);
    let mut score = vec![0; 60];
    score[7..51].fill(40);
    let read = fastq("trim", &g, &score);
    let mut o = options();
    o.trim.quality_windows = true;
    o.trim.window_size = 1;
    o.trim.minimum_passing = 1;
    let r = assemble(&[read], &o).unwrap();
    assert_eq!((r.reads[0].trim_start, r.reads[0].trim_end), (7, 51));
    let p = &r.contigs[0].placements[0];
    assert_eq!(p.calls[0].call_index, 7);
    assert_eq!(p.calls[0].stored_call_index, 7);
    assert_eq!(r.contigs[0].sequence, &g[7..51]);
}
#[test]
fn rich_source_and_nested_unknown_fields_are_rejected_before_expensive_work() {
    let mut read = plain("a", &genome(30));
    read.molecule.original_genbank = Some("x".repeat(16 * 1024 * 1024));
    assert!(matches!(
        batch_revision(&[read], &options()),
        Err(Error::TooLarge)
    ));
    let mut read = plain("a", &genome(30));
    let mut loc = Location::range(0, 1);
    for _ in 0..35 {
        loc = Location::Complement(Box::new(loc));
    }
    read.molecule.features.push(Feature {
        kind: "source".into(),
        location: loc,
        qualifiers: vec![],
        needs_review: false,
    });
    assert!(matches!(
        batch_revision(&[read], &options()),
        Err(Error::TooLarge)
    ));
    let mut value = serde_json::to_value(options()).unwrap();
    value["trim"]["confidence_probability"] = serde_json::json!(0.99);
    assert!(serde_json::from_value::<ContigOptions>(value).is_err());
}
#[test]
fn three_tiled_sanger_reads_with_low_quality_indels_and_reverse_read_reconstruct_observed_template()
{
    let mut x = 17u32;
    let mut g = (0..96)
        .map(|_| {
            x = x.wrapping_mul(1664525).wrapping_add(1013904223);
            b"ACGT"[((x >> 24) & 3) as usize] as char
        })
        .collect::<String>();
    g.replace_range(..1, "A");
    g.replace_range(95.., "C");
    let mut middle = g[22..74].as_bytes().to_vec();
    middle[14] = if middle[14] == b'A' { b'C' } else { b'A' };
    let inserted = *b"ACGT"
        .iter()
        .find(|b| **b != g.as_bytes()[44] && **b != g.as_bytes()[45])
        .unwrap();
    middle.insert(23, inserted);
    let mut scores = vec![40; middle.len()];
    scores[14] = 0;
    scores[23] = 0;
    let rc = Sequence::dna(&g[50..])
        .unwrap()
        .reverse_complement()
        .unwrap();
    let reads = vec![
        fastq("left", &g[..50], &vec![40; 50]),
        fastq("middle", std::str::from_utf8(&middle).unwrap(), &scores),
        fastq("reverse", rc.as_str(), &vec![40; 46]),
    ];
    let r = assemble(&reads, &options()).unwrap();
    let c = best(&r);
    assert!(r.complete, "{:?}", r.limits);
    assert_eq!(c.placements.len(), 3, "{r:#?}");
    assert!(matches_gauge(&c.sequence, &g), "{}", c.sequence);
    assert!(c.complete && c.layout_unique);
    assert!(
        c.columns
            .iter()
            .any(|c| c.state == ConsensusState::SupportedGap)
    );
    assert_eq!(
        c.placements.iter().map(|p| p.calls.len()).sum::<usize>(),
        reads
            .iter()
            .map(|r| r.molecule.sequence.len())
            .sum::<usize>()
    );
}
