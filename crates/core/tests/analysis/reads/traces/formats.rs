use collagene_core::review::reads::traces::*;
use collagene_model::Error;

#[derive(Clone)]
struct Field {
    tag: [u8; 4],
    number: u32,
    kind: u16,
    size: u16,
    count: u32,
    data: Vec<u8>,
}
fn field(tag: &[u8; 4], number: u32, kind: u16, size: u16, data: Vec<u8>) -> Field {
    Field {
        tag: *tag,
        number,
        kind,
        size,
        count: (data.len() / size as usize) as u32,
        data,
    }
}
fn shorts(values: &[i16]) -> Vec<u8> {
    values.iter().flat_map(|n| n.to_be_bytes()).collect()
}
fn abif_fixture(calls: &[u8], quality: bool) -> Vec<u8> {
    let positions: Vec<i16> = if calls.len() == 4 {
        vec![0, 2, 4, 7]
    } else {
        vec![0, 2, 4, 6, 7]
    };
    let mut fields = vec![
        field(b"FWO_", 1, 2, 1, b"TGCA".to_vec()),
        field(b"DATA", 9, 4, 2, shorts(&[-1, 0, 5, 7, 9, 11, 13, 15])),
        field(
            b"DATA",
            10,
            4,
            2,
            shorts(&[100, 101, 102, 103, 104, 105, 106, 107]),
        ),
        field(
            b"DATA",
            11,
            4,
            2,
            shorts(&[200, 201, 202, 203, 204, 205, 206, 207]),
        ),
        field(
            b"DATA",
            12,
            4,
            2,
            shorts(&[300, 301, 302, 303, 304, 305, 306, 307]),
        ),
        field(b"PBAS", 2, 2, 1, calls.to_vec()),
        field(b"PLOC", 2, 4, 2, shorts(&positions)),
        field(b"SMPL", 1, 18, 1, b"\x03syn".to_vec()),
        field(b"RevC", 1, 4, 2, shorts(&[0])),
    ];
    if quality {
        fields.push(field(
            b"PCON",
            2,
            2,
            1,
            (0..calls.len()).map(|i| (i * 10 + 10) as u8).collect(),
        ));
    }
    let mut bytes = vec![0; 128 + 28 * fields.len()];
    bytes[..4].copy_from_slice(b"ABIF");
    bytes[4..6].copy_from_slice(&101u16.to_be_bytes());
    bytes[6..10].copy_from_slice(b"tdir");
    bytes[10..14].copy_from_slice(&1u32.to_be_bytes());
    bytes[14..16].copy_from_slice(&1023u16.to_be_bytes());
    bytes[16..18].copy_from_slice(&28u16.to_be_bytes());
    bytes[22..26].copy_from_slice(&((fields.len() * 28) as u32).to_be_bytes());
    bytes[18..22].copy_from_slice(&(fields.len() as u32).to_be_bytes());
    bytes[26..30].copy_from_slice(&128u32.to_be_bytes());
    for (i, f) in fields.iter().enumerate() {
        let p = 128 + 28 * i;
        bytes[p..p + 4].copy_from_slice(&f.tag);
        bytes[p + 4..p + 8].copy_from_slice(&f.number.to_be_bytes());
        bytes[p + 8..p + 10].copy_from_slice(&f.kind.to_be_bytes());
        bytes[p + 10..p + 12].copy_from_slice(&f.size.to_be_bytes());
        bytes[p + 12..p + 16].copy_from_slice(&f.count.to_be_bytes());
        bytes[p + 16..p + 20].copy_from_slice(&(f.data.len() as u32).to_be_bytes());
        if f.data.len() <= 4 {
            bytes[p + 20..p + 20 + f.data.len()].copy_from_slice(&f.data);
        } else {
            let offset = bytes.len() as u32;
            bytes[p + 20..p + 24].copy_from_slice(&offset.to_be_bytes());
            bytes.extend_from_slice(&f.data);
        }
    }
    bytes
}
fn entry(bytes: &[u8], tag: &[u8; 4], number: u32) -> usize {
    let count = u32::from_be_bytes(bytes[18..22].try_into().unwrap()) as usize;
    (0..count)
        .map(|i| 128 + 28 * i)
        .find(|&p| {
            &bytes[p..p + 4] == tag
                && u32::from_be_bytes(bytes[p + 4..p + 8].try_into().unwrap()) == number
        })
        .unwrap()
}
fn payload(bytes: &[u8], p: usize) -> usize {
    let length = u32::from_be_bytes(bytes[p + 16..p + 20].try_into().unwrap());
    if length <= 4 {
        p + 20
    } else {
        u32::from_be_bytes(bytes[p + 20..p + 24].try_into().unwrap()) as usize
    }
}

#[test]
fn abif_preserves_signed_signals_and_uses_explicit_channel_order() {
    let bytes = abif_fixture(b"ACGTN", true);
    let trace = import_trace(&bytes).unwrap();
    assert_eq!(trace.name, "syn");
    assert_eq!(trace.sequence.as_str(), "ACGTN");
    assert_eq!(
        trace.channels.a,
        vec![300, 301, 302, 303, 304, 305, 306, 307]
    );
    assert_eq!(
        trace.channels.c,
        vec![200, 201, 202, 203, 204, 205, 206, 207]
    );
    assert_eq!(
        trace.channels.g,
        vec![100, 101, 102, 103, 104, 105, 106, 107]
    );
    assert_eq!(trace.channels.t, vec![-1, 0, 5, 7, 9, 11, 13, 15]);
    assert_eq!(trace.peak_positions, vec![0, 2, 4, 6, 7]);
    assert_eq!(
        trace.qualities,
        Some(QualityEvidence::AbifBasecallerScores(vec![
            10, 20, 30, 40, 50
        ]))
    );
    assert_eq!(trace.original_bytes, bytes);
    assert_eq!(trace.source_hash, blake3::hash(&bytes).to_hex().as_str());
    let molecule = trace.molecule().unwrap();
    assert_eq!(molecule.sequence, trace.sequence);
    assert!(molecule.metadata["trace_projection"].contains("Called sequence only"));
}
#[test]
fn abif_inline_values_use_high_order_bytes_without_offset_interpretation() {
    let trace = import_trace(&abif_fixture(b"ACGT", true)).unwrap();
    assert_eq!(trace.sequence.as_str(), "ACGT");
    assert_eq!(trace.name, "syn");
    assert_eq!(
        trace.qualities,
        Some(QualityEvidence::AbifBasecallerScores(vec![10, 20, 30, 40]))
    );
    let no_quality = import_trace(&abif_fixture(b"ACGT", false)).unwrap();
    assert!(no_quality.qualities.is_none());
    assert!(no_quality.warnings.iter().any(|w| w.contains("No PCON")));
}
#[test]
fn abif_rejects_unsupported_orientation_missing_order_and_duplicate_tags() {
    let mut bytes = abif_fixture(b"ACGT", true);
    let p = entry(&bytes, b"RevC", 1);
    bytes[p + 20..p + 22].copy_from_slice(&1u16.to_be_bytes());
    assert!(matches!(import_trace(&bytes), Err(Error::Unsupported(_))));
    let mut bytes = abif_fixture(b"ACGT", true);
    let p = entry(&bytes, b"FWO_", 1);
    bytes[p + 20..p + 24].copy_from_slice(b"TTCA");
    assert!(import_trace(&bytes).is_err());
    let mut bytes = abif_fixture(b"ACGT", true);
    let p = entry(&bytes, b"DATA", 10);
    bytes[p + 4..p + 8].copy_from_slice(&9u32.to_be_bytes());
    assert!(import_trace(&bytes).is_err());
    let mut bytes = abif_fixture(b"ACGT", true);
    bytes[4..6].copy_from_slice(&201u16.to_be_bytes());
    assert!(matches!(import_trace(&bytes), Err(Error::Unsupported(_))));
}
#[test]
fn abif_root_reader_uses_bounded_count_and_offset_not_ignored_writer_hints() {
    let mut bytes = abif_fixture(b"ACGT", true);
    // Vendor permits older libraries to reserve more directory space than entries require.
    bytes[22..26].copy_from_slice(&100000u32.to_be_bytes());
    assert_eq!(import_trace(&bytes).unwrap().sequence.as_str(), "ACGT");
    bytes[26..30].copy_from_slice(&127u32.to_be_bytes());
    assert!(import_trace(&bytes).is_err());
}

#[test]
fn abif_call_editions_are_never_mixed_with_other_peak_or_quality_editions() {
    let mut bytes = abif_fixture(b"ACGT", true);
    for tag in [b"PBAS", b"PLOC", b"PCON"] {
        let p = entry(&bytes, tag, 2);
        bytes[p + 4..p + 8].copy_from_slice(&1u32.to_be_bytes());
    }
    let trace = import_trace(&bytes).unwrap();
    assert_eq!(trace.metadata["abif_call_edition"], "1");
    let mut bytes = abif_fixture(b"ACGT", true);
    let p = entry(&bytes, b"PLOC", 2);
    bytes[p + 4..p + 8].copy_from_slice(&1u32.to_be_bytes());
    assert!(import_trace(&bytes).is_err());
    let mut bytes = abif_fixture(b"ACGT", true);
    let p = entry(&bytes, b"PCON", 2);
    bytes[p + 4..p + 8].copy_from_slice(&1u32.to_be_bytes());
    assert!(import_trace(&bytes).unwrap().qualities.is_none());
}

#[test]
fn abif_validates_entry_dimensions_offsets_and_peak_coordinates() {
    let bytes = abif_fixture(b"ACGTN", true);
    for cutoff in [0, 4, 127, 128, bytes.len() - 1] {
        assert!(import_trace(&bytes[..cutoff]).is_err());
    }
    let mut corrupt = bytes.clone();
    corrupt[18..22].copy_from_slice(&u32::MAX.to_be_bytes());
    assert!(import_trace(&corrupt).is_err());
    let mut corrupt = bytes.clone();
    let p = entry(&corrupt, b"PBAS", 2);
    corrupt[p + 20..p + 24].copy_from_slice(&u32::MAX.to_be_bytes());
    assert!(import_trace(&corrupt).is_err());
    let mut corrupt = bytes.clone();
    let p = entry(&corrupt, b"DATA", 9);
    corrupt[p + 10..p + 12].copy_from_slice(&1u16.to_be_bytes());
    assert!(import_trace(&corrupt).is_err());
    let mut corrupt = bytes.clone();
    let p = entry(&corrupt, b"PCON", 2);
    corrupt[p + 12..p + 16].copy_from_slice(&4u32.to_be_bytes());
    assert!(import_trace(&corrupt).is_err());
    let mut corrupt = bytes.clone();
    let p = entry(&corrupt, b"PLOC", 2);
    let offset = payload(&corrupt, p);
    corrupt[offset..offset + 2].copy_from_slice(&8u16.to_be_bytes());
    assert!(import_trace(&corrupt).is_err());
    let mut corrupt = bytes.clone();
    let p = entry(&corrupt, b"PLOC", 2);
    let offset = payload(&corrupt, p);
    corrupt[offset..offset + 2].copy_from_slice(&(-1i16).to_be_bytes());
    assert!(import_trace(&corrupt).is_err());
}

fn write32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_be_bytes());
}
fn scf_fixture(version: &[u8; 4], precision: u32, code_set: u32, calls: &[u8; 4]) -> Vec<u8> {
    // Explicit hand-calculated second differences, not an encoder that mirrors the reader.
    let raw8: [[u16; 5]; 4] = [
        [250, 15, 248, 251, 3],
        [0, 1, 2, 2, 2],
        [1, 1, 2, 4, 8],
        [255, 0, 3, 0, 252],
    ];
    let raw16: [[u16; 5]; 4] = [
        [65530, 15, 65528, 65531, 3],
        [0, 1, 2, 2, 2],
        [1, 1, 2, 4, 8],
        [65535, 0, 3, 0, 65532],
    ];
    let expected8: [[u16; 5]; 4] = [
        [250, 3, 4, 0, 255],
        [0, 1, 4, 9, 16],
        [1, 3, 7, 15, 31],
        [255, 254, 0, 2, 0],
    ];
    let expected16: [[u16; 5]; 4] = [
        [65530, 3, 4, 0, 65535],
        [0, 1, 4, 9, 16],
        [1, 3, 7, 15, 31],
        [65535, 65534, 0, 2, 0],
    ];
    let major3 = version[0] == b'3';
    let sample_bytes = 20 * precision as usize;
    let base_offset = 128;
    let sample_offset = 176;
    let comment_offset = sample_offset + sample_bytes;
    let comments = b"NAME=synthetic\nMACH=independent fixture\0";
    let private_offset = comment_offset + comments.len();
    let mut bytes = vec![0; private_offset + 3];
    bytes[..4].copy_from_slice(b".scf");
    write32(&mut bytes, 4, 5);
    write32(&mut bytes, 8, sample_offset as u32);
    write32(&mut bytes, 12, 4);
    write32(&mut bytes, 16, 1);
    write32(&mut bytes, 20, 3);
    write32(&mut bytes, 24, base_offset as u32);
    write32(&mut bytes, 28, comments.len() as u32);
    write32(&mut bytes, 32, comment_offset as u32);
    bytes[36..40].copy_from_slice(version);
    write32(&mut bytes, 40, precision);
    write32(&mut bytes, 44, code_set);
    write32(&mut bytes, 48, 3);
    write32(&mut bytes, 52, private_offset as u32);
    let values = if major3 {
        if precision == 1 { raw8 } else { raw16 }
    } else if precision == 1 {
        expected8
    } else {
        expected16
    };
    for (channel, samples) in values.iter().enumerate() {
        for (sample, &value) in samples.iter().enumerate() {
            let index = if major3 {
                channel * 5 + sample
            } else {
                sample * 4 + channel
            };
            let offset = sample_offset + index * precision as usize;
            if precision == 1 {
                bytes[offset] = value as u8;
            } else {
                bytes[offset..offset + 2].copy_from_slice(&value.to_be_bytes());
            }
        }
    }
    for (i, &position) in [0, 1, 3, 4].iter().enumerate() {
        write32(
            &mut bytes,
            base_offset + if major3 { i * 4 } else { i * 12 },
            position,
        );
        for channel in 0..4 {
            bytes[base_offset
                + if major3 {
                    4 * (4 + channel) + i
                } else {
                    i * 12 + 4 + channel
                }] = (channel * 10 + i + 10) as u8;
        }
        bytes[base_offset + if major3 { 4 * 8 + i } else { i * 12 + 8 }] = calls[i];
    }
    bytes[comment_offset..comment_offset + comments.len()].copy_from_slice(comments);
    bytes[private_offset..].copy_from_slice(&[1, 2, 3]);
    bytes
}
#[test]
fn scf_v2_interleaved_and_v3_modular_delta_signals_match_independent_values() {
    for version in [b"2.02", b"3.00"] {
        for precision in [1, 2] {
            let trace = import_trace(&scf_fixture(version, precision, 2, b"ACGT")).unwrap();
            assert_eq!(trace.name, "synthetic");
            assert_eq!(trace.sequence.as_str(), "ACGT");
            assert_eq!(trace.peak_positions, vec![0, 1, 3, 4]);
            assert_eq!(
                trace.channels.a,
                if precision == 1 {
                    vec![250, 3, 4, 0, 255]
                } else {
                    vec![65530, 3, 4, 0, 65535]
                }
            );
            assert_eq!(trace.channels.c, vec![0, 1, 4, 9, 16]);
            assert_eq!(trace.channels.g, vec![1, 3, 7, 15, 31]);
            assert_eq!(
                trace.channels.t,
                if precision == 1 {
                    vec![255, 254, 0, 2, 0]
                } else {
                    vec![65535, 65534, 0, 2, 0]
                }
            );
            assert_eq!(
                trace.qualities,
                Some(QualityEvidence::ScfAccuracy(TraceAccuracy {
                    a: vec![10, 11, 12, 13],
                    c: vec![20, 21, 22, 23],
                    g: vec![30, 31, 32, 33],
                    t: vec![40, 41, 42, 43]
                }))
            );
            assert_eq!(trace.metadata["scf_reported_left_clip_obsolete"], "1");
            assert_eq!(trace.metadata["scf_reported_right_clip_obsolete"], "3");
        }
    }
}
#[test]
fn scf_no_calls_and_accuracy_evidence_are_not_silently_recalibrated() {
    let trace = import_trace(&scf_fixture(b"2.00", 1, 0, b"A-GT")).unwrap();
    assert_eq!(trace.sequence.as_str(), "ANGT");
    assert_eq!(trace.called_symbols, b"A-GT");
    assert!(trace.warnings.iter().any(|w| w.contains("projected to N")));
    let n = import_trace(&scf_fixture(b"3.00", 1, 2, b"ACGN")).unwrap();
    assert_eq!(n.sequence.as_str(), "ACGN");
    assert!(matches!(n.qualities, Some(QualityEvidence::ScfAccuracy(_))));
    assert!(import_trace(&scf_fixture(b"3.00", 1, 0, b"ACGN")).is_err());
    assert!(import_trace(&scf_fixture(b"3.00", 1, 4, b"ACGR")).is_err());
    assert!(import_trace(&scf_fixture(b"3.00", 1, 1, b"ACGT")).is_err());
}
#[test]
fn scf_rejects_bad_versions_overlapping_sections_counts_and_peak_positions() {
    let bytes = scf_fixture(b"3.00", 2, 2, b"ACGT");
    for cutoff in [0, 4, 127, 128, bytes.len() - 1] {
        assert!(import_trace(&bytes[..cutoff]).is_err());
    }
    let mut bad = bytes.clone();
    bad[36..40].copy_from_slice(b"4.00");
    assert!(import_trace(&bad).is_err());
    let mut bad = bytes.clone();
    write32(&mut bad, 4, u32::MAX);
    assert!(import_trace(&bad).is_err());
    let mut bad = bytes.clone();
    write32(&mut bad, 40, 4);
    assert!(import_trace(&bad).is_err());
    let mut bad = bytes.clone();
    write32(&mut bad, 8, 128);
    assert!(import_trace(&bad).is_err());
    let mut bad = bytes.clone();
    write32(&mut bad, 8, u32::MAX);
    assert!(import_trace(&bad).is_err());
    let mut bad = bytes.clone();
    write32(&mut bad, 128, 5);
    assert!(import_trace(&bad).is_err());
    let mut bad = bytes.clone();
    write32(&mut bad, 132, 0);
    write32(&mut bad, 128, 1);
    assert!(import_trace(&bad).is_err());
}
#[test]
fn reverse_complement_transforms_calls_positions_channels_and_qualities_together() {
    for trace in [
        import_trace(&abif_fixture(b"ACGTN", true)).unwrap(),
        import_trace(&scf_fixture(b"3.00", 1, 2, b"ACGN")).unwrap(),
        import_trace(&scf_fixture(b"2.00", 1, 0, b"A-GT")).unwrap(),
    ] {
        let reversed = trace.reverse_complement().unwrap();
        assert_eq!(
            reversed.sequence,
            trace.sequence.reverse_complement().unwrap()
        );
        assert_eq!(
            reversed.channels.a,
            trace.channels.t.iter().rev().copied().collect::<Vec<_>>()
        );
        assert_eq!(
            reversed.peak_positions,
            trace
                .peak_positions
                .iter()
                .rev()
                .map(|&p| (trace.channels.len() - 1) as u32 - p)
                .collect::<Vec<_>>()
        );
        assert_eq!(reversed.original_bytes, trace.original_bytes);
        assert_eq!(reversed.source_hash, trace.source_hash);
        assert_eq!(reversed.reverse_complement().unwrap(), trace);
    }
}
#[test]
fn trace_validation_rejects_inconsistent_projection_and_original_evidence() {
    let mut trace = import_trace(&abif_fixture(b"ACGT", true)).unwrap();
    trace.channels.a.pop();
    assert!(trace.molecule().is_err());
    let mut trace = import_trace(&abif_fixture(b"ACGT", true)).unwrap();
    trace.original_bytes[0] ^= 1;
    assert!(trace.validate().is_err());
    let mut trace = import_trace(&abif_fixture(b"ACGT", true)).unwrap();
    trace.called_symbols[0] = b'C';
    assert!(trace.validate().is_err());
}

#[test]
fn native_attachment_round_trip_reconstructs_original_trace_and_orientation() {
    for trace in [
        import_trace(&abif_fixture(b"ACGTN", true)).unwrap(),
        import_trace(&scf_fixture(b"3.00", 2, 2, b"ACGR"))
            .unwrap()
            .reverse_complement()
            .unwrap(),
    ] {
        let molecule = trace.attached_molecule().unwrap();
        assert!(molecule.metadata.contains_key("trace_original_base64"));
        assert!(!molecule.metadata.contains_key("trace_projection"));
        let native = collagene_core::protocol::save_document(&molecule).unwrap();
        let restored = collagene_core::protocol::load_document(&native).unwrap();
        assert_eq!(trace_from_molecule(&restored).unwrap(), trace);
    }
}
#[test]
fn attached_evidence_never_silently_overrides_edits_or_forged_view_data() {
    let trace = import_trace(&abif_fixture(b"ACGTN", true)).unwrap();
    let attached = trace.attached_molecule().unwrap();
    let edited = attached.replace(0, 1, "T").unwrap();
    assert_eq!(original_trace_bytes(&edited).unwrap(), trace.original_bytes);
    assert!(
        trace_from_molecule(&edited)
            .unwrap_err()
            .to_string()
            .contains("current sequence differs")
    );
    let mut corrupted = attached.clone();
    corrupted
        .metadata
        .insert("trace_source_hash".into(), "forged".into());
    assert!(trace_from_molecule(&corrupted).is_err());
    let mut corrupted = attached.clone();
    corrupted
        .metadata
        .insert("trace_original_base64".into(), "??".into());
    assert!(trace_from_molecule(&corrupted).is_err());
    let mut corrupted = attached.clone();
    corrupted
        .metadata
        .insert("trace_orientation".into(), "unknown".into());
    assert!(trace_from_molecule(&corrupted).is_err());
    let mut modified = trace.clone();
    modified.channels.a[0] += 1;
    assert!(modified.attached_molecule().is_err());
    assert!(trace_from_molecule(&trace.molecule().unwrap()).is_err());
}

#[test]
fn attachment_size_includes_base64_and_native_document_overhead() {
    let mut bytes = abif_fixture(b"ACGT", true);
    // Unused source padding is legal. A 24 MiB attachment expands to 32 MiB of base64
    // before any called sequence or native metadata/header, so it cannot fit .cg2.
    bytes.resize(
        collagene_core::formats::sequence::MAX_IMPORT_BYTES / 4 * 3,
        0,
    );
    let trace = import_trace(&bytes).unwrap();
    assert!(matches!(trace.attached_molecule(), Err(Error::TooLarge)));
}
