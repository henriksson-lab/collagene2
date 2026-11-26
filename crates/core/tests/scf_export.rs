//! Hand-composed v2 source arrays and independent v3 byte expectations.
use collagene_core::{
    scf_export::{self, LossPolicy, Options, SamplePrecision},
    traces::{self, QualityEvidence, Trace, TraceFormat},
};
fn word(bytes: &mut [u8], at: usize, value: u32) {
    bytes[at..at + 4].copy_from_slice(&value.to_be_bytes());
}
fn fixture(wide: bool, extras: bool, code: u32, calls: &[u8; 4]) -> Vec<u8> {
    let channels: [[u16; 8]; 4] = if wide {
        [
            [0, 65535, 0, 32768, 0, 1, 65534, 5],
            [20, 21, 22, 23, 24, 25, 26, 27],
            [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000],
            [65535, 0, 65535, 0, 65535, 0, 65535, 0],
        ]
    } else {
        [
            [0, 255, 0, 128, 0, 1, 254, 5],
            [20, 21, 22, 23, 24, 25, 26, 27],
            [8, 8, 8, 8, 8, 8, 8, 8],
            [255, 0, 255, 0, 255, 0, 255, 0],
        ]
    };
    let precision = if wide { 2 } else { 1 };
    let base_offset = 128 + 32 * precision;
    let comment = if extras {
        b"NAME=actual\nMACH=retained instrument\n\0".as_slice()
    } else {
        b"NAME=actual\n\0".as_slice()
    };
    let comment_offset = base_offset + 48;
    let mut bytes = vec![0; comment_offset + comment.len()];
    bytes[..4].copy_from_slice(b".scf");
    word(&mut bytes, 4, 8);
    word(&mut bytes, 8, 128);
    word(&mut bytes, 12, 4);
    word(&mut bytes, 24, base_offset as u32);
    word(&mut bytes, 28, comment.len() as u32);
    word(&mut bytes, 32, comment_offset as u32);
    bytes[36..40].copy_from_slice(b"2.02");
    word(&mut bytes, 40, precision as u32);
    word(&mut bytes, 44, code);
    let mut at = 128;
    for sample in 0..8 {
        for channel in &channels {
            let value = channel[sample];
            if precision == 1 {
                bytes[at] = value as u8;
            } else {
                bytes[at..at + 2].copy_from_slice(&value.to_be_bytes());
            }
            at += precision;
        }
    }
    for (i, peak) in [0, 2, 4, 7].into_iter().enumerate() {
        let at = base_offset + i * 12;
        word(&mut bytes, at, peak);
        bytes[at + 4..at + 8].copy_from_slice(&[
            1 + i as u8,
            11 + i as u8,
            21 + i as u8,
            31 + i as u8,
        ]);
        bytes[at + 8] = calls[i];
    }
    bytes[comment_offset..].copy_from_slice(comment);
    if extras {
        let at = bytes.len();
        word(&mut bytes, 48, 4);
        word(&mut bytes, 52, at as u32);
        bytes.extend_from_slice(&[99, 88, 77, 66]);
        word(&mut bytes, 16, 2);
        bytes[56] = 9;
        bytes[base_offset + 9] = 7;
        bytes.push(42);
    }
    bytes
}
fn parsed(wide: bool) -> Trace {
    traces::import_trace(&fixture(wide, false, 2, b"ACRN")).unwrap()
}
fn same_evidence(left: &Trace, right: &Trace) {
    assert_eq!(left.sequence, right.sequence);
    assert_eq!(left.called_symbols, right.called_symbols);
    assert_eq!(left.peak_positions, right.peak_positions);
    assert_eq!(left.channels, right.channels);
    assert_eq!(left.qualities, right.qualities);
}
#[test]
fn v3_eight_bit_second_differences_and_channel_quality_layout_match_literal_bytes() {
    let trace = parsed(false);
    let original = trace.original_bytes.clone();
    let out = scf_export::export(&trace, &Options::default()).unwrap();
    assert_eq!(out.sample_size, 1);
    assert!(out.losses.is_empty());
    assert_eq!(&out.bytes[36..40], b"3.00");
    // A-channel original [0,255,0,128,0,1,254,5], second differences
    // modulo256, calculated independently of the production encoder.
    assert_eq!(&out.bytes[128..136], &[0, 255, 2, 127, 0, 129, 252, 10]);
    let b = 160;
    assert_eq!(
        &out.bytes[b..b + 16],
        &[0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 4, 0, 0, 0, 7]
    );
    assert_eq!(
        &out.bytes[b + 16..b + 32],
        &[1, 2, 3, 4, 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34]
    );
    assert_eq!(&out.bytes[b + 32..b + 36], b"ACRN");
    assert_eq!(&out.bytes[b + 36..b + 48], &[0; 12]);
    let restored = traces::import_trace(&out.bytes).unwrap();
    same_evidence(&trace, &restored);
    assert_eq!(trace.original_bytes, original);
}
#[test]
fn unsigned_sixteen_bit_wraps_are_lossless_and_explicit_eight_bit_refuses_rescaling() {
    let trace = parsed(true);
    let out = scf_export::export(&trace, &Options::default()).unwrap();
    assert_eq!(out.sample_size, 2);
    assert_eq!(&out.bytes[128..136], &[0, 0, 255, 255, 0, 2, 127, 255]);
    same_evidence(&trace, &traces::import_trace(&out.bytes).unwrap());
    assert!(
        scf_export::export(
            &trace,
            &Options {
                precision: SamplePrecision::EightBit,
                ..Default::default()
            }
        )
        .is_err()
    );
    let small = parsed(false);
    let wide = scf_export::export(
        &small,
        &Options {
            precision: SamplePrecision::SixteenBit,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(wide.sample_size, 2);
    same_evidence(&small, &traces::import_trace(&wide.bytes).unwrap());
}
#[test]
fn reverse_complement_swaps_signal_and_all_four_accuracy_vectors_and_retains_no_calls() {
    let trace = traces::import_trace(&fixture(false, false, 0, b"AcG-")).unwrap();
    let view = trace.reverse_complement().unwrap();
    let out = scf_export::export(&view, &Options::default()).unwrap();
    let restored = traces::import_trace(&out.bytes).unwrap();
    same_evidence(&view, &restored);
    assert_eq!(restored.called_symbols, b"-CgT");
    assert_eq!(restored.peak_positions, vec![0, 3, 5, 7]);
    assert_eq!(restored.channels.a, vec![0, 255, 0, 255, 0, 255, 0, 255]);
    let Some(QualityEvidence::ScfAccuracy(q)) = restored.qualities else {
        panic!("missing unchanged accuracy vectors")
    };
    assert_eq!(q.a, vec![34, 33, 32, 31]);
    assert_eq!(q.t, vec![4, 3, 2, 1]);
    assert!(matches!(
        restored.format,
        TraceFormat::Scf { code_set: 0, .. }
    ));
    assert_eq!(trace.original_bytes, view.original_bytes);
}
#[test]
fn documentary_omissions_require_explicit_policy_and_leave_original_source_unchanged() {
    let trace = traces::import_trace(&fixture(false, true, 2, b"ACRN")).unwrap();
    let before = trace.clone();
    assert!(scf_export::export(&trace, &Options::default()).is_err());
    let out = scf_export::export(
        &trace,
        &Options {
            loss_policy: LossPolicy::AllowDocumentaryOmissions,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        out.losses
            .iter()
            .map(|l| l.field.as_str())
            .collect::<Vec<_>>(),
        vec![
            "comments",
            "private_data",
            "obsolete_clips",
            "header_reserved",
            "base_reserved",
            "unclassified_container_bytes"
        ]
    );
    same_evidence(&trace, &traces::import_trace(&out.bytes).unwrap());
    assert_eq!(trace, before);
}
#[test]
fn fabricated_signal_missing_quality_changed_calls_and_invalid_dimensions_are_not_exportable() {
    let original = parsed(false);
    let mut changed = original.clone();
    changed.channels.a[0] = 123;
    assert!(scf_export::export(&changed, &Options::default()).is_err());
    let mut negative = original.clone();
    negative.channels.c[0] = -1;
    assert!(scf_export::export(&negative, &Options::default()).is_err());
    let mut missing = original.clone();
    missing.qualities = None;
    assert!(scf_export::export(&missing, &Options::default()).is_err());
    let mut calls = original.clone();
    calls.sequence = collagene_core::Sequence::dna("ACGT").unwrap();
    calls.called_symbols = b"ACGT".to_vec();
    assert!(scf_export::export(&calls, &Options::default()).is_err());
    let mut bad = original;
    bad.channels.t.pop();
    assert!(scf_export::export(&bad, &Options::default()).is_err());
}

#[test]
fn coincident_peaks_and_code_four_are_retained_and_v3_reexport_is_deterministic() {
    let mut raw = fixture(false, false, 4, b"ACNN");
    // Original v2 base section begins after32actualsamplebytes. Retain equal
    // called peaks rather than dropping a base or guessing another peak.
    word(&mut raw, 160 + 2 * 12, 2);
    let trace = traces::import_trace(&raw).unwrap();
    let first = scf_export::export(&trace, &Options::default()).unwrap();
    let imported = traces::import_trace(&first.bytes).unwrap();
    assert_eq!(imported.peak_positions, vec![0, 2, 2, 7]);
    assert!(matches!(
        imported.format,
        TraceFormat::Scf { code_set: 4, .. }
    ));
    same_evidence(&trace, &imported);
    let second = scf_export::export(&imported, &Options::default()).unwrap();
    assert_eq!(first.bytes, second.bytes);
}
