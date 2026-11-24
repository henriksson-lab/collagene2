use collagene_read_alignments::*;
use std::io::{Cursor, Read, Write};
const HEADER: &str = "@HD\tVN:1.6\tSO:unknown\r\n@SQ\tSN:chr1\tLN:1000\r\n@RG\tID:literal\tSM:naive\r\n@CO\tliteral header text\r\n";
fn sam(line: &str) -> SamReader<Cursor<Vec<u8>>> {
    SamReader::new(
        Cursor::new(format!("{HEADER}{line}").into_bytes()),
        Limits::default(),
    )
    .unwrap()
}
fn bam_from_sam(text: &str) -> Vec<u8> {
    let mut r = sam(text);
    let mut w = BamWriter::new(Vec::new(), r.header(), Limits::default()).unwrap();
    while let Some(record) = r.next_record().unwrap() {
        w.write_record(&record).unwrap();
    }
    assert!(r.completion().is_some());
    w.finish().unwrap()
}
fn bgzf(raw: &[u8]) -> Vec<u8> {
    let mut w = BgzfWriter::new(Vec::new());
    w.write_all(raw).unwrap();
    w.finish().unwrap()
}
fn unbgzf(raw: &[u8]) -> Vec<u8> {
    let mut r = BgzfReader::new(Cursor::new(raw), Limits::default()).unwrap();
    let mut b = Vec::new();
    r.read_to_end(&mut b).unwrap();
    assert!(r.source_blake3().is_some());
    b
}
fn literal_bam(record: &[u8]) -> Vec<u8> {
    let text = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n";
    let mut b = b"BAM\x01".to_vec();
    b.extend((text.len() as i32).to_le_bytes());
    b.extend(text);
    b.extend(1i32.to_le_bytes());
    b.extend(5i32.to_le_bytes());
    b.extend(b"chr1\0");
    b.extend(1000i32.to_le_bytes());
    b.extend(record);
    bgzf(&b)
}
// Independently assembled little-endian record: r at chr1:11, 2M1I1M, ACGT, Q0/Q1/Q2/Q93.
fn literal_record() -> Vec<u8> {
    let mut b = Vec::new();
    b.extend(0i32.to_le_bytes());
    b.extend(10i32.to_le_bytes());
    b.extend([2, 60]);
    b.extend(4681u16.to_le_bytes());
    b.extend(3u16.to_le_bytes());
    b.extend(16u16.to_le_bytes());
    b.extend(4i32.to_le_bytes());
    b.extend((-1i32).to_le_bytes());
    b.extend((-1i32).to_le_bytes());
    b.extend(0i32.to_le_bytes());
    b.extend(b"r\0");
    for op in [2u32 << 4, 1 << 4 | 1, 1 << 4] {
        b.extend(op.to_le_bytes());
    }
    b.extend([0x12, 0x48, 0, 1, 2, 93]);
    b.extend(b"XXC\xff");
    b.extend(b"XYs");
    b.extend((-300i16).to_le_bytes());
    b.extend(b"XZI");
    b.extend(u32::MAX.to_le_bytes());
    b.extend(b"ZZZcustom\0");
    b.extend(b"Fff");
    b.extend(0x7fc00001u32.to_le_bytes());
    let mut out = (b.len() as i32).to_le_bytes().to_vec();
    out.extend(b);
    out
}
#[test]
fn literal_sam_exact_evidence_missing_and_all_standard_auxiliary_types() {
    let line = "r\t147\tchr1\t11\t255\t2M1I1M\t=\t31\t-24\tAcgT\t!\"#~\tXA:A:Q\tXI:i:4294967295\tXF:f:0.125\tXZ:Z:unknown note\tXH:H:01AF\tXB:B:c,-128,0,127\tXS:B:S,0,65535\tXU:B:I,0,4294967295\tXR:B:f,1.5,-2\r\n";
    let mut r = sam(line);
    let header = r.header().clone();
    let record = r.next_record().unwrap().unwrap();
    let f = record.fields();
    assert_eq!(f.flags, 147);
    assert_eq!(f.start, Some(10));
    assert_eq!(f.mate_start, Some(30));
    assert_eq!(f.mapping_quality, None);
    assert_eq!(f.sequence.as_deref(), Some("AcgT"));
    assert_eq!(f.qualities.as_deref(), Some([0, 1, 2, 93].as_slice()));
    assert_eq!(f.aux.len(), 9);
    assert_eq!(record.original_bytes(), line.as_bytes());
    assert!(r.completion().is_none());
    assert!(r.next_record().unwrap().is_none());
    let c = r.completion().unwrap();
    assert_eq!(c.records, 1);
    assert_eq!(
        c.source_blake3,
        blake3::hash(format!("{HEADER}{line}").as_bytes())
            .to_hex()
            .to_string()
    );
    let mut w = SamWriter::new(Vec::new(), &header, Limits::default()).unwrap();
    assert!(w.write_record(&record).unwrap().is_empty());
    assert_eq!(w.finish().unwrap(), format!("{HEADER}{line}").as_bytes());
}
#[test]
fn independently_encoded_bam_preserves_integer_widths_float_bits_and_exact_payload() {
    let original = literal_record();
    let bytes = literal_bam(&original);
    let mut r = BamReader::new(Cursor::new(&bytes), Limits::default()).unwrap();
    let header = r.header().clone();
    let record = r.next_record().unwrap().unwrap();
    assert_eq!(record.original_bytes(), original);
    assert_eq!(record.fields().sequence.as_deref(), Some("ACGT"));
    assert_eq!(
        record.fields().aux[0],
        Aux {
            tag: *b"XX",
            wire_type: b'C',
            value: AuxValue::Integer(255)
        }
    );
    assert_eq!(record.fields().aux[1].value, AuxValue::Integer(-300));
    assert_eq!(record.fields().aux[2].value, AuxValue::Integer(4294967295));
    assert_eq!(
        record.fields().aux[4].value,
        AuxValue::FloatBits(0x7fc00001)
    );
    assert!(record.virtual_offset().is_some());
    assert!(r.next_record().unwrap().is_none());
    assert_eq!(
        r.completion().unwrap().source_blake3,
        blake3::hash(&bytes).to_hex().to_string()
    );
    let mut w = BamWriter::new(Vec::new(), &header, Limits::default()).unwrap();
    assert_eq!(
        w.write_record(&record).unwrap(),
        vec![ExportNotice::BgzfFramingNotPreserved]
    );
    assert_eq!(unbgzf(&w.finish().unwrap()), unbgzf(&bytes));
    let mut sam = SamWriter::new(Vec::new(), &header, Limits::default()).unwrap();
    assert!(sam.write_record(&record).is_err());
    assert!(sam.finish().is_err());
}
#[test]
fn cross_format_quality_encoding_unknown_tags_and_loss_notices() {
    let text =
        "r\t0\tchr1\t1\t30\t4M\t*\t0\t0\tacgt\t!\"#~\tZX:i:-2147483648\tZY:Z:not interpreted\n";
    let bytes = bam_from_sam(text);
    let mut r = BamReader::new(Cursor::new(&bytes), Limits::default()).unwrap();
    let h = r.header().clone();
    let record = r.next_record().unwrap().unwrap();
    assert_eq!(record.fields().sequence.as_deref(), Some("ACGT"));
    assert_eq!(
        record.fields().qualities.as_deref(),
        Some([0, 1, 2, 93].as_slice())
    );
    let mut w = SamWriter::new(Vec::new(), &h, Limits::default()).unwrap();
    assert!(
        w.write_record(&record)
            .unwrap()
            .contains(&ExportNotice::AuxIntegerWidthNotPreserved)
    );
    let output = w.finish().unwrap();
    let mut round = SamReader::new(Cursor::new(output), Limits::default()).unwrap();
    assert_eq!(
        round.next_record().unwrap().unwrap().fields(),
        record.fields()
    );
}
#[cfg(feature = "interop-tests")]
#[test]
fn noodles_independently_reads_exported_bam_and_sam_fixture() {
    let text = "r\t0\tchr1\t11\t60\t2M1I1M\t*\t0\t0\tACGT\t!\"#~\tNM:i:1\n";
    let bytes = bam_from_sam(text);
    let mut r = noodles_bam::io::Reader::new(Cursor::new(&bytes));
    let header = r.read_header().unwrap();
    let mut record = noodles_bam::Record::default();
    assert!(r.read_record(&mut record).unwrap() > 0);
    assert_eq!(record.alignment_start().unwrap().unwrap().get(), 11);
    assert_eq!(record.sequence().len(), 4);
    let mut source = BamReader::new(Cursor::new(&bytes), Limits::default()).unwrap();
    let mut output = SamWriter::new(Vec::new(), source.header(), Limits::default()).unwrap();
    output
        .write_record(&source.next_record().unwrap().unwrap())
        .unwrap();
    let mut sr = noodles_sam::io::Reader::new(Cursor::new(output.finish().unwrap()));
    let sh = sr.read_header().unwrap();
    let mut rec = noodles_sam::alignment::RecordBuf::default();
    assert!(sr.read_record_buf(&sh, &mut rec).unwrap() > 0);
    assert_eq!(rec.alignment_start().unwrap().get(), 11);
    assert_eq!(header.reference_sequences().len(), 1);
}
#[test]
fn all_cigar_operations_remain_distinct_with_reverse_hard_clip_coordinates() {
    let mut r = sam("r\t2064\tchr1\t11\t30\t2H3S2M1I2=1X1D4N1P2S5H\t*\t0\t0\tACGTACGTACG\t*\n");
    let record = r.next_record().unwrap().unwrap();
    let map = record.coordinate_map(r.header()).unwrap();
    assert_eq!(map.reported_query_length(), 11);
    assert_eq!(map.implied_original_length(), Some(18));
    assert_eq!(map.runs().len(), 11);
    assert_eq!(
        map.query_position(0),
        QueryProjection::SoftClipped {
            implied_original_query: Some(15)
        }
    );
    assert_eq!(
        map.query_position(3),
        QueryProjection::Aligned {
            reference_position: 10,
            unwrapped_position: 10,
            assertion: CigarKind::Match,
            implied_original_query: Some(12)
        }
    );
    assert_eq!(
        map.query_position(5),
        QueryProjection::Inserted {
            reference_boundary: 12,
            implied_original_query: Some(10)
        }
    );
    assert_eq!(
        map.query_position(8),
        QueryProjection::Aligned {
            reference_position: 14,
            unwrapped_position: 14,
            assertion: CigarKind::Difference,
            implied_original_query: Some(7)
        }
    );
    assert_eq!(
        map.query_position(10),
        QueryProjection::SoftClipped {
            implied_original_query: Some(5)
        }
    );
    assert!(map.runs().iter().any(|r| r.kind == CigarKind::Deletion
        && r.reference_start == 15
        && r.reference_end == 16
        && r.query_start == r.query_end));
    assert!(
        map.runs()
            .iter()
            .any(|r| r.kind == CigarKind::Skip && r.reference_start == 16 && r.reference_end == 20)
    );
}
#[test]
fn unmapped_and_missing_sequence_never_invent_mapping_or_calls() {
    let mut r = sam("r\t260\tchr1\t11\t30\t4M\t*\t0\t0\t*\t*\n");
    let rec = r.next_record().unwrap().unwrap();
    assert!(rec.fields().sequence.is_none());
    let map = rec.coordinate_map(r.header()).unwrap();
    assert!(!map.usable_reported_mapping());
    assert_eq!(map.query_position(0), QueryProjection::Unavailable);
    assert!(!map.sequence_present());
    assert_eq!(map.reported_query_length(), 4);
}
#[test]
fn circular_coordinates_wrap_only_under_explicit_header_topology() {
    let source = "@SQ\tSN:circle\tLN:5\tTP:circular\nr\t0\tcircle\t4\t30\t4M\t*\t0\t0\tACGT\t*\n";
    let mut r = SamReader::new(Cursor::new(source), Limits::default()).unwrap();
    let rec = r.next_record().unwrap().unwrap();
    let map = rec.coordinate_map(r.header()).unwrap();
    assert_eq!(
        map.query_position(2),
        QueryProjection::Aligned {
            reference_position: 0,
            unwrapped_position: 5,
            assertion: CigarKind::Match,
            implied_original_query: Some(2)
        }
    );
    let source = source.replace("\tTP:circular", "");
    let mut r = SamReader::new(Cursor::new(source), Limits::default()).unwrap();
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(
        rec.coordinate_map(r.header()).unwrap().query_position(0),
        QueryProjection::Unavailable
    );
}
#[test]
fn malformed_records_are_poisoned_without_completion() {
    for record in [
        "r\t0\tchr1\t1\t30\t3M\t*\t0\t0\tACGT\t*\n",
        "r\t0\tchr1\t1\t30\t4M\t*\t0\t0\tACGT\t!!!\n",
        "r\t0\tchr1\t1\t30\t4M\t*\t0\t0\tACGT\t*\tXX:i:1\tXX:i:2\n",
        "r\t0\tchr1\t1\t30\t1M1H3M\t*\t0\t0\tACGT\t*\n",
        "r\t0\tchr1\t1\t30\t4B\t*\t0\t0\tACGT\t*\n",
    ] {
        let mut r = sam(record);
        assert!(r.next_record().is_err(), "{record}");
        assert!(matches!(r.next_record(), Err(Error::Poisoned)));
        assert!(r.completion().is_none());
    }
}
#[test]
fn bgzf_partial_headers_oversized_isize_crc_and_final_eof_are_strict() {
    let source = bgzf(b"literal bytes");
    for count in [1, 2, 11, 17, source.len() - 1, source.len() - 28] {
        let mut r = BgzfReader::new(Cursor::new(&source[..count]), Limits::default()).unwrap();
        let mut sink = Vec::new();
        assert!(r.read_to_end(&mut sink).is_err(), "truncated {count}");
        assert!(r.source_blake3().is_none());
    }
    let size = u16::from_le_bytes([source[16], source[17]]) as usize + 1;
    for corrupt in [size - 8, size - 4] {
        let mut bytes = source.clone();
        if corrupt == size - 4 {
            bytes[corrupt..corrupt + 4].copy_from_slice(&65537u32.to_le_bytes());
        } else {
            bytes[corrupt] ^= 1;
        }
        let mut r = BgzfReader::new(Cursor::new(bytes), Limits::default()).unwrap();
        assert!(r.read_to_end(&mut Vec::new()).is_err());
        assert!(r.source_blake3().is_none());
    }
}
#[test]
fn embedded_empty_bgzf_blocks_do_not_signal_eof() {
    let mut bytes = bgzf(b"one");
    bytes.extend(bgzf(b"two"));
    assert_eq!(unbgzf(&bytes), b"onetwo");
}
#[test]
fn record_and_stream_caps_fail_before_false_completeness() {
    let l = Limits {
        max_records: 1,
        ..Limits::default()
    };
    let source =
        format!("{HEADER}r\t4\t*\t0\t0\t*\t*\t0\t0\tA\t*\ns\t4\t*\t0\t0\t*\t*\t0\t0\tT\t*\n");
    let mut r = SamReader::new(Cursor::new(source), l).unwrap();
    assert!(r.next_record().unwrap().is_some());
    assert!(r.next_record().is_err());
    assert!(r.completion().is_none());
    let bytes = literal_bam(&literal_record());
    let l = Limits {
        max_record_bytes: 32,
        ..Limits::default()
    };
    let mut r = BamReader::new(Cursor::new(bytes), l).unwrap();
    assert!(r.next_record().is_err());
    assert!(r.completion().is_none());
}
#[test]
fn malformed_declared_bam_lengths_and_reference_identity_are_rejected() {
    let mut raw = unbgzf(&literal_bam(&literal_record()));
    raw[4..8].copy_from_slice(&i32::MAX.to_le_bytes());
    assert!(BamReader::new(Cursor::new(bgzf(&raw)), Limits::default()).is_err());
    let mut raw = literal_record();
    raw[4..8].copy_from_slice(&2i32.to_le_bytes());
    let mut r = BamReader::new(Cursor::new(literal_bam(&raw)), Limits::default()).unwrap();
    assert!(r.next_record().is_err());
    let mut r = sam("r\t0\tchr1\t1\t30\t4M\t*\t0\t0\tACGT\t*\n");
    let rec = r.next_record().unwrap().unwrap();
    let other = Header::from_sam(b"@SQ\tSN:chr1\tLN:1001\n", &Limits::default()).unwrap();
    assert!(rec.coordinate_map(&other).is_err());
    let mut w = BamWriter::new(Vec::new(), &other, Limits::default()).unwrap();
    assert!(w.write_record(&rec).is_err());
    assert!(w.finish().is_err());
}
#[test]
fn cram_bcf_and_ordinary_gzip_are_not_misidentified_as_sam_bam() {
    assert!(matches!(detect(b"CRAM\x03\0"), Err(Error::Unsupported(_))));
    assert!(detect(b"BCF\x02\x02").is_err());
    assert!(detect(b"##fileformat=VCFv4.3").is_err());
    assert!(
        BamReader::new(
            Cursor::new([31, 139, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            Limits::default()
        )
        .is_err()
    );
}
#[test]
fn virtual_offsets_and_cross_block_records_remain_bounded() {
    let seq = "A".repeat(100_000);
    let bytes = bam_from_sam(&format!(
        "r\t0\tchr1\t1\t30\t100000M\t*\t0\t0\t{seq}\t*\ns\t4\t*\t0\t0\t*\t*\t0\t0\tA\t*\n"
    ));
    let mut r = BamReader::new(Cursor::new(bytes), Limits::default()).unwrap();
    let a = r.next_record().unwrap().unwrap();
    let b = r.next_record().unwrap().unwrap();
    assert_eq!(a.fields().sequence.as_ref().unwrap().len(), 100000);
    assert!(b.virtual_offset().unwrap().compressed > 0);
    assert!(
        b.virtual_offset().unwrap().packed().unwrap()
            > a.virtual_offset().unwrap().packed().unwrap()
    );
    assert!(r.next_record().unwrap().is_none());
}
#[test]
fn long_cigar_cg_transport_is_effective_geometry_without_losing_original_aux() {
    let n = 65536;
    let seq = "A".repeat(n);
    let cigar = "1M".repeat(n);
    let bytes = bam_from_sam(&format!(
        "long\t0\tchr1\t1\t30\t{cigar}\t*\t0\t0\t{seq}\t*\n"
    ));
    let mut r = BamReader::new(Cursor::new(bytes), Limits::default()).unwrap();
    let h = r.header().clone();
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(rec.fields().cigar.len(), n);
    assert!(rec.fields().aux.iter().any(|a| a.tag == *b"CG"));
    let mut w = SamWriter::new(Vec::new(), &h, Limits::default()).unwrap();
    assert!(
        w.write_record(&rec)
            .unwrap()
            .contains(&ExportNotice::LongCigarTransportConverted)
    );
    let output = w.finish().unwrap();
    let mut sr = SamReader::new(Cursor::new(output), Limits::default()).unwrap();
    let record = sr.next_record().unwrap().unwrap();
    assert_eq!(record.fields().cigar.len(), n);
    assert!(record.fields().aux.iter().all(|a| a.tag != *b"CG"));
}
#[test]
fn unsupported_sam_bases_require_explicit_conversion_instead_of_silent_n() {
    let mut r = sam("r\t0\tchr1\t1\t30\t2M\t*\t0\t0\tAU\t*\n");
    let mut w = BamWriter::new(Vec::new(), r.header(), Limits::default()).unwrap();
    assert!(w.write_record(&r.next_record().unwrap().unwrap()).is_err());
    assert!(w.finish().is_err());
}
#[test]
fn unterminated_sam_last_line_is_preserved_and_cannot_be_silently_concatenated() {
    let text = "r\t4\t*\t0\t0\t*\t*\t0\t0\tA\t*";
    let mut r = sam(text);
    let mut w = SamWriter::new(Vec::new(), r.header(), Limits::default()).unwrap();
    let rec = r.next_record().unwrap().unwrap();
    w.write_record(&rec).unwrap();
    assert_eq!(w.finish().unwrap(), format!("{HEADER}{text}").as_bytes());
}

#[test]
fn writer_cumulative_and_stricter_record_policies_cannot_finish_partial_exports() {
    let mut r = sam("r\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n");
    let rec = r.next_record().unwrap().unwrap();
    let limits = Limits {
        max_records: 1,
        ..Limits::default()
    };
    let mut w = BamWriter::new(Vec::new(), r.header(), limits.clone()).unwrap();
    w.write_record(&rec).unwrap();
    assert!(w.write_record(&rec).is_err());
    assert!(w.finish().is_err());
    let mut w = SamWriter::new(Vec::new(), r.header(), limits).unwrap();
    w.write_record(&rec).unwrap();
    assert!(w.write_record(&rec).is_err());
    assert!(w.finish().is_err());
    let limits = Limits {
        max_uncompressed_bytes: (HEADER.len() + rec.original_bytes().len() - 1) as u64,
        ..Limits::default()
    };
    let mut w = SamWriter::new(Vec::new(), r.header(), limits).unwrap();
    assert!(w.write_record(&rec).is_err());
    assert!(w.finish().is_err());
    let mut w = BgzfWriter::with_limit(Vec::new(), 27);
    w.write_all(b"x").unwrap();
    assert!(w.finish().is_err());
    let binary = bam_from_sam("r\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n");
    let mut reader = BamReader::new(Cursor::new(binary), Limits::default()).unwrap();
    let rec = reader.next_record().unwrap().unwrap();
    let limits = Limits {
        max_sequence_bases: 2,
        ..Limits::default()
    };
    let mut w = BamWriter::new(Vec::new(), reader.header(), limits).unwrap();
    assert!(w.write_record(&rec).is_err());
    assert!(w.finish().is_err());
}
#[test]
fn declared_small_sam_stream_limit_stops_without_reading_a_whole_long_line() {
    struct Count {
        bytes: Cursor<Vec<u8>>,
        count: std::rc::Rc<std::cell::Cell<usize>>,
    }
    impl Read for Count {
        fn read(&mut self, b: &mut [u8]) -> std::io::Result<usize> {
            let n = self.bytes.read(b)?;
            self.count.set(self.count.get() + n);
            Ok(n)
        }
    }
    let count = std::rc::Rc::new(std::cell::Cell::new(0));
    let inner = Count {
        bytes: Cursor::new(vec![b'A'; 100000]),
        count: count.clone(),
    };
    let limits = Limits {
        max_input_bytes: 10,
        ..Limits::default()
    };
    assert!(SamReader::new(std::io::BufReader::with_capacity(16, inner), limits).is_err());
    assert!(count.get() <= 16);
}
#[test]
fn bam_bin_handles_unmapped_and_large_reference_positions_explicitly() {
    let b = bam_from_sam("u\t4\t*\t0\t0\t*\t*\t0\t0\tA\t*\n");
    let mut r = BamReader::new(Cursor::new(b), Limits::default()).unwrap();
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(
        u16::from_le_bytes(rec.original_bytes()[14..16].try_into().unwrap()),
        4680
    );
}
#[test]
fn absent_dictionary_is_not_silently_invented_for_bam_conversion() {
    let mut r = SamReader::new(
        Cursor::new("r\t0\tunknown\t1\t30\t1M\t*\t0\t0\tA\t*\n"),
        Limits::default(),
    )
    .unwrap();
    let rec = r.next_record().unwrap().unwrap();
    assert!(!rec.coordinate_map(r.header()).unwrap().dictionary_bound());
    let mut w = BamWriter::new(Vec::new(), r.header(), Limits::default()).unwrap();
    assert!(w.write_record(&rec).is_err());
    assert!(w.finish().is_err());
}
#[test]
fn header_dictionary_alias_collisions_and_unbound_references_are_rejected() {
    for h in [
        b"@SQ\tSN:*bad\tLN:20\n".as_slice(),
        b"@SQ\tSN:a\tLN:20\tAN:b\n@SQ\tSN:b\tLN:30\n",
        b"@HD\tSO:coordinate\n",
        b"@SQ\tSN:a\tLN:20\tLN:21\n",
        b"@PG\tPN:aligner\n",
    ] {
        assert!(Header::from_sam(h, &Limits::default()).is_err());
    }
    let mut r = sam("r\t0\tmissing\t1\t30\t1M\t*\t0\t0\tA\t*\n");
    assert!(r.next_record().is_err());
}
#[test]
fn large_sam_is_consumed_recordwise_without_whole_file_buffers() {
    struct Repeated {
        header: Cursor<Vec<u8>>,
        line: Vec<u8>,
        offset: usize,
        left: u64,
        bytes: u64,
        largest_request: usize,
    }
    impl Read for Repeated {
        fn read(&mut self, b: &mut [u8]) -> std::io::Result<usize> {
            self.largest_request = self.largest_request.max(b.len());
            let n = self.header.read(b)?;
            if n > 0 {
                self.bytes += n as u64;
                return Ok(n);
            }
            if self.left == 0 {
                return Ok(0);
            }
            let n = b.len().min(self.line.len() - self.offset);
            b[..n].copy_from_slice(&self.line[self.offset..self.offset + n]);
            self.offset += n;
            if self.offset == self.line.len() {
                self.offset = 0;
                self.left -= 1;
            }
            self.bytes += n as u64;
            Ok(n)
        }
    }
    let line = format!(
        "stream\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\n",
        "A".repeat(100),
        "I".repeat(100)
    )
    .into_bytes();
    let count = 150000u64;
    let total = HEADER.len() as u64 + line.len() as u64 * count;
    assert!(total > 32 << 20);
    let source = Repeated {
        header: Cursor::new(HEADER.as_bytes().to_vec()),
        line,
        offset: 0,
        left: count,
        bytes: 0,
        largest_request: 0,
    };
    let mut r = SamReader::new(
        std::io::BufReader::with_capacity(4096, source),
        Limits::default(),
    )
    .unwrap();
    let mut n = 0;
    while let Some(record) = r.next_record().unwrap() {
        assert_eq!(record.fields().sequence.as_ref().unwrap().len(), 100);
        n += 1;
    }
    assert_eq!(n, count);
    let receipt = r.completion().unwrap();
    assert_eq!(receipt.records, count);
    assert_eq!(receipt.source_bytes, total);
}
#[test]
fn empty_aux_arrays_roundtrip_without_inventing_elements() {
    let b = bam_from_sam("r\t4\t*\t0\t0\t*\t*\t0\t0\tA\t*\tZZ:B:i\n");
    let mut r = BamReader::new(Cursor::new(b), Limits::default()).unwrap();
    let h = r.header().clone();
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(
        rec.fields().aux[0].value,
        AuxValue::Array {
            element_type: b'i',
            values: vec![]
        }
    );
    let mut w = SamWriter::new(Vec::new(), &h, Limits::default()).unwrap();
    w.write_record(&rec).unwrap();
    let mut r = SamReader::new(Cursor::new(w.finish().unwrap()), Limits::default()).unwrap();
    assert_eq!(
        r.next_record().unwrap().unwrap().fields().aux[0].value,
        rec.fields().aux[0].value
    );
}
#[test]
fn solitary_bam_q9_cannot_silently_become_sam_missing_quality() {
    let mut b = Vec::new();
    b.extend(0i32.to_le_bytes());
    b.extend(0i32.to_le_bytes());
    b.extend([2, 30]);
    b.extend(4681u16.to_le_bytes());
    b.extend(1u16.to_le_bytes());
    b.extend(0u16.to_le_bytes());
    b.extend(1i32.to_le_bytes());
    b.extend((-1i32).to_le_bytes());
    b.extend((-1i32).to_le_bytes());
    b.extend(0i32.to_le_bytes());
    b.extend(b"q\0");
    b.extend(16u32.to_le_bytes());
    b.extend([0x10, 9]);
    let mut record = (b.len() as i32).to_le_bytes().to_vec();
    record.extend(b);
    let mut r = BamReader::new(Cursor::new(literal_bam(&record)), Limits::default()).unwrap();
    let h = r.header().clone();
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(rec.fields().qualities.as_deref(), Some([9u8].as_slice()));
    let mut w = SamWriter::new(Vec::new(), &h, Limits::default()).unwrap();
    assert!(w.write_record(&rec).is_err());
    assert!(w.finish().is_err());
}
#[test]
fn all_bam_sequence_nibbles_preserve_reference_equal_without_inventing_called_dna() {
    let mut b = Vec::new();
    b.extend(0i32.to_le_bytes());
    b.extend(0i32.to_le_bytes());
    b.extend([2, 30]);
    b.extend(4681u16.to_le_bytes());
    b.extend(1u16.to_le_bytes());
    b.extend(0u16.to_le_bytes());
    b.extend(16i32.to_le_bytes());
    b.extend((-1i32).to_le_bytes());
    b.extend((-1i32).to_le_bytes());
    b.extend(0i32.to_le_bytes());
    b.extend(b"n\0");
    b.extend((16u32 << 4).to_le_bytes());
    b.extend([0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef]);
    b.extend([255; 16]);
    let mut record = (b.len() as i32).to_le_bytes().to_vec();
    record.extend(b);
    let mut r = BamReader::new(Cursor::new(literal_bam(&record)), Limits::default()).unwrap();
    let h = r.header().clone();
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(rec.fields().sequence.as_deref(), Some("=ACMGRSVTWYHKDBN"));
    assert_eq!(rec.reported_symbol(0), Some(SequenceSymbol::ReferenceEqual));
    for (i, b) in b"ACMGRSVTWYHKDBN".iter().enumerate() {
        assert_eq!(rec.reported_symbol(i + 1), Some(SequenceSymbol::Iupac(*b)));
    }
    assert!(rec.dna_calls().is_err());
    let mut w = SamWriter::new(Vec::new(), &h, Limits::default()).unwrap();
    w.write_record(&rec).unwrap();
    let output = w.finish().unwrap();
    let mut sam = SamReader::new(Cursor::new(output), Limits::default()).unwrap();
    let h = sam.header().clone();
    let rec2 = sam.next_record().unwrap().unwrap();
    assert_eq!(rec2.fields().sequence, rec.fields().sequence);
    assert!(rec2.dna_calls().is_err());
    let mut w = BamWriter::new(Vec::new(), &h, Limits::default()).unwrap();
    w.write_record(&rec2).unwrap();
    let mut reread = BamReader::new(Cursor::new(w.finish().unwrap()), Limits::default()).unwrap();
    assert_eq!(
        reread.next_record().unwrap().unwrap().fields().sequence,
        rec.fields().sequence
    );
}
#[test]
fn bam_missing_sequence_and_unresolved_sam_calls_have_no_shared_dna_conversion() {
    let mut b = Vec::new();
    b.extend((-1i32).to_le_bytes());
    b.extend((-1i32).to_le_bytes());
    b.extend([2, 255]);
    b.extend(4680u16.to_le_bytes());
    b.extend(0u16.to_le_bytes());
    b.extend(4u16.to_le_bytes());
    b.extend(0i32.to_le_bytes());
    b.extend((-1i32).to_le_bytes());
    b.extend((-1i32).to_le_bytes());
    b.extend(0i32.to_le_bytes());
    b.extend(b"m\0");
    let mut record = (b.len() as i32).to_le_bytes().to_vec();
    record.extend(b);
    let mut r = BamReader::new(Cursor::new(literal_bam(&record)), Limits::default()).unwrap();
    let rec = r.next_record().unwrap().unwrap();
    assert!(rec.fields().sequence.is_none());
    assert!(rec.fields().qualities.is_none());
    assert_eq!(rec.reported_symbol(0), None);
    assert!(rec.dna_calls().is_err());
    let mut r = sam("r\t4\t*\t0\t0\t*\t*\t0\t0\tA.U\t*\n");
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(
        rec.reported_symbol(1),
        Some(SequenceSymbol::Unresolved(b'.'))
    );
    assert_eq!(
        rec.reported_symbol(2),
        Some(SequenceSymbol::Unresolved(b'U'))
    );
    assert!(rec.dna_calls().is_err());
    let mut r = sam("r\t4\t*\t0\t0\t*\t*\t0\t0\tacgtrn\t*\n");
    let rec = r.next_record().unwrap().unwrap();
    let view = rec.dna_calls().unwrap();
    assert_eq!(view.as_str(), "acgtrn");
    assert!(view.case_normalization_required);
}
#[test]
fn header_reference_claims_cannot_be_substituted_with_same_name_and_length() {
    let mut r = sam("r\t0\tchr1\t1\t30\t1M\t*\t0\t0\tA\t*\n");
    let rec = r.next_record().unwrap().unwrap();
    let foreign = Header::from_sam(
        HEADER
            .replace("LN:1000", "LN:1000\tM5:01234567890123456789012345678901")
            .as_bytes(),
        &Limits::default(),
    )
    .unwrap();
    assert!(rec.coordinate_map(&foreign).is_err());
    let mut writer = BamWriter::new(Vec::new(), &foreign, Limits::default()).unwrap();
    assert!(writer.write_record(&rec).is_err());
    assert!(writer.finish().is_err());
}
#[test]
fn partial_bam_failure_progress_retains_emitted_count_but_no_completion_hash() {
    let mut bytes = literal_bam(&literal_record());
    bytes.truncate(bytes.len() - 28);
    let mut r = BamReader::new(Cursor::new(bytes), Limits::default()).unwrap();
    assert_eq!(r.progress().state, StreamState::Reading);
    assert!(r.next_record().unwrap().is_some());
    assert_eq!(r.progress().records_emitted, 1);
    assert!(r.next_record().is_err());
    assert_eq!(r.progress().state, StreamState::Failed);
    assert_eq!(r.progress().records_emitted, 1);
    assert!(r.completion().is_none());
}
#[test]
fn io_write_failure_poisons_both_export_paths() {
    struct Fail {
        left: usize,
    }
    impl Write for Fail {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            if self.left == 0 {
                return Err(std::io::Error::other("injected sink failure"));
            }
            let n = self.left.min(b.len());
            self.left -= n;
            Ok(n)
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let source = format!(
        "r\t4\t*\t0\t0\t100000M\t*\t0\t0\t{}\t*\n",
        "A".repeat(100000)
    );
    let mut r = sam(&source);
    let rec = r.next_record().unwrap().unwrap();
    let mut w = BamWriter::new(Fail { left: 10 }, r.header(), Limits::default()).unwrap();
    assert!(w.write_record(&rec).is_err());
    assert!(matches!(w.write_record(&rec), Err(Error::Poisoned)));
    assert!(w.finish().is_err());
    let mut w = SamWriter::new(
        Fail {
            left: HEADER.len() + 2,
        },
        r.header(),
        Limits::default(),
    )
    .unwrap();
    assert!(w.write_record(&rec).is_err());
    assert!(matches!(w.write_record(&rec), Err(Error::Poisoned)));
    assert!(w.finish().is_err());
}
#[test]
fn bam_edge_coordinate_cannot_overflow_sam_one_based_field() {
    let mut original = literal_record();
    original.truncate(original.len() - 7); // Remove the independent NaN fixture: only coordinate representability must fail.
    let size = (original.len() - 4) as i32;
    original[..4].copy_from_slice(&size.to_le_bytes());
    original[8..12].copy_from_slice(&i32::MAX.to_le_bytes());
    let mut r = BamReader::new(Cursor::new(literal_bam(&original)), Limits::default()).unwrap();
    let h = r.header().clone();
    let rec = r.next_record().unwrap().unwrap();
    assert_eq!(rec.fields().start, Some(i32::MAX as u32));
    let mut w = SamWriter::new(Vec::new(), &h, Limits::default()).unwrap();
    assert!(w.write_record(&rec).is_err());
    assert!(w.finish().is_err());
}
#[test]
fn valid_sam_read_names_are_not_confused_with_binary_format_magic() {
    for name in ["CRAM", "CRAM_read", "BCF", "BCF123", "##fileformat=VCF"] {
        let mut r = sam(&format!("{name}\t4\t*\t0\t0\t*\t*\t0\t0\tA\t*\n"));
        let rec = r.next_record().unwrap().unwrap();
        assert_eq!(rec.fields().name.as_deref(), Some(name));
        assert!(r.next_record().unwrap().is_none());
        assert!(r.completion().is_some());
    }
}
#[test]
fn header_line_caps_are_identical_for_constructors_readers_and_writers() {
    let loose = Limits::default();
    let bytes = b"@CO\t12345\n";
    let tight = Limits {
        max_record_bytes: 8,
        max_header_bytes: 100,
        ..Limits::default()
    };
    assert!(Header::from_sam(bytes, &tight).is_err());
    assert!(SamReader::new(Cursor::new(bytes), tight.clone()).is_err());
    let header = Header::from_sam(bytes, &loose).unwrap();
    assert!(SamWriter::new(Vec::new(), &header, tight.clone()).is_err());
    assert!(BamWriter::new(Vec::new(), &header, tight).is_err());
}
#[test]
fn tighter_bam_header_policy_checks_encoded_overhead_before_writing() {
    let header = Header::from_sam(b"@SQ\tSN:chr1\tLN:1000\n", &Limits::default()).unwrap();
    let tight = Limits {
        max_header_bytes: header.sam_bytes().len(),
        ..Limits::default()
    };
    let mut sink = Vec::new();
    assert!(BamWriter::new(&mut sink, &header, tight).is_err());
    assert!(sink.is_empty());
}
#[test]
fn compact_bam_array_expansion_is_capped_before_any_sam_record_is_written() {
    let mut original = literal_record();
    original.truncate(original.len() - 7);
    original.extend(b"XBBC");
    original.extend(200i32.to_le_bytes());
    original.extend([255u8; 200]);
    let size = (original.len() - 4) as i32;
    original[..4].copy_from_slice(&size.to_le_bytes());
    assert!(original.len() < 400);
    let mut r = BamReader::new(Cursor::new(literal_bam(&original)), Limits::default()).unwrap();
    let h = r.header().clone();
    let rec = r.next_record().unwrap().unwrap();
    let mut sink = Vec::new();
    let limits = Limits {
        max_record_bytes: 400,
        ..Limits::default()
    };
    let mut writer = SamWriter::new(&mut sink, &h, limits).unwrap();
    assert!(matches!(writer.write_record(&rec), Err(Error::Limit(_))));
    assert!(matches!(writer.write_record(&rec), Err(Error::Poisoned)));
    assert!(writer.finish().is_err());
    assert_eq!(sink, h.sam_bytes());
    let mut sink = Vec::new();
    let limits = Limits {
        max_uncompressed_bytes: h.sam_bytes().len() as u64 + 400,
        ..Limits::default()
    };
    let mut writer = SamWriter::new(&mut sink, &h, limits).unwrap();
    assert!(matches!(writer.write_record(&rec), Err(Error::Limit(_))));
    assert!(writer.finish().is_err());
    assert_eq!(sink, h.sam_bytes());
}
#[test]
fn bam_fixed_field_expansion_respects_a_tighter_output_record_cap() {
    let mut reader = sam("r\t4\t*\t0\t0\t*\t*\t0\t0\tA\t*\n");
    let rec = reader.next_record().unwrap().unwrap();
    assert!(rec.original_bytes().len() < 36);
    let mut sink = Vec::new();
    let limits = Limits {
        max_record_bytes: 36,
        ..Limits::default()
    };
    let mut writer = BamWriter::new(&mut sink, reader.header(), limits).unwrap();
    assert!(matches!(writer.write_record(&rec), Err(Error::Limit(_))));
    assert!(writer.finish().is_err());
    assert!(sink.is_empty());
}
