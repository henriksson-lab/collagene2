//! Literal BAI framing plus an independent CIGAR/span oracle.
use collagene_core::review::reads::{
    indexed_bam::IndexedBam,
    store::{Collection, Format, ImportLimits},
};
use collagene_genomic_index::{Format as IndexFormat, Limits, QueryLimits};
use collagene_read_alignments::{BamWriter, Limits as ReadLimits, SamReader};
use std::io::Cursor;

fn source() -> (Vec<u8>, Collection) {
    // Reported spans: [1,5), [5,12) with D/N, [20,21) for zero reference
    // consumption; a placed unmapped record is excluded from mapped selection.
    let sam = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\na\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII\nb\t0\tchr1\t6\t60\t1M2D3N1M\t*\t0\t0\tAT\tII\nc\t0\tchr1\t21\t60\t2I\t*\t0\t0\tGG\tII\nu\t4\tchr1\t22\t0\t*\t*\t0\t0\tA\tI\n";
    let mut reader = SamReader::new(Cursor::new(sam.as_slice()), ReadLimits::default()).unwrap();
    let mut writer = BamWriter::new(Vec::new(), reader.header(), ReadLimits::default()).unwrap();
    while let Some(record) = reader.next_record().unwrap() {
        writer.write_record(&record).unwrap();
    }
    let bytes = writer.finish().unwrap();
    let collection = Collection::import(&bytes, Format::Bam, &ImportLimits::default()).unwrap();
    (bytes, collection)
}
fn bai(source_bytes: usize, has_chunk: bool) -> Vec<u8> {
    // One reference, root bin 0 covering all positions, [0, source_end)
    // candidate chunk; optional linear index omitted. No library index writer.
    let mut b = b"BAI\x01".to_vec();
    b.extend(1i32.to_le_bytes());
    b.extend(i32::from(has_chunk).to_le_bytes());
    if has_chunk {
        b.extend(0u32.to_le_bytes());
        b.extend(1i32.to_le_bytes());
        b.extend(0u64.to_le_bytes());
        b.extend(((source_bytes as u64) << 16).to_le_bytes());
    }
    b.extend(0i32.to_le_bytes());
    b
}
#[test]
fn query_checks_actual_spans_and_cigar_reference_consumption() {
    let (bytes, collection) = source();
    let indexed = IndexedBam::bind(
        &collection,
        &bai(bytes.len(), true),
        IndexFormat::Bai,
        &Limits::default(),
    )
    .unwrap();
    for (a, b, expected) in [
        (0, 1, vec![]),
        (1, 2, vec![0]),
        (4, 6, vec![0, 1]),
        (7, 8, vec![1]),
        (10, 11, vec![1]),
        (20, 21, vec![2]),
        (21, 22, vec![]),
    ] {
        let selected = indexed
            .select(&collection, "chr1", a, b, &QueryLimits::default(), 100, 100)
            .unwrap();
        assert_eq!(selected.record_indices, expected);
        assert_eq!(selected.records_audited, 4);
        assert_eq!(selected.source_blake3, collection.source_blake3());
    }
    assert!(
        indexed
            .select(&collection, "chr1", 5, 5, &QueryLimits::default(), 100, 100)
            .is_err()
    );
    assert!(
        indexed
            .select(
                &collection,
                "chr1",
                0,
                101,
                &QueryLimits::default(),
                100,
                100
            )
            .is_err()
    );
    assert!(
        indexed
            .select(
                &collection,
                "other",
                0,
                10,
                &QueryLimits::default(),
                100,
                100
            )
            .is_err()
    );
    assert!(
        indexed
            .select(&collection, "chr1", 0, 100, &QueryLimits::default(), 3, 100)
            .is_err()
    );
    assert!(
        indexed
            .select(&collection, "chr1", 0, 100, &QueryLimits::default(), 100, 1)
            .is_err()
    );
}
#[test]
fn an_empty_or_out_of_source_index_cannot_claim_no_matching_records() {
    let (bytes, collection) = source();
    let missing = IndexedBam::bind(
        &collection,
        &bai(bytes.len(), false),
        IndexFormat::Bai,
        &Limits::default(),
    )
    .unwrap();
    assert!(
        missing
            .select(&collection, "chr1", 1, 2, &QueryLimits::default(), 100, 100)
            .unwrap_err()
            .to_string()
            .contains("omits")
    );
    let oversized = IndexedBam::bind(
        &collection,
        &bai(bytes.len() + 1, true),
        IndexFormat::Bai,
        &Limits::default(),
    )
    .unwrap();
    assert!(
        oversized
            .select(&collection, "chr1", 1, 2, &QueryLimits::default(), 100, 100)
            .is_err()
    );
    let sam = Collection::import(
        b"@SQ\tSN:chr1\tLN:100\n",
        Format::Sam,
        &ImportLimits::default(),
    )
    .unwrap();
    assert!(
        IndexedBam::bind(
            &sam,
            &bai(bytes.len(), true),
            IndexFormat::Bai,
            &Limits::default()
        )
        .is_err()
    );
    assert!(
        missing
            .select(&sam, "chr1", 1, 2, &QueryLimits::default(), 100, 100)
            .is_err()
    );
}

#[test]
fn bgzf_csi_uses_exact_transport_hash_and_checked_crc_eof_limits() {
    use std::io::Write;
    let (bytes, collection) = source();
    // Literal decoded CSI framing: min_shift14/depth5, no opaque auxiliary
    // bytes, one reference/root bin/one candidate chunk, no pseudo-bin.
    let mut csi = b"CSI\x01".to_vec();
    for n in [14i32, 5, 0, 1, 1] {
        csi.extend(n.to_le_bytes());
    }
    csi.extend(0u32.to_le_bytes()); // root bin
    csi.extend(0u64.to_le_bytes()); // loffset, not an absence assertion
    csi.extend(1i32.to_le_bytes());
    csi.extend(0u64.to_le_bytes());
    csi.extend(((bytes.len() as u64) << 16).to_le_bytes());
    let mut writer = collagene_read_alignments::BgzfWriter::new(Vec::new());
    writer.write_all(&csi).unwrap();
    let compressed = writer.finish().unwrap();
    let indexed = IndexedBam::bind(
        &collection,
        &compressed,
        IndexFormat::Csi,
        &Limits::default(),
    )
    .unwrap();
    assert_eq!(
        indexed.index_blake3(),
        blake3::hash(&compressed).to_hex().as_str()
    );
    assert_eq!(
        indexed
            .select(&collection, "chr1", 7, 8, &QueryLimits::default(), 100, 100)
            .unwrap()
            .record_indices,
        vec![1]
    );
    let mut bad_crc = compressed.clone();
    bad_crc[18] ^= 128;
    assert!(IndexedBam::bind(&collection, &bad_crc, IndexFormat::Csi, &Limits::default()).is_err());
    let mut trailing = compressed.clone();
    trailing.push(0);
    assert!(
        IndexedBam::bind(&collection, &trailing, IndexFormat::Csi, &Limits::default()).is_err()
    );
    assert!(
        IndexedBam::bind(
            &collection,
            &compressed[..compressed.len() - 1],
            IndexFormat::Csi,
            &Limits::default()
        )
        .is_err()
    );
    let small = Limits {
        source_bytes: compressed.len() - 1,
        ..Limits::default()
    };
    assert!(IndexedBam::bind(&collection, &compressed, IndexFormat::Csi, &small).is_err());
}

#[test]
fn fresh_bai_uses_actual_offsets_bins_and_sort_order_not_header_assertions() {
    use collagene_core::review::reads::indexed_bam::build_bai;
    let (_, collection) = source();
    let bai = build_bai(&collection, 1024).unwrap();
    // All four placed records are within the first 16kb bin, with contiguous
    // source offsets, hence one merged chunk in bin4681. No optional metadata.
    assert_eq!(bai.len(), 40);
    assert_eq!(u32::from_le_bytes(bai[12..16].try_into().unwrap()), 4681);
    assert!(build_bai(&collection, 39).is_err());
    let indexed =
        IndexedBam::bind(&collection, &bai, IndexFormat::Bai, &Limits::default()).unwrap();
    assert_eq!(
        indexed
            .select(
                &collection,
                "chr1",
                0,
                100,
                &QueryLimits::default(),
                100,
                100
            )
            .unwrap()
            .record_indices,
        vec![0, 1, 2]
    );

    // Header claims sorting; actual coordinates go backwards. No BAI emitted.
    let sam = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\na\t0\tchr1\t9\t60\t1M\t*\t0\t0\tA\tI\nb\t0\tchr1\t2\t60\t1M\t*\t0\t0\tT\tI\n";
    let mut reader = SamReader::new(Cursor::new(sam.as_slice()), ReadLimits::default()).unwrap();
    let mut writer = BamWriter::new(Vec::new(), reader.header(), ReadLimits::default()).unwrap();
    while let Some(record) = reader.next_record().unwrap() {
        writer.write_record(&record).unwrap();
    }
    let bytes = writer.finish().unwrap();
    let unsorted = Collection::import(&bytes, Format::Bam, &ImportLimits::default()).unwrap();
    assert!(
        build_bai(&unsorted, 1024)
            .unwrap_err()
            .to_string()
            .contains("not coordinate sorted")
    );
}
