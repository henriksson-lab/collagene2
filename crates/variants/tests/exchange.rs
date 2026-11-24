use collagene_variants::{Collection, ExportFormat, Format, Limits};
const VCF:&[u8]=b"##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Declared frequencies\">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Declared structural type\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth assertion\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\nchr1\t7\tv1\tA\tC,G\t01.50\tPASS\tAF=0.25,.\tGT:DP\t0|2:9\t./.:.\nchr1\t20\tbreakend\tG\tG]chr1:30]\t.\t.\tSVTYPE=BND\tGT\t1/0\t.\n";
#[test]
fn lexical_original_typed_headers_and_phased_missing_multiallelic_and_breakend_are_retained() {
    let c = Collection::import(VCF, Format::Vcf, &Limits::default()).unwrap();
    assert_eq!(c.len(), 2);
    assert_eq!(c.original_bytes(), VCF);
    assert_eq!(c.source_blake3(), blake3::hash(VCF).to_hex().as_str());
    assert_eq!(c.header().sample_names().len(), 2);
    let rows = c.page(0, 100).unwrap();
    assert_eq!(rows[0].position_1based, Some(7));
    assert!(rows[0].vcf.contains("0|2:9\t./.:."));
    assert!(rows[0].vcf.contains("C,G"));
    assert!(rows[1].vcf.contains("G]chr1:30]"));
    assert!(rows[1].vcf.ends_with("\t1/0\t.\n"));
    let normalized_vcf = c.export(ExportFormat::Vcf, 1 << 20).unwrap();
    assert!(normalized_vcf.bytes.ends_with(b"\t1/0\t.\n"));
    assert_eq!(
        c.export(ExportFormat::Original, VCF.len()).unwrap().bytes,
        VCF
    );
    let bcf = c.export(ExportFormat::Bcf, 1 << 20).unwrap();
    assert!(!bcf.notices.is_empty());
    let d = Collection::import(&bcf.bytes, Format::Bcf, &Limits::default()).unwrap();
    assert_eq!(d.len(), 2);
    assert_eq!(d.page(0, 100).unwrap(), rows);
    assert_eq!(
        d.export(ExportFormat::Original, 1 << 20).unwrap().bytes,
        bcf.bytes
    );
}
#[test]
fn limits_and_errors_never_publish_partial_collection_or_export() {
    let c = Collection::import(VCF, Format::Vcf, &Limits::default()).unwrap();
    assert!(c.page(0, 101).is_err());
    assert!(c.page(3, 1).is_err());
    assert!(c.export(ExportFormat::Vcf, 10).is_err());
    assert!(c.export(ExportFormat::Bcf, 10).is_err());
    assert_eq!(c.original_bytes(), VCF);
    let mut l = Limits::default();
    l.records = 1;
    assert!(Collection::import(VCF, Format::Vcf, &l).is_err());
    l = Limits::default();
    l.samples = 1;
    assert!(Collection::import(VCF, Format::Vcf, &l).is_err());
    assert!(Collection::import(b"CRAM\x03\x00", Format::Bcf, &Limits::default()).is_err());
    assert!(Collection::import(VCF, Format::Bcf, &Limits::default()).is_err());
}
#[test]
fn bcf_lengths_and_bgzf_declared_size_crc_and_truncation_are_checked() {
    let c = Collection::import(VCF, Format::Vcf, &Limits::default()).unwrap();
    let mut b = c.export(ExportFormat::Bcf, 1 << 20).unwrap().bytes;
    let n = b.len();
    b[n - 4..].copy_from_slice(&65537u32.to_le_bytes());
    assert!(Collection::import(&b, Format::Bcf, &Limits::default()).is_err());
    let mut b = c.export(ExportFormat::Bcf, 1 << 20).unwrap().bytes;
    b[18] ^= 128;
    assert!(Collection::import(&b, Format::Bcf, &Limits::default()).is_err());
    for raw in [
        b"BCF\x02\x02\xff\xff\xff\xff".as_slice(),
        b"\x1f\x8b\x08\x04".as_slice(),
    ] {
        assert!(Collection::import(raw, Format::Bcf, &Limits::default()).is_err());
    }
}

#[test]
fn independently_framed_empty_bcf_dictionary_is_complete_evidence() {
    // BCF2.2 magic, little-endian header byte count including NUL, VCF header.
    // No writer from either this crate or Noodles constructs this fixture.
    let header = b"##fileformat=VCFv4.3\n##contig=<ID=chrLiteral,length=19>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\0";
    let mut bytes = b"BCF\x02\x02".to_vec();
    bytes.extend_from_slice(&(header.len() as u32).to_le_bytes());
    bytes.extend_from_slice(header);
    let c = Collection::import(&bytes, Format::Bcf, &Limits::default()).unwrap();
    assert!(c.is_empty());
    assert!(c.header().contigs().contains_key("chrLiteral"));
    assert_eq!(
        c.export(ExportFormat::Original, bytes.len()).unwrap().bytes,
        bytes
    );
    bytes.push(0); // Partial next record length cannot masquerade as clean EOF.
    assert!(Collection::import(&bytes, Format::Bcf, &Limits::default()).is_err());
}

#[test]
fn dictionary_idx_resource_preflight_ignores_quoted_description_text() {
    for (kind, fields) in [
        ("INFO", "ID=X,Number=1,Type=Integer,Description=\"x\""),
        ("FORMAT", "ID=X,Number=1,Type=Integer,Description=\"x\""),
        ("FILTER", "ID=X,Description=\"x\""),
        ("contig", "ID=chr,length=9"),
    ] {
        for index in ["4096", "4294967295", "18446744073709551615"] {
            let input = format!(
                "##fileformat=VCFv4.3\n##{kind}=<{fields},IDX={index}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            );
            let err = Collection::import(input.as_bytes(), Format::Vcf, &Limits::default())
                .err()
                .unwrap();
            assert!(err.to_string().contains("IDX"), "{err}");
            // Exact same malicious header in literal BCF framing must fail
            // before Noodles creates its header/string maps.
            let mut bcf = b"BCF\x02\x02".to_vec();
            bcf.extend_from_slice(&((input.len() + 1) as u32).to_le_bytes());
            bcf.extend_from_slice(input.as_bytes());
            bcf.push(0);
            assert!(Collection::import(&bcf, Format::Bcf, &Limits::default()).is_err());
        }
    }
    let benign = b"##fileformat=VCFv4.3\n##FILTER=<ID=X,Description=\"quoted,IDX=4294967295 and \\\"IDX=999999\\\"\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let c = Collection::import(benign, Format::Vcf, &Limits::default()).unwrap();
    assert_eq!(c.original_bytes(), benign);
    assert!(c.export(ExportFormat::Bcf, 1 << 20).is_ok());
}

#[test]
fn literal_extended_bcf_lengths_and_zero_alleles_are_rejected_before_raw_indexing() {
    let header=b"##fileformat=VCFv4.3\n##contig=<ID=chr,length=9>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\0";
    for zero_alleles in [false, true] {
        let mut shared = vec![0u8; 24];
        shared[8..12].copy_from_slice(&1i32.to_le_bytes());
        shared[12..16].copy_from_slice(&0x7f800001u32.to_le_bytes());
        if !zero_alleles {
            shared[18..20].copy_from_slice(&1u16.to_le_bytes());
        }
        // Extended typed string count 2^31-1, with no actual string payload.
        shared.extend_from_slice(&[0xf7, 0x13, 0xff, 0xff, 0xff, 0x7f]);
        let mut bytes = b"BCF\x02\x02".to_vec();
        bytes.extend_from_slice(&(header.len() as u32).to_le_bytes());
        bytes.extend_from_slice(header);
        bytes.extend_from_slice(&(shared.len() as u32).to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&shared);
        assert!(Collection::import(&bytes, Format::Bcf, &Limits::default()).is_err());
    }
}

#[test]
fn wholly_missing_info_conversion_reports_codec_limitation_without_dropping_field() {
    let raw=b"##fileformat=VCFv4.3\n##contig=<ID=chr,length=9>\n##INFO=<ID=X,Number=1,Type=Integer,Description=\"x\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr\t1\t.\tA\tC\t.\t.\tX=.\n";
    let c = Collection::import(raw, Format::Vcf, &Limits::default()).unwrap();
    assert!(
        c.export(ExportFormat::Vcf, 1 << 20)
            .unwrap()
            .bytes
            .ends_with(b"X=.\n")
    );
    assert!(
        c.export(ExportFormat::Bcf, 1 << 20)
            .unwrap_err()
            .to_string()
            .contains("missing INFO")
    );
    assert_eq!(
        c.export(ExportFormat::Original, raw.len()).unwrap().bytes,
        raw
    );
}

#[test]
fn bcf_missing_gt_is_not_a_called_allele_and_codec_limits_fail_atomically() {
    let header = "##fileformat=VCFv4.3\n##contig=<ID=chr,length=9>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n";
    // '.' is one missing GT field, './.' is explicitly two missing alleles;
    // neither permits inventing a call or duplicating ploidy from another sample.
    for (a, b) in [("0|2:9", ".:."), (".:.", "./.:."), ("1:7", "0/2:.")] {
        let raw = format!("{header}chr\t1\t.\tA\tC,G\t.\t.\t.\tGT:DP\t{a}\t{b}\n");
        let original = Collection::import(raw.as_bytes(), Format::Vcf, &Limits::default()).unwrap();
        let expected = original.page(0, 100).unwrap();
        let encoded = original.export(ExportFormat::Bcf, 1 << 20).unwrap();
        let reopened = Collection::import(&encoded.bytes, Format::Bcf, &Limits::default()).unwrap();
        assert_eq!(reopened.page(0, 100).unwrap(), expected);
        assert_eq!(original.original_bytes(), raw.as_bytes());
    }
    // Both are valid source evidence, but the pinned BCF codec cannot safely
    // encode these GT shapes. Exact source and normalized VCF remain available.
    for (alt, a, b) in [
        ("C".to_owned(), "0/1", "0/0/1"),
        (vec!["C"; 63].join(","), "63", "."),
    ] {
        let raw = format!("{header}chr\t1\t.\tA\t{alt}\t.\t.\t.\tGT\t{a}\t{b}\n");
        let source = Collection::import(raw.as_bytes(), Format::Vcf, &Limits::default()).unwrap();
        assert!(source.export(ExportFormat::Bcf, 1 << 20).is_err());
        assert!(source.export(ExportFormat::Vcf, 1 << 20).is_ok());
        assert_eq!(source.original_bytes(), raw.as_bytes());
    }
}
