//! Format/source-integrity fixtures independent of genomic analysis.
use collagene_core::{genome::*, *};
fn build(raw: &[u8], info: BuildInfo, step: usize) -> (Vec<EncodedBlock>, FinishedReference) {
    let mut builder = FastaBuilder::begin(info, IngestOptions::default()).unwrap();
    let mut blocks = vec![];
    for bytes in raw.chunks(step) {
        blocks.extend(builder.push(bytes).unwrap());
    }
    let (last, finished) = builder.finish().unwrap();
    blocks.extend(last);
    (blocks, finished)
}
#[test]
fn iupac_lowercase_split_headers_and_odd_padding_round_trip() {
    let raw = b">alpha header\r\nacgtryswkmbdhvn\r\n>a2\nA";
    let (blocks, f) = build(raw, BuildInfo::default(), 1);
    let verifier = BlockVerifier::begin(&f.manifest).unwrap();
    assert_eq!(
        verifier
            .verify("alpha", 0, &blocks[0].encoded)
            .unwrap()
            .dna(),
        b"ACGTRYSWKMBDHVN"
    );
    assert_eq!(blocks[1].encoded, vec![0x10]);
    assert_eq!(
        f.manifest.raw_source_hash,
        blake3::hash(raw).to_hex().to_string()
    );
    let mut semantic = SemanticValidator::begin(&f.manifest).unwrap();
    for b in &blocks {
        semantic
            .push(
                &verifier
                    .verify(&b.contig_id, b.descriptor.index, &b.encoded)
                    .unwrap(),
            )
            .unwrap();
    }
    assert_eq!(semantic.finish().unwrap().total_bases, 16);
    assert!(!f.manifest.expected_complete);
}
#[test]
fn canonical_content_is_format_order_and_label_independent_but_roles_remain_distinct() {
    let (_, a) = build(b">a first\nacgtryn\n>b\nTTCC\n", BuildInfo::default(), 3);
    let mut info = BuildInfo::default();
    info.name = "Different imported label".into();
    let (_, b) = build(b">b different\n ttcc \n>a second\nACGTRY N\n", info, 11);
    assert_eq!(a.manifest.content_revision, b.manifest.content_revision);
    assert_ne!(a.manifest_hash, b.manifest_hash);
    let (_, one) = build(b">a\nAAAA\n", BuildInfo::default(), 4);
    let (_, two) = build(b">a\nAAAA\n>b\nAAAA\n", BuildInfo::default(), 4);
    assert_ne!(one.manifest.content_revision, two.manifest.content_revision);
}
#[test]
fn expected_roles_missing_components_and_contradictory_catalogue_are_explicit() {
    let info = BuildInfo {
        name: "A name is not an attestation".into(),
        assembly_accession: Some("declared.1".into()),
        version: Some("v1".into()),
        provenance: vec![],
        expected_components: vec![
            ExpectedComponent {
                id: "chr1".into(),
                source_id: "NC_fixture.1".into(),
                accession: Some("NC_fixture.1".into()),
                length: 4,
                sequence_hash: Some(blake3::hash(b"ACGT").to_hex().to_string()),
                topology: Topology::Linear,
            },
            ExpectedComponent {
                id: "chr2".into(),
                source_id: "NC_fixture.2".into(),
                accession: Some("NC_fixture.2".into()),
                length: 4,
                sequence_hash: None,
                topology: Topology::Circular,
            },
        ],
    };
    let (_, partial) = build(b">NC_fixture.1\nACGT\n", info.clone(), 3);
    assert_eq!(partial.manifest.contigs[0].id, "chr1");
    assert_eq!(partial.manifest.missing_expected, vec!["chr2"]);
    assert!(!partial.manifest.expected_complete);
    let (_, full) = build(
        b">NC_fixture.1\nACGT\n>NC_fixture.2\nTGCA\n",
        info.clone(),
        7,
    );
    assert!(full.manifest.expected_complete);
    let mut bad = FastaBuilder::begin(info, IngestOptions::default()).unwrap();
    bad.push(b">NC_fixture.1\nAAAA\n").unwrap();
    assert!(bad.finish().is_err());
}
#[test]
fn framing_and_hostile_manifest_counts_are_bounded_before_decode() {
    let (_, f) = build(b">a\nACGTN\n", BuildInfo::default(), 3);
    let bytes = encode_manifest(&f.manifest).unwrap();
    let tail = encode_tail(&f).unwrap();
    let parsed = parse_tail(&HEADER, &tail, f.file_length).unwrap();
    assert_eq!(decode_manifest(&bytes, &parsed).unwrap(), f);
    let mut bad = tail;
    bad[..8].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(parse_tail(&HEADER, &bad, f.file_length).is_err());
    assert!(parse_tail(&[0; 8], &tail, f.file_length).is_err());
    let malformed = vec![1, 255, 255, 255, 255, 255, 255, 255, 255, 255, 1];
    let forged = TailInfo {
        manifest_length: malformed.len() as u64,
        manifest_hash: blake3::hash(&malformed).to_hex().to_string(),
        manifest_offset: 8,
        file_length: 8 + malformed.len() as u64 + 48,
    };
    assert!(decode_manifest(&malformed, &forged).is_err());
}
#[test]
fn missing_duplicate_corrupt_and_reordered_blocks_never_complete_verification() {
    let dna = "ACGT".repeat(BLOCK_BASES / 4 + 3);
    let (blocks, f) = build(format!(">a\n{dna}").as_bytes(), BuildInfo::default(), 8191);
    let verifier = BlockVerifier::begin(&f.manifest).unwrap();
    let mut corrupt = blocks[0].encoded.clone();
    corrupt[0] ^= 1;
    assert!(verifier.verify("a", 0, &corrupt).is_err());
    let mut missing = SemanticValidator::begin(&f.manifest).unwrap();
    missing
        .push(&verifier.verify("a", 0, &blocks[0].encoded).unwrap())
        .unwrap();
    assert!(missing.finish().is_err());
    let mut duplicate = SemanticValidator::begin(&f.manifest).unwrap();
    let first = verifier.verify("a", 0, &blocks[0].encoded).unwrap();
    duplicate.push(&first).unwrap();
    assert!(duplicate.push(&first).is_err());
    assert!(duplicate.finish().is_err());
    let mut reordered = SemanticValidator::begin(&f.manifest).unwrap();
    assert!(
        reordered
            .push(&verifier.verify("a", 1, &blocks[1].encoded).unwrap())
            .is_err()
    );
}
#[test]
fn poisoned_ingest_and_portable_links_do_not_restore_authority() {
    let mut options = IngestOptions::default();
    options.max_push_bytes = 3;
    let mut builder = FastaBuilder::begin(BuildInfo::default(), options).unwrap();
    assert!(builder.push(b">toolong\n").is_err());
    assert!(builder.push(b">a\n").is_err());
    assert!(builder.finish().is_err());
    let (_, f) = build(b">a\nACGT\n", BuildInfo::default(), 4);
    let link = reference_link(&f, "logical-source", Some("local-cache-1".into())).unwrap();
    let (_, attachment) = make_attachment(&link).unwrap();
    assert_eq!(load_link(&attachment.bytes).unwrap(), link);
    let mut json = serde_json::to_value(link).unwrap();
    json["semantic_verified"] = serde_json::json!(true);
    assert!(load_link(&serde_json::to_vec(&json).unwrap()).is_err());
    let position = GlobalPosition {
        content_revision: f.manifest.content_revision.clone(),
        contig_id: "a".into(),
        position: u64::from(u32::MAX) + 7,
    };
    assert_eq!(
        postcard::from_bytes::<GlobalPosition>(&postcard::to_allocvec(&position).unwrap()).unwrap(),
        position
    );
    assert!(position.validate_for(&f.manifest).is_err());
}
