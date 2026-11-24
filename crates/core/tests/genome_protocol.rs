//! Bounded genome worker lifecycle, private live authority and global source bindings.
use collagene_core::{
    Error, Topology,
    crispr::{GuideOptions, GuideQuery},
    genome::*,
    genome_crispr::*,
    genome_session::*,
    protocol::{self, AnalysisJob as J, AnalysisResult as R, AnalysisSession, Request, Response},
};
const SPACER: &str = "ACGTACGTACGTACGTACGT";
fn run(s: &mut AnalysisSession, revision: &str, job: J) -> R {
    let b = protocol::encode_request(&Request {
        id: 801,
        revision: revision.into(),
        job,
    })
    .unwrap();
    assert_eq!(&b[4..6], &21u16.to_le_bytes());
    protocol::decode_response(&s.process_request(&b).unwrap())
        .unwrap()
        .result
}
fn reject(s: &mut AnalysisSession, revision: &str, job: J) {
    let b = protocol::encode_request(&Request {
        id: 802,
        revision: revision.into(),
        job,
    })
    .unwrap();
    assert!(s.process_request(&b).is_err());
}
fn fixture() -> (Vec<u8>, BuildInfo) {
    let mut dna = vec![b'A'; BLOCK_BASES + 100];
    let start = BLOCK_BASES - 10;
    dna[start..start + 20].copy_from_slice(SPACER.as_bytes());
    dna[start + 20..start + 23].copy_from_slice(b"AGG");
    let mut bytes = b">chr\r\n".to_vec();
    bytes.extend(&dna);
    bytes.push(b'\n');
    let mut build = BuildInfo::default();
    build.expected_components = vec![ExpectedComponent {
        id: "chr".into(),
        source_id: "chr".into(),
        accession: Some("test.1".into()),
        length: dna.len() as u64,
        sequence_hash: None,
        topology: Topology::Linear,
    }];
    (bytes, build)
}
fn ingest(s: &mut AnalysisSession) -> (FinishedReference, Vec<EncodedBlock>) {
    let (bytes, build) = fixture();
    let options = IngestOptions::default();
    let rev = fasta_revision(&build, &options, bytes.len() as u64).unwrap();
    let R::GenomeImportStarted(start) = run(
        s,
        &rev,
        J::BeginGenomeFasta {
            build,
            options,
            expected_source_bytes: bytes.len() as u64,
        },
    ) else {
        panic!()
    };
    let mut blocks = vec![];
    let mut at = 0;
    for chunk in bytes.chunks(17001) {
        let R::GenomeImportProgress(p) = run(
            s,
            &start.import_id,
            J::PushGenomeFasta {
                import_id: start.import_id.clone(),
                offset: at,
                bytes: chunk.to_vec(),
            },
        ) else {
            panic!()
        };
        at = p.received;
        blocks.extend(p.blocks);
    }
    let R::GenomeFastaFinished(done) = run(
        s,
        &start.import_id,
        J::FinishGenomeFasta {
            import_id: start.import_id.clone(),
        },
    ) else {
        panic!()
    };
    blocks.extend(done.blocks);
    (done.reference, blocks)
}
fn verify(
    s: &mut AnalysisSession,
    f: &FinishedReference,
    blocks: &[EncodedBlock],
) -> RegisteredGenome {
    let manifest = encode_manifest(&f.manifest).unwrap();
    let tail = encode_tail(f).unwrap().to_vec();
    let rev = native_revision(&HEADER, &tail, f.file_length, &manifest).unwrap();
    let R::GenomeImportStarted(start) = run(
        s,
        &rev,
        J::BeginGenomeNative {
            header: HEADER.to_vec(),
            tail,
            file_length: f.file_length,
            manifest,
        },
    ) else {
        panic!()
    };
    reject(
        s,
        &start.import_id,
        J::FinishGenomeNative {
            import_id: start.import_id.clone(),
        },
    );
    if blocks.len() > 1 {
        reject(
            s,
            &start.import_id,
            J::PushGenomeNativeBlock {
                import_id: start.import_id.clone(),
                block: blocks[1].clone(),
            },
        );
    }
    for b in blocks {
        run(
            s,
            &start.import_id,
            J::PushGenomeNativeBlock {
                import_id: start.import_id.clone(),
                block: b.clone(),
            },
        );
    }
    let R::GenomeRegistered(r) = run(
        s,
        &start.import_id,
        J::FinishGenomeNative {
            import_id: start.import_id.clone(),
        },
    ) else {
        panic!()
    };
    assert!(r.semantic_verified);
    assert!(r.receipt.is_some());
    r
}
fn requested(blocks: &[EncodedBlock], requests: &[BlockRequest]) -> Vec<EncodedBlock> {
    requests
        .iter()
        .map(|r| {
            blocks
                .iter()
                .find(|b| b.contig_id == r.contig_id && b.descriptor.index == r.index)
                .unwrap()
                .clone()
        })
        .collect()
}
#[test]
fn private_readback_global_windows_and_whole_manifest_search_require_exact_ordered_blocks() {
    let mut s = AnalysisSession::default();
    let (f, blocks) = ingest(&mut s);
    let verified = verify(&mut s, &f, &blocks);
    let R::GenomeRegistered(registered) = run(
        &mut s,
        &f.manifest_hash,
        J::RegisterGenome {
            manifest: f.manifest.clone(),
        },
    ) else {
        panic!()
    };
    assert!(!registered.semantic_verified);
    assert!(registered.receipt.is_none());
    assert_ne!(registered.genome.id, verified.genome.id);
    let g = registered.genome;
    let start = BLOCK_BASES as u64 - 10;
    let options = GenomeGuideOptions {
        source_role: "reference".into(),
        interval: GenomicInterval {
            contig_id: "chr".into(),
            start,
            end: start + 1,
        },
        guide: GuideOptions::default(),
    };
    let req = interval_required_blocks(&f.manifest, &options).unwrap();
    let R::GenomeIntervalGuides(guides) = run(
        &mut s,
        &g.manifest_hash,
        J::GenomeIntervalGuides {
            genome: g.clone(),
            options: options.clone(),
            blocks: requested(&blocks, &req),
        },
    ) else {
        panic!()
    };
    guides
        .report
        .validate_bindings_for(&f.manifest, &options)
        .unwrap();
    let candidate = guides
        .report
        .candidates
        .iter()
        .find(|c| c.guide_sequence == SPACER && !c.reverse)
        .unwrap();
    assert!(candidate.guide_positions.contains(&(BLOCK_BASES as u64)));
    let options = GenomeIntervalReadOptions {
        source_role: "reference".into(),
        interval: GenomicInterval {
            contig_id: "chr".into(),
            start,
            end: start + 23,
        },
    };
    let req = read_interval_required_blocks(&f.manifest, &options).unwrap();
    let R::GenomeIntervalRead(read) = run(
        &mut s,
        &g.manifest_hash,
        J::ReadGenomeInterval {
            genome: g.clone(),
            options: options.clone(),
            blocks: requested(&blocks, &req),
        },
    ) else {
        panic!()
    };
    assert_eq!(read.report.sequence, format!("{SPACER}AGG"));
    read.report
        .validate_bindings_for(&f.manifest, &options)
        .unwrap();
    let query = GuideQuery {
        id: "query".into(),
        sequence: SPACER.into(),
    };
    let mut options = GenomeSearchOptions::default();
    options.max_mismatches = 0;
    options.intended = Some(GenomeIntended {
        reference_role: "reference".into(),
        manifest_hash: f.manifest_hash.clone(),
        content_revision: f.manifest.content_revision.clone(),
        contig_id: "chr".into(),
        guide_positions: candidate.guide_positions.clone(),
        pam_positions: candidate.pam_positions.clone(),
        reverse: false,
    });
    let R::GenomeSearchState(mut state) = run(
        &mut s,
        &g.manifest_hash,
        J::StartGenomeSearch {
            genome: g.clone(),
            query: query.clone(),
            options: options.clone(),
        },
    ) else {
        panic!()
    };
    assert!(!state.progress.source_integrity_verified);
    assert!(!state.progress.complete);
    while !state.progress.finished {
        let data = requested(&blocks, &state.progress.next_blocks);
        let R::GenomeSearchState(next) = run(
            &mut s,
            &state.job_id,
            J::AdvanceGenomeSearch {
                job_id: state.job_id.clone(),
                blocks: data,
            },
        ) else {
            panic!()
        };
        state = next;
    }
    state
        .progress
        .validate_bindings_for(&f.manifest, &query, &options)
        .unwrap();
    assert!(
        state.progress.source_integrity_verified
            && state.progress.exhaustive_for_manifest
            && state.progress.complete
    );
    let R::GenomeSearchPage(page) = run(
        &mut s,
        &state.job_id,
        J::PageGenomeSearch {
            job_id: state.job_id.clone(),
            offset: 0,
            limit: 100,
        },
    ) else {
        panic!()
    };
    assert_eq!(page.hits.len(), 1);
    assert!(page.hits[0].is_intended);
    run(
        &mut s,
        &state.job_id,
        J::CancelGenomeJob {
            id: state.job_id.clone(),
        },
    );
    reject(
        &mut s,
        &state.job_id,
        J::PageGenomeSearch {
            job_id: state.job_id.clone(),
            offset: 0,
            limit: 100,
        },
    );
    let mut other = AnalysisSession::default();
    reject(
        &mut other,
        &g.manifest_hash,
        J::ReadGenomeInterval {
            genome: g.clone(),
            options: options_for_read(),
            blocks: blocks.clone(),
        },
    );
}
fn options_for_read() -> GenomeIntervalReadOptions {
    GenomeIntervalReadOptions {
        source_role: "reference".into(),
        interval: GenomicInterval {
            contig_id: "chr".into(),
            start: 0,
            end: 5,
        },
    }
}
#[test]
fn declared_source_length_and_failed_semantic_hashes_never_publish_partial_authority() {
    let mut s = AnalysisSession::default();
    let bytes = b">tiny\nACGT\n";
    let build = BuildInfo::default();
    let options = IngestOptions::default();
    let rev = fasta_revision(&build, &options, bytes.len() as u64).unwrap();
    assert_ne!(
        rev,
        fasta_revision(&build, &options, bytes.len() as u64 + 1).unwrap()
    );
    let R::GenomeImportStarted(start) = run(
        &mut s,
        &rev,
        J::BeginGenomeFasta {
            build,
            options,
            expected_source_bytes: bytes.len() as u64,
        },
    ) else {
        panic!()
    };
    let n = bytes.len() - 2;
    run(
        &mut s,
        &start.import_id,
        J::PushGenomeFasta {
            import_id: start.import_id.clone(),
            offset: 0,
            bytes: bytes[..n].to_vec(),
        },
    );
    reject(
        &mut s,
        &start.import_id,
        J::FinishGenomeFasta {
            import_id: start.import_id.clone(),
        },
    );
    reject(
        &mut s,
        &start.import_id,
        J::PushGenomeFasta {
            import_id: start.import_id.clone(),
            offset: n as u64 + 1,
            bytes: bytes[n..].to_vec(),
        },
    );
    run(
        &mut s,
        &start.import_id,
        J::PushGenomeFasta {
            import_id: start.import_id.clone(),
            offset: n as u64,
            bytes: bytes[n..].to_vec(),
        },
    );
    let R::GenomeFastaFinished(f) = run(
        &mut s,
        &start.import_id,
        J::FinishGenomeFasta {
            import_id: start.import_id.clone(),
        },
    ) else {
        panic!()
    };
    assert_eq!(f.reference.manifest.source_bytes, bytes.len() as u64);
    let (f, mut blocks) = ingest(&mut s);
    let R::GenomeRegistered(r) = run(
        &mut s,
        &f.manifest_hash,
        J::RegisterGenome {
            manifest: f.manifest.clone(),
        },
    ) else {
        panic!()
    };
    let query = GuideQuery {
        id: "q".into(),
        sequence: SPACER.into(),
    };
    let R::GenomeSearchState(state) = run(
        &mut s,
        &f.manifest_hash,
        J::StartGenomeSearch {
            genome: r.genome,
            query,
            options: GenomeSearchOptions::default(),
        },
    ) else {
        panic!()
    };
    blocks[0].encoded[0] ^= 1;
    reject(
        &mut s,
        &state.job_id,
        J::AdvanceGenomeSearch {
            job_id: state.job_id.clone(),
            blocks: requested(&blocks, &state.progress.next_blocks),
        },
    );
    let R::GenomeSearchPage(page) = run(
        &mut s,
        &state.job_id,
        J::PageGenomeSearch {
            job_id: state.job_id.clone(),
            offset: 0,
            limit: 1,
        },
    ) else {
        panic!()
    };
    assert!(page.progress.finished);
    assert!(
        !page.progress.complete
            && !page.progress.exhaustive_for_manifest
            && !page.progress.source_integrity_verified
    );
}
fn raw(magic: &[u8; 4], version: u16, p: &[u8]) -> Vec<u8> {
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(p).as_bytes());
    b.extend(p);
    b
}
#[test]
fn genome_tags_require_schema_nineteen_and_never_expand_the_ordinary_message_budget() {
    for (magic, start, end) in [(b"CGQ\0", 53u32, 66u32), (b"CGR\0", 51u32, 59u32)] {
        for tag in start..=end {
            let prefix = postcard::to_allocvec(&(803u32, "session", tag)).unwrap();
            for version in 1..19 {
                let b = raw(magic, version, &prefix);
                let e = if magic == b"CGQ\0" {
                    AnalysisSession::default().process_request(&b).unwrap_err()
                } else {
                    protocol::decode_response(&b).unwrap_err()
                };
                assert!(e.to_string().contains("schema 19"), "{e}");
            }
        }
        let mut p = postcard::to_allocvec(&(803u32, "session", start)).unwrap();
        p.resize(collagene_core::formats::MAX_IMPORT_BYTES + 1, 0);
        let b = raw(magic, 19, &p);
        let e = if magic == b"CGQ\0" {
            AnalysisSession::default().process_request(&b).unwrap_err()
        } else {
            protocol::decode_response(&b).unwrap_err()
        };
        assert!(matches!(e, Error::TooLarge));
    }
}
#[test]
fn reference_links_are_typed_bounded_producer_links_and_require_worker_nineteen() {
    use collagene_core::workspace::Workspace;
    let mut s = AnalysisSession::default();
    let (f, _) = ingest(&mut s);
    let link = reference_link(&f, "reference", Some("opaque-browser-cache-key".into())).unwrap();
    let (key, attachment) = make_attachment(&link).unwrap();
    let molecule = collagene_core::Molecule::dna("source", "ACGT", Topology::Linear).unwrap();
    let mut workspace = Workspace::from_documents("links", vec![molecule.clone()], 0).unwrap();
    workspace.attachments.insert(key, attachment);
    workspace.validate().unwrap();
    let req = Request {
        id: 804,
        revision: molecule.revision(),
        job: J::SaveWorkspace {
            workspace: workspace.clone(),
        },
    };
    let payload = postcard::to_allocvec(&req).unwrap();
    assert!(
        protocol::process_request(&raw(b"CGQ\0", 18, &payload))
            .unwrap_err()
            .to_string()
            .contains("schema 19")
    );
    let response = Response {
        id: 804,
        revision: molecule.revision(),
        result: R::WorkspaceLoaded(workspace),
    };
    assert!(
        protocol::decode_response(&raw(
            b"CGR\0",
            18,
            &postcard::to_allocvec(&response).unwrap()
        ))
        .unwrap_err()
        .to_string()
        .contains("schema 19")
    );
}
