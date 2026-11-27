use collagene_core::{
    Alphabet, alignment_antibody as a,
    alignment_formats::{self as f, AlignmentDocument, AlignmentFormat},
    antibody_numbering::{self as n, domain_candidates as d},
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
};
const V: &str = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC";
const J: &str = "WGQGTLVTVSS";
fn inputs() -> (AlignmentDocument, a::Source, d::Options, d::Limits) {
    let seq = format!("--{}---{}ARDYGMDV--{J}--", &V[..25], &V[25..]);
    let doc = f::import_alignment(
        format!(">a\n{seq}\n>b\n{seq}\n").as_bytes(),
        AlignmentFormat::AlignedFasta,
        Alphabet::Protein,
        "Authored row fixture",
    )
    .unwrap();
    let pin = a::source(&doc, &doc.rows[0].id).unwrap();
    let o = d::Options {
        v_scoring: n::v_mapping::Scoring {
            observed_reward: 5,
            frequency_bonus: 8,
            unseen_penalty: -6,
            gap_open: -12,
            gap_extend: -2,
            minimum_score: 300,
        },
        min_cdr3_residues: 2,
        max_cdr3_residues: 60,
        max_fr4_mismatches: 1,
    };
    let l = d::Limits {
        v_mapping: n::v_mapping::Limits {
            max_cells: n::profile::MAX_CELLS,
            max_endpoints_per_locus: n::profile::MAX_ENDPOINTS,
            max_trace_states: n::profile::MAX_TRACE_STATES,
        },
        max_scan_windows: d::MAX_SCAN_WINDOWS,
        max_candidates: 4,
    };
    (doc, pin, o, l)
}
#[test]
fn literal_gap_projection_preserves_source_labels_and_equal_rows_never_retarget() {
    let (doc, pin, o, l) = inputs();
    let r = a::analyze(&doc, &pin, &o, &l, &selection(), &mut || false).unwrap();
    r.validate_for(&doc, &pin, &o, &l, &selection()).unwrap();
    assert_eq!(r.analysis.tracks.len(), 1);
    assert_eq!(r.at_column(0, 30).unwrap().source_index, 25);
    assert_eq!(
        r.at_column(0, 30)
            .unwrap()
            .assignment
            .as_ref()
            .unwrap()
            .position
            .to_string(),
        "27"
    );
    assert_eq!(
        r.at_column(0, V.len() + 5)
            .unwrap()
            .assignment
            .as_ref()
            .unwrap()
            .position
            .to_string(),
        "105"
    );
    assert_eq!(
        r.at_column(0, V.len() + 15)
            .unwrap()
            .assignment
            .as_ref()
            .unwrap()
            .position
            .to_string(),
        "118"
    );
    for gap in [
        0,
        1,
        27,
        28,
        29,
        V.len() + 13,
        V.len() + 14,
        doc.column_count() - 1,
    ] {
        assert!(r.at_column(0, gap).is_none());
    }
    let foreign = a::source(&doc, &doc.rows[1].id).unwrap();
    assert_ne!(
        a::request_revision(&doc, &pin, &o, &l, &selection()).unwrap(),
        a::request_revision(&doc, &foreign, &o, &l, &selection()).unwrap()
    );
    assert!(
        r.validate_source(&doc, &foreign, &o, &l, &selection())
            .is_err()
    );
    let mut reordered = doc.clone();
    reordered.rows.swap(0, 1);
    reordered.validate().unwrap();
    assert!(
        r.validate_source(&reordered, &pin, &o, &l, &selection())
            .is_err()
    );
    let mut changed = o.clone();
    changed.max_fr4_mismatches = 0;
    assert!(
        r.validate_source(&doc, &pin, &changed, &l, &selection())
            .is_err()
    );
    let mut shifted = r.clone();
    for col in &mut shifted.residue_columns {
        *col += 1;
    }
    shifted.validate_binding(&r.context_revision).unwrap(); // Shape alone is not source authentication.
    assert!(
        shifted
            .validate_source(&doc, &pin, &o, &l, &selection())
            .is_err()
    );
}
#[test]
fn worker_session_and_current_binary_selectors_keep_native_row_evidence() {
    let (doc, pin, o, l) = inputs();
    let revision = a::request_revision(&doc, &pin, &o, &l, &selection()).unwrap();
    let q = Request {
        id: 37,
        revision: revision.clone(),
        job: AnalysisJob::AnalyzeAlignmentAntibody {
            alignment: doc.clone(),
            source: pin.clone(),
            options: o.clone(),
            limits: l.clone(),
            selection: selection(),
        },
    };
    let bytes = protocol::encode_request(&q).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 196);
    assert_eq!(&bytes[4..6], &protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    let mut session = protocol::AnalysisSession::default();
    for reply in [
        protocol::process_request(&bytes).unwrap(),
        session.process_request(&bytes).unwrap(),
    ] {
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
        assert_eq!(tag, 186);
        let decoded = protocol::decode_response(&reply).unwrap();
        let AnalysisResult::AlignmentAntibody(r) = decoded.result else {
            panic!("alignment antibody")
        };
        r.validate_source(&doc, &pin, &o, &l, &selection()).unwrap();
        assert_eq!(decoded.revision, revision);
    }
    let mut old = bytes;
    old[4..6].copy_from_slice(&118u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    assert!(session.process_request(&old).is_err());
    let mut stale = q;
    stale.revision = "f".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
}
#[test]
fn malformed_projection_and_source_context_fail_even_with_rewritten_packet_checksums() {
    let (doc, pin, o, l) = inputs();
    let r = a::analyze(&doc, &pin, &o, &l, &selection(), &mut || false).unwrap();
    let refuse = |value: a::Report| {
        let body = postcard::to_allocvec(&Response {
            id: 1,
            revision: r.context_revision.clone(),
            result: AnalysisResult::AlignmentAntibody(value),
        })
        .unwrap();
        let mut wire = b"CGR\0".to_vec();
        wire.extend(protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
        wire.extend((body.len() as u32).to_le_bytes());
        wire.extend(blake3::hash(&body).as_bytes());
        wire.extend(body);
        assert!(protocol::decode_response(&wire).is_err());
    };
    let mut bad = r.clone();
    bad.residue_columns[1] = bad.residue_columns[0];
    refuse(bad);
    let mut bad = r.clone();
    bad.source.row_id = doc.rows[1].id.clone();
    refuse(bad);
    let mut bad = r.clone();
    bad.column_count = usize::MAX;
    refuse(bad);
    let mut bad = r.clone();
    bad.engine.push('x');
    refuse(bad);
    let mut bad = r.clone();
    bad.analysis.tracks[0].residues[0].assignment = None;
    refuse(bad);
}
#[test]
fn alphabet_empty_row_source_size_cancellation_and_counts_only_scope_are_explicit() {
    let (doc, pin, o, mut l) = inputs();
    l.max_candidates = 0;
    let r = a::analyze(&doc, &pin, &o, &l, &selection(), &mut || false).unwrap();
    assert!(r.analysis.tracks.is_empty() && r.analysis.evidence.pair_retention_truncated);
    r.validate_source(&doc, &pin, &o, &l, &selection()).unwrap();
    assert!(a::analyze(&doc, &pin, &o, &l, &selection(), &mut || true).is_err());
    let dna = f::import_alignment(
        b">a\nACGT\n",
        AlignmentFormat::AlignedFasta,
        Alphabet::Dna,
        "DNA",
    )
    .unwrap();
    assert!(a::source(&dna, &dna.rows[0].id).is_err());
    let empty = f::import_alignment(
        b">a\n----\n",
        AlignmentFormat::AlignedFasta,
        Alphabet::Protein,
        "Gaps",
    )
    .unwrap();
    assert!(a::source(&empty, &empty.rows[0].id).is_err());
    let big = f::import_alignment(
        format!(">a\n{}\n", "A".repeat(2001)).as_bytes(),
        AlignmentFormat::AlignedFasta,
        Alphabet::Protein,
        "Large",
    )
    .unwrap();
    assert!(a::source(&big, &big.rows[0].id).is_err());
    assert!(a::source(&doc, "unknown").is_err());
}

fn selection() -> collagene_core::antibody_numbering::domain_tracks::Selection {
    collagene_core::antibody_numbering::domain_tracks::Selection {
        numbering_scheme: collagene_core::antibody_numbering::domain_tracks::NumberingScheme::Imgt,
        region_definition:
            collagene_core::antibody_numbering::domain_tracks::RegionDefinition::Imgt,
    }
}
