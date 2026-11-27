use collagene_core::{
    Alphabet, Molecule, Sequence, Topology, protein_liabilities as native,
    protein_liability_annotation,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
};
fn inputs() -> (Molecule, native::Options) {
    (
        Molecule {
            name: "Protein fixture".into(),
            sequence: Sequence::new(Alphabet::Protein, "NATMWX").unwrap(),
            topology: Topology::Linear,
            features: vec![],
            metadata: Default::default(),
            original_genbank: None,
        },
        native::Options {
            max_hits: 100,
            provided_cdrs: None,
        },
    )
}
#[test]
fn native_worker_and_retained_session_share_current_binding_and_candidate_annotation_rules() {
    let (m, o) = inputs();
    let q = Request {
        id: 37,
        revision: native::request_revision(&m, &o).unwrap(),
        job: AnalysisJob::ScreenProteinLiabilities {
            molecule: m.clone(),
            options: o.clone(),
        },
    };
    let bytes = protocol::encode_request(&q).unwrap();
    assert_eq!(&bytes[4..6], &protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 194);
    let mut session = protocol::AnalysisSession::default();
    for reply in [
        protocol::process_request(&bytes).unwrap(),
        session.process_request(&bytes).unwrap(),
    ] {
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
        assert_eq!(tag, 184);
        let r = protocol::decode_response(&reply).unwrap();
        let AnalysisResult::ProteinLiabilities(r) = r.result else {
            panic!("liability result")
        };
        r.validate_for(&m, &o).unwrap();
        let h = r
            .hits
            .iter()
            .find(|h| h.rule == native::RuleId::NGlycosylationSequon)
            .unwrap();
        let f = protein_liability_annotation::candidate(&m, &o, h).unwrap();
        assert!(f.needs_review);
        assert_eq!(f.kind, "misc_feature");
        assert_eq!(f.location.spans(), vec![(0, 3, false)]);
        assert!(
            f.qualifiers
                .iter()
                .any(|(k, v)| k == "note" && v.as_ref().unwrap().contains("not a measured PTM"))
        );
        let mut wrong = h.clone();
        wrong.residues = "NPT".into();
        assert!(protein_liability_annotation::candidate(&m, &o, &wrong).is_err());
    }
    let mut old = bytes;
    old[4..6].copy_from_slice(&116u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    assert!(session.process_request(&old).is_err());
    let mut stale = q;
    stale.revision = "f".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
}
#[test]
fn rewritten_transport_checksums_do_not_admit_false_counts_certainty_or_source_scope() {
    let (m, o) = inputs();
    let report = native::screen(&m, &o, &mut || false).unwrap();
    let expected_revision = report.context_revision.clone();
    let refuse = |r: native::Report| {
        let response = Response {
            id: 37,
            revision: expected_revision.clone(),
            result: AnalysisResult::ProteinLiabilities(r),
        };
        let body = postcard::to_allocvec(&response).unwrap();
        let mut out = b"CGR\0".to_vec();
        out.extend_from_slice(&protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
        out.extend_from_slice(&(body.len() as u32).to_le_bytes());
        out.extend_from_slice(blake3::hash(&body).as_bytes());
        out.extend_from_slice(&body);
        assert!(protocol::decode_response(&out).is_err());
    };
    let mut r = report.clone();
    r.counts[0].certain += 1;
    refuse(r);
    let mut r = report.clone();
    r.hits[0].certainty = native::Certainty::PossiblePattern;
    refuse(r);
    let mut r = report.clone();
    r.hits[0].cdr = native::CdrMembership::OutsideProvided;
    refuse(r);
    let mut r = report.clone();
    r.catalog_revision = "f".repeat(64);
    refuse(r);
    let mut r = report.clone();
    r.hits_truncated = true;
    refuse(r);
    let mut r = report;
    r.options.max_hits = 0;
    refuse(r);
}
