use collagene_core::{
    Alphabet, Molecule, Sequence, Topology,
    antibody_numbering::{self as n, domain_candidates as d, domain_tracks as t},
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
};
fn inputs() -> (Molecule, d::Options, d::Limits) {
    let m=Molecule {name:"Conditional track fixture".into(),sequence:Sequence::new(Alphabet::Protein,"EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDYGMDVWGQGTLVTVSS").unwrap(),topology:Topology::Linear,features:vec![],metadata:Default::default(),original_genbank:None};
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
    (m, o, l)
}
#[test]
fn worker_sessions_use_one_current_source_options_bound_binary_analysis() {
    let (m, o, l) = inputs();
    let revision = t::request_revision(&m, &o, &l, &selection()).unwrap();
    let q = Request {
        id: 37,
        revision: revision.clone(),
        job: AnalysisJob::AnalyzeAntibody {
            molecule: m.clone(),
            options: o.clone(),
            limits: l.clone(),
            selection: selection(),
        },
    };
    let wire = protocol::encode_request(&q).unwrap();
    assert_eq!(&wire[4..6], &protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    assert_eq!(tag, 195);
    let mut session = protocol::AnalysisSession::default();
    for bytes in [
        protocol::process_request(&wire).unwrap(),
        session.process_request(&wire).unwrap(),
    ] {
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(tag, 185);
        let r = protocol::decode_response(&bytes).unwrap();
        assert_eq!(r.revision, revision);
        let AnalysisResult::Antibody(r) = r.result else {
            panic!("antibody response")
        };
        assert_eq!(r.tracks.len(), 1);
        r.validate_binding(&revision).unwrap();
        r.validate_for(&m, &o, &l, &selection()).unwrap();
    }
    let mut old = wire.clone();
    old[4..6].copy_from_slice(&117u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    assert!(session.process_request(&old).is_err());
    let mut stale = q;
    stale.revision = "f".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    let mut cap = l;
    cap.max_candidates = 0;
    assert_ne!(
        revision,
        t::request_revision(&m, &o, &cap, &selection()).unwrap()
    );
}
#[test]
fn rewritten_packet_checksums_cannot_hide_inconsistent_tracks_or_context() {
    let (m, o, l) = inputs();
    let r = t::analyze(&m, &o, &l, &selection(), &mut || false).unwrap();
    let refuse = |value: t::Report| {
        let body = postcard::to_allocvec(&Response {
            id: 37,
            revision: r.context_revision.clone(),
            result: AnalysisResult::Antibody(value),
        })
        .unwrap();
        let mut wire = b"CGR\0".to_vec();
        wire.extend_from_slice(&protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
        wire.extend_from_slice(&(body.len() as u32).to_le_bytes());
        wire.extend_from_slice(blake3::hash(&body).as_bytes());
        wire.extend_from_slice(&body);
        assert!(protocol::decode_response(&wire).is_err());
    };
    let mut bad = r.clone();
    bad.tracks[0].residues[0].assignment = None;
    refuse(bad);
    let mut bad = r.clone();
    bad.tracks[0].residues[0].residue = 'W';
    refuse(bad);
    let mut bad = r.clone();
    bad.tracks[0].unassigned_residues += 1;
    refuse(bad);
    let mut bad = r.clone();
    bad.tracks[0].candidate_id.push('x');
    refuse(bad);
    let mut bad = r.clone();
    bad.evidence.options.max_fr4_mismatches = 0;
    refuse(bad);
    let mut bad = r.clone();
    bad.evidence.v_evidence.corpus_revision = "f".repeat(64);
    refuse(bad);
    let mut bad = r.clone();
    bad.tracks[0].source_end = usize::MAX;
    refuse(bad);
    let mut bad = r.clone();
    bad.evidence.candidates.clear();
    refuse(bad);
    let mut bad = r.clone();
    bad.evidence.pair_retention_truncated = true;
    refuse(bad);
    let mut bad = r.clone();
    bad.evidence.qualified_pairs_in_retained_v = 0;
    refuse(bad);
    let mut bad = r.clone();
    bad.evidence.v_retention_truncated = true;
    refuse(bad);
    let mut bad = r.clone();
    bad.engine = "antibody-domain-tracks-1".into();
    refuse(bad);
}

fn selection() -> collagene_core::antibody_numbering::domain_tracks::Selection {
    collagene_core::antibody_numbering::domain_tracks::Selection {
        numbering_scheme: collagene_core::antibody_numbering::domain_tracks::NumberingScheme::Imgt,
        region_definition:
            collagene_core::antibody_numbering::domain_tracks::RegionDefinition::Imgt,
    }
}

#[test]
fn scheme_and_definition_choices_bind_worker_requests_and_preserve_independent_regions() {
    let (m, o, l) = inputs();
    let mut c = selection();
    c.numbering_scheme = t::NumberingScheme::Kabat1991;
    let q = Request {
        id: 42,
        revision: t::request_revision(&m, &o, &l, &c).unwrap(),
        job: AnalysisJob::AnalyzeAntibody {
            molecule: m.clone(),
            options: o.clone(),
            limits: l.clone(),
            selection: c.clone(),
        },
    };
    let reply = protocol::decode_response(
        &protocol::process_request(&protocol::encode_request(&q).unwrap()).unwrap(),
    )
    .unwrap();
    let AnalysisResult::Antibody(r) = reply.result else {
        panic!("numbered response")
    };
    assert_eq!(r.selection, c);
    assert_eq!(r.tracks[0].unassigned_residues, 8);
    assert_eq!(r.tracks[0].unassigned_regions, 0);
    r.validate_for(&m, &o, &l, &c).unwrap();
    let mut stale = q;
    let AnalysisJob::AnalyzeAntibody { selection, .. } = &mut stale.job else {
        unreachable!()
    };
    selection.region_definition = t::RegionDefinition::Kabat1991;
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
}

#[test]
fn canonical_region_selection_survives_worker_transport_with_independent_imgt_labels() {
    let (m, o, l) = inputs();
    let mut c = selection();
    c.region_definition = t::RegionDefinition::ChothiaCanonical;
    let q = Request {
        id: 43,
        revision: t::request_revision(&m, &o, &l, &c).unwrap(),
        job: AnalysisJob::AnalyzeAntibody {
            molecule: m.clone(),
            options: o.clone(),
            limits: l.clone(),
            selection: c.clone(),
        },
    };
    let wire = protocol::encode_request(&q).unwrap();
    let mut session = protocol::AnalysisSession::default();
    let reply = protocol::decode_response(&session.process_request(&wire).unwrap()).unwrap();
    let AnalysisResult::Antibody(r) = reply.result else {
        panic!("antibody response")
    };
    assert_eq!(r.selection, c);
    assert_eq!(
        r.tracks[0].cdrs[0].status,
        t::CdrStatus::Mapped { start: 25, end: 32 }
    );
    assert_eq!(r.tracks[0].unassigned_residues, 0);
    assert_eq!(r.tracks[0].unassigned_regions, 8);
    r.validate_for(&m, &o, &l, &c).unwrap();
    let mut old = wire;
    old[4..6].copy_from_slice(&120u16.to_le_bytes());
    assert!(session.process_request(&old).is_err());
}
