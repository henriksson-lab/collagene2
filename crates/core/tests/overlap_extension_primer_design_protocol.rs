//! Advisory auto-primer proposals bind original DNA and independent edit intent.
use collagene_core::{
    Molecule, Topology,
    assembly_edit_intent::Edit,
    overlap_extension_primer_design::{self, Options, Request as DesignRequest},
    overlap_primer_design::Fragment,
    primers::DesignOptions,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    transfer::{self, Candidate},
    verification::{self, Status},
};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const EXPECTED: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
fn fixture() -> (Molecule, DesignRequest) {
    (
        Molecule::dna("Immutable original", ORIGINAL, Topology::Linear).unwrap(),
        DesignRequest {
            fragment: Fragment {
                role: "insert".into(),
                input_id: "original".into(),
                start: 0,
                end: 40,
                reverse: false,
            },
            edits: vec![Edit {
                role: "insert".into(),
                start: 16,
                end: 20,
                expected: "GGAC".into(),
                replacement: "TT".into(),
            }],
            options: Options {
                primer_options: DesignOptions {
                    min_length: 8,
                    max_length: 8,
                    min_gc_fraction: 0.,
                    max_gc_fraction: 1.,
                    max_tm_difference_c: 100.,
                    max_homopolymer: 200,
                    max_three_prime_complementarity: 200,
                    ..Default::default()
                },
                homology_flank: 8,
                min_paired_bases: 8,
                ..Default::default()
            },
        },
    )
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut out = magic.to_vec();
    out.extend(version.to_le_bytes());
    out.extend((payload.len() as u32).to_le_bytes());
    out.extend(blake3::hash(&payload).as_bytes());
    out.extend(payload);
    out
}
fn message(source: Molecule, request: DesignRequest, revision: String) -> Request {
    Request {
        id: 164,
        revision,
        job: AnalysisJob::DesignOverlapExtensionPrimers { source, request },
    }
}
#[test]
fn appended164_156_proposal_roundtrip_keeps_independent_expected_and_original_stock() {
    let (source, request) = fixture();
    let before = serde_json::to_value(&source).unwrap();
    let context = overlap_extension_primer_design::design_revision(&source, &request).unwrap();
    assert_ne!(context, source.revision());
    let r = message(source.clone(), request, context.clone());
    let bytes = protocol::encode_request(&r).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 164);
    assert!(protocol::process_request(&wire(b"CGQ\0", 82, &r)).is_err());
    let reply = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
    assert_eq!(tag, 156);
    let response = protocol::decode_response(&reply).unwrap();
    assert_eq!(response.revision, context);
    assert!(protocol::decode_response(&wire(b"CGR\0", 82, &response)).is_err());
    let AnalysisResult::OverlapExtensionPrimersDesigned(report) = response.result else {
        panic!("primer proposal");
    };
    report.validate_binding(&context).unwrap();
    assert_eq!(report.source_revision, source.revision());
    assert_eq!(report.expected.sequence.as_str(), EXPECTED);
    assert_eq!(report.designs.len(), 3);
    assert_eq!(
        report.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.verification.findings
    );
    assert_eq!(
        report
            .verification
            .products
            .iter()
            .map(|m| m.sequence.len())
            .collect::<Vec<_>>(),
        vec![26, 30, 38, 38]
    );
    assert_eq!(report.plan.inputs[0].molecule.sequence.as_str(), ORIGINAL);
    assert_eq!(report.plan.inputs[0].revision, source.revision());
    assert_eq!(
        report.preparation.left.reverse.ordered_sequence(),
        "GGTCTGAAAATACGTAGC"
    );
    assert_eq!(
        report.preparation.right.forward.ordered_sequence(),
        "GCTACGTATTTTCAGACC"
    );
    let candidate = Candidate::new(report.plan.clone());
    assert_eq!(protocol::CANDIDATE_SCHEMA_VERSION, 18);
    assert!(transfer::load_candidate(&transfer::save_candidate(&candidate).unwrap()).is_ok());
    let verified = verification::verify(&candidate.plan, &report.expected);
    assert_eq!(verified.status, Status::ConsistentUnderAssumptions);
    assert_eq!(
        verified.products.last().unwrap().sequence.as_str(),
        EXPECTED
    );
    assert_eq!(serde_json::to_value(&source).unwrap(), before);
}
#[test]
fn outer_context_and_report_binding_include_source_metadata_edit_literals_and_search_options() {
    let (source, request) = fixture();
    let context = overlap_extension_primer_design::design_revision(&source, &request).unwrap();
    assert!(
        protocol::process_request(
            &protocol::encode_request(&message(source.clone(), request.clone(), source.revision()))
                .unwrap()
        )
        .is_err()
    );
    let mut changed_source = source.clone();
    changed_source.name.push_str(" changed source event");
    let mut changed_options = request.clone();
    changed_options.options.max_candidates = 99;
    let mut changed_edit = request.clone();
    changed_edit.edits[0].replacement = "TA".into();
    for (source, request) in [
        (changed_source, request.clone()),
        (source.clone(), changed_options),
        (source.clone(), changed_edit),
    ] {
        assert!(
            protocol::process_request(
                &protocol::encode_request(&message(source, request, context.clone())).unwrap()
            )
            .is_err()
        );
    }
    let report = overlap_extension_primer_design::design(&source, &request).unwrap();
    let mut changed = report.clone();
    changed.preparation.extension.overlap -= 1;
    assert!(changed.validate_binding(&context).is_err());
    let mut changed = report;
    changed.source_revision = "0".repeat(64);
    assert!(changed.validate_binding(&context).is_err());
    let mut invalid = request;
    invalid.edits[0].expected = "AAAA".into();
    assert!(overlap_extension_primer_design::design(&source, &invalid).is_err());
    let mut oversized = invalid;
    oversized.edits[0].expected = "GGAC".into();
    oversized.options.max_candidates = usize::MAX;
    assert!(protocol::encode_request(&message(source, oversized, context)).is_err());
}
#[test]
fn multiple_edits_and_reverse_original_stock_remain_explicit_in_binary_proposals() {
    let (source, mut request) = fixture();
    request.edits.push(Edit {
        role: "insert".into(),
        start: 23,
        end: 24,
        expected: "A".into(),
        replacement: "G".into(),
    });
    let expected = "AGTCCGATGCTACGTATTTTCGGACCTGTAGATCGTCA";
    let reverse_source = Molecule::dna(
        "Reverse original stock",
        "TGACGATCTACAGGTCTGAAGTCCTACGTAGCATCGGACT",
        Topology::Linear,
    )
    .unwrap();
    for (source, reverse) in [(source, false), (reverse_source, true)] {
        request.fragment.reverse = reverse;
        let revision = overlap_extension_primer_design::design_revision(&source, &request).unwrap();
        let response = protocol::decode_response(
            &protocol::process_request(
                &protocol::encode_request(&message(
                    source.clone(),
                    request.clone(),
                    revision.clone(),
                ))
                .unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
        let AnalysisResult::OverlapExtensionPrimersDesigned(report) = response.result else {
            panic!("multiple-edit proposal");
        };
        report.validate_binding(&revision).unwrap();
        let stored_expected = if reverse {
            "TGACGATCTACAGGTCCGAAAATACGTAGCATCGGACT"
        } else {
            expected
        };
        assert_eq!(report.expected.sequence.as_str(), stored_expected);
        assert_eq!(
            report.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            report.verification.findings
        );
        assert_eq!(
            report
                .verification
                .products
                .last()
                .unwrap()
                .sequence
                .as_str(),
            stored_expected
        );
        assert_eq!(report.plan.inputs[0].revision, source.revision());
        assert_eq!(report.request.edits.len(), 2);
        assert_eq!(report.request.fragment.reverse, reverse);
    }
}
