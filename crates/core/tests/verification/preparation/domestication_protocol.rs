//! End-to-end binary proposal/reference/preparation transport and version guards.
use collagene_core::{
    cloning::{
        golden_gate::domestication::{self as design, preparation as physical},
        restriction::catalog as catalog_digest,
    },
    planning::plan::PlanInput,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::codon_design,
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
fn fixture() -> (Molecule, design::Request, design::Selection) {
    let source = Molecule::dna(
        "Original template",
        "AGTCCGATGGTCTCGTAGGACTTCAGACCTGTA",
        Topology::Circular,
    )
    .unwrap();
    let request = design::Request {
        source_revision: source.revision(),
        selection: catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &["bsai".into()],
        )
        .unwrap(),
        scope: vec![design::SourceRange {
            start: 0,
            end: source.sequence.len(),
        }],
        retained_sites: vec![],
        coding_features: vec![],
        protected: vec![],
        options: codon_design::DesignOptions::default(),
        limits: design::Limits::default(),
    };
    let patches = vec![design::Patch {
        position: 10,
        expected: "T".into(),
        replacement: "C".into(),
    }];
    let selection = design::Selection {
        context_revision: design::propose_revision(&source, &request).unwrap(),
        candidate_id: design::candidate_id(&request, &patches).unwrap(),
        patches,
    };
    (source, request, selection)
}
fn selector(bytes: &[u8]) -> u32 {
    postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..])
        .unwrap()
        .0
        .2
}
fn exchange(job: AnalysisJob, revision: String, job_tag: u32, result_tag: u32) -> AnalysisResult {
    let encoded = protocol::request::Request::encode(&Request {
        id: 610,
        revision,
        job,
    })
    .unwrap();
    assert_eq!(selector(&encoded), job_tag);
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&60u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    let response = protocol::execution::process_request(&encoded).unwrap();
    assert_eq!(selector(&response), result_tag);
    let mut old = response.clone();
    old[4..6].copy_from_slice(&60u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    protocol::response::Response::decode(&response)
        .unwrap()
        .result
}
#[test]
fn inspection_proposal_check_and_original_template_preparation_keep_append_only_selectors() {
    let (source, request, selection) = fixture();
    let AnalysisResult::DomesticationInspected(inspection) = exchange(
        AnalysisJob::InspectDomestication {
            source: source.clone(),
            request: request.clone(),
        },
        design::inspect_revision(&source, &request).unwrap(),
        137,
        129,
    ) else {
        panic!("inspection required")
    };
    assert_eq!(inspection.sites.len(), 1);
    let AnalysisResult::DomesticationProposed(proposal) = exchange(
        AnalysisJob::ProposeDomestication {
            source: source.clone(),
            request: request.clone(),
        },
        design::propose_revision(&source, &request).unwrap(),
        138,
        130,
    ) else {
        panic!("proposal required")
    };
    assert!(
        proposal.candidate.is_none(),
        "No coding feature was selected for synonymous search."
    );
    let AnalysisResult::DomesticationChecked { report, target } = exchange(
        AnalysisJob::CheckDomestication {
            source: source.clone(),
            request: request.clone(),
            selection: selection.clone(),
        },
        design::check_revision(&source, &request, &selection).unwrap(),
        139,
        131,
    ) else {
        panic!("checked literal reference required")
    };
    assert_eq!(
        target.sequence.as_str(),
        "AGTCCGATGGCCTCGTAGGACTTCAGACCTGTA"
    );
    assert_eq!(report.target_revision, target.revision());
    let input = PlanInput {
        id: physical::INPUT_ID.into(),
        revision: source.revision(),
        molecule: source,
    };
    let request = physical::Request {
        origin: None,
        design: request,
        selection,
        route: physical::Route::InversePcr {
            edit_start: 10,
            edit_end: 11,
            forward_length: 8,
            reverse_length: 8,
            reverse_tail_bases: 0,
        },
    };
    let revision = physical::revision(&input, &request, None).unwrap();
    let AnalysisResult::DomesticationPrepared(prepared) = exchange(
        AnalysisJob::PrepareDomestication {
            workspace: None,
            input: input.clone(),
            request: request.clone(),
        },
        revision,
        140,
        132,
    ) else {
        panic!("actual original-source preparation required")
    };
    assert_eq!(prepared.target, target);
    assert_eq!(prepared.plan.inputs[0].molecule, input.molecule);
    assert_eq!(
        prepared.verification.status,
        Status::ConsistentUnderAssumptions
    );
    let stale = protocol::request::Request::encode(&Request {
        id: 611,
        revision: "0".repeat(64),
        job: AnalysisJob::PrepareDomestication {
            workspace: None,
            input,
            request,
        },
    })
    .unwrap();
    assert!(protocol::execution::process_request(&stale).is_err());
}
