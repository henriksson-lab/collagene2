//! The same actual DNA junction has different explicit code-1/code-2 translation.
use collagene_core::{
    cloning::assembly::junctions::{self as assembly_junctions, Frame, Options},
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
fn fixture() -> (DesignPlan, Molecule) {
    let inputs = [("left", "AAAT"), ("right", "TGATTT")]
        .into_iter()
        .map(|(id, dna)| {
            let molecule = Molecule::dna(id, dna, Topology::Linear).unwrap();
            PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule,
            }
        })
        .collect();
    let target = Molecule::dna("independent target", "AAATGATTT", Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "literal stop-codon junction".into(),
        target_revision: target.revision(),
        inputs,
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("left".into()), Source::Input("right".into())],
            overlaps: vec![1],
            circular: false,
        }],
    };
    (plan, target)
}
fn options(code: u8) -> Options {
    Options {
        genetic_code: code,
        frames: vec![Frame {
            reverse: false,
            frame: 0,
        }],
        include_cds: false,
        ..Default::default()
    }
}
fn request(plan: DesignPlan, target: Molecule, options: Options) -> Vec<u8> {
    let revision = assembly_junctions::inspection_revision(&plan, &target, &options).unwrap();
    protocol::request::Request::encode(&Request {
        id: 996,
        revision,
        job: AnalysisJob::InspectAssemblyJunctions {
            plan,
            target,
            options,
        },
    })
    .unwrap()
}
#[test]
fn ordinary_translation_uses_actual_junction_dna_and_explicit_genetic_code() {
    let (plan, target) = fixture();
    for (code, protein) in [(1, "K*F"), (2, "KWF")] {
        let bytes = request(plan.clone(), target.clone(), options(code));
        let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
        assert_eq!(tag, 94);
        let response = protocol::execution::process_request(&bytes).unwrap();
        let decoded = protocol::response::Response::decode(&response).unwrap();
        let AnalysisResult::AssemblyJunctionInspection(report) = decoded.result else {
            panic!("expected junction report")
        };
        assert_eq!(
            report.verification.status,
            Status::ConsistentUnderAssumptions
        );
        assert_eq!(report.junctions.len(), 1);
        assert_eq!(report.junctions[0].homology[0].start, 3);
        assert_eq!(report.junctions[0].homology[0].end, 4);
        assert_eq!(report.junctions[0].window.dna, "AAATGATTT");
        assert_eq!(report.junctions[0].frames[0].translation.protein, protein);
        assert_eq!(
            report.junctions[0].frames[0]
                .crossing_codon
                .as_ref()
                .unwrap()
                .dna,
            "TGA"
        );
        for old in [1u16, 32] {
            let mut downgraded = bytes.clone();
            downgraded[4..6].copy_from_slice(&old.to_le_bytes());
            assert!(protocol::execution::process_request(&downgraded).is_err());
            let mut downgraded = response.clone();
            downgraded[4..6].copy_from_slice(&old.to_le_bytes());
            assert!(protocol::response::Response::decode(&downgraded).is_err());
        }
    }
}
#[test]
fn target_and_code_changes_invalidate_context_but_a_fresh_wrong_target_keeps_actual_dna_evidence() {
    let (plan, target) = fixture();
    let old = assembly_junctions::inspection_revision(&plan, &target, &options(1)).unwrap();
    let changed = protocol::request::Request::encode(&Request {
        id: 996,
        revision: old,
        job: AnalysisJob::InspectAssemblyJunctions {
            plan: plan.clone(),
            target: target.clone(),
            options: options(2),
        },
    })
    .unwrap();
    assert!(protocol::execution::process_request(&changed).is_err());
    let wrong = target.replace(8, 9, "G").unwrap();
    let mut wrong_plan = plan;
    wrong_plan.target_revision = wrong.revision();
    let response =
        protocol::execution::process_request(&request(wrong_plan, wrong, options(1))).unwrap();
    let AnalysisResult::AssemblyJunctionInspection(report) =
        protocol::response::Response::decode(&response)
            .unwrap()
            .result
    else {
        panic!("expected report")
    };
    assert_eq!(report.verification.status, Status::Invalid);
    assert_eq!(report.junctions[0].window.dna, "AAATGATTT");
}
