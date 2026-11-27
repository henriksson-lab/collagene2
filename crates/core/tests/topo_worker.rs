use collagene_core::{
    Molecule, Topology,
    protocol::{self, AnalysisJob as J, AnalysisResult as R, Request},
    topo as t, topo_worker as w,
};
use t::physical::{self as p, projection as v};
fn fixture(four: bool) -> (Vec<t::Input>, t::Request) {
    let inputs = vec![t::Input {
        role: "explicit".into(),
        molecule: Molecule::dna("actual source", "ACGT", Topology::Linear).unwrap(),
    }];
    let count = if four { 4 } else { 2 };
    let strands = (0..count)
        .map(|i| t::PreparedStrand {
            id: format!("strand{i}"),
            segments: vec![t::SourceInterval {
                source_role: "explicit".into(),
                start: 0,
                length: 4,
                reverse: i % 2 == 1,
            }],
            five_prime: p::FivePrime::Hydroxyl,
            three_prime: p::ThreePrime::Hydroxyl,
        })
        .collect();
    let pairing = (0..count / 2)
        .map(|i| p::PairSpan {
            strand_a: i * 2,
            start_a: 0,
            strand_b: i * 2 + 1,
            start_b: 0,
            length: 4,
        })
        .collect();
    (
        inputs,
        t::Request {
            recipe: t::GenerationRecipe {
                schema_version: 1,
                materials: vec![t::MaterialRecipe::Prepared {
                    declaration_id: "declared strands".into(),
                    assumption: "Source correspondence and OH termini are explicit hypotheses"
                        .into(),
                    strands,
                    pairing,
                }],
            },
            steps: vec![],
            engine_limits: p::Limits::default(),
            projection_options: v::Options::default(),
        },
    )
}
fn run(revision: String, job: J) -> (R, u32) {
    let request = protocol::encode_request(&Request {
        id: 3,
        revision: revision.clone(),
        job,
    })
    .unwrap();
    let bytes = protocol::process_request(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    let response = protocol::decode_response(&bytes).unwrap();
    assert_eq!(response.revision, revision);
    (response.result, tag)
}
#[test]
fn streamed_revisions_bind_fresh_simulation_projection_and_explicit_target_policy() {
    let (inputs, request) = fixture(false);
    let batch = w::simulation_revision(&inputs, &request).unwrap();
    assert_eq!(batch, t::batch_revision(&inputs, &request).unwrap());
    let (report, tag) = run(
        batch.clone(),
        J::SimulateTopo {
            inputs: inputs.clone(),
            request: request.clone(),
        },
    );
    assert_eq!(tag, 74);
    let R::TopoSimulation(report) = report else {
        panic!("expected simulation")
    };
    assert_eq!(report.physical.source_bases, 8);
    let selection = t::Selection {
        batch_revision: batch,
        component: 0,
        reference_seed: p::Call {
            strand: 0,
            index: 0,
        },
    };
    let revision = w::projection_revision(&inputs, &request, &selection).unwrap();
    assert_eq!(
        revision,
        t::product_revision(&inputs, &request, &selection).unwrap()
    );
    let (product, tag) = run(
        revision,
        J::ProjectTopoProduct {
            inputs: inputs.clone(),
            request: request.clone(),
            selection: selection.clone(),
        },
    );
    assert_eq!(tag, 76);
    let R::TopoProjection(product) = product else {
        panic!("expected projection")
    };
    assert_eq!(product.reference.unwrap().sequence.as_str(), "ACGT");
    let target = Molecule::dna("independent target", "ACGT", Topology::Linear).unwrap();
    let policy = t::AcceptancePolicy::default();
    let revision = w::target_revision(&inputs, &request, &selection, &target, &policy).unwrap();
    let (assessment, tag) = run(
        revision.clone(),
        J::CheckTopoTarget {
            inputs: inputs.clone(),
            request: request.clone(),
            selection: selection.clone(),
            target: target.clone(),
            policy: policy.clone(),
        },
    );
    assert_eq!(tag, 77);
    let R::TopoTarget(assessment) = assessment else {
        panic!("expected target check")
    };
    assert_eq!(assessment.batch_revision, revision);
    assert!(matches!(
        assessment.outcome,
        t::TargetOutcome::AcceptedUnderDeclaredGraph
    ));
    let mut changed = inputs;
    changed[0]
        .molecule
        .metadata
        .insert("changed".into(), "source identity".into());
    assert!(w::projection_revision(&changed, &request, &selection).is_err());
}
#[test]
fn both_declared_orientation_graphs_and_new_schema_guards_remain_explicit() {
    let (inputs, request) = fixture(true);
    let revision = w::orientations_revision(&inputs, &request, [0, 1], [2, 3], false).unwrap();
    let (report, tag) = run(
        revision.clone(),
        J::InspectTopoOrientations {
            inputs: inputs.clone(),
            request: request.clone(),
            vector: [0, 1],
            insert: [2, 3],
            directional: false,
        },
    );
    assert_eq!(tag, 75);
    let R::TopoOrientations(report) = report else {
        panic!("expected orientations")
    };
    assert_eq!(report.reference.candidates.len(), 2);
    assert_eq!(report.candidate_requests.len(), 2);
    assert!(!report.reference.biochemical_mixture_exhaustive);
    assert_eq!(
        revision,
        t::inspect_orientations(&inputs, &request, [0, 1], [2, 3], false)
            .unwrap()
            .batch_revision
    );
    for tag in 82u32..=85 {
        let payload = postcard::to_allocvec(&(3u32, "binding", tag)).unwrap();
        let mut bytes = b"CGQ\0".to_vec();
        bytes.extend(25u16.to_le_bytes());
        bytes.extend((payload.len() as u32).to_le_bytes());
        bytes.extend(blake3::hash(&payload).as_bytes());
        bytes.extend(payload);
        assert!(matches!(
            protocol::process_request(&bytes),
            Err(collagene_core::Error::Unsupported(_))
        ));
    }
}
