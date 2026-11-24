//! Literal append-only transport checks; readonly annotation DTOs are not replay authority.
use collagene_core::{
    Error, Feature, Location, Molecule, Topology,
    protocol::{self, AnalysisJob as J, AnalysisResult as R, Request},
    topo as t, topo_worker as w,
};
use t::physical as p;

fn fixture() -> (Vec<t::Input>, t::Request, t::Selection) {
    let mut source = Molecule::dna("literal source", "ACGT", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Range {
            start: 0,
            end: 2,
            before: false,
            after: false,
        },
        qualifiers: vec![("label".into(), Some("actual top annotation".into()))],
        needs_review: false,
    });
    let inputs = vec![t::Input {
        role: "source".into(),
        molecule: source,
    }];
    let strands = (0..2)
        .map(|index| t::PreparedStrand {
            id: format!("strand{index}"),
            segments: vec![t::SourceInterval {
                source_role: "source".into(),
                start: 0,
                length: 4,
                reverse: index == 1,
            }],
            five_prime: p::FivePrime::Hydroxyl,
            three_prime: p::ThreePrime::Hydroxyl,
        })
        .collect();
    let request = t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![t::MaterialRecipe::Prepared {
                declaration_id: "literal paired interval".into(),
                assumption: "Explicit antiparallel reference correspondence and OH end hypothesis"
                    .into(),
                strands,
                pairing: vec![p::PairSpan {
                    strand_a: 0,
                    start_a: 0,
                    strand_b: 1,
                    start_b: 0,
                    length: 4,
                }],
            }],
        },
        transfers: vec![],
        engine_limits: Default::default(),
        projection_options: Default::default(),
    };
    let selection = t::Selection {
        batch_revision: w::simulation_revision(&inputs, &request).unwrap(),
        component: 0,
        reference_seed: p::Call {
            strand: 0,
            index: 0,
        },
    };
    (inputs, request, selection)
}
fn run(revision: String, job: J, job_tag: u32, result_tag: u32) -> R {
    let encoded = protocol::encode_request(&Request {
        id: 17,
        revision: revision.clone(),
        job,
    })
    .unwrap();
    assert_eq!(
        u16::from_le_bytes([encoded[4], encoded[5]]),
        protocol::WORKER_SCHEMA_VERSION
    );
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, job_tag);
    let mut old = encoded.clone();
    old[4..6].copy_from_slice(&26u16.to_le_bytes());
    assert!(matches!(
        protocol::process_request(&old),
        Err(Error::Unsupported(_))
    ));
    let bytes = protocol::process_request(&encoded).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, result_tag);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&26u16.to_le_bytes());
    assert!(matches!(
        protocol::decode_response(&old),
        Err(Error::Unsupported(_))
    ));
    let response = protocol::decode_response(&bytes).unwrap();
    assert_eq!(response.revision, revision);
    response.result
}
#[test]
fn annotated_jobs_append_without_relabeling_physical_reports_or_restoring_authority() {
    let (inputs, request, selection) = fixture();
    let revision = w::annotated_projection_revision(&inputs, &request, &selection).unwrap();
    assert_eq!(
        revision,
        t::annotated_product_revision(&inputs, &request, &selection).unwrap()
    );
    let R::AnnotatedTopoProjection(product) = run(
        revision.clone(),
        J::ProjectAnnotatedTopo {
            inputs: inputs.clone(),
            request: request.clone(),
            selection: selection.clone(),
        },
        86,
        78,
    ) else {
        panic!("expected annotated projection")
    };
    let reference = product.reference.as_ref().unwrap();
    assert_eq!(reference.sequence.as_str(), "ACGT");
    assert_eq!(reference.features.len(), 1);
    assert_eq!(reference.features[0].label(), "actual top annotation");
    assert!(
        product
            .physical
            .reference
            .as_ref()
            .unwrap()
            .features
            .is_empty()
    );
    assert_eq!(
        product.annotations.as_ref().unwrap().sources[0].projected,
        1
    );
    assert_eq!(
        reference
            .metadata
            .get("source:legacy:chemistry_unmapped")
            .map(String::as_str),
        Some("true")
    );
    let mut forged = product.clone();
    forged.reference.as_mut().unwrap().features.clear();
    assert!(w::validate_annotated_projection(&forged, &revision).is_err());
    let target = Molecule::dna("independently supplied target", "ACGT", Topology::Linear).unwrap();
    let policy = t::AcceptancePolicy::default();
    let revision =
        w::annotated_target_revision(&inputs, &request, &selection, &target, &policy).unwrap();
    let R::AnnotatedTopoTarget(assessment) = run(
        revision,
        J::CheckAnnotatedTopoTarget {
            inputs,
            request,
            selection,
            target,
            policy,
        },
        87,
        79,
    ) else {
        panic!("expected annotated target")
    };
    assert!(matches!(
        assessment.outcome,
        t::TargetOutcome::AcceptedUnderDeclaredGraph
    ));
    assert_eq!(
        assessment
            .product
            .reference
            .as_ref()
            .unwrap()
            .features
            .len(),
        1
    );
}
