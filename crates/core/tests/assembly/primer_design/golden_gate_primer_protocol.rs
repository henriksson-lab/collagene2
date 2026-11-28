use collagene_core::{
    cloning::{
        assembly::recipe::{FragmentSpec, Preparation},
        golden_gate::{
            GoldenGateEnzyme,
            primers::{self as planner, Fragment, Request as DesignRequest},
        },
        ligation::recipe as ligation_recipe,
        restriction::catalog as catalog_digest,
    },
    planning::plan::PlanInput,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
};
use collagene_model::{Molecule, Topology};
fn input() -> PlanInput {
    let molecule = Molecule::dna("original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    PlanInput {
        id: "source".into(),
        revision: molecule.revision(),
        molecule,
    }
}
#[test]
fn type_iis_planning_roundtrips_and_rejects_pre_40_envelopes_and_stale_context() {
    let inputs = vec![input()];
    let target = Molecule::dna(
        "independent intended circle",
        "ACGACCCCGGGGTAAC",
        Topology::Circular,
    )
    .unwrap();
    let request = DesignRequest {
        name: "literal Type IIS proposal".into(),
        fragments: vec![Fragment {
            role: "backbone".into(),
            input_id: "source".into(),
            start: 0,
            end: 16,
            reverse: false,
            forward_length: 8,
            reverse_length: 8,
            right_overhang: "ACGA".into(),
        }],
        enzyme: GoldenGateEnzyme::BsaI,
        clamp: "AAAA".into(),
        spacer: "A".into(),
    };
    let revision = planner::preview_revision(&inputs, &target, &request).unwrap();
    let mut job = Request {
        id: 401,
        revision: revision.clone(),
        job: AnalysisJob::DesignGoldenGatePrimers {
            inputs,
            target,
            request,
        },
    };
    let bytes = protocol::request::Request::encode(&job).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 101);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&39u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    let response = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!(tag, 93);
    let decoded = protocol::response::Response::decode(&response).unwrap();
    let AnalysisResult::GoldenGatePrimersDesigned(report) = decoded.result else {
        panic!("expected source-derived proposal")
    };
    report.validate_binding(&revision).unwrap();
    assert_eq!(report.fragments.len(), 1);
    assert!(matches!(
        &report.fragments[0].preparation,
        Preparation::PcrCatalogDigest {
            fragment_index: 1,
            ..
        }
    ));
    let mut old = response;
    old[4..6].copy_from_slice(&39u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    if let AnalysisJob::DesignGoldenGatePrimers { request, .. } = &mut job.job {
        request.fragments[0].right_overhang = "TTAA".into();
    }
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&job).unwrap())
            .is_err()
    );
}
#[test]
fn new_preparation_cannot_hide_inside_an_old_recipe_envelope() {
    let inputs = vec![input()];
    let target = Molecule::dna("independent", "ACGACCCCGGGG", Topology::Linear).unwrap();
    let recipe = ligation_recipe::Recipe {
        standard: None,
        name: "literal PCR digest".into(),
        fragments: vec![ligation_recipe::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "insert".into(),
                input_id: "source".into(),
                reverse: false,
                preparation: Preparation::PcrCatalogDigest {
                    conditions: None,
                    forward: Primer::new("F", "ACGACCCC", "GGTCTCA").unwrap(),
                    reverse: Primer::new("R", "GTTACCCC", "GGTCTCA").unwrap(),
                    selection: catalog_digest::select(
                        &catalog_digest::builtin_catalog().unwrap(),
                        &["bsai".into()],
                    )
                    .unwrap(),
                    fragment_index: 1,
                    repair: None,
                    template_strand: collagene_core::chemistry::Strand::Top,
                },
            },
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    };
    let revision = ligation_recipe::preview_revision(
        &inputs,
        &target,
        &recipe,
        collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    let job = Request {
        id: 402,
        revision,
        job: AnalysisJob::PreviewLigationRecipe {
            policy: collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
            inputs,
            target,
            recipe,
        },
    };
    let mut bytes = protocol::request::Request::encode(&job).unwrap();
    bytes[4..6].copy_from_slice(&39u16.to_le_bytes());
    assert!(protocol::execution::process_request(&bytes).is_err());
}
