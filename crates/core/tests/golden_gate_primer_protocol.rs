use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::{FragmentSpec, Preparation},
    catalog_digest,
    golden_gate::GoldenGateEnzyme,
    golden_gate_primer_design::{self as planner, Fragment, Request as DesignRequest},
    ligation_recipe,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::PlanInput,
};
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
    let bytes = protocol::encode_request(&job).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 101);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&39u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    let response = protocol::process_request(&bytes).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
    assert_eq!(tag, 93);
    let decoded = protocol::decode_response(&response).unwrap();
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
    assert!(protocol::decode_response(&old).is_err());
    if let AnalysisJob::DesignGoldenGatePrimers { request, .. } = &mut job.job {
        request.fragments[0].right_overhang = "TTAA".into();
    }
    assert!(protocol::process_request(&protocol::encode_request(&job).unwrap()).is_err());
}
#[test]
fn new_preparation_cannot_hide_inside_an_old_recipe_envelope() {
    let inputs = vec![input()];
    let target = Molecule::dna("independent", "ACGACCCCGGGG", Topology::Linear).unwrap();
    let recipe = ligation_recipe::Recipe {
        name: "literal PCR digest".into(),
        fragments: vec![ligation_recipe::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "insert".into(),
                input_id: "source".into(),
                reverse: false,
                preparation: Preparation::PcrCatalogDigest {
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
        collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    let job = Request {
        id: 402,
        revision,
        job: AnalysisJob::PreviewLigationRecipe {
            policy: collagene_core::assembly_recipe::CleavagePolicy::AssumedComplete,
            inputs,
            target,
            recipe,
        },
    };
    let mut bytes = protocol::encode_request(&job).unwrap();
    bytes[4..6].copy_from_slice(&39u16.to_le_bytes());
    assert!(protocol::process_request(&bytes).is_err());
}
