//! Generic saved-material TOPO worker input is a fresh replay request, never vector activation authority.
use collagene_core::{
    analysis::Primer,
    material_derivation as native,
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    run_lineage::EventBinding,
    run_material::MaterialOrigin,
    topo::{
        self as t,
        physical::{Enzyme, FivePrime, ThreePrime},
    },
    topo_preparation as prep, topo_run_preparation as saved, verification,
    workspace::Workspace,
};
use collagene_model::{Molecule, Topology};
fn strand(role: &str, length: usize) -> t::PreparedStrand {
    t::PreparedStrand {
        id: role.into(),
        segments: vec![t::SourceInterval {
            source_role: role.into(),
            start: 0,
            length,
            reverse: false,
        }],
        five_prime: FivePrime::Phosphate,
        three_prime: ThreePrime::Phosphotyrosyl {
            enzyme: Enzyme::VacciniaTypeIBTyr274,
            activation_id: role.into(),
        },
    }
}
#[test]
fn derived_insert_template_crosses_current_binary_worker_and_preserves_explicit_vector_geometry() {
    let source = Molecule::dna("Original", "CACCATGACG", Topology::Linear).unwrap();
    let w = Workspace::from_documents("Target-free source", vec![source.clone()], 0).unwrap();
    let d = &w.documents[0];
    let event = EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    };
    let forward = Primer::new("F", "CACCA", "").unwrap();
    let reverse = Primer::new("R", "CGTCA", "").unwrap();
    let preview = verification::simulate_pcr_documents(
        &source,
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap()
    .remove(0);
    let prepared = native::prepare_install(
        &w,
        &native::Request {
            sources: vec![collagene_core::material_derivation::BoundSource {
                role: collagene_core::material_derivation::SourceRole::Template,
                event: event,
                origin: native::SourceOrigin::OriginalReference,
            }],
            derivation: native::Derivation::ExactPcr {
                forward,
                reverse,
                selection: native::PcrSelection {
                    selected_index: 0,
                    expected_preview_revision: preview.molecule.revision(),
                    geometry: preview.geometry,
                },
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            name: "Reusable PCR insert template".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let origin = MaterialOrigin::Derived(
        native::capsule(&prepared.workspace, &prepared.receipt.product).unwrap(),
    );
    let template = prepared.workspace.revisions[&prepared.receipt.product.revision].clone();
    let inputs = vec![
        t::Input {
            role: "vt".into(),
            molecule: Molecule::dna("Declared vector top", "AAGGGGCATCCCTT", Topology::Linear)
                .unwrap(),
        },
        t::Input {
            role: "vb".into(),
            molecule: Molecule::dna(
                "Declared vector bottom",
                "GGTGAAGGGATGCCCCTT",
                Topology::Linear,
            )
            .unwrap(),
        },
        t::Input {
            role: "template".into(),
            molecule: template,
        },
    ];
    let declaration = prep::PreparationRequest {
        mode: prep::PreparationMode::DirectionalCacc,
        source: prep::PreparationSource::Original {
            vector: prep::PreparedVectorDeclaration {
                declaration_id: "Explicit activation".into(),
                assumption: "Synthetic vector activation hypothesis; no supplier stock authority"
                    .into(),
                strands: vec![strand("vt", 14), strand("vb", 18)],
                pairing: vec![t::physical::PairSpan {
                    strand_a: 0,
                    start_a: 0,
                    strand_b: 1,
                    start_b: 4,
                    length: 14,
                }],
            },
            insert: prep::ExactPcrDeclaration {
                material_id: "New TOPO insert PCR".into(),
                template_role: "template".into(),
                forward: Primer::new("F", "CACCA", "").unwrap(),
                reverse: Primer::new("R", "CGTCA", "").unwrap(),
                tail: None,
            },
        },
        engine_limits: Default::default(),
        projection_options: Default::default(),
    };
    let compiled = prep::compile(&declaration).unwrap();
    let task = saved::Task::Orientations {
        request: compiled.request,
        vector: compiled.vector_strands,
        insert: compiled.insert_strands,
        directional: compiled.directional,
    };
    let context = saved::revision(
        &inputs,
        Some(&saved::BoundTemplate {
            role: "template".into(),
            origin: origin.clone(),
        }),
        &task,
    )
    .unwrap();
    let q = Request {
        id: 169,
        revision: context.clone(),
        job: AnalysisJob::InspectTopoRun {
            inputs: inputs.clone(),
            template: Some(saved::BoundTemplate {
                role: "template".into(),
                origin: origin.clone(),
            }),
            task: task.clone(),
        },
    };
    let tag = postcard::take_from_bytes::<u32>(&postcard::to_allocvec(&q.job).unwrap())
        .unwrap()
        .0;
    assert_eq!(tag, 164);
    let encoded = protocol::request::Request::encode(&q).unwrap();
    assert_eq!(
        &encoded[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let encoded = protocol::session::AnalysisSession::default()
        .process_request(&encoded)
        .unwrap();
    assert_eq!(
        &encoded[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::response::Response::decode(&encoded).unwrap();
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(&response.result).unwrap())
            .unwrap()
            .0,
        154
    );
    let AnalysisResult::TopoRunInspection(report) = response.result else {
        panic!("fresh TOPO inspection")
    };
    report
        .validate_for(
            &inputs,
            Some(&saved::BoundTemplate {
                role: "template".into(),
                origin: origin.clone(),
            }),
            &task,
            &context,
        )
        .unwrap();
    let saved::Reply::Orientations(inspection) = report.reply else {
        panic!("orientations")
    };
    assert_eq!(inspection.reference.candidates.len(), 2);
    assert!(!inspection.reference.biochemical_mixture_exhaustive);
    let mut wrong = q.clone();
    wrong.revision = "0".repeat(64);
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&wrong).unwrap())
            .is_err()
    );
    let mut wrong = q;
    let AnalysisJob::InspectTopoRun { template, .. } = &mut wrong.job else {
        unreachable!()
    };
    template.as_mut().unwrap().role = "vt".into();
    assert!(protocol::request::Request::encode(&wrong).is_err());
}
