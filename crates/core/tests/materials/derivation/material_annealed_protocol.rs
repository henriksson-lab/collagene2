//! Current binary inspection and native two-role installation use the same source binding.
use collagene_core::{
    chemistry::AnnealOptions,
    documents::workspace::{self, Workspace},
    material::{derivation as native, lineage::EventBinding},
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Molecule, Topology};
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> Workspace {
    Workspace::from_documents(
        "Actual supplied oligos",
        vec![
            Molecule::dna("Top", "AATTCGTACT", Topology::Linear).unwrap(),
            Molecule::dna("Bottom", "AGTACG", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap()
}
fn tag(x: &impl serde::Serialize) -> u32 {
    postcard::take_from_bytes::<u32>(&postcard::to_allocvec(x).unwrap())
        .unwrap()
        .0
}
fn wire(r: &Response) -> Vec<u8> {
    let p = postcard::to_allocvec(r).unwrap();
    let mut out = b"CGR\0".to_vec();
    out.extend(protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    out.extend((p.len() as u32).to_le_bytes());
    out.extend(blake3::hash(&p).as_bytes());
    out.extend(p);
    out
}
#[test]
fn exact_two_role_preview_install_archive_and_corrupt_reply_refusal() {
    let w = fixture();
    let top = event(&w, 0);
    let bottom = event(&w, 1);
    let options = AnnealOptions {
        min_paired_bases: 6,
        ..Default::default()
    };
    let context = native::annealed_oligos_preview_revision(&w, &top, &bottom, &options).unwrap();
    let job = AnalysisJob::InspectMaterialAnnealedOligos {
        workspace: w.clone(),
        top: top.clone(),
        bottom: bottom.clone(),
        options: options.clone(),
    };
    assert_eq!(tag(&job), 190);
    let response = protocol::session::AnalysisSession::default()
        .process_request(
            &protocol::request::Request::encode(&Request {
                id: 1,
                revision: context.clone(),
                job,
            })
            .unwrap(),
        )
        .unwrap();
    let response = protocol::response::Response::decode(&response).unwrap();
    assert_eq!(tag(&response.result), 180);
    let AnalysisResult::MaterialAnnealedOligosPreview(preview) = response.result else {
        panic!("annealing preview")
    };
    preview.validate_binding(&context).unwrap();
    assert_eq!(preview.products[0].sequence.as_str(), "AATTCGTACT");
    let selection = preview
        .selection(&preview.preview.annealing.candidates[0].id)
        .unwrap();
    let request = native::Request {
        sources: vec![
            native::BoundSource {
                role: native::SourceRole::TopOligo,
                event: top.clone(),
                origin: native::SourceOrigin::OriginalReference,
            },
            native::BoundSource {
                role: native::SourceRole::BottomOligo,
                event: bottom.clone(),
                origin: native::SourceOrigin::OriginalReference,
            },
        ],
        derivation: native::Derivation::AnnealedOligos {
            options: options.clone(),
            selection,
        },
        name: "Binary adapter".into(),
        save_primers: false,
    };
    let revision = native::revision(&w, &request).unwrap();
    let job = AnalysisJob::InstallMaterialDerivation {
        workspace: w.clone(),
        request,
    };
    let result = protocol::execution::process_request(
        &protocol::request::Request::encode(&Request {
            id: 2,
            revision,
            job,
        })
        .unwrap(),
    )
    .unwrap();
    let result = protocol::response::Response::decode(&result).unwrap();
    let AnalysisResult::MaterialDerivationInstalled(installed) = result.result else {
        panic!("installed")
    };
    let reopened = workspace::load_workspace(&installed.archive).unwrap();
    installed.receipt.validate_archive(&reopened).unwrap();
    assert_eq!(
        native::reconstruct(&reopened, &installed.receipt.product)
            .unwrap()
            .tracked()
            .duplex
            .bottom
            .as_str(),
        "AGTACG"
    );
    let mut forged = Response {
        id: 3,
        revision: context.clone(),
        result: AnalysisResult::MaterialAnnealedOligosPreview(preview),
    };
    let AnalysisResult::MaterialAnnealedOligosPreview(report) = &mut forged.result else {
        unreachable!()
    };
    report.sources.swap(0, 1);
    assert!(protocol::response::Response::decode(&wire(&forged)).is_err());
    let changed = Request {
        id: 4,
        revision: context,
        job: AnalysisJob::InspectMaterialAnnealedOligos {
            workspace: w,
            top: bottom,
            bottom: top,
            options,
        },
    };
    assert!(
        protocol::execution::process_request(
            &protocol::request::Request::encode(&changed).unwrap()
        )
        .is_err()
    );
}
