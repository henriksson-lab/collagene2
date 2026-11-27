//! Binary admission and replay for explicit target-free terminal actions.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    chemistry::{EndTreatment, Strand},
    material_derivation as native,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
    run_lineage::EventBinding,
    run_material::MaterialOrigin,
    tailing::TailingBase,
    workspace::{self, Workspace},
};
fn event(w: &Workspace) -> EventBinding {
    let d = w
        .documents
        .iter()
        .find(|d| d.id == w.active_document)
        .unwrap();
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn actual_pcr() -> (Workspace, EventBinding, native::SourceOrigin) {
    let source = Molecule::dna("Literal source", "AAACCC", Topology::Linear).unwrap();
    let w = Workspace::from_documents("PCR source", vec![source.clone()], 0).unwrap();
    let f = Primer::new("Initial F", "AAA", "").unwrap();
    let r = Primer::new("Initial R", "GGG", "").unwrap();
    let product = native::preview_pcr(
        &source,
        &native::SourceOrigin::OriginalReference,
        &f,
        &r,
        Strand::Top,
    )
    .unwrap()
    .remove(0);
    let request = native::Request {
        sources: vec![collagene_core::material_derivation::BoundSource {
            role: collagene_core::material_derivation::SourceRole::Template,
            event: event(&w),
            origin: native::SourceOrigin::OriginalReference,
        }],
        derivation: native::Derivation::ExactPcr {
            forward: f,
            reverse: r,
            template_strand: Strand::Top,
            selection: native::PcrSelection {
                selected_index: 0,
                expected_preview_revision: product.molecule.revision(),
                geometry: product.geometry,
            },
        },
        name: "Actual unmodified PCR".into(),
        save_primers: true,
    };
    let installed = native::prepare_install(&w, &request).unwrap();
    let origin = native::SourceOrigin::Saved(MaterialOrigin::Derived(
        native::capsule(&installed.workspace, &installed.receipt.product).unwrap(),
    ));
    (installed.workspace, installed.receipt.product, origin)
}
fn tag(v: &impl serde::Serialize) -> u32 {
    postcard::take_from_bytes::<u32>(&postcard::to_allocvec(v).unwrap())
        .unwrap()
        .0
}
fn call(job: AnalysisJob, revision: String) -> Response {
    let q = Request {
        id: 31,
        revision: revision.clone(),
        job,
    };
    let encoded = protocol::encode_request(&q).unwrap();
    assert_eq!(
        &encoded[4..6],
        &protocol::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let response = protocol::AnalysisSession::default()
        .process_request(&encoded)
        .unwrap();
    let response = protocol::decode_response(&response).unwrap();
    assert_eq!(response.revision, revision);
    response
}
fn wire(r: &Response) -> Vec<u8> {
    let p = postcard::to_allocvec(r).unwrap();
    let mut b = b"CGR\0".to_vec();
    b.extend(protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&p).as_bytes());
    b.extend(p);
    b
}
#[test]
fn actual_pcr_tail_then_phosphorylation_replay_through_current_binary_without_target_or_new_primers()
 {
    let (mut frame, mut source, mut origin) = actual_pcr();
    for action in [
        native::TerminalAction::Tail3Prime {
            base: TailingBase::T,
        },
        native::TerminalAction::TreatEnds {
            treatment: EndTreatment::Phosphorylate,
        },
    ] {
        let before = workspace::save_workspace(&frame).unwrap();
        let context =
            native::terminal_processing_preview_revision(&frame, &source, &origin, &action)
                .unwrap();
        let job = AnalysisJob::InspectMaterialTerminalProcessing {
            workspace: frame.clone(),
            source: source.clone(),
            origin: origin.clone(),
            action: action.clone(),
        };
        assert_eq!(tag(&job), 189);
        let r = call(job, context.clone());
        assert_eq!(tag(&r.result), 179);
        let AnalysisResult::MaterialTerminalProcessingPreview(preview) = r.result else {
            panic!("terminal preview")
        };
        preview.validate_binding(&context).unwrap();
        assert_eq!(preview.product.sequence.as_str(), "AAACCCT");
        assert!(!preview.incomplete);
        let request = native::Request {
            sources: vec![collagene_core::material_derivation::BoundSource {
                role: collagene_core::material_derivation::SourceRole::Template,
                event: source,
                origin: origin,
            }],
            derivation: native::Derivation::TerminalProcessing { action },
            name: "Explicit terminal descendant".into(),
            save_primers: false,
        };
        let context = native::revision(&frame, &request).unwrap();
        let r = call(
            AnalysisJob::InstallMaterialDerivation {
                workspace: frame.clone(),
                request,
            },
            context,
        );
        let AnalysisResult::MaterialDerivationInstalled(installed) = r.result else {
            panic!("install")
        };
        assert!(installed.receipt.primer_attachment.is_none());
        let reopened = workspace::load_workspace(&installed.archive).unwrap();
        installed.receipt.validate_archive(&reopened).unwrap();
        let fresh = native::reconstruct(&reopened, &installed.receipt.product).unwrap();
        assert_eq!(fresh.reference().sequence.as_str(), "AAACCCT");
        assert!(!fresh.incomplete());
        assert_eq!(workspace::save_workspace(&frame).unwrap(), before);
        source = installed.receipt.product;
        origin = native::SourceOrigin::Saved(MaterialOrigin::Derived(
            native::capsule(&reopened, &source).unwrap(),
        ));
        frame = reopened;
    }
}
#[test]
fn changed_action_source_and_documentary_hash_refuse_without_any_installation() {
    let (frame, source, origin) = actual_pcr();
    let before = workspace::save_workspace(&frame).unwrap();
    let action = native::TerminalAction::Tail3Prime {
        base: TailingBase::A,
    };
    let context =
        native::terminal_processing_preview_revision(&frame, &source, &origin, &action).unwrap();
    let changed = AnalysisJob::InspectMaterialTerminalProcessing {
        workspace: frame.clone(),
        source: source.clone(),
        origin: origin.clone(),
        action: native::TerminalAction::Tail3Prime {
            base: TailingBase::G,
        },
    };
    assert!(
        protocol::process_request(
            &protocol::encode_request(&Request {
                id: 41,
                revision: context.clone(),
                job: changed
            })
            .unwrap()
        )
        .is_err()
    );
    let mut source_alias = source.clone();
    source_alias.event_id = "0".repeat(64);
    assert!(
        protocol::encode_request(&Request {
            id: 42,
            revision: context.clone(),
            job: AnalysisJob::InspectMaterialTerminalProcessing {
                workspace: frame.clone(),
                source: source_alias,
                origin: origin.clone(),
                action: action.clone()
            }
        })
        .is_err()
    );
    let preview = native::preview_terminal_processing(&frame, &source, &origin, &action).unwrap();
    let mut response = Response {
        id: 43,
        revision: context,
        result: AnalysisResult::MaterialTerminalProcessingPreview(preview),
    };
    let AnalysisResult::MaterialTerminalProcessingPreview(report) = &mut response.result else {
        unreachable!()
    };
    report.source_sequence_blake3 = "0".repeat(64);
    assert!(protocol::decode_response(&wire(&response)).is_err());
    let request = native::Request {
        sources: vec![collagene_core::material_derivation::BoundSource {
            role: collagene_core::material_derivation::SourceRole::Template,
            event: source,
            origin: origin,
        }],
        derivation: native::Derivation::TerminalProcessing { action },
        name: "No invented primer order".into(),
        save_primers: true,
    };
    assert!(native::preflight(&frame, &request).is_err());
    assert_eq!(workspace::save_workspace(&frame).unwrap(), before);
}
