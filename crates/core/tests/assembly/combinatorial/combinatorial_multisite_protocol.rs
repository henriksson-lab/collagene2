//! Nested method semantics are gated independently of the unchanged plan/report formats.
use collagene_core::{
    cloning::recombination::{
        multisite::{self as adapter, multisite},
        reference::{Family, Limits, Reaction, SiteVariant, SourceRecord},
    },
    documents::workspace::{self, Workspace},
    material::lineage::EventBinding,
    planning::combinatorial::{
        self as design, Bin, Choice, ChoiceSelection, ConstructOverrides, Design, Method,
        SelectedConstruct, TargetDeclaration, batch,
    },
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
        session::AnalysisSession,
    },
};
use collagene_model::{Molecule, Topology};
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> (Workspace, Design, SelectedConstruct) {
    let mut profile = multisite::Profile {
        id: "synthetic-library-LR".into(),
        edition: "literal-v1".into(),
        sources: vec![],
        variants: vec![],
        limitations: vec!["Synthetic reference fixture, not a commercial profile".into()],
    };
    for (id, family, dna) in [
        ("L", Family::L, "ACGGGTACAAATCC"),
        ("R", Family::R, "TGGAGTACAAAGCT"),
    ] {
        profile.sources.push(SourceRecord {
            id: id.into(),
            url: "https://example.org/literal-source".into(),
            raw_sha256: "1".repeat(64),
            vector_name: "synthetic literal reference".into(),
            vector_length: 14,
            start_1based: 1,
            end_1based: 14,
            reverse_extraction: false,
            map_url: "https://example.org/literal-map".into(),
            map_sha256: "2".repeat(64),
        });
        profile.variants.push(SiteVariant {
            id: id.into(),
            family,
            subtype: 1,
            sequence: dna.into(),
            overlap_start: 4,
            source_id: id.into(),
        });
    }
    let w = Workspace::from_documents(
        "Nested-method compatibility fixture",
        vec![
            Molecule::dna("L source", "ACGGGTACAAATCC", Topology::Linear).unwrap(),
            Molecule::dna("R source", "TGGAGTACAAAGCT", Topology::Linear).unwrap(),
            Molecule::dna(
                "Independent literal intended reference",
                "ACGGGTACAAAGCT",
                Topology::Linear,
            )
            .unwrap(),
        ],
        2,
    )
    .unwrap();
    let bins = (0..2)
        .map(|i| Bin {
            id: format!("bin-{i}"),
            label: format!("Parent {i}"),
            choices: vec![Choice {
                anneal_bottom: None,
                generation: None,
                id: format!("choice-{i}"),
                label: format!("Original parent {i}"),
                source: binding(&w, i),
                selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    start: 0,
                    end: 14,
                    reverse: false,
                },
                preparation_reverse: false,
                preparation: None,
                treatment: None,
            }],
        })
        .collect();
    let d = Design {
        material_bindings: vec![],
        id: "reference-library".into(),
        name: "Two-parent reference library".into(),
        circular: false,
        bins,
        method: Method::MultiSiteReference {
            request: adapter::Request {
                reaction: multisite::Request {
                    profile,
                    reaction: Reaction::LR,
                    specificities: vec![1],
                    limits: Limits::default(),
                },
                annotation_limits:
                    collagene_core::cloning::recombination::AnnotationLimits::default(),
            },
            selection: None,
        },
    };
    let row = SelectedConstruct {
        id: "declared-reference".into(),
        name: "Independent source-bound reference".into(),
        choices: (0..2)
            .map(|i| ChoiceSelection {
                bin_id: format!("bin-{i}"),
                choice_id: format!("choice-{i}"),
            })
            .collect(),
        target: TargetDeclaration::Existing(binding(&w, 2)),
        overrides: None,
    };
    (w, d, row)
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
fn request(w: Workspace, d: Design, row: SelectedConstruct) -> Request {
    let matrix = batch::Request {
        design: Some(d),
        rows: vec![batch::RowRequest::Ready(row)],
    };
    Request {
        id: 128,
        revision: batch::begin_revision(&w, &matrix).unwrap(),
        job: AnalysisJob::BeginCombinatorialBatch {
            workspace: w,
            matrix,
        },
    }
}
#[test]
fn multisite_design_and_explicit_override_use_current_batch_encoding() {
    assert!(protocol::schema::CANDIDATE_SCHEMA_VERSION >= 16);
    assert!(workspace::PLAN_SCHEMA_VERSION >= 14);
    assert!(protocol::schema::WORKER_SCHEMA_VERSION >= 69);
    for in_override in [false, true] {
        let (w, mut d, mut row) = fixture();
        if in_override {
            let method = std::mem::replace(
                &mut d.method,
                Method::Ligation {
                    allow_nicked: false,
                },
            );
            row.overrides = Some(ConstructOverrides {
                preparations: vec![],
                method: Some(method),
                circular: None,
            });
        }
        let req = request(w, d, row);
        assert!(
            AnalysisSession::default()
                .process_request(&wire(b"CGQ\0", 68, &req))
                .is_err()
        );
        let bytes = AnalysisSession::default()
            .process_request(&protocol::request::Request::encode(&req).unwrap())
            .unwrap();
        assert!(matches!(
            protocol::response::Response::decode(&bytes).unwrap().result,
            AnalysisResult::CombinatorialBatchStarted(_)
        ));
    }
}

#[test]
fn job155_returns_bound_inspection147_before_component_choice_and_is_never_worker68() {
    let (w, d, row) = fixture();
    let revision = design::inspect_multisite_revision(&d, &row, &w).unwrap();
    let req = Request {
        id: 155,
        revision: revision.clone(),
        job: AnalysisJob::InspectCombinatorialMultiSite {
            workspace: w,
            design: d,
            construct: row,
        },
    };
    let encoded = protocol::request::Request::encode(&req).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!(tag, 150);
    let bytes = protocol::execution::process_request(&encoded).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 142);
    let response = protocol::response::Response::decode(&bytes).unwrap();
    let AnalysisResult::CombinatorialMultiSiteInspection(inspection) = response.result else {
        panic!("expected source-pinned inspection")
    };
    inspection.validate_binding(&revision).unwrap();
    assert_eq!(inspection.source_revision.len(), 64);
    assert_eq!(inspection.batch.reference.candidates.len(), 1);
    assert_eq!(inspection.batch.documents.len(), 2);
    assert!(
        inspection
            .batch
            .documents
            .iter()
            .any(|d| d.molecule.sequence.as_str() == "ACGGGTACAAAGCT")
    );
    assert!(protocol::execution::process_request(&wire(b"CGQ\0", 68, &req)).is_err());
    let mut old = bytes;
    old[4..6].copy_from_slice(&68u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    let mut stale = req;
    stale.revision = "0".repeat(64);
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&stale).unwrap())
            .is_err()
    );
    let stale = Response {
        id: 155,
        revision: "0".repeat(64),
        result: AnalysisResult::CombinatorialMultiSiteInspection(inspection),
    };
    assert!(
        protocol::response::Response::decode(&wire(
            b"CGR\0",
            protocol::schema::WORKER_SCHEMA_VERSION,
            &stale
        ))
        .is_err()
    );
}
fn captured(method: bool) -> batch::Commit {
    let (w, mut d, row) = fixture();
    if !method {
        d.method = Method::Ligation {
            allow_nicked: false,
        };
    }
    let matrix = batch::Request {
        design: Some(d),
        rows: vec![batch::RowRequest::Ready(row)],
    };
    let mut session = batch::Session::begin(w, matrix).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    session.into_capture(&context).unwrap()
}
#[test]
fn nested_method_in_immutable_batch_workspace_is_gated_on_save_load_and_capture_reply() {
    for method in [true, false] {
        let saved = captured(method);
        let record = batch::load(&saved.workspace, &saved.attachment_hash).unwrap();
        assert_eq!(record.rows.len(), 1);
        let save = Request {
            id: 1,
            revision: saved
                .workspace
                .head(&saved.workspace.active_document)
                .unwrap()
                .revision(),
            job: AnalysisJob::SaveWorkspace {
                workspace: saved.workspace.clone(),
            },
        };
        let bytes = wire(b"CGQ\0", 68, &save);
        if method {
            assert!(protocol::execution::process_request(&bytes).is_err());
        } else {
            // This is a newly computed science20 record, not historical68 data.
            assert!(protocol::execution::process_request(&bytes).is_err());
            protocol::execution::process_request(
                &protocol::request::Request::encode(&save).unwrap(),
            )
            .unwrap();
        }
        let archive = workspace::save_workspace(&saved.workspace).unwrap();
        let load = Request {
            id: 2,
            revision: blake3::hash(&archive).to_hex().to_string(),
            job: AnalysisJob::LoadWorkspace { bytes: archive },
        };
        if method {
            assert!(protocol::execution::process_request(&wire(b"CGQ\0", 68, &load)).is_err());
        } else {
            assert!(protocol::execution::process_request(&wire(b"CGQ\0", 68, &load)).is_err());
            protocol::execution::process_request(
                &protocol::request::Request::encode(&load).unwrap(),
            )
            .unwrap();
        }
        let response = Response {
            id: 3,
            revision: saved.capture_revision.clone(),
            result: AnalysisResult::CombinatorialBatchCaptured(saved),
        };
        let current = wire(b"CGR\0", protocol::schema::WORKER_SCHEMA_VERSION, &response);
        protocol::response::Response::decode(&current).unwrap();
        let old = wire(b"CGR\0", 68, &response);
        if method {
            assert!(protocol::response::Response::decode(&old).is_err());
        } else {
            // The method is old; the actual captured report is science20.
            assert!(protocol::response::Response::decode(&old).is_err());
        }
    }
}
