use collagene_core::{
    cloning::restriction::{
        catalog as catalog_digest,
        inspection::InspectionOptions,
        setups::{CleavageMode, Recipe},
    },
    documents::workspace::Workspace,
    material::lineage::EventBinding,
    output::gel::{Calibration, CalibrationPoint, LaneInput, Preparation},
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
};
use collagene_model::{Molecule, Topology};

fn request(pinned: bool) -> Request {
    let molecule = Molecule::dna("AluI template", "AAAAGCTTTT", Topology::Linear).unwrap();
    let frame = Workspace::from_documents("source", vec![molecule.clone()], 0).unwrap();
    let doc = &frame.documents[0];
    let recipe = Recipe {
        origin: collagene_core::material::derivation::SourceOrigin::OriginalReference,
        source: EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        },
        selection: catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &["alui".into()],
        )
        .unwrap(),
        options: InspectionOptions::default(),
        cleavage_mode: CleavageMode::AssumedCompleteMethylationUnassessed,
    };
    let lane = LaneInput {
        label: "Pinned AluI".into(),
        molecule: molecule.clone(),
        preparation: if pinned {
            Preparation::CatalogDigest { recipe }
        } else {
            Preparation::Uncut
        },
        mass_ng: Some(10.),
    };
    let calibration = Calibration {
        label: "Literal illustrative ladder".into(),
        points: vec![
            CalibrationPoint {
                size_bp: 1.,
                migration: 0.9,
            },
            CalibrationPoint {
                size_bp: 1000.,
                migration: 0.1,
            },
        ],
        resolution: 0.01,
    };
    Request {
        id: 7,
        revision: molecule.revision(),
        job: if pinned {
            AnalysisJob::PredictPinnedGel {
                molecule,
                lanes: vec![lane],
                calibration,
            }
        } else {
            AnalysisJob::PredictGel {
                molecule,
                lanes: vec![lane],
                calibration,
            }
        },
    }
}
fn tag(bytes: &[u8]) -> u32 {
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    tag
}

#[test]
fn pinned_worker_replays_nonstarter_geometry_in_current_reports() {
    let req = request(true);
    let wire = protocol::request::Request::encode(&req).unwrap();
    assert_eq!(tag(&wire), 122);
    let output = protocol::execution::process_request(&wire).unwrap();
    assert_eq!(tag(&output), 114);
    let response = protocol::response::Response::decode(&output).unwrap();
    assert_eq!((response.id, response.revision), (req.id, req.revision));
    let AnalysisResult::PinnedGel(report) = response.result else {
        panic!("expected pinned gel")
    };
    let lane = &report.lanes[0];
    assert!(lane.error.is_none(), "{:?}", lane.error);
    assert_eq!(
        lane.fragments
            .iter()
            .map(|f| (f.top_nt, f.bottom_nt))
            .collect::<Vec<_>>(),
        vec![(5, 5), (5, 5)]
    );
    assert_eq!(lane.bands.len(), 1);
    assert_eq!(lane.bands[0].fragment_indexes, vec![0, 1]);
    assert_eq!(
        lane.fragments
            .iter()
            .map(|f| f.amount_ng.unwrap())
            .sum::<f64>(),
        10.
    );
    assert!(lane.warnings.iter().any(|w| w.contains("methylation")));
}
#[test]
fn pinned_gel_rejects_unsupported_request_and_result_versions() {
    let mut new_wire = protocol::request::Request::encode(&request(true)).unwrap();
    let mut response = protocol::execution::process_request(&new_wire).unwrap();
    new_wire[4..6].copy_from_slice(&54u16.to_le_bytes());
    response[4..6].copy_from_slice(&54u16.to_le_bytes());
    assert!(protocol::execution::process_request(&new_wire).is_err());
    assert!(protocol::response::Response::decode(&response).is_err());
}

#[test]
fn stale_pins_and_unpinned_operation_cannot_relabel_catalog_digest() {
    let mut req = request(true);
    if let AnalysisJob::PredictPinnedGel { lanes, .. } = &mut req.job {
        let Preparation::CatalogDigest { recipe } = &mut lanes[0].preparation else {
            unreachable!()
        };
        recipe.source.revision = "a".repeat(64);
    }
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&req).unwrap())
            .is_err()
    );
    let mut req = request(true);
    let AnalysisJob::PredictPinnedGel {
        molecule,
        lanes,
        calibration,
    } = req.job
    else {
        unreachable!()
    };
    req.job = AnalysisJob::PredictGel {
        molecule,
        lanes,
        calibration,
    };
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&req).unwrap())
            .is_err()
    );
}

#[test]
fn pinned_gels_keep_the_ordinary_binary_message_budget() {
    let maximum = collagene_core::formats::sequence::MAX_IMPORT_BYTES;
    let mut wire = vec![0; 42 + maximum + 1];
    wire[4..6].copy_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
    wire[6..10].copy_from_slice(&((maximum + 1) as u32).to_le_bytes());
    for (magic, selector, response) in [(b"CGQ\0", 124u32, false), (b"CGR\0", 116u32, true)] {
        wire[..4].copy_from_slice(magic);
        let prefix = postcard::to_allocvec(&(7u32, "source", selector)).unwrap();
        wire[42..42 + prefix.len()].copy_from_slice(&prefix);
        let rejected = if response {
            protocol::response::Response::decode(&wire).map(|_| ())
        } else {
            protocol::execution::process_request(&wire).map(|_| ())
        };
        assert!(matches!(rejected, Err(collagene_model::Error::TooLarge)));
    }
}
