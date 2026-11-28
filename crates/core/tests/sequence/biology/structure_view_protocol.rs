use collagene_core::{
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    structure_view::{
        self as native,
        coordinates::{
            Format,
            scene::{Alternates, Chains, Selection},
        },
    },
    workspace::{self, Workspace},
};
use collagene_model::{Molecule, Topology};
const PDB:&[u8]=b"ATOM      1  CA AALA A  -2B     11.104  13.207   2.100  0.60 20.00           C  \nHETATM    2 CA    CA Z  17      12.000  13.000   2.000  1.00 10.00          CA2+\nCONECT    1    2\nEND\n";
fn request(selection: Option<Selection>) -> Request {
    Request {
        id: 17,
        revision: native::revision(PDB, Format::Pdb, &selection).unwrap(),
        job: AnalysisJob::PrepareStructure {
            bytes: PDB.to_vec(),
            format: Format::Pdb,
            selection,
        },
    }
}
#[test]
fn current_worker_and_session_prepare_source_indexed_coordinates_without_sequence_material() {
    let q = request(None);
    let bytes = protocol::request::Request::encode(&q).unwrap();
    assert_eq!(
        &bytes[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let ((_, _, selector), _) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(selector, 192);
    let mut session = protocol::session::AnalysisSession::default();
    for reply in [
        protocol::execution::process_request(&bytes).unwrap(),
        session.process_request(&bytes).unwrap(),
    ] {
        let ((_, _, selector), _) =
            postcard::take_from_bytes::<(u32, &str, u32)>(&reply[42..]).unwrap();
        assert_eq!(selector, 182);
        let r = protocol::response::Response::decode(&reply).unwrap();
        let AnalysisResult::StructurePrepared(p) = r.result else {
            panic!("structure projection")
        };
        p.validate_binding(&q.revision).unwrap();
        assert_eq!(p.source_hash, blake3::hash(PDB).to_hex().as_str());
        assert_eq!(p.atoms.len(), 2);
        assert_eq!(p.atoms[0].position, [11.104, 13.207, 2.1]);
        assert_eq!(p.scene.edges[0].atoms, [0, 1]);
        assert_eq!(p.atoms[1].formal_charge, Some(2));
    }
    let q = request(Some(Selection {
        model: 1,
        chains: Chains::Author(Some("Z".into())),
        alternates: Alternates::All,
        hydrogens: true,
    }));
    let r = protocol::response::Response::decode(
        &protocol::execution::process_request(&protocol::request::Request::encode(&q).unwrap())
            .unwrap(),
    )
    .unwrap();
    let AnalysisResult::StructurePrepared(p) = r.result else {
        panic!()
    };
    assert_eq!(p.scene.points.len(), 1);
    assert_eq!(p.scene.points[0].atom_index, 1);
    assert_eq!(p.scene.unresolved_connections, 1);
    let mut old = bytes;
    old[4..6].copy_from_slice(&114u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    assert!(session.process_request(&old).is_err());
}
#[test]
fn source_options_geometry_and_scene_index_tampering_fail_transport_validation() {
    let q = request(None);
    let report = native::prepare(PDB, Format::Pdb, None).unwrap();
    let refuse = |p: native::Prepared| {
        let response = Response {
            id: 17,
            revision: q.revision.clone(),
            result: AnalysisResult::StructurePrepared(p),
        };
        let mut bytes = b"CGR\0".to_vec();
        bytes.extend_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
        let body = postcard::to_allocvec(&response).unwrap();
        bytes.extend_from_slice(&(body.len() as u32).to_le_bytes());
        bytes.extend_from_slice(blake3::hash(&body).as_bytes());
        bytes.extend_from_slice(&body);
        assert!(protocol::response::Response::decode(&bytes).is_err());
    };
    let mut p = report.clone();
    p.source_hash = "a".repeat(64);
    refuse(p);
    let mut p = report.clone();
    p.scene.points[0].atom_index = 999;
    refuse(p);
    let mut p = report.clone();
    p.scene.points[0].position[0] += 0.25;
    refuse(p);
    let mut p = report.clone();
    p.scene.points.pop();
    refuse(p);
    let mut p = report.clone();
    p.scene.selection.model = 2;
    refuse(p);
    let mut p = report.clone();
    p.scene.edges[0].atoms = [0, 0];
    refuse(p);
    let mut p = report;
    p.notes.atoms_without_element = 1;
    refuse(p);
    let mut q = request(None);
    q.revision = "b".repeat(64);
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&q).unwrap())
            .is_err()
    );
    let mut changed = PDB.to_vec();
    changed.push(b'\n');
    assert_ne!(
        native::revision(PDB, Format::Pdb, &None).unwrap(),
        native::revision(&changed, Format::Pdb, &None).unwrap()
    );
}
#[test]
fn retained_original_roundtrips_workspace_and_other_attachment_roles_are_not_replaced() {
    let mut frame = Workspace::from_documents(
        "Coordinate evidence",
        vec![Molecule::dna("Unrelated sequence", "ACGT", Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    let before = frame.revisions.clone();
    let hash = native::retain(
        &mut frame.attachments,
        "fixture.pdb",
        PDB.to_vec(),
        Format::Pdb,
    )
    .unwrap();
    assert_eq!(
        native::retain(
            &mut frame.attachments,
            "renamed.pdb",
            PDB.to_vec(),
            Format::Pdb
        )
        .unwrap(),
        hash
    );
    assert_eq!(frame.attachments[&hash].name, "fixture.pdb");
    let loaded = workspace::load_workspace(&workspace::save_workspace(&frame).unwrap()).unwrap();
    assert_eq!(loaded.attachments[&hash].bytes, PDB);
    assert_eq!(loaded.revisions, before);
    assert_eq!(
        native::source_format(&loaded.attachments[&hash].kind),
        Some(Format::Pdb)
    );
    frame.attachments.get_mut(&hash).unwrap().kind = "unrelated-evidence".into();
    assert!(
        native::retain(
            &mut frame.attachments,
            "fixture.pdb",
            PDB.to_vec(),
            Format::Pdb
        )
        .is_err()
    );
    assert_eq!(frame.attachments[&hash].kind, "unrelated-evidence");
    assert!(
        native::retain(
            &mut frame.attachments,
            "bad\nname",
            PDB.to_vec(),
            Format::Pdb
        )
        .is_err()
    );
}
