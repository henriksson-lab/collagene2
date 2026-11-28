//! Key-bound pool intent cannot masquerade as index-only historical formats.
use collagene_core::{
    cloning::golden_gate::{GoldenGateOptions, prepared::Selection},
    documents::workspace::Workspace,
    material::{lineage::EventBinding, outputs as cloning_outputs},
    planning::combinatorial::{
        self as design, Bin, Choice, ChoiceSelection, ConstructOverrides, Design, Method,
        PoolSelection, PoolSubstrate, PoolSubstrateKind, SelectedConstruct, TargetDeclaration,
        batch, export, release,
    },
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
        session::AnalysisSession,
    },
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
use std::io::{Cursor, Read};
const A: &str = "ATGACCCCGGGGTTTCCAGA";
fn fixture() -> (Workspace, Design, SelectedConstruct) {
    let source = Molecule::dna(
        "whole original vector",
        &format!("ATGCGGTCTCA{A}ATGATGAGACCGCAT"),
        Topology::Circular,
    )
    .unwrap();
    let w = Workspace::from_documents("selected reaction format fixture", vec![source], 0).unwrap();
    let doc = &w.documents[0];
    let choice = Choice {
        anneal_bottom: None,
        generation: Some(
            collagene_core::chemistry::tracked_duplex::GenerationAssumptions {
                host: None,
                dam: collagene_core::chemistry::methylation::MethylationState::Unmodified,
                dcm: collagene_core::chemistry::methylation::MethylationState::Unmodified,
                cpg: collagene_core::chemistry::methylation::MethylationState::Unmodified,
            },
        ),
        id: "vector-choice".into(),
        label: "literal original contribution".into(),
        source: EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        },
        selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
            start: 11,
            end: 11 + A.len(),
            reverse: false,
        },
        preparation_reverse: false,
        preparation: None,
        treatment: None,
    };
    let mut d = Design {
        material_bindings: vec![],
        id: "pool-wire".into(),
        name: "one chosen whole-source product".into(),
        circular: true,
        bins: vec![Bin {
            id: "vector".into(),
            label: "vector".into(),
            choices: vec![choice],
        }],
        method: Method::GoldenGatePool {
            options: GoldenGateOptions {
                max_fragments: 2,
                ..Default::default()
            },
            substrates: vec![PoolSubstrate {
                bin_id: "vector".into(),
                kind: PoolSubstrateKind::WholeSource,
            }],
            selection: None,
        },
    };
    let row = SelectedConstruct {
        standard: None,
        id: "selected-circle".into(),
        name: "Independent selected circle".into(),
        choices: vec![ChoiceSelection {
            bin_id: "vector".into(),
            choice_id: "vector-choice".into(),
        }],
        target: TargetDeclaration::OriginalRegions {
            substitutions: vec![],
            edits: vec![],
            junction_edits: vec![],
        },
        overrides: None,
    };
    let report = design::inspect_pool(&d, &row, &w).unwrap();
    let (index, c) = report
        .pool
        .as_ref()
        .unwrap()
        .candidates
        .iter()
        .enumerate()
        .find(|(_, c)| c.restriction_site_free && c.duplex.top.as_str() == A)
        .expect("literal one-fragment circle");
    let selection = Selection {
        context_revision: report.context_revision,
        candidate_index: index,
        directed_sequence_key: c.directed_sequence_key.clone(),
        physical_sequence_key: c.physical_sequence_key.clone(),
    };
    let resolved = design::resolve_pool(&d, &row, &w).unwrap();
    let design::Recipe::GoldenGatePool {
        source_revision, ..
    } = resolved.recipe
    else {
        unreachable!()
    };
    let Method::GoldenGatePool { selection: pin, .. } = &mut d.method else {
        unreachable!()
    };
    *pin = Some(PoolSelection {
        selection,
        source_revision,
    });
    (w, d, row)
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    let payload = postcard::to_allocvec(value).unwrap();
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((payload.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&payload).as_bytes());
    b.extend(payload);
    b
}

fn released() -> release::Commit {
    let (w, d, row) = fixture();
    let mut analysis = batch::Session::begin(
        w,
        batch::Request {
            design: Some(d),
            rows: vec![batch::RowRequest::Ready(row)],
        },
    )
    .unwrap();
    let context = analysis.summary().context_revision;
    assert_eq!(
        analysis.advance(&context, 0).unwrap().row.unwrap().status,
        Some(Status::ConsistentUnderAssumptions)
    );
    let saved = analysis.into_capture(&context).unwrap();
    let request = release::Request {
        batch_hash: saved.attachment_hash,
        rows: vec![release::Selection {
            construct_id: "selected-circle".into(),
            options: cloning_outputs::Options {
                name: "Released literal circle".into(),
                ..Default::default()
            },
        }],
    };
    let mut session = release::ReleaseSession::begin(saved.workspace, request).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    assert_eq!(session.inspect(0).unwrap().product.sequence.as_str(), A);
    session.into_commit(&context).unwrap()
}

#[test]
fn pool_method_cannot_hide_in_per_construct_overrides_under_worker57_or58() {
    let (w, mut d, mut row) = fixture();
    let pool = std::mem::replace(
        &mut d.method,
        Method::Ligation {
            allow_nicked: false,
        },
    );
    row.overrides = Some(ConstructOverrides {
        preparations: vec![],
        method: Some(pool),
        circular: None,
    });
    let matrix = batch::Request {
        design: Some(d),
        rows: vec![batch::RowRequest::Ready(row)],
    };
    let request = Request {
        id: 128,
        revision: batch::begin_revision(&w, &matrix).unwrap(),
        job: AnalysisJob::BeginCombinatorialBatch {
            workspace: w,
            matrix,
        },
    };
    for version in [57, 58] {
        assert!(
            AnalysisSession::default()
                .process_request(&wire(b"CGQ\0", version, &request))
                .is_err()
        );
    }
    let reply = AnalysisSession::default()
        .process_request(&protocol::request::Request::encode(&request).unwrap())
        .unwrap();
    let result = protocol::response::Response::decode(&reply).unwrap();
    assert!(matches!(
        result.result,
        AnalysisResult::CombinatorialBatchStarted(_)
    ));
}

#[test]
fn historical_batch_selected_pool_releases_and_exports_protocols_through_binary_worker() {
    let committed = released();
    let options = export::Options {
        genbank: true,
        fasta: true,
        primer_csv: false,
        primer_tsv: false,
        protocol_html: true,
        protocol_pdf: false,
        construct_ids: None,
    };
    let request = Request {
        id: 135,
        revision: export::revision(&committed.workspace, &committed.manifest_hash, &options)
            .unwrap(),
        job: AnalysisJob::ExportCombinatorialOutputs {
            workspace: committed.workspace,
            manifest_hash: committed.manifest_hash,
            options,
        },
    };
    let reply = protocol::execution::process_request(
        &protocol::request::Request::encode(&request).unwrap(),
    )
    .unwrap();
    let response = protocol::response::Response::decode(&reply).unwrap();
    let AnalysisResult::CombinatorialOutputsArchive(archive) = response.result else {
        panic!("missing archive")
    };
    archive.validate().unwrap();
    assert!(archive.bytes.starts_with(b"PK"));
    let mut zip = zip::ZipArchive::new(Cursor::new(&archive.bytes)).unwrap();
    let mut html = String::new();
    let mut fasta = String::new();
    for i in 0..zip.len() {
        let mut entry = zip.by_index(i).unwrap();
        if entry.name().ends_with(".html") {
            entry.read_to_string(&mut html).unwrap();
        } else if entry.name().ends_with(".fasta") {
            entry.read_to_string(&mut fasta).unwrap();
        }
    }
    assert!(html.contains("Golden Gate"));
    assert!(html.contains("candidate"));
    assert!(fasta.contains(A));
}
