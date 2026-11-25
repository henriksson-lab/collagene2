//! Key-bound pool intent cannot masquerade as index-only historical formats.
use collagene_core::{
    Molecule, Topology, cloning_outputs, combinatorial_batch as batch,
    combinatorial_design::{
        self as design, Bin, Choice, ChoiceSelection, ConstructOverrides, Design, Method,
        PoolSelection, PoolSubstrate, PoolSubstrateKind, SelectedConstruct, TargetDeclaration,
    },
    combinatorial_export as export, combinatorial_release as release,
    golden_gate::GoldenGateOptions,
    prepared_golden_gate::Selection,
    protocol::{self, AnalysisJob, AnalysisResult, AnalysisSession, Request, Response},
    run_lineage::{self, EventBinding, RunRecord},
    transfer::{self, Candidate},
    verification::{self, DesignPlan, Operation, Status},
    workspace::{self, Workspace},
};
use std::io::{Cursor, Read, Write};
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
        id: "vector-choice".into(),
        label: "literal original contribution".into(),
        source: EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        },
        start: 11,
        end: 11 + A.len(),
        reverse: false,
        preparation: None,
        treatment: None,
    };
    let mut d = Design {
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
        id: "selected-circle".into(),
        name: "Independent selected circle".into(),
        choices: vec![ChoiceSelection {
            bin_id: "vector".into(),
            choice_id: "vector-choice".into(),
        }],
        target: TargetDeclaration::OriginalRegions {
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
fn index_only(plan: &DesignPlan) -> DesignPlan {
    let mut old = plan.clone();
    for step in &mut old.steps {
        if let Operation::SelectedGoldenGate {
            templates,
            options,
            candidate_index,
            ..
        } = step
        {
            *step = Operation::GoldenGate {
                templates: templates.clone(),
                options: options.clone(),
                candidate_index: *candidate_index,
            };
        }
    }
    old
}
fn plan_schema(archive: &[u8], schema: u32) -> Vec<u8> {
    let mut input = zip::ZipArchive::new(Cursor::new(archive)).unwrap();
    let mut output = zip::ZipWriter::new(Cursor::new(Vec::new()));
    for i in 0..input.len() {
        let mut entry = input.by_index(i).unwrap();
        let name = entry.name().to_owned();
        let mut bytes = Vec::new();
        entry.read_to_end(&mut bytes).unwrap();
        if name == "manifest.json" {
            let mut value: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
            value["plan_schema"] = schema.into();
            bytes = serde_json::to_vec(&value).unwrap();
        }
        output
            .start_file(name, zip::write::SimpleFileOptions::default())
            .unwrap();
        output.write_all(&bytes).unwrap();
    }
    output.finish().unwrap().into_inner()
}
fn rekey(record: &mut RunRecord) {
    record.plan_revision = verification::plan_revision(&record.plan);
    record.run_id = blake3::hash(
        &serde_json::to_vec(&(
            "cloning-run-v1",
            record.schema_version,
            record.plan_schema,
            &record.plan_revision,
            &record.plan,
            &record.target,
            &record.inputs,
            record
                .products
                .iter()
                .map(|p| (&p.step_index, &p.revision))
                .collect::<Vec<_>>(),
            &record.assessment,
        ))
        .unwrap(),
    )
    .to_hex()
    .to_string();
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
fn selected_opcode18_requires_candidate12_and_worker59_without_rejecting_index_only_history() {
    let (w, d, row) = fixture();
    let compiled = design::compile(&d, &row, &w).unwrap();
    assert!(protocol::WORKER_SCHEMA_VERSION >= 59);
    assert_eq!(protocol::CANDIDATE_SCHEMA_VERSION, 12);
    assert_eq!(workspace::PLAN_SCHEMA_VERSION, 10);
    assert_eq!(
        postcard::to_allocvec(compiled.plan.steps.last().unwrap()).unwrap()[0],
        18
    );
    let selected = Candidate::new(compiled.plan.clone());
    let encoded = transfer::save_candidate(&selected).unwrap();
    transfer::load_candidate(&encoded).unwrap();
    assert!(transfer::load_candidate(&wire(b"CGP\0", 11, &selected)).is_err());
    let old = Candidate::new(index_only(&compiled.plan));
    transfer::load_candidate(&wire(b"CGP\0", 11, &old)).unwrap();
    for version in [57, 58] {
        let request = Request {
            id: 18,
            revision: compiled.resolved.target.revision(),
            job: AnalysisJob::Verify {
                plan: compiled.plan.clone(),
                target: compiled.resolved.target.clone(),
            },
        };
        assert!(protocol::process_request(&wire(b"CGQ\0", version, &request)).is_err());
    }
    let request = Request {
        id: 6,
        revision: compiled.resolved.target.revision(),
        job: AnalysisJob::Verify {
            plan: old.plan,
            target: compiled.resolved.target,
        },
    };
    let response = protocol::decode_response(
        &protocol::process_request(&wire(b"CGQ\0", 58, &request)).unwrap(),
    )
    .unwrap();
    assert!(matches!(response.result, AnalysisResult::Verification(_)));
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
        .process_request(&protocol::encode_request(&request).unwrap())
        .unwrap();
    let result = protocol::decode_response(&reply).unwrap();
    assert!(matches!(
        result.result,
        AnalysisResult::CombinatorialBatchStarted(_)
    ));
    let mut ordinary = request.clone();
    let AnalysisJob::BeginCombinatorialBatch { workspace, matrix } = &mut ordinary.job else {
        unreachable!()
    };
    let batch::RowRequest::Ready(row) = &mut matrix.rows[0] else {
        unreachable!()
    };
    row.overrides.as_mut().unwrap().method = None;
    ordinary.revision = batch::begin_revision(workspace, matrix).unwrap();
    for version in [57, 58] {
        let reply = AnalysisSession::default()
            .process_request(&wire(b"CGQ\0", version, &ordinary))
            .unwrap();
        protocol::decode_response(&reply).unwrap();
    }
}
#[test]
fn science15_reports_cannot_claim_worker58() {
    let (w, d, row) = fixture();
    let compiled = design::compile(&d, &row, &w).unwrap();
    let report = verification::verify(&compiled.plan, &compiled.resolved.target);
    assert!(report.engine_version.ends_with("/science-15"));
    assert_eq!(report.status, Status::ConsistentUnderAssumptions);
    let response = Response {
        id: 15,
        revision: compiled.resolved.target.revision(),
        result: AnalysisResult::Verification(report),
    };
    protocol::decode_response(&wire(b"CGR\0", 59, &response)).unwrap();
    assert!(protocol::decode_response(&wire(b"CGR\0", 58, &response)).is_err());
}
#[test]
fn workspace_and_rekeyed_runrecord_cannot_downcast_selected_semantics_to_plan9() {
    let (mut w, d, row) = fixture();
    let compiled = design::compile(&d, &row, &w).unwrap();
    w.candidates.push(Candidate::new(compiled.plan.clone()));
    let archive = workspace::save_workspace(&w).unwrap();
    workspace::load_workspace(&archive).unwrap();
    assert!(workspace::load_workspace(&plan_schema(&archive, 9)).is_err());
    w.candidates[0] = Candidate::new(index_only(&compiled.plan));
    let old = workspace::save_workspace(&w).unwrap();
    workspace::load_workspace(&plan_schema(&old, 9)).unwrap();
    let committed = released();
    let manifest = release::load(&committed.workspace, &committed.manifest_hash).unwrap();
    let mut record =
        run_lineage::load_run(&committed.workspace, &manifest.rows[0].output.run.run_id).unwrap();
    assert_eq!(record.plan_schema, 10);
    record.plan_schema = 9;
    record.assessment.engine_version = "0.1.0/science-14".into();
    rekey(&mut record);
    assert!(
        run_lineage::load_run_record(&serde_json::to_vec(&record).unwrap()).is_err(),
        "valid rekey does not hide selected operation"
    );
    record.plan = index_only(&record.plan);
    rekey(&mut record);
    run_lineage::load_run_record(&serde_json::to_vec(&record).unwrap()).unwrap();
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
    let reply = protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap();
    let response = protocol::decode_response(&reply).unwrap();
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
