//! Independent ZIP decoding and literal sequence/order/protocol expectations.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_recipe::Preparation,
    chemistry::EndTreatment,
    cloning_outputs, combinatorial_batch as batch,
    combinatorial_design::{
        Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct, TargetDeclaration,
    },
    combinatorial_export::{self as export, Options},
    combinatorial_release::{self as release, ReleaseSession, Request, Selection},
    formats,
    history::DocumentHistory,
    run_lineage::{self, EventBinding},
    workspace::{Workspace, WorkspaceAttachment},
};
use serde_json::Value;
use std::{
    collections::BTreeMap,
    io::{Cursor, Read},
};

fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn released() -> release::Commit {
    let w = Workspace::from_documents(
        "archive fixture",
        vec![
            Molecule::dna("left", "TTGCAATCGTCGATGCTAGC", Topology::Linear).unwrap(),
            Molecule::dna("right", "ATCGTCGATGCTAGCGGTAC", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let d = Design {
        id: "archive-design".into(),
        name: "Bulk design β".into(),
        circular: false,
        method: Method::Ligation {
            allow_nicked: false,
        },
        bins: (0..2)
            .map(|i| {
                let (f, r) = if i == 0 {
                    ("TTGCAATC", "GCTAGCAT")
                } else {
                    ("ATCGTCGA", "GTACCGCT")
                };
                Bin {
                    id: format!("bin-{i}"),
                    label: format!("Bin {i}"),
                    choices: vec![Choice {
                        id: format!("choice-{i}"),
                        label: format!("Choice {i}"),
                        source: binding(&w, i),
                        start: 0,
                        end: 20,
                        reverse: false,
                        preparation: Some(Preparation::ExactPcr {
                            forward: Primer::new(format!("F α{i}"), f, "").unwrap(),
                            reverse: Primer::new(format!("R β{i}"), r, "").unwrap(),
                        }),
                        treatment: Some(EndTreatment::Phosphorylate),
                    }],
                }
            })
            .collect(),
    };
    let choices = (0..2)
        .map(|i| ChoiceSelection {
            bin_id: format!("bin-{i}"),
            choice_id: format!("choice-{i}"),
        })
        .collect::<Vec<_>>();
    let rows = ["..", "second"]
        .into_iter()
        .enumerate()
        .map(|(i, id)| {
            batch::RowRequest::Ready(SelectedConstruct {
                id: id.into(),
                name: format!("Cloning Ω construct{i}"),
                choices: choices.clone(),
                target: TargetDeclaration::OriginalRegions {
                    junction_edits: vec![],
                },
                overrides: None,
            })
        })
        .collect();
    let mut session = batch::Session::begin(
        w,
        batch::Request {
            design: Some(d),
            rows,
        },
    )
    .unwrap();
    let context = session.summary().context_revision;
    for i in 0..2 {
        session.advance(&context, i).unwrap();
    }
    let saved = session.into_capture(&context).unwrap();
    let request = Request {
        batch_hash: saved.attachment_hash,
        rows: ["..", "second"]
            .into_iter()
            .enumerate()
            .map(|(i, id)| Selection {
                construct_id: id.into(),
                options: cloning_outputs::Options {
                    name: format!("../Released β{i}"),
                    show_prepared: false,
                    save_primers: true,
                    allow_incomplete: false,
                },
            })
            .collect(),
    };
    let mut session = ReleaseSession::begin(saved.workspace, request).unwrap();
    let context = session.summary().context_revision;
    for i in 0..2 {
        session.advance(&context, i).unwrap();
    }
    session.into_commit(&context).unwrap()
}
fn unpack(bytes: &[u8]) -> BTreeMap<String, Vec<u8>> {
    let mut zip = zip::ZipArchive::new(Cursor::new(bytes)).unwrap();
    let mut entries = BTreeMap::new();
    for i in 0..zip.len() {
        let mut file = zip.by_index(i).unwrap();
        let name = file.name().to_string();
        assert!(file.enclosed_name().is_some());
        assert!(!name.starts_with('/') && !name.contains("../") && !name.contains('\\'));
        let mut bytes = vec![];
        file.read_to_end(&mut bytes).unwrap();
        assert!(entries.insert(name, bytes).is_none());
    }
    entries
}
fn all() -> Options {
    Options {
        genbank: true,
        fasta: true,
        primer_csv: true,
        primer_tsv: true,
        protocol_html: true,
        protocol_pdf: true,
        construct_ids: None,
    }
}

#[test]
fn all_formats_preserve_literal_dna_oligos_utf8_and_exact_lineage_under_hashed_paths() {
    let release = released();
    let options = all();
    let archive = export::export(&release.workspace, &release.manifest_hash, &options).unwrap();
    archive.validate().unwrap();
    assert_eq!(
        archive.context_revision,
        export::revision(&release.workspace, &release.manifest_hash, &options).unwrap()
    );
    let files = unpack(&archive.bytes);
    assert_eq!(files.len(), 13);
    assert_eq!(archive.files.len(), files.len());
    let manifest: Value = serde_json::from_slice(&files["manifest.json"]).unwrap();
    assert_eq!(manifest["release_manifest_hash"], release.manifest_hash);
    assert_eq!(manifest["rows"].as_array().unwrap().len(), 2);
    assert_eq!(manifest["files"].as_array().unwrap().len(), 12);
    assert_eq!(manifest["rows"][0]["released"]["construct_id"], "..");
    assert_eq!(
        manifest["rows"][0]["fresh_status"],
        "ConsistentUnderAssumptions"
    );
    let released = release::load(&release.workspace, &release.manifest_hash).unwrap();
    assert_eq!(
        manifest["rows"][0]["released"]["inputs"],
        serde_json::to_value(&released.rows[0].inputs).unwrap()
    );
    assert_eq!(
        manifest["rows"][0]["released"]["output"]["construct"],
        serde_json::to_value(&released.rows[0].output.construct).unwrap()
    );
    for summary in &archive.files {
        let bytes = &files[&summary.path];
        assert_eq!(summary.bytes, bytes.len());
        assert_eq!(summary.blake3, blake3::hash(bytes).to_hex().as_str());
        if summary.path.ends_with(".gb") {
            let text = std::str::from_utf8(bytes).unwrap();
            let m = formats::import_genbank(text).unwrap();
            assert_eq!(
                m.molecules[0].sequence.as_str(),
                "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC"
            );
        }
        if summary.path.ends_with(".fasta") {
            let m = formats::import_fasta(std::str::from_utf8(bytes).unwrap()).unwrap();
            assert_eq!(
                m.molecules[0].sequence.as_str(),
                "TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC"
            );
            assert!(!summary.warnings.is_empty());
        }
        if summary.path.ends_with(".csv") || summary.path.ends_with(".tsv") {
            let mut reader = csv::ReaderBuilder::new()
                .delimiter(if summary.path.ends_with(".tsv") {
                    b'\t'
                } else {
                    b','
                })
                .from_reader(bytes.as_slice());
            assert_eq!(reader.headers().unwrap().get(1), Some("Ordered sequence"));
            let records = reader.records().map(Result::unwrap).collect::<Vec<_>>();
            assert_eq!(records.len(), 4);
            let mut calls = records
                .iter()
                .map(|r| r.get(1).unwrap())
                .collect::<Vec<_>>();
            calls.sort();
            let mut expected = vec!["TTGCAATC", "GCTAGCAT", "ATCGTCGA", "GTACCGCT"];
            expected.sort();
            assert_eq!(calls, expected);
            assert!(records.iter().any(|r| r.get(0).unwrap().contains('α')));
        }
        if summary.path.ends_with(".html") {
            let html = std::str::from_utf8(bytes).unwrap();
            assert!(html.contains("Cloning Ω construct"));
            assert!(html.contains("ConsistentUnderAssumptions"));
            assert!(html.contains("TTGCAATC"));
            assert!(html.contains("Assumptions and limits"));
        }
        if summary.path.ends_with(".pdf") {
            assert!(bytes.starts_with(b"%PDF-"));
        }
    }
}

#[test]
fn subset_and_sequence_only_exchange_are_explicitly_historical_and_do_not_rebase_edited_heads() {
    let release = released();
    let manifest = release::load(&release.workspace, &release.manifest_hash).unwrap();
    let old = manifest.rows[0].output.construct.clone();
    let mut changed = release.workspace.clone();
    let mut histories = DocumentHistory::restore_all(&changed).unwrap();
    histories
        .iter_mut()
        .find(|d| d.id == old.document_id)
        .unwrap()
        .commit(
            Molecule::dna("later edited head", "AAAA", Topology::Linear).unwrap(),
            "Later edit",
        )
        .unwrap();
    changed.documents.clear();
    changed.history.clear();
    for history in histories {
        history.contribute(&mut changed).unwrap();
    }
    changed.validate().unwrap();
    let options = Options {
        genbank: false,
        fasta: true,
        primer_csv: false,
        primer_tsv: false,
        protocol_html: false,
        protocol_pdf: false,
        construct_ids: Some(vec!["..".into()]),
    };
    let archive = export::export(&changed, &release.manifest_hash, &options).unwrap();
    let files = unpack(&archive.bytes);
    assert_eq!(files.len(), 2);
    let metadata: Value = serde_json::from_slice(&files["manifest.json"]).unwrap();
    assert!(metadata["rows"][0]["fresh_status"].is_null());
    assert_eq!(
        metadata["rows"][0]["released"]["output"]["construct"],
        serde_json::to_value(&old).unwrap()
    );
    let fasta = files
        .iter()
        .find(|(name, _)| name.ends_with(".fasta"))
        .unwrap()
        .1;
    assert!(
        std::str::from_utf8(fasta)
            .unwrap()
            .contains("TTGCAATCGTCGATGCTAGCATCGTCGATGCTAGCGGTAC")
    );
    assert_eq!(
        changed.head(&old.document_id).unwrap().sequence.as_str(),
        "AAAA"
    );
    assert_ne!(
        archive.context_revision,
        export::revision(&release.workspace, &release.manifest_hash, &options).unwrap()
    );
    assert!(
        archive
            .warnings
            .iter()
            .any(|w| w.contains("Historical sequence"))
    );
}

#[test]
fn unsafe_catalogs_tampered_bytes_wrong_bindings_and_an_earlier_preparation_never_become_archives()
{
    let release = released();
    let before = run_lineage::workspace_fingerprint(&release.workspace).unwrap();
    let mut archive = export::export(&release.workspace, &release.manifest_hash, &all()).unwrap();
    archive.bytes[20] ^= 1;
    assert!(archive.validate().is_err());
    let mut archive = export::export(
        &release.workspace,
        &release.manifest_hash,
        &Options {
            protocol_pdf: false,
            ..all()
        },
    )
    .unwrap();
    archive.files[0].path = "../escape.gb".into();
    assert!(archive.validate().is_err());
    let mut options = all();
    options.construct_ids = Some(vec!["missing".into()]);
    assert!(export::export(&release.workspace, &release.manifest_hash, &options).is_err());
    assert!(export::export(&release.workspace, &"0".repeat(64), &all()).is_err());
    let mut changed = release.workspace.clone();
    let mut manifest = release::load(&changed, &release.manifest_hash).unwrap();
    let prefix = &manifest.rows[0].output.run.products[0];
    manifest.rows[0].output.construct = EventBinding {
        document_id: prefix.document_id.clone(),
        event_id: prefix.event_id.clone(),
        revision: prefix.revision.clone(),
    };
    let bytes = serde_json::to_vec(&manifest).unwrap();
    let hash = blake3::hash(&bytes).to_hex().to_string();
    changed.attachments.insert(
        hash.clone(),
        WorkspaceAttachment {
            kind: release::KIND.into(),
            name: "tampered prefix".into(),
            media_type: "application/json".into(),
            bytes,
        },
    );
    assert!(export::export(&changed, &hash, &all()).is_err());
    assert_eq!(
        run_lineage::workspace_fingerprint(&release.workspace).unwrap(),
        before
    );
}
