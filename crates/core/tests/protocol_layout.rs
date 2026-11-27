use collagene_core::{
    Molecule, Topology,
    assembly_recipe::Preparation,
    chemistry::Terminus,
    combinatorial_batch as batch,
    combinatorial_design::{
        Bin, Choice, ChoiceSelection, Design, Method, SelectedConstruct, TargetDeclaration,
    },
    combinatorial_release as release,
    history::DocumentHistory,
    protocol_layout as p,
    run_lineage::{self, EventBinding},
    workspace::{self, Workspace},
};
#[path = "common/pdf.rs"]
mod pdf_readback;
use pdf_readback::{check_xref, extracted};
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> (Workspace, String) {
    let w = Workspace::from_documents(
        "Actual native batch",
        vec![
            Molecule::dna("Left", "AACC", Topology::Linear).unwrap(),
            Molecule::dna("Right", "TTAA", Topology::Linear).unwrap(),
            Molecule::dna("Unused but declared choice", "GGGG", Topology::Linear).unwrap(),
            Molecule::dna("Unrelated private", "CCGG", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let prepare = Some(Preparation::DeclaredDuplex {
        left: Terminus::CLEAVED,
        right: Terminus::CLEAVED,
    });
    let choice = |id: &str, index| Choice {
        generation: None,
        id: id.into(),
        label: id.into(),
        source: event(&w, index),
        start: 0,
        end: 4,
        reverse: false,
        preparation: prepare.clone(),
        treatment: None,
    };
    let design = Design {
        id: "native-layout-batch".into(),
        name: "Actual intended eight bases".into(),
        circular: false,
        method: Method::Ligation {
            allow_nicked: false,
        },
        bins: vec![
            Bin {
                id: "left".into(),
                label: "Left".into(),
                choices: vec![choice("left", 0), choice("unused", 2)],
            },
            Bin {
                id: "right".into(),
                label: "Right".into(),
                choices: vec![choice("right", 1)],
            },
        ],
    };
    let rows = ["alpha", "beta"]
        .into_iter()
        .map(|id| {
            batch::RowRequest::Ready(SelectedConstruct {
                id: id.into(),
                name: format!("Construct {id}"),
                choices: vec![
                    ChoiceSelection {
                        bin_id: "left".into(),
                        choice_id: "left".into(),
                    },
                    ChoiceSelection {
                        bin_id: "right".into(),
                        choice_id: "right".into(),
                    },
                ],
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
            design: Some(design),
            rows,
        },
    )
    .unwrap();
    let context = session.summary().context_revision;
    for i in 0..2 {
        let step = session.advance(&context, i).unwrap();
        assert_eq!(
            step.row.unwrap().status,
            Some(collagene_core::verification::Status::ConsistentUnderAssumptions)
        );
    }
    let commit = session.into_capture(&context).unwrap();
    (commit.workspace, commit.attachment_hash)
}
fn planned(w: &Workspace, h: &str, id: &str) -> p::Source {
    let b = batch::load(w, h).unwrap();
    let row = b.rows.iter().find(|r| r.construct_id == id).unwrap();
    p::Source::Planned {
        batch_hash: h.into(),
        construct_id: id.into(),
        declaration_revision: row.declaration_revision.clone().unwrap(),
    }
}
fn request(w: &Workspace, h: &str, plate: p::Plate) -> p::Request {
    let (rows, columns) = plate.dimensions();
    p::Request {
        name: "Reviewed layout".into(),
        plate,
        assignments: vec![
            p::Assignment {
                id: "first".into(),
                well: p::Well { row: 0, column: 0 },
                label: "alpha, \"one\" <review>".into(),
                source: planned(w, h, "alpha"),
            },
            p::Assignment {
                id: "replicate".into(),
                well: p::Well {
                    row: rows - 1,
                    column: columns - 1,
                },
                label: "Explicit replicate".into(),
                source: planned(w, h, "alpha"),
            },
        ],
    }
}
fn released(w: &Workspace, h: &str) -> (Workspace, String) {
    let r = release::Request {
        batch_hash: h.into(),
        rows: ["alpha", "beta"]
            .into_iter()
            .map(|id| release::Selection {
                construct_id: id.into(),
                options: collagene_core::cloning_outputs::Options {
                    name: format!("Released {id}"),
                    show_prepared: false,
                    save_primers: false,
                    allow_incomplete: false,
                },
            })
            .collect(),
    };
    let mut session = release::ReleaseSession::begin(w.clone(), r).unwrap();
    let ctx = session.summary().context_revision;
    for i in 0..2 {
        session.advance(&ctx, i).unwrap();
    }
    let c = session.finish_to(w, &ctx).unwrap();
    (c.workspace, c.manifest_hash)
}
#[test]
fn boundary_wells_replicates_capture_restart_and_all_three_exports_preserve_exact_ids() {
    let (w, h) = fixture();
    let r = request(&w, &h, p::Plate::Wells96);
    let before = serde_json::to_value(&w).unwrap();
    let preview = p::preview(&w, &r).unwrap();
    preview
        .validate_binding(&p::preview_revision(&w, &r).unwrap())
        .unwrap();
    assert!(preview.resolved.iter().all(|r| r.missing_quantities));
    let c = p::capture(&w, &r).unwrap();
    c.validate_binding(&preview.context_revision).unwrap();
    assert_eq!(before, serde_json::to_value(&w).unwrap());
    assert_eq!(
        serde_json::to_value(&c.workspace.documents).unwrap(),
        before["documents"]
    );
    assert_eq!(
        serde_json::to_value(&c.workspace.history).unwrap(),
        before["history"]
    );
    assert_eq!(c.workspace.active_document, w.active_document);
    let restored =
        workspace::load_workspace(&workspace::save_workspace(&c.workspace).unwrap()).unwrap();
    let l = p::load(&restored, &c.attachment_hash).unwrap();
    assert_eq!(l.request.assignments[0].well.label(), "A1");
    assert_eq!(l.request.assignments[1].well.label(), "H12");
    assert_eq!(l.request.assignments.len(), 2);
    for format in [
        p::ExportFormat::Csv,
        p::ExportFormat::Html,
        p::ExportFormat::Pdf,
    ] {
        let out = p::export(&restored, &c.attachment_hash, format).unwrap();
        out.validate_binding(&p::export_revision(&restored, &c.attachment_hash, format).unwrap())
            .unwrap();
        match format {
            p::ExportFormat::Csv => {
                let text = String::from_utf8(out.bytes).unwrap();
                assert!(text.contains("\"alpha, \"\"one\"\" <review>\""));
                assert!(
                    text.contains(&h)
                        && text.contains("H12")
                        && text.contains("Missing; no quantities inferred")
                );
            }
            p::ExportFormat::Html => {
                let text = String::from_utf8(out.bytes).unwrap();
                assert!(text.contains("&lt;review&gt;") && !text.contains("<review>"));
                assert!(text.contains(&h) && text.contains(&l.context_revision));
            }
            p::ExportFormat::Pdf => {
                check_xref(&out.bytes);
                let text = extracted(&out.bytes);
                assert!(text.contains("PLATE LAYOUT"));
                assert!(!text.contains("SEQUENCE COMPARISON"));
                assert!(text.contains(&h) && text.contains(&l.context_revision));
                assert!(
                    text.contains("H12")
                        && text.contains("quantities")
                        && text.contains("Layout")
                        && text.contains("Context")
                );
            }
        }
    }
}
#[test]
fn three_eighty_four_last_well_and_empty_layout_are_explicit_not_automatic_fill() {
    let (w, h) = fixture();
    let mut r = request(&w, &h, p::Plate::Wells384);
    let l = p::preview(&w, &r).unwrap();
    assert_eq!(l.request.assignments[1].well.label(), "P24");
    r.assignments.clear();
    let c = p::capture(&w, &r).unwrap();
    assert!(c.layout.resolved.is_empty());
    assert!(
        p::load(&c.workspace, &c.attachment_hash)
            .unwrap()
            .request
            .assignments
            .is_empty()
    );
}
#[test]
fn wrong_declarations_construct_ids_and_released_equal_dna_event_aliases_refuse() {
    let (w, h) = fixture();
    let mut r = request(&w, &h, p::Plate::Wells96);
    if let p::Source::Planned {
        declaration_revision,
        ..
    } = &mut r.assignments[0].source
    {
        *declaration_revision = "0".repeat(64);
    }
    assert!(p::preview(&w, &r).is_err());
    let mut r = request(&w, &h, p::Plate::Wells96);
    if let p::Source::Planned { construct_id, .. } = &mut r.assignments[0].source {
        *construct_id = "beta".into();
    }
    assert!(p::preview(&w, &r).is_err());
    let (mut w, mh) = released(&w, &h);
    let manifest = release::load(&w, &mh).unwrap();
    let row = &manifest.rows[0];
    let product = row.output.construct.clone();
    let alias =
        DocumentHistory::new("equal-output-alias", w.revisions[&product.revision].clone()).unwrap();
    let alias_id = alias.id.clone();
    alias.contribute(&mut w).unwrap();
    let d = w.documents.iter().find(|d| d.id == alias_id).unwrap();
    let alias = EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    };
    assert_eq!(alias.revision, product.revision);
    let mut r = request(&w, &h, p::Plate::Wells96);
    r.assignments.truncate(1);
    r.assignments[0].source = p::Source::Released {
        manifest_hash: mh,
        construct_id: row.construct_id.clone(),
        run_id: row.output.run.run_id.clone(),
        product,
    };
    assert!(p::preview(&w, &r).is_ok());
    if let p::Source::Released { product, .. } = &mut r.assignments[0].source {
        *product = alias;
    }
    assert!(p::preview(&w, &r).is_err());
}
#[test]
fn released_assignment_closes_whole_indivisible_manifest_without_unrelated_workspace() {
    let (w, h) = fixture();
    let (w, mh) = released(&w, &h);
    let m = release::load(&w, &mh).unwrap();
    let mut r = request(&w, &h, p::Plate::Wells96);
    r.assignments.truncate(1);
    r.assignments[0].source = p::Source::Released {
        manifest_hash: mh.clone(),
        construct_id: m.rows[0].construct_id.clone(),
        run_id: m.rows[0].output.run.run_id.clone(),
        product: m.rows[0].output.construct.clone(),
    };
    let c = p::capture(&w, &r).unwrap();
    let deps = p::source_bindings(&c.workspace, &c.layout).unwrap();
    assert!(deps.attachments.contains(&mh) && deps.attachments.contains(&h));
    for row in &m.rows {
        assert!(deps.runs.contains(&row.output.run.run_id));
        assert!(deps.events.contains(&row.output.construct));
        assert!(deps.attachments.contains(&row.output.run.attachment_hash));
    }
    assert!(
        deps.events.contains(&event(&w, 2)),
        "Unselected choice is actual content of the indivisible design archive"
    );
    assert!(!deps.events.contains(&event(&w, 3)));
    assert_eq!(
        c.layout.resolved[0].product.as_ref(),
        Some(&m.rows[0].output.construct)
    );
}
#[test]
fn duplicate_coordinates_ids_label_bounds_and_oversize_decode_do_not_capture() {
    let (w, h) = fixture();
    let r = request(&w, &h, p::Plate::Wells96);
    let before = run_lineage::workspace_fingerprint(&w).unwrap();
    let mut duplicate = r.clone();
    duplicate.assignments[1].well = duplicate.assignments[0].well;
    assert!(p::capture(&w, &duplicate).is_err());
    let mut duplicate = r.clone();
    duplicate.assignments[1].id = duplicate.assignments[0].id.clone();
    assert!(p::capture(&w, &duplicate).is_err());
    let mut out = r.clone();
    out.assignments[0].well.row = 8;
    assert!(p::preview(&w, &out).is_err());
    let mut out = r.clone();
    out.assignments[0].label = "x".repeat(65);
    assert!(p::preview(&w, &out).is_err());
    let mut oversized = r.clone();
    oversized.assignments = vec![r.assignments[0].clone(); 385];
    assert!(p::preview(&w, &oversized).is_err());
    assert!(
        serde_json::from_slice::<p::Request>(&serde_json::to_vec(&oversized).unwrap()).is_err()
    );
    assert!(p::load_record(&vec![b' '; p::MAX_BYTES + 1]).is_err());
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), before);
}
#[test]
fn stale_preview_cannot_be_new_capture_and_archive_layout_survives_later_heads() {
    let (w, h) = fixture();
    let r = request(&w, &h, p::Plate::Wells96);
    let c = p::capture(&w, &r).unwrap();
    let mut docs = DocumentHistory::restore_all(&c.workspace).unwrap();
    docs[0]
        .commit(
            Molecule::dna("Later stock", "TCCC", Topology::Linear).unwrap(),
            "Later edit",
        )
        .unwrap();
    let mut later = c.workspace.clone();
    later.documents.clear();
    later.history.clear();
    for doc in docs {
        doc.contribute(&mut later).unwrap();
    }
    later.validate().unwrap();
    assert!(p::load(&later, &c.attachment_hash).is_ok());
    let next = p::preview_revision(&later, &r).unwrap();
    assert_ne!(next, c.context_revision);
    assert!(c.layout.validate_binding(&next).is_err());
    let mut forged = c.clone();
    forged.workspace = later;
    assert!(forged.validate_binding(&c.context_revision).is_err());
    let mut output = p::export(&c.workspace, &c.attachment_hash, p::ExportFormat::Csv).unwrap();
    output.bytes[0] ^= 1;
    assert!(output.validate_binding(&output.context_revision).is_err());
}

#[test]
fn native_layout_count_limit_is_enforced_without_replacing_old_captures() {
    let (mut w, h) = fixture();
    for i in 0..p::MAX_LAYOUTS {
        let mut r = request(&w, &h, p::Plate::Wells96);
        r.name = format!("Captured plate {i}");
        w = p::capture(&w, &r).unwrap().workspace;
    }
    assert_eq!(p::list(&w).unwrap().len(), 128);
    let fingerprint = run_lineage::workspace_fingerprint(&w).unwrap();
    let r = request(&w, &h, p::Plate::Wells96);
    assert!(p::capture(&w, &r).is_err());
    assert_eq!(run_lineage::workspace_fingerprint(&w).unwrap(), fingerprint);
}
#[test]
fn decoded_rows_and_resolved_receipts_cannot_change_native_bindings() {
    let (w, h) = fixture();
    let r = request(&w, &h, p::Plate::Wells96);
    let layout = p::preview(&w, &r).unwrap();
    let mut changed = layout.clone();
    changed.resolved[0].name = "Unrelated name assertion".into();
    assert!(changed.validate_for_workspace(&w).is_err());
    let mut raw = serde_json::to_value(&layout).unwrap();
    raw["resolved"] = serde_json::Value::Array(vec![raw["resolved"][0].clone(); 385]);
    assert!(serde_json::from_value::<p::Layout>(raw).is_err());
}

#[test]
fn planning_an_invalid_recorded_construct_does_not_approve_or_release_it() {
    let (w, h) = fixture();
    let mut request = batch::load(&w, &h).unwrap().request;
    let batch::RowRequest::Ready(row) = &mut request.rows[1] else {
        panic!("ready declaration")
    };
    row.target = TargetDeclaration::Existing(event(&w, 3));
    let mut session = batch::Session::begin(w.clone(), request).unwrap();
    let ctx = session.summary().context_revision;
    for i in 0..2 {
        session.advance(&ctx, i).unwrap();
    }
    let saved = session.into_capture(&ctx).unwrap();
    let before = serde_json::to_value(&saved.workspace.documents).unwrap();
    let request = p::Request {
        name: "Documentary invalid row".into(),
        plate: p::Plate::Wells96,
        assignments: vec![p::Assignment {
            id: "unapproved".into(),
            well: p::Well { row: 0, column: 0 },
            label: String::new(),
            source: planned(&saved.workspace, &saved.attachment_hash, "beta"),
        }],
    };
    let c = p::capture(&saved.workspace, &request).unwrap();
    assert_eq!(
        c.layout.resolved[0].recorded_status,
        Some(collagene_core::verification::Status::Invalid)
    );
    assert!(c.layout.resolved[0].run_id.is_none() && c.layout.resolved[0].product.is_none());
    assert_eq!(
        serde_json::to_value(&c.workspace.documents).unwrap(),
        before
    );
    assert_eq!(release::list(&c.workspace).unwrap().len(), 0);
}
