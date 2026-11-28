use collagene_core::{
    chemistry::EndTreatment,
    run_lineage::*,
    verification::*,
    workspace::{self, RevisionEvent, Workspace},
    *,
};
use collagene_model::{Molecule, Topology};
fn binding(w: &Workspace, n: usize) -> EventBinding {
    let d = &w.documents[n];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> (Workspace, DesignPlan) {
    let s = Molecule::dna("source", "GATATCAAAA", Topology::Circular).unwrap();
    let t = Molecule::dna("independent target", "GATATCAAAA", Topology::Circular).unwrap();
    let p = DesignPlan {
        name: "digest PNK ligation".into(),
        target_revision: t.revision(),
        inputs: vec![PlanInput {
            id: "s".into(),
            revision: s.revision(),
            molecule: s.clone(),
        }],
        steps: vec![
            Operation::Digest {
                template: Source::Input("s".into()),
                enzymes: vec!["EcoRV".into()],
                fragment_index: 0,
            },
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Phosphorylate,
            },
            Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                }],
                circular: true,
                allow_nicked: false,
                options: Default::default(),
            },
        ],
    };
    (
        Workspace::from_documents("test", vec![s.clone(), s, t], 2).unwrap(),
        p,
    )
}
fn capture(w: &Workspace, p: &DesignPlan) -> RunCommit {
    record_run(
        w,
        p,
        binding(w, 2),
        vec![InputEventBinding {
            input_id: "s".into(),
            event: binding(w, 1),
        }],
        false,
    )
    .unwrap()
}
fn identity(r: &RunRecord) -> String {
    let tuple = (
        "cloning-run-v1",
        r.schema_version,
        r.plan_schema,
        &r.plan_revision,
        &r.plan,
        &r.target,
        &r.inputs,
        r.products
            .iter()
            .map(|p| (&p.step_index, &p.revision))
            .collect::<Vec<_>>(),
        &r.assessment,
    );
    blake3::hash(&serde_json::to_vec(&tuple).unwrap())
        .to_hex()
        .to_string()
}
#[test]
fn retained_graph_uses_exact_events_keeps_sources_and_target_and_is_idempotent() {
    let (w, p) = fixture();
    assert_eq!(w.documents[0].head, w.documents[1].head);
    let before = serde_json::to_value(&w).unwrap();
    let c = capture(&w, &p);
    c.validate().unwrap();
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
    assert_eq!(c.workspace.active_document, w.active_document);
    assert_eq!(
        serde_json::to_value(&c.workspace.documents[..3]).unwrap(),
        serde_json::to_value(&w.documents).unwrap()
    );
    assert_eq!(c.receipt.products.len(), 3);
    let r = load_run(&c.workspace, &c.receipt.run_id).unwrap();
    for (i, b) in r.products.iter().enumerate() {
        let e = c
            .workspace
            .history
            .iter()
            .find(|e| e.id == b.event_id)
            .unwrap();
        assert_eq!(
            e.parents,
            vec![if i == 0 {
                binding(&w, 1).event_id
            } else {
                r.products[i - 1].event_id.clone()
            }]
        );
        assert!(!e.parents.contains(&binding(&w, 2).event_id));
    }
    let again = capture(&c.workspace, &p);
    assert!(again.receipt.already_recorded);
    assert_eq!(again.receipt.run_id, c.receipt.run_id);
    assert_eq!(
        serde_json::to_value(&again.workspace).unwrap(),
        serde_json::to_value(&c.workspace).unwrap()
    );
    let replay = replay_run(&c.workspace, &c.receipt.run_id).unwrap();
    assert!(replay.products_match && replay.status_matches && replay.engine_matches);
    replay.validate().unwrap();
    assert_eq!(list_runs(&c.workspace).unwrap().len(), 1);
}
#[test]
fn exact_same_revision_event_choice_changes_run_and_parent() {
    let (w, p) = fixture();
    let a = capture(&w, &p);
    let b = record_run(
        &w,
        &p,
        binding(&w, 2),
        vec![InputEventBinding {
            input_id: "s".into(),
            event: binding(&w, 0),
        }],
        false,
    )
    .unwrap();
    assert_ne!(a.receipt.run_id, b.receipt.run_id);
    assert_ne!(
        a.receipt.products[0].event_id,
        b.receipt.products[0].event_id
    );
}
#[test]
fn wrong_bindings_wrong_target_partial_and_forward_graphs_are_rejected() {
    let (w, p) = fixture();
    let mut wrong = binding(&w, 1);
    wrong.document_id = "document-1".into();
    assert!(
        record_run(
            &w,
            &p,
            binding(&w, 2),
            vec![InputEventBinding {
                input_id: "s".into(),
                event: wrong
            }],
            false
        )
        .is_err()
    );
    assert!(
        record_run(
            &w,
            &p,
            binding(&w, 1),
            vec![InputEventBinding {
                input_id: "s".into(),
                event: binding(&w, 1)
            }],
            true
        )
        .is_err()
    );
    let mut broken = p.clone();
    broken.steps.push(Operation::Digest {
        template: Source::Step(2),
        enzymes: vec!["does-not-exist".into()],
        fragment_index: 0,
    });
    assert!(
        record_run(
            &w,
            &broken,
            binding(&w, 2),
            vec![InputEventBinding {
                input_id: "s".into(),
                event: binding(&w, 1)
            }],
            true
        )
        .is_err()
    );
    broken.steps = vec![Operation::TreatEnds {
        template: Source::Step(0),
        treatment: EndTreatment::Phosphorylate,
    }];
    assert!(
        record_run(
            &w,
            &broken,
            binding(&w, 2),
            vec![InputEventBinding {
                input_id: "s".into(),
                event: binding(&w, 1)
            }],
            true
        )
        .is_err()
    );
}
#[test]
fn explicit_tail_operation_is_preserved_and_replayed() {
    let (_, mut p) = fixture();
    let target = Molecule::dna("tail target", "ATCAAAAGATA", Topology::Linear).unwrap();
    let w = Workspace::from_documents(
        "tail",
        vec![
            p.inputs[0].molecule.clone(),
            p.inputs[0].molecule.clone(),
            target.clone(),
        ],
        2,
    )
    .unwrap();
    p.target_revision = target.revision();
    p.steps.truncate(1);
    p.steps.push(Operation::Tail3Prime {
        template: Source::Step(0),
        base: collagene_core::tailing::TailingBase::A,
    });
    let c = capture(&w, &p);
    let r = replay_run(&c.workspace, &c.receipt.run_id).unwrap();
    assert!(r.products_match);
    assert_eq!(
        r.current_report.products.last().unwrap().sequence.as_str(),
        "ATCAAAAGATA"
    );
}
#[test]
fn source_arity_bound_is_explicit_and_input_is_unchanged() {
    let (w, mut p) = fixture();
    let before = workspace_fingerprint(&w).unwrap();
    p.steps = vec![Operation::Overlap {
        fragments: vec![Source::Input("s".into()); 129],
        overlaps: vec![0; 128],
        circular: false,
    }];
    assert!(
        record_run(
            &w,
            &p,
            binding(&w, 2),
            vec![InputEventBinding {
                input_id: "s".into(),
                event: binding(&w, 1)
            }],
            true
        )
        .is_err()
    );
    assert_eq!(workspace_fingerprint(&w).unwrap(), before);
}
#[test]
fn capacity_failure_is_transactional_and_ids_avoid_collisions() {
    let (w, p) = fixture();
    let c = capture(&w, &p);
    let mut collision = w.clone();
    let id = c.receipt.products[0].document_id.clone();
    let m = collision.head("document-1").unwrap().clone();
    let e = RevisionEvent::new(id.clone(), m.revision(), vec![], "unrelated").unwrap();
    let mut d = collision.documents[0].clone();
    d.id = id;
    d.head_event = Some(e.id.clone());
    collision.documents.push(d);
    collision.history.push(e);
    let c2 = capture(&collision, &p);
    assert_ne!(
        c2.receipt.products[0].document_id,
        c.receipt.products[0].document_id
    );
    let mut full = w.clone();
    for n in 3..workspace::MAX_DOCUMENTS {
        let id = format!("capacity-{n}");
        let mut d = full.documents[0].clone();
        let e = RevisionEvent::new(id.clone(), d.head.clone(), vec![], "import").unwrap();
        d.id = id;
        d.head_event = Some(e.id.clone());
        full.documents.push(d);
        full.history.push(e);
    }
    full.validate().unwrap();
    let before = workspace_fingerprint(&full).unwrap();
    assert!(
        record_run(
            &full,
            &p,
            binding(&full, 2),
            vec![InputEventBinding {
                input_id: "s".into(),
                event: binding(&full, 1)
            }],
            false
        )
        .is_err()
    );
    assert_eq!(workspace_fingerprint(&full).unwrap(), before);
}
#[test]
fn archive_roundtrip_structural_crossrefs_and_opaque_attachment_tamper() {
    let (w, p) = fixture();
    let c = capture(&w, &p);
    let bytes = workspace::save_workspace(&c.workspace).unwrap();
    let loaded = workspace::load_workspace(&bytes).unwrap();
    assert_eq!(list_runs(&loaded).unwrap()[0].run_id, c.receipt.run_id);
    let mut wrong = c.workspace.clone();
    let h = c.receipt.attachment_hash.clone();
    let mut attachment = wrong.attachments.remove(&h).unwrap();
    let mut raw: serde_json::Value = serde_json::from_slice(&attachment.bytes).unwrap();
    raw["products"][0]["event_id"] =
        "event-ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff".into();
    attachment.bytes = serde_json::to_vec(&raw).unwrap();
    let h = blake3::hash(&attachment.bytes).to_hex().to_string();
    wrong.attachments.insert(h, attachment);
    assert!(wrong.validate().is_err());
    let mut wrong = c.workspace.clone();
    wrong
        .attachments
        .get_mut(&c.receipt.attachment_hash)
        .unwrap()
        .bytes[0] = b'[';
    assert!(wrong.validate().is_err());
}
#[test]
fn record_unknown_nested_fields_are_rejected() {
    let (w, p) = fixture();
    let c = capture(&w, &p);
    let a = &c.workspace.attachments[&c.receipt.attachment_hash];
    let mut v: serde_json::Value = serde_json::from_slice(&a.bytes).unwrap();
    v["plan"]["inputs"][0]["injected"] = true.into();
    assert!(load_run_record(&serde_json::to_vec(&v).unwrap()).is_err());
}
#[test]
fn fingerprint_binds_navigation_and_attachment_metadata() {
    let (w, p) = fixture();
    let c = capture(&w, &p);
    let fp = workspace_fingerprint(&c.workspace).unwrap();
    let mut edited = c.workspace.clone();
    edited
        .attachments
        .get_mut(&c.receipt.attachment_hash)
        .unwrap()
        .name = "renamed.json".into();
    assert_ne!(workspace_fingerprint(&edited).unwrap(), fp);
    edited = c.workspace.clone();
    edited.active_document = "document-1".into();
    assert_ne!(workspace_fingerprint(&edited).unwrap(), fp);
    edited = c.workspace;
    edited.revisions.values_mut().next().unwrap().name = "forged".into();
    assert!(workspace_fingerprint(&edited).is_err());
}
#[test]
fn historical_engines_are_retained_but_current_replay_reports_the_difference() {
    let (w, p) = fixture();
    let mut c = capture(&w, &p);
    let mut a = c
        .workspace
        .attachments
        .remove(&c.receipt.attachment_hash)
        .unwrap();
    let mut r = load_run_record(&a.bytes).unwrap();
    let old_id = r.run_id.clone();
    let graph_hash = c
        .workspace
        .attachments
        .iter()
        .find_map(|(hash, attachment)| {
            if attachment.kind != physical_run_evidence::KIND {
                return None;
            }
            let graph: physical_run_evidence::Evidence =
                serde_json::from_slice(&attachment.bytes).unwrap();
            (graph.run_id == old_id).then(|| hash.clone())
        })
        .unwrap();
    let mut graph_attachment = c.workspace.attachments.remove(&graph_hash).unwrap();
    let mut graph: physical_run_evidence::Evidence =
        serde_json::from_slice(&graph_attachment.bytes).unwrap();
    r.assessment.engine_version = "0.0.1/science-previous".into();
    // Historical scientific stamps bind the nested report, its context, the
    // full-evidence product receipt, immutable revision and owning native events.
    for step in &mut graph.steps {
        let PhysicalProcessingEvidence::LigationOutcomes(evidence) = &mut step.evidence else {
            continue;
        };
        evidence.report.scientific_engine = r.assessment.engine_version.clone();
        evidence.report.context_revision = blake3::hash(
            &serde_json::to_vec(&(
                &evidence.report.engine_version,
                &evidence.report.scientific_engine,
                &evidence.report.material_engine,
                &evidence.report.sources,
                &evidence.report.request,
            ))
            .unwrap(),
        )
        .to_hex()
        .to_string();
        evidence.validate_documentary().unwrap();
        let product = r
            .products
            .iter_mut()
            .find(|p| p.step_index == step.step_index)
            .unwrap();
        let mut molecule = c.workspace.revisions[&product.revision].clone();
        molecule.metadata.insert(
            "source:ligation_outcomes:context".into(),
            evidence.report.context_revision.clone(),
        );
        molecule.metadata.insert(
            "source:ligation_outcomes:evidence_hash".into(),
            blake3::hash(&serde_json::to_vec(&*evidence).unwrap())
                .to_hex()
                .to_string(),
        );
        product.revision = molecule.revision();
        c.workspace
            .revisions
            .insert(product.revision.clone(), molecule);
    }
    r.run_id = identity(&r);
    assert_ne!(old_id, r.run_id);
    let mut previous = std::collections::BTreeMap::new();
    for (i, binding) in r.products.iter_mut().enumerate() {
        let event = c
            .workspace
            .history
            .iter_mut()
            .find(|e| e.id == binding.event_id)
            .unwrap();
        let old = event.id.clone();
        let parents = event
            .parents
            .iter()
            .map(|id| previous.get(id).cloned().unwrap_or_else(|| id.clone()))
            .collect();
        let replacement = RevisionEvent::new(
            event.document_id.clone(),
            binding.revision.clone(),
            parents,
            format!("Cloning run {} step {}", r.run_id, i + 1),
        )
        .unwrap();
        let document = c
            .workspace
            .documents
            .iter_mut()
            .find(|d| d.id == binding.document_id)
            .unwrap();
        document.head = binding.revision.clone();
        document.head_event = Some(replacement.id.clone());
        binding.event_id = replacement.id.clone();
        previous.insert(old, replacement.id.clone());
        *event = replacement;
    }
    graph.run_id = r.run_id.clone();
    graph.engine_version = r.assessment.engine_version.clone();
    graph.products = r.products.clone();
    graph_attachment.name = format!("{}.physical-processing-graphs.json", r.run_id);
    graph_attachment.bytes = serde_json::to_vec(&graph).unwrap();
    c.workspace.attachments.insert(
        blake3::hash(&graph_attachment.bytes).to_hex().to_string(),
        graph_attachment,
    );
    a.name = format!("{}.cloning-run.json", r.run_id);
    a.bytes = save_run_record(&r).unwrap();
    c.workspace
        .attachments
        .insert(blake3::hash(&a.bytes).to_hex().to_string(), a);
    c.workspace.validate().unwrap();
    let restarted =
        workspace::load_workspace(&workspace::save_workspace(&c.workspace).unwrap()).unwrap();
    let retained = physical_run_evidence::load(&restarted, &r.run_id)
        .unwrap()
        .unwrap();
    assert_eq!(retained.products, r.products);
    assert_eq!(retained.engine_version, r.assessment.engine_version);
    let replay = replay_run(&restarted, &r.run_id).unwrap();
    assert!(!replay.engine_matches);
    assert!(replay.status_matches);
    // Current native evidence carries its actual current scientific stamp. Equal
    // reference DNA cannot substitute for the historical evidence-bearing revision.
    assert!(!replay.products_match);
    assert_eq!(replay.current_report.products.len(), r.products.len());
    for (current, historical) in replay.current_report.products.iter().zip(&r.products) {
        assert_eq!(
            current.sequence,
            restarted.revisions[&historical.revision].sequence
        );
        assert_eq!(
            current.topology,
            restarted.revisions[&historical.revision].topology
        );
    }
    assert_eq!(replay.current_report.engine_version, ENGINE_VERSION);
    let last = r.products.last().unwrap();
    assert!(
        run_material::reconstruct(
            &restarted,
            &run_material::Selection {
                run_id: r.run_id,
                step_index: last.step_index,
                event: EventBinding {
                    document_id: last.document_id.clone(),
                    event_id: last.event_id.clone(),
                    revision: last.revision.clone()
                }
            }
        )
        .is_err()
    );
}

#[test]
fn matched_profile_limit_requires_explicit_incomplete_capture_and_keeps_both_roles() {
    let a = Molecule::dna("a", "AACGTCGATGCTAGCA", Topology::Linear).unwrap();
    let b = Molecule::dna("b", "ACGTCGATGCTAGCAT", Topology::Linear).unwrap();
    let target = Molecule::dna("target", "AACGTCGATGCTAGCAT", Topology::Linear).unwrap();
    let w = Workspace::from_documents(
        "limited profile",
        vec![a.clone(), b.clone(), target.clone()],
        2,
    )
    .unwrap();
    let p = DesignPlan {
        name: "profile repeat budget".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "a".into(),
                revision: a.revision(),
                molecule: a,
            },
            PlanInput {
                id: "b".into(),
                revision: b.revision(),
                molecule: b,
            },
        ],
        steps: vec![Operation::ProfiledOverlap {
            fragments: vec![
                OrientedSource {
                    source: Source::Input("a".into()),
                    reverse: false,
                },
                OrientedSource {
                    source: Source::Input("b".into()),
                    reverse: false,
                },
            ],
            overlaps: vec![15],
            circular: false,
            options: collagene_core::overlap_profiles::ProfileOptions {
                max_comparisons: 15,
                ..Default::default()
            },
        }],
    };
    let inputs = vec![
        InputEventBinding {
            input_id: "a".into(),
            event: binding(&w, 0),
        },
        InputEventBinding {
            input_id: "b".into(),
            event: binding(&w, 1),
        },
    ];
    let report = verify(&p, &target);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert_eq!(report.products.len(), 1);
    assert!(record_run(&w, &p, binding(&w, 2), inputs.clone(), false).is_err());
    let c = record_run(&w, &p, binding(&w, 2), inputs, true).unwrap();
    let e = c
        .workspace
        .history
        .iter()
        .find(|e| e.id == c.receipt.products[0].event_id)
        .unwrap();
    assert_eq!(
        e.parents,
        vec![binding(&w, 0).event_id, binding(&w, 1).event_id]
    );
    assert!(
        replay_run(&c.workspace, &c.receipt.run_id)
            .unwrap()
            .products_match
    );
}
#[test]
fn native_loaded_product_metadata_never_grants_physical_input_authority() {
    let (w, p) = fixture();
    let c = capture(&w, &p);
    let loaded =
        workspace::load_workspace(&workspace::save_workspace(&c.workspace).unwrap()).unwrap();
    let product = &c.receipt.products[0];
    let molecule = loaded.revisions[&product.revision].clone();
    let target = Molecule::dna(
        "independent imported intermediate",
        molecule.sequence.as_str(),
        Topology::Linear,
    )
    .unwrap();
    let isolated = Workspace::from_documents(
        "imported product",
        vec![molecule.clone(), target.clone()],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "metadata cannot authorize physical treatment".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "s".into(),
            revision: molecule.revision(),
            molecule,
        }],
        steps: vec![Operation::TreatEnds {
            template: Source::Input("s".into()),
            treatment: EndTreatment::Phosphorylate,
        }],
    };
    assert!(
        record_run(
            &isolated,
            &plan,
            binding(&isolated, 1),
            vec![InputEventBinding {
                input_id: "s".into(),
                event: binding(&isolated, 0)
            }],
            true
        )
        .is_err()
    );
}
#[test]
fn retained_historical_source_event_remains_exact_after_head_changes() {
    let (mut w, p) = fixture();
    let old = binding(&w, 1);
    let m = Molecule::dna("new source", "CCCCAAAA", Topology::Linear).unwrap();
    let event = RevisionEvent::new(
        old.document_id.clone(),
        m.revision(),
        vec![old.event_id.clone()],
        "edit after design",
    )
    .unwrap();
    let d = &mut w.documents[1];
    d.undo.push(d.head.clone());
    d.undo_events.push(old.event_id.clone());
    d.head = m.revision();
    d.head_event = Some(event.id.clone());
    w.revisions.insert(m.revision(), m);
    w.history.push(event);
    w.validate().unwrap();
    let c = record_run(
        &w,
        &p,
        binding(&w, 2),
        vec![InputEventBinding {
            input_id: "s".into(),
            event: old.clone(),
        }],
        false,
    )
    .unwrap();
    let e = c
        .workspace
        .history
        .iter()
        .find(|e| e.id == c.receipt.products[0].event_id)
        .unwrap();
    assert_eq!(e.parents, vec![old.event_id]);
    assert_eq!(
        serde_json::to_value(&w.documents[1]).unwrap(),
        serde_json::to_value(&c.workspace.documents[1]).unwrap()
    );
}
#[test]
fn native_outer_plan_schema_cannot_hide_retained_new_operations() {
    use std::io::{Cursor, Read, Write};
    let (w, p) = fixture();
    let c = capture(&w, &p);
    let bytes = workspace::save_workspace(&c.workspace).unwrap();
    let mut zip = zip::ZipArchive::new(Cursor::new(bytes)).unwrap();
    let mut out = zip::ZipWriter::new(Cursor::new(Vec::new()));
    for i in 0..zip.len() {
        let mut f = zip.by_index(i).unwrap();
        let name = f.name().to_string();
        let mut data = Vec::new();
        f.read_to_end(&mut data).unwrap();
        if name == "manifest.json" {
            let mut value: serde_json::Value = serde_json::from_slice(&data).unwrap();
            value["plan_schema"] = 2.into();
            data = serde_json::to_vec(&value).unwrap();
        }
        out.start_file(
            name,
            zip::write::SimpleFileOptions::default()
                .compression_method(zip::CompressionMethod::Stored),
        )
        .unwrap();
        out.write_all(&data).unwrap();
    }
    assert!(workspace::load_workspace(&out.finish().unwrap().into_inner()).is_err());
}

#[test]
fn replay_transport_validation_binds_stored_product_comparisons_and_workspace() {
    let (w, p) = fixture();
    let c = capture(&w, &p);
    let report = replay_run(&c.workspace, &c.receipt.run_id).unwrap();
    report.validate_for_workspace(&c.workspace).unwrap();
    let mut forged = report.clone();
    forged.products_match = false;
    assert!(forged.validate_for_workspace(&c.workspace).is_err());
    let mut changed = c.workspace.clone();
    changed.active_document = "document-1".into();
    assert!(report.validate_for_workspace(&changed).is_err());
}
