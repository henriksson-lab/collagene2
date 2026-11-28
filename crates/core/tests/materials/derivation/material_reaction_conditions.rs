//! Target-free native digest choices survive archive/replay without creating authority.
use collagene_core::{
    chemistry::EndRepairMethod,
    cloning::restriction::{
        catalog,
        conditions::{ChosenFormulation, EvidenceStatus, ReactionConditions},
    },
    documents::workspace::{self, Workspace},
    material::{derivation as d, lineage::EventBinding, run::MaterialOrigin},
};
use collagene_model::{Feature, Location, Molecule, Topology};
fn original() -> Workspace {
    let mut molecule = Molecule::dna("Original stock", "AAGAATTCTT", Topology::Linear).unwrap();
    molecule.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(7, 10),
        qualifiers: vec![("label".into(), Some("Retained source annotation".into()))],
        needs_review: false,
    });
    Workspace::from_documents("Exact material conditions", vec![molecule], 0).unwrap()
}
fn event(w: &Workspace) -> EventBinding {
    let doc = &w.documents[0];
    EventBinding {
        document_id: doc.id.clone(),
        event_id: doc.head_event.clone().unwrap(),
        revision: doc.head.clone(),
    }
}
fn options(formulation: Option<&str>) -> ReactionConditions {
    ReactionConditions {
        schema_version: 1,
        formulations: vec![ChosenFormulation {
            enzyme_id: "ecori".into(),
            formulation_id: formulation.map(str::to_owned),
        }],
        buffer: "rCutSmart".into(),
        temperature_milli_c: 37000,
        incubation_minutes: 60,
        minimum_activity_percent: 75,
        activity_context: Some("NEB unit-assay relative activity".into()),
        temperature_context: Some("vendor recommended incubation".into()),
    }
}
fn request(w: &Workspace, conditions: Option<ReactionConditions>) -> d::Request {
    d::Request {
        sources: vec![d::BoundSource {
            role: d::SourceRole::Template,
            event: event(w),
            origin: d::SourceOrigin::OriginalReference,
        }],
        derivation: d::Derivation::CatalogDigest {
            selection: catalog::select(&catalog::builtin_catalog().unwrap(), &["ecori".into()])
                .unwrap(),
            conditions,
            options: Default::default(),
            fragment_index: 1,
        },
        name: "Literal selected right fragment".into(),
        save_primers: false,
    }
}
#[test]
fn documentary_constraints_survive_restart_and_do_not_withhold_actual_fragment_geometry() {
    let w = original();
    let before = workspace::save_workspace(&w).unwrap();
    let mut physical = None;
    for (choice, status) in [
        (None, None),
        (Some(options(None)), Some(EvidenceStatus::Unknown)),
        (
            Some(options(Some("neb-r0101"))),
            Some(EvidenceStatus::BelowThreshold),
        ),
        (
            Some(options(Some("neb-r3101"))),
            Some(EvidenceStatus::Known),
        ),
    ] {
        let request = request(&w, choice.clone());
        let installed = d::prepare_install(&w, &request).unwrap();
        let restarted =
            workspace::load_workspace(&workspace::save_workspace(&installed.workspace).unwrap())
                .unwrap();
        let fresh = d::reconstruct(&restarted, &installed.receipt.product).unwrap();
        assert_eq!(fresh.reference().sequence.as_str(), "AATTCTT");
        assert_eq!(
            fresh.reference().features[0].location.spans(),
            vec![(4, 7, false)]
        );
        let actual = serde_json::to_value(fresh.tracked()).unwrap();
        if let Some(expected) = &physical {
            assert_eq!(&actual, expected);
        } else {
            physical = Some(actual);
        }
        assert_eq!(
            fresh.reaction_conditions(),
            installed.receipt.reaction_conditions
        );
        assert_eq!(
            fresh.incomplete(),
            status.is_some_and(|s| s != EvidenceStatus::Known)
        );
        let record =
            d::load_record(&restarted.attachments[&installed.receipt.record_hash].bytes).unwrap();
        assert_eq!(record.schema_version, 5);
        let d::Derivation::CatalogDigest { conditions, .. } = record.request.derivation else {
            panic!("digest")
        };
        assert_eq!(conditions, choice);
        if let Some(status) = status {
            let evidence = &fresh.reaction_conditions()[0];
            assert_eq!(evidence.report.documentary_status, status);
            assert_eq!(evidence.source_revision, event(&w).revision);
            assert_eq!(evidence.report.options, choice.unwrap());
            evidence.validate_historical().unwrap();
            let mut forged = installed.receipt.clone();
            forged.reaction_conditions.clear();
            assert!(forged.validate_archive(&restarted).is_err());
        } else {
            assert!(fresh.reaction_conditions().is_empty());
        }
        let capsule = d::capsule(&restarted, &installed.receipt.product).unwrap();
        assert_eq!(capsule.schema_version, 5);
        let replay = d::replay(&capsule, fresh.reference()).unwrap();
        assert_eq!(replay.reaction_conditions(), fresh.reaction_conditions());
    }
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
#[test]
fn replayed_ancestry_keeps_fresh_documentary_evidence_through_end_repair() {
    let w = original();
    let installed = d::prepare_install(&w, &request(&w, Some(options(None)))).unwrap();
    let source = installed.receipt.product.clone();
    let origin = MaterialOrigin::Derived(d::capsule(&installed.workspace, &source).unwrap());
    let request = d::Request {
        sources: vec![d::BoundSource {
            role: d::SourceRole::Template,
            event: source,
            origin: d::SourceOrigin::Saved(origin.clone()),
        }],
        derivation: d::Derivation::EndRepair {
            method: EndRepairMethod::FillFivePrimeOverhangs,
        },
        name: "Filled conditional fragment".into(),
        save_primers: false,
    };
    let repaired = d::prepare_install(&installed.workspace, &request).unwrap();
    let fresh = d::reconstruct(&repaired.workspace, &repaired.receipt.product).unwrap();
    assert_eq!(fresh.reference().sequence.as_str(), "AATTCTT");
    assert!(fresh.incomplete());
    assert_eq!(
        fresh.reaction_conditions(),
        installed.receipt.reaction_conditions
    );
    let replay = origin
        .replay(&installed.workspace.revisions[&installed.receipt.product.revision])
        .unwrap();
    assert_eq!(replay.reaction_conditions(), fresh.reaction_conditions());
}
#[test]
fn current_nullable_schema_and_bounded_documentary_transport_are_strict() {
    let w = original();
    let request = request(&w, Some(options(Some("neb-r3101"))));
    let installed = d::prepare_install(&w, &request).unwrap();
    let mut json = serde_json::to_value(&request).unwrap();
    json["derivation"]["CatalogDigest"]
        .as_object_mut()
        .unwrap()
        .remove("conditions");
    assert!(serde_json::from_value::<d::Request>(json).is_err());
    let reports = &installed.receipt.reaction_conditions;
    d::reaction_conditions::validate(reports).unwrap();
    let mut duplicates = reports.to_vec();
    duplicates.push(reports[0].clone());
    assert!(d::reaction_conditions::validate(&duplicates).is_err());
    let oversized = vec![reports[0].clone(); d::reaction_conditions::MAX_CONTEXTS + 1];
    assert!(matches!(
        d::reaction_conditions::validate(&oversized),
        Err(collagene_model::Error::TooLarge)
    ));
    let mut old_magic = installed.workspace.attachments[&installed.receipt.record_hash]
        .bytes
        .clone();
    old_magic[7] = 2;
    assert!(d::load_record(&old_magic).is_err());
    let mut old_capsule = d::capsule(&installed.workspace, &installed.receipt.product).unwrap();
    old_capsule.schema_version = 2;
    assert!(old_capsule.validate_documentary().is_err());
    let mut mismatched = request.clone();
    let d::Derivation::CatalogDigest {
        conditions: Some(options),
        ..
    } = &mut mismatched.derivation
    else {
        unreachable!()
    };
    options.formulations[0].enzyme_id = "bamhi".into();
    assert!(d::prepare_install(&w, &mismatched).is_err());
}
