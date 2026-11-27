//! Literal TA/GC source alternatives over the sole tailing recipe compiler.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_intent::Selector,
    assembly_recipe::{Preparation, RunPreparation},
    chemistry::{EndTreatment, Strand},
    cloning_outputs, combinatorial_batch as batch,
    combinatorial_design::{
        self as design, Bin, Choice, ChoiceSelection, Design, LibraryMaterialBinding, Method,
        SelectedConstruct, TailingAction, TargetDeclaration,
    },
    combinatorial_release as release, material_derivation as derived,
    run_lineage::{self, EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    tailing::TailingBase,
    tailing_recipe::{Action, Profile, ProfileStatus},
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
    workspace::{self, Workspace},
};
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn input(m: Molecule) -> PlanInput {
    PlanInput {
        id: "stock".into(),
        revision: m.revision(),
        molecule: m,
    }
}
fn pcr(f: &str, r: &str) -> Preparation {
    Preparation::ExactPcr {
        forward: Primer::new("F", f, "").unwrap(),
        reverse: Primer::new("R", r, "").unwrap(),
        template_strand: Strand::Top,
    }
}
fn method(profile: Profile, left: Action, right: Action) -> Method {
    Method::Tailing {
        profile,
        actions: vec![
            TailingAction {
                bin_id: "left".into(),
                action: left,
            },
            TailingAction {
                bin_id: "right".into(),
                action: right,
            },
        ],
        allow_nicked: false,
    }
}
fn fixture(profile: Profile) -> (Workspace, Design, SelectedConstruct) {
    let target = match profile {
        Profile::Ta => "AAACCCTGGGTTTA",
        Profile::Gc => "AAACCCCGGGTTTG",
    };
    let w = Workspace::from_documents(
        "Explicit tailed library",
        vec![
            Molecule::dna("Left stock", "AAACCC", Topology::Linear).unwrap(),
            Molecule::dna("Right stock", "GGGTTT", Topology::Linear).unwrap(),
            Molecule::dna("Independent scarred circle", target, Topology::Circular).unwrap(),
            Molecule::dna(
                "Wrong independent circle",
                "AAACCCTGGGTTTC",
                Topology::Circular,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let (l, r) = match profile {
        Profile::Ta => (TailingBase::T, TailingBase::A),
        Profile::Gc => (TailingBase::C, TailingBase::G),
    };
    let d = Design {
        id: "tail-library".into(),
        name: "Explicit complementary tails".into(),
        circular: true,
        material_bindings: vec![],
        method: method(
            profile,
            Action::AddTail { base: l },
            Action::AddTail { base: r },
        ),
        bins: (0..2)
            .map(|i| Bin {
                id: if i == 0 { "left" } else { "right" }.into(),
                label: format!("Contribution {i}"),
                choices: vec![Choice {
                    id: format!("choice-{i}"),
                    label: format!("Stock {i}"),
                    source: event(&w, i),
                    selector: Selector::ReferenceRange {
                        start: 0,
                        end: 6,
                        reverse: false,
                    },
                    preparation_reverse: false,
                    preparation: Some(if i == 0 {
                        pcr("AAA", "GGG")
                    } else {
                        pcr("GGG", "AAA")
                    }),
                    treatment: Some(EndTreatment::Phosphorylate),
                    generation: None,
                }],
            })
            .collect(),
    };
    let c = SelectedConstruct {
        id: "intended".into(),
        name: "Literal scarred circle".into(),
        choices: vec![
            ChoiceSelection {
                bin_id: "left".into(),
                choice_id: "choice-0".into(),
            },
            ChoiceSelection {
                bin_id: "right".into(),
                choice_id: "choice-1".into(),
            },
        ],
        target: TargetDeclaration::Existing(event(&w, 2)),
        overrides: None,
    };
    (w, d, c)
}
fn circle_eq(a: &str, b: &str) -> bool {
    a.len() == b.len() && format!("{a}{a}").contains(b)
}
#[test]
fn original_ta_and_gc_preserve_literal_scars_and_finite_outcomes() {
    for profile in [Profile::Ta, Profile::Gc] {
        let (w, d, c) = fixture(profile);
        let before = serde_json::to_vec(&w).unwrap();
        let e = design::evaluate(&d, &c, &w).unwrap();
        let design::Preview::Tailing(p) = &e.preview else {
            panic!("sole tailing preview")
        };
        assert_eq!(
            p.profile_status,
            ProfileStatus::Matched,
            "{:?}",
            p.profile_findings
        );
        assert_eq!(
            e.verification().status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            e.verification().findings
        );
        let expected = match profile {
            Profile::Ta => "AAACCCTGGGTTTA",
            Profile::Gc => "AAACCCCGGGTTTG",
        };
        assert!(circle_eq(
            e.verification().products.last().unwrap().sequence.as_str(),
            expected
        ));
        assert_eq!(
            e.plan()
                .steps
                .iter()
                .filter(|op| matches!(op, Operation::Tail3Prime { .. }))
                .count(),
            2
        );
        assert!(matches!(
            e.plan().steps.last(),
            Some(Operation::Ligate { .. })
        ));
        let outcomes = e
            .verification()
            .physical_processing_reports
            .iter()
            .find_map(|step| {
                if let verification::PhysicalProcessingEvidence::LigationOutcomes(evidence) =
                    &step.evidence
                {
                    Some(&evidence.report)
                } else {
                    None
                }
            })
            .expect("actual finite ligation paths retained");
        assert!(outcomes.paths_checked > 0);
        assert_eq!(outcomes.request.circular, true);
        outcomes
            .validate_binding(&outcomes.context_revision)
            .unwrap();
        assert_eq!(before, serde_json::to_vec(&w).unwrap());
    }
}
#[test]
fn explicit_actions_are_complete_unique_profile_bound_and_bin_bound() {
    let (_, d, _) = fixture(Profile::Ta);
    let Method::Tailing { actions, .. } = &d.method else {
        unreachable!()
    };
    for bad_actions in [
        vec![],
        vec![actions[0].clone(), actions[0].clone()],
        vec![
            actions[0].clone(),
            TailingAction {
                bin_id: "foreign".into(),
                action: Action::KeepExisting,
            },
        ],
        vec![
            actions[0].clone(),
            TailingAction {
                bin_id: "right".into(),
                action: Action::AddTail {
                    base: TailingBase::C,
                },
            },
        ],
    ] {
        let mut bad = d.clone();
        bad.method = Method::Tailing {
            profile: Profile::Ta,
            actions: bad_actions,
            allow_nicked: false,
        };
        assert!(design::validate_design(&bad).is_err());
    }
}
fn saved_tailed_vector(w: Workspace) -> (Workspace, LibraryMaterialBinding) {
    let original = w.revisions[&event(&w, 0).revision].clone();
    let target =
        Molecule::dna("Independent historical T tail", "AAACCCT", Topology::Linear).unwrap();
    let mut w = w;
    collagene_core::history::DocumentHistory::new("historical-tail-target", target.clone())
        .unwrap()
        .contribute(&mut w)
        .unwrap();
    let plan = DesignPlan {
        name: "Actual saved T-tailed vector".into(),
        target_revision: target.revision(),
        inputs: vec![input(original)],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("stock".into()),
                forward: Primer::new("F", "AAA", "").unwrap(),
                reverse: Primer::new("R", "GGG", "").unwrap(),
                template_strand: Strand::Top,
            },
            Operation::Tail3Prime {
                template: Source::Step(0),
                base: TailingBase::T,
            },
            Operation::TreatEnds {
                template: Source::Step(1),
                treatment: EndTreatment::Phosphorylate,
            },
        ],
    };
    let commit = run_lineage::record_run(
        &w,
        &plan,
        event(&w, w.documents.len() - 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        false,
    )
    .unwrap();
    let p = &commit.receipt.products[2];
    let ev = EventBinding {
        document_id: p.document_id.clone(),
        event_id: p.event_id.clone(),
        revision: p.revision.clone(),
    };
    let origin = MaterialOrigin::Run(
        run_material::capsule(
            &commit.workspace,
            &run_material::Selection {
                run_id: commit.receipt.run_id,
                step_index: 2,
                event: ev.clone(),
            },
        )
        .unwrap(),
    );
    (
        commit.workspace,
        LibraryMaterialBinding { event: ev, origin },
    )
}
#[test]
fn saved_run_keep_and_derived_add_use_exact_owned_material_after_archive_restart() {
    let (w, mut d, c) = fixture(Profile::Ta);
    let (w, vector) = saved_tailed_vector(w);
    let right = event(&w, 1);
    let stock = &w.revisions[&right.revision];
    let f = Primer::new("Derived F", "GGG", "").unwrap();
    let r = Primer::new("Derived R", "AAA", "").unwrap();
    let docs = verification::simulate_pcr_documents(stock, &f, &r, Strand::Top).unwrap();
    assert_eq!(docs.len(), 1);
    let prepared = derived::prepare_install(
        &w,
        &derived::Request {
            sources: vec![collagene_core::material_derivation::BoundSource {
                role: collagene_core::material_derivation::SourceRole::Template,
                event: right,
                origin: derived::SourceOrigin::OriginalReference,
            }],
            derivation: derived::Derivation::ExactPcr {
                forward: f,
                reverse: r,
                template_strand: Strand::Top,
                selection: derived::PcrSelection {
                    selected_index: 0,
                    expected_preview_revision: docs[0].molecule.revision(),
                    geometry: docs[0].geometry.clone(),
                },
            },
            name: "Actual reusable insert".into(),
            save_primers: true,
        },
    )
    .unwrap();
    let right = prepared.receipt.product.clone();
    let right_origin = derived::capsule(&prepared.workspace, &right).unwrap();
    let w = workspace::load_workspace(&workspace::save_workspace(&prepared.workspace).unwrap())
        .unwrap();
    d.bins[0].choices[0].source = vector.event.clone();
    d.bins[0].choices[0].preparation = Some(Preparation::FromRun {
        origin: match &vector.origin {
            MaterialOrigin::Run(o) => o.clone(),
            _ => unreachable!(),
        },
        preparation: RunPreparation::WholeLinear,
    });
    d.bins[0].choices[0].treatment = None;
    d.bins[1].choices[0].source = right.clone();
    d.bins[1].choices[0].preparation = Some(Preparation::FromDerived {
        origin: right_origin.clone(),
        preparation: RunPreparation::WholeLinear,
    });
    d.material_bindings = vec![
        vector,
        LibraryMaterialBinding {
            event: right,
            origin: MaterialOrigin::Derived(right_origin),
        },
    ];
    d.method = method(
        Profile::Ta,
        Action::KeepExisting,
        Action::AddTail {
            base: TailingBase::A,
        },
    );
    let e = design::evaluate(&d, &c, &w).unwrap();
    assert_eq!(
        e.verification().status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        e.verification().findings
    );
    let design::Preview::Tailing(p) = &e.preview else {
        unreachable!()
    };
    assert_eq!(p.profile_status, ProfileStatus::Matched);
    assert!(p.compiled.rows[0].tail_step.is_none());
    assert!(p.compiled.rows[1].tail_step.is_some());
    assert!(circle_eq(
        e.verification().products.last().unwrap().sequence.as_str(),
        "AAACCCTGGGTTTA"
    ));
    let mut alias_workspace = w.clone();
    let same = alias_workspace.revisions[&d.bins[0].choices[0].source.revision].clone();
    collagene_core::history::DocumentHistory::new("equal-dna-other-owner", same)
        .unwrap()
        .contribute(&mut alias_workspace)
        .unwrap();
    let mut forged = d.clone();
    forged.bins[0].choices[0].source = event(&alias_workspace, alias_workspace.documents.len() - 1);
    assert!(design::evaluate(&forged, &c, &alias_workspace).is_err());
}
#[test]
fn batch_wrong_target_is_local_and_release_remains_atomic_and_replayable() {
    let (w, d, c) = fixture(Profile::Ta);
    let mut wrong = c.clone();
    wrong.id = "wrong".into();
    wrong.target = TargetDeclaration::Existing(event(&w, 3));
    let request = batch::Request {
        design: Some(d),
        rows: vec![batch::RowRequest::Ready(c), batch::RowRequest::Ready(wrong)],
    };
    let mut session = batch::Session::begin(w.clone(), request).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    session.advance(&context, 1).unwrap();
    let record = session.finish_record().unwrap();
    assert_eq!(
        record.rows[0].status,
        Some(Status::ConsistentUnderAssumptions)
    );
    assert_eq!(record.rows[1].status, Some(Status::Invalid));
    assert_eq!(
        record.rows[0].tailing_profile.as_ref().unwrap().status,
        ProfileStatus::Matched
    );
    assert_eq!(
        record.rows[1].tailing_profile.as_ref().unwrap().status,
        ProfileStatus::Matched
    );
    let mut tampered = record.clone();
    tampered.rows[0].tailing_profile.as_mut().unwrap().profile = Profile::Gc;
    assert!(tampered.validate().is_err());
    let commit = session.into_capture(&context).unwrap();
    let w =
        workspace::load_workspace(&workspace::save_workspace(&commit.workspace).unwrap()).unwrap();
    assert_eq!(
        batch::replay(&w, &commit.attachment_hash, 0)
            .unwrap()
            .current
            .status,
        record.rows[0].status
    );
    let options = cloning_outputs::Options {
        name: "Released TA clone".into(),
        show_prepared: false,
        save_primers: true,
        allow_incomplete: false,
    };
    let mut release = release::ReleaseSession::begin(
        w.clone(),
        release::Request {
            batch_hash: commit.attachment_hash.clone(),
            rows: vec![
                release::Selection {
                    construct_id: "intended".into(),
                    options: options.clone(),
                },
                release::Selection {
                    construct_id: "wrong".into(),
                    options: options.clone(),
                },
            ],
        },
    )
    .unwrap();
    let context = release.summary().context_revision;
    release.advance(&context, 0).unwrap();
    assert!(release.advance(&context, 1).is_err());
    assert!(release.finish_to(&w, &context).is_err());
    let mut release = release::ReleaseSession::begin(
        w.clone(),
        release::Request {
            batch_hash: commit.attachment_hash,
            rows: vec![release::Selection {
                construct_id: "intended".into(),
                options,
            }],
        },
    )
    .unwrap();
    let context = release.summary().context_revision;
    release.advance(&context, 0).unwrap();
    let commit = release.finish_to(&w, &context).unwrap();
    commit.validate().unwrap();
    assert_eq!(
        run_lineage::list_runs(&commit.workspace).unwrap().len(),
        run_lineage::list_runs(&w).unwrap().len() + 1
    );
    let manifest = release::load(&commit.workspace, &commit.manifest_hash).unwrap();
    assert!(manifest.rows[0].output.primer_attachment.is_some());
    let final_event = &manifest.rows[0].output.construct;
    assert!(circle_eq(
        commit.workspace.revisions[&final_event.revision]
            .sequence
            .as_str(),
        "AAACCCTGGGTTTA"
    ));
}

#[test]
fn compatible_blunt_ligation_does_not_claim_ta_geometry() {
    let (mut w, mut d, mut c) = fixture(Profile::Ta);
    let target = Molecule::dna(
        "Independent blunt product",
        "AAACCCGGGTTT",
        Topology::Linear,
    )
    .unwrap();
    collagene_core::history::DocumentHistory::new("blunt-target", target)
        .unwrap()
        .contribute(&mut w)
        .unwrap();
    c.target = TargetDeclaration::Existing(event(&w, w.documents.len() - 1));
    d.circular = false;
    d.method = method(Profile::Ta, Action::KeepExisting, Action::KeepExisting);
    for bin in &mut d.bins {
        bin.choices[0].preparation = Some(Preparation::DeclaredDuplex {
            left: collagene_core::chemistry::Terminus::CLEAVED,
            right: collagene_core::chemistry::Terminus::CLEAVED,
        });
        bin.choices[0].treatment = None;
    }
    let e = design::evaluate(&d, &c, &w).unwrap();
    assert_eq!(e.verification().status, Status::ConsistentUnderAssumptions);
    let design::Preview::Tailing(p) = &e.preview else {
        unreachable!()
    };
    assert_eq!(p.profile_status, ProfileStatus::Incompatible);
    assert!(!p.profile_findings.is_empty());
    assert!(
        e.plan()
            .steps
            .iter()
            .all(|op| !matches!(op, Operation::Tail3Prime { .. }))
    );
}

#[test]
fn unscarred_independent_original_regions_are_not_replaced_by_tailed_product() {
    let (w, d, mut c) = fixture(Profile::Ta);
    c.target = TargetDeclaration::OriginalRegions {
        junction_edits: vec![],
        substitutions: vec![],
        edits: vec![],
    };
    let e = design::evaluate(&d, &c, &w).unwrap();
    assert_eq!(e.resolved.target.sequence.as_str(), "AAACCCGGGTTT");
    assert_eq!(e.verification().status, Status::Invalid);
    let design::Preview::Tailing(p) = &e.preview else {
        unreachable!()
    };
    assert_eq!(p.profile_status, ProfileStatus::Matched);
    assert!(circle_eq(
        e.verification().products.last().unwrap().sequence.as_str(),
        "AAACCCTGGGTTTA"
    ));
}
