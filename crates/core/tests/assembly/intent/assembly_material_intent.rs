//! Literal source-strand intent is independent of the later physical recipe.
use collagene_core::{
    chemistry::{AnnealOptions, Strand, Terminus, tracked_duplex as modified_duplex},
    cloning::{
        assembly::{
            edit_intent as edit,
            intent::{self, Contribution, ContributionPosition, MaterialBinding, Selector},
            junction_edit as junction, substitution_intent as substitution,
        },
        restriction::catalog as catalog_digest,
    },
    documents::{
        history::DocumentHistory,
        workspace::{self, Workspace},
    },
    material::{
        derivation,
        lineage::{EventBinding, InputEventBinding},
        outputs as cloning_outputs,
        run::{self as run_material, MaterialOrigin},
    },
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    verification::{self, evidence::Status},
};
use collagene_model::{Feature, Location, Molecule, Sequence, Topology};
const TOP: &str = "AATTGCTACGGAATTCTTACG";
const BOTTOM: &str = "AATTCGTAAGAATTCCGTAGC";
const ACTUAL: &str = "AATTCGTAAG";
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn annealed(
    top: &str,
    bottom: &str,
    min_pair: usize,
    capped: bool,
) -> (Workspace, EventBinding, MaterialOrigin) {
    let mut t = Molecule::dna("Explicit top oligo", top, Topology::Linear).unwrap();
    if top == TOP {
        t.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::range(15, 21),
            qualifiers: vec![("label".into(), Some("Mapped source top".into()))],
            needs_review: false,
        });
    }
    let b = Molecule::dna("Explicit bottom oligo", bottom, Topology::Linear).unwrap();
    let target = Molecule::dna("Independent annealing top", top, Topology::Linear).unwrap();
    let plan = DesignPlan {
        name: "Native oligo annealing".into(),
        target_revision: target.revision(),
        inputs: vec![
            PlanInput {
                id: "top".into(),
                revision: t.revision(),
                molecule: t.clone(),
            },
            PlanInput {
                id: "bottom".into(),
                revision: b.revision(),
                molecule: b.clone(),
            },
        ],
        steps: vec![Operation::AnnealOligos {
            top: Source::Input("top".into()),
            bottom: Source::Input("bottom".into()),
            options: AnnealOptions {
                min_paired_bases: min_pair,
                max_candidates: if capped { 1 } else { 100 },
                top_chemistry: Terminus::CLEAVED,
                bottom_chemistry: Terminus::CLEAVED,
                ..Default::default()
            },
            candidate_index: 0,
        }],
    };
    let w = Workspace::from_documents(
        "Original oligos and independent target",
        vec![t, b, target],
        2,
    )
    .unwrap();
    let commit = cloning_outputs::finalize(
        &w,
        &plan,
        event(&w, 2),
        vec![
            InputEventBinding {
                input_id: "top".into(),
                event: event(&w, 0),
            },
            InputEventBinding {
                input_id: "bottom".into(),
                event: event(&w, 1),
            },
        ],
        &cloning_outputs::Options {
            name: "Saved actual duplex".into(),
            allow_incomplete: capped,
            save_primers: false,
            show_prepared: false,
        },
    )
    .unwrap();
    let pin = commit.receipt.construct;
    let origin = MaterialOrigin::Run(
        collagene_core::material::capsule::capsule(
            &commit.workspace,
            &run_material::Selection {
                run_id: commit.receipt.run.run_id,
                step_index: 0,
                event: pin.clone(),
            },
        )
        .unwrap(),
    );
    let w =
        workspace::load_workspace(&workspace::save_workspace(&commit.workspace).unwrap()).unwrap();
    (w, pin, origin)
}
fn source(w: &Workspace, pin: &EventBinding) -> Molecule {
    w.revisions[&pin.revision].clone()
}
fn inputs(w: &Workspace, pin: &EventBinding) -> Vec<PlanInput> {
    let m = source(w, pin);
    vec![PlanInput {
        id: "saved".into(),
        revision: m.revision(),
        molecule: m,
    }]
}
fn request(
    w: &Workspace,
    pin: &EventBinding,
    origin: &MaterialOrigin,
    start: usize,
    end: usize,
    strand: Strand,
) -> intent::Request {
    let fresh = origin.replay(&source(w, pin)).unwrap();
    let state = modified_duplex::summary(fresh.tracked())
        .unwrap()
        .state_revision;
    intent::Request {
        name: "Independent actual-strand reference".into(),
        fragments: vec![Contribution {
            role: "insert".into(),
            input_id: "saved".into(),
            selector: Selector::MaterialStrandRange {
                strand,
                start,
                end,
                expected_state_revision: state,
            },
        }],
        circular: false,
        material_bindings: vec![MaterialBinding {
            input_id: "saved".into(),
            event: pin.clone(),
            origin: origin.clone(),
        }],
    }
}
fn eco() -> catalog_digest::CatalogSelection {
    let c = catalog_digest::builtin_catalog().unwrap();
    let e = c
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap();
    catalog_digest::select(&c, &[e.id.clone()]).unwrap()
}
#[test]
fn actual_bottom_calls_extrude_beyond_displayed_top_and_map_only_real_paired_coordinates() {
    let (w, pin, origin) = annealed(TOP, BOTTOM, 17, false);
    let ins = inputs(&w, &pin);
    let before = workspace::save_workspace(&w).unwrap();
    let r = request(&w, &pin, &origin, 0, 10, Strand::Bottom);
    let report = intent::build(&ins, &r, Some(&w)).unwrap();
    assert_eq!(report.target.sequence.as_str(), ACTUAL);
    assert_ne!(
        report.target.sequence.as_str(),
        Sequence::dna("AATTCTTACG")
            .unwrap()
            .reverse_complement()
            .unwrap()
            .as_str()
    );
    let p = &report.fragments[0];
    assert_eq!(p.original_sequence.as_str(), ACTUAL);
    assert_eq!(
        (
            p.source_domain.as_ref().unwrap().start,
            p.source_domain.as_ref().unwrap().end
        ),
        (4, 25)
    );
    assert_eq!(
        (
            p.top_mappings[0].source_start,
            p.top_mappings[0].source_end,
            p.top_mappings[0].destination_start,
            p.top_mappings[0].reverse
        ),
        (15, 21, 4, true)
    );
    assert_eq!(
        p.source_position(0).unwrap(),
        ContributionPosition::Material {
            strand: Strand::Bottom,
            offset: 0
        }
    );
    let feature = report
        .target
        .features
        .iter()
        .find(|f| {
            f.qualifiers
                .iter()
                .any(|(_, v)| v.as_deref() == Some("Mapped source top"))
        })
        .unwrap();
    assert_eq!(
        feature.location,
        Location::Complement(Box::new(Location::range(4, 10)))
    );
    assert_eq!(
        feature.location.extract(&report.target.sequence).unwrap(),
        "CTTACG"
    );
    let descriptors = intent::material_source_descriptors(&report.target).unwrap();
    assert_eq!(descriptors.len(), 1);
    assert_eq!(descriptors[0].event, pin);
    assert_eq!(descriptors[0].selection, origin.selection());
    let metadata = serde_json::to_string(&report.target.metadata).unwrap();
    assert!(!metadata.contains(BOTTOM));
    assert!(!metadata.contains("record_bytes"));
    assert!(collagene_core::chemistry::Duplex::from_molecule(&report.target).is_err());
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    report
        .validate_binding(&intent::intent_revision(&ins, &r, Some(&w)).unwrap())
        .unwrap();
}
#[test]
fn independent_actual_bottom_target_verifies_only_after_native_reconstruction_digest_and_reverse_ligation()
 {
    let (mut w, pin, origin) = annealed(TOP, BOTTOM, 17, false);
    let ins = inputs(&w, &pin);
    let mut r = request(&w, &pin, &origin, 0, 10, Strand::Bottom);
    r.circular = true;
    let target = intent::build(&ins, &r, Some(&w)).unwrap().target;
    assert_eq!(target.sequence.as_str(), ACTUAL);
    let old_docs = serde_json::to_vec(&w.documents).unwrap();
    DocumentHistory::new("independent-bottom-target", target.clone())
        .unwrap()
        .contribute(&mut w)
        .unwrap();
    let target_pin = event(&w, w.documents.len() - 1);
    let MaterialOrigin::Run(capsule) = origin else {
        panic!("fixture run")
    };
    let plan = DesignPlan {
        name: "Actual sticky bottom assembly".into(),
        target_revision: target.revision(),
        inputs: ins,
        steps: vec![
            Operation::ReconstructRunProduct {
                template: Source::Input("saved".into()),
                origin: capsule,
            },
            Operation::CatalogDigest {
                template: Source::Step(0),
                selection: eco(),
                fragment_index: 1,
            },
            Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(1),
                    reverse: true,
                }],
                circular: true,
                allow_nicked: false,
                options: Default::default(),
            },
        ],
    };
    let v = verification::engine::verify(&plan, &target);
    assert_eq!(v.status, Status::ConsistentUnderAssumptions);
    assert_eq!(v.products.last().unwrap().sequence.as_str(), ACTUAL);
    let wrong =
        Molecule::dna("Wrong independent target", "AATTCGTAAA", Topology::Circular).unwrap();
    let mut wrong_plan = plan.clone();
    wrong_plan.target_revision = wrong.revision();
    assert_eq!(
        verification::engine::verify(&wrong_plan, &wrong).status,
        Status::Invalid
    );
    assert!(verification::output::export_protocol_html(&plan, &target).contains("AATTCGTAAG"));
    assert!(
        collagene_core::output::pdf::export_protocol_pdf(&plan, &target)
            .unwrap()
            .bytes
            .starts_with(b"%PDF")
    );
    let saved = cloning_outputs::finalize(
        &w,
        &plan,
        target_pin,
        vec![InputEventBinding {
            input_id: "saved".into(),
            event: pin,
        }],
        &cloning_outputs::Options {
            name: "Actual reverse circle".into(),
            save_primers: false,
            show_prepared: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    assert_eq!(
        serde_json::to_vec(&saved.workspace.documents[..w.documents.len() - 1]).unwrap(),
        old_docs
    );
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap()).unwrap();
    let selection = run_material::Selection {
        run_id: saved.receipt.run.run_id,
        step_index: 2,
        event: saved.receipt.construct,
    };
    let fresh = run_material::reconstruct(&reopened, &selection).unwrap();
    assert_eq!(fresh.tracked().duplex.top.as_str(), ACTUAL);
}
#[test]
fn missing_workspace_wrong_event_state_or_strand_bounds_cannot_substitute_equal_dna() {
    let (mut w, pin, origin) = annealed(TOP, BOTTOM, 17, false);
    let ins = inputs(&w, &pin);
    let r = request(&w, &pin, &origin, 0, 10, Strand::Bottom);
    assert!(intent::build(&ins, &r, None).is_err());
    let mut bad = r.clone();
    if let Selector::MaterialStrandRange {
        expected_state_revision,
        ..
    } = &mut bad.fragments[0].selector
    {
        *expected_state_revision = "0".repeat(64)
    }
    assert!(intent::build(&ins, &bad, Some(&w)).is_err());
    let mut bad = r.clone();
    if let Selector::MaterialStrandRange { end, .. } = &mut bad.fragments[0].selector {
        *end = 22
    }
    assert!(intent::build(&ins, &bad, Some(&w)).is_err());
    let m = source(&w, &pin);
    DocumentHistory::new("equal-dna-alias", m)
        .unwrap()
        .contribute(&mut w)
        .unwrap();
    let mut bad = r.clone();
    bad.material_bindings[0].event = event(&w, w.documents.len() - 1);
    assert!(intent::build(&ins, &bad, Some(&w)).is_err());
    assert!(r.fragments[0].reference_fragment().is_err());
}
#[test]
fn literal_trims_substitutions_and_indels_keep_original_actual_strand_receipt() {
    let (w, pin, origin) = annealed(TOP, BOTTOM, 17, false);
    let ins = inputs(&w, &pin);
    let r = request(&w, &pin, &origin, 0, 10, Strand::Bottom);
    let trims = vec![junction::Edit {
        junction: 0,
        remove_left: 1,
        remove_right: 0,
        inserted: "GG".into(),
    }];
    // For one linear contribution junction0 explicitly names its outgoing end;
    // it trims its suffix and appends the literal insert (no incoming row).
    let junction = junction::build_intent(&ins, &r, &trims, Some(&w)).unwrap();
    assert_eq!(junction.fragments[0].original_sequence.as_str(), ACTUAL);
    assert_eq!(junction.target.sequence.as_str(), "AATTCGTAAGG");
    assert_eq!(
        (
            junction.fragments[0].retained_start,
            junction.fragments[0].retained_end
        ),
        (0, 9)
    );
    let sr = substitution::Request {
        base: r,
        junction_edits: vec![],
        substitutions: vec![substitution::Substitution {
            role: "insert".into(),
            offset: 4,
            expected: "C".into(),
            replacement: "G".into(),
        }],
    };
    let sub = substitution::build(&ins, &sr, Some(&w)).unwrap();
    assert_eq!(sub.target.sequence.as_str(), "AATTGGTAAG");
    assert_eq!(
        sub.applied[0].source_position,
        ContributionPosition::Material {
            strand: Strand::Bottom,
            offset: 4
        }
    );
    let er = edit::Request {
        base: sr,
        edits: vec![edit::Edit {
            role: "insert".into(),
            start: 1,
            end: 3,
            expected: "AT".into(),
            replacement: "CCC".into(),
        }],
    };
    let edited = edit::build(&ins, &er, Some(&w)).unwrap();
    assert_eq!(edited.target.sequence.as_str(), "ACCCTGGTAAG");
    assert_eq!(edited.fragments[0].original_sequence.as_str(), ACTUAL);
    edited
        .validate_sources(&ins, &edited.context_revision, Some(&w))
        .unwrap();
    let mut changed = er;
    changed.edits[0].expected = "TT".into();
    assert!(edit::build(&ins, &changed, Some(&w)).is_err());
}
#[test]
fn derived_catalog_fragment_uses_actual_new_native_source_domain_after_restart() {
    let (w, pin, origin) = annealed(TOP, BOTTOM, 17, false);
    let requested = derivation::Request {
        sources: vec![collagene_core::material::derivation::BoundSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: pin,
            origin: derivation::SourceOrigin::Saved(origin),
        }],
        derivation: derivation::Derivation::CatalogDigest {
            selection: eco(),
            options: Default::default(),
            fragment_index: 1,
        },
        name: "Native derived fragment".into(),
        save_primers: false,
    };
    let installed = derivation::prepare_install(&w, &requested).unwrap();
    let pin = installed.receipt.product;
    let w = workspace::load_workspace(&workspace::save_workspace(&installed.workspace).unwrap())
        .unwrap();
    let origin = MaterialOrigin::Derived(derivation::capsule(&w, &pin).unwrap());
    let r = request(&w, &pin, &origin, 0, 10, Strand::Bottom);
    let report = intent::build(&inputs(&w, &pin), &r, Some(&w)).unwrap();
    assert_eq!(report.target.sequence.as_str(), ACTUAL);
    assert_eq!(report.fragments[0].source_domain.as_ref().unwrap().start, 4);
    assert_eq!(
        (
            report.fragments[0].top_mappings[0].source_start,
            report.fragments[0].top_mappings[0].source_end
        ),
        (4, 10)
    );
    assert_eq!(
        intent::material_source_descriptors(&report.target).unwrap()[0].selection,
        origin.selection()
    );
}
#[test]
fn unequal_thirteen_seventeen_strands_preserve_bottom_only_calls_without_widening_top() {
    let (w, pin, origin) = annealed("GCTACGGAATTCT", "AAAAAGAATTCCGTAGC", 13, false);
    let fresh = origin.replay(&source(&w, &pin)).unwrap();
    assert_eq!(
        (
            fresh.tracked().duplex.top.len(),
            fresh.tracked().duplex.bottom.len()
        ),
        (13, 17)
    );
    let r = request(&w, &pin, &origin, 0, 17, Strand::Bottom);
    let result = intent::build(&inputs(&w, &pin), &r, Some(&w)).unwrap();
    assert_eq!(result.target.sequence.as_str(), "AAAAAGAATTCCGTAGC");
    assert_eq!(
        (
            result.fragments[0].top_mappings[0].source_start,
            result.fragments[0].top_mappings[0].source_end,
            result.fragments[0].top_mappings[0].destination_start
        ),
        (0, 13, 4)
    );
    let mut bad = r;
    if let Selector::MaterialStrandRange { strand, .. } = &mut bad.fragments[0].selector {
        *strand = Strand::Top
    }
    assert!(intent::build(&inputs(&w, &pin), &bad, Some(&w)).is_err());
}
#[test]
fn incomplete_native_origin_is_preserved_as_source_evidence_without_inventing_isolation() {
    // The real17bp overlap is retained before the candidate cap stops searching
    // the final allowed offset (min16), making incompleteness native, not metadata.
    let (w, pin, origin) = annealed(TOP, BOTTOM, 16, true);
    assert!(origin.replay(&source(&w, &pin)).unwrap().incomplete());
    let r = request(&w, &pin, &origin, 0, 10, Strand::Bottom);
    let report = intent::build(&inputs(&w, &pin), &r, Some(&w)).unwrap();
    assert_eq!(report.target.sequence.as_str(), ACTUAL);
    assert!(report.fragments[0].incomplete);
    assert!(intent::material_source_descriptors(&report.target).unwrap()[0].incomplete);
    assert!(report.warnings.iter().any(|w| w.contains("incomplete")));
}
#[test]
fn digest_selector_uses_stored_bottom_offsets_and_only_actual_circles_wrap() {
    let domain = intent::MaterialDomain {
        top_length: 21,
        bottom_start: 4,
        bottom_length: 21,
        topology: Topology::Linear,
        state_revision: "a".repeat(64),
    };
    let fragment = intent::DigestGeometry {
        top_start: 11,
        top_length: 10,
        bottom_start: 15,
        bottom_length: 10,
    };
    let bottom = intent::digest_material_selector(&domain, &fragment, true).unwrap();
    assert!(matches!(
        bottom,
        Selector::MaterialStrandRange {
            strand: Strand::Bottom,
            start: 0,
            end: 10,
            ..
        }
    ));
    let forward = intent::digest_material_selector(&domain, &fragment, false).unwrap();
    assert!(matches!(
        forward,
        Selector::MaterialStrandRange {
            strand: Strand::Top,
            start: 11,
            end: 21,
            ..
        }
    ));
    let mut bad = fragment;
    bad.bottom_start = 20;
    assert!(intent::digest_material_selector(&domain, &bad, true).is_err());
    let circle = intent::MaterialDomain {
        top_length: 10,
        bottom_start: 0,
        bottom_length: 10,
        topology: Topology::Circular,
        state_revision: "b".repeat(64),
    };
    let wrapped = intent::DigestGeometry {
        top_start: 8,
        top_length: 5,
        bottom_start: 8,
        bottom_length: 5,
    };
    assert!(intent::digest_material_selector(&circle, &wrapped, false).is_ok());
    let mut linear = circle;
    linear.topology = Topology::Linear;
    assert!(intent::digest_material_selector(&linear, &wrapped, false).is_err());
}
