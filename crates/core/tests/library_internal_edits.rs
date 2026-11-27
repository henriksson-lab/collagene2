//! Independent row targets, generic native source ownership and durable release.
use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    assembly_edit_intent::Edit,
    assembly_edited_primer_design as primers,
    assembly_recipe::{Preparation, RunPreparation},
    assembly_substitution_intent::Substitution,
    cloning_outputs, combinatorial_batch as batch, combinatorial_design as library,
    combinatorial_release as release, material_derivation as derivation,
    overlap_extension_primer_design::Options,
    overlap_profiles::ProfileOptions,
    run_lineage::{self, EventBinding},
    run_material::MaterialOrigin,
    verification::{self, Status},
    workspace::{self, Workspace},
};
const VECTOR: &str = "CACGTTCGTAGGACCA";
const INSERT: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const EDITED: &str = "CACGTTCGTAGGACCAAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn options() -> Options {
    let mut o = Options {
        homology_flank: 8,
        ..Default::default()
    };
    o.primer_options.min_length = 8;
    o.primer_options.max_length = 8;
    o.primer_options.min_gc_fraction = 0.;
    o.primer_options.max_gc_fraction = 1.;
    o.primer_options.max_tm_difference_c = 100.;
    o.primer_options.max_homopolymer = 200;
    o.primer_options.max_three_prime_complementarity = 200;
    o
}
fn fixture() -> (Workspace, library::Design, library::SelectedConstruct) {
    let w = Workspace::from_documents(
        "Saved library source",
        vec![
            Molecule::dna("vector", VECTOR, Topology::Linear).unwrap(),
            Molecule::dna("original insert", INSERT, Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let forward = Primer::new("Original F", "AGTCCGAT", "").unwrap();
    let reverse = Primer::new("Original R", "TGACGATC", "").unwrap();
    let documents = verification::simulate_pcr_documents(
        &w.revisions[&event(&w, 1).revision],
        &forward,
        &reverse,
        collagene_core::chemistry::Strand::Top,
    )
    .unwrap();
    assert_eq!(documents.len(), 1);
    let request = derivation::Request {
        source: event(&w, 1),
        origin: derivation::SourceOrigin::OriginalReference,
        derivation: derivation::Derivation::ExactPcr {
            forward: forward.clone(),
            reverse: reverse.clone(),
            selection: derivation::PcrSelection {
                selected_index: 0,
                expected_preview_revision: documents[0].molecule.revision(),
                geometry: documents[0].geometry.clone(),
            },
            template_strand: collagene_core::chemistry::Strand::Top,
        },
        name: "Actual saved insert PCR".into(),
        save_primers: true,
    };
    let installation = derivation::prepare_install(&w, &request).unwrap();
    let product = installation.receipt.product;
    let w = workspace::load_workspace(&workspace::save_workspace(&installation.workspace).unwrap())
        .unwrap();
    let capsule = derivation::capsule(&w, &product).unwrap();
    let origin = MaterialOrigin::Derived(capsule.clone());
    let design = library::Design {
        material_bindings: vec![library::LibraryMaterialBinding {
            event: product.clone(),
            origin: origin.clone(),
        }],
        id: "internal-library".into(),
        name: "Independent row edits".into(),
        circular: false,
        method: library::Method::Overlap {
            finishing: collagene_core::assembly_recipe::OverlapFinishing::ReferenceOnly,
            overlaps: vec![8],
            options: ProfileOptions::default(),
        },
        bins: vec![
            library::Bin {
                id: "vector".into(),
                label: "Vector".into(),
                choices: vec![library::Choice {
                    id: "vector-choice".into(),
                    label: "Vector".into(),
                    source: event(&w, 0),
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: 0,
                        end: 16,
                        reverse: false,
                    },
                    preparation_reverse: false,
                    generation: None,
                    treatment: None,
                    preparation: Some(Preparation::ExactPcr {
                        forward: Primer::new("VF", "CACGTTCG", "").unwrap(),
                        reverse: Primer::new("VR", "TGGTCCTA", "").unwrap(),
                        template_strand: collagene_core::chemistry::Strand::Top,
                    }),
                }],
            },
            library::Bin {
                id: "insert".into(),
                label: "Insert".into(),
                choices: vec![library::Choice {
                    id: "insert-choice".into(),
                    label: "Saved insert".into(),
                    source: product,
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: 0,
                        end: 40,
                        reverse: false,
                    },
                    preparation_reverse: false,
                    generation: None,
                    treatment: None,
                    preparation: Some(Preparation::FromDerived {
                        origin: capsule,
                        preparation: RunPreparation::ExactPcr {
                            forward,
                            reverse,
                            template_strand: collagene_core::chemistry::Strand::Top,
                        },
                    }),
                }],
            },
        ],
    };
    let row = library::SelectedConstruct {
        id: "edited-row".into(),
        name: "Literal edited construct".into(),
        overrides: None,
        choices: vec![
            library::ChoiceSelection {
                bin_id: "vector".into(),
                choice_id: "vector-choice".into(),
            },
            library::ChoiceSelection {
                bin_id: "insert".into(),
                choice_id: "insert-choice".into(),
            },
        ],
        target: library::TargetDeclaration::OriginalRegions {
            junction_edits: vec![],
            substitutions: vec![],
            edits: vec![Edit {
                role: "fragment-insert".into(),
                start: 16,
                end: 20,
                expected: "GGAC".into(),
                replacement: "TT".into(),
            }],
        },
    };
    (w, design, row)
}
fn proposal(
    w: &Workspace,
    d: &library::Design,
    row: &library::SelectedConstruct,
) -> primers::Report {
    let resolved = library::resolve(d, row, w).unwrap();
    let library::IndependentIntent::Edits(intent) = resolved.intent.unwrap() else {
        panic!("actual general intent")
    };
    assert_eq!(
        intent.target.sequence.as_str(),
        EDITED,
        "independent literal before primer proposal"
    );
    let request = primers::Request {
        intent: intent.request,
        overlap: 8,
        options: options(),
        profile_options: ProfileOptions::default(),
        finishing: collagene_core::assembly_recipe::OverlapFinishing::ReferenceOnly,
    };
    let p = primers::design(&resolved.inputs, &request, Some(w)).unwrap();
    p.validate_sources(&resolved.inputs, &p.context_revision, Some(w))
        .unwrap();
    assert_eq!(
        p.preview.verification.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        p.preview.verification.findings
    );
    assert_eq!(
        p.preview
            .verification
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        EDITED
    );
    p
}
fn apply(row: &mut library::SelectedConstruct, p: &primers::Report) {
    row.overrides = Some(library::ConstructOverrides {
        preparations: p
            .recipe
            .fragments
            .iter()
            .map(|f| library::PreparationOverride {
                preparation_reverse: false,
                bin_id: if f.role == "fragment-vector" {
                    "vector"
                } else {
                    assert_eq!(f.role, "fragment-insert");
                    "insert"
                }
                .into(),
                preparation: Some(f.preparation.clone()),
                treatment: None,
                generation: f.generation.clone(),
            })
            .collect(),
        method: None,
        circular: None,
    });
}
#[test]
fn derived_source_row_proposal_batch_restart_release_and_native_replay_preserve_literal_target() {
    let (w, d, mut row) = fixture();
    let heads = serde_json::to_vec(&w.documents).unwrap();
    let history = serde_json::to_vec(&w.history).unwrap();
    let p = proposal(&w, &d, &row);
    assert!(matches!(
        &p.recipe.fragments[1].preparation,
        Preparation::FromDerived {
            preparation: RunPreparation::OverlapExtensionPcr { .. },
            ..
        }
    ));
    assert_eq!(
        p.rows[1].template_strand,
        collagene_core::chemistry::Strand::Top
    );
    assert_eq!(p.rows[1].search.len(), 3);
    apply(&mut row, &p);
    let mut stale = row.clone();
    stale.id = "stale-row".into();
    let library::TargetDeclaration::OriginalRegions { edits, .. } = &mut stale.target else {
        unreachable!()
    };
    edits[0].expected = "AAAA".into();
    let r = batch::Request {
        design: Some(d),
        rows: vec![
            batch::RowRequest::Ready(row),
            batch::RowRequest::Ready(stale),
        ],
    };
    let mut session = batch::Session::begin(w.clone(), r).unwrap();
    let context = session.summary().context_revision;
    let good = session.advance(&context, 0).unwrap().row.unwrap();
    assert_eq!(good.status, Some(Status::ConsistentUnderAssumptions));
    assert_eq!(
        session.advance(&context, 1).unwrap().row.unwrap().state,
        batch::RowState::Error
    );
    let captured = session.into_capture(&context).unwrap();
    assert_eq!(
        serde_json::to_vec(&captured.workspace.documents).unwrap(),
        heads
    );
    assert_eq!(
        serde_json::to_vec(&captured.workspace.history).unwrap(),
        history
    );
    let restarted =
        workspace::load_workspace(&workspace::save_workspace(&captured.workspace).unwrap())
            .unwrap();
    let replay = batch::replay(&restarted, &captured.attachment_hash, 0).unwrap();
    assert_eq!(replay.current.status, good.status);
    assert_eq!(
        replay.current.declaration_revision,
        good.declaration_revision
    );
    let request = release::Request {
        batch_hash: captured.attachment_hash,
        rows: vec![release::Selection {
            construct_id: "edited-row".into(),
            options: cloning_outputs::Options {
                name: "Released literal edit".into(),
                show_prepared: false,
                save_primers: true,
                allow_incomplete: false,
            },
        }],
    };
    let mut release = release::ReleaseSession::begin(restarted, request).unwrap();
    let context = release.summary().context_revision;
    release.advance(&context, 0).unwrap();
    let inspected = release.inspect(0).unwrap();
    assert_eq!(inspected.product.sequence.as_str(), EDITED);
    let committed = release.into_commit(&context).unwrap();
    let manifest = release::load(&committed.workspace, &committed.manifest_hash).unwrap();
    assert!(manifest.rows[0].output.primer_attachment.is_some());
    let replay =
        run_lineage::replay_run(&committed.workspace, &manifest.rows[0].output.run.run_id).unwrap();
    assert!(replay.products_match);
    assert_eq!(
        replay.current_report.status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(
        replay
            .current_report
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        EDITED
    );
    for d in &w.documents {
        assert_eq!(
            committed
                .workspace
                .documents
                .iter()
                .find(|n| n.id == d.id)
                .unwrap()
                .head,
            d.head
        );
    }
}
#[test]
fn saved_proposal_requires_actual_workspace_and_full_source_and_refuses_documentary_tampering() {
    let (w, d, row) = fixture();
    let resolved = library::resolve(&d, &row, &w).unwrap();
    let library::IndependentIntent::Edits(intent) = resolved.intent.unwrap() else {
        unreachable!()
    };
    let request = primers::Request {
        intent: intent.request,
        overlap: 8,
        options: options(),
        profile_options: ProfileOptions::default(),
        finishing: collagene_core::assembly_recipe::OverlapFinishing::ReferenceOnly,
    };
    assert!(
        primers::source_preflight(
            &resolved
                .inputs
                .iter()
                .map(|i| (i.id.as_str(), &i.molecule))
                .collect::<Vec<_>>(),
            &request,
            None
        )
        .is_err()
    );
    let foreign = Workspace::from_documents(
        "unrelated",
        vec![Molecule::dna("same letters", INSERT, Topology::Linear).unwrap()],
        0,
    )
    .unwrap();
    assert!(primers::preflight(&resolved.inputs, &request, Some(&foreign)).is_err());
    let mut changed = resolved.inputs.clone();
    changed[1].molecule.name.push('!');
    changed[1].revision = changed[1].molecule.revision();
    assert!(primers::preflight(&changed, &request, Some(&w)).is_err());
    let p = primers::design(&resolved.inputs, &request, Some(&w)).unwrap();
    let mut tampered = p.clone();
    tampered.recipe.fragments[1].reverse = true;
    assert!(
        tampered
            .validate_sources(&resolved.inputs, &p.context_revision, Some(&w))
            .is_err()
    );
}
#[test]
fn substitutions_and_indels_have_independent_oriented_offsets_and_conflicts_are_not_reinterpreted()
{
    let (w, d, mut row) = fixture();
    let library::TargetDeclaration::OriginalRegions { substitutions, .. } = &mut row.target else {
        unreachable!()
    };
    substitutions.push(Substitution {
        role: "fragment-insert".into(),
        offset: 23,
        expected: "A".into(),
        replacement: "G".into(),
    });
    let r = library::resolve(&d, &row, &w).unwrap();
    assert_eq!(
        r.target.sequence.as_str(),
        "CACGTTCGTAGGACCAAGTCCGATGCTACGTATTTTCGGACCTGTAGATCGTCA"
    );
    let library::IndependentIntent::Edits(report) = r.intent.unwrap() else {
        unreachable!()
    };
    report
        .validate_sources(&r.inputs, &report.context_revision, Some(&w))
        .unwrap();
    let library::TargetDeclaration::OriginalRegions { substitutions, .. } = &mut row.target else {
        unreachable!()
    };
    substitutions[0].offset = 18;
    substitutions[0].expected = "A".into();
    assert!(
        library::resolve(&d, &row, &w).is_err(),
        "overlapping original coordinates must not become output offsets"
    );
}
#[test]
fn source_event_alias_and_declared_literal_limits_refuse_before_preparation() {
    let (w, mut d, row) = fixture();
    // The original stock is equal DNA but is not the native Derived product event.
    d.bins[1].choices[0].source = event(&w, 1);
    assert!(library::resolve(&d, &row, &w).is_err());
    let (w, d, mut row) = fixture();
    let library::TargetDeclaration::OriginalRegions { edits, .. } = &mut row.target else {
        unreachable!()
    };
    edits[0].replacement = "A".repeat(100_001);
    assert!(library::validate_construct(&d, &row).is_err());
    assert!(
        batch::Session::begin(
            w,
            batch::Request {
                design: Some(d),
                rows: vec![batch::RowRequest::Ready(row)]
            }
        )
        .is_err()
    );
}
