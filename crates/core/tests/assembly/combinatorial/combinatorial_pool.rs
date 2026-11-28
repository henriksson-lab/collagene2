//! Literal whole-reaction substrates over independent source/target declarations.
use collagene_core::{
    cloning::{
        assembly::recipe::Preparation,
        golden_gate::{GoldenGateOptions, prepared::Selection, reaction as golden_gate_reaction},
    },
    documents::workspace::Workspace,
    material::{lineage::EventBinding, outputs as cloning_outputs},
    planning::{
        combinatorial::{
            self as design, Bin, Choice, ChoiceSelection, Design, Method, PoolSelection,
            PoolSubstrate, PoolSubstrateKind, SelectedConstruct, TargetDeclaration, batch, release,
        },
        plan::Operation,
    },
    sequence::analysis::Primer,
    verification::evidence::Status,
};
use collagene_model::{Molecule, Topology};
const A: &str = "ATGACCCCGGGGTTTCCAGA";
const B: &str = "ATGAGTTAGCGACCTTGGAA";
const AB: &str = "ATGACCCCGGGGTTTCCAGAATGAGTTAGCGACCTTGGAA";
fn binding(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture() -> (Workspace, Design, SelectedConstruct) {
    let w = Workspace::from_documents(
        "literal reaction library",
        vec![
            Molecule::dna(
                "whole vector",
                &format!("ATGCGGTCTCA{A}ATGATGAGACCGCAT"),
                Topology::Circular,
            )
            .unwrap(),
            Molecule::dna("insert template", B, Topology::Linear).unwrap(),
            Molecule::dna("uncut carryover", "ACGTACGTACGT", Topology::Circular).unwrap(),
            Molecule::dna("independent intended AB", AB, Topology::Circular).unwrap(),
            Molecule::dna("independent wrong target", "AACCGGTT", Topology::Circular).unwrap(),
            Molecule::dna(
                "whole vector",
                &format!("ATGCGGTCTCA{A}ATGATGAGACCGCAT"),
                Topology::Circular,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let ids = ["backbone", "insert", "carry"];
    let bins = ids
        .iter()
        .enumerate()
        .map(|(i, id)| Bin {
            id: (*id).into(),
            label: (*id).into(),
            choices: vec![Choice {
                anneal_bottom: None,
                generation: Some(
                    collagene_core::chemistry::tracked_duplex::GenerationAssumptions {
                        dam: collagene_core::chemistry::methylation::MethylationState::Unmodified,
                        dcm: collagene_core::chemistry::methylation::MethylationState::Unmodified,
                        cpg: collagene_core::chemistry::methylation::MethylationState::Unmodified,
                    },
                ),
                id: format!("{id}-choice"),
                label: (*id).into(),
                source: binding(&w, i),
                selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    start: if i == 0 { 11 } else { 0 },
                    end: if i == 0 {
                        11 + A.len()
                    } else if i == 1 {
                        B.len()
                    } else {
                        12
                    },
                    reverse: false,
                },
                preparation_reverse: false,
                preparation: if i == 1 {
                    Some(Preparation::ExactPcr {
                        forward: Primer::new("manual forward", "ATGAGTTA", "ATGCGGTCTCA").unwrap(),
                        reverse: Primer::new("manual reverse", "TTCCAAGG", "ATGCGGTCTCATCAT")
                            .unwrap(),
                        template_strand: collagene_core::chemistry::Strand::Top,
                    })
                } else {
                    None
                },
                treatment: None,
            }],
        })
        .collect();
    let d = Design {
        material_bindings: vec![],
        id: "reaction-library".into(),
        name: "Whole-material library".into(),
        circular: true,
        bins,
        method: Method::GoldenGatePool {
            options: GoldenGateOptions {
                max_fragments: 3,
                ..Default::default()
            },
            substrates: ids
                .iter()
                .enumerate()
                .map(|(i, id)| PoolSubstrate {
                    bin_id: (*id).into(),
                    kind: if i == 1 {
                        PoolSubstrateKind::Prepared
                    } else {
                        PoolSubstrateKind::WholeSource
                    },
                })
                .collect(),
            selection: None,
        },
    };
    let row = SelectedConstruct {
        id: "chosen-ab".into(),
        name: "Chosen AB".into(),
        choices: ids
            .iter()
            .map(|id| ChoiceSelection {
                bin_id: (*id).into(),
                choice_id: format!("{id}-choice"),
            })
            .collect(),
        target: TargetDeclaration::Existing(binding(&w, 3)),
        overrides: None,
    };
    (w, d, row)
}
fn pin_literal_ab(w: &Workspace, d: &mut Design, row: &SelectedConstruct) {
    let report = design::inspect_pool(d, row, w).unwrap();
    let pool = report.pool.as_ref().unwrap();
    // Literal fixture DNA chooses the test's user action, never an engine target match.
    let (index, c) = pool
        .candidates
        .iter()
        .enumerate()
        .find(|(_, c)| c.restriction_site_free && c.duplex.top.as_str() == AB)
        .expect("literal AB candidate");
    let pin = Selection {
        context_revision: report.context_revision,
        candidate_index: index,
        directed_sequence_key: c.directed_sequence_key.clone(),
        physical_sequence_key: c.physical_sequence_key.clone(),
    };
    let resolved = design::resolve_pool(d, row, w).unwrap();
    let design::Recipe::GoldenGatePool {
        source_revision, ..
    } = resolved.recipe
    else {
        unreachable!()
    };
    let Method::GoldenGatePool { selection, .. } = &mut d.method else {
        unreachable!()
    };
    *selection = Some(PoolSelection {
        selection: pin,
        source_revision,
    });
}
#[test]
fn all_digest_material_and_carryover_remain_visible_without_implicit_selection() {
    let (w, d, row) = fixture();
    let snapshot = serde_json::to_vec(&w).unwrap();
    let resolved = design::resolve_pool(&d, &row, &w).unwrap();
    assert_eq!(resolved.target.sequence.as_str(), AB);
    assert_eq!(
        resolved.bindings[0].event.event_id,
        w.documents[0].head_event.clone().unwrap()
    );
    let design::Recipe::GoldenGatePool {
        request, selection, ..
    } = &resolved.recipe
    else {
        unreachable!()
    };
    assert!(selection.is_none());
    assert!(matches!(
        request.rows[0].substrate,
        golden_gate_reaction::Substrate::WholeSource { origin: None }
    ));
    let report = design::inspect_pool(&d, &row, &w).unwrap();
    let pool = report.pool.unwrap();
    assert_eq!(pool.fragments.len(), 6); // two vector, three complete amplicon, uncut circle
    assert_eq!(pool.carryover_sources, vec![2]);
    assert!(pool.candidates.iter().any(|c| !c.restriction_site_free));
    assert!(design::compile(&d, &row, &w).is_err());
    assert!(design::evaluate(&d, &row, &w).is_err());
    assert_eq!(serde_json::to_vec(&w).unwrap(), snapshot);
}
#[test]
fn explicit_native_pin_preserves_preparation_and_wrong_target_cannot_choose_a_product() {
    let (w, mut d, row) = fixture();
    pin_literal_ab(&w, &mut d, &row);
    let compiled = design::compile(&d, &row, &w).unwrap();
    assert_eq!(compiled.plan.steps.len(), 5);
    assert!(matches!(
        compiled.plan.steps[0],
        Operation::DeclareSourceModifications { .. }
    ));
    assert!(matches!(compiled.plan.steps[1], Operation::Pcr { .. }));
    assert!(matches!(
        compiled.plan.steps[4],
        Operation::SelectedGoldenGate { .. }
    ));
    let evaluated = design::evaluate(&d, &row, &w).unwrap();
    assert_eq!(
        evaluated.verification().status,
        Status::ConsistentUnderAssumptions
    );
    assert_eq!(
        evaluated
            .verification()
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        AB
    );
    assert!(evaluated.preparation_complete());
    let mut wrong = row.clone();
    wrong.target = TargetDeclaration::Existing(binding(&w, 4));
    assert!(
        design::compile(&d, &wrong, &w).is_err(),
        "target change invalidates old discovery context"
    );
    // Reconfirm the SAME literal candidate against a wrong independently declared target.
    let Method::GoldenGatePool { selection, .. } = &mut d.method else {
        unreachable!()
    };
    *selection = None;
    pin_literal_ab(&w, &mut d, &wrong);
    let failed = design::evaluate(&d, &wrong, &w).unwrap();
    assert_eq!(failed.verification().status, Status::Invalid);
    assert_eq!(
        failed
            .verification()
            .products
            .last()
            .unwrap()
            .sequence
            .as_str(),
        AB
    );
}
#[test]
fn changed_search_context_or_forged_candidate_keys_never_gain_authority() {
    let (w, mut d, row) = fixture();
    pin_literal_ab(&w, &mut d, &row);
    let original = d.clone();
    let Method::GoldenGatePool { options, .. } = &mut d.method else {
        unreachable!()
    };
    options.max_states += 1;
    assert!(design::compile(&d, &row, &w).is_err());
    d = original;
    let Method::GoldenGatePool { selection, .. } = &mut d.method else {
        unreachable!()
    };
    selection.as_mut().unwrap().selection.physical_sequence_key = "0".repeat(64);
    let invalid = design::evaluate(&d, &row, &w).unwrap();
    assert_eq!(invalid.verification().status, Status::Invalid);
    assert_eq!(
        invalid.verification().products.len(),
        4,
        "only three declared-generation prefixes and PCR exist; forged pool candidate was withheld"
    );
}
#[test]
fn bounded_shapes_and_whole_source_chemistry_are_explicit() {
    let (w, mut d, row) = fixture();
    let original = d.clone();
    let Method::GoldenGatePool { options, .. } = &mut d.method else {
        unreachable!()
    };
    options.circular = false;
    assert!(d.validate_construct(&row).is_err());
    d = original.clone();
    let Method::GoldenGatePool { substrates, .. } = &mut d.method else {
        unreachable!()
    };
    substrates[0].bin_id = "foreign-bin".into();
    assert!(d.validate_construct(&row).is_err());
    d = original;
    d.bins[0].choices[0].treatment = Some(collagene_core::chemistry::EndTreatment::Phosphorylate);
    assert!(design::resolve_pool(&d, &row, &w).is_err());
}
#[test]
fn equal_dna_event_rebinding_and_changed_original_geometry_require_new_consent() {
    let (w, mut d, row) = fixture();
    pin_literal_ab(&w, &mut d, &row);
    let original = d.clone();
    let original_revision = d.bins[0].choices[0].source.revision.clone();
    d.bins[0].choices[0].source = binding(&w, 5);
    assert_eq!(d.bins[0].choices[0].source.revision, original_revision);
    let resolved = design::resolve_pool(&d, &row, &w).unwrap();
    assert_eq!(resolved.bindings[0].event.document_id, w.documents[5].id);
    assert!(
        design::compile(&d, &row, &w).is_err(),
        "equal letters do not authorize a new source event"
    );
    d = original;
    let collagene_core::cloning::assembly::intent::Selector::ReferenceRange { start, .. } =
        &mut d.bins[0].choices[0].selector
    else {
        unreachable!()
    };
    *start += 1;
    assert!(
        design::compile(&d, &row, &w).is_err(),
        "reaction substrate is unchanged but original design intent changed"
    );
}
#[test]
fn selected_pool_crosses_historical_batch_into_fresh_atomic_release() {
    let (w, mut d, row) = fixture();
    pin_literal_ab(&w, &mut d, &row);
    let mut batch = batch::Session::begin(
        w,
        batch::Request {
            design: Some(d),
            rows: vec![batch::RowRequest::Ready(row.clone())],
        },
    )
    .unwrap();
    let context = batch.summary().context_revision;
    let progress = batch.advance(&context, 0).unwrap();
    assert_eq!(
        progress.row.unwrap().status,
        Some(Status::ConsistentUnderAssumptions)
    );
    let historical = batch.into_capture(&context).unwrap();
    let request = release::Request {
        batch_hash: historical.attachment_hash,
        rows: vec![release::Selection {
            construct_id: row.id,
            options: cloning_outputs::Options {
                name: "Released AB".into(),
                save_primers: true,
                ..Default::default()
            },
        }],
    };
    let mut session = release::ReleaseSession::begin(historical.workspace, request).unwrap();
    let context = session.summary().context_revision;
    session.advance(&context, 0).unwrap();
    let inspection = session.inspect(0).unwrap();
    assert_eq!(inspection.product.sequence.as_str(), AB);
    assert_eq!(inspection.product.name, "Released AB");
    assert!(inspection.row.output.primer_attachment.is_some());
    let commit = session.into_commit(&context).unwrap();
    let manifest = release::load(&commit.workspace, &commit.manifest_hash).unwrap();
    manifest
        .validate_current_declarations(&commit.workspace)
        .unwrap();
    assert_eq!(manifest.rows.len(), 1);
}
