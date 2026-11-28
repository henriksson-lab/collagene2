use super::*;
use crate::assembly::draft::assembly_role_ids::Role;
use crate::cloning::preparation::cloning_methods_model::Treatment;
use collagene_core::{
    cloning::ligation::tailing::TailingBase, documents::workspace::Workspace,
    planning::plan::PlanInput,
};

fn captured(name: &str, dna: &str) -> Arc<CapturedSource> {
    let molecule = Molecule::dna(name, dna, Topology::Linear).unwrap();
    let workspace = Workspace::from_documents(name, vec![molecule.clone()], 0).unwrap();
    let document = &workspace.documents[0];
    Arc::new(CapturedSource {
        source: SourceInfo {
            binding: collagene_core::material::lineage::EventBinding {
                document_id: document.id.clone(),
                event_id: document.head_event.clone().unwrap(),
                revision: document.head.clone(),
            },
            name: molecule.name.clone(),
            length: molecule.sequence.len(),
            alphabet: molecule.sequence.alphabet(),
            topology: molecule.topology,
        },
        molecule,
    })
}
fn draft() -> Draft {
    let mut draft = Draft::default();
    draft.methods.mode = Mode::Tailing;
    draft.methods.source_roles = Some(vec![
        Role {
            input_id: "pinned-vector".into(),
            role: "vector-core".into(),
        },
        Role {
            input_id: "pinned-insert".into(),
            role: "insert-core".into(),
        },
    ]);
    for (slot, source) in draft.slots.iter_mut().zip([
        captured("Original vector", "AACC"),
        captured("Original insert", "GGTT"),
    ]) {
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                end: value,
                ..
            } = &mut slot.selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = source.molecule.sequence.len();
        }
        slot.captured = Some(source);
    }
    draft.sync_rows();
    draft
}

#[test]
fn two_source_tailing_roundtrip_preserves_pins_roles_preparations_and_literal_intent() {
    let mut draft = draft();
    {
        let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
            reverse: value,
            ..
        } = &mut draft.slots[1].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = true;
    }
    draft.methods.rows[0].preparation = OverlapPreparation::Digest;
    draft.methods.rows[0].enzymes = "EcoRV".into();
    draft.methods.rows[0].fragment = "2".into();
    draft.methods.rows[1].forward = "GG".into();
    draft.methods.rows[1].reverse_primer = "AA".into();
    draft.methods.tailing_profile = collagene_core::cloning::ligation::tailing_recipe::Profile::Gc;
    draft.methods.vector.action =
        collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
            base: TailingBase::C,
        };
    draft.methods.insert.action =
        collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
            base: TailingBase::G,
        };
    draft.methods.insert.treatment = Treatment::Phosphorylate;
    draft.junction_edits.push(assembly_junction_edit::Edit {
        junction: 0,
        remove_left: 0,
        remove_right: 0,
        inserted: "T".into(),
    });
    draft.sync_rows();
    let bytes = encode_raw(&draft).unwrap();
    let restored = decode(&bytes).unwrap();
    assert!(restored == draft);
    assert_eq!(restored.methods.rows[0].document, 0);
    assert_eq!(restored.methods.rows[1].document, 1);
    assert!(restored.methods.rows[1].reverse);
    assert_eq!(restored.junction_edits[0].inserted, "T");
    assert_eq!(restored.methods.rows[0].enzymes, "EcoRV");
}

#[test]
fn linear_terminal_intent_is_tailing_only_and_has_no_downstream_trim() {
    let mut draft = draft();
    draft.circular = false;
    draft.junction_edits = vec![assembly_junction_edit::Edit {
        junction: 1,
        remove_left: 1,
        remove_right: 0,
        inserted: "A".into(),
    }];
    draft.sync_rows();
    draft.validate().unwrap();
    assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
    draft.junction_edits[0].remove_right = 1;
    assert!(draft.validate().is_err());
    draft.junction_edits[0].remove_right = 0;
    for mode in [Mode::Overlap, Mode::GoldenGate, Mode::MultiSiteReference] {
        draft.methods.mode = mode;
        assert!(draft.validate().is_err());
    }
}

#[test]
fn moving_tailing_sources_keeps_each_sources_tail_chemistry_and_preparation() {
    let mut draft = draft();
    {
        let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
            reverse: value,
            ..
        } = &mut draft.slots[0].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = true;
    }
    draft.methods.rows[0].preparation = OverlapPreparation::Digest;
    draft.methods.rows[0].enzymes = "EcoRV".into();
    draft.methods.rows[1].forward = "GG".into();
    draft.methods.tailing_profile = collagene_core::cloning::ligation::tailing_recipe::Profile::Gc;
    draft.methods.vector.action =
        collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
            base: TailingBase::C,
        };
    draft.methods.vector.treatment = Treatment::Dephosphorylate;
    draft.methods.insert.action =
        collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
            base: TailingBase::G,
        };
    draft.methods.insert.treatment = Treatment::Phosphorylate;
    draft.sync_rows();
    let original = draft.clone();
    draft.slots.swap(0, 1);
    draft.methods.swap_rows(0, 1);
    draft.sync_rows();
    draft.validate().unwrap();
    assert!(draft.slots[0].captured == original.slots[1].captured);
    assert!(draft.slots[1].captured == original.slots[0].captured);
    assert!(
        draft.methods.vector.action
            == (collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                base: TailingBase::G
            })
    );
    assert!(draft.methods.vector.treatment == Treatment::Phosphorylate);
    assert!(
        draft.methods.insert.action
            == (collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                base: TailingBase::C
            })
    );
    assert!(draft.methods.insert.treatment == Treatment::Dephosphorylate);
    assert_eq!(draft.methods.rows[0].forward, "GG");
    assert_eq!(draft.methods.rows[1].enzymes, "EcoRV");
    assert!(!draft.methods.rows[0].reverse && draft.methods.rows[1].reverse);
    assert!(
        draft.methods.source_roles.as_ref().unwrap()[0]
            == original.methods.source_roles.as_ref().unwrap()[1]
    );
    draft.methods.swap_rows(0, 0);
    assert!(
        draft.methods.vector.action
            == (collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                base: TailingBase::G
            })
    );
    draft.slots.swap(1, 0);
    draft.methods.swap_rows(1, 0);
    draft.sync_rows();
    assert!(draft == original);
}

#[test]
fn kept_ends_and_added_tail_choices_follow_exact_sources_on_reorder() {
    use collagene_core::cloning::ligation::tailing_recipe::{Action, Profile};
    let mut draft = draft();
    draft.methods.tailing_profile = Profile::Ta;
    draft.methods.vector.action = Action::KeepExisting;
    draft.methods.vector.treatment = Treatment::Dephosphorylate;
    draft.methods.insert.action = Action::AddTail {
        base: TailingBase::A,
    };
    draft.methods.insert.treatment = Treatment::Phosphorylate;
    let before = draft.clone();
    draft.slots.swap(0, 1);
    draft.methods.swap_rows(0, 1);
    draft.sync_rows();
    assert!(draft.methods.vector.action == before.methods.insert.action);
    assert!(draft.methods.insert.action == Action::KeepExisting);
    assert!(draft.methods.vector.treatment == Treatment::Phosphorylate);
    assert!(draft.methods.insert.treatment == Treatment::Dephosphorylate);
    assert!(draft.methods.tailing_profile == Profile::Ta);
    assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
    draft.slots.swap(0, 1);
    draft.methods.swap_rows(0, 1);
    draft.sync_rows();
    assert!(draft == before);
}

#[test]
fn tailing_rejects_wrong_count_source_drift_and_unsupported_schema_but_retains_library_settings() {
    let mut invalid = draft();
    invalid.methods.rows[1].document = 0;
    assert!(invalid.validate().is_err());
    invalid.sync_rows();
    invalid.methods.rows[1].reverse = true;
    assert!(invalid.validate().is_err());
    invalid.sync_rows();
    invalid.version = CURRENT_VERSION - 1;
    assert!(invalid.validate().is_err());
    invalid.version = CURRENT_VERSION;
    invalid.slots.push(Slot::empty("Third source"));
    invalid.sync_rows();
    assert!(invalid.validate().is_err());

    let mut library =
        crate::assembly::draft::assembly_alternatives::enable(&Draft::default(), "ta-exclusion")
            .unwrap();
    let saved = &mut library.alternatives.as_mut().unwrap().rows[0].edits.methods;
    saved.mode = Mode::Tailing;
    library.methods.mode = Mode::Tailing;
    assert!(library.validate().is_ok());
    let stored = encode(&library).unwrap();
    let restored = decode(&stored).unwrap();
    assert!(
        restored.alternatives.as_ref().unwrap().rows[0]
            .edits
            .methods
            .mode
            == Mode::Tailing
    );
}

#[test]
fn circular_target_uses_only_operator_declared_junction_bases() {
    let mut draft = draft();
    draft.junction_edits = vec![
        assembly_junction_edit::Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "T".into(),
        },
        assembly_junction_edit::Edit {
            junction: 1,
            remove_left: 0,
            remove_right: 0,
            inserted: "A".into(),
        },
    ];
    draft.validate().unwrap();
    let roles = draft.methods.role_ids(2).unwrap();
    let inputs = draft
        .slots
        .iter()
        .zip(&roles)
        .map(|(slot, role)| {
            let captured = slot.captured.as_ref().unwrap();
            PlanInput {
                id: role.input_id.clone(),
                revision: captured.source.binding.revision.clone(),
                molecule: captured.molecule.clone(),
            }
        })
        .collect::<Vec<_>>();
    let base = collagene_core::cloning::assembly::intent::Request {
        material_bindings: Vec::new(),
        name: draft.name.clone(),
        circular: true,
        fragments: draft
            .slots
            .iter()
            .zip(&roles)
            .map(
                |(slot, role)| collagene_core::cloning::assembly::intent::Contribution {
                    role: role.role.clone(),
                    input_id: role.input_id.clone(),
                    selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                        start: slot.selector.reference_interval().unwrap().0,
                        end: slot.selector.reference_interval().unwrap().1,
                        reverse: slot.selector.reference_interval().unwrap().2,
                    },
                },
            )
            .collect(),
    };
    let target =
        assembly_junction_edit::build_intent(&inputs, &base, &draft.junction_edits, None).unwrap();
    assert_eq!(target.target.sequence.as_str(), "AACCTGGTTA");
    draft.methods.tailing_profile = collagene_core::cloning::ligation::tailing_recipe::Profile::Gc;
    draft.methods.vector.action =
        collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
            base: TailingBase::C,
        };
    draft.methods.insert.action =
        collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
            base: TailingBase::G,
        };
    draft.sync_rows();
    let unchanged =
        assembly_junction_edit::build_intent(&inputs, &base, &draft.junction_edits, None).unwrap();
    assert_eq!(target.context_revision, unchanged.context_revision);
    assert_eq!(target.target, unchanged.target);
}
