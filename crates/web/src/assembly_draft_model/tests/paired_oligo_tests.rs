use super::*;
use crate::assembly_role_ids::Role;
use collagene_core::workspace::Workspace;

fn captured(name: &str, dna: &str) -> Arc<CapturedSource> {
    captured_topology(name, dna, Topology::Linear)
}
fn captured_topology(name: &str, dna: &str, topology: Topology) -> Arc<CapturedSource> {
    let molecule = Molecule::dna(name, dna, topology).unwrap();
    let workspace = Workspace::from_documents(name, vec![molecule.clone()], 0).unwrap();
    let document = &workspace.documents[0];
    Arc::new(CapturedSource {
        source: SourceInfo {
            binding: collagene_core::run_lineage::EventBinding {
                document_id: document.id.clone(),
                event_id: document.head_event.clone().unwrap(),
                revision: document.head.clone(),
            },
            name: molecule.name.clone(),
            length: molecule.sequence.len(),
            alphabet: Alphabet::Dna,
            topology,
        },
        molecule,
    })
}
fn paired() -> Draft {
    let mut draft = Draft::default();
    draft.methods.mode = Mode::Ligation;
    draft.methods.source_roles = Some(vec![
        Role {
            input_id: "backbone-parent".into(),
            role: "backbone".into(),
        },
        Role {
            input_id: "adaptor-top".into(),
            role: "adaptor".into(),
        },
    ]);
    for (slot, source) in draft.slots.iter_mut().zip([
        captured("Backbone", "AAACCC"),
        captured("Top oligo", "AAGGCC"),
    ]) {
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut slot.selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = source.molecule.sequence.len();
        }
        slot.captured = Some(source);
    }
    draft.methods.rows[1].preparation = OverlapPreparation::AnnealedOligos;
    draft.methods.rows[1].anneal_options.min_paired_bases = 4;
    draft.slots[1].anneal_bottom = Some(captured("Bottom oligo", "GGCCTT"));
    draft.sync_rows();
    draft.methods.rows[1].anneal_candidate = Some(0);
    draft.validate().unwrap();
    draft
}

#[test]
fn ligation_intent_literals_and_linear_suffixes_are_independent_of_preparation() {
    let mut draft = paired();
    let preparation = draft.methods.clone();
    draft.junction_edits = vec![
        assembly_junction_edit::Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "G".into(),
        },
        assembly_junction_edit::Edit {
            junction: 1,
            remove_left: 0,
            remove_right: 0,
            inserted: "C".into(),
        },
    ];
    draft.validate().unwrap();
    assert!(draft.methods == preparation);
    assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
    draft.circular = false;
    draft.sync_rows();
    draft.validate().unwrap();
    draft.junction_edits[1].remove_right = 1;
    assert!(draft.validate().is_err());
    draft.junction_edits[1].remove_right = 0;
    draft.junction_edits[1].inserted = "A".repeat(201);
    assert!(draft.validate().is_err());
}

#[test]
fn swapped_ordered_roles_reset_intent_to_the_literal_old_bottom_and_clear_stale_pins() {
    use collagene_core::{
        chemistry::{FivePrime, Terminus, ThreePrime},
        verification::PlanInput,
    };
    let mut draft = paired();
    draft.slots[0].captured = Some(captured_topology(
        "Backbone",
        "GGAATTCTTACG",
        Topology::Circular,
    ));
    {
        let collagene_core::assembly_intent::Selector::ReferenceRange { start: value, .. } =
            &mut draft.slots[0].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = 2;
    }
    {
        let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
            &mut draft.slots[0].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = 14;
    }
    draft.slots[1].captured = Some(captured("Top", "AATTGCTACG"));
    draft.slots[1].anneal_bottom = Some(captured("Bottom", "AATTCGTAGC"));
    {
        let collagene_core::assembly_intent::Selector::ReferenceRange { start: value, .. } =
            &mut draft.slots[1].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = 1;
    }
    {
        let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
            &mut draft.slots[1].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = 9;
    }
    {
        let collagene_core::assembly_intent::Selector::ReferenceRange { reverse: value, .. } =
            &mut draft.slots[1].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = true;
    }
    draft.slots[1].digest_geometry = Some(DigestGeometry {
        material: None,
        top_start: 1,
        top_length: 8,
        bottom_start: 1,
        bottom_length: 8,
    });
    draft.methods.rows[1].anneal_options.min_paired_bases = 6;
    draft.methods.rows[1].anneal_options.top_chemistry = Terminus {
        five_prime: FivePrime::Hydroxyl,
        three_prime: ThreePrime::Blocked,
    };
    draft.methods.rows[1].anneal_options.bottom_chemistry = Terminus::CLEAVED;
    draft.methods.rows[1].treatment = crate::cloning_methods_model::Treatment::Phosphorylate;
    draft.junction_edits = vec![assembly_junction_edit::Edit {
        junction: 0,
        remove_left: 1,
        remove_right: 0,
        inserted: "T".into(),
    }];
    draft.sync_rows();
    draft.methods.rows[1].anneal_candidate = Some(0);
    draft.validate().unwrap();
    let before = draft.clone();
    let swapped = draft.swap_oligo_roles(1).unwrap();
    assert!(draft == before);
    assert!(swapped.slots[1].captured == before.slots[1].anneal_bottom);
    assert!(swapped.slots[1].anneal_bottom == before.slots[1].captured);
    assert_eq!(
        (
            swapped.slots[1].selector.reference_interval().unwrap().0,
            swapped.slots[1].selector.reference_interval().unwrap().1,
            swapped.slots[1].selector.reference_interval().unwrap().2
        ),
        (0, 10, false)
    );
    assert_eq!(
        swapped.methods.rows[1].anneal_options.top_chemistry,
        Terminus::CLEAVED
    );
    assert_eq!(
        swapped.methods.rows[1].anneal_options.bottom_chemistry,
        before.methods.rows[1].anneal_options.top_chemistry
    );
    assert!(swapped.methods.rows[1].treatment == before.methods.rows[1].treatment);
    assert!(swapped.methods.rows[0] == before.methods.rows[0]);
    assert_eq!(swapped.methods.source_roles, before.methods.source_roles);
    assert_eq!(swapped.methods.rows[1].anneal_candidate, None);
    assert_ne!(
        swapped.methods.rows[1].anneal_source_context,
        before.methods.rows[1].anneal_source_context
    );
    assert!(swapped.slots[1].digest_geometry.is_none() && swapped.junction_edits.is_empty());
    assert!(decode(&encode_raw(&swapped).unwrap()).unwrap() == swapped);
    let ids = swapped.methods.all_input_ids().unwrap();
    let sources = [
        swapped.slots[0].captured.as_ref().unwrap(),
        swapped.slots[1].captured.as_ref().unwrap(),
        swapped.slots[1].anneal_bottom.as_ref().unwrap(),
    ];
    let inputs = sources
        .into_iter()
        .zip(&ids)
        .map(|(capture, id)| PlanInput {
            id: id.clone(),
            revision: capture.source.binding.revision.clone(),
            molecule: capture.molecule.clone(),
        })
        .collect::<Vec<_>>();
    let roles = swapped.methods.role_ids(2).unwrap();
    let request = collagene_core::assembly_intent::Request {
        material_bindings: Vec::new(),
        name: "Independently reversed adaptor".into(),
        circular: true,
        fragments: swapped
            .slots
            .iter()
            .zip(&roles)
            .map(
                |(slot, role)| collagene_core::assembly_intent::Contribution {
                    role: role.role.clone(),
                    input_id: role.input_id.clone(),
                    selector: collagene_core::assembly_intent::Selector::ReferenceRange {
                        start: slot.selector.reference_interval().unwrap().0,
                        end: slot.selector.reference_interval().unwrap().1,
                        reverse: slot.selector.reference_interval().unwrap().2,
                    },
                },
            )
            .collect(),
    };
    let intent = collagene_core::assembly_intent::build(&inputs, &request, None).unwrap();
    assert_eq!(intent.target.sequence.as_str(), "AATTCTTACGGGAATTCGTAGC");
}

#[test]
fn swapping_one_pair_preserves_unrelated_four_slot_junction_intent() {
    let mut draft = paired();
    for (label, dna, id) in [
        ("Third", "GTAC", "third-parent"),
        ("Fourth", "CAAT", "fourth-parent"),
    ] {
        let source = captured(label, dna);
        let mut slot = Slot::empty(label);
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut slot.selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = source.molecule.sequence.len();
        }
        slot.captured = Some(source);
        draft.slots.push(slot);
        draft.methods.rows.push(OverlapInput::default());
        draft.methods.source_roles.as_mut().unwrap().push(Role {
            input_id: id.into(),
            role: label.to_lowercase(),
        });
    }
    draft.methods.rows[2].preparation = OverlapPreparation::ExactPcr;
    draft.methods.rows[2].forward = "GT".into();
    draft.methods.rows[2].reverse_primer = "GT".into();
    draft.junction_edits = vec![
        assembly_junction_edit::Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 0,
            inserted: "A".into(),
        },
        assembly_junction_edit::Edit {
            junction: 1,
            remove_left: 0,
            remove_right: 0,
            inserted: "C".into(),
        },
        assembly_junction_edit::Edit {
            junction: 2,
            remove_left: 1,
            remove_right: 1,
            inserted: "GG".into(),
        },
        assembly_junction_edit::Edit {
            junction: 3,
            remove_left: 1,
            remove_right: 0,
            inserted: "T".into(),
        },
    ];
    draft.sync_rows();
    draft.methods.rows[1].anneal_candidate = Some(0);
    let before = draft.clone();
    let swapped = draft.swap_oligo_roles(1).unwrap();
    assert_eq!(swapped.junction_edits, before.junction_edits[2..]);
    assert!(swapped.slots[0] == before.slots[0] && swapped.slots[2..] == before.slots[2..]);
    assert!(
        swapped.methods.rows[0] == before.methods.rows[0]
            && swapped.methods.rows[2..] == before.methods.rows[2..]
    );
    assert_eq!(swapped.methods.source_roles, before.methods.source_roles);
    assert_eq!(swapped.methods.rows[1].anneal_candidate, None);
    assert!(swapped.slots[1].captured == before.slots[1].anneal_bottom);
    assert!(decode(&encode_raw(&swapped).unwrap()).unwrap() == swapped);
    assert!(draft == before);
    // The first circular contribution has the closing boundary as its
    // incoming neighbor; it must not erase the middle two declarations.
    draft.slots.swap(0, 1);
    draft.methods.swap_rows(0, 1);
    draft.sync_rows();
    let first = draft.swap_oligo_roles(0).unwrap();
    assert_eq!(first.junction_edits, draft.junction_edits[1..3]);
}

#[test]
fn swapping_originals_and_full_termini_matches_asymmetric_physical_reversal() {
    use collagene_core::chemistry::{self, FivePrime, Terminus, ThreePrime};
    let mut draft = paired();
    draft.slots[1].captured = Some(captured("Top", "AAGCTACG"));
    {
        let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
            &mut draft.slots[1].selector
        else {
            panic!("fixture requires reference selector")
        };
        *value = 8;
    }
    draft.slots[1].anneal_bottom = Some(captured("Bottom", "TTTCGTAGC"));
    draft.methods.rows[1].anneal_options.min_paired_bases = 6;
    draft.methods.rows[1].anneal_options.top_chemistry = Terminus {
        five_prime: FivePrime::Unknown,
        three_prime: ThreePrime::Blocked,
    };
    draft.methods.rows[1].anneal_options.bottom_chemistry = Terminus::CLEAVED;
    draft.sync_rows();
    let original = chemistry::anneal_oligos(
        &draft.slots[1].captured.as_ref().unwrap().molecule.sequence,
        &draft.slots[1]
            .anneal_bottom
            .as_ref()
            .unwrap()
            .molecule
            .sequence,
        &draft.methods.rows[1].anneal_options,
    )
    .unwrap();
    let swapped = draft.swap_oligo_roles(1).unwrap();
    let fresh = chemistry::anneal_oligos(
        &swapped.slots[1]
            .captured
            .as_ref()
            .unwrap()
            .molecule
            .sequence,
        &swapped.slots[1]
            .anneal_bottom
            .as_ref()
            .unwrap()
            .molecule
            .sequence,
        &swapped.methods.rows[1].anneal_options,
    )
    .unwrap();
    assert!(original.complete && fresh.complete && !original.candidates.is_empty());
    assert_eq!(original.candidates.len(), fresh.candidates.len());
    for candidate in &original.candidates {
        let reversed = candidate.duplex.reverse().unwrap();
        assert!(
            fresh
                .candidates
                .iter()
                .any(|candidate| candidate.duplex == reversed)
        );
    }
    let asymmetric = original
        .candidates
        .iter()
        .find(|candidate| candidate.duplex.bottom_start == 2)
        .unwrap();
    let (left, right) = asymmetric.duplex.ends().unwrap();
    assert_eq!(left.overhang.len(), 2);
    assert_eq!(right.overhang.len(), 3);
    assert!(asymmetric.duplex.nicks.is_empty());
    let before = serde_json::to_value(&draft).unwrap();
    assert!(draft.swap_oligo_roles(0).is_err());
    assert!(draft.swap_oligo_roles(99).is_err());
    assert_eq!(serde_json::to_value(&draft).unwrap(), before);
}

#[test]
fn paired_draft_flattens_three_exact_originals_and_roundtrips_current_schema() {
    let draft = paired();
    assert_eq!(draft.version, CURRENT_VERSION);
    let expected = format!("anneal-bottom-{}", blake3::hash(b"adaptor-top").to_hex());
    assert_eq!(draft.methods.bottom_input_id(1).unwrap(), expected);
    assert_eq!(
        draft.methods.all_input_ids().unwrap(),
        vec![
            "backbone-parent".to_string(),
            "adaptor-top".into(),
            expected
        ]
    );
    assert_eq!(draft.methods.source_indices().unwrap(), vec![0, 1, 2]);
    assert_eq!(
        crate::assembly_alternatives::all_captures(&draft).count(),
        3,
        "Workspace restore must validate the bottom event as well as primary contributions"
    );
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    assert!(restored == draft);
    assert_eq!(restored.methods.rows[1].anneal_candidate, Some(0));
    let mut legacy = draft.clone();
    legacy.version = CURRENT_VERSION - 1;
    assert!(legacy.validate().is_err());
    assert!(
        crate::assembly_alternatives::for_storage(&draft)
            .unwrap()
            .version
            == CURRENT_VERSION
    );
}

#[test]
fn moving_source_slots_keeps_pair_events_options_candidate_and_stable_identity() {
    let mut draft = paired();
    let before = draft.clone();
    let bottom_id = draft.methods.bottom_input_id(1).unwrap();
    draft.slots.swap(0, 1);
    draft.methods.swap_rows(0, 1);
    draft.sync_rows();
    draft.validate().unwrap();
    assert!(draft.slots[0].anneal_bottom == before.slots[1].anneal_bottom);
    assert_eq!(draft.methods.bottom_input_id(0).unwrap(), bottom_id);
    assert_eq!(draft.methods.rows[0].anneal_candidate, Some(0));
    assert_eq!(
        draft.methods.rows[0].anneal_source_context,
        before.methods.rows[1].anneal_source_context
    );
    assert_eq!(draft.methods.source_indices().unwrap(), vec![0, 1, 2]);
    draft.slots.swap(1, 0);
    draft.methods.swap_rows(1, 0);
    draft.sync_rows();
    assert!(draft == before);
}

#[test]
fn equal_dna_different_bottom_event_and_changed_options_clear_candidate() {
    let mut draft = paired();
    let prior = draft.methods.rows[1].anneal_source_context.clone();
    let bottom = Arc::make_mut(draft.slots[1].anneal_bottom.as_mut().unwrap());
    bottom.source.binding.event_id = "another-original-bottom-event".into();
    assert!(draft.validate().is_err());
    draft.sync_rows();
    assert_eq!(draft.methods.rows[1].anneal_candidate, None);
    assert_ne!(draft.methods.rows[1].anneal_source_context, prior);
    draft.methods.rows[1].anneal_candidate = Some(0);
    draft.methods.rows[1].anneal_options.min_paired_bases = 5;
    draft.sync_rows();
    assert_eq!(draft.methods.rows[1].anneal_candidate, None);
    draft.slots[1].anneal_bottom = None;
    draft.sync_rows();
    assert!(draft.methods.source_indices().is_err());
    assert_eq!(draft.methods.rows[1].anneal_bottom_document, None);
}

#[test]
fn collisions_non_ligation_libraries_and_unbounded_oligos_are_rejected() {
    let mut draft = paired();
    let secondary = draft.methods.bottom_input_id(1).unwrap();
    draft.methods.source_roles.as_mut().unwrap()[0].input_id = secondary;
    assert!(draft.methods.all_input_ids().is_err());
    let mut wrong = paired();
    wrong.methods.mode = Mode::Overlap;
    assert!(wrong.validate().is_err());
    let mut library = crate::assembly_alternatives::enable(&Draft::default(), "pair-ban").unwrap();
    library.slots[0].anneal_bottom = paired().slots[1].anneal_bottom.clone();
    library.sync_rows();
    assert!(library.validate().unwrap_err().contains("libraries"));
    let mut too_long = paired();
    too_long.slots[1].anneal_bottom = Some(captured(
        "Oversized oligo",
        &"A".repeat(collagene_core::chemistry::MAX_ANNEAL_OLIGO_BASES + 1),
    ));
    too_long.sync_rows();
    assert!(too_long.validate().is_err());
}
