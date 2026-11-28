//! Literal source/bond fixtures; no second recognition or digestion implementation.
use super::*;
use crate::chemistry::methylation::MethylationState;

fn declared() -> SourceInput {
    SourceInput::DeclaredOriginal {
        generation: GenerationAssumptions {
            host: None,
            dam: MethylationState::Unmodified,
            dcm: MethylationState::Unmodified,
            cpg: MethylationState::Unmodified,
        },
    }
}
fn request(source: &Molecule, enzyme: &str, start: usize, length: usize, reverse: bool) -> Request {
    Request {
        ivt: ivt::Request {
            source_revision: source.revision(),
            profile: ivt::profiles().unwrap().remove(0).reference,
            transcript: ivt::Transcript {
                start,
                length,
                coding_strand: if reverse {
                    ivt::CodingStrand::Reverse
                } else {
                    ivt::CodingStrand::Forward
                },
            },
            cap: ivt::CapChoice::CleanCapAgN7113,
            initiator: ivt::InitiatorRequirement::NoAdditionalRequirement,
            linearization: None,
            poly_a: None,
            provenance: Provenance {
                statement: "Explicit source-bound IVT interval, no physical preparation asserted"
                    .into(),
                reference: None,
            },
        },
        source: declared(),
        selection: catalog::select(&catalog::builtin_catalog().unwrap(), &[enzyme.into()]).unwrap(),
        conditions: None,
        policy: Policy {
            allowed_end_kinds: vec![ivt::EndKind::Blunt, ivt::EndKind::FivePrimeOverhang],
            require_known_conditions: false,
            require_resolved_ivt_checks: false,
            require_complete_origin: true,
        },
        limits: Limits {
            max_comparisons: 1_000_000,
            max_sites_per_enzyme: 64,
            max_report_bytes: MAX_REPORT_BYTES,
        },
    }
}
fn source(site: &str, circular: bool) -> Molecule {
    Molecule::dna(
        "Literal T7 source",
        &format!("TAATACGACTCACTATAAGTT{site}TT"),
        if circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )
    .unwrap()
}
fn reverse(source: &Molecule) -> Molecule {
    let mut reversed = source.clone();
    reversed.sequence = source.sequence.reverse_complement().unwrap();
    reversed
}
fn rotate(source: &Molecule, offset: usize) -> Molecule {
    let text = source.sequence.as_str();
    let mut rotated = source.clone();
    rotated.sequence =
        collagene_model::Sequence::dna(&format!("{}{}", &text[offset..], &text[..offset])).unwrap();
    rotated
}

#[test]
fn literal_blunt_and_staggered_bonds_map_to_forward_and_reverse_runoff_without_editing_source() {
    for (site, enzyme, length, coding, template, kind) in [
        ("CCCGGG", "smai", 7, 7, 7, ivt::EndKind::Blunt),
        ("GAATTC", "ecori", 9, 5, 9, ivt::EndKind::FivePrimeOverhang),
    ] {
        for backwards in [false, true] {
            let original = source(site, false);
            let molecule = if backwards {
                reverse(&original)
            } else {
                original
            };
            let before = molecule.clone();
            let start = if backwards {
                molecule.sequence.len() - 18
            } else {
                17
            };
            let q = request(&molecule, enzyme, start, length, backwards);
            let report = assess(&molecule, &q).unwrap();
            let candidate = &report.candidates[0];
            assert_eq!(
                candidate.state,
                CandidateState::EligibleUnderDeclarations,
                "{:?}",
                candidate.issues
            );
            let mapping = candidate.mapping.as_ref().unwrap();
            assert_eq!(
                (mapping.coding_cut_offset, mapping.template_cut_offset),
                (coding, template)
            );
            let ivt = candidate.ivt.as_ref().unwrap();
            assert_eq!(ivt.linearization.as_ref().unwrap().kind, kind);
            assert_eq!(
                ivt.request
                    .linearization
                    .as_ref()
                    .unwrap()
                    .template_cut_offset,
                length
            );
            assert!(candidate.retained_fragment.is_some());
            assert_eq!(
                candidate.fragments.len(),
                2,
                "linear single cut retains both fragment geometries"
            );
            assert_eq!(molecule, before);
            report.validate_for(&molecule, &q).unwrap();
        }
    }
}

#[test]
fn circular_origin_crossing_preserves_exact_stagger_and_template_fragment() {
    let original = source("GAATTC", true);
    // Original EcoRI bonds22/26 become27/2 after rotating29bases by24.
    // +1 moves17->22; the lifted coding/template offsets remain5/9.
    let molecule = rotate(&original, 24);
    let q = request(&molecule, "ecori", 22, 9, false);
    let report = assess(&molecule, &q).unwrap();
    let candidate = &report.candidates[0];
    assert_eq!(
        candidate.state,
        CandidateState::EligibleUnderDeclarations,
        "{:?}",
        candidate.issues
    );
    let mapping = candidate.mapping.as_ref().unwrap();
    assert_eq!((mapping.top_bond, mapping.bottom_bond), (27, 2));
    assert_eq!(
        (
            mapping.coding_cut_offset,
            mapping.template_cut_offset,
            mapping.signed_template_minus_coding
        ),
        (5, 9, 4)
    );
    assert_eq!(candidate.fragments.len(), 1);
    assert_eq!(candidate.retained_fragment, Some(0));
    assert_eq!(candidate.fragments[0].geometry.topology, Topology::Linear);
    assert_eq!(candidate.ivt.as_ref().unwrap().transcript_spans.len(), 2);
    let reverse_source = reverse(&molecule);
    let reverse_start = molecule.sequence.len() - 1 - 22;
    let reverse_q = request(&reverse_source, "ecori", reverse_start, 9, true);
    let reversed = assess(&reverse_source, &reverse_q).unwrap();
    assert_eq!(
        reversed.candidates[0].state,
        CandidateState::EligibleUnderDeclarations
    );
    assert_eq!(
        reversed.candidates[0]
            .mapping
            .as_ref()
            .unwrap()
            .template_cut_offset,
        9
    );
}

#[test]
fn a_second_turn_is_not_invented_to_match_a_short_runoff_declaration() {
    let molecule = Molecule::dna(
        "coordinate-only fixture",
        &"A".repeat(29),
        Topology::Circular,
    )
    .unwrap();
    let selected = catalog::select(&catalog::builtin_catalog().unwrap(), &["ecori".into()])
        .unwrap()
        .enzymes
        .remove(0);
    // Unit input is a literal pair of already-established bonds, not a claim
    // that this coordinate-only DNA contains EcoRI. Independent modulo lifting
    // would incorrectly turn template offset32 into3 and hide the crossed +1.
    let site = tracked::CheckedSite {
        enzyme: "EcoRI".into(),
        recognition_start: 4,
        reverse: false,
        top_cut: Some(5),
        bottom_cut: Some(9),
        outcome: Outcome::PermittedUnderProfile,
        evidence: vec![],
        reason: "coordinate-only arithmetic fixture".into(),
    };
    let transcript = ivt::Transcript {
        start: 6,
        length: 3,
        coding_strand: ivt::CodingStrand::Forward,
    };
    let projected = project(
        molecule.sequence.len(),
        molecule.topology,
        &transcript,
        &selected,
        &site,
        0,
    )
    .unwrap();
    assert_eq!(
        (projected.coding_cut_offset, projected.template_cut_offset),
        (28, 32)
    );
    assert_ne!(projected.template_cut_offset, 3);
    let mut malformed = site;
    malformed.top_cut = Some(usize::MAX);
    assert!(
        project(
            molecule.sequence.len(),
            molecule.topology,
            &transcript,
            &selected,
            &malformed,
            0,
        )
        .is_none(),
        "untrusted documentary coordinates are bounded before signed arithmetic"
    );
}

#[test]
fn promoter_disruption_wrong_runoff_multicuts_and_unknown_cleavage_remain_explicit() {
    let molecule = source("CCCGGG", true);
    // HinfI recognizes GACTC inside the literal consensus promoter at6.
    let promoter = assess(&molecule, &request(&molecule, "hinfi", 17, 7, false)).unwrap();
    assert_eq!(promoter.candidates[0].state, CandidateState::Rejected);
    assert!(
        promoter.candidates[0]
            .issues
            .iter()
            .any(|i| i.kind == IssueKind::PromoterOrTranscriptNotRetained)
    );
    let wrong = assess(&molecule, &request(&molecule, "smai", 17, 6, false)).unwrap();
    assert_eq!(wrong.candidates[0].state, CandidateState::Rejected);
    assert!(
        wrong.candidates[0]
            .issues
            .iter()
            .any(|i| i.kind == IssueKind::RunoffMismatch)
    );
    let multiple = source("CCCGGGTTCCCGGG", true);
    let multi = assess(&multiple, &request(&multiple, "smai", 17, 7, false)).unwrap();
    assert_eq!(multi.candidates[0].state, CandidateState::Rejected);
    assert_eq!(multi.candidates[0].cleavage.sites.len(), 2);
    assert_eq!(multi.candidates[0].fragments.len(), 2);
    assert!(
        multi.candidates[0]
            .issues
            .iter()
            .any(|i| i.kind == IssueKind::MultiplePermittedCuts)
    );
    let mut unknown = request(&molecule, "smai", 17, 7, false);
    let SourceInput::DeclaredOriginal { generation } = &mut unknown.source else {
        unreachable!()
    };
    generation.cpg = MethylationState::Unknown;
    let uncertain = assess(&molecule, &unknown).unwrap();
    assert_eq!(uncertain.candidates[0].state, CandidateState::Unresolved);
    assert!(!uncertain.candidates[0].cleavage.products_determined);
    assert!(uncertain.candidates[0].fragments.is_empty());
    assert!(uncertain.candidates[0].ivt.is_none());
    let SourceInput::DeclaredOriginal { generation } = &mut unknown.source else {
        unreachable!()
    };
    generation.cpg = MethylationState::FullyModified;
    let blocked = assess(&molecule, &unknown).unwrap();
    assert_eq!(blocked.candidates[0].state, CandidateState::Rejected);
    assert!(
        blocked.candidates[0]
            .issues
            .iter()
            .any(|i| i.kind == IssueKind::NoPermittedCut)
    );
    assert_eq!(blocked.candidates[0].fragments.len(), 1);
    assert_eq!(
        blocked.candidates[0].fragments[0].geometry.topology,
        Topology::Circular
    );
    assert!(
        blocked.candidates[0].ivt.is_none(),
        "blocked cleavage preserves its intact circular survivor, not a linearized template"
    );
}

#[test]
fn explicit_formulation_conditions_never_auto_select_hf_or_turn_geometry_into_success() {
    let molecule = source("GAATTC", false);
    let mut q = request(&molecule, "ecori", 17, 9, false);
    q.policy.require_known_conditions = true;
    let absent = assess(&molecule, &q).unwrap();
    assert_eq!(absent.candidates[0].state, CandidateState::Unresolved);
    for (formulation, state) in [
        (None, CandidateState::Unresolved),
        (Some("neb-r0101"), CandidateState::Rejected),
        (Some("neb-r3101"), CandidateState::EligibleUnderDeclarations),
    ] {
        q.conditions = Some(conditions::ReactionConditions {
            schema_version: 1,
            formulations: vec![conditions::ChosenFormulation {
                enzyme_id: "ecori".into(),
                formulation_id: formulation.map(str::to_owned),
            }],
            buffer: "rCutSmart".into(),
            temperature_milli_c: 37000,
            incubation_minutes: 60,
            minimum_activity_percent: 75,
            activity_context: Some("NEB unit-assay relative activity".into()),
            temperature_context: Some("vendor recommended incubation".into()),
        });
        let report = assess(&molecule, &q).unwrap();
        assert_eq!(report.candidates[0].state, state);
        assert_eq!(
            report.candidates[0]
                .mapping
                .as_ref()
                .unwrap()
                .template_cut_offset,
            9
        );
        assert_eq!(
            report.candidates[0]
                .conditions
                .as_ref()
                .unwrap()
                .report
                .options
                .formulations[0]
                .formulation_id
                .as_deref(),
            formulation
        );
        assert_eq!(
            report.candidates[0]
                .ivt
                .as_ref()
                .unwrap()
                .request
                .source_revision,
            molecule.revision()
        );
    }
}

#[test]
fn entire_installed_catalogue_is_admitted_and_shared_limits_do_not_declare_unsearched_absence() {
    let molecule = source("CCCGGG", true);
    let catalog = catalog::builtin_catalog().unwrap();
    let ids: Vec<_> = catalog
        .edition()
        .enzymes
        .iter()
        .map(|e| e.id.clone())
        .collect();
    assert!(ids.len() <= MAX_ENZYMES);
    let mut q = request(&molecule, "smai", 17, 7, false);
    q.selection = catalog::select(&catalog, &ids).unwrap();
    let all = assess(&molecule, &q).unwrap();
    assert_eq!(
        all.candidates
            .iter()
            .map(|c| c.enzyme_id.clone())
            .collect::<Vec<_>>(),
        ids
    );
    q.selection = catalog::select(&catalog, &["smai".into(), "ecori".into()]).unwrap();
    q.limits.max_comparisons = 1;
    let limited = assess(&molecule, &q).unwrap();
    assert_eq!(limited.comparisons, 1);
    assert_eq!(limited.candidates.len(), 2);
    for candidate in &limited.candidates {
        assert_eq!(candidate.state, CandidateState::Unresolved);
        assert!(!candidate.cleavage.complete);
        assert!(
            !candidate
                .issues
                .iter()
                .any(|i| i.kind == IssueKind::NoPermittedCut)
        );
    }
    q.limits.max_report_bytes = 1024;
    assert!(matches!(assess(&molecule, &q), Err(Error::TooLarge)));
}

#[test]
fn retained_report_and_request_pins_are_not_reusable_science_after_tampering() {
    let molecule = source("CCCGGG", false);
    let q = request(&molecule, "smai", 17, 7, false);
    let report = assess(&molecule, &q).unwrap();
    report.validate_for(&molecule, &q).unwrap();
    let mut stale = q.clone();
    stale.selection.catalog_content_hash = "0".repeat(64);
    assert!(assess(&molecule, &stale).is_err());
    let mut changed = report.clone();
    changed.candidates[0]
        .mapping
        .as_mut()
        .unwrap()
        .template_cut_offset += 1;
    assert!(changed.validate_binding(&report.context_revision).is_err());
    let mut foreign = molecule;
    foreign.name = "Another same-letter source".into();
    assert!(report.validate_for(&foreign, &q).is_err());
}

#[test]
fn saved_material_replays_the_existing_origin_instead_of_redeclaring_its_generation() {
    use crate::{
        documents::workspace::Workspace,
        material::{
            capsule,
            lineage::{EventBinding, InputEventBinding},
            outputs, run,
        },
        sequence::analysis::Primer,
    };
    let original = source("CCCGGG", false);
    let mut target = original.clone();
    target.name = "Independent unchanged PCR template".into();
    let plan = DesignPlan {
        standard: None,
        name: "Literal original-source PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: original.revision(),
            molecule: original.clone(),
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("literal first8", "TAATACGA", "").unwrap(),
            reverse: Primer::new("literal reverse last8", "AACCCGGG", "").unwrap(),
            template_strand: crate::chemistry::Strand::Top,
        }],
    };
    let workspace =
        Workspace::from_documents("Saved IVT template", vec![original, target], 1).unwrap();
    let event = |index: usize| {
        let doc = &workspace.documents[index];
        EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        }
    };
    let commit = outputs::finalize(
        &workspace,
        &plan,
        event(1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: event(0),
        }],
        &outputs::Options {
            name: "Saved actual PCR template".into(),
            ..Default::default()
        },
    )
    .unwrap();
    let saved = commit
        .workspace
        .head(&commit.receipt.construct.document_id)
        .unwrap();
    let origin = capsule::capsule(
        &commit.workspace,
        &run::Selection {
            run_id: commit.receipt.run.run_id.clone(),
            step_index: 0,
            event: commit.receipt.construct.clone(),
        },
    )
    .unwrap();
    let mut q = request(saved, "smai", 17, 7, false);
    q.source = SourceInput::Saved {
        origin: MaterialOrigin::Run(origin),
    };
    let report = assess(saved, &q).unwrap();
    assert_eq!(
        report.candidates[0].state,
        CandidateState::EligibleUnderDeclarations
    );
    assert!(
        report.source_modifications.generation.is_empty(),
        "fresh PCR source does not gain an invented original-host generation"
    );
    assert!(!report.source_history_incomplete);
    q.source = declared();
    assert!(
        assess(saved, &q).is_err(),
        "saved material cannot be relabeled as an original reference generation"
    );
}
