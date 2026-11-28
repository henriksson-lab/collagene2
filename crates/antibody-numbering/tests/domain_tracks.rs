use collagene_antibody_numbering::{
    Error, Insertion, Position, domain_candidates as d, domain_tracks::*, profile, v_mapping as v,
};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
const V: &str = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC";
const J: &str = "WGQGTLVTVSS";
fn source(sequence: &str) -> Molecule {
    Molecule {
        name: "Authored conditional track fixture".into(),
        sequence: Sequence::new(Alphabet::Protein, sequence).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn options() -> d::Options {
    d::Options {
        references: collagene_antibody_numbering::v_mapping::ReferenceSelection::human_all(),
        v_scoring: v::Scoring {
            observed_reward: 5,
            frequency_bonus: 8,
            unseen_penalty: -6,
            gap_open: -12,
            gap_extend: -2,
            minimum_score: 300,
        },
        min_cdr3_residues: 2,
        max_cdr3_residues: 60,
        max_fr4_mismatches: 1,
    }
}
fn limits() -> d::Limits {
    d::Limits {
        v_mapping: v::Limits {
            max_cells: profile::MAX_CELLS,
            max_endpoints_per_locus: profile::MAX_ENDPOINTS,
            max_trace_states: profile::MAX_TRACE_STATES,
        },
        max_scan_windows: d::MAX_SCAN_WINDOWS,
        max_candidates: d::MAX_CANDIDATES,
    }
}
fn label(s: &str) -> Position {
    s.parse().unwrap()
}
#[test]
fn region_definition_has_literal_boundaries_and_preserves_exact_insertion_labels() {
    for (range, region) in [
        (1..=26, Region::Fr1),
        (27..=38, Region::Cdr1),
        (39..=55, Region::Fr2),
        (56..=65, Region::Cdr2),
        (66..=104, Region::Fr3),
        (105..=117, Region::Cdr3),
        (118..=128, Region::Fr4),
    ] {
        for n in range {
            assert_eq!(imgt_region(&label(&n.to_string())).unwrap(), region);
        }
    }
    for s in ["111.1", "111.10", "112.1", "112.10"] {
        let p = label(s);
        assert_eq!(imgt_region(&p).unwrap(), Region::Cdr3);
        assert_eq!(p.to_string(), s);
    }
    for s in ["129", "104.1", "111A", "112AA", "118.1"] {
        assert!(imgt_region(&label(s)).is_err());
    }
    assert!(serde_json::from_str::<NumberingScheme>("\"Kabat\"").is_err());
    assert!(serde_json::from_str::<RegionDefinition>("\"Chothia\"").is_err());
}
#[test]
fn authored_track_preserves_source_and_distinguishes_framework_loop_and_j_evidence() {
    let m = source(&format!("{V}ARDYGMDV{J}"));
    let r = analyze(&m, &options(), &limits(), &selection(), &mut || false).unwrap();
    assert_eq!(r.tracks.len(), 1);
    let t = &r.tracks[0];
    assert_eq!(t.candidate_id, r.evidence.candidates[0].id);
    assert_eq!(t.source_end, m.sequence.len());
    assert_eq!(
        t.residues.len(),
        t.source_end - t.first_supported_source_index
    );
    for row in &t.residues {
        assert_eq!(
            row.residue,
            m.sequence.as_str().as_bytes()[row.source_index] as char
        );
    }
    for (index, position, region, evidence) in [
        (
            21,
            "23",
            Region::Fr1,
            AssignmentEvidence::InvariantVReference,
        ),
        (
            25,
            "27",
            Region::Cdr1,
            AssignmentEvidence::ConditionalLoopRule,
        ),
        (
            32,
            "38",
            Region::Cdr1,
            AssignmentEvidence::ConditionalLoopRule,
        ),
        (
            50,
            "56",
            Region::Cdr2,
            AssignmentEvidence::ConditionalLoopRule,
        ),
        (
            57,
            "65",
            Region::Cdr2,
            AssignmentEvidence::ConditionalLoopRule,
        ),
        (
            V.len(),
            "105",
            Region::Cdr3,
            AssignmentEvidence::ConditionalLoopRule,
        ),
        (
            V.len() + 8,
            "118",
            Region::Fr4,
            AssignmentEvidence::JFr4Reference,
        ),
    ] {
        let row = t.residues.iter().find(|r| r.source_index == index).unwrap();
        assert_eq!(
            row.assignment,
            Some(Assignment {
                position: label(position),
                evidence
            })
        );
        assert_eq!(row.region, Some(region));
    }
    assert_eq!(
        t.unassigned_residues,
        t.residues.iter().filter(|r| r.assignment.is_none()).count()
    );
    r.validate_for(&m, &options(), &limits(), &selection())
        .unwrap();
}
#[test]
fn tied_observed_paths_keep_unresolved_rows_without_a_representative_fallback() {
    let panel: serde_json::Value =
        serde_json::from_str(include_str!("../data/rcsb-assessment/panel.json")).unwrap();
    let case = panel["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "1HZH:2")
        .unwrap();
    let m = source(case["sequence"].as_str().unwrap());
    let r = analyze(&m, &options(), &limits(), &selection(), &mut || false).unwrap();
    let t = r
        .tracks
        .iter()
        .find(|t| t.locus == v::Locus::Kappa)
        .unwrap();
    let c = r
        .evidence
        .candidates
        .iter()
        .find(|c| c.id == t.candidate_id)
        .unwrap();
    assert!(matches!(c.cdr1, d::LoopAssignment::Unresolved { .. }));
    assert!(t.unassigned_residues > 0);
    let p = r
        .evidence
        .v_evidence
        .mappings
        .iter()
        .find(|m| m.locus == t.locus)
        .unwrap()
        .projections
        .iter()
        .find(|p| p.end == c.v_endpoint)
        .unwrap();
    assert!(p.projection.is_none());
    for row in &t.residues {
        if let Some(a) = &row.assignment {
            assert_ne!(
                row.region,
                Some(Region::Cdr1),
                "unresolved loop must not inherit reference storage labels"
            );
            if a.evidence == AssignmentEvidence::InvariantVReference {
                assert!(
                    p.invariant_correspondences
                        .as_ref()
                        .unwrap()
                        .iter()
                        .any(|m| m.source_index == row.source_index
                            && m.reference_position == a.position)
                );
            }
        }
    }
}
#[test]
fn unsupported_loop_is_unlabeled_and_does_not_remove_the_remaining_track() {
    let m = source(&format!(
        "{}{}{}ARDYGMDV{J}",
        &V[..27],
        "O".repeat(20),
        &V[27..]
    ));
    let r = analyze(&m, &options(), &limits(), &selection(), &mut || false).unwrap();
    let c = r
        .evidence
        .candidates
        .iter()
        .find(|c| c.locus == v::Locus::Heavy)
        .unwrap();
    assert!(matches!(c.cdr1, d::LoopAssignment::Unresolved { .. }));
    let t = r.tracks.iter().find(|t| t.candidate_id == c.id).unwrap();
    assert!(t.unassigned_residues >= 20);
    assert!(t.residues.iter().all(|r| r.region != Some(Region::Cdr1)));
    assert!(t.residues.iter().any(|r| r.region == Some(Region::Cdr3)));
}
#[test]
fn tandem_tracks_keep_absolute_offsets_and_long_cdr3_labels_in_source_order() {
    let block = format!("{V}{}{J}", "A".repeat(40));
    let linker = "GGGGSGGGGS";
    let m = source(&format!("{block}{linker}{block}"));
    let r = analyze(&m, &options(), &limits(), &selection(), &mut || false).unwrap();
    for offset in [0, block.len() + linker.len()] {
        let t = r
            .tracks
            .iter()
            .find(|t| t.source_end == offset + block.len() && t.locus == v::Locus::Heavy)
            .unwrap();
        assert!(t.residues.iter().all(|r| r.source_index >= offset));
        let loop_rows = t
            .residues
            .iter()
            .filter(|r| r.region == Some(Region::Cdr3))
            .collect::<Vec<_>>();
        assert_eq!(loop_rows.len(), 40);
        assert_eq!(loop_rows[0].source_index, offset + V.len());
        let labels = loop_rows
            .iter()
            .map(|r| r.assignment.as_ref().unwrap().position.to_string())
            .collect::<Vec<_>>();
        assert!(labels.contains(&"111.10".into()) && labels.contains(&"112.10".into()));
        assert!(
            labels.iter().position(|s| s == "112.10") < labels.iter().position(|s| s == "112.1")
        );
        assert!(loop_rows.iter().any(|r| matches!(
            r.assignment.as_ref().unwrap().position.insertion,
            Insertion::Decimal(_)
        )));
    }
}
#[test]
fn truncated_n_terminus_stays_unresolved_and_track_construction_is_cancellable() {
    let m = source(&format!("{}ARDYGMDV{J}", &V[8..]));
    let r = analyze(&m, &options(), &limits(), &selection(), &mut || false).unwrap();
    let t = r
        .tracks
        .iter()
        .find(|t| t.locus == v::Locus::Heavy)
        .unwrap();
    assert_eq!(t.invariant_position_one, None);
    assert!(
        t.residues
            .iter()
            .filter_map(|r| r.assignment.as_ref())
            .all(|a| a.position.number != 1)
    );
    assert!(
        t.residues
            .iter()
            .filter_map(|r| r.assignment.as_ref())
            .any(|a| a.position.number == 23)
    );
    let mut candidate_calls = 0;
    d::find_candidates(&m, &options(), &limits(), &mut || {
        candidate_calls += 1;
        false
    })
    .unwrap();
    let mut track_calls = 0;
    assert_eq!(
        analyze(&m, &options(), &limits(), &selection(), &mut || {
            track_calls += 1;
            track_calls > candidate_calls + 5
        }),
        Err(Error::Cancelled)
    );
    assert!(track_calls > candidate_calls);
}
#[test]
fn reports_bind_source_options_context_rows_and_resource_scope() {
    let m = source(&format!("{V}ARDYGMDV{J}"));
    let r = analyze(&m, &options(), &limits(), &selection(), &mut || false).unwrap();
    let mut forged = r.clone();
    forged.tracks[0].residues[0].assignment = None;
    assert!(
        forged
            .validate_for(&m, &options(), &limits(), &selection())
            .is_err()
    );
    let mut changed = m.clone();
    changed.name.push('x');
    assert!(
        r.validate_for(&changed, &options(), &limits(), &selection())
            .is_err()
    );
    let mut altered = options();
    altered.max_fr4_mismatches = 0;
    assert!(
        r.validate_for(&m, &altered, &limits(), &selection())
            .is_err()
    );
    let mut forged = r.clone();
    forged.engine.push('x');
    assert_eq!(
        forged.validate_for(&m, &options(), &limits(), &selection()),
        Err(Error::Binding)
    );
    assert_eq!(
        analyze(&m, &options(), &limits(), &selection(), &mut || true),
        Err(Error::Cancelled)
    );
    let mut cap = limits();
    cap.max_candidates = 0;
    let counts = analyze(&m, &options(), &cap, &selection(), &mut || false).unwrap();
    assert!(counts.tracks.is_empty() && counts.evidence.pair_retention_truncated);
    assert!(counts.evidence.qualified_pairs_in_retained_v > 0);
}

fn selection() -> collagene_antibody_numbering::domain_tracks::Selection {
    collagene_antibody_numbering::domain_tracks::Selection {
        numbering_scheme: collagene_antibody_numbering::domain_tracks::NumberingScheme::Imgt,
        region_definition: collagene_antibody_numbering::domain_tracks::RegionDefinition::Imgt,
    }
}
