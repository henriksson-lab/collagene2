use collagene_antibody_numbering::{Error, profile, v_mapping::*};
use collagene_model::{Alphabet, Molecule, PROTEIN_RESIDUES, Sequence, Topology};
// OGRDB IGHV3-23*01, G00060 v10: labeled V-prefix with gaps removed.
// CC0 data provenance: data/ogrdb-human-vj/ATTRIBUTION.md.
const HEAVY: &str = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC";
fn source(s: &str) -> Molecule {
    Molecule {
        name: "V-prefix correspondence fixture".into(),
        sequence: Sequence::new(Alphabet::Protein, s).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn scoring() -> Scoring {
    Scoring {
        observed_reward: 5,
        frequency_bonus: 8,
        unseen_penalty: -6,
        gap_open: -12,
        gap_extend: -2,
        minimum_score: 300,
    }
}
fn limits() -> Limits {
    Limits {
        max_cells: profile::MAX_CELLS,
        max_endpoints_per_locus: profile::MAX_ENDPOINTS,
        max_trace_states: profile::MAX_TRACE_STATES,
    }
}
#[test]
fn distinct_proteins_and_explicit_missing_columns_preserve_actual_reference_evidence() {
    let models = models(&ReferenceSelection::human_all(), &scoring(), &mut || false).unwrap();
    assert_eq!(
        models
            .iter()
            .map(|m| (
                m.training_alleles.len(),
                m.omitted_alleles.len(),
                m.distinct_proteins,
                m.profile.states.len()
            ))
            .collect::<Vec<_>>(),
        vec![(192, 2, 140, 100), (72, 2, 60, 94), (84, 2, 72, 97)]
    );
    assert_eq!(
        models[0]
            .omitted_positions
            .iter()
            .map(|p| p.number)
            .collect::<Vec<_>>(),
        [10, 32, 33, 73]
    );
    assert_eq!(
        models[1]
            .omitted_positions
            .iter()
            .map(|p| p.number)
            .collect::<Vec<_>>(),
        [58, 59, 60, 61, 62, 63, 64, 73, 81, 82]
    );
    assert_eq!(
        models[2]
            .omitted_positions
            .iter()
            .map(|p| p.number)
            .collect::<Vec<_>>(),
        [10, 32, 33, 34, 60, 61, 73]
    );
    // Literal reference-column facts computed separately from the Rust builder.
    let w = PROTEIN_RESIDUES.iter().position(|&a| a == b'W').unwrap();
    let c = PROTEIN_RESIDUES.iter().position(|&a| a == b'C').unwrap();
    for (model, expected_w) in models.iter().zip([139, 60, 72]) {
        let i = model.positions.iter().position(|p| p.number == 41).unwrap();
        assert_eq!(model.counts[i].residues[w], expected_w);
        assert_eq!(
            model.counts[i].residues[c],
            usize::from(model.locus == Locus::Heavy)
        );
        for anchor in [23, 104] {
            let i = model
                .positions
                .iter()
                .position(|p| p.number == anchor)
                .unwrap();
            assert!(model.profile.states[i].required);
            assert_eq!(model.counts[i].residues[c], model.distinct_proteins);
            assert_eq!(model.profile.states[i].emissions[c], Some(13));
            assert_eq!(
                model.profile.states[i]
                    .emissions
                    .iter()
                    .filter(|e| e.is_some())
                    .count(),
                1
            );
        }
        assert!(
            model
                .counts
                .iter()
                .all(|c| c.gaps + c.residues.iter().sum::<usize>() == model.distinct_proteins)
        );
    }
    // Heavy position27: G137/Q1/W2, all140 distinct sequences represented once.
    let i = models[0]
        .positions
        .iter()
        .position(|p| p.number == 27)
        .unwrap();
    let g = PROTEIN_RESIDUES.iter().position(|&a| a == b'G').unwrap();
    let q = PROTEIN_RESIDUES.iter().position(|&a| a == b'Q').unwrap();
    assert_eq!(models[0].counts[i].residues[g], 137);
    assert_eq!(models[0].counts[i].residues[q], 1);
    assert_eq!(models[0].profile.states[i].emissions[g], Some(12));
    assert_eq!(models[0].profile.states[i].emissions[q], Some(5));
    assert_eq!(models[0].profile.states[i].emissions[c], Some(-6));
    models[0].validate().unwrap();
    let mut forged = models[0].clone();
    forged.positions[0].number = 2;
    assert_eq!(forged.validate(), Err(Error::Binding));
}
#[test]
fn embedded_reference_query_keeps_original_source_indices_and_all_locus_evidence() {
    let prefix = "MHHHHHH*";
    let m = source(&format!("{prefix}{HEAVY}*GGGG"));
    let report = map_query(
        &m,
        &ReferenceSelection::human_all(),
        &scoring(),
        &limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(
        report
            .mappings
            .iter()
            .map(|m| m.locus.clone())
            .collect::<Vec<_>>(),
        [Locus::Heavy, Locus::Kappa, Locus::Lambda]
    );
    assert_eq!(report.cells, 294 * (m.sequence.len() + 1));
    let heavy = &report.mappings[0];
    assert_eq!(
        heavy.alignment.endpoints[0].upper.end,
        prefix.len() + HEAVY.len()
    );
    let projection = heavy.projections[0]
        .projection
        .as_ref()
        .expect("literal reference-prefix mapping should have a unique scored path");
    for position in [23, 104] {
        let matched = projection
            .matches
            .iter()
            .find(|r| r.reference_position.number == position)
            .unwrap();
        assert_eq!(matched.residue, 'C');
        assert_eq!(m.sequence.as_str().as_bytes()[matched.source_index], b'C');
    }
    assert!(
        projection
            .matches
            .iter()
            .all(|m| m.source_index >= prefix.len() && m.source_index < prefix.len() + HEAVY.len())
    );
    assert!(
        projection
            .matches
            .windows(2)
            .all(|w| w[0].source_index < w[1].source_index)
    );
    assert_eq!(
        projection.matches.len() + projection.unassigned_residues.len(),
        heavy.alignment.endpoints[0].upper.end - heavy.alignment.endpoints[0].upper.start
    );
    report
        .validate_for(&m, &ReferenceSelection::human_all(), &scoring(), &limits())
        .unwrap();
}
#[test]
fn tandem_endpoints_are_retained_and_ambiguity_and_score_filters_remain_explicit() {
    let repeated = map_query(
        &source(&format!("{HEAVY}*{HEAVY}")),
        &ReferenceSelection::human_all(),
        &scoring(),
        &limits(),
        &mut || false,
    )
    .unwrap();
    let heavy = &repeated.mappings[0];
    assert_eq!(heavy.alignment.qualified_endpoints, 2);
    assert!(!heavy.alignment.truncated);
    assert_eq!(
        heavy.projections.iter().map(|p| p.end).collect::<Vec<_>>(),
        [HEAVY.len(), 2 * HEAVY.len() + 1]
    );
    assert!(
        heavy
            .projections
            .iter()
            .all(|p| p.projection.is_some() && p.projection_issue.is_none())
    );
    let mut capped = limits();
    capped.max_endpoints_per_locus = 1;
    let limited = map_query(
        &source(&format!("{HEAVY}*{HEAVY}")),
        &ReferenceSelection::human_all(),
        &scoring(),
        &capped,
        &mut || false,
    )
    .unwrap();
    assert_eq!(limited.mappings[0].alignment.qualified_endpoints, 2);
    assert!(limited.mappings[0].alignment.truncated);
    assert_eq!(limited.mappings[0].projections.len(), 1);
    capped.max_endpoints_per_locus = 0;
    let counted = map_query(
        &source(HEAVY),
        &ReferenceSelection::human_all(),
        &scoring(),
        &capped,
        &mut || false,
    )
    .unwrap();
    assert_eq!(counted.mappings[0].alignment.qualified_endpoints, 1);
    assert!(counted.mappings[0].projections.is_empty());
    assert!(counted.mappings[0].alignment.truncated);
    let ambiguous = map_query(
        &source(&format!("X*{HEAVY}")),
        &ReferenceSelection::human_all(),
        &scoring(),
        &limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(
        ambiguous.mappings[0].projections[0].projection_issue,
        Some(ProjectionIssue::AmbiguousSource)
    );
    let mut high = scoring();
    high.minimum_score = i32::MAX;
    let filtered = map_query(
        &source(HEAVY),
        &ReferenceSelection::human_all(),
        &high,
        &limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        filtered
            .mappings
            .iter()
            .all(|m| m.alignment.qualified_endpoints == 0 && m.projections.is_empty())
    );
    assert!(filtered.mappings[0].alignment.reachable_endpoints > 0);
    let absent = map_query(
        &source("AAACAAA"),
        &ReferenceSelection::human_all(),
        &scoring(),
        &limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        absent
            .mappings
            .iter()
            .all(|m| m.alignment.reachable_endpoints == 0)
    );
}
#[test]
fn reports_reconstruct_models_source_and_limits_and_jobs_cancel_without_partial_success() {
    let m = source(HEAVY);
    let report = map_query(
        &m,
        &ReferenceSelection::human_all(),
        &scoring(),
        &limits(),
        &mut || false,
    )
    .unwrap();
    let round: Report = serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    let actual_models =
        models(&ReferenceSelection::human_all(), &scoring(), &mut || false).unwrap();
    for (mapping, model) in report.mappings.iter().zip(&actual_models) {
        assert_eq!(mapping.training_alleles, model.training_alleles);
        assert_eq!(mapping.omitted_alleles, model.omitted_alleles);
    }
    let mut altered = round.clone();
    altered.mappings[0].training_alleles.clear();
    assert!(
        altered
            .validate_for(&m, &ReferenceSelection::human_all(), &scoring(), &limits())
            .is_err()
    );
    round
        .validate_for(&m, &ReferenceSelection::human_all(), &scoring(), &limits())
        .unwrap();
    let mut forged = report.clone();
    forged.mappings[0].model_revision = "forged".into();
    assert_eq!(
        forged.validate_for(&m, &ReferenceSelection::human_all(), &scoring(), &limits()),
        Err(Error::Binding)
    );
    let mut renamed = m.clone();
    renamed.name = "Changed source metadata".into();
    assert_eq!(
        report.validate_for(
            &renamed,
            &ReferenceSelection::human_all(),
            &scoring(),
            &limits()
        ),
        Err(Error::Binding)
    );
    assert!(
        map_query(
            &m,
            &ReferenceSelection::human_all(),
            &scoring(),
            &Limits {
                max_cells: report.cells - 1,
                max_endpoints_per_locus: profile::MAX_ENDPOINTS,
                max_trace_states: profile::MAX_TRACE_STATES,
            },
            &mut || false
        )
        .is_err()
    );
    assert!(
        map_query(
            &m,
            &ReferenceSelection::human_all(),
            &scoring(),
            &Limits {
                max_cells: 0,
                max_endpoints_per_locus: profile::MAX_ENDPOINTS,
                max_trace_states: profile::MAX_TRACE_STATES,
            },
            &mut || false
        )
        .is_err()
    );
    let mut invalid = scoring();
    invalid.frequency_bonus = 501;
    assert!(models(&ReferenceSelection::human_all(), &invalid, &mut || false).is_err());
    let mut calls = 0;
    assert_eq!(
        map_query(
            &m,
            &ReferenceSelection::human_all(),
            &scoring(),
            &limits(),
            &mut || {
                calls += 1;
                calls == 10
            }
        )
        .unwrap_err(),
        Error::Cancelled
    );
    assert!(calls >= 10);
    let dna = Molecule::dna("DNA", "ACGT", Topology::Linear).unwrap();
    assert!(
        map_query(
            &dna,
            &ReferenceSelection::human_all(),
            &scoring(),
            &limits(),
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn query_insertions_remain_unassigned_and_forged_correspondences_fail_reconstruction() {
    // O is an explicit protein residue, absent from this reference corpus. A
    // long interior run tests path insertions without guessing CDR label rules.
    let split = 45;
    let extra = "O".repeat(20);
    let m = source(&format!("{}{}{}", &HEAVY[..split], extra, &HEAVY[split..]));
    let report = map_query(
        &m,
        &ReferenceSelection::human_all(),
        &scoring(),
        &limits(),
        &mut || false,
    )
    .unwrap();
    let projection = report.mappings[0].projections[0]
        .projection
        .as_ref()
        .unwrap();
    assert_eq!(
        projection.unassigned_residues,
        (split..split + 20).collect::<Vec<_>>()
    );
    assert!(
        !projection
            .matches
            .iter()
            .any(|m| (split..split + 20).contains(&m.source_index))
    );
    let mut forged = report.clone();
    forged.mappings[0].projections[0]
        .projection
        .as_mut()
        .unwrap()
        .matches[0]
        .source_index += 1;
    assert_eq!(
        forged.validate_for(&m, &ReferenceSelection::human_all(), &scoring(), &limits()),
        Err(Error::Binding)
    );
    let mut forged = report.clone();
    forged.mappings[0].projections[0]
        .projection
        .as_mut()
        .unwrap()
        .unassigned_residues
        .clear();
    assert_eq!(
        forged.validate_for(&m, &ReferenceSelection::human_all(), &scoring(), &limits()),
        Err(Error::Binding)
    );
}

#[test]
fn fixed_terminal_anchor_survives_tied_internal_correspondences_without_projecting_them() {
    let panel: serde_json::Value =
        serde_json::from_str(include_str!("../data/rcsb-assessment/panel.json")).unwrap();
    let row = panel["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "1HZH:1")
        .unwrap();
    let m = source(row["sequence"].as_str().unwrap());
    let r = map_query(
        &m,
        &ReferenceSelection::human_all(),
        &scoring(),
        &limits(),
        &mut || false,
    )
    .unwrap();
    let heavy = &r.mappings[0];
    let end = heavy.projections.iter().find(|p| p.end == 96).unwrap();
    assert_eq!(
        end.projection_issue,
        Some(ProjectionIssue::TiedOptimalPaths)
    );
    assert!(end.projection.is_none());
    let anchor = end
        .invariant_correspondences
        .as_ref()
        .unwrap()
        .iter()
        .find(|p| p.reference_position.number == 104)
        .unwrap();
    assert_eq!(anchor.source_index, 95);
    assert_eq!(anchor.residue, 'C');
    assert_eq!(anchor.reference_position.number, 104);
    let model = models(&ReferenceSelection::human_all(), &scoring(), &mut || false)
        .unwrap()
        .remove(0);
    let last = model.profile.states.last().unwrap();
    assert!(last.required && !last.allow_insert_after);
    assert_eq!(model.positions.last().unwrap().number, 104);
    let mut forged = r.clone();
    forged.mappings[0]
        .projections
        .iter_mut()
        .find(|p| p.end == 96)
        .unwrap()
        .invariant_correspondences
        .as_mut()
        .unwrap()
        .iter_mut()
        .find(|p| p.reference_position.number == 104)
        .unwrap()
        .source_index += 1;
    assert_eq!(
        forged.validate_for(&m, &ReferenceSelection::human_all(), &scoring(), &limits()),
        Err(Error::Binding)
    );
}
