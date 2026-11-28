//! Synthetic conditional routing fixtures, not a species-classifier calibration.
use collagene_antibody_numbering::{
    domain_candidates as d, domain_tracks as t, profile,
    references::{self, Corpus, SourcePolicy, Species},
    v_mapping::{self, Locus, ReferenceSelection},
};
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
const V: &str = "QIQLQQSGPELVKPGASVKISCKASGYTFTDYYINWVKQRPGQGLEWIGWIYPGSGNTKYNEKFKGKATLTVDTSSSTAYMQLSSLTSEDSAVYFC";
const J: &str = "WGTGTTVTVSS";
fn reference() -> ReferenceSelection {
    ReferenceSelection {
        species: Species::MouseC57Bl6,
        loci: vec![Locus::Heavy],
    }
}
fn scoring() -> v_mapping::Scoring {
    v_mapping::Scoring {
        observed_reward: 5,
        frequency_bonus: 8,
        unseen_penalty: -6,
        gap_open: -12,
        gap_extend: -2,
        minimum_score: 300,
    }
}
fn options() -> d::Options {
    d::Options {
        references: reference(),
        v_scoring: scoring(),
        min_cdr3_residues: 2,
        max_cdr3_residues: 64,
        max_fr4_mismatches: 0,
    }
}
fn limits() -> d::Limits {
    d::Limits {
        v_mapping: v_mapping::Limits {
            max_cells: profile::MAX_CELLS,
            max_endpoints_per_locus: profile::MAX_ENDPOINTS,
            max_trace_states: profile::MAX_TRACE_STATES,
        },
        max_scan_windows: d::MAX_SCAN_WINDOWS,
        max_candidates: d::MAX_CANDIDATES,
    }
}
fn source() -> Molecule {
    Molecule {
        name: "Authored mouse V plus ARDY plus observed J FR4".into(),
        sequence: Sequence::new(Alphabet::Protein, &format!("{V}ARDY{J}")).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn selection() -> t::Selection {
    t::Selection {
        numbering_scheme: t::NumberingScheme::Imgt,
        region_definition: t::RegionDefinition::Imgt,
    }
}
#[test]
fn pinned_mouse_facts_keep_atypical_functional_alleles_and_explicit_j_label_origin() {
    let corpus = references::load(&Species::MouseC57Bl6).unwrap();
    assert!(std::ptr::eq(
        corpus,
        references::load(&Species::MouseC57Bl6).unwrap()
    ));
    assert_eq!(
        (
            corpus.sources().len(),
            corpus.v_entries().len(),
            corpus.j_entries().len()
        ),
        (1, 95, 4)
    );
    assert_eq!(corpus.sources()[0].set_ref, "OGRDB:C57BL/6 IGH.5");
    assert_eq!(corpus.sources()[0].excluded.len(), 15);
    let SourcePolicy::MousePinnedColumns {
        atypical_anchors,
        j_phase_policy,
        ..
    } = corpus.source_policy()
    else {
        panic!("explicit mouse provenance")
    };
    assert_eq!(atypical_anchors.len(), 4);
    assert!(j_phase_policy.contains("differs"));
    let v = corpus
        .v_entries()
        .iter()
        .find(|v| v.label == "IGHV-24BS")
        .unwrap();
    assert_eq!(v.aligned_protein.replace('-', ""), V);
    for j in corpus.j_entries() {
        assert_eq!(
            j.fr4_positions,
            (118..=128).map(|n| n.to_string()).collect::<Vec<_>>()
        );
        assert_eq!(j.anchor_offset % 3, usize::from(j.provider_codon_frame - 1));
        assert_ne!(
            j.anchor_offset % 3,
            usize::from((4 - j.provider_codon_frame) % 3)
        );
        assert_eq!(j.trailing_nucleotides, "G");
    }
    assert_eq!(
        corpus
            .j_entries()
            .iter()
            .find(|j| j.label == "IGHJ-32C2")
            .unwrap()
            .fr4_protein,
        J
    );
    assert_ne!(
        corpus.revision().unwrap(),
        references::load(&Species::Human)
            .unwrap()
            .revision()
            .unwrap()
    );
}
#[test]
fn model_admits_92_canonical_mouse_references_without_repair_or_species_relabeling() {
    let corpus = references::load(&Species::MouseC57Bl6).unwrap();
    let models = v_mapping::models(&reference(), &scoring(), &mut || false).unwrap();
    assert_eq!(models.len(), 1);
    let model = &models[0];
    assert_eq!(model.species, Species::MouseC57Bl6);
    assert_eq!(model.locus, Locus::Heavy);
    assert_eq!(
        (
            model.training_alleles.len(),
            model.omitted_alleles.len(),
            model.distinct_proteins
        ),
        (92, 3, 92)
    );
    for name in ["IGHV-ZPQB", "IGHV-IQE5", "IGHV-OSZZ"] {
        let entry = corpus.v_entries().iter().find(|v| v.label == name).unwrap();
        assert!(model.omitted_alleles.contains(&entry.allele_ref));
        assert!(!model.training_alleles.contains(&entry.allele_ref));
    }
    let w41 = corpus
        .v_entries()
        .iter()
        .find(|v| v.label == "IGHV-IQE5")
        .unwrap();
    assert!(w41.non_cysteine_anchors.is_empty());
    assert_eq!(w41.aligned_protein.as_bytes()[40], b'C');
    model.validate().unwrap();
    let human =
        v_mapping::models(&ReferenceSelection::human_all(), &scoring(), &mut || false).unwrap();
    let subset = v_mapping::models(
        &ReferenceSelection {
            species: Species::Human,
            loci: vec![Locus::Kappa],
        },
        &scoring(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(subset[0], human[1]);
}
#[test]
fn unsupported_mouse_loci_and_unadmitted_scheme_regions_are_explicit_errors() {
    for loci in [
        vec![],
        vec![Locus::Kappa],
        vec![Locus::Heavy, Locus::Lambda],
        vec![Locus::Heavy, Locus::Heavy],
    ] {
        assert!(
            ReferenceSelection {
                species: Species::MouseC57Bl6,
                loci
            }
            .validate()
            .is_err()
        );
    }
    for chosen in [
        t::Selection {
            numbering_scheme: t::NumberingScheme::Kabat1991,
            region_definition: t::RegionDefinition::Imgt,
        },
        t::Selection {
            numbering_scheme: t::NumberingScheme::Imgt,
            region_definition: t::RegionDefinition::North2011,
        },
    ] {
        assert!(t::request_revision(&source(), &options(), &limits(), &chosen).is_err());
    }
}
#[test]
fn synthetic_mouse_heavy_route_and_current_bindings_do_not_certify_species_or_other_schemes() {
    let molecule = source();
    let options = options();
    let limits = limits();
    let selected = selection();
    let request = t::request_revision(&molecule, &options, &limits, &selected).unwrap();
    let result = t::analyze(&molecule, &options, &limits, &selected, &mut || false).unwrap();
    result.validate_binding(&request).unwrap();
    assert!(!result.tracks.is_empty());
    assert!(
        result
            .tracks
            .iter()
            .all(|track| track.locus == Locus::Heavy)
    );
    assert_eq!(result.evidence.v_evidence.selection, reference());
    let admission = &result.evidence.v_evidence.mappings[0];
    assert_eq!(
        (
            admission.training_alleles.len(),
            admission.omitted_alleles.len()
        ),
        (92, 3)
    );
    let models = v_mapping::models(&reference(), &scoring(), &mut || false).unwrap();
    assert_eq!(admission.training_alleles, models[0].training_alleles);
    assert_eq!(admission.omitted_alleles, models[0].omitted_alleles);
    let mut removed = result.clone();
    removed.evidence.v_evidence.mappings[0]
        .omitted_alleles
        .clear();
    assert!(removed.validate_binding(&request).is_err());
    let mut human = options.clone();
    human.references = ReferenceSelection::human_all();
    assert_ne!(
        request,
        t::request_revision(&molecule, &human, &limits, &selected).unwrap()
    );
    let mut forged = result;
    forged.evidence.options.references = ReferenceSelection::human_all();
    assert!(forged.validate_binding(&request).is_err());
}
#[test]
fn corpus_species_policy_and_authored_fr4_labels_cannot_be_silently_rebound() {
    let corpus = references::load(&Species::MouseC57Bl6).unwrap();
    let mut value = serde_json::to_value(corpus).unwrap();
    value["species"] = "Human".into();
    assert!(
        serde_json::from_value::<Corpus>(value)
            .unwrap()
            .revision()
            .is_err()
    );
    let mut value = serde_json::to_value(corpus).unwrap();
    value["j_entries"][0]["fr4_positions"][0] = "119".into();
    assert!(
        serde_json::from_value::<Corpus>(value)
            .unwrap()
            .revision()
            .is_err()
    );
    let mut value = serde_json::to_value(corpus).unwrap();
    value["v_entries"][0]["aligned_nucleotides"] = "é".repeat(156).into();
    assert!(
        serde_json::from_value::<Corpus>(value)
            .unwrap()
            .revision()
            .is_err()
    );
}
