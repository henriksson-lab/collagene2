use collagene_core::{
    cloning::{
        assembly::standards,
        pcr::{collection_binding, primers::ThermoConditions},
    },
    construct_qc,
    protocol::{
        execution::process_request,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult, WORKER_SCHEMA_VERSION},
    },
};
use collagene_model::{Feature, Location, Molecule, Sequence, Topology};
use collagene_oligo_binding::{BudgetLimits, Collection, Oligo, Options};

#[test]
fn public_perfect_duplex_model_rejects_non_dna_before_dinucleotide_lookup() {
    let conditions = ThermoConditions::default();
    for sequence in ["AU", "ACGU", "NN", "acgt", "é", "A"] {
        assert!(
            conditions.perfect_duplex_tm(sequence).is_err(),
            "{sequence}"
        );
    }
    assert!(conditions.perfect_duplex_tm(&"A".repeat(201)).is_err());
    let estimate = conditions.perfect_duplex_tm("ACGTACGT").unwrap();
    estimate.validate_for_conditions(&conditions).unwrap();
    let mut changed = conditions;
    changed.sodium_molar = 0.1;
    assert!(estimate.validate_for_conditions(&changed).is_err());
}

fn roundtrip(revision: String, job: AnalysisJob) -> Response {
    let request = Request {
        id: 91,
        revision: revision.clone(),
        job,
    };
    let bytes = request.encode().unwrap();
    assert_eq!(u16::from_le_bytes([bytes[4], bytes[5]]), 133);
    let output = process_request(&bytes).unwrap();
    assert_eq!(
        Response::admitted_frame_limit(&output, output.len()).unwrap(),
        collagene_core::formats::sequence::MAX_IMPORT_BYTES
    );
    let response = Response::decode(&output).unwrap();
    assert_eq!(response.id, 91);
    assert_eq!(response.revision, revision);
    let mut old = bytes;
    old[4..6].copy_from_slice(&(WORKER_SCHEMA_VERSION - 1).to_le_bytes());
    assert!(process_request(&old).is_err());
    response
}
#[test]
fn qc_worker_retains_declared_annotation_scope_and_independent_rule_results() {
    let mut m = Molecule::dna("coding", "ATGAAATAA", Topology::Linear).unwrap();
    m.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, 9),
        qualifiers: vec![],
        needs_review: false,
    });
    let mut options = construct_qc::Options::default();
    let unreviewed = construct_qc::assess(&m, &options).unwrap();
    assert!(
        unreviewed
            .findings
            .iter()
            .any(|f| f.rule == construct_qc::Rule::AnnotationBasis
                && f.verdict == construct_qc::Verdict::Indeterminate)
    );
    options
        .annotations
        .push(construct_qc::AnnotationAssumption {
            feature_index: 0,
            source_revision: m.revision(),
            annotation_revision: construct_qc::annotation_revision(&m.features[0]).unwrap(),
            rationale: "Declared CDS coordinates, not a biological identity assertion".into(),
        });
    let revision = construct_qc::revision(&m, &options).unwrap();
    let response = roundtrip(
        revision.clone(),
        AnalysisJob::AnalyzeConstructQc {
            molecule: m.clone(),
            options: options.clone(),
        },
    );
    let AnalysisResult::ConstructQc(report) = response.result else {
        panic!("QC report")
    };
    report.validate_for(&m, &options).unwrap();
    assert_eq!(report.statement, construct_qc::STATEMENT);
    assert!(report.findings.iter().any(
        |f| f.rule == construct_qc::Rule::CdsFrame && f.verdict == construct_qc::Verdict::Pass
    ));
    let mut other = options;
    other.genetic_code = 2;
    assert!(report.validate_for(&m, &other).is_err());
    assert!(
        Request {
            id: 1,
            revision: "0".repeat(64),
            job: AnalysisJob::AnalyzeConstructQc {
                molecule: m,
                options: other
            }
        }
        .encode()
        .is_err()
    );
}
fn standard() -> standards::Request {
    let profile = standards::catalogue().unwrap().remove(0);
    let slot = profile.definition.grammar.slot("noncoding").unwrap();
    let sequence = format!("{}AACCGGTT{}", slot.left_addition, slot.right_addition);
    standards::Request {
        profile_id: profile.profile_id,
        profile_identity: profile.profile_identity,
        parts: vec![standards::SourcePart {
            id: "part".into(),
            molecule: Molecule::dna("part", &sequence, Topology::Linear).unwrap(),
            start: 0,
            length: sequence.len(),
            reverse: false,
            slot: "noncoding".into(),
            body_start: 22,
            body_end: 30,
        }],
        limits: Default::default(),
    }
}
#[test]
fn standard_worker_pins_profile_source_crop_and_full_part_format_without_assembly_claim() {
    let request = standard();
    let revision = standards::revision(&request).unwrap();
    let response = roundtrip(
        revision,
        AnalysisJob::InspectAssemblyStandard {
            request: request.clone(),
        },
    );
    let AnalysisResult::AssemblyStandard(report) = response.result else {
        panic!("standard report")
    };
    report.validate_for(&request).unwrap();
    assert!(report.assessment.complete);
    assert_eq!(
        report.assessment.outcome,
        collagene_assembly_standards::assessment::Outcome::Conforms
    );
    assert_eq!(
        report.sources[0].map_interval(22, 30).unwrap(),
        vec![(22, 30)]
    );
    let mut stale = request;
    stale.profile_identity = "0".repeat(64);
    assert!(standards::assess(&stale).is_err());
    let wrapped = standards::SourceBinding {
        id: "wrapped".into(),
        revision: "a".repeat(64),
        start: 8,
        length: 5,
        reverse: true,
        source_length: 10,
        topology: Topology::Circular,
    };
    assert_eq!(wrapped.map_interval(0, 4).unwrap(), vec![(9, 10), (0, 3)]);
    let mut malformed = standard();
    malformed.parts[0].molecule.features.push(Feature {
        kind: "CDS".into(),
        location: Location::range(0, usize::MAX),
        qualifiers: vec![],
        needs_review: false,
    });
    assert!(standards::assess(&malformed).is_err());
}
fn binding() -> collection_binding::Request {
    collection_binding::Request {
        source: Molecule::dna("template", "ACGTACGTACGTACGT", Topology::Circular).unwrap(),
        collection: Collection {
            oligos: vec![Oligo {
                id: "p".into(),
                name: "primer".into(),
                ordered: Sequence::dna("GGACGTACGT").unwrap(),
                declared_tail_length: 2,
            }],
        },
        options: Options {
            minimum_annealing_length: 8,
            maximum_extra_unpaired: 0,
            maximum_mismatches: 0,
            minimum_score: 0,
            ..Default::default()
        },
        budget: BudgetLimits::default(),
        conditions: ThermoConditions::default(),
    }
}
#[test]
fn collection_worker_keeps_distinct_footprint_and_hypothetical_whole_oligo_tm() {
    let request = binding();
    let revision = collection_binding::revision(&request).unwrap();
    let response = roundtrip(
        revision,
        AnalysisJob::SearchOligoCollection {
            request: request.clone(),
        },
    );
    let AnalysisResult::OligoCollectionBinding(report) = response.result else {
        panic!("collection report")
    };
    report.validate_for(&request).unwrap();
    assert!(!report.search.hits.is_empty());
    let row = &report.thermodynamics[0];
    let collection_binding::TmValue::Available(a) = &row.annealing else {
        panic!("perfect footprint Tm")
    };
    let collection_binding::TmValue::Available(b) = &row.whole_oligo else {
        panic!("hypothetical whole oligo Tm")
    };
    assert!((a.tm_c - b.tm_c).abs() > 0.01);
    let mut no_work = request.clone();
    no_work.budget.work = 0;
    assert!(collection_binding::assess(&no_work).is_err());
    no_work.budget.work = collagene_oligo_binding::Budget::admission_cost(
        &no_work.source.sequence,
        &no_work.collection,
    )
    .unwrap()
    .work;
    let absent = collection_binding::assess(&no_work).unwrap();
    assert!(absent.search.hits.is_empty());
    assert!(absent.search.oligos.iter().all(|p| p.status
        != collagene_oligo_binding::Status::Complete
        && p.endpoints_examined == 0));
    assert_eq!(
        absent.search.oligos[0].status,
        collagene_oligo_binding::Status::Limited(collagene_oligo_binding::Limit::Work)
    );
    assert!(report.validate_for(&no_work).is_err());
    let mut mismatched = request;
    mismatched.collection.oligos[0].ordered = Sequence::dna("GGACGAACGT").unwrap();
    mismatched.options.maximum_mismatches = 1;
    let uncertain = collection_binding::assess(&mismatched).unwrap();
    assert!(!uncertain.search.hits.is_empty());
    assert!(uncertain.thermodynamics.iter().all(|t| matches!(
        &t.annealing,
        collection_binding::TmValue::MismatchAvailable(_)
    )));
    let collection_binding::TmValue::MismatchAvailable(value) =
        &uncertain.thermodynamics[0].annealing
    else {
        panic!("single internal mismatch estimate")
    };
    // Independently sum canonical AC/CG/GA/AA/AC/CG/GT plus A/T ends,
    // then replace GA/CA and rotated TA/AA using the published table facts.
    assert!((value.estimate.delta_h_kcal_per_mol - -40.0).abs() < 1e-10);
    assert!((value.estimate.delta_s_cal_per_mol_kelvin_at_1m - -110.3).abs() < 1e-10);
    let salt_s = -110.3 + 0.368 * 7.0 * 0.05_f64.ln();
    let expected = -40_000.0 / (salt_s + 1.987 * (500e-9_f64 / 4.0).ln()) - 273.15;
    assert!((value.estimate.tm_c - expected).abs() < 1e-10);
    assert!(!value.estimate.self_complementary);
    uncertain.validate_for(&mismatched).unwrap();
    let mut forged = uncertain;
    forged.thermodynamics[0].annealing = collection_binding::TmValue::Available(a.clone());
    assert!(forged.validate_binding(&forged.context_revision).is_err());
}
