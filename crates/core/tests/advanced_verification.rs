//! Parent integration audit: external targets, hand-derived products and historical envelopes.
use collagene_core::{
    analysis::Primer,
    chemistry::{Duplex, DuplexSummary},
    golden_gate::GoldenGateOptions,
    mutagenesis::PcrOptions,
    primers::BindingOptions,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    transfer::{self, Candidate},
    verification::{self, DesignPlan, Operation, PlanInput, Source, Status},
    *,
};
fn dna(name: &str, s: &str, t: Topology) -> Molecule {
    Molecule::dna(name, s, t).unwrap()
}
fn feature(label: &str, location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![
            ("label".into(), Some(label.into())),
            ("translation".into(), Some("SOURCE_ASSERTION".into())),
            ("codon_start".into(), Some("1".into())),
        ],
        needs_review: false,
    }
}
fn named<'a>(m: &'a Molecule, label: &str) -> &'a Feature {
    m.features
        .iter()
        .find(|f| f.label() == label)
        .unwrap_or_else(|| panic!("missing {label}: {:?}", m.features))
}
fn source() -> Molecule {
    let mut source = dna(
        "original",
        "AGTCCGATGCTACGTAGGACTTCAGACCTGTA",
        Topology::Linear,
    );
    source
        .features
        .push(feature("changed forward", Location::range(0, 8)));
    source
        .features
        .push(feature("intact", Location::range(8, 24)));
    source.features.push(feature(
        "changed reverse",
        Location::Complement(Box::new(Location::range(24, 32))),
    ));
    source
}
fn pcr_options() -> PcrOptions {
    PcrOptions {
        binding: BindingOptions {
            max_mismatches: 1,
            min_three_prime_exact: 2,
            ..Default::default()
        },
        ..Default::default()
    }
}
fn primers() -> (Primer, Primer) {
    (
        Primer::new("F", "AGACCGAT", "AA").unwrap(),
        Primer::new("R", "TCCAGGTC", "CC").unwrap(),
    )
}
fn plan(target: &Molecule, inputs: Vec<(&str, Molecule)>, steps: Vec<Operation>) -> DesignPlan {
    DesignPlan {
        name: "Independent integration audit".into(),
        target_revision: target.revision(),
        inputs: inputs
            .into_iter()
            .map(|(id, molecule)| PlanInput {
                id: id.into(),
                revision: molecule.revision(),
                molecule,
            })
            .collect(),
        steps,
    }
}
fn mismatch_step(template: Source, options: PcrOptions) -> Operation {
    let (f, r) = primers();
    Operation::PcrWithMismatches {
        template,
        forward: f,
        reverse: r,
        options,
    }
}
fn mutant_target() -> Molecule {
    dna(
        "externally specified target",
        "AAAGACCGATGCTACGTAGGACTTCAGACCTGGAGG",
        Topology::Linear,
    )
}
#[test]
fn typed_mismatch_reconstruction_changes_both_gene_strands_and_preserves_unchanged_assertions() {
    let source = source();
    let target = mutant_target();
    let plan = plan(
        &target,
        vec![("template", source.clone())],
        vec![mismatch_step(
            Source::Input("template".into()),
            pcr_options(),
        )],
    );
    let report = verification::verify(&plan, &target);
    assert_eq!(
        report.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        report.findings
    );
    assert_eq!(report.products.len(), 1);
    let product = &report.products[0];
    assert_eq!(product.sequence, target.sequence);
    assert_eq!(
        named(product, "changed forward").location,
        Location::range(2, 10)
    );
    assert_eq!(
        named(product, "changed reverse")
            .location
            .extract(&product.sequence)
            .unwrap(),
        "TCCAGGTC"
    );
    for label in ["changed forward", "changed reverse"] {
        let f = named(product, label);
        assert!(f.needs_review);
        assert!(
            !f.qualifiers
                .iter()
                .any(|(k, _)| k == "translation" || k == "codon_start")
        );
        assert!(f.qualifiers.iter().any(
            |(k, v)| k == "note" && v.as_ref().is_some_and(|v| v.contains("SOURCE_ASSERTION"))
        ));
    }
    assert!(!named(product, "intact").needs_review);
    assert!(
        named(product, "intact")
            .qualifiers
            .contains(&("translation".into(), Some("SOURCE_ASSERTION".into())))
    );
    let physical: DuplexSummary = serde_json::from_str(&product.metadata["duplex_model"]).unwrap();
    assert_eq!(physical.top_length, target.sequence.len());
    assert_eq!(physical.bottom_length, target.sequence.len());
    assert_eq!(physical.bottom_start, 0);
    assert!(physical.nicks.is_empty());
    assert_eq!(
        named(&source, "changed forward").location,
        Location::range(0, 8)
    );
    assert!(!named(&source, "changed forward").needs_review);
}
#[test]
fn wrong_independent_target_and_strict_three_prime_policy_do_not_inherit_success() {
    let source = source();
    let wrong = dna(
        "wrong target",
        "AAAGTCCGATGCTACGTAGGACTTCAGACCTGTAGG",
        Topology::Linear,
    );
    let p = plan(
        &wrong,
        vec![("template", source.clone())],
        vec![mismatch_step(
            Source::Input("template".into()),
            pcr_options(),
        )],
    );
    let r = verification::verify(&p, &wrong);
    assert_eq!(r.status, Status::Invalid);
    assert_eq!(r.products.len(), 1);
    assert!(r.findings.iter().any(|s| s.contains("does not match")));
    let target = mutant_target();
    let strict = PcrOptions {
        binding: BindingOptions {
            min_three_prime_exact: 8,
            ..pcr_options().binding
        },
        ..pcr_options()
    };
    let p = plan(
        &target,
        vec![("template", source)],
        vec![mismatch_step(Source::Input("template".into()), strict)],
    );
    let r = verification::verify(&p, &target);
    assert_eq!(r.status, Status::Invalid, "{:?}", r.findings);
    assert!(r.products.is_empty());
}
#[test]
fn truncated_or_overlapping_mismatch_models_cannot_become_verified_duplexes() {
    let source = source();
    let target = mutant_target();
    let options = PcrOptions {
        binding: BindingOptions {
            max_work: 8,
            ..pcr_options().binding
        },
        ..pcr_options()
    };
    let p = plan(
        &target,
        vec![("template", source)],
        vec![mismatch_step(Source::Input("template".into()), options)],
    );
    let r = verification::verify(&p, &target);
    assert_eq!(r.status, Status::Incomplete, "{:?}", r.findings);
    let source = dna("overlap", "ACGTAA", Topology::Linear);
    let target = dna("overlap sequence preview", "ATGTAA", Topology::Linear);
    let p = plan(
        &target,
        vec![("template", source)],
        vec![Operation::PcrWithMismatches {
            template: Source::Input("template".into()),
            forward: Primer::new("F", "ATGT", "").unwrap(),
            reverse: Primer::new("R", "TTAC", "").unwrap(),
            options: PcrOptions {
                binding: BindingOptions {
                    max_mismatches: 1,
                    min_three_prime_exact: 1,
                    ..Default::default()
                },
                ..Default::default()
            },
        }],
    );
    let r = verification::verify(&p, &target);
    assert_eq!(r.status, Status::Incomplete, "{:?}", r.findings);
    assert!(r.products.is_empty());
}
#[test]
fn imported_physical_metadata_never_grants_mismatch_pcr_substrate_authority() {
    let mut source = source();
    let target = mutant_target();
    source.metadata.insert("duplex_model".into(), "{}".into());
    let p = plan(
        &target,
        vec![("template", source)],
        vec![mismatch_step(
            Source::Input("template".into()),
            pcr_options(),
        )],
    );
    let r = verification::verify(&p, &target);
    assert_eq!(r.status, Status::Incomplete);
    assert!(r.findings.iter().any(|f| f.contains("physical evidence")));
}
#[test]
fn mismatch_pcr_accepts_an_actual_intact_preceding_pcr_intermediate() {
    let source = source();
    let target = mutant_target();
    let p = plan(
        &target,
        vec![("template", source)],
        vec![
            Operation::Pcr {
                template: Source::Input("template".into()),
                forward: Primer::new("exact F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("exact R", "TACAGGTC", "").unwrap(),
            },
            mismatch_step(Source::Step(0), pcr_options()),
        ],
    );
    let r = verification::verify(&p, &target);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    assert_eq!(r.products.len(), 2);
}
fn donor(left: &str, payload: &str, right: &str) -> Molecule {
    dna(
        "independent donor",
        &format!("GGTCTCA{left}{payload}{right}AGAGACC"),
        Topology::Linear,
    )
}
fn gate_inputs() -> Vec<Molecule> {
    vec![donor("AATG", "CCCC", "GCTT"), donor("GCTT", "TTTT", "AATG")]
}
fn gate_target() -> Molecule {
    dna("gate target", "AATGCCCCGCTTTTTT", Topology::Circular)
}
fn gate_options() -> GoldenGateOptions {
    GoldenGateOptions {
        max_fragments: 4,
        ..Default::default()
    }
}
fn gate_plan(
    target: &Molecule,
    inputs: Vec<Molecule>,
    options: GoldenGateOptions,
    index: usize,
) -> DesignPlan {
    plan(
        target,
        inputs
            .into_iter()
            .enumerate()
            .map(|(i, m)| (if i == 0 { "a" } else { "b" }, m))
            .collect(),
        vec![Operation::GoldenGate {
            templates: vec![Source::Input("a".into()), Source::Input("b".into())],
            options,
            candidate_index: index,
        }],
    )
}
#[test]
fn golden_gate_selected_candidate_is_reconstructed_with_annotations_and_external_target() {
    let mut inputs = gate_inputs();
    inputs[0]
        .features
        .push(feature("payload", Location::range(11, 15)));
    let target = gate_target();
    let p = gate_plan(&target, inputs.clone(), gate_options(), 0);
    let r = verification::verify(&p, &target);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    assert_eq!(r.products.len(), 1);
    let product = &r.products[0];
    assert_eq!(product.sequence, target.sequence);
    assert_eq!(product.topology, Topology::Circular);
    assert_eq!(named(product, "payload").location, Location::range(4, 8));
    assert_eq!(
        named(product, "payload")
            .location
            .extract(&product.sequence)
            .unwrap(),
        "CCCC"
    );
    assert!(!named(product, "payload").needs_review);
    assert!(!r.exhaustive_alternative_search);
    assert!(
        r.findings
            .iter()
            .any(|s| s.contains("one-copy") && s.contains("uniqueness"))
    );
    let wrong = dna("wrong gate target", "AATGCCCTGCTTTTTT", Topology::Circular);
    let p = gate_plan(&wrong, inputs, gate_options(), 0);
    let r = verification::verify(&p, &wrong);
    assert_eq!(r.status, Status::Invalid);
    assert_eq!(r.products.len(), 1);
}
#[test]
fn golden_gate_recutting_and_unavailable_candidate_are_conclusive_failures() {
    let inputs = vec![
        donor("AATG", "CCCCG", "GTCT"),
        donor("GTCT", "CAAAA", "AATG"),
    ];
    let target = dna(
        "unstable cycling target",
        "AATGCCCCGGTCTCAAAA",
        Topology::Circular,
    );
    let p = gate_plan(&target, inputs, gate_options(), 0);
    let r = verification::verify(&p, &target);
    assert_eq!(r.status, Status::Invalid, "{:?}", r.findings);
    assert!(r.findings.iter().any(|s| s.contains("retains recognition")));
    let target = gate_target();
    let p = gate_plan(&target, gate_inputs(), gate_options(), usize::MAX);
    let r = verification::verify(&p, &target);
    assert_eq!(r.status, Status::Invalid);
    assert!(r.products.is_empty());
}
#[test]
fn golden_gate_capped_search_retains_known_product_but_not_verified_uniqueness() {
    let target = gate_target();
    let options = GoldenGateOptions {
        max_products: 1,
        ..gate_options()
    };
    let p = gate_plan(&target, gate_inputs(), options, 0);
    let r = verification::verify(&p, &target);
    assert_eq!(r.status, Status::Incomplete, "{:?}", r.findings);
    assert_eq!(r.products.len(), 1);
    assert_eq!(r.products[0].sequence, target.sequence);
    assert!(
        r.findings
            .iter()
            .any(|s| s.contains("enumeration incomplete"))
    );
    let p = gate_plan(
        &target,
        gate_inputs(),
        GoldenGateOptions {
            max_states: 1,
            ..gate_options()
        },
        0,
    );
    let r = verification::verify(&p, &target);
    assert_eq!(r.status, Status::Incomplete);
    assert!(r.products.is_empty());
}
#[test]
fn golden_gate_accepts_reconstructed_staggered_digest_steps_without_imported_metadata_authority() {
    let target = gate_target();
    let inputs = gate_inputs();
    let p = plan(
        &target,
        vec![("a", inputs[0].clone()), ("b", inputs[1].clone())],
        vec![
            Operation::Digest {
                template: Source::Input("a".into()),
                enzymes: vec!["BsaI".into()],
                fragment_index: 1,
            },
            Operation::Digest {
                template: Source::Input("b".into()),
                enzymes: vec!["BsaI".into()],
                fragment_index: 1,
            },
            Operation::GoldenGate {
                templates: vec![Source::Step(0), Source::Step(1)],
                options: gate_options(),
                candidate_index: 0,
            },
        ],
    );
    let r = verification::verify(&p, &target);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    assert_eq!(r.products.len(), 3);
    assert_eq!(r.products[2].sequence, target.sequence);
    let mut imported = gate_inputs();
    imported[0]
        .metadata
        .insert("duplex_model".into(), "{}".into());
    let p = gate_plan(&target, imported, gate_options(), 0);
    let r = verification::verify(&p, &target);
    assert_eq!(r.status, Status::Incomplete);
}
fn worker(source: &Molecule, job: AnalysisJob) -> protocol::Response {
    let request = Request {
        id: 912,
        revision: source.revision(),
        job,
    };
    let bytes = protocol::encode_request(&request).unwrap();
    assert_eq!(&bytes[4..6], &21u16.to_le_bytes());
    let response = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    assert_eq!(response.id, 912);
    assert_eq!(response.revision, source.revision());
    response
}
#[test]
fn new_worker_jobs_preserve_results_annotation_review_and_request_identity() {
    let source = source();
    let (f, r) = primers();
    let response = worker(
        &source,
        AnalysisJob::MismatchPcrPreview {
            molecule: source.clone(),
            forward: f,
            reverse: r,
            options: pcr_options(),
        },
    );
    let AnalysisResult::MismatchPcrPreview(preview) = response.result else {
        panic!("wrong worker result")
    };
    assert_eq!(preview.report.products.len(), 1);
    assert_eq!(preview.documents.len(), 1);
    assert_eq!(preview.documents[0].sequence, mutant_target().sequence);
    assert!(named(&preview.documents[0], "changed reverse").needs_review);
    assert!(!preview.documents[0].metadata.contains_key("duplex_model"));
    let inputs = gate_inputs();
    let response = worker(
        &inputs[0],
        AnalysisJob::GoldenGate {
            molecule: inputs[0].clone(),
            templates: inputs.clone(),
            options: gate_options(),
        },
    );
    let AnalysisResult::GoldenGate(report) = response.result else {
        panic!("wrong gate result")
    };
    assert_eq!(report.stable_physical_product_count, 1);
    assert!(!report.exhaustive_mixture_search);
}
#[test]
fn old_envelopes_cannot_smuggle_new_steps_or_new_worker_results() {
    let source = source();
    let target = mutant_target();
    let p = plan(
        &target,
        vec![("template", source.clone())],
        vec![mismatch_step(
            Source::Input("template".into()),
            pcr_options(),
        )],
    );
    for p in [
        p,
        gate_plan(&gate_target(), gate_inputs(), gate_options(), 0),
    ] {
        let candidate = Candidate::new(p.clone());
        let current = transfer::save_candidate(&candidate).unwrap();
        assert_eq!(&current[4..6], &6u16.to_le_bytes());
        for old in [1u16, 2] {
            let mut bytes = current.clone();
            bytes[4..6].copy_from_slice(&old.to_le_bytes());
            assert!(matches!(
                transfer::load_candidate(&bytes),
                Err(Error::Unsupported(_))
            ));
        }
        for (label, job) in [
            (
                "verify",
                AnalysisJob::Verify {
                    target: target.clone(),
                    plan: p.clone(),
                },
            ),
            (
                "PDF",
                AnalysisJob::ExportProtocolPdf {
                    target: target.clone(),
                    plan: p,
                },
            ),
        ] {
            let request = Request {
                id: 81,
                revision: target.revision(),
                job,
            };
            let mut bytes = protocol::encode_request(&request).unwrap();
            bytes[4..6].copy_from_slice(&4u16.to_le_bytes());
            assert!(
                matches!(
                    protocol::process_request(&bytes),
                    Err(Error::Unsupported(_))
                ),
                "{label} must also enforce the schema of its nested plan operations"
            );
        }
    }
    let (f, r) = primers();
    let request = Request {
        id: 82,
        revision: source.revision(),
        job: AnalysisJob::MismatchPcr {
            molecule: source,
            forward: f,
            reverse: r,
            options: pcr_options(),
        },
    };
    let mut bytes = protocol::encode_request(&request).unwrap();
    let mut reply = protocol::process_request(&bytes).unwrap();
    bytes[4..6].copy_from_slice(&4u16.to_le_bytes());
    assert!(matches!(
        protocol::process_request(&bytes),
        Err(Error::Unsupported(_))
    ));
    reply[4..6].copy_from_slice(&4u16.to_le_bytes());
    assert!(matches!(
        protocol::decode_response(&reply),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn a_known_selected_product_mismatching_target_is_invalid_even_when_other_paths_are_unsearched() {
    let wrong = dna(
        "wrong target under capped search",
        "AATGCCCTGCTTTTTT",
        Topology::Circular,
    );
    let p = gate_plan(
        &wrong,
        gate_inputs(),
        GoldenGateOptions {
            max_products: 1,
            ..gate_options()
        },
        0,
    );
    let r = verification::verify(&p, &wrong);
    assert_eq!(r.products.len(), 1);
    assert_eq!(r.products[0].sequence.as_str(), "AATGCCCCGCTTTTTT");
    assert_ne!(r.products[0].sequence, wrong.sequence);
    assert_eq!(
        r.status,
        Status::Invalid,
        "A capped search cannot make this explicitly selected, fully reconstructed wrong candidate unknown: {:?}",
        r.findings
    );
    assert!(r.findings.iter().any(|s| s.contains("does not match")));
}

#[test]
fn physical_template_mismatch_and_partial_step_failure_remain_explicit_not_flattened() {
    let source = source();
    let (f, r) = primers();
    let wrong = dna(
        "different paired substrate",
        "TGACCGATGCTACGTAGGACTTCAGACCTGTA",
        Topology::Linear,
    );
    let substrate = Duplex::from_molecule(&wrong).unwrap();
    assert!(
        collagene_core::mutagenesis::simulate_pcr_reconstructed(
            &source,
            &substrate,
            &f,
            &r,
            &pcr_options()
        )
        .is_err()
    );
    let substrate = Duplex::from_molecule(&source).unwrap();
    assert!(
        collagene_core::mutagenesis::simulate_pcr_reconstructed(
            &source,
            &substrate,
            &f,
            &r,
            &pcr_options()
        )
        .unwrap()
        .enumeration_complete
    );
    let target = mutant_target();
    let p = plan(
        &target,
        vec![("donor", donor("AATG", "CCCC", "GCTT"))],
        vec![
            Operation::Digest {
                template: Source::Input("donor".into()),
                enzymes: vec!["BsaI".into()],
                fragment_index: 1,
            },
            mismatch_step(Source::Step(0), pcr_options()),
        ],
    );
    let report = verification::verify(&p, &target);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert_eq!(report.products.len(), 1);
    assert_eq!(report.products[0].sequence.as_str(), "AATGCCCC");
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.contains("staggered/nicked"))
    );
    assert!(
        !report.findings.iter().any(|f| f.contains("does not match")),
        "A partial intermediate is not the final product"
    );
}
