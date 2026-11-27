use collagene_core::{
    Feature, Location, Molecule, Topology,
    analysis::Primer,
    chemistry::{AnnealOptions, EndTreatment, FivePrime, Terminus, ThreePrime},
    mutagenesis::PcrOptions,
    oligo_formats::{self, Delimiter},
    plan_oligos::*,
    verification::{DesignPlan, Operation, PlanInput, Source},
};
fn fixture() -> (DesignPlan, Molecule) {
    let source = Molecule::dna("source", "ACGTACGTACGTACGT", Topology::Linear).unwrap();
    let target = Molecule::dna("independent target", "TTTTCCCC", Topology::Linear).unwrap();
    let policy = PcrOptions::default();
    let plan = DesignPlan {
        name: "Declared orders, not verified PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "s".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("s".into()),
                forward: Primer::new("=forward", "ACGT", "GG").unwrap(),
                reverse: Primer::new("same name", "TTAA", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::PcrWithMismatches {
                template: Source::Step(0),
                forward: Primer::new("renamed", "GACGT", "G").unwrap(),
                reverse: Primer::new("same name", "CCGG", "").unwrap(),
                options: policy,
            },
            Operation::TreatEnds {
                template: Source::Step(1),
                treatment: EndTreatment::Phosphorylate,
            },
        ],
    };
    (plan, target)
}
#[test]
fn exact_order_dedup_preserves_each_name_domain_boundary_and_mismatch_policy() {
    let (plan, target) = fixture();
    let before = serde_json::to_value(&plan).unwrap();
    let report = build(&plan, &target, &Options::default()).unwrap();
    assert!(report.complete);
    assert!(report.unavailable.is_empty());
    assert_eq!(report.oligos.records.len(), 3);
    assert_eq!(report.uses.len(), 4);
    assert_eq!(report.uses[0].oligo_id, report.uses[2].oligo_id);
    assert_ne!(report.uses[1].oligo_id, report.uses[3].oligo_id);
    assert_eq!(report.uses[0].primer_name, "=forward");
    assert_eq!(report.uses[2].primer_name, "renamed");
    assert_eq!(report.uses[0].primer.as_ref().unwrap().tail, "GG");
    assert_eq!(report.uses[2].primer.as_ref().unwrap().tail, "G");
    assert!(report.uses[2].mismatch_options.is_some());
    assert!(report.uses[0].mismatch_options.is_none());
    assert_eq!(
        report.oligos.records[0].raw_ordered_sequence.as_deref(),
        Some("GGACGT")
    );
    assert_eq!(
        report.oligos.records[0].raw_annealing_sequence.as_deref(),
        Some("")
    );
    assert!(
        report
            .oligos
            .records
            .iter()
            .all(|r| r.modification.is_none() && r.quantity.is_none())
    );
    assert!(report.warnings.iter().any(|w| w.contains("partitions")));
    assert_eq!(serde_json::to_value(&plan).unwrap(), before);
}
#[test]
fn spreadsheet_safe_csv_tsv_and_non_dedup_use_exact_declared_dna() {
    let (plan, target) = fixture();
    let csv = export(&plan, &target, &Options::default(), Delimiter::Comma).unwrap();
    assert!(csv.export.text.contains("'=forward"));
    assert!(csv.export.text.contains("GGACGT"));
    csv.validate_binding(
        &export_revision(&plan, &target, &Options::default(), Delimiter::Comma).unwrap(),
    )
    .unwrap();
    let tsv = export(
        &plan,
        &target,
        &Options { deduplicate: false },
        Delimiter::Tab,
    )
    .unwrap();
    assert_eq!(tsv.report.oligos.records.len(), 4);
    assert_ne!(csv.context_revision, tsv.context_revision);
    assert_eq!(tsv.report.uses.len(), 4);
    assert!(tsv.export.text.contains('\t'));
    let value = serde_json::to_vec(&tsv).unwrap();
    let decoded: Exported = serde_json::from_slice(&value).unwrap();
    decoded.validate_binding(&tsv.context_revision).unwrap();
}
#[test]
fn wrong_target_stale_input_and_native_annotation_depth_are_rejected_before_hashing() {
    let (mut plan, target) = fixture();
    let other = Molecule::dna(
        "another independent target",
        target.sequence.as_str(),
        target.topology,
    )
    .unwrap();
    assert!(build(&plan, &other, &Options::default()).is_err());
    plan.inputs[0].revision = "0".repeat(64);
    assert!(build(&plan, &target, &Options::default()).is_err());
    let (mut plan, target) = fixture();
    let mut location = Location::range(0, 1);
    for _ in 0..128 {
        location = Location::Complement(Box::new(location));
    }
    plan.inputs[0].molecule.features.push(Feature {
        kind: "misc_feature".into(),
        location,
        qualifiers: vec![],
        needs_review: false,
    });
    assert!(preflight(&plan, &target, &Options::default()).is_err());
}
#[test]
fn starting_anneal_termini_are_not_order_modifications_and_steps_are_not_synthesized() {
    let (mut plan, target) = fixture();
    let mut original =
        Molecule::dna("declared original ss oligo", "AATT", Topology::Linear).unwrap();
    original
        .metadata
        .insert("genbank.molecule_type".into(), "ss-DNA".into());
    plan.inputs.push(PlanInput {
        id: "oligo".into(),
        revision: original.revision(),
        molecule: original,
    });
    let mut options = AnnealOptions::default();
    options.bottom_chemistry = Terminus {
        five_prime: FivePrime::Phosphate,
        three_prime: ThreePrime::Hydroxyl,
    };
    plan.steps.push(Operation::AnnealOligos {
        top: Source::Step(0),
        bottom: Source::Input("oligo".into()),
        options,
        candidate_index: 0,
    });
    let report = build(&plan, &target, &Options::default()).unwrap();
    assert!(report.complete);
    assert_eq!(report.uses.len(), 5);
    let usage = report.uses.last().unwrap();
    assert_eq!(usage.kind, "anneal-bottom");
    assert_eq!(
        usage.declared_terminus.unwrap().five_prime,
        FivePrime::Phosphate
    );
    assert!(usage.primer.is_none());
    assert!(
        report
            .oligos
            .records
            .iter()
            .all(|r| r.modification.is_none())
    );
    assert!(report.warnings.iter().any(|w| w.contains("derived Step")));
}
#[test]
fn opaque_original_modification_evidence_makes_order_incomplete_without_flattening() {
    use base64::{Engine, engine::general_purpose::STANDARD};
    let (mut plan, target) = fixture();
    let mapping = oligo_formats::OligoImportOptions {
        columns: oligo_formats::ColumnMapping {
            name: 0,
            ordered_sequence: Some(1),
            modification: Some(2),
            ..Default::default()
        },
        ..Default::default()
    };
    let table = oligo_formats::import_oligos(
        "Name,Sequence,Modification\nmodified,ACGT,opaque supplied modification\n",
        &mapping,
    )
    .unwrap();
    let record = &table.records[0];
    assert!(record.primer.is_none());
    let mut original = Molecule::dna("modified original", "ACGT", Topology::Linear).unwrap();
    original.metadata.insert(
        "oligo:original_record_base64".into(),
        record.original_record_base64.clone(),
    );
    original.metadata.insert(
        "oligo:original_record_blake3".into(),
        record.original_record_blake3.clone(),
    );
    original.metadata.insert(
        "oligo:source_mapping".into(),
        serde_json::to_string(&mapping).unwrap(),
    );
    original.metadata.insert(
        "oligo:source_fields".into(),
        serde_json::to_string(&record.original_fields).unwrap(),
    );
    assert!(
        !STANDARD
            .decode(&record.original_record_base64)
            .unwrap()
            .is_empty()
    );
    plan.inputs.push(PlanInput {
        id: "modified".into(),
        revision: original.revision(),
        molecule: original,
    });
    plan.steps.push(Operation::AnnealOligos {
        top: Source::Input("modified".into()),
        bottom: Source::Step(0),
        options: AnnealOptions::default(),
        candidate_index: 0,
    });
    let report = build(&plan, &target, &Options::default()).unwrap();
    assert!(!report.complete);
    assert!(!report.unavailable.is_empty());
    assert_eq!(report.uses.len(), 4);
    assert!(
        report
            .oligos
            .records
            .iter()
            .all(|r| r.modification.is_none())
    );
    assert!(report.unavailable.iter().any(|w| w.contains("opaque")));
}
