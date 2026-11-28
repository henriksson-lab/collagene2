//! Documentary choices remain explicit constraints independent of sequence reconstruction.
use collagene_core::{
    chemistry::{Duplex, tracked_duplex::TrackedDuplex},
    cloning::{
        assembly::recipe::{self as assembly, FragmentSpec, Preparation},
        ligation::recipe as ligation,
        restriction::{
            catalog,
            conditions::{
                ChosenFormulation, EvidenceStatus, ReactionConditions, source::SourceReport,
            },
        },
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification::{engine, evidence::Status, output::export_protocol_html},
};
use collagene_model::{Molecule, Topology};
#[path = "../../common/pdf.rs"]
mod pdf_readback;
fn choices(enzyme: &str, formulation: Option<&str>) -> ReactionConditions {
    ReactionConditions {
        schema_version: 1,
        formulations: vec![ChosenFormulation {
            enzyme_id: enzyme.into(),
            formulation_id: formulation.map(str::to_owned),
        }],
        buffer: "rCutSmart".into(),
        temperature_milli_c: 37000,
        incubation_minutes: 60,
        minimum_activity_percent: 75,
        activity_context: Some("NEB unit-assay relative activity".into()),
        temperature_context: Some("vendor recommended incubation".into()),
    }
}
fn fixture(conditions: Option<ReactionConditions>) -> (DesignPlan, Molecule) {
    let source = Molecule::dna("original", "AAGAATTCTT", Topology::Linear).unwrap();
    let target = Molecule::dna("independent reference", "AATTCTT", Topology::Linear).unwrap();
    let selection =
        catalog::select(&catalog::builtin_catalog().unwrap(), &["ecori".into()]).unwrap();
    (
        DesignPlan {
            name: "documentary digest".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "source".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps: vec![Operation::CatalogDigest {
                template: Source::Input("source".into()),
                selection,
                conditions,
                fragment_index: 1,
            }],
        },
        target,
    )
}
#[test]
fn declared_unknown_and_unmet_conditions_keep_products_and_independent_comparison() {
    for (formulation, status, outcome) in [
        (
            Some("neb-r0101"),
            EvidenceStatus::BelowThreshold,
            Status::Incomplete,
        ),
        (
            Some("neb-r3101"),
            EvidenceStatus::Known,
            Status::ConsistentUnderAssumptions,
        ),
        (None, EvidenceStatus::Unknown, Status::Incomplete),
    ] {
        let (plan, target) = fixture(Some(choices("ecori", formulation)));
        let report = engine::verify(&plan, &target);
        assert_eq!(report.status, outcome, "{:?}", report.findings);
        assert_eq!(report.products.len(), 1);
        assert_eq!(report.products[0].sequence.as_str(), "AATTCTT");
        assert!(
            report
                .findings
                .iter()
                .any(|f| f.contains("matches the separately bound target"))
        );
        let evidence = report.step_reports[0].reaction_conditions.as_ref().unwrap();
        assert_eq!(evidence.source_revision, plan.inputs[0].revision);
        assert_eq!(evidence.report.documentary_status, status);
        report.step_reports[0].validate_documentary().unwrap();
        let html = export_protocol_html(&plan, &target);
        assert!(html.contains("Documentary reaction-condition evidence:"));
        assert!(
            html.contains("Star activity and incubation-duration compatibility remain unassessed")
        );
        assert!(html.contains("Documentary source provenance:"));
        assert!(html.contains("incubation_minutes&quot;:60"));
        let pdf = collagene_core::output::pdf::export_protocol_pdf(&plan, &target).unwrap();
        assert_eq!(pdf.status, outcome);
        pdf_readback::check_xref(&pdf.bytes);
        let text = pdf_readback::extracted(&pdf.bytes)
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ");
        assert!(text.contains("Documentary reaction-condition evidence"));
        assert!(text.contains("neb-r") || formulation.is_none());
    }
}
#[test]
fn no_choices_are_unassessed_and_current_plan_requires_explicit_nullable_field() {
    let (plan, target) = fixture(None);
    let report = engine::verify(&plan, &target);
    assert!(report.step_reports.is_empty());
    assert_eq!(report.status, Status::ConsistentUnderAssumptions);
    assert!(
        export_protocol_html(&plan, &target)
            .contains("Documentary reaction conditions unassessed:")
    );
    let encoded = postcard::to_allocvec(&plan).unwrap();
    let decoded: DesignPlan = postcard::from_bytes(&encoded).unwrap();
    assert_eq!(
        serde_json::to_value(&decoded).unwrap(),
        serde_json::to_value(&plan).unwrap()
    );
    let mut json = serde_json::to_value(&plan).unwrap();
    json["steps"][0]["CatalogDigest"]
        .as_object_mut()
        .unwrap()
        .remove("conditions");
    assert!(serde_json::from_value::<DesignPlan>(json).is_err());
}
#[test]
fn source_material_options_and_retained_evidence_are_bound_without_promoting_history() {
    let (plan, _) = fixture(None);
    let source = &plan.inputs[0].molecule;
    let material =
        TrackedDuplex::from_unknown_duplex(&Duplex::from_molecule(source).unwrap()).unwrap();
    let Operation::CatalogDigest { selection, .. } = &plan.steps[0] else {
        unreachable!()
    };
    let c = catalog::builtin_catalog().unwrap();
    let options = choices("ecori", Some("neb-r3101"));
    let evidence = SourceReport::assess(&c, source, &material, selection, &options).unwrap();
    evidence
        .validate_for(&c, source, &material, selection, &options)
        .unwrap();
    let mut changed = options.clone();
    changed.incubation_minutes = 61;
    assert!(
        evidence
            .validate_for(&c, source, &material, selection, &changed)
            .is_err()
    );
    let mut tampered = evidence.clone();
    tampered.report.documentary_status = EvidenceStatus::Unknown;
    assert!(tampered.validate_historical().is_err());
    let mut wrong = source.clone();
    wrong.name = "different source identity".into();
    assert!(
        evidence
            .validate_for(&c, &wrong, &material, selection, &options)
            .is_err()
    );
    let foreign = Molecule::dna("foreign", "AAGAATTCTA", Topology::Linear).unwrap();
    assert!(SourceReport::assess(&c, &foreign, &material, selection, &options).is_err());
}
#[test]
fn pcr_digest_conditions_bind_fresh_amplicon_and_html_escapes_declared_context() {
    let source = Molecule::dna("original", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    let target = Molecule::dna("independent", "ACGACCCCGGGG", Topology::Linear).unwrap();
    let mut options = choices("bsai", None);
    options.activity_context = Some("<declared & unresolved>".into());
    let selection =
        catalog::select(&catalog::builtin_catalog().unwrap(), &["bsai".into()]).unwrap();
    let recipe = ligation::Recipe {
        name: "PCR digest".into(),
        fragments: vec![ligation::Fragment {
            source: FragmentSpec {
                generation: None,
                role: "insert".into(),
                input_id: "s".into(),
                reverse: false,
                preparation: Preparation::PcrCatalogDigest {
                    forward: Primer::new("F", "ACGACCCC", "GGTCTCA").unwrap(),
                    reverse: Primer::new("R", "GTTACCCC", "GGTCTCA").unwrap(),
                    template_strand: collagene_core::chemistry::Strand::Top,
                    selection,
                    conditions: Some(options.clone()),
                    fragment_index: 1,
                    repair: None,
                },
            },
            treatment: None,
        }],
        circular: false,
        allow_nicked: false,
    };
    let input = PlanInput {
        id: "s".into(),
        revision: source.revision(),
        molecule: source,
    };
    let preview = ligation::preview(
        &[input],
        &target,
        &recipe,
        assembly::CleavagePolicy::AssumedComplete,
    )
    .unwrap();
    assert_eq!(
        preview.verification.status,
        Status::Incomplete,
        "{:?}",
        preview.verification.findings
    );
    assert_eq!(preview.verification.products.len(), 3);
    let evidence = preview
        .verification
        .step_reports
        .iter()
        .find_map(|r| r.reaction_conditions.as_ref())
        .unwrap();
    assert_eq!(
        evidence.source_revision,
        preview.verification.products[0].revision()
    );
    assert_ne!(
        evidence.source_revision,
        preview.compiled.plan.inputs[0].revision
    );
    assert_eq!(evidence.report.options, options);
    let html = export_protocol_html(&preview.compiled.plan, &target);
    assert!(html.contains("&lt;declared &amp; unresolved&gt;"));
    assert!(!html.contains("<declared & unresolved>"));
    preview.validate_binding(&preview.context_revision).unwrap();
}

#[test]
fn checked_digest_preserves_documentary_failure_separately_from_positional_cleavage() {
    use collagene_core::chemistry::{
        methylation::MethylationState as M,
        tracked_duplex::{GenerationAssumptions, checked_inspection_options},
    };
    let (mut plan, target) = fixture(Some(choices("ecori", Some("neb-r0101"))));
    let Operation::CatalogDigest {
        selection,
        conditions,
        fragment_index,
        ..
    } = plan.steps.remove(0)
    else {
        unreachable!()
    };
    plan.steps = vec![
        Operation::DeclareSourceModifications {
            template: Source::Input("source".into()),
            assumptions: GenerationAssumptions {
                dam: M::Unmodified,
                dcm: M::Unmodified,
                cpg: M::Unmodified,
            },
        },
        Operation::CheckedDigest {
            template: Source::Step(0),
            selection,
            conditions,
            fragment_index,
            options: checked_inspection_options(),
        },
    ];
    let report = engine::verify(&plan, &target);
    assert_eq!(report.status, Status::Incomplete, "{:?}", report.findings);
    assert_eq!(report.products.len(), 2);
    assert_eq!(report.products[1].sequence.as_str(), "AATTCTT");
    let evidence = report.step_reports[0].reaction_conditions.as_ref().unwrap();
    assert_eq!(report.step_reports[0].step_index, 1);
    assert_eq!(
        evidence.report.documentary_status,
        EvidenceStatus::BelowThreshold
    );
    assert_eq!(evidence.source_revision, report.products[0].revision());
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.contains("matches the separately bound target"))
    );
}
