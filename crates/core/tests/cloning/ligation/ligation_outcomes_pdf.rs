//! Independently extracted PDF text and raw embedded streams, no external tools.
use collagene_core::{
    chemistry::EndTreatment,
    cloning::ligation::outcomes::Options,
    output::pdf,
    planning::plan::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    sequence::analysis::Primer,
    verification,
};
use collagene_model::{Molecule, Topology};
#[path = "../../common/pdf.rs"]
mod pdf_readback;
use pdf_readback::{check_xref, extracted};

fn fixture(options: Options, nicked: bool) -> (DesignPlan, Molecule) {
    let vector = Molecule::dna("Literal vector", "AGTCCG", Topology::Linear).unwrap();
    let insert = Molecule::dna("Literal insert", "GACTGA", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent intended product",
        "AGTCCGGACTGA",
        Topology::Circular,
    )
    .unwrap();
    let input = |id: &str, m: Molecule| PlanInput {
        id: id.into(),
        revision: m.revision(),
        molecule: m,
    };
    let pcr = |id: &str, f: &str, r: &str| Operation::Pcr {
        template: Source::Input(id.into()),
        forward: Primer::new("F", f, "").unwrap(),
        reverse: Primer::new("R", r, "").unwrap(),
        template_strand: collagene_core::chemistry::Strand::Top,
    };
    let treatment = if nicked {
        EndTreatment::Dephosphorylate
    } else {
        EndTreatment::Phosphorylate
    };
    let plan = DesignPlan {
        name: "Intended path with actual competing orientation".into(),
        target_revision: target.revision(),
        inputs: vec![input("vector", vector), input("insert", insert)],
        steps: vec![
            pcr("vector", "AGT", "CGG"),
            Operation::TreatEnds {
                template: Source::Step(0),
                // At least one covalent strand must seal at each blunt junction.
                // An all-OH mixture cannot be rescued by allow_nicked.
                treatment: EndTreatment::Phosphorylate,
            },
            pcr("insert", "GAC", "TCA"),
            Operation::TreatEnds {
                template: Source::Step(2),
                treatment,
            },
            Operation::Ligate {
                fragments: vec![
                    OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    },
                    OrientedSource {
                        source: Source::Step(3),
                        reverse: false,
                    },
                ],
                circular: true,
                allow_nicked: nicked,
                options,
            },
        ],
    };
    (plan, target)
}
fn streams(bytes: &[u8]) -> Vec<&[u8]> {
    let mut result = vec![];
    let mut start = 0;
    while let Some(relative) = bytes[start..].windows(8).position(|s| s == b"\nstream\n") {
        let content = start + relative + 8;
        let dictionary = String::from_utf8_lossy(&bytes[start..content]);
        let length = dictionary
            .rsplit("/Length ")
            .next()
            .and_then(|s| s.split_whitespace().next())
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap();
        assert!(content + length <= bytes.len());
        result.push(&bytes[content..content + length]);
        start = content + length;
    }
    result
}
#[test]
fn pdf_prints_both_actual_sequences_and_embeds_full_bound_evidence_and_fasta() {
    let (plan, target) = fixture(Options::default(), false);
    let fresh = verification::engine::verify(&plan, &target);
    assert_eq!(
        fresh.status,
        verification::evidence::Status::ConsistentUnderAssumptions,
        "{:?}",
        fresh.findings
    );
    assert_eq!(fresh.products.len(), 5);
    let evidence = fresh
        .physical_processing_reports
        .iter()
        .find_map(|p| match &p.evidence {
            verification::evidence::PhysicalProcessingEvidence::LigationOutcomes(e) => Some(e),
            _ => None,
        })
        .expect("fresh alternative evidence");
    assert_eq!(evidence.report.intended.sequence.as_str(), "AGTCCGGACTGA");
    assert_eq!(evidence.report.alternatives.len(), 1);
    assert_eq!(
        evidence.report.alternatives[0].sequence.as_str(),
        "AGTCCGTCAGTC"
    );
    let output = pdf::export_protocol_pdf(&plan, &target).unwrap();
    // PDF line wrapping is presentation; preserve word order while checking
    // the scientific statements independently of line boundaries.
    let text = extracted(&output.bytes)
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ");
    for required in [
        "Intended path:",
        "Alternative 1:",
        "AGTCCGGACTGA",
        "AGTCCGTCAGTC",
        "Circular",
        "No probabilities are assigned",
        "one copy",
        "Retained nicks",
        "Actual modification state",
        "enumeration complete true",
    ] {
        assert!(text.contains(required), "missing {required}: {text}")
    }
    let objects = streams(&output.bytes);
    let embedded = objects
        .iter()
        .find_map(|s| {
            serde_json::from_slice::<verification::evidence::LigationOutcomeEvidence>(s).ok()
        })
        .expect("full embedded JSON evidence");
    assert_eq!(
        serde_json::to_value(&embedded).unwrap(),
        serde_json::to_value(evidence).unwrap()
    );
    let fasta = objects
        .iter()
        .find(|s| s.starts_with(b">step_5_"))
        .expect("embedded strand FASTA");
    let fasta = std::str::from_utf8(fasta).unwrap();
    assert!(fasta.contains(">step_5_intended_0_top"));
    assert!(fasta.contains(">step_5_alternative_1_bottom"));
    assert!(fasta.contains("AGTCCGGACTGA"));
    assert!(fasta.contains("AGTCCGTCAGTC"));
    let syntax = String::from_utf8_lossy(&output.bytes);
    assert!(syntax.contains("/EmbeddedFiles"));
    assert!(syntax.contains("ligation-step-0005.json"));
    assert!(syntax.contains("ligation-step-0005.fasta"));
    assert!(!syntax.contains("/OpenAction"));
    assert!(!syntax.contains("/JavaScript"));
    check_xref(&output.bytes);
}
#[test]
fn bounded_enumeration_and_real_nicks_remain_explicit_in_pdf() {
    let mut options = Options::default();
    options.max_products = 1;
    let (plan, target) = fixture(options, true);
    let fresh = verification::engine::verify(&plan, &target);
    assert_eq!(
        fresh.status,
        verification::evidence::Status::ConsistentUnderAssumptions,
        "{:?}",
        fresh.findings
    );
    assert_eq!(fresh.products.len(), 5);
    let evidence = fresh
        .physical_processing_reports
        .iter()
        .find_map(|p| match &p.evidence {
            verification::evidence::PhysicalProcessingEvidence::LigationOutcomes(e) => Some(e),
            _ => None,
        })
        .unwrap();
    assert!(!evidence.report.intended.duplex.nicks.is_empty());
    assert!(!evidence.report.enumeration_complete);
    let output = pdf::export_protocol_pdf(&plan, &target).unwrap();
    // PDF line wrapping is presentation; preserve word order while checking
    // the scientific statements independently of line boundaries.
    let text = extracted(&output.bytes)
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ");
    assert!(text.contains("enumeration complete false"));
    assert!(text.contains("Products"));
    assert!(text.contains("Retained nicks"));
    assert!(text.contains("unresolved"));
    assert!(text.contains("not absent alternatives"));
    let embedded = streams(&output.bytes)
        .into_iter()
        .find_map(|s| {
            serde_json::from_slice::<verification::evidence::LigationOutcomeEvidence>(s).ok()
        })
        .unwrap();
    assert_eq!(
        embedded.report.intended.duplex.nicks,
        evidence.report.intended.duplex.nicks
    );
    assert_eq!(
        embedded.report.withheld_products,
        evidence.report.withheld_products
    );
}
