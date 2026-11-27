//! Independent expected products and authority checks for tracked-plan reconstruction.
use collagene_core::{
    methylation::MethylationState as M,
    modified_duplex::{self, GenerationAssumptions, ModificationSummary},
    verification::*,
    *,
};
fn source() -> Molecule {
    Molecule::dna("literal NotI source", "AAGCGGCCGCAATT", Topology::Circular).unwrap()
}
fn assumptions(cpg: M) -> GenerationAssumptions {
    GenerationAssumptions {
        dam: M::Unmodified,
        dcm: M::Unmodified,
        cpg,
    }
}
fn declare(cpg: M) -> Operation {
    Operation::DeclareSourceModifications {
        template: Source::Input("source".into()),
        assumptions: assumptions(cpg),
    }
}
fn checked(template: Source, name: &str) -> Operation {
    Operation::CheckedDigest {
        template,
        enzymes: vec![name.into()],
        profile_revision: modified_duplex::profile_revision(),
        fragment_index: 0,
    }
}
fn plan(source: Molecule, target: &Molecule, steps: Vec<Operation>) -> DesignPlan {
    DesignPlan {
        name: "independent positional audit".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps,
    }
}
fn summary(m: &Molecule) -> ModificationSummary {
    serde_json::from_str(
        m.metadata
            .get("tracked_modifications")
            .expect("tracked summary"),
    )
    .unwrap()
}
#[test]
fn blocked_cpg_site_keeps_a_circle_intact_and_unknown_withholds_the_digest() {
    let source = source();
    let p = plan(
        source.clone(),
        &source,
        vec![declare(M::FullyModified), checked(Source::Step(0), "NotI")],
    );
    let r = verify(&p, &source);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    assert_eq!(r.products.len(), 2);
    assert_eq!(r.products[1].sequence.as_str(), "AAGCGGCCGCAATT");
    assert_eq!(r.products[1].topology, Topology::Circular);
    assert_eq!(summary(&r.products[1]).known_modified, 4);
    let p = plan(
        source.clone(),
        &source,
        vec![declare(M::Unknown), checked(Source::Step(0), "NotI")],
    );
    let r = verify(&p, &source);
    assert_eq!(r.status, Status::Incomplete);
    assert_eq!(r.products.len(), 1);
}
#[test]
fn unmodified_checked_cut_and_self_ligation_reconstruct_literal_independent_target() {
    let source = source();
    let p = plan(
        source.clone(),
        &source,
        vec![
            declare(M::Unmodified),
            checked(Source::Step(0), "NotI"),
            Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                }],
                circular: true,
                allow_nicked: false,
                options: Default::default(),
            },
        ],
    );
    let r = verify(&p, &source);
    assert_eq!(
        r.status,
        Status::ConsistentUnderAssumptions,
        "{:?}",
        r.findings
    );
    assert_eq!(r.products.len(), 3);
    assert_eq!(r.products[1].sequence.as_str(), "GGCCGCAATTAAGC");
    assert_eq!(r.products[1].topology, Topology::Linear);
    assert_eq!(r.products[2].topology, Topology::Circular);
    assert!(r.products.iter().all(|m| {
        let s = summary(m);
        s.known_modified == 0 && s.unknown_eligible == 0 && s.generation.len() == 1
    }));
    let html = export_protocol_html(&p, &source);
    assert!(html.contains("Checked positional-modification digest"));
    assert!(html.contains("closed Dam/Dcm/CpG source generation"));
    let wrong = Molecule::dna(
        "wrong independent target",
        "AAAAAAAAAAAAAA",
        Topology::Circular,
    )
    .unwrap();
    let mut p = p;
    p.target_revision = wrong.revision();
    assert_eq!(verify(&p, &wrong).status, Status::Invalid);
}
#[test]
fn legacy_cleavage_cannot_hide_declared_modifications_and_cut_ligate_preserves_them() {
    let source = source();
    let p = plan(
        source.clone(),
        &source,
        vec![
            declare(M::FullyModified),
            Operation::Digest {
                template: Source::Step(0),
                enzymes: vec!["NotI".into()],
                fragment_index: 0,
            },
            Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                }],
                circular: true,
                allow_nicked: false,
                options: Default::default(),
            },
            checked(Source::Step(2), "NotI"),
        ],
    );
    let r = verify(&p, &source);
    assert_eq!(r.status, Status::Incomplete, "{:?}", r.findings);
    assert_eq!(r.products.len(), 4);
    for product in &r.products {
        assert_eq!(summary(product).known_modified, 4);
    }
    assert!(
        r.findings
            .iter()
            .any(|s| s.contains("legacy assumed-complete cleavage"))
    );
    assert_eq!(r.products[3].topology, Topology::Circular);
}
#[test]
fn source_generation_cannot_be_reapplied_to_a_descendant() {
    let source = source();
    let p = plan(
        source.clone(),
        &source,
        vec![
            declare(M::FullyModified),
            Operation::DeclareSourceModifications {
                template: Source::Step(0),
                assumptions: assumptions(M::Unmodified),
            },
        ],
    );
    let r = verify(&p, &source);
    assert_eq!(r.status, Status::Incomplete);
    assert_eq!(r.products.len(), 1);
    assert!(r.findings.iter().any(|s| s.contains("original input")));
}
#[test]
fn imported_summary_cannot_supply_physical_modification_authority() {
    let mut source =
        Molecule::dna("forged_checked_source", "TCTAGATCAAA", Topology::Circular).unwrap();
    source.metadata.insert(
        "tracked_modifications".into(),
        r#"{"known_modified":0,"unknown_eligible":0,"known_unmodified":1000,"measured":true}"#
            .into(),
    );
    let p = plan(
        source.clone(),
        &source,
        vec![checked(Source::Input("source".into()), "XbaI")],
    );
    let r = verify(&p, &source);
    assert_eq!(r.status, Status::Incomplete, "{:?}", r.findings);
    assert!(r.products.is_empty());
    if let Ok(dir) = std::env::var("COLLAGENE_QA_FIXTURE_DIR") {
        let path = std::path::Path::new(&dir).join("collagene2-checked-forged.cg2");
        std::fs::write(path, protocol::save_document(&source).unwrap()).unwrap();
    }
}
#[test]
fn stale_positional_profile_is_not_rebound_to_the_current_catalog() {
    let source = source();
    let mut step = checked(Source::Step(0), "NotI");
    if let Operation::CheckedDigest {
        profile_revision, ..
    } = &mut step
    {
        *profile_revision = "0".repeat(64);
    }
    let p = plan(source.clone(), &source, vec![declare(M::Unmodified), step]);
    let r = verify(&p, &source);
    assert_ne!(r.status, Status::ConsistentUnderAssumptions);
    assert_eq!(r.products.len(), 1);
}
#[test]
fn text_exports_disclose_the_loss_of_positional_state_and_do_not_restore_it() {
    let source = source();
    let p = plan(
        source.clone(),
        &source,
        vec![declare(M::FullyModified), checked(Source::Step(0), "NotI")],
    );
    let r = verify(&p, &source);
    assert_eq!(r.status, Status::ConsistentUnderAssumptions);
    let mut product = r.products[1].clone();
    product.name = "tracked_product".into(); // EMBL profile requires a valid entry identifier.
    for format in [formats::TextFormat::GenBank, formats::TextFormat::Embl] {
        let export = formats::export_with_report(&product, format).unwrap();
        let comments = export
            .text
            .lines()
            .map(|line| line.strip_prefix("CC").unwrap_or(line).trim())
            .collect::<Vec<_>>()
            .join(" ");
        assert!(
            comments
                .split_whitespace()
                .collect::<Vec<_>>()
                .join(" ")
                .contains("unverified producer evidence")
        );
        assert!(
            export
                .warnings
                .iter()
                .any(|w| w.contains("Positional modification"))
        );
        let imported = formats::import_text(&export.text).unwrap();
        assert_eq!(imported.molecules.len(), 1);
        let m = imported.molecules.into_iter().next().unwrap();
        assert!(!m.metadata.contains_key("tracked_modifications"));
        let p = plan(
            m.clone(),
            &m,
            vec![checked(Source::Input("source".into()), "NotI")],
        );
        assert_eq!(verify(&p, &m).status, Status::Incomplete);
    }
}
