use collagene_core::{
    analysis::{self, Enzyme},
    chemistry,
    restriction_inspection::*,
};
use collagene_model::{Feature, Location, Molecule, Sequence, Topology};
fn enzyme(name: &str) -> Enzyme {
    analysis::starter_enzymes()
        .into_iter()
        .find(|e| e.name == name)
        .unwrap()
}
fn dna(s: &str, t: Topology) -> Molecule {
    Molecule::dna("independent", s, t).unwrap()
}
fn run(m: &Molecule, es: &[Enzyme]) -> RestrictionInspection {
    inspect("source", m, es, &InspectionOptions::default()).unwrap()
}
#[test]
fn edge_type_iis_preserves_valid_interior_and_outside_recognition_without_products() {
    let m = dna("AAAGGTCTCTAGTCCCCCAAAAGGTCTC", Topology::Linear);
    let es = [enzyme("BsaI")];
    assert!(analysis::restriction_sites(&m, &es[0]).is_err());
    let r = run(&m, &es);
    assert!(r.scan_complete);
    assert_eq!(r.sites.len(), 2);
    assert!(r.sites[0].top_cut.is_some() && r.sites[0].bottom_cut.is_some());
    let edge = &r.sites[1];
    assert_eq!(edge.recognition_start, 22);
    assert_eq!(edge.recognition_length, 6);
    assert!(
        edge.diagnostics
            .contains(&SiteDiagnostic::CutsOutsideSubstrate)
    );
    assert!(edge.bottom_cut.is_none());
    assert!(!r.products_determined && r.fragments.is_empty() && r.ref_top_lengths.is_none());
    assert!(product("source", &m, &es, &InspectionOptions::default(), 0).is_err());
}
#[test]
fn exact_strand_products_have_distinct_top_lengths_and_real_end_chemistry() {
    let m = dna("TTTGAATTCAAA", Topology::Linear);
    let es = [enzyme("EcoRI")];
    let r = run(&m, &es);
    assert!(r.scan_complete && r.products_determined);
    assert_eq!(r.ref_top_lengths.as_ref().unwrap(), &vec![8, 4]);
    assert_eq!(r.fragments.len(), 2);
    assert_eq!(
        (
            r.fragments[0].duplex.top_length,
            r.fragments[0].duplex.bottom_length
        ),
        (4, 8)
    );
    let a = product("source", &m, &es, &r.options, 0).unwrap();
    let b = product("source", &m, &es, &r.options, 1).unwrap();
    assert_eq!(a.molecule.sequence.as_str(), "TTTG");
    assert_eq!(b.molecule.sequence.as_str(), "AATTCAAA");
    assert_ne!(a.batch_revision, b.batch_revision);
    assert_eq!(a.inspection_revision, r.batch_revision);
    a.validate_bindings_for("source", &m, &es, &r.options, 0)
        .unwrap();
    let (_, right) = a.summary.duplex.ends.unwrap();
    assert_eq!(right.overhang, "AATT");
    assert_eq!(right.chemistry, chemistry::Terminus::CLEAVED);
    assert!(a.molecule.original_genbank.is_none());
}
#[test]
fn every_origin_and_reverse_type_iis_geometry_matches_literal_coordinate_oracle() {
    let template = "AAAGGTCTCTAGTCCCCCAA";
    let es = [enzyme("BsaI")];
    for reverse in [false, true] {
        let s = if reverse {
            Sequence::dna(template)
                .unwrap()
                .reverse_complement()
                .unwrap()
                .as_str()
                .to_owned()
        } else {
            template.into()
        };
        for origin in 0..s.len() {
            let rotated = format!("{}{}", &s[origin..], &s[..origin]);
            let m = dna(&rotated, Topology::Circular);
            let r = run(&m, &es);
            assert!(r.scan_complete && r.products_determined);
            assert_eq!(r.sites.len(), 1);
            let site = &r.sites[0];
            let original_start = if reverse { s.len() - 3 - 6 } else { 3 };
            let start = (original_start + s.len() - origin) % s.len();
            assert_eq!(site.recognition_start, start);
            assert_eq!(site.reverse, reverse);
            let (a, b) = if reverse {
                (start as i64 - 5, start as i64 - 1)
            } else {
                (start as i64 + 7, start as i64 + 11)
            };
            assert_eq!(
                (site.top_cut, site.bottom_cut),
                (
                    Some(a.rem_euclid(s.len() as i64) as usize),
                    Some(b.rem_euclid(s.len() as i64) as usize)
                )
            );
            let positions: Vec<_> = site.spans.iter().flat_map(|p| p.start..p.end).collect();
            assert_eq!(
                positions,
                (0..6).map(|j| (start + j) % s.len()).collect::<Vec<_>>()
            );
            let p = product("source", &m, &es, &r.options, 0).unwrap();
            assert_eq!(p.molecule.topology, Topology::Linear);
            assert_eq!(p.molecule.sequence.len(), s.len());
        }
    }
}
#[test]
fn possible_recognition_and_ambiguous_substrate_never_assert_physical_products() {
    for s in ["TTTGAANTCAAA", "NTTGAATTCAAA"] {
        let m = dna(s, Topology::Linear);
        let r = run(&m, &[enzyme("EcoRI")]);
        assert!(r.scan_complete);
        assert!(!r.products_determined);
        if s.contains("GAANTC") {
            assert!(!r.sites[0].certain);
            assert!(
                r.sites[0]
                    .diagnostics
                    .contains(&SiteDiagnostic::PossibleRecognition)
            );
            assert!(r.ref_top_lengths.is_none());
        } else {
            assert!(r.sites[0].certain);
            assert_eq!(r.ref_top_lengths, Some(vec![8, 4]));
        }
    }
    let m = dna("AAAACCCC", Topology::Circular);
    let r = run(&m, &[enzyme("EcoRI")]);
    assert!(r.products_determined);
    assert_eq!(r.fragments[0].duplex.topology, Topology::Circular);
    assert!(r.fragments[0].duplex.ends.is_none());
}
#[test]
fn comparison_site_output_and_fragment_limits_remain_explicit() {
    let m = dna("TTTGAATTCAAATTTGAATTCAAA", Topology::Linear);
    let es = [enzyme("EcoRI")];
    for o in [
        InspectionOptions {
            max_comparisons: 1,
            ..Default::default()
        },
        InspectionOptions {
            max_sites: 1,
            ..Default::default()
        },
        InspectionOptions {
            max_fragments: 1,
            ..Default::default()
        },
    ] {
        let r = inspect("source", &m, &es, &o).unwrap();
        assert!(!r.products_determined && r.fragments.is_empty() && !r.limits.is_empty());
        r.validate().unwrap();
        assert!(product("source", &m, &es, &o, 0).is_err());
    }
    let long = dna(&"GAATTC".repeat(100), Topology::Circular);
    let o = InspectionOptions {
        max_report_bytes: 4096,
        ..Default::default()
    };
    let r = inspect("source", &long, &es, &o).unwrap();
    assert!(r.limits.contains(&InspectionLimit::ReportBytes));
    assert!(!r.products_determined);
    assert!(serde_json::to_vec(&r).unwrap().len() <= 4096);
}
#[test]
fn multi_enzyme_duplicate_and_conflicting_cuts_do_not_invent_unique_fragments() {
    let m = dna("AAAGGTCTCTAGTCCCCCAA", Topology::Circular);
    let a = enzyme("BsaI");
    let mut b = a.clone();
    b.name = "Other declared geometry".into();
    let r = run(&m, &[a.clone(), b.clone()]);
    assert_eq!(r.sites.len(), 2);
    assert!(r.products_determined);
    assert_eq!(r.fragments.len(), 1);
    b.bottom_cut = 12;
    let r = run(&m, &[a, b]);
    assert!(r.scan_complete && !r.products_determined);
    assert!(
        r.diagnostics
            .iter()
            .any(|s| s.contains("conflicting bottom"))
    );
}
#[test]
fn annotation_projection_keeps_clipped_review_and_unresolved_source_evidence() {
    let mut m = dna("TTTGAATTCAAA", Topology::Linear);
    m.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Range {
            start: 1,
            end: 7,
            before: false,
            after: false,
        },
        qualifiers: vec![("label".into(), Some("cut feature".into()))],
        needs_review: false,
    });
    m.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Unresolved("REMOTE.1:1..3".into()),
        qualifiers: vec![],
        needs_review: false,
    });
    let p = product(
        "source",
        &m,
        &[enzyme("EcoRI")],
        &InspectionOptions::default(),
        0,
    )
    .unwrap();
    assert!(p.molecule.features.iter().any(|f| f.needs_review));
    assert!(!p.projection_findings.is_empty());
    assert!(
        p.molecule
            .metadata
            .keys()
            .any(|k| k.starts_with("annotation_mapping:unmapped:"))
    );
    assert!(p.molecule.original_genbank.is_none());
    p.validate().unwrap();
}
#[test]
fn historical_shape_and_exact_current_binding_reject_tamper_without_replaying() {
    let m = dna("TTTGAATTCAAA", Topology::Linear);
    let es = [enzyme("EcoRI")];
    let o = InspectionOptions::default();
    let r = run(&m, &es);
    r.validate_bindings_for("source", &m, &es, &o).unwrap();
    assert!(r.validate_bindings_for("different", &m, &es, &o).is_err());
    let mut wrong = r.clone();
    wrong.sites[0].bottom_cut = Some(3);
    assert!(wrong.validate().is_err());
    let mut wrong = r.clone();
    wrong.fragments[0].duplex.top_length += 1;
    assert!(wrong.validate().is_err());
    let p = product("source", &m, &es, &o, 0).unwrap();
    let mut wrong = p.clone();
    wrong.molecule.metadata.insert(
        "restriction_inspection:source_revision".into(),
        "0".repeat(64),
    );
    assert!(wrong.validate().is_err());
    assert!(p.validate_bindings_for("source", &m, &es, &o, 1).is_err());
    let mut es2 = es.clone();
    es2[0].bottom_cut = 4;
    assert_ne!(
        batch_revision("source", &m, &es, &o).unwrap(),
        batch_revision("source", &m, &es2, &o).unwrap()
    );
    assert!(r.validate_bindings_for("source", &m, &es2, &o).is_err());
}

#[test]
fn adversarial_historical_shapes_are_rejected_before_descriptor_or_value_cloning() {
    let m = dna("TTTGAATTCAAA", Topology::Linear);
    let es = [enzyme("EcoRI")];
    let report = run(&m, &es);
    let mut oversized = report.clone();
    oversized.enzymes = vec![report.enzymes[0].clone(); 10_001];
    assert!(oversized.validate().is_err());
    let mut oversized = report.clone();
    oversized.enzymes[0].recognition = "A".repeat(4096);
    assert!(oversized.validate().is_err());
    let mut oversized = report.clone();
    oversized.sites[0].spans = vec![ReferenceSpan { start: 0, end: 1 }; 10_001];
    assert!(oversized.validate().is_err());
    let mut oversized = report.clone();
    oversized.sites[0].diagnostics = vec![SiteDiagnostic::PossibleRecognition; 10_001];
    assert!(oversized.validate().is_err());
    let mut oversized = product("source", &m, &es, &InspectionOptions::default(), 0).unwrap();
    oversized.enzymes = vec![report.enzymes[0].clone(); 10_001];
    assert!(oversized.validate().is_err());
}

#[test]
fn palindromic_asymmetric_declared_cleavage_cannot_choose_a_strand_orientation() {
    let m = dna("TTTGAATTCAAAGGATCCAAA", Topology::Linear);
    let mut ambiguous = enzyme("EcoRI");
    ambiguous.name = "Declared asymmetric palindrome".into();
    ambiguous.bottom_cut = 4;
    let es = [ambiguous, enzyme("BamHI")];
    let r = run(&m, &es);
    assert!(r.scan_complete);
    assert_eq!(r.sites.len(), 2);
    let s = &r.sites[0];
    assert!(s.certain);
    assert_eq!((s.top_boundary, s.bottom_boundary), (4, 7));
    assert_eq!((s.top_cut, s.bottom_cut), (None, None));
    assert_eq!(s.diagnostics, vec![SiteDiagnostic::AmbiguousStrandGeometry]);
    assert!(r.sites[1].top_cut.is_some() && r.sites[1].bottom_cut.is_some());
    assert!(!r.products_determined && r.fragments.is_empty() && r.ref_top_lengths.is_none());
    assert!(product("source", &m, &es, &r.options, 0).is_err());
    r.validate_bindings_for("source", &m, &es, &r.options)
        .unwrap();
}

#[test]
fn complete_fragment_set_must_conserve_both_source_strands() {
    let m = dna("TTTGAATTCAAA", Topology::Linear);
    let mut r = run(&m, &[enzyme("EcoRI")]);
    r.fragments[0].duplex.bottom_length += 1;
    assert!(r.validate().is_err());
}

#[test]
fn native_recipe_and_annotated_text_roundtrip_retain_only_unverified_producer_claims() {
    use collagene_core::formats::{self, TextFormat};
    let raw = "LOCUS       source_record             12 bp    DNA     linear   SYN 01-JAN-2026\nDEFINITION  Independent source with reference coordinates.\nACCESSION   AB123456\nVERSION     AB123456.2\nREFERENCE   1  (bases 1 to 12)\n  AUTHORS   Independent,A.\n  TITLE     Source-only reference\n  JOURNAL   Unpublished\nFEATURES             Location/Qualifiers\n     misc_feature    2..7\n                     /label=\"clipped annotation\"\nORIGIN\n        1 tttgaattca aa\n//\n";
    let m = formats::import_genbank(raw).unwrap().molecules.remove(0);
    let p = product(
        "source",
        &m,
        &[enzyme("EcoRI")],
        &InspectionOptions::default(),
        0,
    )
    .unwrap();
    let recipe = unverified_recipe(&p.molecule).unwrap().unwrap();
    assert_eq!(recipe.schema_version, 1);
    assert_eq!(recipe.enzymes[0].recognition, "GAATTC");
    assert_eq!(
        (recipe.enzymes[0].top_cut, recipe.enzymes[0].bottom_cut),
        (1, 5)
    );
    assert_eq!(recipe.fragment_index, 0);
    assert_eq!(recipe.source_role, "source");
    assert_eq!(recipe.source_revision, m.revision());
    assert_eq!(recipe.batch_revision, p.batch_revision);
    let native = collagene_core::protocol::save_document(&p.molecule).unwrap();
    let restored = collagene_core::protocol::load_document(&native).unwrap();
    assert_eq!(
        restored.metadata[RECIPE_KEY],
        p.molecule.metadata[RECIPE_KEY]
    );
    assert!(chemistry::Duplex::from_molecule(&restored).is_err());
    for format in [TextFormat::GenBank, TextFormat::Embl] {
        let exported = formats::export_with_report(&p.molecule, format).unwrap();
        assert!(
            exported
                .warnings
                .iter()
                .any(|w| w.contains("unverified producer"))
        );
        assert!(!exported.text.contains("Source-only reference"));
        let imported = formats::import_text(&exported.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(imported.sequence, p.molecule.sequence);
        assert_eq!(imported.features.len(), p.molecule.features.len());
        assert!(imported.features.iter().any(|f| f.needs_review));
        assert_eq!(
            imported.metadata["source:restriction_inspection:recipe"],
            p.molecule.metadata[RECIPE_KEY]
        );
        assert_eq!(
            imported.metadata["source:restriction_inspection:duplex_summary"],
            p.molecule.metadata["duplex_model"]
        );
        assert!(!imported.metadata.contains_key(RECIPE_KEY));
        assert!(!imported.metadata.contains_key("duplex_model"));
        assert_eq!(
            imported.metadata["source:legacy:chemistry_unmapped"],
            "true"
        );
        assert!(chemistry::Duplex::from_molecule(&imported).is_err());
    }
    let mut tampered = p.clone();
    let mut value: serde_json::Value =
        serde_json::from_str(&tampered.molecule.metadata[RECIPE_KEY]).unwrap();
    value["fragment_index"] = serde_json::json!(1);
    tampered
        .molecule
        .metadata
        .insert(RECIPE_KEY.into(), value.to_string());
    assert!(tampered.validate().is_err());
    // A matching sentence inside an ordinary feature note is not the exact COMMENT marker.
    let mut ordinary = m.clone();
    ordinary.features[0].qualifiers.push((
        "note".into(),
        Some("Collagene physical product v1: unverified sequence-only evidence".into()),
    ));
    let text = formats::export_genbank(&ordinary).unwrap();
    let imported = formats::import_genbank(&text).unwrap().molecules.remove(0);
    assert!(
        !imported
            .metadata
            .contains_key("source:legacy:chemistry_unmapped")
    );
    assert!(chemistry::Duplex::from_molecule(&imported).is_ok());
}
