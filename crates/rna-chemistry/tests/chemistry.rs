//! Authored stoichiometric oracles and published neutral monomer checks.
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
use collagene_rna_chemistry::*;
fn source(sequence: &str, topology: Topology) -> Molecule {
    Molecule {
        name: "RNA fixture".into(),
        sequence: Sequence::new(Alphabet::Rna, sequence).unwrap(),
        topology,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn evidence() -> Provenance {
    Provenance {
        statement: "Explicit authored chemical model; not measured state".into(),
        reference: None,
    }
}
fn declaration(m: &Molecule) -> Declaration {
    Declaration {
        version: DECLARATION_VERSION,
        source_revision: source_revision(m).unwrap(),
        coverage: Coverage::ClosedModel,
        ends: if m.topology == Topology::Circular {
            Ends::Circular
        } else {
            Ends::Linear {
                five: FivePrime::Hydroxyl,
                three: ThreePrime::Hydroxyl,
            }
        },
        sites: vec![],
        links: vec![],
        annotation_reviews: vec![],
        provenance: evidence(),
    }
}
fn site(index: usize, nucleoside: Nucleoside, sugar: Sugar) -> Site {
    Site {
        index,
        nucleoside,
        sugar,
        provenance: evidence(),
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{a} != {b}");
}
#[test]
fn neutral_canonical_nucleosides_match_independently_published_formulas_and_masses() {
    // ChEBI 16335,17562,16750,16704. Values are published rounded quantities,
    // not implementation-generated expected results.
    for (base, formula, mono, avg) in [
        ("A", "C10H13N5O4", 267.09675, 267.245),
        ("C", "C9H13N3O5", 243.08552, 243.219),
        ("G", "C10H13N5O5", 283.09167, 283.244),
        ("U", "C9H12N2O6", 244.06954, 244.203),
    ] {
        let m = source(base, Topology::Linear);
        let d = declaration(&m);
        let r = analyze(&m, &d).unwrap();
        assert_eq!(source_revision(&m).unwrap(), m.revision());
        r.validate_binding(&r.context_revision).unwrap();
        r.validate_for(&m, &d).unwrap();
        assert!(r.identity_complete);
        assert_eq!(r.formula.unwrap().hill(), formula);
        let mass = r.mass.unwrap();
        close(mass.monoisotopic_da, mono, 0.000006);
        close(mass.average_da, avg, 0.000001);
    }
}
#[test]
fn every_short_linear_and_circular_reference_matches_a_nucleotide_condensation_oracle() {
    // Independent route: sum free 5'-monophosphate nucleotides, then remove
    // waters for polymerization and remove the unused terminal phosphate on
    // the linear 5'-OH species. Production instead sums nucleosides and adds bonds.
    let free = [
        [10_i64, 14, 5, 7, 1],
        [9, 14, 3, 8, 1],
        [10, 14, 5, 8, 1],
        [9, 13, 2, 9, 1],
    ];
    for len in 1..=4 {
        for word in 0..4usize.pow(len as u32) {
            let mut sequence = String::new();
            let mut n = word;
            let mut total = [0_i64; 5];
            for _ in 0..len {
                let i = n % 4;
                n /= 4;
                sequence.push(b"ACGU"[i] as char);
                for k in 0..5 {
                    total[k] += free[i][k];
                }
            }
            for topology in [Topology::Linear, Topology::Circular] {
                let mut expected = total;
                let waters = if topology == Topology::Circular {
                    len
                } else {
                    len - 1
                };
                expected[1] -= 2 * waters as i64;
                expected[3] -= waters as i64;
                if topology == Topology::Linear {
                    expected[1] -= 1;
                    expected[3] -= 3;
                    expected[4] -= 1;
                }
                let m = source(&sequence, topology);
                let r = analyze(&m, &declaration(&m)).unwrap();
                let f = r.formula.unwrap();
                assert_eq!(
                    [f.c, f.h, f.n, f.o, f.p],
                    expected.map(|n| n as u64),
                    "{sequence} {topology:?}"
                );
            }
        }
    }
    let m = source("AC", Topology::Linear);
    assert_eq!(
        analyze(&m, &declaration(&m))
            .unwrap()
            .formula
            .unwrap()
            .hill(),
        "C19H25N8O11P"
    );
}
#[test]
fn modified_nucleosides_remain_distinct_even_when_isobaric() {
    for (base, kind, formula, mono) in [
        ("A", Nucleoside::N6Methyladenosine, "C11H15N5O4", 281.11240),
        ("C", Nucleoside::FiveMethylcytidine, "C10H15N3O5", 257.10117),
        ("U", Nucleoside::FiveMethyluridine, "C10H14N2O6", 258.08519),
        ("U", Nucleoside::Pseudouridine, "C9H12N2O6", 244.06954),
        ("A", Nucleoside::Inosine, "C10H12N4O5", 268.08077),
    ] {
        let m = source(base, Topology::Linear);
        let plain = analyze(&m, &declaration(&m)).unwrap();
        let mut d = declaration(&m);
        d.sites.push(site(0, kind, Sugar::Ribose));
        let r = analyze(&m, &d).unwrap();
        assert_eq!(r.formula.unwrap().hill(), formula);
        close(r.mass.as_ref().unwrap().monoisotopic_da, mono, 0.000006);
        assert_ne!(r.context_revision, plain.context_revision);
        if base == "U" && formula == "C9H12N2O6" {
            assert_eq!(r.formula, plain.formula);
            assert_ne!(r.declaration, plain.declaration);
        }
    }
    let m = source("C", Topology::Linear);
    let mut a = declaration(&m);
    a.sites
        .push(site(0, Nucleoside::FiveMethylcytidine, Sugar::Ribose));
    let mut b = declaration(&m);
    b.sites
        .push(site(0, Nucleoside::Canonical, Sugar::TwoPrimeOMethyl));
    let ra = analyze(&m, &a).unwrap();
    let rb = analyze(&m, &b).unwrap();
    assert_eq!(ra.formula, rb.formula);
    assert_ne!(ra.context_revision, rb.context_revision);
}
#[test]
fn sugar_and_backbone_changes_use_atom_substitutions_with_separate_stereo_identity() {
    let m = source("C", Topology::Linear);
    for (sugar, formula) in [
        (Sugar::TwoPrimeOMethyl, "C10H15N3O5"),
        (Sugar::TwoPrimeDeoxy, "C9H13N3O4"),
        (Sugar::TwoPrimeFluoro, "C9H12FN3O4"),
    ] {
        let mut d = declaration(&m);
        d.sites.push(site(0, Nucleoside::Canonical, sugar));
        assert_eq!(analyze(&m, &d).unwrap().formula.unwrap().hill(), formula);
    }
    let m = source("U", Topology::Linear);
    let mut d = declaration(&m);
    d.sites
        .push(site(0, Nucleoside::FiveMethyluridine, Sugar::TwoPrimeDeoxy));
    assert_eq!(
        analyze(&m, &d).unwrap().formula.unwrap().hill(),
        "C10H14N2O5"
    );
    let m = source("AC", Topology::Linear);
    let plain = analyze(&m, &declaration(&m)).unwrap();
    let mut reports = vec![];
    for stereo in [Stereo::Rp, Stereo::Sp, Stereo::Unspecified] {
        let mut d = declaration(&m);
        d.links.push(Link {
            after: 0,
            linkage: Linkage::Phosphorothioate {
                stereo: stereo.clone(),
            },
            provenance: evidence(),
        });
        let r = analyze(&m, &d).unwrap();
        assert_eq!(r.formula.unwrap().hill(), "C19H25N8O10PS");
        assert!(r.formula_complete);
        assert_eq!(r.identity_complete, stereo != Stereo::Unspecified);
        close(
            r.mass.as_ref().unwrap().monoisotopic_da - plain.mass.as_ref().unwrap().monoisotopic_da,
            15.97715655483,
            1e-9,
        );
        reports.push(r);
    }
    assert_eq!(reports[0].formula, reports[1].formula);
    assert_ne!(reports[0].context_revision, reports[1].context_revision);
    assert_eq!(reports[2].unspecified_stereocenters, vec![0]);
}
#[test]
fn terminal_phosphates_and_circle_seams_account_for_complete_condensation() {
    let m = source("A", Topology::Linear);
    let mut d = declaration(&m);
    d.ends = Ends::Linear {
        five: FivePrime::Monophosphate,
        three: ThreePrime::Hydroxyl,
    };
    let amp = analyze(&m, &d).unwrap();
    assert_eq!(amp.formula.unwrap().hill(), "C10H14N5O7P");
    close(amp.mass.unwrap().monoisotopic_da, 347.06308, 0.000006);
    d.ends = Ends::Linear {
        five: FivePrime::Triphosphate,
        three: ThreePrime::Hydroxyl,
    };
    let atp = analyze(&m, &d).unwrap();
    assert_eq!(atp.formula.unwrap().hill(), "C10H16N5O13P3");
    close(atp.mass.unwrap().monoisotopic_da, 506.99575, 0.000006);
    d.ends = Ends::Linear {
        five: FivePrime::Hydroxyl,
        three: ThreePrime::TwoThreeCyclicPhosphate,
    };
    let terminal = analyze(&m, &d).unwrap();
    let ring = source("A", Topology::Circular);
    let circle = analyze(&ring, &declaration(&ring)).unwrap();
    assert_eq!(terminal.formula, circle.formula);
    assert_ne!(terminal.context_revision, circle.context_revision);
    d.sites
        .push(site(0, Nucleoside::Canonical, Sugar::TwoPrimeOMethyl));
    assert!(analyze(&m, &d).is_err());
    d.sites[0].sugar = Sugar::Unknown {
        label: "2' chemistry unresolved".into(),
    };
    assert!(!analyze(&m, &d).unwrap().formula_complete);
    let mut d = declaration(&ring);
    d.links.push(Link {
        after: 0,
        linkage: Linkage::Phosphorothioate { stereo: Stereo::Rp },
        provenance: evidence(),
    });
    assert_eq!(
        analyze(&ring, &d).unwrap().formula.unwrap().hill(),
        "C10H12N5O5PS"
    );
}
#[test]
fn partial_unknown_ambiguous_and_retained_foreign_chemistry_never_become_unmodified_mass() {
    let m = source("AC", Topology::Linear);
    let mut declarations = vec![];
    let mut d = declaration(&m);
    d.coverage = Coverage::Partial;
    declarations.push(d);
    let mut d = declaration(&m);
    d.sites.push(site(
        0,
        Nucleoside::Unknown {
            label: "m6A spelling alone is not a typed declaration".into(),
        },
        Sugar::Ribose,
    ));
    declarations.push(d);
    let mut d = declaration(&m);
    d.links.push(Link {
        after: 0,
        linkage: Linkage::Unknown {
            label: "unsupported backbone".into(),
        },
        provenance: evidence(),
    });
    declarations.push(d);
    let mut d = declaration(&m);
    d.ends = Ends::Linear {
        five: FivePrime::Unknown {
            label: "cap identity/protonation not modeled".into(),
        },
        three: ThreePrime::Hydroxyl,
    };
    declarations.push(d);
    for d in declarations {
        let r = analyze(&m, &d).unwrap();
        assert!(!r.formula_complete);
        assert!(!r.identity_complete);
        assert!(r.formula.is_none() && r.mass.is_none());
        assert!(!r.issues.is_empty());
        r.validate_binding(&r.context_revision).unwrap();
        r.validate_for(&m, &d).unwrap();
    }
    let m = source(&"N".repeat(120), Topology::Linear);
    let r = analyze(&m, &declaration(&m)).unwrap();
    assert_eq!(r.issues[0].kind, IssueKind::AmbiguousReference);
    assert_eq!(r.issues[0].count, 120);
    assert_eq!(r.issues[0].positions.len(), 100);
    assert!(r.issues[0].positions_truncated);
    assert!(r.mass.is_none());
    let mut m = source("AC", Topology::Linear);
    m.features.push(Feature {
        kind: "modified_base".into(),
        location: Location::range(0, 1),
        qualifiers: vec![("mod_base".into(), Some("m6A".into()))],
        needs_review: false,
    });
    m.metadata.insert(
        "rna:modifications".into(),
        "retained source notation".into(),
    );
    let r = analyze(&m, &declaration(&m)).unwrap();
    assert!(r.mass.is_none());
    assert_eq!(r.retained_evidence.len(), 2);
    assert!(
        r.issues
            .iter()
            .any(|i| i.kind == IssueKind::RetainedSourceEvidence)
    );
}
#[test]
fn invalid_parents_indices_versions_topologies_and_budgets_are_refused_before_model_claims() {
    let m = source("AC", Topology::Linear);
    let mut d = declaration(&m);
    d.sites
        .push(site(0, Nucleoside::FiveMethylcytidine, Sugar::Ribose));
    assert!(analyze(&m, &d).is_err());
    d = declaration(&m);
    d.sites.push(site(0, Nucleoside::Canonical, Sugar::Ribose));
    d.sites
        .push(site(0, Nucleoside::Canonical, Sugar::TwoPrimeOMethyl));
    assert!(analyze(&m, &d).is_err());
    d = declaration(&m);
    d.links.push(Link {
        after: 1,
        linkage: Linkage::Phosphodiester,
        provenance: evidence(),
    });
    assert!(analyze(&m, &d).is_err());
    for version in [0, 1, 3] {
        let mut d = declaration(&m);
        d.version = version;
        assert!(analyze(&m, &d).is_err());
    }
    let mut d = declaration(&m);
    d.ends = Ends::Circular;
    assert!(analyze(&m, &d).is_err());
    d = declaration(&m);
    d.source_revision = "f".repeat(64);
    assert_eq!(analyze(&m, &d).unwrap_err(), Error::Binding);
    let dna = Molecule::dna("DNA", "AC", Topology::Linear).unwrap();
    assert!(source_revision(&dna).is_err());
    let too_long = source(&"A".repeat(MAX_BASES + 1), Topology::Linear);
    assert!(source_revision(&too_long).is_err());
    let mut value = serde_json::to_value(declaration(&m)).unwrap();
    value.as_object_mut().unwrap().remove("coverage");
    assert!(serde_json::from_value::<Declaration>(value).is_err());
}
#[test]
fn report_validation_reconstructs_source_declarations_catalog_formula_and_quantities() {
    let m = source("ACGU", Topology::Linear);
    let d = declaration(&m);
    let r = analyze(&m, &d).unwrap();
    let roundtrip: Report = serde_json::from_slice(&serde_json::to_vec(&r).unwrap()).unwrap();
    roundtrip.validate_for(&m, &d).unwrap();
    roundtrip.validate_binding(&r.context_revision).unwrap();
    assert_eq!(request_revision(&m, &d).unwrap(), r.context_revision);
    let mut changed = r.clone();
    changed.engine = "rna-chemistry-0".into();
    assert!(changed.validate_for(&m, &d).is_err());
    assert!(changed.validate_binding(&r.context_revision).is_err());
    let mut changed = r.clone();
    changed.formula.as_mut().unwrap().p += 1;
    assert!(changed.validate_for(&m, &d).is_err());
    assert!(changed.validate_binding(&r.context_revision).is_err());
    let mut changed = r.clone();
    changed.mass.as_mut().unwrap().monoisotopic_da += 1.0;
    assert!(changed.validate_for(&m, &d).is_err());
    assert!(changed.validate_binding(&r.context_revision).is_err());
    let mut changed = r.clone();
    changed.identity_complete = false;
    assert!(changed.validate_for(&m, &d).is_err());
    assert!(changed.validate_binding(&r.context_revision).is_err());
    let mut changed = r.clone();
    changed.catalog_revision = "0".repeat(64);
    assert!(changed.validate_for(&m, &d).is_err());
    assert!(changed.validate_binding(&r.context_revision).is_err());
    let mut changed = d.clone();
    changed.provenance.statement = "Another explicit source hypothesis".into();
    assert!(r.validate_for(&m, &changed).is_err());
    let changed = source("AGGU", Topology::Linear);
    assert!(r.validate_for(&changed, &d).is_err());
}

fn annotated(base: &str, code: &str) -> Molecule {
    let mut m = source(base, Topology::Linear);
    m.features.push(Feature {
        kind: "modified_base".into(),
        location: Location::range(0, 1),
        qualifiers: vec![("mod_base".into(), Some(code.into()))],
        needs_review: false,
    });
    m
}
fn reviewed(m: &Molecule) -> Declaration {
    let c = annotation_candidates(m).unwrap().remove(0);
    let mut d = declaration(m);
    d.sites.push(c.suggested_site.unwrap());
    d.annotation_reviews.push(FeatureReview {
        feature_index: c.feature_index,
        expected_feature_revision: c.feature_revision,
        site_entry: 0,
    });
    d
}
#[test]
fn insdc_annotation_review_constrains_chemistry_without_erasing_original_evidence() {
    // INSDC 11.4 codes, checked against independently specified monomer formulas.
    for (base, code, formula) in [
        ("A", "m6a", "C11H15N5O4"),
        ("C", "m5c", "C10H15N3O5"),
        ("U", "m5u", "C10H14N2O6"),
        ("U", "p", "C9H12N2O6"),
        ("A", "i", "C10H12N4O5"),
        ("C", "cm", "C10H15N3O5"),
        ("G", "gm", "C11H15N5O5"),
        ("U", "um", "C10H14N2O6"),
        ("U", "fm", "C10H14N2O6"),
        ("U", "tm", "C11H16N2O6"),
    ] {
        let m = annotated(base, code);
        let before = m.clone();
        assert!(analyze(&m, &declaration(&m)).unwrap().mass.is_none());
        let d = reviewed(&m);
        let r = analyze(&m, &d).unwrap();
        assert_eq!(r.formula.unwrap().hill(), formula, "{code}");
        assert_eq!(
            r.retained_evidence,
            vec![RetainedEvidence::Feature { index: 0 }]
        );
        assert!(r.unreviewed_evidence.is_empty());
        r.validate_binding(&r.context_revision).unwrap();
        r.validate_for(&m, &d).unwrap();
        assert_eq!(m, before);
    }
    // A base-only annotation permits explicitly declared extra sugar chemistry.
    let m = annotated("A", "m6a");
    let mut d = reviewed(&m);
    d.sites[0].sugar = Sugar::TwoPrimeOMethyl;
    assert_eq!(
        analyze(&m, &d).unwrap().formula.unwrap().hill(),
        "C12H17N5O4"
    );
    // A sugar-only annotation permits an explicitly declared base modification.
    let m = annotated("C", "cm");
    let mut d = reviewed(&m);
    d.sites[0].nucleoside = Nucleoside::FiveMethylcytidine;
    assert_eq!(
        analyze(&m, &d).unwrap().formula.unwrap().hill(),
        "C11H17N3O5"
    );
}
#[test]
fn review_cannot_hide_unknown_metadata_uncertain_locations_or_conflicting_chemistry() {
    for (base, code) in [("A", "OTHER"), ("U", "t"), ("U", "m1f"), ("C", "m6a")] {
        let m = annotated(base, code);
        let c = annotation_candidates(&m).unwrap().remove(0);
        assert!(c.suggested_site.is_none() && c.unsupported_reason.is_some());
        assert_eq!(c.raw_code.as_deref(), Some(code));
        let mut d = declaration(&m);
        d.sites.push(site(0, Nucleoside::Canonical, Sugar::Ribose));
        d.annotation_reviews.push(FeatureReview {
            feature_index: 0,
            expected_feature_revision: c.feature_revision,
            site_entry: 0,
        });
        assert!(analyze(&m, &d).is_err());
    }
    let mut m = annotated("A", "m6a");
    m.metadata
        .insert("rna:modifications".into(), "unresolved cap".into());
    let r = analyze(&m, &reviewed(&m)).unwrap();
    assert!(r.mass.is_none());
    assert_eq!(r.retained_evidence.len(), 2);
    assert_eq!(
        r.unreviewed_evidence,
        vec![RetainedEvidence::Metadata {
            key: "rna:modifications".into()
        }]
    );
    for location in [
        Location::Complement(Box::new(Location::range(0, 1))),
        Location::Range {
            start: 0,
            end: 1,
            before: true,
            after: false,
        },
        Location::Join(vec![Location::range(0, 1), Location::range(1, 2)]),
    ] {
        let mut m = annotated("AA", "m6a");
        m.features[0].location = location;
        assert!(
            annotation_candidates(&m).unwrap()[0]
                .suggested_site
                .is_none()
        );
    }
    let mut m = annotated("A", "m6a");
    m.features[0].needs_review = true;
    assert!(
        annotation_candidates(&m).unwrap()[0]
            .suggested_site
            .is_none()
    );
    for code in ["cm", "fm", "tm"] {
        let m = annotated(if code == "cm" { "C" } else { "U" }, code);
        let mut d = reviewed(&m);
        d.sites[0].sugar = Sugar::Ribose;
        assert!(analyze(&m, &d).is_err());
    }
    let m = annotated("U", "fm");
    let mut d = reviewed(&m);
    d.sites[0].nucleoside = Nucleoside::Canonical;
    assert!(analyze(&m, &d).is_err());
}
#[test]
fn review_is_bound_to_exact_feature_and_declared_site_and_is_never_implicit() {
    let m = annotated("AA", "m6a");
    let d = reviewed(&m);
    let mut variants = vec![];
    let mut changed = d.clone();
    changed.annotation_reviews[0].feature_index = 1;
    variants.push(changed);
    let mut changed = d.clone();
    changed.annotation_reviews[0].site_entry = 1;
    variants.push(changed);
    let mut changed = d.clone();
    changed.annotation_reviews[0].expected_feature_revision = "f".repeat(64);
    variants.push(changed);
    let mut changed = d.clone();
    changed
        .annotation_reviews
        .push(changed.annotation_reviews[0].clone());
    variants.push(changed);
    let mut changed = d.clone();
    changed.sites[0].index = 1;
    variants.push(changed);
    let mut changed = d.clone();
    changed.sites[0].nucleoside = Nucleoside::Canonical;
    variants.push(changed);
    for changed in variants {
        assert!(analyze(&m, &changed).is_err());
    }
    let mut changed = m.clone();
    changed.features[0]
        .qualifiers
        .push(("note".into(), Some("source corrected".into())));
    let mut rebound = d.clone();
    rebound.source_revision = source_revision(&changed).unwrap();
    assert!(analyze(&changed, &rebound).is_err()); // Source rebinding cannot repair a stale feature review.
    let mut no_review = d.clone();
    no_review.annotation_reviews.clear();
    assert!(analyze(&m, &no_review).unwrap().mass.is_none());
    let mut partial = d.clone();
    partial.coverage = Coverage::Partial;
    assert!(analyze(&m, &partial).unwrap().mass.is_none());
    let mut duplicated = m.clone();
    duplicated.features.push(duplicated.features[0].clone());
    let only_first = reviewed(&duplicated);
    let r = analyze(&duplicated, &only_first).unwrap();
    assert!(r.mass.is_none());
    assert_eq!(
        r.unreviewed_evidence,
        vec![RetainedEvidence::Feature { index: 1 }]
    );
}
