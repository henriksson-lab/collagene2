use collagene_core::{
    analysis,
    catalog_digest::*,
    chemistry, formats,
    restriction_history::{self, SourceBinding},
    restriction_inspection,
    workspace::{RevisionEvent, Workspace},
    *,
};
fn selected(c: &catalog::Catalog, ids: &[&str]) -> CatalogSelection {
    select(c, &ids.iter().map(|s| s.to_string()).collect::<Vec<_>>()).unwrap()
}
fn dna(s: &str, t: Topology) -> Molecule {
    Molecule::dna("independent catalog fixture", s, t).unwrap()
}
fn binding(w: &Workspace, i: usize) -> SourceBinding {
    let d = &w.documents[i];
    SourceBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn frame() -> Workspace {
    let raw = "LOCUS       catalog_source            12 bp    DNA     circular SYN 01-JAN-2026\nACCESSION   ZZ123456\nVERSION     ZZ123456.1\nFEATURES             Location/Qualifiers\n     misc_feature    join(11..12,1..2)\n                     /label=\"origin\"\nORIGIN\n        1 tttgaattca aa\n//\n";
    Workspace::from_documents(
        "catalog native fixture",
        formats::import_text(raw).unwrap().molecules,
        0,
    )
    .unwrap()
}
#[test]
fn original_twelve_use_identical_science_without_replacing_algorithms() {
    let c = builtin_catalog().unwrap();
    for e in analysis::starter_enzymes() {
        let s = selected(&c, &[&e.name.to_ascii_lowercase()]);
        assert_eq!(s.enzymes[0].enzyme.recognition, e.recognition);
        assert_eq!(
            (s.enzymes[0].enzyme.top_cut, s.enzymes[0].enzyme.bottom_cut),
            (e.top_cut, e.bottom_cut)
        );
        let m = dna(
            &format!("AAA{}AAAAAAAAAAAAAAAA", e.recognition),
            Topology::Linear,
        );
        let o = InspectionOptions::default();
        let legacy = restriction_inspection::inspect("source", &m, &[e.clone()], &o).unwrap();
        let wrapped = inspect(&c, &s, "source", &m, &o).unwrap();
        assert_eq!(
            serde_json::to_value(&legacy).unwrap(),
            serde_json::to_value(&wrapped.inspection).unwrap()
        );
        assert_ne!(wrapped.batch_revision, legacy.batch_revision);
        wrapped
            .validate_bindings_for(&c, &s, "source", &m, &o)
            .unwrap();
        if legacy.products_determined {
            let p = product(&c, &s, "source", &m, &o, 0).unwrap();
            let old = restriction_inspection::product("source", &m, &[e], &o, 0).unwrap();
            assert_eq!(p.product.molecule.sequence, old.molecule.sequence);
            assert_eq!(p.product.batch_revision, old.batch_revision);
            assert_eq!(unverified_selection(&p.product.molecule).unwrap(), Some(s));
        }
    }
}
#[test]
fn added_blunt_and_degenerate_identity_have_literal_cuts_and_real_ends() {
    let c = builtin_catalog().unwrap();
    let o = InspectionOptions::default();
    let s = selected(&c, &["alui"]);
    let m = dna("TTTAGCTAAA", Topology::Linear);
    let r = inspect(&c, &s, "source", &m, &o).unwrap();
    assert_eq!(
        (
            r.inspection.sites[0].top_cut,
            r.inspection.sites[0].bottom_cut
        ),
        (Some(5), Some(5))
    );
    let a = product(&c, &s, "source", &m, &o, 0).unwrap();
    let b = product(&c, &s, "source", &m, &o, 1).unwrap();
    assert_eq!(a.product.molecule.sequence.as_str(), "TTTAG");
    assert_eq!(b.product.molecule.sequence.as_str(), "CTAAA");
    assert_eq!(
        a.product.summary.duplex.ends.unwrap().1.polarity,
        chemistry::Polarity::Blunt
    );
    let s = selected(&c, &["hinfi"]);
    let m = dna("AAAGACTCAAA", Topology::Linear);
    let r = inspect(&c, &s, "source", &m, &o).unwrap();
    assert!(r.inspection.sites[0].certain);
    assert_eq!(
        (
            r.inspection.sites[0].top_cut,
            r.inspection.sites[0].bottom_cut
        ),
        (Some(4), Some(7))
    );
    let p = product(&c, &s, "source", &m, &o, 0).unwrap();
    // Left fragment protrudes on the bottom strand: reverse-complement of
    // reference ACT is AGT in that strand's own 5′→3′ direction.
    assert_eq!(p.product.summary.duplex.ends.unwrap().1.overhang, "AGT");
    let right = product(&c, &s, "source", &m, &o, 1).unwrap();
    assert_eq!(right.product.summary.duplex.ends.unwrap().0.overhang, "ACT");
}
#[test]
fn same_cut_identities_are_distinct_and_selection_never_expands_formulations_or_aliases() {
    let c = builtin_catalog().unwrap();
    let m = dna("AAACCGGTTT", Topology::Linear);
    let o = InspectionOptions::default();
    let a = selected(&c, &["hpaii"]);
    let b = selected(&c, &["mspi"]);
    let pa = product(&c, &a, "source", &m, &o, 0).unwrap();
    let pb = product(&c, &b, "source", &m, &o, 0).unwrap();
    assert_eq!(pa.product.molecule.sequence, pb.product.molecule.sequence);
    assert_ne!(pa.batch_revision, pb.batch_revision);
    assert_ne!(
        pa.product.molecule.revision(),
        pb.product.molecule.revision()
    );
    for id in ["EcoRI", "R3101", "EcoRI-HF", "not-an-id"] {
        assert!(select(&c, &[id.into()]).is_err());
    }
    assert!(select(&c, &["ecori".into(), "ecori".into()]).is_err());
    let mut e = c.edition().clone();
    for en in &mut e.enzymes[..2] {
        en.aliases.push(catalog::Alias {
            name: "shared alias".into(),
            shared: true,
            citations: vec![catalog::Citation {
                source: "neb-recognition".into(),
                locator: "independent synthetic ambiguous alias fixture".into(),
            }],
        });
    }
    let c = catalog::Catalog::from_edition(catalog::seal(e).unwrap()).unwrap();
    assert!(matches!(
        c.lookup_enzyme("shared alias"),
        catalog::Lookup::Ambiguous(_)
    ));
    assert!(select(&c, &["shared alias".into()]).is_err());
}
#[test]
fn ordered_identity_pin_rejects_forgery_changed_editions_and_stale_sources() {
    let c = builtin_catalog().unwrap();
    let s = selected(&c, &["ecori", "bamhi"]);
    let reverse = selected(&c, &["bamhi", "ecori"]);
    let m = dna("TTTGAATTCAAAGGATCCAAA", Topology::Linear);
    let o = InspectionOptions::default();
    assert_ne!(
        batch_revision(&c, &s, "source", &m, &o).unwrap(),
        batch_revision(&c, &reverse, "source", &m, &o).unwrap()
    );
    let r = inspect(&c, &s, "source", &m, &o).unwrap();
    let mut forged = s.clone();
    forged.enzymes[0].enzyme.top_cut += 1;
    assert!(source_preflight(&c, &forged, "source", &m, &o).is_err());
    let mut e = c.edition().clone();
    e.edition = "independent future fixture edition".into();
    let next = catalog::Catalog::from_edition(catalog::seal(e).unwrap()).unwrap();
    assert!(s.validate_for(&next).is_err());
    let new = selected(&next, &["ecori", "bamhi"]);
    assert_ne!(
        r.batch_revision,
        batch_revision(&next, &new, "source", &m, &o).unwrap()
    );
    assert!(
        r.validate_bindings_for(&next, &new, "source", &m, &o)
            .is_err()
    );
    assert!(
        r.validate_bindings_for(&c, &s, "different role", &m, &o)
            .is_err()
    );
    let mut renamed = m.clone();
    renamed.name = "another exact source revision".into();
    assert!(
        r.validate_bindings_for(&c, &s, "source", &renamed, &o)
            .is_err()
    );
    let mut p = product(&c, &s, "source", &m, &o, 0).unwrap();
    p.product.molecule.metadata.remove(RECIPE_KEY);
    p.product.molecule.metadata.insert(
        IMPORTED_RECIPE_KEY.into(),
        serde_json::to_string(&s).unwrap(),
    );
    assert_eq!(unverified_selection(&p.product.molecule).unwrap(), Some(s));
    assert!(p.validate().is_err());
}
#[test]
fn uncertainty_caps_and_physical_producer_claims_never_gain_catalog_authority() {
    let c = builtin_catalog().unwrap();
    let s = selected(&c, &["ecori"]);
    let o = InspectionOptions::default();
    let m = dna("TTTGAANTCAAA", Topology::Linear);
    let r = inspect(&c, &s, "source", &m, &o).unwrap();
    assert!(r.inspection.scan_complete && !r.inspection.products_determined);
    assert!(product(&c, &s, "source", &m, &o, 0).is_err());
    let mut m = dna("TTTGAATTCAAA", Topology::Linear);
    let low = InspectionOptions {
        max_comparisons: 1,
        ..o.clone()
    };
    assert!(
        !inspect(&c, &s, "source", &m, &low)
            .unwrap()
            .inspection
            .scan_complete
    );
    assert!(product(&c, &s, "source", &m, &low, 0).is_err());
    m.metadata.insert(
        "duplex_model".into(),
        "unverified producer claims are not an intact duplex".into(),
    );
    assert!(
        !inspect(&c, &s, "source", &m, &o)
            .unwrap()
            .inspection
            .products_determined
    );
    let mut excessive = s.clone();
    excessive.enzymes = vec![s.enzymes[0].clone(); 65];
    assert!(excessive.validate().is_err());
}
#[test]
fn native_catalog_install_adds_metadata_before_event_and_preserves_only_exact_parent() {
    let c = builtin_catalog().unwrap();
    let s = selected(&c, &["ecori"]);
    let w = frame();
    let b = binding(&w, 0);
    let o = InspectionOptions::default();
    let original = serde_json::to_value(&w).unwrap();
    let out = install(&c, &s, &w, &b, &o, 0).unwrap();
    let decoded = workspace::load_workspace(&out.installed.archive).unwrap();
    out.validate_archive_bindings_for(&c, &s, &w, &b, &o, 0, &decoded)
        .unwrap();
    out.validate_archive_for(&c, &s, &w, &b, &o, 0, &decoded)
        .unwrap();
    assert_eq!(serde_json::to_value(&w).unwrap(), original);
    let m = decoded
        .revisions
        .get(&out.installed.product_revision)
        .unwrap();
    assert_eq!(m.sequence.as_str(), "AATTCAAATTTG");
    assert_eq!(unverified_selection(m).unwrap(), Some(s));
    assert!(m.original_genbank.is_none());
    assert!(!m.features.is_empty());
    assert!(chemistry::Duplex::from_molecule(m).is_err());
    assert!(
        decoded
            .revisions
            .get(&b.revision)
            .unwrap()
            .original_genbank
            .is_some()
    );
    let event = decoded
        .history
        .iter()
        .find(|e| e.id == out.installed.event_id)
        .unwrap();
    assert_eq!(event.parents, vec![b.event_id]);
    assert_eq!(event.revision, m.revision());
}
#[test]
fn exact_current_source_event_and_identical_content_roles_are_required() {
    let c = builtin_catalog().unwrap();
    let s = selected(&c, &["ecori"]);
    let o = InspectionOptions::default();
    let w = frame();
    let original = w.revisions[&w.documents[0].head].clone();
    let mut twins =
        Workspace::from_documents("two roles", vec![original.clone(), original], 0).unwrap();
    let a = binding(&twins, 0);
    let b = binding(&twins, 1);
    assert_eq!(a.revision, b.revision);
    assert_ne!(
        restriction_history::source_role(&a).unwrap(),
        restriction_history::source_role(&b).unwrap()
    );
    assert!(prepare_install(&c, &s, &twins, &b, &o, 0).is_err());
    let rev_a = installation_revision(&c, &s, &twins, &a, &o, 0).unwrap();
    twins.active_document = b.document_id.clone();
    let rev_b = installation_revision(&c, &s, &twins, &b, &o, 0).unwrap();
    assert_ne!(rev_a, rev_b);
    let next = RevisionEvent::new(
        b.document_id.clone(),
        b.revision.clone(),
        vec![b.event_id.clone()],
        "same DNA, new current event",
    )
    .unwrap();
    twins.documents[1].head_event = Some(next.id.clone());
    twins.history.push(next);
    twins.validate().unwrap();
    assert!(prepare_install(&c, &s, &twins, &b, &o, 0).is_err());
}
#[test]
fn cheap_archive_request_crosscheck_rejects_mislabeled_fragment_without_replaying_science() {
    let c = builtin_catalog().unwrap();
    let s = selected(&c, &["ecori"]);
    let m = dna("TTTGAATTCAAA", Topology::Linear);
    let w = Workspace::from_documents("wrong fragment fixture", vec![m], 0).unwrap();
    let b = binding(&w, 0);
    let o = InspectionOptions::default();
    let mut wrong = install(&c, &s, &w, &b, &o, 1).unwrap();
    let decoded = workspace::load_workspace(&wrong.installed.archive).unwrap();
    let revision = installation_revision(&c, &s, &w, &b, &o, 0).unwrap();
    wrong.batch_revision = revision.clone();
    wrong.installed.batch_revision = revision;
    wrong.validate_archive(&decoded).unwrap();
    assert!(
        wrong
            .validate_archive_bindings_for(&c, &s, &w, &b, &o, 0, &decoded)
            .is_err()
    );
    assert!(
        wrong
            .validate_archive_for(&c, &s, &w, &b, &o, 0, &decoded)
            .is_err()
    );
    let different = InspectionOptions {
        max_sites: 99,
        ..o.clone()
    };
    assert!(
        wrong
            .validate_archive_bindings_for(&c, &s, &w, &b, &different, 0, &decoded)
            .is_err()
    );
}
