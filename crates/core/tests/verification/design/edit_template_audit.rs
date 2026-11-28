//! Literal independent edited references, supplied-donor mutations and coordinate walks.
use collagene_core::{crispr, edit_template::*, *};
use collagene_model::{Error, Feature, Location, Molecule, Topology};
fn mol(s: &str, circle: bool) -> Molecule {
    Molecule::dna(
        "source",
        s,
        if circle {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )
    .unwrap()
}
fn opts(
    start: usize,
    delete: usize,
    replacement: &str,
    left: usize,
    right: usize,
) -> EditTemplateOptions {
    EditTemplateOptions {
        edit_start: start,
        delete_length: delete,
        replacement: replacement.into(),
        left_arm: left,
        right_arm: right,
        reverse_donor: false,
        max_cut_distance: 100,
    }
}
fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => panic!(),
        })
        .collect()
}
#[test]
fn literal_linear_substitution_insertion_deletion_and_reverse_donors() {
    let s = mol("AACCGGTTAACC", false);
    for (o, target, donor) in [
        (opts(4, 4, "AG", 3, 3), "AACCAGAACC", "ACCAGAAC"),
        (opts(4, 4, "", 3, 3), "AACCAACC", "ACCAAC"),
        (opts(4, 0, "AG", 3, 3), "AACCAGGGTTAACC", "ACCAGGGT"),
    ] {
        for reverse in [false, true] {
            let mut o = o.clone();
            o.reverse_donor = reverse;
            let t = mol(target, false);
            let r = design(&s, &t, &o, None).unwrap();
            assert!(r.check.consistent);
            assert_eq!(r.edited_reference.sequence.as_str(), target);
            assert_eq!(
                r.donor.sequence.as_str(),
                if reverse { rc(donor) } else { donor.into() }
            );
            assert_eq!(r.left.source_positions, vec![1, 2, 3]);
            assert_eq!(
                r.right.source_positions,
                (4 + o.delete_length..7 + o.delete_length).collect::<Vec<_>>()
            );
            replay(&s, &t, &o, None, &r).unwrap();
            assert_eq!(s.sequence.as_str(), "AACCGGTTAACC");
        }
    }
}
#[test]
fn origin_crossing_edit_has_explicit_new_origin_and_nonoverlapping_homology() {
    let s = mol("AACCGGTT", true);
    let t = mol("ACCGGGA", true);
    let mut o = opts(6, 3, "GA", 2, 2);
    let r = design(&s, &t, &o, None).unwrap();
    assert!(r.check.consistent);
    assert_eq!(r.edited_reference.sequence.as_str(), "GAACCGG");
    assert_eq!(r.donor.sequence.as_str(), "GGGAAC");
    assert_eq!(r.deleted_positions, vec![6, 7, 0]);
    assert_eq!(r.left.source_positions, vec![4, 5]);
    assert_eq!(r.right.source_positions, vec![1, 2]);
    o.reverse_donor = true;
    assert_eq!(
        design(&s, &t, &o, None).unwrap().donor.sequence.as_str(),
        "GTTCCC"
    );
    o.left_arm = 4;
    assert!(source_preflight(&s, &t, &o, None).is_err());
}
#[test]
fn independent_target_never_supplies_generated_bases_and_foreign_reports_cannot_pass() {
    let s = mol("AACCGGTTAACC", false);
    let good = mol("AACCAGAACC", false);
    let wrong = mol("AAAAAAAAAA", false);
    let o = opts(4, 4, "AG", 3, 3);
    let a = design(&s, &good, &o, None).unwrap();
    let b = design(&s, &wrong, &o, None).unwrap();
    assert!(a.check.consistent);
    assert!(!b.check.consistent && !b.check.target_matches);
    assert_eq!(a.donor.sequence, b.donor.sequence);
    assert_eq!(a.edited_reference.sequence, b.edited_reference.sequence);
    let mut edited = a.clone();
    edited.check.target_matches = false;
    assert!(edited.validate().is_err());
    assert!(a.validate_bindings_for(&s, &wrong, &o, None).is_err());
    let mut forged = a.clone();
    forged.donor.name = "forged producer name".into();
    forged.validate_bindings_for(&s, &good, &o, None).unwrap();
    assert!(replay(&s, &good, &o, None, &forged).is_err());
}
#[test]
fn independent_checker_detects_each_arm_replacement_length_and_orientation_error() {
    let s = mol("AACCGGTTAACC", false);
    let t = mol("AACCAGAACC", false);
    let o = opts(4, 4, "AG", 3, 3);
    let r = design(&s, &t, &o, None).unwrap();
    for i in 0..r.donor.sequence.len() {
        let mut bytes = r.donor.sequence.as_str().as_bytes().to_vec();
        bytes[i] = if bytes[i] == b'A' { b'C' } else { b'A' };
        let text = String::from_utf8(bytes).unwrap();
        let c = check(&s, &t, &o, None, &text).unwrap();
        assert!(!c.consistent);
        assert_eq!(c.left_matches, i >= 3);
        assert_eq!(c.replacement_matches, !(3..5).contains(&i));
        assert_eq!(c.right_matches, i < 5);
        assert!(c.target_matches);
    }
    assert!(!check(&s, &t, &o, None, "ACCAGAACA").unwrap().consistent);
    let mut reverse = o.clone();
    reverse.reverse_donor = true;
    assert!(
        !check(&s, &t, &reverse, None, r.donor.sequence.as_str())
            .unwrap()
            .consistent
    );
}
#[test]
fn selected_guide_is_revalidated_and_nominal_distance_policy_is_explicit() {
    let s = mol(
        "AAAAAACGTCAGTACGATCGTACGATTGGTTACCGGATCCGATAGCTAGCATCGATCGGATC",
        false,
    );
    let t = mol(
        "AAAAAACGTCAGTACGATCGTACAATTGGTTACCGGATCCGATAGCTAGCATCGATCGGATC",
        false,
    );
    let go = crispr::GuideOptions::default();
    let report = crispr::design_guides(&s, &go).unwrap();
    let candidate = report
        .candidates
        .iter()
        .find(|c| c.guide_positions == (6..26).collect::<Vec<_>>() && !c.reverse)
        .unwrap()
        .clone();
    assert!(candidate.passes_filters);
    let guide = GuideBinding {
        options: go,
        candidate,
    };
    let mut o = opts(23, 1, "A", 4, 4);
    o.max_cut_distance = 0;
    let r = design(&s, &t, &o, Some(&guide)).unwrap();
    assert!(r.check.consistent);
    assert_eq!(r.check.guide_cut_distance, Some(0));
    let mut forged = guide.clone();
    forged.candidate.nominal_top_cut = Some(0);
    assert!(design(&s, &t, &o, Some(&forged)).is_err());
    let mut moved = o.clone();
    moved.edit_start = 40;
    let wrong_target = mol(
        "AAAAAACGTCAGTACGATCGTACGATTGGTTACCGGATCCGAATAGCTAGCATCGATCGGATC",
        false,
    );
    let moved = design(&s, &wrong_target, &moved, Some(&guide)).unwrap();
    assert_eq!(moved.check.guide_within_distance, Some(false));
    assert!(!moved.check.consistent);
}
#[test]
fn source_limits_stale_bindings_and_unresolved_inputs_are_errors_before_hashing() {
    let s = mol("AACCGGTTAACC", false);
    let t = mol("AACCAGAACC", false);
    let o = opts(4, 4, "AG", 3, 3);
    let r = design(&s, &t, &o, None).unwrap();
    let mut stale = s.clone();
    stale.name.push('!');
    assert!(r.validate_bindings_for(&stale, &t, &o, None).is_err());
    let mut huge = s.clone();
    huge.original_genbank = Some("x".repeat(8 * 1024 * 1024));
    assert!(matches!(
        batch_revision(&huge, &t, &o, None),
        Err(Error::TooLarge)
    ));
    let uncertain = mol("AACCGNTTAACC", false);
    assert!(design(&uncertain, &t, &o, None).is_err());
    let mut bad = o.clone();
    bad.edit_start = usize::MAX;
    assert!(design(&s, &t, &bad, None).is_err());
    let mut bad = o.clone();
    bad.left_arm = 0;
    assert!(design(&s, &t, &bad, None).is_err());
    let unchanged = mol("AACCGGTTAACC", false);
    let no_op = opts(4, 1, "G", 3, 3);
    let r = design(&s, &unchanged, &no_op, None).unwrap();
    assert!(!r.check.changes_sequence && !r.check.consistent);
}
#[test]
fn source_annotation_mapping_and_original_evidence_stay_separate_from_new_dna() {
    let mut s = mol("AACCGGTTAACC", false);
    s.original_genbank = Some("source evidence only".into());
    s.metadata
        .insert("source_quality".into(), "not copied".into());
    s.features.push(Feature {
        kind: "gene".into(),
        location: Location::range(1, 4),
        qualifiers: vec![("label".into(), Some("left_gene".into()))],
        needs_review: false,
    });
    let t = mol("AACCAGAACC", false);
    let mut o = opts(4, 4, "AG", 3, 3);
    let r = design(&s, &t, &o, None).unwrap();
    let f = r
        .donor
        .features
        .iter()
        .find(|f| f.label() == "left_gene")
        .unwrap();
    assert_eq!(f.location, Location::range(0, 3));
    assert!(r.donor.original_genbank.is_none() && !r.donor.metadata.contains_key("source_quality"));
    assert_eq!(s.original_genbank.as_deref(), Some("source evidence only"));
    o.reverse_donor = true;
    let r = design(&s, &t, &o, None).unwrap();
    let f = r
        .donor
        .features
        .iter()
        .find(|f| f.label() == "left_gene")
        .unwrap();
    assert_eq!(
        f.location,
        Location::Complement(Box::new(Location::range(5, 8)))
    );
}
#[test]
fn native_and_text_exports_preserve_only_declared_provenance_without_check_authority() {
    let s = mol("AACCGGTTAACC", false);
    let t = mol("AACCAGAACC", false);
    let o = opts(4, 4, "AG", 3, 3);
    let r = design(&s, &t, &o, None).unwrap();
    let native = protocol::load_document(&protocol::save_document(&r.donor).unwrap()).unwrap();
    assert_eq!(native, r.donor);
    for format in [formats::TextFormat::GenBank, formats::TextFormat::Embl] {
        let exported = formats::export_with_report(&r.donor, format).unwrap();
        assert!(exported.text.contains("edit_template:source_revision"));
        assert!(
            exported
                .warnings
                .iter()
                .any(|w| w.contains("unverified producer"))
        );
        let imported = formats::import_text(&exported.text)
            .unwrap()
            .molecules
            .remove(0);
        assert!(!imported.metadata.contains_key("edit_template:options"));
        assert!(
            imported
                .metadata
                .keys()
                .all(|key| !key.contains("consistent"))
        );
    }
    let mut oversized = r.donor.clone();
    oversized.metadata.insert(
        "edit_template:guide".into(),
        format!("{}\u{0007}", "x".repeat(20_000)),
    );
    let exported = formats::export_with_report(&oversized, formats::TextFormat::GenBank).unwrap();
    assert!(
        exported
            .text
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ")
            .contains("truncated; retain native")
    );
    assert!(!exported.text.contains('\u{0007}'));
}
#[test]
fn short_binary_sources_cover_every_edit_boundary_deletion_wrap_and_donor_orientation() {
    for n in 3..=5 {
        for bits in 0..1usize << n {
            let source = (0..n)
                .map(|i| if bits & (1 << i) == 0 { 'A' } else { 'C' })
                .collect::<String>();
            for circle in [false, true] {
                let s = mol(&source, circle);
                for start in 0..n {
                    for delete in 0..=n - 2 {
                        if !circle && (start == 0 || start + delete >= n) {
                            continue;
                        }
                        for insert in ["", "A", "C", "AC"] {
                            if delete == 0 && insert.is_empty() {
                                continue;
                            }
                            let mut sequence = String::new();
                            if circle && start + delete > n {
                                sequence.push_str(insert);
                                let first = (start + delete) % n;
                                sequence.extend(
                                    (0..n - delete)
                                        .map(|i| source.as_bytes()[(first + i) % n] as char),
                                );
                            } else {
                                for (i, b) in source.bytes().enumerate() {
                                    if i == start {
                                        sequence.push_str(insert);
                                    }
                                    if i < start || i >= start + delete {
                                        sequence.push(b as char);
                                    }
                                }
                            }
                            let t = mol(&sequence, circle);
                            let left = source.as_bytes()[(start + n - 1) % n] as char;
                            let right = source.as_bytes()[(start + delete) % n] as char;
                            let donor = format!("{left}{insert}{right}");
                            let unchanged = if !circle {
                                sequence == source
                            } else if sequence.len() != n {
                                false
                            } else {
                                (0..n).any(|origin| {
                                    format!("{}{}", &sequence[origin..], &sequence[..origin])
                                        == source
                                })
                            };
                            for reverse in [false, true] {
                                let mut o = opts(start, delete, insert, 1, 1);
                                o.reverse_donor = reverse;
                                let r = design(&s, &t, &o, None).unwrap();
                                assert_eq!(
                                    r.edited_reference.sequence.as_str(),
                                    sequence,
                                    "source{source} circle{circle} start{start} delete{delete} insert{insert}"
                                );
                                assert_eq!(
                                    r.donor.sequence.as_str(),
                                    if reverse { rc(&donor) } else { donor.clone() }
                                );
                                assert!(
                                    r.check.left_matches
                                        && r.check.right_matches
                                        && r.check.replacement_matches
                                        && r.check.target_matches
                                );
                                assert_eq!(r.check.changes_sequence, !unchanged);
                                assert_eq!(r.check.consistent, !unchanged);
                            }
                        }
                    }
                }
            }
        }
    }
}
