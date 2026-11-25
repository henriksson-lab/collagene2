use collagene_model::{Molecule, Sequence, Topology};
use collagene_recombination::multisite::{self, Profile, Request};
use collagene_recombination::{Family, Input, Limits, Reaction, SiteVariant, SourceRecord};

fn add(
    profile: &mut Profile,
    id: &str,
    family: Family,
    subtype: u8,
    sequence: String,
    overlap_start: usize,
) {
    profile.sources.push(SourceRecord {
        id: id.into(),
        url: "https://example.org/synthetic-test-record".into(),
        raw_sha256: "1".repeat(64),
        vector_name: "independently declared synthetic fixture, not a manufacturer profile".into(),
        vector_length: sequence.len(),
        start_1based: 1,
        end_1based: sequence.len(),
        reverse_extraction: false,
        map_url: "https://example.org/synthetic-test-map".into(),
        map_sha256: "2".repeat(64),
    });
    profile.variants.push(SiteVariant {
        id: id.into(),
        family,
        subtype,
        sequence,
        overlap_start,
        source_id: id.into(),
    });
}
fn empty_profile() -> Profile {
    Profile {
        id: "synthetic-four-specificities".into(),
        edition: "literal-fixture-v1".into(),
        sources: vec![],
        variants: vec![],
        limitations: vec![
            "Synthetic reference bookkeeping fixture; no biochemical applicability claim".into(),
        ],
    }
}
fn core(s: u8) -> &'static str {
    match s {
        1 => "GTACAAA",
        2 => "GTACAAG",
        3 => "GTACAAC",
        4 => "GTACAAT",
        _ => panic!("fixture specificity"),
    }
}
fn site(f: Family, s: u8) -> String {
    match f {
        Family::L => format!("ACGG{}TCC", core(s)),
        Family::R => format!("TGGA{}GCT", core(s)),
        _ => panic!("fixture family"),
    }
}
fn fixture() -> (Vec<Input>, Request) {
    let mut profile = empty_profile();
    for s in 1..=4 {
        for (f, prefix) in [(Family::L, "L"), (Family::R, "R")] {
            add(&mut profile, &format!("{prefix}{s}"), f, s, site(f, s), 4);
        }
    }
    let sequences = [
        format!("{}CCCCA{}AAAAC", site(Family::L, 4), site(Family::R, 1)),
        format!("{}GGGGT{}TTTTC", site(Family::L, 1), site(Family::L, 2)),
        format!("{}AAAAG{}CCCCT", site(Family::R, 2), site(Family::L, 3)),
        format!("{}TTTTA{}GGGGA", site(Family::R, 4), site(Family::R, 3)),
    ];
    let inputs = sequences
        .iter()
        .enumerate()
        .map(|(i, s)| Input {
            role: format!("parent-{i}"),
            molecule: Molecule::dna(format!("parent {i}"), s, Topology::Circular).unwrap(),
        })
        .collect();
    (
        inputs,
        Request {
            profile,
            reaction: Reaction::LR,
            specificities: vec![4, 1, 2, 3],
            limits: Limits::default(),
        },
    )
}
const MAIN: &str = concat!(
    "GTACAATTCCCCCCATGGA",
    "GTACAAATCCGGGGTACGG",
    "GTACAAGGCTAAAAGACGG",
    "GTACAACGCTGGGGATGGA"
);
const BACKBONES: &str = concat!(
    "GTACAAAGCTAAAACACGG",
    "GTACAATGCTTTTTATGGA",
    "GTACAACTCCCCCCTTGGA",
    "GTACAAGTCCTTTTCACGG"
);
fn rc(s: &str) -> String {
    Sequence::dna(s)
        .unwrap()
        .reverse_complement()
        .unwrap()
        .as_str()
        .into()
}
fn circular_eq(actual: &str, expected: &str) -> bool {
    actual.len() == expected.len()
        && (format!("{actual}{actual}").contains(expected)
            || format!("{}{}", rc(actual), rc(actual)).contains(expected))
}
fn conservation(inputs: &[Input], candidate: &multisite::Candidate) {
    for (i, input) in inputs.iter().enumerate() {
        let mut counts = vec![0u8; input.molecule.sequence.len()];
        for m in candidate
            .products
            .iter()
            .flat_map(|p| &p.segments)
            .filter(|m| m.input == i)
        {
            for count in &mut counts[m.source_start..m.source_end] {
                *count += 1;
            }
        }
        assert!(
            counts.iter().all(|count| *count == 1),
            "each original reference coordinate, including byproducts, occurs exactly once"
        );
    }
}
#[test]
fn four_parent_exchange_has_literal_main_and_all_reciprocal_reference_bases() {
    let (inputs, request) = fixture();
    let report = multisite::recombine(&inputs, &request).unwrap();
    assert!(report.complete);
    assert_eq!(report.candidates.len(), 1);
    let c = &report.candidates[0];
    assert_eq!(c.sites.len(), 8);
    assert_eq!(c.pairs, vec![[0, 1], [2, 3], [4, 5], [6, 7]]);
    assert_eq!(c.products.len(), 2);
    assert_eq!(
        c.products
            .iter()
            .map(|p| p.molecule.sequence.len())
            .sum::<usize>(),
        152
    );
    assert_eq!(MAIN.len(), 76);
    assert!(
        c.products
            .iter()
            .any(|p| circular_eq(p.molecule.sequence.as_str(), MAIN))
    );
    let backbone = BACKBONES.replace(' ', "");
    assert_eq!(backbone.len(), 76);
    assert!(
        c.products
            .iter()
            .any(|p| circular_eq(p.molecule.sequence.as_str(), &backbone))
    );
    assert!(
        c.products
            .iter()
            .all(|p| p.molecule.topology == Topology::Circular)
    );
    assert_eq!(
        c.products.iter().map(|p| p.junctions.len()).sum::<usize>(),
        8
    );
    conservation(&inputs, c);
    multisite::audit_candidate(&inputs, &request.profile, c).unwrap();
    multisite::replay(&inputs, &request, &report).unwrap();
    assert!(report.validate_binding(&"0".repeat(64)).is_err());
    let mut edited = inputs.clone();
    edited[0].molecule.name = "changed exact source snapshot".into();
    assert_ne!(
        multisite::revision(&edited, &request).unwrap(),
        report.batch_revision
    );
    assert!(multisite::replay(&edited, &request, &report).is_err());
}
#[test]
fn reverse_parents_origin_spanning_sites_and_unreacted_fifth_parent_are_preserved() {
    let (mut inputs, request) = fixture();
    let s = inputs[0].molecule.sequence.as_str();
    let rotated = format!("{}{}", &s[7..], &s[..7]);
    inputs[0].molecule = Molecule::dna("origin shifted", &rotated, Topology::Circular).unwrap();
    inputs[1].molecule = Molecule::dna(
        "reverse entry",
        &rc(inputs[1].molecule.sequence.as_str()),
        Topology::Circular,
    )
    .unwrap();
    inputs.push(Input {
        role: "unreacted".into(),
        molecule: Molecule::dna(
            "unreacted circular parent",
            "CCCGGGTTTAAACC",
            Topology::Circular,
        )
        .unwrap(),
    });
    let report = multisite::recombine(&inputs, &request).unwrap();
    assert!(report.complete);
    assert_eq!(report.candidates.len(), 1);
    let c = &report.candidates[0];
    assert_eq!(c.products.len(), 3);
    assert!(c.sites.iter().any(|h| h.reverse));
    assert!(
        c.sites
            .iter()
            .any(|h| h.start + h.length > inputs[h.input].molecule.sequence.len())
    );
    assert!(
        c.products
            .iter()
            .any(|p| circular_eq(p.molecule.sequence.as_str(), MAIN))
    );
    assert!(
        c.products
            .iter()
            .any(|p| p.molecule.sequence.as_str() == "CCCGGGTTTAAACC")
    );
    conservation(&inputs, c);
}
fn nested_fixture(parent: &str) -> (Vec<Input>, Request) {
    let mut profile = empty_profile();
    let b = format!("A{}C", core(3));
    let l = format!("GG{b}TT");
    let p = format!("CC{l}AA");
    add(&mut profile, "B3", Family::B, 3, b, 1);
    add(&mut profile, "L3", Family::L, 3, l.clone(), 3);
    add(&mut profile, "P3", Family::P, 3, p.clone(), 5);
    add(
        &mut profile,
        "R3",
        Family::R,
        3,
        format!("TTT{}GGG", core(3)),
        3,
    );
    let first = if parent == "P" { p } else { l };
    let inputs = vec![
        Input {
            role: "nested".into(),
            molecule: Molecule::dna("nested context", &first, Topology::Linear).unwrap(),
        },
        Input {
            role: "R".into(),
            molecule: Molecule::dna("cognate R", &format!("TTT{}GGG", core(3)), Topology::Linear)
                .unwrap(),
        },
    ];
    (
        inputs,
        Request {
            profile,
            reaction: Reaction::LR,
            specificities: vec![3],
            limits: Limits::default(),
        },
    )
}
#[test]
fn all_family_maximal_context_classification_never_turns_p_donor_into_l_entry() {
    let (inputs, request) = nested_fixture("P");
    let report = multisite::recombine(&inputs, &request).unwrap();
    assert!(report.complete);
    for f in [Family::B, Family::L, Family::P] {
        assert!(report.sites.iter().any(|h| h.input == 0 && h.family == f));
    }
    assert!(report.candidates.is_empty());
    let (inputs, mut request) = nested_fixture("L");
    let mut alias = request
        .profile
        .variants
        .iter()
        .find(|v| v.id == "L3")
        .unwrap()
        .clone();
    alias.id = "Z-L3-identical-context".into();
    request.profile.variants.push(alias);
    let report = multisite::recombine(&inputs, &request).unwrap();
    assert!(report.complete);
    assert_eq!(report.candidates.len(), 1);
    assert_eq!(report.candidates[0].sites[0].variant_id, "L3");
    assert!(
        report
            .warnings
            .iter()
            .any(|s| s.contains("enumerated once"))
    );
    conservation(&inputs, &report.candidates[0]);
}
#[test]
fn equal_maximal_conflicting_families_and_ambiguous_sources_are_not_fresh_bindings() {
    let (inputs, mut request) = nested_fixture("L");
    let mut conflict = request
        .profile
        .variants
        .iter()
        .find(|v| v.id == "L3")
        .unwrap()
        .clone();
    conflict.id = "P3-equal-context".into();
    conflict.family = Family::P;
    request.profile.variants.push(conflict);
    let report = multisite::recombine(&inputs, &request).unwrap();
    assert!(!report.complete);
    assert!(report.candidates.is_empty());
    assert!(!report.unsupported.is_empty());
    let (mut inputs, request) = nested_fixture("L");
    let mut bytes = inputs[0].molecule.sequence.as_str().as_bytes().to_vec();
    bytes[0] = b'N';
    inputs[0].molecule.sequence = Sequence::dna(std::str::from_utf8(&bytes).unwrap()).unwrap();
    let report = multisite::recombine(&inputs, &request).unwrap();
    assert!(!report.complete);
    assert!(report.candidates.is_empty());
}
#[test]
fn caps_stop_global_scanning_or_whole_candidate_before_partial_reference_outputs() {
    let (inputs, request) = fixture();
    for mode in 0..5 {
        let mut r = request.clone();
        match mode {
            0 => r.limits.max_comparisons = 1,
            1 => r.limits.max_sites = 1,
            2 => r.limits.max_combinations = 1,
            3 => r.limits.max_components_per_candidate = 1,
            _ => r.limits.max_product_bases = 151,
        };
        let report = multisite::recombine(&inputs, &r).unwrap();
        assert!(!report.complete, "mode {mode}");
        assert!(report.candidates.is_empty(), "mode {mode}");
        assert!(!report.limits.is_empty());
    }
    let mut r = request.clone();
    r.limits.max_report_bytes = 1024;
    assert!(multisite::recombine(&inputs, &r).is_err());
    r.limits.max_comparisons = 1;
    assert!(
        multisite::recombine(&inputs, &r).is_err(),
        "incomplete scan still obeys the declared report-byte cap"
    );
    r.limits.max_comparisons = request.limits.max_comparisons;
    r.reaction = Reaction::BP;
    assert!(
        multisite::recombine(&inputs, &r).is_err(),
        "missing cognate groups still obey the declared report-byte cap"
    );
}
#[test]
fn noncognate_profiles_duplicate_specificities_and_overlapping_requested_sites_are_rejected_or_withheld()
 {
    let (inputs, mut request) = fixture();
    request.specificities.push(4);
    assert!(multisite::preflight(&inputs, &request).is_err());
    request.specificities.pop();
    request.profile.variants[1]
        .sequence
        .replace_range(4..11, "AAAAAAA");
    assert!(request.profile.validate().is_err());
    // Same source coordinates cannot serve two different explicitly requested pairings.
    let (inputs, mut request) = nested_fixture("L");
    for original in request.profile.variants.clone() {
        let mut v = original;
        v.id = format!("{}-specificity4", v.id);
        v.subtype = 4;
        request.profile.variants.push(v);
    }
    request.specificities = vec![3, 4];
    let report = multisite::recombine(&inputs, &request).unwrap();
    assert!(report.complete);
    assert!(report.candidates.is_empty());
}

#[test]
fn bp_specificity_sixteen_swaps_both_reference_tails_without_an_lr_assumption() {
    let mut profile = empty_profile();
    add(
        &mut profile,
        "B16",
        Family::B,
        16,
        "ACGACGTCGATTA".into(),
        3,
    );
    add(
        &mut profile,
        "P16",
        Family::P,
        16,
        "GGTACGTCGACCA".into(),
        3,
    );
    let inputs = vec![
        Input {
            role: "B".into(),
            molecule: Molecule::dna("B parent", "ACGACGTCGATTA", Topology::Linear).unwrap(),
        },
        Input {
            role: "P".into(),
            molecule: Molecule::dna("P parent", "GGTACGTCGACCA", Topology::Linear).unwrap(),
        },
    ];
    let request = Request {
        profile,
        reaction: Reaction::BP,
        specificities: vec![16],
        limits: Limits::default(),
    };
    let report = multisite::recombine(&inputs, &request).unwrap();
    assert!(report.complete);
    assert_eq!(report.candidates.len(), 1);
    let c = &report.candidates[0];
    let sequences = c
        .products
        .iter()
        .map(|p| p.molecule.sequence.as_str())
        .collect::<Vec<_>>();
    assert!(sequences.contains(&"ACGACGTCGACCA"));
    // The first remaining unmatched port is B's right end, so this product
    // traverses B[3..13] backwards and then P[0..3] backwards. Its directed
    // reference is the exact opposite strand of GGTACGTCGATTA.
    assert!(sequences.contains(&"TAATCGACGTACC"));
    let reciprocal = c
        .products
        .iter()
        .find(|p| p.molecule.sequence.as_str() == "TAATCGACGTACC")
        .unwrap();
    assert_eq!(
        reciprocal
            .segments
            .iter()
            .map(|m| (
                m.input,
                m.source_start,
                m.source_end,
                m.reverse,
                m.output_start,
                m.output_end,
            ))
            .collect::<Vec<_>>(),
        vec![(0, 3, 13, true, 0, 10), (1, 0, 3, true, 10, 13),]
    );
    assert!(
        c.products
            .iter()
            .all(|p| p.molecule.topology == Topology::Linear)
    );
    let families = c
        .products
        .iter()
        .flat_map(|p| &p.junctions)
        .map(|j| j.family)
        .collect::<Vec<_>>();
    assert!(families.contains(&Family::L));
    assert!(families.contains(&Family::R));
    conservation(&inputs, c);
}
#[test]
fn candidate_cap_keeps_a_whole_alternative_and_never_claims_a_unique_complete_search() {
    let (mut inputs, mut request) = nested_fixture("L");
    let original = inputs[0].molecule.sequence.as_str();
    inputs[0].molecule = Molecule::dna(
        "two separated L contexts",
        &format!("{original}CCCC{original}"),
        Topology::Linear,
    )
    .unwrap();
    request.limits.max_candidates = 1;
    let report = multisite::recombine(&inputs, &request).unwrap();
    assert!(!report.complete);
    assert_eq!(report.candidates.len(), 1);
    assert!(report.limits.iter().any(|s| s.contains("Candidate limit")));
    conservation(&inputs, &report.candidates[0]);
}
#[test]
fn public_map_audit_rejects_unbounded_map_shapes_before_reconstructing_dna() {
    let (inputs, request) = fixture();
    let report = multisite::recombine(&inputs, &request).unwrap();
    let mut c = report.candidates[0].clone();
    let map = c.products[0].segments[0].clone();
    c.products[0].segments = vec![map; 75];
    assert!(multisite::audit_candidate(&inputs, &request.profile, &c).is_err());
}
