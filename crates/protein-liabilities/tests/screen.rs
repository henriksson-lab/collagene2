use collagene_model::{Alphabet, Molecule, Sequence, Topology};
use collagene_protein_liabilities::*;
fn source(s: &str) -> Molecule {
    Molecule {
        name: "Authored protein".into(),
        sequence: Sequence::new(Alphabet::Protein, s).unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    }
}
fn options() -> Options {
    Options {
        max_hits: MAX_HITS,
        provided_cdrs: None,
    }
}
fn cdr(m: &Molecule, start: usize, end: usize) -> ProvidedCdrs {
    ProvidedCdrs {
        source_revision: m.revision(),
        definition_label: "Explicit authored CDR context, no numbering claim".into(),
        provenance: "Coordinate fixture".into(),
        regions: vec![CdrRegion {
            domain_id: "variable-a".into(),
            kind: CdrKind::Cdr1,
            start,
            end,
        }],
    }
}
fn count(r: &Report, id: RuleId) -> (usize, usize) {
    let c = r.counts.iter().find(|c| c.rule == id).unwrap();
    (c.certain, c.possible)
}
#[test]
fn canonical_literature_motifs_and_explicit_exclusions_are_distinct_from_measured_ptms() {
    use RuleId::*;
    for pair in ["NG", "NS", "NN", "NT", "NH"] {
        let m = source(pair);
        let r = screen(&m, &options(), &mut || false).unwrap();
        assert_eq!(count(&r, CanonicalAsnDeamidation), (1, 0), "{pair}");
    }
    for pair in ["NA", "ND", "NF", "NP", "NQ"] {
        let m = source(pair);
        let r = screen(&m, &options(), &mut || false).unwrap();
        assert_eq!(count(&r, CanonicalAsnDeamidation), (0, 0), "{pair}");
    }
    for pair in ["DG", "DS", "DD", "DT", "DH"] {
        let m = source(pair);
        let r = screen(&m, &options(), &mut || false).unwrap();
        assert_eq!(count(&r, CanonicalAspIsomerization), (1, 0), "{pair}");
    }
    for pair in ["DN", "DA", "DR", "DY"] {
        let m = source(pair);
        let r = screen(&m, &options(), &mut || false).unwrap();
        assert_eq!(count(&r, CanonicalAspIsomerization), (0, 0), "{pair}");
    }
    for (word, n) in [
        ("NAS", 1),
        ("NAT", 1),
        ("NPS", 0),
        ("NPT", 0),
        ("NSS", 1),
        ("NTS", 1),
        ("NJS", 1),
        ("NUS", 1),
        ("NOS", 1),
        ("NA", 0),
        ("N*S", 0),
    ] {
        let r = screen(&source(word), &options(), &mut || false).unwrap();
        assert_eq!(count(&r, NGlycosylationSequon), (n, 0), "{word}");
    }
    let m = source("MWM*");
    let r = screen(&m, &options(), &mut || false).unwrap();
    assert_eq!(count(&r, MethionineOxidation), (2, 0));
    assert_eq!(count(&r, TryptophanOxidation), (1, 0));
    assert_eq!(r.stop_symbols, 1);
    assert!(r.cdr_rule_unassessed);
    assert_eq!(r.source_revision, m.revision());
    r.validate_for(&m, &options()).unwrap();
}
// Independent concrete-string oracle: no production bit masks, profiles or matcher.
fn concrete(s: &str) -> std::collections::BTreeSet<(usize, u8)> {
    let b = s.as_bytes();
    let mut out = std::collections::BTreeSet::new();
    for i in 0..b.len() {
        if i + 2 <= b.len() {
            let pair = &s[i..i + 2];
            if ["NG", "NS", "NN", "NT", "NH"].contains(&pair) {
                out.insert((i, 0));
            }
            if ["DG", "DS", "DD", "DT", "DH"].contains(&pair) {
                out.insert((i, 1));
            }
        }
        if i + 3 <= b.len() && b[i] == b'N' && b[i + 1] != b'P' && [b'S', b'T'].contains(&b[i + 2])
        {
            out.insert((i, 2));
        }
        if b[i] == b'M' {
            out.insert((i, 3));
        }
        if b[i] == b'W' {
            out.insert((i, 4));
        }
    }
    out
}
fn expansions(s: &str) -> Vec<String> {
    let mut rows = vec![String::new()];
    for c in s.chars() {
        let choices = match c {
            'B' => "DN".to_owned(),
            'Z' => "EQ".into(),
            'J' => "IL".into(),
            'X' => "ACDEFGHIKLMNPQRSTVWYOU".into(),
            c => c.to_string(),
        };
        rows = rows
            .into_iter()
            .flat_map(|prefix| {
                choices
                    .chars()
                    .map(move |c| format!("{prefix}{c}"))
                    .collect::<Vec<_>>()
            })
            .collect();
    }
    rows
}
#[test]
fn ambiguous_patterns_equal_union_and_intersection_of_concrete_residue_assignments() {
    for len in 1..=3 {
        for mut word in 0..6_usize.pow(len) {
            let mut s = String::new();
            for _ in 0..len {
                s.push(b"NPBSTX"[word % 6] as char);
                word /= 6;
            }
            let variants = expansions(&s);
            let mut union = std::collections::BTreeSet::new();
            let mut intersection = concrete(&variants[0]);
            for v in variants {
                let matched = concrete(&v);
                union.extend(matched.clone());
                intersection = intersection.intersection(&matched).copied().collect();
            }
            let r = screen(&source(&s), &options(), &mut || false).unwrap();
            let actual = r
                .hits
                .iter()
                .map(|h| {
                    let id = match h.rule {
                        RuleId::CanonicalAsnDeamidation => 0,
                        RuleId::CanonicalAspIsomerization => 1,
                        RuleId::NGlycosylationSequon => 2,
                        RuleId::MethionineOxidation => 3,
                        RuleId::TryptophanOxidation => 4,
                        _ => panic!("CDR rule without context"),
                    };
                    assert_eq!(
                        h.certainty == Certainty::CertainPattern,
                        intersection.contains(&(h.start, id)),
                        "{s}"
                    );
                    (h.start, id)
                })
                .collect();
            assert_eq!(union, actual, "{s}");
        }
    }
    let r = screen(&source("BZJXUO"), &options(), &mut || false).unwrap();
    assert_eq!(r.ambiguous_residues, 4);
}
#[test]
fn provided_cdrs_are_exact_source_context_not_automatic_numbering_or_global_absence() {
    let m = source("NAANH");
    let mut o = options();
    o.provided_cdrs = Some(cdr(&m, 0, 1));
    let r = screen(&m, &o, &mut || false).unwrap();
    assert_eq!(count(&r, RuleId::ProvidedCdrAsparagineReview), (1, 0));
    assert_eq!(count(&r, RuleId::CanonicalAsnDeamidation), (1, 0));
    assert!(!r.cdr_rule_unassessed);
    assert!(
        r.hits
            .iter()
            .any(|h| h.start == 3 && h.cdr == CdrMembership::OutsideProvided)
    );
    let m2 = source("NGAAN");
    let mut o = options();
    o.provided_cdrs = Some(cdr(&m2, 0, 1));
    let r = screen(&m2, &o, &mut || false).unwrap();
    let h = r
        .hits
        .iter()
        .find(|h| h.rule == RuleId::CanonicalAsnDeamidation)
        .unwrap();
    assert_eq!((h.start, h.end), (0, 2));
    assert_eq!(h.cdr, CdrMembership::Inside { region_index: 0 });
    let mut empty = o.clone();
    empty.provided_cdrs.as_mut().unwrap().regions.clear();
    assert!(
        screen(&m2, &empty, &mut || false)
            .unwrap()
            .cdr_rule_unassessed
    );
    let mut changed = m2.clone();
    changed.name = "Different source annotation".into();
    assert_eq!(
        screen(&changed, &o, &mut || false).unwrap_err(),
        Error::Binding
    );
}
#[test]
fn exact_counts_survive_retention_limits_overlapping_patterns_and_cancellation() {
    let m = source("NNSTMWNAT");
    let full = screen(&m, &options(), &mut || false).unwrap();
    for max_hits in [0, 1, 2, full.hits.len(), full.hits.len() + 1] {
        let o = Options {
            max_hits,
            provided_cdrs: None,
        };
        let r = screen(&m, &o, &mut || false).unwrap();
        assert_eq!(r.counts, full.counts);
        assert_eq!(r.hits, full.hits[..max_hits.min(full.hits.len())]);
        assert_eq!(r.hits_truncated, max_hits < full.hits.len());
        r.validate_for(&m, &o).unwrap();
    }
    assert_eq!(
        screen(&m, &options(), &mut || true).unwrap_err(),
        Error::Cancelled
    );
    let m = source(&"A".repeat(4096));
    let mut polls = 0;
    assert_eq!(
        screen(&m, &options(), &mut || {
            polls += 1;
            polls == 2
        })
        .unwrap_err(),
        Error::Cancelled
    );
}
#[test]
fn malformed_context_source_types_and_resource_limits_refuse_before_a_screen_is_claimed() {
    let m = source("N*AST");
    let mut o = options();
    o.provided_cdrs = Some(cdr(&m, 0, 3));
    assert!(screen(&m, &o, &mut || false).is_err());
    let m = source("NAAST");
    for (start, end) in [(0, 0), (2, 1), (0, 6), (5, 6)] {
        let mut o = options();
        o.provided_cdrs = Some(cdr(&m, start, end));
        assert!(screen(&m, &o, &mut || false).is_err());
    }
    let mut o = options();
    let mut c = cdr(&m, 0, 3);
    let mut other = c.regions[0].clone();
    other.kind = CdrKind::Cdr2;
    other.start = 2;
    other.end = 4;
    c.regions.push(other);
    o.provided_cdrs = Some(c);
    assert!(screen(&m, &o, &mut || false).is_err());
    let mut c = cdr(&m, 0, 1);
    let mut other = c.regions[0].clone();
    other.start = 2;
    other.end = 3;
    c.regions.push(other);
    o.provided_cdrs = Some(c);
    assert!(screen(&m, &o, &mut || false).is_err());
    let mut o = options();
    o.max_hits = MAX_HITS + 1;
    assert!(screen(&m, &o, &mut || false).is_err());
    assert!(source_revision(&Molecule::dna("DNA", "AC", Topology::Linear).unwrap()).is_err());
    assert!(source_revision(&source(&"A".repeat(MAX_RESIDUES + 1))).is_err());
}
#[test]
fn reports_are_reconstructed_against_source_rules_options_and_complete_counts() {
    let m = source("NGNATMW");
    let o = options();
    let r = screen(&m, &o, &mut || false).unwrap();
    let roundtrip: Report = serde_json::from_slice(&serde_json::to_vec(&r).unwrap()).unwrap();
    roundtrip.validate_for(&m, &o).unwrap();
    assert_eq!(r.context_revision, request_revision(&m, &o).unwrap());
    let mut changed = r.clone();
    changed.hits[0].certainty = Certainty::PossiblePattern;
    assert!(changed.validate_for(&m, &o).is_err());
    let mut changed = r.clone();
    changed.counts[0].certain += 1;
    assert!(changed.validate_for(&m, &o).is_err());
    let mut changed = r.clone();
    changed.catalog_revision = "f".repeat(64);
    assert!(changed.validate_for(&m, &o).is_err());
    let mut changed = o.clone();
    changed.max_hits = 0;
    assert!(r.validate_for(&m, &changed).is_err());
    assert!(r.validate_for(&source("NGNATMWA"), &o).is_err());
}
