//! Literal source/CDS/site fixtures, independent of search-produced expectations.
use collagene_core::{
    cloning::{
        golden_gate::domestication::{self as d, Limits, Patch, Request, Selection, SourceRange},
        restriction::catalog as catalog_digest,
    },
    sequence::codon_design::DesignOptions,
};
use collagene_model::{Feature, Location, Molecule, Sequence, Topology};
const CODING: &str = "ATGGGTCTCGAATAA";
fn range(start: usize, end: usize) -> Location {
    Location::Range {
        start,
        end,
        before: false,
        after: false,
    }
}
fn feature(location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![
            ("label".into(), Some("literal gene".into())),
            (
                "note".into(),
                Some("Preserve this original qualifier exactly".into()),
            ),
        ],
        needs_review: false,
    }
}
fn source(text: &str, topology: Topology, location: Option<Location>) -> Molecule {
    let mut s = Molecule::dna("literal original", text, topology).unwrap();
    if let Some(location) = location {
        s.features.push(feature(location));
    }
    s
}
fn request(s: &Molecule, enzyme: &str, selected: &[usize]) -> Request {
    Request {
        source_revision: s.revision(),
        selection: catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &[enzyme.into()],
        )
        .unwrap(),
        scope: vec![SourceRange {
            start: 0,
            end: s.sequence.len(),
        }],
        retained_sites: vec![],
        coding_features: selected
            .iter()
            .map(|i| d::feature_pin(s, *i).unwrap())
            .collect(),
        protected: vec![],
        options: DesignOptions::default(),
        limits: Limits::default(),
    }
}
fn selection(s: &Molecule, r: &Request, patches: &[(usize, &str, &str)]) -> Selection {
    let patches = patches
        .iter()
        .map(|(position, a, b)| Patch {
            position: *position,
            expected: (*a).into(),
            replacement: (*b).into(),
        })
        .collect::<Vec<_>>();
    Selection {
        context_revision: d::propose_revision(s, r).unwrap(),
        candidate_id: d::candidate_id(r, &patches).unwrap(),
        patches,
    }
}
#[test]
fn literal_synonymous_site_removal_preserves_complete_source_annotation_and_private_target() {
    let s = source(CODING, Topology::Linear, Some(range(0, 15)));
    let r = request(&s, "bsai", &[0]);
    let before = serde_json::to_vec(&s).unwrap();
    let inspected = d::inspect(&s, &r).unwrap();
    assert_eq!(inspected.sites.len(), 1);
    assert_eq!(inspected.sites[0].pin.start, 3);
    assert_eq!(inspected.features[0].protein.as_deref(), Some("MGLE"));
    let picked = selection(&s, &r, &[(5, "T", "C")]);
    let proof = d::check(&s, &r, &picked).unwrap();
    assert_eq!(proof.target().sequence.as_str(), "ATGGGCCTCGAATAA");
    assert!(proof.report().sites_after.is_empty());
    assert_eq!(proof.report().features[0].before_protein, "MGLE");
    assert_eq!(proof.report().features[0].after_protein, "MGLE");
    assert!(proof.report().features[0].has_terminal_stop);
    assert_eq!(proof.target().features[0].location, s.features[0].location);
    assert_eq!(
        proof.target().features[0].qualifiers,
        s.features[0].qualifiers
    );
    assert!(proof.target().features[0].needs_review);
    assert_eq!(proof.report().target_revision, proof.target().revision());
    proof
        .report()
        .validate_binding(&d::check_revision(&s, &r, &picked).unwrap())
        .unwrap();
    assert_eq!(serde_json::to_vec(&s).unwrap(), before);
    let proposal = d::propose(&s, &r).unwrap();
    proposal
        .validate_binding(&d::propose_revision(&s, &r).unwrap())
        .unwrap();
    assert!(
        proposal
            .searches
            .iter()
            .all(|search| search.work_used <= search.reserved_work)
    );
    assert!(
        proposal
            .searches
            .iter()
            .map(|search| search.work_used)
            .sum::<u64>()
            <= r.limits.max_work
    );
    let candidate = proposal
        .candidate
        .as_ref()
        .expect("one-base synonymous alternative exists");
    assert_eq!(candidate.patches.len(), 1);
    assert!(
        !proposal.search_complete,
        "independent per-feature search is not a joint global optimum proof"
    );
    assert!(
        d::check(&s, &r, candidate)
            .unwrap()
            .report()
            .sites_after
            .is_empty()
    );
}
#[test]
fn reverse_compound_and_origin_spanning_cds_keep_their_true_source_coordinates() {
    let reverse = Sequence::dna(CODING).unwrap().reverse_complement().unwrap();
    let s = source(
        reverse.as_str(),
        Topology::Linear,
        Some(Location::Complement(Box::new(range(0, 15)))),
    );
    let r = request(&s, "bsai", &[0]);
    let proof = d::check(&s, &r, &selection(&s, &r, &[(9, "A", "G")])).unwrap();
    assert_eq!(proof.target().sequence.as_str(), "TTATTCGAGGCCCAT");
    assert_eq!(proof.report().features[0].after_protein, "MGLE");
    let joined = Location::Join(vec![range(0, 12), range(15, 18)]);
    let s = source("ATGGGTCTCGAACCCTAA", Topology::Linear, Some(joined.clone()));
    let r = request(&s, "bsai", &[0]);
    let proof = d::check(&s, &r, &selection(&s, &r, &[(5, "T", "C")])).unwrap();
    assert_eq!(proof.target().sequence.as_str(), "ATGGGCCTCGAACCCTAA");
    assert_eq!(proof.target().features[0].location, joined);
    assert_eq!(proof.report().features[0].after_protein, "MGLE");
    let joined = Location::Join(vec![range(9, 15), range(0, 9)]);
    let s = source("CTCGAATAAATGGGT", Topology::Circular, Some(joined.clone()));
    let mut r = request(&s, "bsai", &[0]);
    r.scope = vec![SourceRange { start: 12, end: 18 }];
    let inspected = d::inspect(&s, &r).unwrap();
    assert_eq!(inspected.sites[0].pin.start, 12);
    let proof = d::check(&s, &r, &selection(&s, &r, &[(14, "T", "C")])).unwrap();
    assert_eq!(proof.target().sequence.as_str(), "CTCGAATAAATGGGC");
    assert_eq!(proof.target().features[0].location, joined);
    assert_eq!(proof.report().features[0].after_protein, "MGLE");
}
#[test]
fn phase_prefixes_alternative_initiators_and_terminal_stops_are_independently_preserved() {
    for (prefix, phase) in [("A", "2"), ("CC", "3")] {
        let mut s = source(
            &format!("{prefix}{CODING}"),
            Topology::Linear,
            Some(range(0, prefix.len() + 15)),
        );
        s.features[0]
            .qualifiers
            .push(("codon_start".into(), Some(phase.into())));
        let r = request(&s, "bsai", &[0]);
        let proof = d::check(&s, &r, &selection(&s, &r, &[(prefix.len() + 5, "T", "C")])).unwrap();
        assert_eq!(proof.report().features[0].after_protein, "MGLE");
        let mut edits = vec![
            (0, &s.sequence.as_str()[..1], "G"),
            (prefix.len() + 5, "T", "C"),
        ];
        if s.sequence.as_str().starts_with('G') {
            edits[0].2 = "A";
        }
        assert!(
            d::check(&s, &r, &selection(&s, &r, &edits)).is_err(),
            "incomplete leading codon remains literal"
        );
    }
    let mut s = source("CTGGGTCTCGAATAA", Topology::Linear, Some(range(0, 15)));
    s.features[0]
        .qualifiers
        .push(("transl_table".into(), Some("11".into())));
    let r = request(&s, "bsai", &[0]);
    assert_eq!(
        d::check(&s, &r, &selection(&s, &r, &[(5, "T", "C")]))
            .unwrap()
            .report()
            .features[0]
            .after_protein,
        "MGLE"
    );
    assert!(
        d::check(&s, &r, &selection(&s, &r, &[(2, "G", "A"), (5, "T", "C")])).is_err(),
        "ordinary leucine synonym cannot erase recognized initiation"
    );
    assert!(
        d::check(&s, &r, &selection(&s, &r, &[(5, "T", "C"), (12, "T", "C")])).is_err(),
        "terminal stop cannot become glutamine"
    );
}
#[test]
fn overlapping_or_unsupported_cds_are_never_silently_dropped() {
    let mut s = source(CODING, Topology::Linear, Some(range(0, 15)));
    s.features.push(feature(range(4, 13)));
    let r = request(&s, "bsai", &[0]);
    assert!(
        d::check(&s, &r, &selection(&s, &r, &[(5, "T", "C")])).is_err(),
        "GGT→GGC is silent in selected frame but GTC→GCC changes overlapping frame"
    );
    assert!(
        d::propose(&s, &r).unwrap().candidate.is_none(),
        "overlap is conservatively locked during independent search"
    );
    let mut s = source(CODING, Topology::Linear, Some(range(0, 15)));
    s.features[0]
        .qualifiers
        .push(("exception".into(), Some("custom biology".into())));
    let r = request(&s, "bsai", &[0]);
    assert!(!d::inspect(&s, &r).unwrap().features[0].supported);
    assert!(d::check(&s, &r, &selection(&s, &r, &[(5, "T", "C")])).is_err());
    for location in [
        Location::Range {
            start: 10,
            end: 15,
            before: true,
            after: false,
        },
        Location::Join(vec![
            range(10, 15),
            Location::Unresolved("remote:unknown".into()),
        ]),
    ] {
        let mut s = source(CODING, Topology::Linear, Some(range(0, 15)));
        s.features.push(feature(location));
        let r = request(&s, "bsai", &[0]);
        assert!(
            d::check(&s, &r, &selection(&s, &r, &[(5, "T", "C")])).is_err(),
            "unknown coverage can include a patch outside nominal known spans"
        );
        assert!(d::propose(&s, &r).unwrap().candidate.is_none());
    }
}
#[test]
fn full_source_rescan_preserves_wanted_sites_and_blocks_new_reverse_or_origin_sites() {
    let s = source("GGTCTCAACCGGTCTC", Topology::Linear, None);
    let mut r = request(&s, "bsai", &[]);
    r.scope = vec![SourceRange { start: 0, end: 6 }];
    let sites = d::inspect(&s, &r).unwrap().sites;
    r.retained_sites = vec![
        sites
            .iter()
            .find(|v| v.pin.start == 10)
            .unwrap()
            .pin
            .clone(),
    ];
    let proof = d::check(&s, &r, &selection(&s, &r, &[(2, "T", "A")])).unwrap();
    assert_eq!(proof.report().sites_after.len(), 1);
    assert_eq!(proof.report().sites_after[0].pin.start, 10);
    assert!(d::check(&s, &r, &selection(&s, &r, &[(2, "T", "A"), (12, "T", "A")])).is_err());
    r.protected = vec![SourceRange { start: 2, end: 3 }];
    assert!(d::check(&s, &r, &selection(&s, &r, &[(2, "T", "A")])).is_err());
    let s = source("GGTCTCAAAAGAGACT", Topology::Linear, None);
    let r = request(&s, "bsai", &[]);
    assert!(
        d::check(&s, &r, &selection(&s, &r, &[(2, "T", "A"), (15, "T", "C")])).is_err(),
        "new reverse GAGACC is found away from old site"
    );
    let s = source("CTCAAAGGTCTCGGA", Topology::Circular, None);
    let r = request(&s, "bsai", &[]);
    assert!(
        d::check(&s, &r, &selection(&s, &r, &[(8, "T", "A"), (14, "A", "T")])).is_err(),
        "new BsaI motif spans the source origin"
    );
}
#[test]
fn complete_codon_permission_can_require_an_extra_change_outside_the_unwanted_motif() {
    let s = source("ATGAAAGCTTCTGAATAA", Topology::Linear, Some(range(0, 18)));
    let mut r = request(&s, "hindiii", &[0]);
    // AAGCTT spans4..10; only the first base of Ser TCT lies inside it.
    // Earlier motif bases are protected. Removing this site requires Ser AGT,
    // changing T9→A AND C10→G, the latter outside the recognition span.
    r.scope = vec![SourceRange { start: 4, end: 10 }];
    r.protected = vec![SourceRange { start: 0, end: 9 }];
    let proposed = d::propose(&s, &r).unwrap();
    let candidate = proposed
        .candidate
        .as_ref()
        .expect("complete-codon search must admit AGT");
    assert_eq!(
        candidate.patches,
        vec![
            Patch {
                position: 9,
                expected: "T".into(),
                replacement: "A".into()
            },
            Patch {
                position: 10,
                expected: "C".into(),
                replacement: "G".into()
            }
        ]
    );
    let proof = d::check(&s, &r, candidate).unwrap();
    assert_eq!(proof.target().sequence.as_str(), "ATGAAAGCTAGTGAATAA");
    assert_eq!(proof.report().features[0].after_protein, "MKASE");
}
#[test]
fn distant_substitutions_preserve_intervening_annotations_and_never_create_a_fake_spanning_edit() {
    let mut s = source(
        &format!("GGTCTC{}GGTCTC", "A".repeat(120_000)),
        Topology::Linear,
        None,
    );
    s.features.push(Feature {
        kind: "misc_feature".into(),
        location: range(50, 100),
        qualifiers: vec![("label".into(), Some("unchanged between edits".into()))],
        needs_review: false,
    });
    let r = request(&s, "bsai", &[]);
    let proof = d::check(
        &s,
        &r,
        &selection(&s, &r, &[(2, "T", "A"), (120_008, "T", "A")]),
    )
    .unwrap();
    assert_eq!(proof.target().features[0], s.features[0]);
    assert_eq!(proof.target().sequence.len(), s.sequence.len());
    assert!(proof.report().sites_after.is_empty());
}
#[test]
fn stale_source_feature_catalog_candidate_and_limits_fail_before_any_target_is_accepted() {
    let s = source(CODING, Topology::Linear, Some(range(0, 15)));
    let r = request(&s, "bsai", &[0]);
    let picked = selection(&s, &r, &[(5, "T", "C")]);
    let mut wrong = picked.clone();
    wrong.candidate_id = "0".repeat(64);
    assert!(d::check(&s, &r, &wrong).is_err());
    let mut wrong = picked.clone();
    wrong.patches[0].expected = "G".into();
    wrong.candidate_id = d::candidate_id(&r, &wrong.patches).unwrap();
    assert!(d::check(&s, &r, &wrong).is_err());
    let mut changed = s.clone();
    changed.features[0]
        .qualifiers
        .push(("note".into(), Some("changed feature".into())));
    let mut stale = r.clone();
    stale.source_revision = changed.revision();
    assert!(d::preflight(&changed, &stale).is_err());
    let mut stale = r.clone();
    stale.selection.catalog_content_hash = "0".repeat(64);
    assert!(d::preflight(&s, &stale).is_err());
    let mut bounded = r.clone();
    bounded.limits.max_work = 1;
    assert!(d::preflight(&s, &bounded).is_err());
    let proof = d::check(&s, &r, &picked).unwrap();
    let mut historical = proof.report().clone();
    historical.target_revision = "0".repeat(64);
    // Shape-only imported receipts are never a constructor of the private proof.
    historical
        .validate_binding(proof.context_revision())
        .unwrap();
    assert_ne!(
        historical.target_revision,
        d::check(&s, &r, &picked).unwrap().target().revision()
    );
    let mut malformed = proof.report().clone();
    malformed.candidate_id = "0".repeat(64);
    assert!(
        malformed
            .validate_binding(proof.context_revision())
            .is_err()
    );
}

#[test]
fn successful_pass_refunds_its_unused_ceiling_so_two_cds_sites_share_one_real_budget() {
    let mut s = source(
        &format!("{CODING}AACCAACCAA{CODING}"),
        Topology::Linear,
        Some(range(0, 15)),
    );
    s.features.push(feature(range(25, 40)));
    let mut r = request(&s, "bsai", &[0, 1]);
    r.limits.max_work = 2_000_000;
    r.options.max_work = 2_000_000;
    let before = d::inspect(&s, &r).unwrap();
    assert_eq!(before.sites.len(), 2);
    let proposed = d::propose(&s, &r).unwrap();
    assert_eq!(
        proposed.searches.len(),
        2,
        "both selected source CDSs must actually be visited"
    );
    assert!(proposed.searches.iter().all(|s| s.search_complete));
    assert!(
        proposed
            .searches
            .iter()
            .map(|s| s.reserved_work)
            .sum::<u64>()
            > r.limits.max_work,
        "sequential pass ceilings may exceed aggregate budget when earlier unused work is refunded"
    );
    assert!(proposed.searches.iter().map(|s| s.work_used).sum::<u64>() < r.limits.max_work);
    let candidate = proposed
        .candidate
        .as_ref()
        .expect("both independent sites have a bounded synonymous solution");
    let checked = d::check(&s, &r, candidate).unwrap();
    assert!(checked.report().sites_after.is_empty());
    assert_eq!(checked.report().features.len(), 2);
    assert!(
        checked
            .report()
            .features
            .iter()
            .all(|f| f.before_protein == "MGLE"
                && f.after_protein == "MGLE"
                && f.has_terminal_stop)
    );
    assert_eq!(
        checked.target().features[0].location,
        s.features[0].location
    );
    assert_eq!(
        checked.target().features[1].location,
        s.features[1].location
    );
    proposed
        .validate_binding(&d::propose_revision(&s, &r).unwrap())
        .unwrap();
}
