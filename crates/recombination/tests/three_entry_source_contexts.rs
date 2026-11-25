use collagene_model::{Molecule, Sequence, Topology};
use collagene_recombination::{Family, Input, Limits, Reaction, multisite};
fn rc(s: &str) -> String {
    Sequence::dna(s)
        .unwrap()
        .reverse_complement()
        .unwrap()
        .as_str()
        .into()
}
#[test]
fn finite_profile_pins_exact_contexts_source_spans_orientations_and_nested_p3_competitor() {
    let p = multisite::three_entry_source_contexts().unwrap();
    assert_eq!(p.id, "tol2kit-228-383-302-vectorii-lr-reference");
    assert_eq!(p.edition, "2026-11-25.1");
    assert_eq!(p.variants.len(), 9);
    assert_eq!(p.sources.len(), 9);
    assert_eq!(
        p.variants.iter().map(|v| v.sequence.len()).sum::<usize>(),
        1188
    );
    let facts = [
        (Family::L, 4, 593, 688, false, 82, "GTATAGA"),
        (Family::R, 1, 856, 1012, false, 7, "GTACAAA"),
        (Family::L, 1, 569, 668, false, 83, "GTACAAA"),
        (Family::L, 2, 1452, 1551, true, 83, "GTACAAG"),
        (Family::R, 2, 591, 747, true, 7, "GTACAAG"),
        (Family::L, 3, 947, 1042, true, 82, "GTATAAT"),
        (Family::R, 4, 37, 161, false, 7, "GTATAGA"),
        (Family::R, 3, 2064, 2188, true, 7, "GTATAAT"),
        (Family::P, 3, 2746, 2977, true, 82, "GTATAAT"),
    ];
    for (family, subtype, start, end, reverse, offset, core) in facts {
        let v = p
            .variants
            .iter()
            .find(|v| v.family == family && v.subtype == subtype)
            .unwrap();
        let source = p.sources.iter().find(|s| s.id == v.source_id).unwrap();
        assert_eq!(
            (
                source.start_1based,
                source.end_1based,
                source.reverse_extraction
            ),
            (start, end, reverse)
        );
        assert_eq!(v.sequence.len(), end - start + 1);
        assert_eq!(v.overlap_start, offset);
        assert_eq!(&v.sequence[offset..offset + 7], core);
        assert_eq!(source.raw_sha256.len(), 64);
        assert_eq!(source.map_sha256.len(), 64);
        assert!(source.url.starts_with("https://"));
        assert!(source.map_url.starts_with("https://"));
    }
    let sources = [
        (
            "p5E-MCS",
            2810,
            "9a4cfe981fd870a56783a0af1641cd59f884f66b1eb0abcad1b756ca0a8c8e39",
            "https://tol2kitkwan.genetics.utah.edu/images/0/06/228.p5E-MCS.zip",
        ),
        (
            "pME-EGFP",
            3327,
            "bd4677c48120be4c1655ccb6444032d0e489882a71d215b00db5abfebf6be3cf",
            "https://tol2kitkwan.genetics.utah.edu/images/e/e9/383.pME-EGFP.zip",
        ),
        (
            "p3E-polyA",
            2838,
            "2b7d36c1ce28408040f22fbf05fcdcdaf27732134d7be59f9fdaad372398066b",
            "https://tol2kitkwan.genetics.utah.edu/images/f/f7/302.p3E-polyA.zip",
        ),
        (
            "pDEST R4-R3 Vector II",
            4555,
            "e91c6ae899edcfd029be1c796177f4e09466eb50e36a346edc239228858da3ea",
            "https://documents.thermofisher.com/TFS-Assets/vectors/pdestr4-r3_vector_ii_seq.txt",
        ),
        (
            "pDONR P2R-P3",
            4773,
            "16b796fec73d1e06b38f5feeec059aa20f7e149d503ca75e0fc605222c1b10d9",
            "https://documents.thermofisher.com/TFS-Assets/vectors/pdonrp2rp3_seq.txt",
        ),
    ];
    for source in &p.sources {
        let fact = sources.iter().find(|f| f.0 == source.vector_name).unwrap();
        assert_eq!(source.vector_length, fact.1);
        assert_eq!(source.raw_sha256, fact.2);
        assert_eq!(source.url, fact.3);
    }
    let l3 = p
        .variants
        .iter()
        .find(|v| v.family == Family::L && v.subtype == 3)
        .unwrap();
    assert_eq!(
        l3.sequence,
        "AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATAAGCAATGCTTTTTTATAATGCCAACTTTGTATAATAAAGTTG"
    );
    let p3 = p
        .variants
        .iter()
        .find(|v| v.family == Family::P && v.subtype == 3)
        .unwrap();
    assert_eq!(&p3.sequence[..96], l3.sequence);
    assert_eq!(p3.sequence.len(), 232);
    let original = p.revision().unwrap();
    let mut changed = p.clone();
    changed.sources[0].raw_sha256 = "0".repeat(64);
    assert_ne!(changed.revision().unwrap(), original);
    let mut changed = p;
    changed.sources[0].end_1based -= 1;
    assert!(changed.validate().is_err());
}
#[test]
fn four_synthetic_parents_using_real_context_facts_make_independently_frozen_four_scar_reference() {
    let profile = multisite::three_entry_source_contexts().unwrap();
    let site = |family, subtype| {
        profile
            .variants
            .iter()
            .find(|v| v.family == family && v.subtype == subtype)
            .unwrap()
            .sequence
            .clone()
    };
    // Functional contexts plus declared synthetic bodies only: these are not the
    // deposited plasmids. L2/R2/L3/R3 occur in reverse canonical orientation.
    let dna = [
        format!("{}CCCCA{}AAAAC", site(Family::L, 4), site(Family::R, 1)),
        format!(
            "{}GGGGT{}TTTTC",
            site(Family::L, 1),
            rc(&site(Family::L, 2))
        ),
        format!(
            "{}AAAAG{}CCCCT",
            rc(&site(Family::R, 2)),
            rc(&site(Family::L, 3))
        ),
        format!(
            "{}TTTTA{}GGGGA",
            site(Family::R, 4),
            rc(&site(Family::R, 3))
        ),
    ];
    let inputs = dna
        .iter()
        .enumerate()
        .map(|(i, s)| Input {
            role: format!("synthetic-parent-{i}"),
            molecule: Molecule::dna(
                format!("Synthetic context parent{i}"),
                s,
                Topology::Circular,
            )
            .unwrap(),
        })
        .collect::<Vec<_>>();
    // Independently frozen before invoking the graph: core4 + L4 suffix/body/R1
    // prefix; core1 + L1 suffix/body/RC(L2) prefix; RC(core2) + RC(R2)
    // suffix/body/RC(L3) prefix; RC(core3) + RC(R3) suffix/body/R4 prefix.
    const EXPECTED: &str = "GTATAGAAAAGTTGCCCCACAAGTTTGTACAAAAAAGCAGGCTGGGGTACCCAGCTTTCTTGTACAAAGTGGAAAAGCAACTTTATTATACATAGTTGGGGGACAACTTT";
    assert_eq!(EXPECTED.len(), 110);
    assert_eq!(
        inputs
            .iter()
            .map(|i| i.molecule.sequence.len())
            .sum::<usize>(),
        996
    );
    let request = multisite::Request {
        profile,
        reaction: Reaction::LR,
        specificities: vec![4, 1, 2, 3],
        limits: Limits::default(),
    };
    let report = multisite::recombine(&inputs, &request).unwrap();
    assert!(report.complete);
    assert_eq!(report.candidates.len(), 1);
    assert!(!report.sites.iter().any(|h| h.family == Family::P));
    let candidate = &report.candidates[0];
    assert_eq!(candidate.products.len(), 2);
    let main = candidate
        .products
        .iter()
        .find(|p| {
            let s = p.molecule.sequence.as_str();
            s.len() == EXPECTED.len()
                && (format!("{s}{s}").contains(EXPECTED)
                    || format!("{}{}", rc(s), rc(s)).contains(EXPECTED))
        })
        .expect("independent110bp synthetic intended reference");
    assert_eq!(main.junctions.len(), 4);
    assert!(main.junctions.iter().all(|j| j.family == Family::B));
    assert!(
        candidate
            .products
            .iter()
            .any(|p| p.molecule.sequence.len() == 886
                && p.junctions.iter().all(|j| j.family == Family::P))
    );
    multisite::audit_candidate(&inputs, &request.profile, candidate).unwrap();
    for (input, parent) in inputs.iter().enumerate() {
        let mut counts = vec![0u8; parent.molecule.sequence.len()];
        for map in candidate
            .products
            .iter()
            .flat_map(|p| &p.segments)
            .filter(|m| m.input == input)
        {
            for n in &mut counts[map.source_start..map.source_end] {
                *n += 1;
            }
        }
        assert!(counts.iter().all(|n| *n == 1));
    }
}
