//! Original-coordinate chemical oracle independent of implementation mappings.
use collagene_core::{
    chemistry::{
        Strand,
        methylation::MethylationState as M,
        tracked_duplex::{
            CallState, Chemical, GenerationAssumptions, StrandModifications, TrackedDuplex, digest,
            ligate,
        },
    },
    sequence::analysis::Enzyme,
};
use collagene_model::{Molecule, Topology};
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
fn rotate(s: &str, n: usize) -> String {
    format!("{}{}", &s[n..], &s[..n])
}
fn state(s: &StrandModifications, c: Chemical, i: usize) -> CallState {
    s.spans
        .iter()
        .find(|s| s.chemical == c && s.start <= i && i < s.end)
        .map(|s| s.state)
        .unwrap_or(match c {
            Chemical::N6A => s.n6a_default,
            Chemical::C5C => s.c5c_default,
        })
}
fn expected(seq: &str, c: Chemical, i: usize) -> CallState {
    let n = seq.len();
    for (pattern, offset, chemical) in [
        ("GATC", 1, Chemical::N6A),
        ("CCAGG", 1, Chemical::C5C),
        ("CCTGG", 1, Chemical::C5C),
        ("CG", 0, Chemical::C5C),
    ] {
        if c != chemical {
            continue;
        }
        for start in 0..n {
            if (start + offset) % n == i
                && (0..pattern.len())
                    .all(|j| seq.as_bytes()[(start + j) % n] == pattern.as_bytes()[j])
            {
                return CallState::Modified;
            }
        }
    }
    CallState::Unmodified
}
fn assert_mapped(m: &StrandModifications, source: &str, offset: usize) {
    for i in 0..source.len() {
        let old = (offset + i) % source.len();
        let c = match source.as_bytes()[old] {
            b'A' => Chemical::N6A,
            b'C' => Chemical::C5C,
            _ => continue,
        };
        assert_eq!(
            state(m, c, i),
            expected(source, c, old),
            "old coordinate {old}, new {i}, offset {offset}"
        );
    }
}
#[test]
fn all_origins_and_orientations_preserve_actual_strand_chemical_positions_through_digest_and_closure()
 {
    let base = "AAGGATCCAACCGGCCAAGCGGCCGCTTAACCAGGTT";
    let enzyme = Enzyme {
        name: "BamHI".into(),
        recognition: "GGATCC".into(),
        top_cut: 1,
        bottom_cut: 5,
    };
    for reverse in [false, true] {
        let original = if reverse { rc(base) } else { base.into() };
        for origin in 0..original.len() {
            let text = rotate(&original, origin);
            let bottom = rc(&text);
            let n = text.len();
            let m = Molecule::dna("independent ledger", &text, Topology::Circular).unwrap();
            let tracked = TrackedDuplex::from_input(
                &m,
                &GenerationAssumptions {
                    dam: M::FullyModified,
                    dcm: M::FullyModified,
                    cpg: M::FullyModified,
                },
            )
            .unwrap();
            assert_mapped(&tracked.top_modifications, &text, 0);
            assert_mapped(&tracked.bottom_modifications, &bottom, 0);
            let reverse_twice = tracked.reverse().unwrap().reverse().unwrap();
            assert_eq!(reverse_twice, tracked);
            let site = (0..n)
                .find(|&start| (0..6).all(|j| text.as_bytes()[(start + j) % n] == b"GGATCC"[j]))
                .unwrap();
            let top = (site + 1) % n;
            let bottom_cut = top + 4;
            let fs = digest(&tracked, std::slice::from_ref(&enzyme)).unwrap();
            assert_eq!(fs.len(), 1);
            let f = &fs[0];
            assert_eq!((f.top_start, f.bottom_start), (top, bottom_cut as i64));
            assert_eq!(f.duplex.duplex.top.as_str(), rotate(&text, top));
            let bi = (n - bottom_cut % n) % n;
            assert_eq!(f.duplex.duplex.bottom.as_str(), rotate(&bottom, bi));
            assert_mapped(&f.duplex.top_modifications, &text, top);
            assert_mapped(&f.duplex.bottom_modifications, &bottom, bi);
            let closed = ligate(std::slice::from_ref(&f.duplex), true, false).unwrap();
            assert_eq!(closed.duplex.topology, Topology::Circular);
            assert_eq!(closed.duplex.top.as_str(), rotate(&text, top));
            assert_eq!(
                closed.duplex.bottom.as_str(),
                rc(closed.duplex.top.as_str())
            );
            assert_mapped(&closed.top_modifications, &text, top);
            assert_mapped(&closed.bottom_modifications, &bottom, (n - top) % n);
            assert_eq!(
                closed.summary().unwrap().known_modified,
                tracked.summary().unwrap().known_modified
            );
        }
    }
}
#[test]
fn reversal_swaps_actual_five_prime_ledgers_without_complementing_their_chemistry() {
    let m = Molecule::dna(
        "strand chemical ledger",
        "AAGATCAACGCCAGGTT",
        Topology::Linear,
    )
    .unwrap();
    let t = TrackedDuplex::from_input(
        &m,
        &GenerationAssumptions {
            dam: M::FullyModified,
            dcm: M::FullyModified,
            cpg: M::FullyModified,
        },
    )
    .unwrap();
    let r = t.reverse().unwrap();
    assert_eq!(r.duplex.top, t.duplex.bottom);
    assert_eq!(r.duplex.bottom, t.duplex.top);
    for c in [Chemical::N6A, Chemical::C5C] {
        for i in 0..t.duplex.top.len() {
            assert_eq!(
                state(&r.top_modifications, c, i),
                state(&t.bottom_modifications, c, i)
            );
            assert_eq!(
                state(&r.bottom_modifications, c, i),
                state(&t.top_modifications, c, i)
            );
        }
    }
    let _strand = Strand::Bottom;
}
#[test]
fn fill_in_keeps_original_dam_calls_and_marks_each_added_adenine_unmodified() {
    let m = Molecule::dna(
        "literal repair oracle",
        "GAATTCAAAAGATC",
        Topology::Circular,
    )
    .unwrap();
    let t = TrackedDuplex::from_input(
        &m,
        &GenerationAssumptions {
            dam: M::FullyModified,
            dcm: M::Unmodified,
            cpg: M::Unmodified,
        },
    )
    .unwrap();
    let enzyme = Enzyme {
        name: "EcoRI".into(),
        recognition: "GAATTC".into(),
        top_cut: 1,
        bottom_cut: 5,
    };
    let f = digest(&t, &[enzyme]).unwrap().remove(0).duplex;
    assert_eq!(f.duplex.top.as_str(), "AATTCAAAAGATCG");
    assert_eq!(
        state(&f.top_modifications, Chemical::N6A, 10),
        CallState::Modified
    );
    assert_eq!(
        state(&f.bottom_modifications, Chemical::N6A, 6),
        CallState::Modified
    );
    let repaired = f
        .repair(collagene_core::chemistry::EndRepairMethod::FillFivePrimeOverhangs)
        .unwrap();
    assert_eq!(repaired.duplex.top.as_str(), "AATTCAAAAGATCGAATT");
    assert_eq!(
        repaired.duplex.bottom,
        repaired.duplex.top.reverse_complement().unwrap()
    );
    assert_eq!(
        state(&repaired.top_modifications, Chemical::N6A, 10),
        CallState::Modified
    );
    assert_eq!(
        state(&repaired.bottom_modifications, Chemical::N6A, 6),
        CallState::Modified
    );
    for s in [&repaired.top_modifications, &repaired.bottom_modifications] {
        for i in [14, 15] {
            assert_eq!(state(s, Chemical::N6A, i), CallState::Unmodified);
        }
    }
    assert_eq!(repaired.summary().unwrap().known_modified, 2);
    assert_eq!(
        repaired.summary().unwrap().known_unmodified,
        f.summary().unwrap().known_unmodified + 4
    );
    let closed = ligate(&[repaired], true, false).unwrap();
    assert_eq!(closed.duplex.top.as_str(), "AATTCAAAAGATCGAATT");
    assert_eq!(closed.summary().unwrap().known_modified, 2);
}
