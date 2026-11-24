//! Independently assembled little-endian fixtures, no index generator/tool.
use collagene_genomic_index::{Association, Format, Index, Limits, QueryLimits, VirtualOffset};
fn u32_(b: &mut Vec<u8>, n: u32) {
    b.extend_from_slice(&n.to_le_bytes());
}
fn u64_(b: &mut Vec<u8>, n: u64) {
    b.extend_from_slice(&n.to_le_bytes());
}
fn bin(b: &mut Vec<u8>, id: u32, csi: bool, pairs: &[(u64, u64)]) {
    u32_(b, id);
    if csi {
        u64_(b, 0);
    }
    u32_(b, pairs.len() as u32);
    for &(a, z) in pairs {
        u64_(b, a);
        u64_(b, z);
    }
}
fn bai() -> Vec<u8> {
    let mut b = b"BAI\x01".to_vec();
    u32_(&mut b, 1);
    u32_(&mut b, 4);
    // Deliberately unsorted bins; metadata counters are not virtual offsets.
    bin(&mut b, 4682, false, &[(0x18000, 0x30000)]);
    bin(&mut b, 37450, false, &[(0x10000, 0x40000), (7, u64::MAX)]);
    bin(&mut b, 0, false, &[(0x30000, 0x40000)]);
    bin(&mut b, 4681, false, &[(0x10000, 0x20000)]);
    u32_(&mut b, 2);
    u64_(&mut b, 0x10000);
    u64_(&mut b, 0x18000);
    u64_(&mut b, 11);
    b
}
fn pairs(index: &Index, start: u64, end: u64) -> Vec<(u64, u64)> {
    index
        .query(0, start, end, &QueryLimits::default())
        .unwrap()
        .iter()
        .map(|c| (c.start().packed(), c.end().packed()))
        .collect()
}
#[test]
fn literal_bai_boundaries_merge_candidates_and_separate_metadata() {
    let bytes = bai();
    let index = Index::parse(&bytes, Format::Bai, &Limits::default()).unwrap();
    assert_eq!(index.original_bytes(), bytes);
    assert_eq!(index.reference_count(), 1);
    assert_eq!(index.coordinate_limit(), 1 << 29);
    assert_eq!(
        index.references()[0]
            .bins()
            .iter()
            .map(|b| b.id())
            .collect::<Vec<_>>(),
        [0, 4681, 4682]
    );
    assert_eq!(
        index.references()[0].metadata().unwrap().placed_unmapped,
        u64::MAX
    );
    assert_eq!(index.unplaced_count(), Some(11));
    assert_eq!(
        pairs(&index, 0, 1),
        [(0x10000, 0x20000), (0x30000, 0x40000)]
    );
    assert_eq!(pairs(&index, 16383, 16385), [(0x10000, 0x40000)]);
    assert_eq!(pairs(&index, 16384, 16385), [(0x18000, 0x40000)]);
    assert!(index.query(0, 5, 5, &QueryLimits::default()).is_err());
    assert!(
        index
            .query(0, 0, (1 << 29) + 1, &QueryLimits::default())
            .is_err()
    );
    assert!(index.query(1, 0, 1, &QueryLimits::default()).is_err());
}
#[test]
fn literal_csi_generalized_coordinates_and_auxiliary_bytes() {
    let mut b = b"CSI\x01".to_vec();
    u32_(&mut b, 4);
    u32_(&mut b, 2);
    u32_(&mut b, 3);
    b.extend_from_slice(&[0, 255, 17]);
    u32_(&mut b, 1);
    u32_(&mut b, 3);
    bin(&mut b, 9, true, &[(0x10000, 0x20000)]);
    bin(&mut b, 10, true, &[(0x20000, 0x30000)]);
    bin(&mut b, 74, true, &[(0x10000, 0x30000), (100, 2)]);
    let i = Index::parse(&b, Format::Csi, &Limits::default()).unwrap();
    assert_eq!(
        (i.min_shift(), i.depth(), i.coordinate_limit()),
        (4, 2, 1024)
    );
    assert_eq!(i.auxiliary_bytes(), [0, 255, 17]);
    assert_eq!(pairs(&i, 15, 17), [(0x10000, 0x30000)]);
    assert_eq!(pairs(&i, 16, 17), [(0x20000, 0x30000)]);
    assert!(pairs(&i, 32, 33).is_empty());
    assert_eq!(i.references()[0].metadata().unwrap().mapped, 100);
}
#[test]
fn explicit_association_pins_source_identity_but_is_not_index_authentication() {
    let i = Index::parse(&bai(), Format::Bai, &Limits::default()).unwrap();
    let association = Association {
        source_blake3: "a".repeat(64),
        source_bytes: 4,
        reference_names: vec!["declared_reference".into()],
    };
    let b = i.bind(&association).unwrap();
    assert_eq!(b.association(), &association);
    assert_eq!(
        b.query(&"a".repeat(64), 0, 0, 1, &QueryLimits::default())
            .unwrap()
            .len(),
        2
    );
    assert!(
        b.query(&"b".repeat(64), 0, 0, 1, &QueryLimits::default())
            .is_err()
    );
    let mut bad = association.clone();
    bad.source_bytes = 3;
    assert!(i.bind(&bad).is_err());
    bad = association.clone();
    bad.reference_names.clear();
    assert!(i.bind(&bad).is_err());
    let o = VirtualOffset::new((1 << 48) - 1, 65535).unwrap();
    assert_eq!(o.raw(), u64::MAX);
    assert!(VirtualOffset::new(1 << 48, 0).is_err());
}
#[test]
fn truncation_trailing_sparse_caps_and_duplicate_bins_reject_atomically() {
    let bytes = bai();
    for n in 0..bytes.len() {
        if n == bytes.len() - 8 {
            continue;
        }
        assert!(
            Index::parse(&bytes[..n], Format::Bai, &Limits::default()).is_err(),
            "prefix {n}"
        );
    }
    let mut extra = bytes.clone();
    extra.push(0);
    assert!(Index::parse(&extra, Format::Bai, &Limits::default()).is_err());
    let mut l = Limits::default();
    l.chunks = 1;
    assert!(Index::parse(&bytes, Format::Bai, &l).is_err());
    l = Limits::default();
    l.retained_bytes = 10;
    assert!(Index::parse(&bytes, Format::Bai, &l).is_err());
    let mut duplicate = b"BAI\x01".to_vec();
    u32_(&mut duplicate, 1);
    u32_(&mut duplicate, 2);
    bin(&mut duplicate, 0, false, &[]);
    bin(&mut duplicate, 0, false, &[]);
    u32_(&mut duplicate, 0);
    assert!(Index::parse(&duplicate, Format::Bai, &Limits::default()).is_err());
    let mut reserved = b"BAI\x01".to_vec();
    u32_(&mut reserved, 1);
    u32_(&mut reserved, 1);
    bin(&mut reserved, 37449, false, &[]);
    u32_(&mut reserved, 0);
    assert!(Index::parse(&reserved, Format::Bai, &Limits::default()).is_err());
    let i = Index::parse(&bytes, Format::Bai, &Limits::default()).unwrap();
    assert!(
        i.query(
            0,
            0,
            1,
            &QueryLimits {
                chunks: 1,
                work: 1000
            }
        )
        .is_err()
    );
    assert!(
        i.query(
            0,
            0,
            1,
            &QueryLimits {
                chunks: 100,
                work: 0
            }
        )
        .is_err()
    );
}
#[test]
fn csi_shift_depth_and_signed_lengths_are_checked_without_unsafe_shifts() {
    for (shift, depth) in [(u32::MAX, 0), (64, 0), (34, 10), (0, 11)] {
        let mut b = b"CSI\x01".to_vec();
        u32_(&mut b, shift);
        u32_(&mut b, depth);
        u32_(&mut b, 0);
        u32_(&mut b, 0);
        assert!(Index::parse(&b, Format::Csi, &Limits::default()).is_err());
    }
    let mut b = b"CSI\x01".to_vec();
    u32_(&mut b, 63);
    u32_(&mut b, 0);
    u32_(&mut b, 0);
    u32_(&mut b, 1);
    u32_(&mut b, 1);
    bin(&mut b, 0, true, &[(1, 2)]);
    let i = Index::parse(&b, Format::Csi, &Limits::default()).unwrap();
    assert_eq!(i.coordinate_limit(), 1u64 << 63);
    assert_eq!(pairs(&i, (1u64 << 63) - 1, 1u64 << 63), [(1, 2)]);
}
