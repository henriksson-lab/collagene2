//! One length-framed, domain-separated content identity representation.
use super::*;
fn bytes(hash: &mut blake3::Hasher, value: &[u8]) {
    hash.update(&(value.len() as u64).to_le_bytes());
    hash.update(value);
}
fn number(hash: &mut blake3::Hasher, value: usize) {
    hash.update(&(value as u64).to_le_bytes());
}
pub(super) fn source(source: &Sequence, topology: Topology) -> String {
    let mut h = blake3::Hasher::new();
    bytes(&mut h, b"oligo-binding-source-1");
    bytes(&mut h, source.as_str().as_bytes());
    h.update(&[match topology {
        Topology::Linear => 0,
        Topology::Circular => 1,
    }]);
    h.finalize().to_hex().to_string()
}
pub(super) fn collection(collection: &Collection) -> String {
    let mut h = blake3::Hasher::new();
    bytes(&mut h, b"oligo-binding-collection-1");
    number(&mut h, collection.oligos.len());
    for o in &collection.oligos {
        bytes(&mut h, o.id.as_bytes());
        bytes(&mut h, o.name.as_bytes());
        bytes(&mut h, o.ordered.as_str().as_bytes());
        number(&mut h, o.declared_tail_length);
    }
    h.finalize().to_hex().to_string()
}
pub(super) fn options(options: &Options) -> String {
    let mut h = blake3::Hasher::new();
    bytes(&mut h, MODEL.as_bytes());
    for n in [
        options.minimum_annealing_length,
        options.maximum_extra_unpaired,
        options.maximum_mismatches,
        options.minimum_three_prime_pairing,
    ] {
        number(&mut h, n);
    }
    for n in [
        options.match_score,
        options.mismatch_score,
        options.minimum_score,
    ] {
        h.update(&n.to_le_bytes());
    }
    h.finalize().to_hex().to_string()
}
