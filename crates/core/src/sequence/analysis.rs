use collagene_model::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub fn base_mask(base: u8) -> u8 {
    match base {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Match {
    pub start: usize,
    pub length: usize,
    pub reverse: bool,
    pub certain: bool,
}

/// Borrowed coordinates of a real selected source strand. This is neither a
/// replacement document nor a serialized material authority.
#[derive(Clone, Copy)]
pub struct TemplateView<'a> {
    pub sequence: &'a Sequence,
    pub topology: Topology,
    pub template_strand: crate::chemistry::Strand,
}
impl<'a> TemplateView<'a> {
    pub fn reference(source: &'a Molecule) -> Self {
        Self {
            sequence: &source.sequence,
            topology: source.topology,
            template_strand: crate::chemistry::Strand::Top,
        }
    }
}
/// Ambiguous template bases give possible matches, never definitive binding evidence.
pub fn find_matches(mol: &Molecule, query: &str, both_strands: bool) -> Result<Vec<Match>, Error> {
    find_matches_view(TemplateView::reference(mol), query, both_strands)
}
/// Shared limit for one exact search or all four exact-PCR binding searches.
/// Concrete inputs use linear KMP; ambiguous comparisons consume this budget.
pub const MAX_EXACT_SEARCH_WORK: usize = 200_000_000;
const EXACT_WORK_ERROR: &str =
    "exact search/PCR binding exceeded its 200-million work budget; shorten or refine the inputs";

/// One actual-work ledger for exact binding searches and PCR materialization.
/// The same ledger can span independent operations; failed searches keep charges.
/// It is native execution state and is never accepted from serialized evidence.
pub struct ExactSearchBudget {
    work: super::matching::Work,
}
impl Default for ExactSearchBudget {
    fn default() -> Self {
        Self {
            work: super::matching::Work::new(MAX_EXACT_SEARCH_WORK),
        }
    }
}
impl ExactSearchBudget {
    pub fn new(limit: usize) -> Result<Self, Error> {
        if limit == 0 || limit > MAX_EXACT_SEARCH_WORK {
            return Err(Error::TooLarge);
        }
        Ok(Self {
            work: super::matching::Work::new(limit),
        })
    }
    pub fn used(&self) -> usize {
        self.work.used()
    }
    pub fn remaining(&self) -> usize {
        self.work.remaining()
    }
    pub(crate) fn charge(&mut self, amount: usize) -> Result<(), Error> {
        if self.work.take(amount) {
            Ok(())
        } else {
            Err(Error::Unsupported(EXACT_WORK_ERROR.into()))
        }
    }
}

/// Same exact/IUPAC search with a caller-owned actual-work ledger.
pub fn find_matches_with_budget(
    mol: &Molecule,
    query: &str,
    both_strands: bool,
    budget: &mut ExactSearchBudget,
) -> Result<Vec<Match>, Error> {
    find_with_work(
        TemplateView::reference(mol),
        query,
        both_strands,
        &mut budget.work,
    )
}

pub fn find_matches_view(
    mol: TemplateView<'_>,
    query: &str,
    both_strands: bool,
) -> Result<Vec<Match>, Error> {
    let mut work = super::matching::Work::new(MAX_EXACT_SEARCH_WORK);
    find_with_work(mol, query, both_strands, &mut work)
}
fn find_with_work(
    mol: TemplateView<'_>,
    query: &str,
    both_strands: bool,
    work: &mut super::matching::Work,
) -> Result<Vec<Match>, Error> {
    if mol.sequence.alphabet() == Alphabet::Protein {
        return Err(Error::Unsupported("nucleotide search on protein".into()));
    }
    if !work.take(query.len()) {
        return Err(Error::Unsupported(EXACT_WORK_ERROR.into()));
    }
    let pattern = Sequence::new(mol.sequence.alphabet(), query)?;
    if pattern.len() > mol.sequence.len() {
        return Ok(vec![]);
    }
    let mut patterns = vec![(pattern, false)];
    if both_strands {
        if !work.take(patterns[0].0.len()) {
            return Err(Error::Unsupported(EXACT_WORK_ERROR.into()));
        }
        let rc = patterns[0].0.reverse_complement()?;
        if rc != patterns[0].0 {
            patterns.push((rc, true));
        }
    }
    let mut hits = vec![];
    for (pattern, reverse) in patterns {
        let mut hit_limit = false;
        let complete = super::matching::scan(
            mol.sequence,
            &pattern,
            mol.topology == Topology::Circular,
            work,
            |placement| {
                if hits.len() >= 100_000 {
                    hit_limit = true;
                    return false;
                }
                hits.push(Match {
                    start: placement.start,
                    length: pattern.len(),
                    reverse,
                    certain: placement.certain,
                });
                true
            },
        )?;
        if !complete {
            return Err(Error::Unsupported(if hit_limit {
                "search has over 100,000 matches; use a longer query".into()
            } else {
                EXACT_WORK_ERROR.into()
            }));
        }
    }
    hits.sort_by_key(|h| (h.start, h.reverse));
    Ok(hits)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Enzyme {
    pub name: String,
    pub recognition: String,
    pub top_cut: i32,
    pub bottom_cut: i32,
}

/// A small factual starter set, not a supplier-condition or methylation database.
pub fn starter_enzymes() -> Vec<Enzyme> {
    [
        ("EcoRI", "GAATTC", 1, 5),
        ("BamHI", "GGATCC", 1, 5),
        ("HindIII", "AAGCTT", 1, 5),
        ("EcoRV", "GATATC", 3, 3),
        ("SmaI", "CCCGGG", 3, 3),
        ("NotI", "GCGGCCGC", 2, 6),
        ("XbaI", "TCTAGA", 1, 5),
        ("BsaI", "GGTCTC", 7, 11),
        ("PstI", "CTGCAG", 5, 1),
        ("BsmBI", "CGTCTC", 7, 11),
        ("BbsI", "GAAGAC", 8, 12),
        ("SapI", "GCTCTTC", 8, 11),
    ]
    .into_iter()
    .map(|(n, r, t, b)| Enzyme {
        name: n.into(),
        recognition: r.into(),
        top_cut: t,
        bottom_cut: b,
    })
    .collect()
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CutSite {
    pub enzyme: String,
    pub recognition_start: usize,
    pub reverse: bool,
    pub top_cut: usize,
    pub bottom_cut: usize,
    pub certain: bool,
}

pub fn restriction_sites(mol: &Molecule, enzyme: &Enzyme) -> Result<Vec<CutSite>, Error> {
    restriction_sites_with_budget(mol, enzyme, &mut ExactSearchBudget::default())
}
/// Same native recognition/cut kernel with a caller-owned actual matching ledger.
/// Failed searches retain their charges and do not publish partial cut lists.
pub fn restriction_sites_with_budget(
    mol: &Molecule,
    enzyme: &Enzyme,
    budget: &mut ExactSearchBudget,
) -> Result<Vec<CutSite>, Error> {
    if mol.sequence.alphabet() != Alphabet::Dna {
        return Err(Error::Unsupported(
            "restriction analysis requires DNA".into(),
        ));
    }
    let recognition = Sequence::dna(&enzyme.recognition)?;
    let n = mol.sequence.len() as i64;
    let m = recognition.len() as i64;
    let mut sites = vec![];
    for hit in find_matches_with_budget(mol, &enzyme.recognition, true, budget)? {
        let offset = hit.start as i64;
        let (top, bottom) = if hit.reverse {
            (
                offset + m - i64::from(enzyme.bottom_cut),
                offset + m - i64::from(enzyme.top_cut),
            )
        } else {
            (
                offset + i64::from(enzyme.top_cut),
                offset + i64::from(enzyme.bottom_cut),
            )
        };
        if mol.topology == Topology::Linear && (top < 0 || bottom < 0 || top > n || bottom > n) {
            return Err(Error::Unsupported(format!(
                "{} recognition site at base {} lacks sequence for both modeled strand cuts",
                enzyme.name,
                hit.start + 1
            )));
        }
        let pos = |v: i64| {
            if mol.topology == Topology::Circular {
                v.rem_euclid(n) as usize
            } else {
                v as usize
            }
        };
        sites.push(CutSite {
            enzyme: enzyme.name.clone(),
            recognition_start: hit.start,
            reverse: hit.reverse,
            top_cut: pos(top),
            bottom_cut: pos(bottom),
            certain: hit.certain,
        });
    }
    Ok(sites)
}

/// Predicted lengths along the reference/top strand. Assumes cleavage; no methylation inference.
pub fn digest_lengths(mol: &Molecule, sites: &[CutSite]) -> Result<Vec<usize>, Error> {
    let n = mol.sequence.len();
    if sites.iter().any(|s| !s.certain) {
        return Err(Error::Unsupported(
            "ambiguous restriction sites cannot establish a complete digest".into(),
        ));
    }
    if sites.iter().any(|s| s.top_cut > n || s.bottom_cut > n) {
        return Err(Error::Coordinates("cut outside molecule".into()));
    }
    let cuts: BTreeSet<usize> = sites
        .iter()
        .map(|s| {
            if mol.topology == Topology::Circular {
                s.top_cut % n
            } else {
                s.top_cut
            }
        })
        .collect();
    if cuts.is_empty() {
        return Ok(vec![n]);
    }
    let mut positions: Vec<_> = cuts.into_iter().collect();
    let mut lengths = vec![];
    if mol.topology == Topology::Linear {
        positions.insert(0, 0);
        positions.push(n);
    }
    for pair in positions.windows(2) {
        if pair[1] > pair[0] {
            lengths.push(pair[1] - pair[0]);
        }
    }
    if mol.topology == Topology::Circular {
        lengths.push(n - positions.last().unwrap() + positions[0]);
    }
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    Ok(lengths)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Primer {
    pub name: String,
    pub annealing: Sequence,
    pub tail: String,
}
impl Primer {
    pub fn new(name: impl Into<String>, annealing: &str, tail: &str) -> Result<Self, Error> {
        let annealing = Sequence::dna(annealing)?;
        let tail = if tail.trim().is_empty() {
            String::new()
        } else {
            Sequence::dna(tail)?.as_str().into()
        };
        Ok(Self {
            name: name.into(),
            annealing,
            tail,
        })
    }
    pub fn ordered_sequence(&self) -> String {
        format!("{}{}", self.tail, self.annealing.as_str())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PcrProduct {
    pub template_strand: crate::chemistry::Strand,
    pub sequence: Sequence,
    pub forward_start: usize,
    pub reverse_start: usize,
    pub crosses_origin: bool,
    pub forward_primer: String,
    pub reverse_primer: String,
}

/// Enumerate exact, non-overlapping binding geometries for both oligos in both roles.
/// Includes possible single-oligo amplification; primer names do not restrict binding.
/// Geometry retains complete selected-frame coordinates and synthesized tails.
#[derive(Debug, Clone)]
pub struct PcrGeometryProduct {
    pub product: PcrProduct,
    pub template_length: usize,
    pub forward_tail_length: usize,
    pub reverse_tail_length: usize,
}
pub fn simulate_pcr(
    mol: &Molecule,
    forward: &Primer,
    reverse: &Primer,
) -> Result<Vec<PcrProduct>, Error> {
    let mut seen = BTreeSet::new();
    Ok(simulate_pcr_geometry(mol, forward, reverse)?
        .into_iter()
        .filter_map(|g| {
            let p = g.product;
            let key = (
                p.forward_start,
                p.reverse_start,
                blake3::hash(p.sequence.as_str().as_bytes())
                    .to_hex()
                    .to_string(),
            );
            seen.insert(key).then_some(p)
        })
        .collect())
}
pub fn simulate_pcr_geometry(
    mol: &Molecule,
    forward: &Primer,
    reverse: &Primer,
) -> Result<Vec<PcrGeometryProduct>, Error> {
    simulate_pcr_geometry_view(TemplateView::reference(mol), forward, reverse)
}
pub fn simulate_pcr_geometry_view(
    mol: TemplateView<'_>,
    forward: &Primer,
    reverse: &Primer,
) -> Result<Vec<PcrGeometryProduct>, Error> {
    simulate_pcr_geometry_view_with_budget(mol, forward, reverse, &mut ExactSearchBudget::default())
}
/// Canonical exact-PCR path with one caller-owned work budget. Charge preprocessing,
/// all four binding scans, pair visits and emitted source/product copies together.
pub fn simulate_pcr_geometry_view_with_budget(
    mol: TemplateView<'_>,
    forward: &Primer,
    reverse: &Primer,
    budget: &mut ExactSearchBudget,
) -> Result<Vec<PcrGeometryProduct>, Error> {
    if mol.sequence.alphabet() != Alphabet::Dna
        || forward.annealing.alphabet() != Alphabet::Dna
        || reverse.annealing.alphabet() != Alphabet::Dna
    {
        return Err(Error::Unsupported(
            "PCR requires DNA template and primers".into(),
        ));
    }
    for primer in [forward, reverse] {
        budget.charge(primer.name.len())?;
        budget.charge(primer.annealing.len())?;
        budget.charge(primer.tail.len())?;
    }
    // Public serialized primer fields are revalidated at the scientific boundary.
    let forward = Primer::new(&forward.name, forward.annealing.as_str(), &forward.tail)?;
    let reverse = Primer::new(&reverse.name, reverse.annealing.as_str(), &reverse.tail)?;
    if forward
        .annealing
        .as_str()
        .bytes()
        .chain(reverse.annealing.as_str().bytes())
        .any(|b| base_mask(b).count_ones() != 1)
    {
        return Err(Error::Unsupported(
            "degenerate primer annealing is not yet modeled".into(),
        ));
    }
    let primers = [forward, reverse];
    let bindings = primers
        .iter()
        .map(|primer| {
            let forward = find_with_work(mol, primer.annealing.as_str(), false, &mut budget.work)?;
            if !budget.work.take(primer.annealing.len()) {
                return Err(Error::Unsupported(EXACT_WORK_ERROR.into()));
            }
            let reverse_sequence = primer.annealing.reverse_complement()?;
            let reverse = find_with_work(mol, reverse_sequence.as_str(), false, &mut budget.work)?;
            if forward.iter().chain(&reverse).any(|h| !h.certain) {
                return Err(Error::Unsupported(
                    "ambiguous template binding prevents complete PCR enumeration".into(),
                ));
            }
            Ok((forward, reverse))
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let possibilities = bindings
        .iter()
        .map(|b| b.0.len())
        .sum::<usize>()
        .saturating_mul(bindings.iter().map(|b| b.1.len()).sum());
    if possibilities > 10_000 {
        return Err(Error::Unsupported("too many possible primer pairs".into()));
    }
    let n = mol.sequence.len();
    let mut products = vec![];
    let mut seen = BTreeSet::new();
    let mut product_bases = 0usize;
    for (fi, forward) in primers.iter().enumerate() {
        for (ri, reverse) in primers.iter().enumerate() {
            budget.charge(reverse.tail.len())?;
            let reverse_tail = if reverse.tail.is_empty() {
                String::new()
            } else {
                Sequence::dna(&reverse.tail)?
                    .reverse_complement()?
                    .as_str()
                    .to_owned()
            };
            for f in &bindings[fi].0 {
                for r in &bindings[ri].1 {
                    budget.charge(1)?;
                    let end = r.start + r.length;
                    let length = if mol.topology == Topology::Circular {
                        let distance = (end + n - f.start) % n;
                        if distance == 0 { n } else { distance }
                    } else if end > f.start {
                        end - f.start
                    } else {
                        continue;
                    };
                    // This baseline excludes overlapping annealing sites and multi-turn products.
                    if length < f.length + r.length || length > n {
                        continue;
                    }
                    let product_length = forward.tail.len() + length + reverse_tail.len();
                    if product_length > MAX_SEQUENCE_LEN
                        || product_bases.saturating_add(product_length) > 32_000_000
                    {
                        return Err(Error::Unsupported(
                            "PCR products exceed the materialization budget; refine the primers"
                                .into(),
                        ));
                    }
                    budget.charge(product_length)?;
                    let template: String = (0..length)
                        .map(|j| mol.sequence.as_str().as_bytes()[(f.start + j) % n] as char)
                        .collect();
                    let sequence =
                        Sequence::dna(&format!("{}{template}{reverse_tail}", forward.tail))?;
                    let key = (
                        f.start,
                        r.start,
                        blake3::hash(sequence.as_str().as_bytes())
                            .to_hex()
                            .to_string(),
                        forward.tail.len(),
                        reverse_tail.len(),
                        length,
                    );
                    if seen.insert(key) {
                        product_bases += product_length;
                        products.push(PcrGeometryProduct {
                            template_length: length,
                            forward_tail_length: forward.tail.len(),
                            reverse_tail_length: reverse_tail.len(),
                            product: PcrProduct {
                                template_strand: mol.template_strand,
                                sequence,
                                forward_start: f.start,
                                reverse_start: r.start,
                                crosses_origin: f.start + length > n,
                                forward_primer: forward.name.clone(),
                                reverse_primer: reverse.name.clone(),
                            },
                        });
                    }
                }
            }
        }
    }
    Ok(products)
}

/// Return the rotation offset if b is a rotation of a, in linear time.
pub fn rotation_offset(a: &str, b: &str) -> Option<usize> {
    if a.len() != b.len() || a.is_empty() {
        return None;
    }
    let pattern = b.as_bytes();
    let n = pattern.len();
    let mut prefix = vec![0; n];
    for i in 1..n {
        let mut j = prefix[i - 1];
        while j > 0 && pattern[i] != pattern[j] {
            j = prefix[j - 1];
        }
        if pattern[i] == pattern[j] {
            j += 1;
        }
        prefix[i] = j;
    }
    let mut j = 0;
    for (i, c) in a.bytes().chain(a.bytes().take(n - 1)).enumerate() {
        while j > 0 && c != pattern[j] {
            j = prefix[j - 1];
        }
        if c == pattern[j] {
            j += 1;
        }
        if j == n {
            return Some(i + 1 - n);
        }
    }
    None
}

#[cfg(test)]
mod restriction_budget_tests {
    use super::*;
    use crate::chemistry::{
        self,
        tracked_duplex::{self, CleavageBudget, TrackedDuplex},
    };

    #[test]
    fn plain_digest_shared_actual_work_preserves_literal_cuts_and_native_modification_projection() {
        let source = Molecule::dna("literal EcoRI", "AAAAGAATTCTTTT", Topology::Linear).unwrap();
        let enzyme = Enzyme {
            name: "EcoRI".into(),
            recognition: "GAATTC".into(),
            top_cut: 1,
            bottom_cut: 5,
        };
        let mut exact = ExactSearchBudget::new(100_000).unwrap();
        let sites = restriction_sites_with_budget(&source, &enzyme, &mut exact).unwrap();
        assert_eq!(sites.len(), 1);
        assert_eq!(
            (
                sites[0].recognition_start,
                sites[0].top_cut,
                sites[0].bottom_cut
            ),
            (4, 5, 9)
        );
        assert!(exact.used() > 0);
        let duplex = chemistry::Duplex::from_molecule(&source).unwrap();
        let mut shared = CleavageBudget::new(100_000).unwrap();
        let fragments =
            chemistry::digest_with_budget(&duplex, &[enzyme.clone()], &mut shared).unwrap();
        assert_eq!(shared.used(), exact.used());
        assert_eq!(
            fragments
                .iter()
                .map(|f| (f.duplex.top.len(), f.duplex.bottom.len()))
                .collect::<Vec<_>>(),
            [(5, 9), (9, 5)]
        );
        let ordinary = chemistry::digest(&duplex, &[enzyme.clone()]).unwrap();
        assert_eq!(
            serde_json::to_value(&fragments).unwrap(),
            serde_json::to_value(&ordinary).unwrap()
        );
        let tracked = TrackedDuplex::from_unknown_duplex(&duplex).unwrap();
        let mut tracked_budget = CleavageBudget::new(100_000).unwrap();
        let retained =
            tracked_duplex::digest_with_budget(&tracked, &[enzyme.clone()], &mut tracked_budget)
                .unwrap();
        let ordinary = tracked_duplex::digest(&tracked, &[enzyme]).unwrap();
        assert_eq!(tracked_budget.used(), exact.used());
        assert_eq!(
            serde_json::to_value(retained).unwrap(),
            serde_json::to_value(ordinary).unwrap()
        );
    }

    #[test]
    fn failed_plain_restriction_queries_keep_spent_work_and_cannot_reset_the_enclosing_budget() {
        let source = Molecule::dna("absent site", "ACGTACGTACGTACGT", Topology::Linear).unwrap();
        let enzyme = Enzyme {
            name: "EcoRI".into(),
            recognition: "GAATTC".into(),
            top_cut: 1,
            bottom_cut: 5,
        };
        let mut exact = ExactSearchBudget::new(100_000).unwrap();
        assert!(
            restriction_sites_with_budget(&source, &enzyme, &mut exact)
                .unwrap()
                .is_empty()
        );
        assert!(exact.used() > 0);
        let duplex = chemistry::Duplex::from_molecule(&source).unwrap();
        let mut shared = CleavageBudget::new(1).unwrap();
        assert!(chemistry::digest_with_budget(&duplex, &[enzyme.clone()], &mut shared).is_err());
        assert_eq!(shared.used(), 1);
        assert_eq!(shared.remaining(), 0);
        assert!(matches!(
            chemistry::digest_with_budget(&duplex, &[enzyme], &mut shared),
            Err(Error::TooLarge)
        ));
        assert_eq!(shared.used(), 1);

        // Native Type-IIS geometry refusal occurs AFTER actual recognition work.
        let edge =
            Molecule::dna("reverse Type-IIS edge", "TTTTGAGACCTTTT", Topology::Linear).unwrap();
        let bsa = Enzyme {
            name: "BsaI".into(),
            recognition: "GGTCTC".into(),
            top_cut: 7,
            bottom_cut: 11,
        };
        let mut attempted = ExactSearchBudget::new(100_000).unwrap();
        assert!(restriction_sites_with_budget(&edge, &bsa, &mut attempted).is_err());
        let mut shared = CleavageBudget::new(100_000).unwrap();
        let duplex = chemistry::Duplex::from_molecule(&edge).unwrap();
        assert!(chemistry::digest_with_budget(&duplex, &[bsa], &mut shared).is_err());
        assert_eq!(shared.used(), attempted.used());
        assert!(shared.used() > 0);
    }
}
