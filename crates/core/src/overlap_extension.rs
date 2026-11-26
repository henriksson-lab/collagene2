//! Selected reciprocal overlap extension of two actual intact linear duplexes.
//!
//! The new pairing is A.top / B.bottom after explicit denaturation/reannealing,
//! not a reference assembly or a ligation. Both paired 3′ ends extend on the
//! opposite strand through the existing five-prime-overhang fill-in kernel.
//! See Ho et al., Gene77 (1989), https://pubmed.ncbi.nlm.nih.gov/2744487/.
//!
//! This geometry has no intended-target argument and creates no private material
//! authority. The TrackedDuplex wrapper constructs the selected pair from
//! freshly resolved tracked sources: it copies only A.top and B.bottom with the private
//! `from_pieces` helper, unions their generation records, then uses tracked
//! `repair(FillFivePrimeOverhangs)`. Copied modification states stay on the active
//! source strands; newly synthesized spans are standard-unmodified. A subsequent
//! ordinary PCR must consume that private result and the supplied outer primers.
//! A serialized Geometry must never be used to restore that private state.
use crate::{
    Error, Topology,
    chemistry::{
        self, Duplex, EndRepairMethod, EndRepairReport, Strand, StrandMapping, Terminus, ThreePrime,
    },
};
use serde::{Deserialize, Serialize};
use std::io::Write;

pub const ENGINE_VERSION: &str = "overlap-extension-1";
/// Aggregate bases in the FOUR original stored strands, checked before cloning.
pub const MAX_SOURCE_BASES: usize = 8 << 20;
pub const MAX_COMPARISONS: usize = 50_000_000;
pub const MAX_CANDIDATES: usize = 1000;
const MAX_REPORT_BYTES: usize = 32 << 20;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Options {
    /// Selected suffix(A.top)/prefix(B.top) pairing length, not a target-derived join.
    pub overlap: usize,
    /// Explicit search/model threshold, not a temperature or efficacy prediction.
    pub min_paired_bases: usize,
    /// Includes validation of the selected pairing and alternative-terminal search.
    pub max_comparisons: usize,
    pub max_candidates: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            overlap: 20,
            min_paired_bases: 12,
            max_comparisons: 20_000_000,
            max_candidates: 100,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct TerminalPairings {
    /// Retained proper terminal overlaps in ascending length order. The explicitly
    /// selected, independently checked overlap is always included.
    pub overlaps: Vec<usize>,
    pub complete: bool,
    /// Unique only within the declared A.top/B.bottom terminal-pair model and
    /// minimum length, never a claim of reaction-mixture or amplification uniqueness.
    pub unique: bool,
    /// May remain present even if the candidate-retention budget stops the search.
    pub first_alternative: Option<usize>,
    pub comparisons: usize,
    pub limits: Vec<String>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SourcedMapping {
    /// 0=A, 1=B. Coordinates address that supplied operand's stored 5′→3′ strand.
    /// An operation wrapper must also retain any original operand orientation.
    pub source: usize,
    pub mapping: StrandMapping,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Geometry {
    pub engine_version: String,
    /// Full duplex GEOMETRY identities, including terminal groups and nicks.
    /// Positional modification identities belong to the tracked wrapper.
    pub source_geometry_revisions: [String; 2],
    pub options: Options,
    pub homology: String,
    /// Actual selected A.top/B.bottom pairing, with both 5′ overhangs retained.
    pub annealed: Duplex,
    /// Fresh native fill-in result. `fill.source` is exactly `annealed`.
    pub fill: EndRepairReport,
    /// Maps actual source strands to the final stored product strands. Synthesized
    /// entries identify the complementary template; they are not copied molecules.
    pub mappings: Vec<SourcedMapping>,
    pub terminal_pairings: TerminalPairings,
    pub warnings: Vec<String>,
}
impl Geometry {
    /// Validate bounded documentary evidence, never source identity or material
    /// authority. The original four source strands are not present here; only
    /// fresh execution against actual tracked sources can establish their origin.
    pub fn validate_documentary(&self) -> Result<(), Error> {
        let o = &self.options;
        let p = &self.terminal_pairings;
        // Check all collection/string bounds before validation or native repair
        // can allocate. The streaming serialized-size check allocates no buffer.
        if self.mappings.len() > 8
            || self.fill.mappings.len() > 8
            || p.overlaps.len() > MAX_CANDIDATES
            || self.warnings.len() > 100
            || self.fill.warnings.len() > 100
            || p.limits.len() > 100
            || self.homology.len() > MAX_SOURCE_BASES / 2
            || [
                &self.fill.top_added_3prime,
                &self.fill.bottom_added_3prime,
                &self.fill.top_removed_3prime,
                &self.fill.bottom_removed_3prime,
            ]
            .iter()
            .any(|s| s.len() > MAX_SOURCE_BASES)
            || self
                .warnings
                .iter()
                .chain(&self.fill.warnings)
                .chain(&p.limits)
                .any(|s| s.len() > 4096)
        {
            return Err(Error::TooLarge);
        }
        for duplex in [&self.annealed, &self.fill.source, &self.fill.duplex] {
            if duplex.top.len() > MAX_SOURCE_BASES
                || duplex.bottom.len() > MAX_SOURCE_BASES
                || !duplex.nicks.is_empty()
            {
                return Err(bad("invalid documentary duplex bounds or nicks"));
            }
        }
        if self.engine_version != ENGINE_VERSION
            || self.source_geometry_revisions.iter().any(|revision| {
                revision.len() != 64
                    || !revision
                        .bytes()
                        .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            })
            || o.min_paired_bases == 0
            || o.overlap < o.min_paired_bases
            || o.max_comparisons == 0
            || o.max_comparisons > MAX_COMPARISONS
            || o.max_candidates == 0
            || o.max_candidates > MAX_CANDIDATES
        {
            return Err(bad("invalid documentary version, revisions or options"));
        }
        serde_json::to_writer(Count(0), self).map_err(|_| Error::TooLarge)?;
        let a = self.annealed.top.len();
        let b = self.annealed.bottom.len();
        let k = o.overlap;
        if k >= a.min(b)
            || a.checked_add(b)
                .and_then(|n| n.checked_mul(2))
                .is_none_or(|n| n > MAX_SOURCE_BASES)
            || self.annealed.topology != Topology::Linear
            || self.annealed.bottom_start != (a - k) as i64
            || self.homology.len() != k
            || self.homology.as_bytes() != &self.annealed.top.as_str().as_bytes()[a - k..]
            || self.annealed.left.three_prime != ThreePrime::Hydroxyl
            || self.annealed.right.three_prime != ThreePrime::Hydroxyl
            || self.fill.source != self.annealed
        {
            return Err(bad(
                "documentary pairing differs from selected strand geometry",
            ));
        }
        self.annealed.validate()?;
        self.fill.source.validate()?;
        self.fill.duplex.validate()?;
        // Native duplex validation establishes the homology's antiparallel
        // complement without inventing either discarded original source strand.
        let replay =
            chemistry::repair_ends(&self.annealed, EndRepairMethod::FillFivePrimeOverhangs)?;
        if replay != self.fill || self.mappings.len() != self.fill.mappings.len() {
            return Err(bad("documentary fill differs from native repair"));
        }
        for (sourced, mapping) in self.mappings.iter().zip(&self.fill.mappings) {
            let source_len = match mapping.source_strand {
                Strand::Top => a,
                Strand::Bottom => b,
            };
            let product_len = match mapping.product_strand {
                Strand::Top => self.fill.duplex.top.len(),
                Strand::Bottom => self.fill.duplex.bottom.len(),
            };
            if sourced.source != usize::from(mapping.source_strand == Strand::Bottom)
                || sourced.mapping != *mapping
                || mapping.source_start >= mapping.source_end
                || mapping.source_end > source_len
                || mapping.product_start >= mapping.product_end
                || mapping.product_end > product_len
                || mapping.source_end - mapping.source_start
                    != mapping.product_end - mapping.product_start
            {
                return Err(bad("invalid documentary sourced strand mapping"));
            }
        }
        let proper = |overlap: usize| overlap >= o.min_paired_bases && overlap < a.min(b);
        let retained_alternative = p.overlaps.iter().any(|overlap| *overlap != k);
        if p.overlaps.is_empty()
            || p.overlaps.len() > o.max_candidates
            || !p.overlaps.contains(&k)
            || p.overlaps.iter().any(|overlap| !proper(*overlap))
            || p.overlaps.windows(2).any(|pair| pair[0] >= pair[1])
            || p.comparisons < k
            || p.comparisons > o.max_comparisons
            || p.complete != p.limits.is_empty()
            || p.unique != (p.complete && p.first_alternative.is_none())
            || retained_alternative && p.first_alternative.is_none()
            || p.first_alternative.is_some_and(|alternative| {
                alternative == k
                    || !proper(alternative)
                    || p.complete && !p.overlaps.contains(&alternative)
                    || p.overlaps
                        .iter()
                        .any(|retained| *retained != k && *retained > alternative)
            })
        {
            return Err(bad("inconsistent documentary terminal-pair evidence"));
        }
        Ok(())
    }
}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("overlap extension: {message}"))
}
struct Count(usize);
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= MAX_REPORT_BYTES)
            .ok_or_else(|| std::io::Error::other("overlap-extension report limit"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
struct Hash(blake3::Hasher);
impl Write for Hash {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0.update(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn revision(source: &Duplex) -> Result<String, Error> {
    let mut hash = Hash(blake3::Hasher::new());
    serde_json::to_writer(&mut hash, source).map_err(|_| Error::TooLarge)?;
    Ok(hash.0.finalize().to_hex().to_string())
}
fn shape(left: &Duplex, right: &Duplex, options: &Options) -> Result<usize, Error> {
    if options.min_paired_bases == 0
        || options.max_comparisons == 0
        || options.max_comparisons > MAX_COMPARISONS
        || options.max_candidates == 0
        || options.max_candidates > MAX_CANDIDATES
    {
        return Err(bad("invalid pairing threshold or bounded search settings"));
    }
    let mut bases = 0usize;
    for source in [left, right] {
        bases = bases
            .checked_add(source.top.len())
            .and_then(|n| n.checked_add(source.bottom.len()))
            .filter(|n| *n <= MAX_SOURCE_BASES)
            .ok_or(Error::TooLarge)?;
        if source.nicks.len() > 10_000 {
            return Err(Error::TooLarge);
        }
        source.validate()?;
        if source.topology != Topology::Linear
            || source.bottom_start != 0
            || source.top.len() != source.bottom.len()
            || !source.nicks.is_empty()
        {
            return Err(bad(
                "both operands must be intact, blunt, nick-free linear duplexes",
            ));
        }
    }
    let a = left.top.len();
    let b = right.top.len();
    if options.overlap < options.min_paired_bases
        || options.overlap == 0
        || options.overlap >= a.min(b)
    {
        return Err(bad(
            "choose a proper nonempty terminal overlap meeting the declared minimum; both strands must have a templated extension",
        ));
    }
    if left.right.three_prime != ThreePrime::Hydroxyl
        || right.left.three_prime != ThreePrime::Hydroxyl
    {
        return Err(bad(
            "the selected A.top and B.bottom extending 3′ ends must have known accessible hydroxyl groups",
        ));
    }
    a.checked_add(b)
        .and_then(|n| n.checked_sub(options.overlap))
        .ok_or(Error::TooLarge)
}
struct Comparisons {
    used: usize,
    maximum: usize,
}
impl Comparisons {
    fn equal(&mut self, a: u8, b: u8) -> Option<bool> {
        if self.used == self.maximum {
            return None;
        }
        self.used += 1;
        Some(a == b)
    }
}
/// Independent KMP border computation and one linear scan. No product candidates
/// are materialized during terminal search, and no sequence alignment is invoked.
fn terminal_pairings(
    left: &[u8],
    right: &[u8],
    options: &Options,
    budget: &mut Comparisons,
) -> TerminalPairings {
    let mut report = TerminalPairings {
        overlaps: vec![options.overlap],
        complete: false,
        unique: false,
        first_alternative: None,
        comparisons: budget.used,
        limits: vec![],
    };
    let length = left.len().min(right.len()) - 1;
    let minimal = length.saturating_sub(1).saturating_add(left.len());
    if budget.maximum - budget.used < minimal {
        report.limits.push(
            "The remaining comparison budget cannot complete even the minimum terminal-pair scan."
                .into(),
        );
        return report;
    }
    let pattern = &right[..length];
    let mut prefix = vec![0usize; length];
    for i in 1..length {
        let mut paired = prefix[i - 1];
        loop {
            match budget.equal(pattern[i], pattern[paired]) {
                Some(true) => {
                    paired += 1;
                    break;
                }
                Some(false) if paired > 0 => paired = prefix[paired - 1],
                Some(false) => break,
                None => {
                    report.comparisons = budget.used;
                    report.limits.push(
                        "Comparison budget exhausted while indexing terminal homology.".into(),
                    );
                    return report;
                }
            }
        }
        prefix[i] = paired;
    }
    let mut paired = 0usize;
    for &base in left {
        if paired == length {
            paired = prefix[paired - 1];
        }
        loop {
            match budget.equal(base, pattern[paired]) {
                Some(true) => {
                    paired += 1;
                    break;
                }
                Some(false) if paired > 0 => paired = prefix[paired - 1],
                Some(false) => break,
                None => {
                    report.comparisons = budget.used;
                    report.limits.push(
                        "Comparison budget exhausted while searching terminal pairings.".into(),
                    );
                    return report;
                }
            }
        }
    }
    report.complete = true;
    while paired >= options.min_paired_bases {
        if paired != options.overlap {
            report.first_alternative.get_or_insert(paired);
            if report.overlaps.len() == options.max_candidates {
                report.complete = false;
                report.limits.push("Terminal-pair candidate retention budget exhausted; further pairings remain unretained.".into());
                break;
            }
            report.overlaps.push(paired);
        }
        paired = prefix[paired - 1];
    }
    report.overlaps.sort_unstable();
    report.unique = report.complete && report.first_alternative.is_none();
    report.comparisons = budget.used;
    report
}

/// Reconstruct one explicitly selected reciprocal extension candidate from actual
/// stored strands. Incomplete/ambiguous terminal search is retained as diagnostic
/// evidence; a caller must not infer a uniquely determined material product from it.
/// Self-pairing, the other denatured strand combination and reaction kinetics are
/// outside this two-operand terminal-pair model, even when `unique` is true.
pub fn geometry(left: &Duplex, right: &Duplex, options: &Options) -> Result<Geometry, Error> {
    let n = shape(left, right, options)?;
    let a = left.top.len();
    let b = right.top.len();
    let k = options.overlap;
    let mut budget = Comparisons {
        used: 0,
        maximum: options.max_comparisons,
    };
    for (&l, &r) in left.top.as_str().as_bytes()[a - k..]
        .iter()
        .zip(right.top.as_str().as_bytes())
    {
        match budget.equal(l, r) {
            Some(true) => {}
            Some(false) => {
                return Err(bad(
                    "the selected A suffix and B prefix do not pair exactly",
                ));
            }
            None => {
                return Err(bad(
                    "the selected terminal pairing cannot be established within the comparison budget",
                ));
            }
        }
    }
    let pairings = terminal_pairings(
        left.top.as_str().as_bytes(),
        right.top.as_str().as_bytes(),
        options,
        &mut budget,
    );
    let annealed = Duplex {
        top: left.top.clone(),
        bottom: right.bottom.clone(),
        bottom_start: (a - k) as i64,
        topology: Topology::Linear,
        left: Terminus {
            five_prime: left.left.five_prime,
            three_prime: right.left.three_prime,
        },
        right: Terminus {
            five_prime: right.right.five_prime,
            three_prime: left.right.three_prime,
        },
        nicks: vec![],
    };
    annealed.validate()?;
    let fill = chemistry::repair_ends(&annealed, EndRepairMethod::FillFivePrimeOverhangs)?;
    if fill.duplex.top.len() != n
        || fill.duplex.bottom.len() != n
        || fill.duplex.bottom_start != 0
        || fill.top_added_3prime.len() != b - k
        || fill.bottom_added_3prime.len() != a - k
        || !fill.top_removed_3prime.is_empty()
        || !fill.bottom_removed_3prime.is_empty()
    {
        return Err(bad(
            "native fill-in does not match the selected reciprocal extension geometry",
        ));
    }
    let mappings = fill
        .mappings
        .iter()
        .map(|mapping| SourcedMapping {
            source: usize::from(mapping.source_strand == Strand::Bottom),
            mapping: mapping.clone(),
        })
        .collect();
    let mut warnings=vec![
        "One explicitly denatured A.top/B.bottom terminal pairing is extended under the standard-unmodified synthesis assumption. The other source strands and reaction-mixture mass balance are not modeled.".into(),
        "Terminal-search uniqueness concerns only this orientation and the recorded minimum pairing length. Self-association, other strand combinations, mismatched annealing, stability, concentrations, cycle kinetics and experimental success are not established.".into(),
        "Extension preserves the copied source strands and their 5′ groups; it requires no 5′ phosphorylation and performs no ligation or nick sealing. Supply an explicit outer primer pair for a subsequent fresh PCR.".into(),
    ];
    warnings.extend(pairings.limits.iter().cloned());
    if !pairings.complete || !pairings.unique {
        warnings.push("The selected geometry is diagnostic: incomplete or competing terminal pairing evidence cannot establish a uniquely determined material product.".into());
    }
    let report = Geometry {
        engine_version: ENGINE_VERSION.into(),
        source_geometry_revisions: [revision(left)?, revision(right)?],
        options: options.clone(),
        homology: left.top.as_str()[a - k..].into(),
        annealed,
        fill,
        mappings,
        terminal_pairings: pairings,
        warnings,
    };
    serde_json::to_writer(Count(0), &report).map_err(|_| Error::TooLarge)?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Molecule, Sequence,
        analysis::{self, Primer},
        chemistry::{FivePrime, Nick},
    };
    const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
    const LEFT: &str = "AGTCCGATGCTACGTATTTTCAGACC";
    const RIGHT: &str = "GCTACGTATTTTCAGACCTGTAGATCGTCA";
    const EXPECTED: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
    fn fixture() -> (Duplex, Duplex, Options) {
        // U(16) + oldGGAC(4) + V(20), independently declared GGAC→TT.
        // Both physical source fragments come from real primer reconstruction;
        // none of the kernel inputs is an intended-target molecule.
        let source = Molecule::dna("Unchanged original", ORIGINAL, Topology::Linear).unwrap();
        let before = serde_json::to_value(&source).unwrap();
        let left = analysis::simulate_pcr_geometry(
            &source,
            &Primer::new("Left outer", "AGTCCGAT", "").unwrap(),
            &Primer::new("Left mutagenic", "TACGTAGC", "GGTCTGAAAA").unwrap(),
        )
        .unwrap();
        let right = analysis::simulate_pcr_geometry(
            &source,
            &Primer::new("Right mutagenic", "TTCAGACC", "GCTACGTATT").unwrap(),
            &Primer::new("Right outer", "TGACGATC", "").unwrap(),
        )
        .unwrap();
        assert_eq!(left.len(), 1);
        assert_eq!(right.len(), 1);
        assert_eq!(left[0].product.sequence.as_str(), LEFT);
        assert_eq!(right[0].product.sequence.as_str(), RIGHT);
        assert_eq!(serde_json::to_value(&source).unwrap(), before);
        (
            Duplex::exact_pcr(&left[0].product.sequence).unwrap(),
            Duplex::exact_pcr(&right[0].product.sequence).unwrap(),
            Options {
                overlap: 18,
                ..Default::default()
            },
        )
    }
    #[test]
    fn actual_two_pcr_products_extend_reciprocally_to_independently_expected_38_base_insert() {
        let (left, right, options) = fixture();
        assert_eq!(left.top.len(), 26);
        assert_eq!(right.top.len(), 30);
        let before = serde_json::to_value((&left, &right)).unwrap();
        let result = geometry(&left, &right, &options).unwrap();
        assert_eq!(result.homology, "GCTACGTATTTTCAGACC");
        assert_eq!(result.annealed.top, left.top);
        assert_eq!(result.annealed.bottom, right.bottom);
        assert_eq!(result.annealed.bottom_start, 8);
        assert_eq!(result.fill.source, result.annealed);
        assert_eq!(result.fill.duplex.top.as_str(), EXPECTED);
        assert_eq!(result.fill.duplex.top.len(), 38);
        assert_eq!(
            result.fill.duplex.bottom,
            Sequence::dna(EXPECTED)
                .unwrap()
                .reverse_complement()
                .unwrap()
        );
        assert_eq!(result.fill.duplex.bottom_start, 0);
        assert_eq!(result.fill.top_added_3prime, "TGTAGATCGTCA");
        assert_eq!(result.fill.bottom_added_3prime, "ATCGGACT");
        assert_eq!(result.fill.top_added_3prime.len(), 12);
        assert_eq!(result.fill.bottom_added_3prime.len(), 8);
        assert!(result.fill.duplex.nicks.is_empty());
        assert!(result.terminal_pairings.complete);
        assert!(result.terminal_pairings.unique);
        assert_eq!(result.terminal_pairings.overlaps, vec![18]);
        assert_eq!(result.mappings.len(), 4);
        let copy_top = result
            .mappings
            .iter()
            .find(|m| !m.mapping.synthesized && m.mapping.product_strand == Strand::Top)
            .unwrap();
        assert_eq!(copy_top.source, 0);
        assert_eq!(
            (
                copy_top.mapping.source_start,
                copy_top.mapping.source_end,
                copy_top.mapping.product_start,
                copy_top.mapping.product_end
            ),
            (0, 26, 0, 26)
        );
        let copy_bottom = result
            .mappings
            .iter()
            .find(|m| !m.mapping.synthesized && m.mapping.product_strand == Strand::Bottom)
            .unwrap();
        assert_eq!(copy_bottom.source, 1);
        assert_eq!(
            (
                copy_bottom.mapping.source_start,
                copy_bottom.mapping.source_end,
                copy_bottom.mapping.product_start,
                copy_bottom.mapping.product_end
            ),
            (0, 30, 0, 30)
        );
        for m in result.mappings.iter().filter(|m| m.mapping.synthesized) {
            assert!(m.mapping.reverse);
            let source = if m.source == 0 { &left } else { &right };
            let template = match m.mapping.source_strand {
                Strand::Top => &source.top,
                Strand::Bottom => &source.bottom,
            };
            let expected =
                Sequence::dna(&template.as_str()[m.mapping.source_start..m.mapping.source_end])
                    .unwrap()
                    .reverse_complement()
                    .unwrap();
            let product = match m.mapping.product_strand {
                Strand::Top => &result.fill.duplex.top,
                Strand::Bottom => &result.fill.duplex.bottom,
            };
            assert_eq!(
                expected.as_str(),
                &product.as_str()[m.mapping.product_start..m.mapping.product_end]
            );
        }
        // Supplied outer primers amplify actual filled DNA, never the independent target.
        let actual = result
            .fill
            .duplex
            .molecule("Actual extended duplex")
            .unwrap();
        let outer = analysis::simulate_pcr_geometry(
            &actual,
            &Primer::new("Outer forward", "AGTCCGAT", "").unwrap(),
            &Primer::new("Outer reverse", "TGACGATC", "").unwrap(),
        )
        .unwrap();
        assert_eq!(outer.len(), 1);
        assert_eq!(outer[0].product.sequence.as_str(), EXPECTED);
        assert_eq!(serde_json::to_value((&left, &right)).unwrap(), before);
        // Operand orientation is explicit: reversing the assembly means reversing
        // and exchanging the actual sources, not reorienting reference strings.
        let reversed = geometry(
            &right.reverse().unwrap(),
            &left.reverse().unwrap(),
            &options,
        )
        .unwrap();
        assert_eq!(reversed.fill.duplex, result.fill.duplex.reverse().unwrap());
        assert_eq!(reversed.fill.top_added_3prime.len(), 8);
        assert_eq!(reversed.fill.bottom_added_3prime.len(), 12);
    }
    #[test]
    fn selected_strand_terminal_groups_are_preserved_without_ligation_or_unused_strand_repair() {
        let (mut left, mut right, options) = fixture();
        left.left.five_prime = FivePrime::Phosphate;
        right.right.five_prime = FivePrime::Unknown;
        // These belong to the OTHER, unused source strands and must not be
        // silently repaired or mistaken for the paired extending3′ endpoints.
        left.left.three_prime = ThreePrime::Blocked;
        right.right.three_prime = ThreePrime::Blocked;
        let before = serde_json::to_value((&left, &right)).unwrap();
        let result = geometry(&left, &right, &options).unwrap();
        assert_eq!(result.fill.duplex.left.five_prime, FivePrime::Phosphate);
        assert_eq!(result.fill.duplex.right.five_prime, FivePrime::Unknown);
        assert_eq!(result.fill.duplex.left.three_prime, ThreePrime::Hydroxyl);
        assert_eq!(result.fill.duplex.right.three_prime, ThreePrime::Hydroxyl);
        assert_eq!(serde_json::to_value((&left, &right)).unwrap(), before);
    }
    #[test]
    fn unsupported_termini_nicks_topologies_staggers_mismatches_and_ambiguity_are_rejected() {
        let (left, right, options) = fixture();
        for blocked in [ThreePrime::Blocked, ThreePrime::Unknown] {
            let mut changed = left.clone();
            changed.right.three_prime = blocked;
            assert!(geometry(&changed, &right, &options).is_err());
            let mut changed = right.clone();
            changed.left.three_prime = blocked;
            assert!(geometry(&left, &changed, &options).is_err());
        }
        let mut nicked = left.clone();
        nicked.nicks.push(Nick {
            strand: Strand::Top,
            position: 3,
        });
        assert!(geometry(&nicked, &right, &options).is_err());
        let mut circular = left.clone();
        circular.topology = Topology::Circular;
        assert!(geometry(&circular, &right, &options).is_err());
        let mut staggered = right.clone();
        staggered.bottom = Sequence::dna(&right.top.as_str()[..right.top.len() - 1])
            .unwrap()
            .reverse_complement()
            .unwrap();
        staggered.validate().unwrap();
        assert!(geometry(&left, &staggered, &options).is_err());
        let mut changed = LEFT.to_owned();
        changed.pop();
        changed.push('A');
        let wrong = Duplex::exact_pcr(&Sequence::dna(&changed).unwrap()).unwrap();
        assert!(geometry(&wrong, &right, &options).is_err());
        let mut ambiguous = left.clone();
        ambiguous.top = Sequence::dna(&format!("N{}", &LEFT[1..])).unwrap();
        ambiguous.bottom = ambiguous.top.reverse_complement().unwrap();
        assert!(geometry(&ambiguous, &right, &options).is_err());
        for overlap in [0, left.top.len().min(right.top.len())] {
            let mut bad = options.clone();
            bad.overlap = overlap;
            assert!(geometry(&left, &right, &bad).is_err());
        }
        let mut bad = options.clone();
        bad.min_paired_bases = options.overlap + 1;
        assert!(geometry(&left, &right, &bad).is_err());
    }
    #[test]
    fn competing_terminal_pairings_and_search_limits_do_not_claim_unique_material() {
        let left = Duplex::exact_pcr(&Sequence::dna("AAAAAA").unwrap()).unwrap();
        let right = left.clone();
        let options = Options {
            overlap: 2,
            min_paired_bases: 2,
            ..Default::default()
        };
        let full = geometry(&left, &right, &options).unwrap();
        assert!(full.terminal_pairings.complete);
        assert!(!full.terminal_pairings.unique);
        assert_eq!(full.terminal_pairings.overlaps, vec![2, 3, 4, 5]);
        assert!(full.terminal_pairings.first_alternative.is_some());
        assert_eq!(full.fill.duplex.top.as_str(), "AAAAAAAAAA");
        let mut limited = options.clone();
        limited.max_candidates = 1;
        let retained = geometry(&left, &right, &limited).unwrap();
        assert!(!retained.terminal_pairings.complete);
        assert!(!retained.terminal_pairings.unique);
        assert_eq!(retained.terminal_pairings.overlaps, vec![2]);
        assert!(retained.terminal_pairings.first_alternative.is_some());
        assert!(!retained.terminal_pairings.limits.is_empty());
        limited = options.clone();
        limited.max_comparisons = options.overlap;
        let unsearched = geometry(&left, &right, &limited).unwrap();
        assert!(!unsearched.terminal_pairings.complete);
        assert!(!unsearched.terminal_pairings.unique);
        assert_eq!(unsearched.terminal_pairings.comparisons, options.overlap);
        assert_eq!(unsearched.terminal_pairings.overlaps, vec![2]);
        limited.max_comparisons = options.overlap - 1;
        assert!(geometry(&left, &right, &limited).is_err());
        for (comparisons, candidates) in [
            (0, 10),
            (MAX_COMPARISONS + 1, 10),
            (100, 0),
            (100, MAX_CANDIDATES + 1),
        ] {
            let mut bad = options.clone();
            bad.max_comparisons = comparisons;
            bad.max_candidates = candidates;
            assert!(geometry(&left, &right, &bad).is_err());
        }
    }
    #[test]
    fn complete_linear_terminal_search_matches_independent_short_word_enumeration() {
        fn word(bits: usize) -> String {
            (0..4)
                .map(|i| if bits & (1 << i) == 0 { 'A' } else { 'C' })
                .collect()
        }
        for l in 0..16 {
            for r in 0..16 {
                let left_text = word(l);
                let right_text = word(r);
                let expected = (1..4)
                    .filter(|&k| left_text[4 - k..] == right_text[..k])
                    .collect::<Vec<_>>();
                let Some(&overlap) = expected.first() else {
                    continue;
                };
                let left = Duplex::exact_pcr(&Sequence::dna(&left_text).unwrap()).unwrap();
                let right = Duplex::exact_pcr(&Sequence::dna(&right_text).unwrap()).unwrap();
                let options = Options {
                    overlap,
                    min_paired_bases: 1,
                    ..Default::default()
                };
                let result = geometry(&left, &right, &options).unwrap();
                assert!(result.terminal_pairings.complete);
                assert_eq!(
                    result.terminal_pairings.overlaps, expected,
                    "{left_text}/{right_text}"
                );
                assert_eq!(result.terminal_pairings.unique, expected.len() == 1);
                assert_eq!(
                    result.fill.duplex.top.as_str(),
                    format!("{left_text}{}", &right_text[overlap..])
                );
            }
        }
    }
}
