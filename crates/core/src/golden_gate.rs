//! Independent finite Golden Gate reaction geometry. No external software or kinetics.
use crate::{
    Error, Molecule, Topology,
    analysis::{self, Enzyme, Match},
    chemistry::{self, Duplex, End, FivePrime, Polarity, ThreePrime},
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub const MODEL: &str = "complete-digest/one-copy-fragment-pool-1";
const MAX_POOL: usize = 128;
const MAX_REPORT_BASES: usize = 16 * 1024 * 1024;
const MAX_RECONSTRUCTED_BASES: usize = 100_000_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum GoldenGateEnzyme {
    BsaI,
    BsmBI,
    BbsI,
    SapI,
}
impl GoldenGateEnzyme {
    /// Factual recognition/cut offsets, measured from the recognition's 5′ boundary.
    pub fn profile(self) -> Enzyme {
        let (name, recognition, top_cut, bottom_cut) = match self {
            Self::BsaI => ("BsaI", "GGTCTC", 7, 11),
            Self::BsmBI => ("BsmBI", "CGTCTC", 7, 11),
            Self::BbsI => ("BbsI", "GAAGAC", 8, 12),
            Self::SapI => ("SapI", "GCTCTTC", 8, 11),
        };
        Enzyme {
            name: name.into(),
            recognition: recognition.into(),
            top_cut,
            bottom_cut,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct GoldenGateOptions {
    pub enzyme: GoldenGateEnzyme,
    pub circular: bool,
    pub both_orientations: bool,
    /// Include one-fragment self-ligation by default so background is not silently hidden.
    pub min_fragments: usize,
    pub max_fragments: usize,
    pub max_products: usize,
    pub max_states: usize,
}
impl Default for GoldenGateOptions {
    fn default() -> Self {
        Self {
            enzyme: GoldenGateEnzyme::BsaI,
            circular: true,
            both_orientations: true,
            min_fragments: 1,
            max_fragments: 16,
            max_products: 128,
            max_states: 50_000,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GoldenGateSource {
    pub source_index: usize,
    pub name: String,
    pub revision: String,
    pub recognition_sites: Vec<Match>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GoldenGateFragment {
    pub source_index: usize,
    pub source_revision: String,
    /// Index from complete digest, in ascending source top-cut order.
    pub fragment_index: usize,
    pub top_start: usize,
    /// Unwrapped source-reference boundary, not relative to top_start.
    pub bottom_start: i64,
    pub duplex: Duplex,
    pub left_end: Option<End>,
    pub right_end: Option<End>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct OrientedFragment {
    /// Index into GoldenGateReport.fragments, not the original input list.
    pub fragment: usize,
    pub reverse: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecutSite {
    pub recognition_start: usize,
    pub reverse: bool,
    /// Circular boundaries are modulo product length; linear boundaries are signed/unclipped.
    pub top_cut: i64,
    pub bottom_cut: i64,
    pub fully_double_stranded: bool,
    pub cleavable_under_model: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GoldenGateCandidate {
    pub path: Vec<OrientedFragment>,
    pub duplex: Duplex,
    pub recut_sites: Vec<RecutSite>,
    /// No complete recognition motif in the top reference, on either orientation.
    /// This does not establish yield, absence of star activity or kinetic stability.
    pub restriction_site_free: bool,
    /// Circular top rotations collapse here, but reverse complements remain separate.
    pub directed_sequence_key: String,
    /// Circle rotation+reverse-complement equivalence, or actual oriented linear duplex geometry.
    /// It does not establish a unique annotation/source history.
    pub physical_sequence_key: String,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GoldenGateReport {
    pub model: String,
    pub options: GoldenGateOptions,
    pub sources: Vec<GoldenGateSource>,
    pub fragments: Vec<GoldenGateFragment>,
    /// Intact uncut circular input indices; not silently claimed to be newly assembled products.
    pub carryover_sources: Vec<usize>,
    pub candidates: Vec<GoldenGateCandidate>,
    /// Complete only for the stated finite one-copy pool/orientation options, not the real mixture.
    pub enumeration_complete: bool,
    pub exhaustive_mixture_search: bool,
    pub states_examined: usize,
    pub truncation_reasons: Vec<String>,
    pub stable_physical_product_count: usize,
    pub warnings: Vec<String>,
}

fn bad(message: &str) -> Error {
    Error::Unsupported(message.into())
}
fn sealable(end: &End) -> bool {
    end.chemistry.five_prime == FivePrime::Phosphate
        && end.chemistry.three_prime == ThreePrime::Hydroxyl
}
fn complementary(a: &str, b: &str) -> bool {
    a.len() == b.len()
        && a.bytes()
            .zip(b.bytes().rev())
            .all(|(x, y)| x == crate::model::complement(y, crate::Alphabet::Dna))
}
fn joins(right: &End, left: &End) -> bool {
    right.polarity == left.polarity
        && right.polarity != Polarity::Blunt
        && complementary(&right.overhang, &left.overhang)
        && sealable(right)
        && sealable(left)
}

/// Inspect recognition sites after the actual junction sequence has been reconstructed.
/// A near-linear-end motif remains visible even when there is insufficient duplex for cleavage.
pub fn assess_product(duplex: &Duplex, enzyme: GoldenGateEnzyme) -> Result<Vec<RecutSite>, Error> {
    duplex.validate()?;
    let profile = enzyme.profile();
    let molecule = Molecule::dna(
        "Golden Gate recut reference",
        duplex.top.as_str(),
        duplex.topology,
    )?;
    let n = duplex.top.len() as i64;
    let paired_start = duplex.bottom_start.max(0);
    let bottom_end = duplex.bottom_start + duplex.bottom.len() as i64;
    let paired_end = bottom_end.min(n);
    let mut sites = Vec::new();
    for hit in analysis::find_matches(&molecule, &profile.recognition, true)? {
        let start = hit.start as i64;
        let m = profile.recognition.len() as i64;
        let (mut top, mut bottom) = if hit.reverse {
            (
                start + m - i64::from(profile.bottom_cut),
                start + m - i64::from(profile.top_cut),
            )
        } else {
            (
                start + i64::from(profile.top_cut),
                start + i64::from(profile.bottom_cut),
            )
        };
        let circular = duplex.topology == Topology::Circular;
        let paired = circular || start >= paired_start && start + m <= paired_end;
        let cleavable = paired
            && if circular {
                (bottom - top).unsigned_abs() < (n as u64)
            } else {
                top > 0 && top < n && bottom > duplex.bottom_start && bottom < bottom_end
            };
        if circular {
            top = top.rem_euclid(n);
            bottom = bottom.rem_euclid(n);
        }
        sites.push(RecutSite {
            recognition_start: hit.start,
            reverse: hit.reverse,
            top_cut: top,
            bottom_cut: bottom,
            fully_double_stranded: paired,
            cleavable_under_model: cleavable,
        });
    }
    Ok(sites)
}

// Booth-style two-candidate elimination: least cyclic rotation in linear work.
fn least_rotation(bytes: &[u8]) -> usize {
    let n = bytes.len();
    let (mut a, mut b, mut k) = (0usize, 1usize, 0usize);
    while a < n && b < n && k < n {
        let (x, y) = (bytes[(a + k) % n], bytes[(b + k) % n]);
        if x == y {
            k += 1;
            continue;
        }
        if x > y {
            a += k + 1;
            if a <= b {
                a = b + 1;
            }
        } else {
            b += k + 1;
            if b <= a {
                b = a + 1;
            }
        }
        k = 0;
    }
    a.min(b).min(n.saturating_sub(1))
}
fn rotation_bytes(sequence: &str) -> Vec<u8> {
    let index = least_rotation(sequence.as_bytes());
    sequence.as_bytes()[index..]
        .iter()
        .chain(&sequence.as_bytes()[..index])
        .copied()
        .collect()
}
fn sequence_keys(duplex: &Duplex) -> Result<(String, String), Error> {
    let directed = if duplex.topology == Topology::Circular {
        rotation_bytes(duplex.top.as_str())
    } else {
        duplex.top.as_str().as_bytes().to_vec()
    };
    let mut directed_hasher = blake3::Hasher::new();
    directed_hasher.update(if duplex.topology == Topology::Circular {
        b"circular"
    } else {
        b"linear"
    });
    directed_hasher.update(&directed);
    let physical = if duplex.topology == Topology::Circular {
        let rc = duplex.top.reverse_complement()?;
        directed.clone().min(rotation_bytes(rc.as_str()))
    } else {
        let a = serde_json::to_vec(duplex).map_err(|e| Error::Format(e.to_string()))?;
        let b = serde_json::to_vec(&duplex.reverse()?).map_err(|e| Error::Format(e.to_string()))?;
        a.min(b)
    };
    Ok((
        directed_hasher.finalize().to_hex().to_string(),
        blake3::hash(&physical).to_hex().to_string(),
    ))
}

#[derive(Clone)]
struct Oriented {
    choice: OrientedFragment,
    duplex: Duplex,
    left: End,
    right: End,
}
struct Search<'a> {
    report: &'a mut GoldenGateReport,
    oriented: Vec<Oriented>,
    edges: Vec<Vec<usize>>,
    path: Vec<usize>,
    used: u128,
    stored_bases: usize,
    reconstructed_bases: usize,
    stopped: bool,
}
impl Search<'_> {
    fn truncate(&mut self, reason: &str, stop: bool) {
        self.report.enumeration_complete = false;
        if !self.report.truncation_reasons.iter().any(|r| r == reason) {
            self.report.truncation_reasons.push(reason.into());
        }
        self.stopped |= stop;
    }
    fn emit(&mut self) -> Result<(), Error> {
        if self.report.candidates.len() >= self.report.options.max_products {
            self.truncate(
                "product output limit reached before completing finite pool enumeration",
                true,
            );
            return Ok(());
        }
        let bases = self
            .path
            .iter()
            .try_fold(0usize, |n, i| {
                n.checked_add(
                    self.oriented[*i].duplex.top.len() + self.oriented[*i].duplex.bottom.len(),
                )
            })
            .ok_or(Error::TooLarge)?;
        self.reconstructed_bases = self
            .reconstructed_bases
            .checked_add(bases)
            .ok_or(Error::TooLarge)?;
        if self.reconstructed_bases > MAX_RECONSTRUCTED_BASES {
            self.truncate("reconstruction work exceeds 100 million strand bases", true);
            return Ok(());
        }
        if self.stored_bases.saturating_add(bases) > MAX_REPORT_BASES {
            self.truncate(
                "stored pool/products exceed the 16 MiB strand-base budget",
                true,
            );
            return Ok(());
        }
        let parts = self
            .path
            .iter()
            .map(|i| self.oriented[*i].duplex.clone())
            .collect::<Vec<_>>();
        let product = chemistry::ligate(&parts, self.report.options.circular, false)?;
        let recut_sites = assess_product(&product, self.report.options.enzyme)?;
        let (directed_sequence_key, physical_sequence_key) = sequence_keys(&product)?;
        self.stored_bases += bases;
        self.report.candidates.push(GoldenGateCandidate {
            path: self.path.iter().map(|i| self.oriented[*i].choice).collect(),
            duplex: product,
            restriction_site_free: recut_sites.is_empty(),
            recut_sites,
            directed_sequence_key,
            physical_sequence_key,
        });
        Ok(())
    }
    fn visit(&mut self, index: usize) -> Result<(), Error> {
        if self.stopped {
            return Ok(());
        }
        if self.report.states_examined >= self.report.options.max_states {
            self.truncate("assembly search-state limit reached", true);
            return Ok(());
        }
        self.report.states_examined += 1;
        let bit = 1u128 << self.oriented[index].choice.fragment;
        self.used |= bit;
        self.path.push(index);
        let first = self.path[0];
        let closed = joins(&self.oriented[index].right, &self.oriented[first].left);
        if self.path.len() >= self.report.options.min_fragments
            && (!self.report.options.circular || closed)
        {
            self.emit()?;
        }
        if !self.stopped {
            let edges = self.edges[index].clone();
            for next in edges {
                let fragment = self.oriented[next].choice.fragment;
                if self.used & (1u128 << fragment) != 0 {
                    continue;
                }
                // A circular cycle is represented once per physical orientation, with its
                // lowest pool index at the origin. No sequence/origin history is discarded.
                if self.report.options.circular && fragment < self.oriented[first].choice.fragment {
                    continue;
                }
                if self.path.len() >= self.report.options.max_fragments {
                    self.truncate(
                        "fragment depth limit omitted an extendable assembly path",
                        false,
                    );
                    continue;
                }
                self.visit(next)?;
                if self.stopped {
                    break;
                }
            }
        }
        self.path.pop();
        self.used &= !bit;
        Ok(())
    }
}

/// Complete digest followed by bounded exact ligation of one copy of each pool element.
/// This is not a kinetic reaction simulator and cannot prove absence of concatemers,
/// partial-digest products, star activity, carryover, host effects or repeated-copy products.
pub fn simulate(
    inputs: &[Molecule],
    options: &GoldenGateOptions,
) -> Result<GoldenGateReport, Error> {
    simulate_model(inputs, None, options)
}

/// Physical graph integration: supply actual reconstructed substrates, never chemistry
/// inferred from imported Molecule metadata. Source top sequence/topology must agree.
pub fn simulate_reconstructed(
    inputs: &[Molecule],
    substrates: &[Duplex],
    options: &GoldenGateOptions,
) -> Result<GoldenGateReport, Error> {
    if inputs.len() != substrates.len() {
        return Err(Error::Format(
            "Golden Gate substrate/source counts differ".into(),
        ));
    }
    simulate_model(inputs, Some(substrates), options)
}
fn simulate_model(
    inputs: &[Molecule],
    substrates: Option<&[Duplex]>,
    options: &GoldenGateOptions,
) -> Result<GoldenGateReport, Error> {
    if inputs.is_empty()
        || inputs.len() > 32
        || options.min_fragments == 0
        || options.min_fragments > options.max_fragments
        || options.max_fragments > 32
        || options.max_products == 0
        || options.max_products > 1000
        || options.max_states == 0
        || options.max_states > 1_000_000
    {
        return Err(bad(
            "Golden Gate needs 1–32 inputs, fragment bounds 1–32, products 1–1000 and search states 1–1,000,000",
        ));
    }
    let mut report=GoldenGateReport{model:MODEL.into(),options:options.clone(),sources:vec![],fragments:vec![],carryover_sources:vec![],candidates:vec![],enumeration_complete:true,exhaustive_mixture_search:false,states_examined:0,truncation_reasons:vec![],stable_physical_product_count:0,warnings:vec![
        "Complete cleavage and sealed, exactly complementary ligation are assumed. Methylation, supplier conditions, fidelity, kinetics, partial digestion and experimental success are not evaluated.".into(),
        "Enumeration uses at most one copy of each complete-digest pool fragment. It cannot establish absence of repeated-copy concatemers or other real-mixture products; uncut circular carryover is listed separately.".into(),
        "Restriction-site-free means absence of the selected complete recognition motif in the reconstructed reference. It is not an empirical stability or yield prediction. Recuttable products remain visible.".into(),
        "Product paths retain source/fragment orientation and coordinates. Annotation projection must use these immutable source records; the duplex alone does not reconstruct annotation history.".into(),
    ]};
    if !options.both_orientations {
        report.warnings.push(
            "Reverse fragment choices are excluded by the requested orientation scope.".into(),
        );
    }
    let profile = options.enzyme.profile();
    let mut stored = 0usize;
    for (source_index, input) in inputs.iter().enumerate() {
        let duplex = if let Some(substrates) = substrates {
            input.validate()?;
            let d = &substrates[source_index];
            d.validate()?;
            if d.top != input.sequence || d.topology != input.topology {
                return Err(Error::Format(
                    "Golden Gate source differs from reconstructed physical substrate".into(),
                ));
            }
            d.clone()
        } else {
            Duplex::from_molecule(input)?
        };
        let recognition_sites = analysis::find_matches(input, &profile.recognition, true)?;
        let revision = input.revision();
        if duplex.topology == Topology::Circular && recognition_sites.is_empty() {
            report.carryover_sources.push(source_index);
        }
        report.sources.push(GoldenGateSource {
            source_index,
            name: input.name.chars().take(1024).collect(),
            revision: revision.clone(),
            recognition_sites,
        });
        let fragments = chemistry::digest(&duplex, std::slice::from_ref(&profile))?;
        if report.fragments.len() + fragments.len() > MAX_POOL {
            return Err(bad(
                "Golden Gate complete digest exceeds the 128-fragment pool limit",
            ));
        }
        for (fragment_index, fragment) in fragments.into_iter().enumerate() {
            stored = stored
                .checked_add(fragment.duplex.top.len() + fragment.duplex.bottom.len())
                .ok_or(Error::TooLarge)?;
            if stored > MAX_REPORT_BASES {
                return Err(Error::TooLarge);
            }
            let ends = if fragment.duplex.topology == Topology::Linear {
                Some(fragment.duplex.ends()?)
            } else {
                None
            };
            report.fragments.push(GoldenGateFragment {
                source_index,
                source_revision: revision.clone(),
                fragment_index,
                top_start: fragment.top_start,
                bottom_start: fragment.bottom_start,
                duplex: fragment.duplex,
                left_end: ends.as_ref().map(|e| e.0.clone()),
                right_end: ends.map(|e| e.1),
            });
        }
    }
    if !report.carryover_sources.is_empty() {
        report.warnings.push("Uncut circular inputs remain potential carryover and are not silently counted as newly assembled candidates.".into());
    }
    let mut oriented = Vec::new();
    for (fragment, record) in report.fragments.iter().enumerate() {
        if record.duplex.topology != Topology::Linear {
            continue;
        }
        for reverse in [false, true] {
            if reverse && !options.both_orientations {
                continue;
            }
            let duplex = if reverse {
                record.duplex.reverse()?
            } else {
                record.duplex.clone()
            };
            let (left, right) = duplex.ends()?;
            if options.circular && (!sealable(&left) || !sealable(&right)) {
                continue;
            }
            oriented.push(Oriented {
                choice: OrientedFragment { fragment, reverse },
                duplex,
                left,
                right,
            });
        }
    }
    let edges = oriented
        .iter()
        .map(|a| {
            oriented
                .iter()
                .enumerate()
                .filter_map(|(i, b)| {
                    (a.choice.fragment != b.choice.fragment && joins(&a.right, &b.left))
                        .then_some(i)
                })
                .collect::<Vec<_>>()
        })
        .collect();
    let mut search = Search {
        report: &mut report,
        oriented,
        edges,
        path: vec![],
        used: 0,
        stored_bases: stored,
        reconstructed_bases: 0,
        stopped: false,
    };
    for start in 0..search.oriented.len() {
        search.visit(start)?;
        if search.stopped {
            break;
        }
    }
    report.stable_physical_product_count = report
        .candidates
        .iter()
        .filter(|c| c.restriction_site_free)
        .map(|c| &c.physical_sequence_key)
        .collect::<BTreeSet<_>>()
        .len();
    if report.candidates.iter().any(|c| !c.restriction_site_free) {
        report.warnings.push("Some ligatable candidates retain or recreate a recognition site, including sites formed across junctions; those candidates are not restriction-site-free products.".into());
    }
    Ok(report)
}
