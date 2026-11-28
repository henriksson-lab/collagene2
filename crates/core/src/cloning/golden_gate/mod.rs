//! Independent finite Golden Gate reaction geometry. No external software or kinetics.
use crate::{
    chemistry::{
        self, Duplex, End, FivePrime, Polarity, ThreePrime,
        methylation::Outcome,
        tracked_duplex::{
            self as modified_duplex, CheckedCleavageEvidence, CheckedDigestOptions, CleavageBudget,
            TrackedDuplex,
        },
    },
    sequence::analysis::{Enzyme, Match},
};
use collagene_model::{Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub const MODEL: &str = "positional-checked-digest/one-copy-fragment-pool-1";
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
    pub outcome: Outcome,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GoldenGateCandidate {
    pub path: Vec<OrientedFragment>,
    pub duplex: Duplex,
    pub recut_sites: Vec<RecutSite>,
    pub cleavage: CheckedCleavageEvidence,
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
    pub profile_revision: String,
    pub source_cleavage: Vec<CheckedCleavageEvidence>,
    pub products_determined: bool,
    /// Current enclosing-plan ledger, including preparative checked cuts and
    /// Golden Gate source/recut scans. Nested historical origin replay uses
    /// separate bounded ledgers under the shared native ancestry limits.
    pub comparisons: usize,
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
            .all(|(x, y)| x == collagene_model::complement(y, collagene_model::Alphabet::Dna))
}
fn joins(right: &End, left: &End) -> bool {
    right.polarity == left.polarity
        && right.polarity != Polarity::Blunt
        && complementary(&right.overhang, &left.overhang)
        && sealable(right)
        && sealable(left)
}

/// Actual positional recut assessment; source chemistry is never inferred from
/// a serialized endpoint summary. Incomplete scans do not prove site absence.
pub fn assess_product(
    source: &TrackedDuplex,
    enzyme: GoldenGateEnzyme,
    budget: &mut CleavageBudget,
) -> Result<(Vec<RecutSite>, CheckedCleavageEvidence), Error> {
    let selection = modified_duplex::checked_selection(&enzyme.profile().name)?;
    let checked = modified_duplex::checked_digest_with_budget(
        source,
        &CheckedDigestOptions {
            selection,
            max_comparisons: budget.remaining().max(1),
            max_sites: 4096,
        },
        budget,
    )?;
    let profile = enzyme.profile();
    let duplex = &source.duplex;
    let n = duplex.top.len() as i64;
    let paired_start = duplex.bottom_start.max(0);
    let bottom_end = duplex.bottom_start + duplex.bottom.len() as i64;
    let paired_end = bottom_end.min(n);
    let mut sites = Vec::new();
    for hit in &checked.sites {
        let start = hit.recognition_start as i64;
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
        let geometric = paired
            && if circular {
                (bottom - top).unsigned_abs() < n as u64
            } else {
                top > 0 && top < n && bottom > duplex.bottom_start && bottom < bottom_end
            };
        if circular {
            top = top.rem_euclid(n);
            bottom = bottom.rem_euclid(n)
        }
        sites.push(RecutSite {
            recognition_start: hit.recognition_start,
            reverse: hit.reverse,
            top_cut: top,
            bottom_cut: bottom,
            fully_double_stranded: paired,
            cleavable_under_model: geometric && hit.outcome == Outcome::PermittedUnderProfile,
            outcome: hit.outcome,
        });
    }
    Ok((sites, checked.evidence()))
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
pub(crate) fn sequence_keys(duplex: &Duplex) -> Result<(String, String), Error> {
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
    tracked_fragments: Vec<TrackedDuplex>,
    selected_index: Option<usize>,
    selected: Option<TrackedDuplex>,
    budget: &'a mut CleavageBudget,
    report_bytes: usize,
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
            .map(|i| {
                let choice = self.oriented[*i].choice;
                let actual = &self.tracked_fragments[choice.fragment];
                if choice.reverse {
                    actual.reverse()
                } else {
                    Ok(actual.clone())
                }
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let product = modified_duplex::ligate(&parts, self.report.options.circular, false)?;
        let (recut_sites, cleavage) =
            assess_product(&product, self.report.options.enzyme, self.budget)?;
        if !cleavage.complete {
            self.truncate("Candidate recut assessment incomplete under shared cleavage limits; site absence is not established",false);
        }
        if cleavage
            .sites
            .iter()
            .any(|s| matches!(s.outcome, Outcome::Unknown | Outcome::Impaired))
        {
            let warning = "A retained candidate has unresolved positional recutting evidence; its native assembly exists conditionally, not as an isolated or methylation-cleared cycling product.";
            if !self.report.warnings.iter().any(|w| w == warning) {
                self.report.warnings.push(warning.into())
            }
        }
        let (directed_sequence_key, physical_sequence_key) = sequence_keys(&product.duplex)?;
        let candidate = GoldenGateCandidate {
            path: self.path.iter().map(|i| self.oriented[*i].choice).collect(),
            duplex: product.duplex.clone(),
            restriction_site_free: cleavage.complete && recut_sites.is_empty(),
            recut_sites,
            cleavage,
            directed_sequence_key,
            physical_sequence_key,
        };
        self.report_bytes = self
            .report_bytes
            .checked_add(
                serde_json::to_vec(&candidate)
                    .map_err(|_| Error::TooLarge)?
                    .len(),
            )
            .filter(|n| *n <= 32 * 1024 * 1024 - 65536)
            .ok_or(Error::TooLarge)?;
        if self.selected_index == Some(self.report.candidates.len()) {
            self.selected = Some(product)
        }
        self.stored_bases += bases;
        self.report.candidates.push(candidate);
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

/// Original references remain chemically unknown; explicitly generated/saved
/// substrates enter the same kernel through simulate_tracked.
pub fn simulate(
    inputs: &[Molecule],
    options: &GoldenGateOptions,
) -> Result<GoldenGateReport, Error> {
    if inputs.is_empty() || inputs.len() > 32 {
        return Err(Error::TooLarge);
    }
    collagene_model::bounds::model_preflight(&inputs.iter().collect::<Vec<_>>())?;
    let substrates = inputs
        .iter()
        .map(|m| TrackedDuplex::from_unknown_duplex(&Duplex::from_molecule(m)?))
        .collect::<Result<Vec<_>, Error>>()?;
    let mut budget = CleavageBudget::new(100_000_000)?;
    Ok(simulate_tracked(inputs, &substrates, options, None, &mut budget)?.report)
}
/// Private native selection from the same checked fragment/candidate pass.
/// The report never restores material and no desired target participates.
#[derive(Debug)]
pub struct FreshGoldenGate {
    pub report: GoldenGateReport,
    pub selected: Option<TrackedDuplex>,
    _private: (),
}
pub fn simulate_tracked(
    inputs: &[Molecule],
    substrates: &[TrackedDuplex],
    options: &GoldenGateOptions,
    selected_index: Option<usize>,
    budget: &mut CleavageBudget,
) -> Result<FreshGoldenGate, Error> {
    if inputs.len() != substrates.len() {
        return Err(bad("Golden Gate source/tracked substrate counts differ"));
    }
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
    collagene_model::bounds::model_preflight(&inputs.iter().collect::<Vec<_>>())?;
    crate::cloning::assembly::recipe::hash(&(inputs, substrates))?;
    let mut report=GoldenGateReport{model:MODEL.into(),options:options.clone(),profile_revision:modified_duplex::profile_revision(),source_cleavage:vec![],products_determined:true,comparisons:budget.used(),sources:vec![],fragments:vec![],carryover_sources:vec![],candidates:vec![],enumeration_complete:true,exhaustive_mixture_search:false,states_examined:0,truncation_reasons:vec![],stable_physical_product_count:0,warnings:vec![
        "Positional checked cleavage uses actual inherited/synthesized modification states under pinned current profiles. Sealed complementary ligation, supplier conditions, fidelity and kinetics remain separate assumptions; this is not experimental success.".into(),
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
    let mut report_storage = serde_json::to_vec(&report)
        .map_err(|_| Error::TooLarge)?
        .len();
    let mut tracked_fragments = Vec::new();
    let selection = modified_duplex::checked_selection(&profile.name)?;
    for (source_index, input) in inputs.iter().enumerate() {
        input.validate()?;
        let tracked = &substrates[source_index];
        tracked.validate()?;
        let duplex = &tracked.duplex;
        if duplex.top != input.sequence || duplex.topology != input.topology {
            return Err(bad(
                "Golden Gate source differs from actual tracked substrate",
            ));
        }
        let checked = modified_duplex::checked_digest_with_budget(
            tracked,
            &CheckedDigestOptions {
                selection: selection.clone(),
                max_comparisons: budget.remaining().max(1),
                max_sites: 4096,
            },
            budget,
        )?;
        let recognition_sites = checked
            .sites
            .iter()
            .map(|site| Match {
                start: site.recognition_start,
                length: profile.recognition.len(),
                reverse: site.reverse,
                certain: true,
            })
            .collect();
        let revision = input.revision();
        let source_record = GoldenGateSource {
            source_index,
            name: input.name.chars().take(1024).collect(),
            revision: revision.clone(),
            recognition_sites,
        };
        let evidence = checked.evidence();
        report_storage = report_storage
            .checked_add(
                serde_json::to_vec(&source_record)
                    .map_err(|_| Error::TooLarge)?
                    .len(),
            )
            .and_then(|n| n.checked_add(serde_json::to_vec(&evidence).ok()?.len()))
            .filter(|n| *n <= 32 * 1024 * 1024 - 65536)
            .ok_or(Error::TooLarge)?;
        report.sources.push(source_record);

        if !checked.complete {
            report.enumeration_complete = false;
            report.truncation_reasons.push(format!(
                "Source {source_index} checked cleavage scan incomplete"
            ))
        }
        if !checked.products_determined {
            report.products_determined = false
        }
        if checked.products_determined
            && checked.fragments.len() == 1
            && checked.fragments[0].duplex.duplex.topology == Topology::Circular
        {
            report.carryover_sources.push(source_index)
        }
        report.source_cleavage.push(evidence);
        let fragments = checked.fragments;
        if report.fragments.len() + fragments.len() > MAX_POOL {
            return Err(bad(
                "Golden Gate complete digest exceeds the 128-fragment pool limit",
            ));
        }
        for (fragment_index, fragment) in fragments.into_iter().enumerate() {
            let actual = fragment.duplex;
            let fragment = chemistry::DigestFragment {
                duplex: actual.duplex.clone(),
                top_start: fragment.top_start,
                bottom_start: fragment.bottom_start,
            };
            tracked_fragments.push(actual);
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
            let fragment_record = GoldenGateFragment {
                source_index,
                source_revision: revision.clone(),
                fragment_index,
                top_start: fragment.top_start,
                bottom_start: fragment.bottom_start,
                duplex: fragment.duplex,
                left_end: ends.as_ref().map(|e| e.0.clone()),
                right_end: ends.map(|e| e.1),
            };
            report_storage = report_storage
                .checked_add(
                    serde_json::to_vec(&fragment_record)
                        .map_err(|_| Error::TooLarge)?
                        .len(),
                )
                .filter(|n| *n <= 32 * 1024 * 1024 - 65536)
                .ok_or(Error::TooLarge)?;
            report.fragments.push(fragment_record);
        }
    }
    report.comparisons = budget.used();
    if !report.products_determined {
        report.fragments.clear();
        // Retain independently established intact circular carryover identities;
        // no fragment subset is enumerated as the requested complete pool.
        report.enumeration_complete = false;
        report.truncation_reasons.push("Requested pool withheld: one or more complete source fragment sets are unresolved; no available subset was substituted".into());
        crate::cloning::assembly::recipe::hash(&report)?;
        return Ok(FreshGoldenGate {
            report,
            selected: None,
            _private: (),
        });
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
        tracked_fragments,
        selected_index,
        selected: None,
        budget,
        report_bytes: 0,
        oriented,
        edges,
        path: vec![],
        used: 0,
        stored_bases: stored,
        reconstructed_bases: 0,
        stopped: false,
    };
    search.report_bytes = serde_json::to_vec(&search.report)
        .map_err(|_| Error::TooLarge)?
        .len();
    for start in 0..search.oriented.len() {
        search.visit(start)?;
        if search.stopped {
            break;
        }
    }
    let selected = search.selected.take();
    drop(search);
    report.comparisons = budget.used();
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
    crate::cloning::assembly::recipe::hash(&report)?;
    Ok(FreshGoldenGate {
        report,
        selected,
        _private: (),
    })
}

pub mod domestication;
pub mod fusion_design;
pub mod fusion_rank;
pub mod junction_edit;
pub mod prepared;
pub mod primers;
pub mod reaction;
