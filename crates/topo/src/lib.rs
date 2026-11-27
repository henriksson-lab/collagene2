#![doc = include_str!("../README.md")]
//! Explicit elementary TOPO strand transfer. Activation is an input chemical
//! assumption, never inferred from a sequence record or supplier metadata.
use collagene_model::{Alphabet, Sequence};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
pub const ENGINE: &str = "topo-explicit-strand-transfer-1";
pub const PRIMARY_MODEL: &str = "vaccinia-tyr274-explicit-graph-ta-cacc-1";
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum FivePrime {
    Hydroxyl,
    Phosphate,
    Blocked,
    Unknown,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Enzyme {
    VacciniaTypeIBTyr274,
    Uncharacterized,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ThreePrime {
    Hydroxyl,
    Phosphate,
    Blocked,
    Unknown,
    Phosphotyrosyl {
        enzyme: Enzyme,
        activation_id: String,
    },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Strand {
    pub role: String,
    pub id: String,
    pub sequence: Sequence,
    pub five_prime: FivePrime,
    pub three_prime: ThreePrime,
}
/// Both intervals use actual 5′→3′ source indices. a[start+k] pairs with
/// b[start+length-1-k]; the strands are antiparallel.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PairSpan {
    pub strand_a: usize,
    pub start_a: usize,
    pub strand_b: usize,
    pub start_b: usize,
    pub length: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Substrates {
    pub strands: Vec<Strand>,
    pub pairing: Vec<PairSpan>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Geometry {
    /// Religation across an actual nick, templated by a continuous opposite backbone.
    NickReseal,
    Flush,
    ThreePrimeTA,
    DirectionalCacc,
}
/// Original source strands only. This consumes the donor's actual terminal3′
/// activation and acceptor's actual terminal5′OH; it is not an end search.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Transfer {
    pub donor: usize,
    pub acceptor: usize,
    pub geometry: Geometry,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Limits {
    pub max_source_bases: usize,
    pub max_strands: usize,
    pub max_pair_spans: usize,
    pub max_steps: usize,
    pub max_work_cells: usize,
    pub max_work_bytes: usize,
    pub max_output_bases: usize,
    pub max_report_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_source_bases: 100_000,
            max_strands: 32,
            max_pair_spans: 1024,
            max_steps: 16,
            max_work_cells: 16_000_000,
            max_work_bytes: 32 * 1024 * 1024,
            max_output_bases: 400_000,
            max_report_bytes: 16 * 1024 * 1024,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Budget {
    SourceBases,
    Strands,
    PairSpans,
    Steps,
    WorkCells,
    WorkBytes,
    OutputBases,
    ReportBytes,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    Invalid(String),
    Limit {
        kind: Budget,
        required: usize,
        limit: usize,
    },
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "invalid TOPO substrate/model: {s}"),
            Self::Limit {
                kind,
                required,
                limit,
            } => write!(f, "{kind:?} budget exceeded: {required} > {limit}"),
        }
    }
}
impl std::error::Error for Error {}
fn bad(s: &str) -> Error {
    Error::Invalid(s.into())
}
fn bounded(n: usize, max: usize, kind: Budget) -> Result<(), Error> {
    if n > max {
        Err(Error::Limit {
            kind,
            required: n,
            limit: max,
        })
    } else {
        Ok(())
    }
}
fn text(s: &str, n: usize) -> bool {
    s.len() <= n && !s.trim().is_empty()
}
fn add(a: usize, b: usize) -> Result<usize, Error> {
    a.checked_add(b).ok_or_else(|| bad("integer overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize, Error> {
    a.checked_mul(b).ok_or_else(|| bad("integer overflow"))
}
fn comp(b: u8) -> Option<u8> {
    match b {
        b'A' => Some(b'T'),
        b'T' => Some(b'A'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        _ => None,
    }
}
fn mask(b: u8) -> u8 {
    match b {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
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
fn complementary_masks(a: u8, b: u8) -> bool {
    let a = mask(a);
    let c = ((a & 1) << 3) | ((a & 8) >> 3) | ((a & 2) << 1) | ((a & 4) >> 1);
    c & mask(b) != 0
}
fn activation(t: &ThreePrime) -> Option<&str> {
    if let ThreePrime::Phosphotyrosyl { activation_id, .. } = t {
        Some(activation_id)
    } else {
        None
    }
}
fn limits(l: &Limits) -> Result<(), Error> {
    for (v, max, kind) in [
        (l.max_source_bases, 1_000_000, Budget::SourceBases),
        (l.max_strands, 64, Budget::Strands),
        (l.max_pair_spans, 4096, Budget::PairSpans),
        (l.max_steps, 64, Budget::Steps),
        (l.max_work_cells, 256_000_000, Budget::WorkCells),
        (l.max_work_bytes, 128 * 1024 * 1024, Budget::WorkBytes),
        (l.max_output_bases, 4_000_000, Budget::OutputBases),
        (l.max_report_bytes, 32 * 1024 * 1024, Budget::ReportBytes),
    ] {
        if v == 0 {
            return Err(bad("zero resource bound"));
        }
        bounded(v, max, kind)?
    }
    Ok(())
}
/// Borrowed length/string/span checks BEFORE hashing, cloning or allocating
/// per-base graph arrays. No recursive user structures are accepted.
pub fn preflight(s: &Substrates, steps: &[Transfer], l: &Limits) -> Result<usize, Error> {
    limits(l)?;
    bounded(s.strands.len(), l.max_strands, Budget::Strands)?;
    bounded(s.pairing.len(), l.max_pair_spans, Budget::PairSpans)?;
    bounded(steps.len(), l.max_steps, Budget::Steps)?;
    if s.strands.is_empty() {
        return Err(bad("no source strands"));
    }
    let mut ids = BTreeSet::new();
    let mut activations = BTreeSet::new();
    let mut total = 0;
    for strand in &s.strands {
        if !text(&strand.role, 256)
            || !text(&strand.id, 256)
            || !ids.insert((&strand.role, &strand.id))
        {
            return Err(bad("source strand identifiers must be bounded and unique"));
        }
        if strand.sequence.alphabet() != Alphabet::Dna {
            return Err(bad("actual DNA strands required"));
        }
        total = add(total, strand.sequence.len())?;
        bounded(total, l.max_source_bases, Budget::SourceBases)?;
        if let Some(id) = activation(&strand.three_prime)
            && (!text(id, 256) || !activations.insert(id))
        {
            return Err(bad("distinct declared activation identities required"));
        }
    }
    for p in &s.pairing {
        if p.strand_a >= s.strands.len()
            || p.strand_b >= s.strands.len()
            || p.strand_a == p.strand_b
            || p.length == 0
            || p.start_a
                .checked_add(p.length)
                .is_none_or(|n| n > s.strands[p.strand_a].sequence.len())
            || p.start_b
                .checked_add(p.length)
                .is_none_or(|n| n > s.strands[p.strand_b].sequence.len())
        {
            return Err(bad("pair span outside actual source strands"));
        }
    }
    for t in steps {
        if t.donor >= s.strands.len()
            || t.acceptor >= s.strands.len()
            || (t.donor == t.acceptor && t.geometry != Geometry::NickReseal)
        {
            return Err(bad("distinct donor/acceptor source strands required"));
        }
    }
    // Reserve graph construction, pairing checks, traversal, union/find and source
    // audit work before allocation. Explicit local geometry adds at most64 cells/step.
    bounded(
        add(mul(total, 64)?, mul(steps.len(), 64)?)?,
        l.max_work_cells,
        Budget::WorkCells,
    )?;
    bounded(
        add(
            add(mul(total, 128)?, mul(s.pairing.len(), 128)?)?,
            add(
                add(mul(s.strands.len(), 4096)?, mul(steps.len(), 2048)?)?,
                8192,
            )?,
        )?,
        l.max_work_bytes,
        Budget::WorkBytes,
    )?;
    bounded(total, l.max_output_bases, Budget::OutputBases)?;
    Ok(total)
}
fn hash<T: Serialize>(v: &T) -> Result<String, Error> {
    struct Sink(blake3::Hasher);
    impl std::io::Write for Sink {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.0.update(b);
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut h = Sink(blake3::Hasher::new());
    serde_json::to_writer(&mut h, v).map_err(|_| bad("hash serialization"))?;
    Ok(h.0.finalize().to_hex().to_string())
}
pub fn batch_revision(s: &Substrates, steps: &[Transfer], l: &Limits) -> Result<String, Error> {
    preflight(s, steps, l)?;
    hash(&(ENGINE, PRIMARY_MODEL, s, steps, l))
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, PartialOrd, Ord)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct Call {
    pub strand: usize,
    pub index: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CallPair {
    pub a: Call,
    pub b: Call,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Outcome {
    Applied,
    BlockedUnderDeclaredModel,
    Unknown,
    Unsupported,
}
/// Fate of the four old complementary calls during an explicitly modeled invasion.
/// The complete-displacement and unpaired-continuation cases are declared graph
/// extensions, not assertions of commercial-kit geometry or efficiency.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum DisplacementKind {
    DuplexAnchoredThreePrimeFlap,
    FreeDisplacedChain,
    UnpairedTailExtension,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Displacement {
    pub kind: DisplacementKind,
    pub calls_five_to_three: Vec<Call>,
    pub five_prime_continuation: Option<Call>,
    pub continuation_pair: Option<Call>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct TransferEvidence {
    pub step: usize,
    pub transfer: Transfer,
    pub outcome: Outcome,
    pub reason: String,
    pub donor_call: Call,
    pub acceptor_call: Call,
    pub removed_pairs: Vec<CallPair>,
    pub added_pairs: Vec<CallPair>,
    pub released_activation: Option<String>,
    pub displacement: Option<Displacement>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SourcePlacement {
    pub source_strand: usize,
    pub start: usize,
    pub end: usize,
    pub output_start: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ChainEnds {
    Linear {
        five_prime: FivePrime,
        three_prime: ThreePrime,
    },
    CovalentlyCircular,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Chain {
    pub sequence: Sequence,
    pub ends: ChainEnds,
    pub sources: Vec<SourcePlacement>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Nick {
    pub three: Call,
    pub five: Call,
    pub three_prime: ThreePrime,
    pub five_prime: FivePrime,
}
/// A source-call run not currently paired. This includes retained displaced
/// terminal flaps, unreacted overhangs and declared ssDNA; none are deleted.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct UnpairedRun {
    pub chain: usize,
    pub start: usize,
    pub end: usize,
    pub terminal_five: bool,
    pub terminal_three: bool,
    pub displaced_by: Vec<usize>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Component {
    pub chains: Vec<usize>,
    pub paired_nucleotides: usize,
    pub unpaired_nucleotides: usize,
    pub covalently_closed: bool,
}
/// Pairing and nick Call coordinates remain immutable SOURCE-strand indices.
/// Chain.sources provides the exact source→output-coordinate translation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProductGraph {
    pub chains: Vec<Chain>,
    pub source_pairing: Vec<PairSpan>,
    pub nicks: Vec<Nick>,
    pub unpaired: Vec<UnpairedRun>,
    pub components: Vec<Component>,
    pub released_activations: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine: String,
    pub model: String,
    pub batch_revision: String,
    pub limits: Limits,
    pub source_revisions: Vec<String>,
    pub steps: Vec<TransferEvidence>,
    pub product: ProductGraph,
    pub source_bases: usize,
    pub work_complete: bool,
    pub all_requested_transfers_applied: bool,
    pub unsupported_or_unknown: bool,
    pub fully_specified_dna: bool,
    pub assumptions: Vec<String>,
}
#[derive(Clone)]
struct Graph {
    concrete: bool,
    offsets: Vec<usize>,
    pairs: Vec<Option<usize>>,
    next: Vec<Option<usize>>,
    prev: Vec<Option<usize>>,
    displaced: Vec<Option<usize>>,
}
impl Graph {
    fn new(s: &Substrates, n: usize) -> Result<Self, Error> {
        let mut offsets = vec![];
        let mut cursor = 0;
        for st in &s.strands {
            offsets.push(cursor);
            cursor += st.sequence.len()
        }
        let mut g = Self {
            concrete: s
                .strands
                .iter()
                .all(|st| st.sequence.as_str().bytes().all(|b| comp(b).is_some())),
            offsets,
            pairs: vec![None; n],
            next: vec![None; n],
            prev: vec![None; n],
            displaced: vec![None; n],
        };
        for (i, st) in s.strands.iter().enumerate() {
            let a = g.offsets[i];
            for k in 1..st.sequence.len() {
                g.next[a + k - 1] = Some(a + k);
                g.prev[a + k] = Some(a + k - 1)
            }
        }
        for p in &s.pairing {
            for k in 0..p.length {
                let a = g.offsets[p.strand_a] + p.start_a + k;
                let b = g.offsets[p.strand_b] + p.start_b + p.length - 1 - k;
                if g.pairs[a].is_some() || g.pairs[b].is_some() {
                    return Err(bad("a source nucleotide has multiple pairing partners"));
                }
                let ca = g.base(s, a);
                let cb = g.base(s, b);
                if !complementary_masks(ca, cb) {
                    return Err(bad(
                        "declared base pairing has no compatible concrete calls",
                    ));
                }
                g.pairs[a] = Some(b);
                g.pairs[b] = Some(a)
            }
        }
        Ok(g)
    }
    fn call(&self, node: usize) -> Call {
        let strand = self.offsets.partition_point(|&s| s <= node) - 1;
        Call {
            strand,
            index: node - self.offsets[strand],
        }
    }
    fn base(&self, s: &Substrates, node: usize) -> u8 {
        let c = self.call(node);
        s.strands[c.strand].sequence.as_str().as_bytes()[c.index]
    }
    fn flat(&self, c: Call) -> usize {
        self.offsets[c.strand] + c.index
    }
    fn pair(&self, a: usize, b: usize) -> CallPair {
        CallPair {
            a: self.call(a),
            b: self.call(b),
        }
    }
    fn fully_concrete(&self, _s: &Substrates) -> bool {
        self.concrete
    }
    fn paired_run(&self, start: usize, len: usize) -> bool {
        if len == 0 {
            return false;
        }
        let Some(mut other) = self.pairs[start] else {
            return false;
        };
        for k in 0..len {
            if self.pairs[start + k] != Some(other) {
                return false;
            }
            if k + 1 < len {
                let Some(previous) = self.prev[other] else {
                    return false;
                };
                other = previous;
            }
        }
        true
    }
}
type PairChanges = (Vec<(usize, usize)>, Vec<(usize, usize)>);
fn eligibility(
    s: &Substrates,
    g: &Graph,
    t: &Transfer,
) -> Result<PairChanges, (Outcome, &'static str)> {
    let donor = &s.strands[t.donor];
    let acceptor = &s.strands[t.acceptor];
    let d = g.offsets[t.donor] + donor.sequence.len() - 1;
    let a = g.offsets[t.acceptor];
    if g.next[d].is_some() || g.prev[a].is_some() {
        return Err((
            Outcome::BlockedUnderDeclaredModel,
            "one requested terminal has already formed a backbone bond",
        ));
    }
    match &donor.three_prime {
        ThreePrime::Phosphotyrosyl {
            enzyme: Enzyme::VacciniaTypeIBTyr274,
            ..
        } => {}
        ThreePrime::Phosphotyrosyl { .. } => {
            return Err((
                Outcome::Unsupported,
                "uncharacterized phosphotyrosyl enzyme",
            ));
        }
        ThreePrime::Unknown => {
            return Err((Outcome::Unknown, "donor terminal activation is unknown"));
        }
        _ => {
            return Err((
                Outcome::BlockedUnderDeclaredModel,
                "no declared phosphotyrosyl donor activation",
            ));
        }
    }
    match acceptor.five_prime {
        FivePrime::Hydroxyl => {}
        FivePrime::Unknown => return Err((Outcome::Unknown, "acceptor5′ chemistry is unknown")),
        _ => {
            return Err((
                Outcome::BlockedUnderDeclaredModel,
                "acceptor is not a free5′OH; phosphatase/ligase not inferred",
            ));
        }
    }
    if !g.fully_concrete(s) {
        return Err((
            Outcome::Unknown,
            "ambiguous DNA/pairing is retained but transfer is unresolved under the concrete-call profile",
        ));
    }
    let ds = donor.sequence.as_str().as_bytes();
    if ds.len() < 5 || !matches!(&ds[ds.len() - 5..], b"CCCTT" | b"TCCTT") {
        return Err((
            Outcome::Unsupported,
            "terminal donor context is outside the supported Vaccinia(Y)CCTT profile",
        ));
    }
    match t.geometry {
        Geometry::NickReseal => {
            let Some(dp) = g.pairs[d] else {
                return Err((Outcome::Unsupported, "nick donor must be paired"));
            };
            let Some(ap) = g.pairs[a] else {
                return Err((Outcome::Unsupported, "nick acceptor must be paired"));
            };
            if !g.paired_run(d - 4, 5) || g.prev[dp] != Some(ap) || g.next[ap] != Some(dp) {
                return Err((
                    Outcome::Unsupported,
                    "nick transfer requires continuous antiparallel template across exactly adjacent termini",
                ));
            }
            Ok((vec![], vec![]))
        }
        Geometry::Flush => {
            let Some(dp) = g.pairs[d] else {
                return Err((Outcome::Unsupported, "flush donor end is not paired"));
            };
            let Some(ap) = g.pairs[a] else {
                return Err((Outcome::Unsupported, "flush acceptor end is not paired"));
            };
            if g.prev[dp].is_some() || g.next[ap].is_some() || !g.paired_run(d - 4, 5) {
                return Err((
                    Outcome::Unsupported,
                    "not two fully flush supported junction faces",
                ));
            }
            Ok((vec![], vec![]))
        }
        Geometry::ThreePrimeTA => {
            if donor.sequence.len() < 5
                || acceptor.sequence.is_empty()
                || g.pairs[d].is_some()
                || !g.paired_run(d - 4, 4)
            {
                return Err((
                    Outcome::Unsupported,
                    "TA donor requires exactly one unpaired3′T with paired upstream motif",
                ));
            }
            let Some(dp) = g.pairs[d - 1] else {
                return Err((Outcome::Unsupported, "TA donor duplex missing"));
            };
            let Some(ap) = g.pairs[a] else {
                return Err((Outcome::Unsupported, "TA acceptor5′ duplex missing"));
            };
            let Some(overhang) = g.next[ap] else {
                return Err((Outcome::Unsupported, "acceptor lacks a3′A overhang"));
            };
            if g.prev[dp].is_some()
                || g.next[overhang].is_some()
                || g.pairs[overhang].is_some()
                || g.base(s, d) != b'T'
                || g.base(s, overhang) != b'A'
            {
                return Err((
                    Outcome::Unsupported,
                    "outside exact complementary1ntTA geometry; mismatch transfer is not ruled out biologically",
                ));
            }
            Ok((vec![], vec![(d, overhang)]))
        }
        Geometry::DirectionalCacc => {
            let Some(dp) = g.pairs[d] else {
                return Err((
                    Outcome::Unsupported,
                    "directional donor3′ is not duplex paired",
                ));
            };
            let c = g.call(dp);
            if c.index != 4
                || g.prev[g.offsets[c.strand]].is_some()
                || !g.paired_run(d - 4, 5)
                || ds.len() < 5
                || acceptor.sequence.len() < 4
                || &acceptor.sequence.as_str().as_bytes()[..4] != b"CACC"
                || &s.strands[c.strand].sequence.as_str().as_bytes()[..4] != b"GGTG"
            {
                return Err((
                    Outcome::Unsupported,
                    "outside declared CACC/GGTG4nt invasion geometry",
                ));
            }
            let mut removed = vec![];
            let mut added = vec![];
            for k in 0..4 {
                let over = g.offsets[c.strand] + 3 - k;
                let ai = a + k;
                let Some(old) = g.pairs[ai] else {
                    return Err((
                        Outcome::Unsupported,
                        "incoming CACC must initially be a blunt duplex",
                    ));
                };
                if g.pairs[over].is_some()
                    || g.base(s, over) != comp(g.base(s, ai)).expect("concrete")
                {
                    return Err((
                        Outcome::Unsupported,
                        "vector overhang is not free/complementary",
                    ));
                }
                removed.push((ai, old));
                added.push((ai, over));
            }
            let old0 = removed[0].1;
            if g.next[old0].is_some()
                || !(1..4).all(|k| g.prev[removed[k - 1].1] == Some(removed[k].1))
            {
                return Err((
                    Outcome::Unsupported,
                    "displacement must retain an actual terminal3′4nt flap",
                ));
            }
            Ok((removed, added))
        }
    }
}
fn pairing_spans(g: &Graph) -> Vec<PairSpan> {
    let mut spans: Vec<PairSpan> = vec![];
    for (node, p) in g.pairs.iter().enumerate() {
        let Some(other) = *p else { continue };
        if node >= other {
            continue;
        }
        let a = g.call(node);
        let b = g.call(other);
        if let Some(last) = spans.last_mut()
            && last.strand_a == a.strand
            && last.strand_b == b.strand
            && last.start_a + last.length == a.index
            && b.index + 1 == last.start_b
        {
            last.start_b = b.index;
            last.length += 1;
            continue;
        }
        spans.push(PairSpan {
            strand_a: a.strand,
            start_a: a.index,
            strand_b: b.strand,
            start_b: b.index,
            length: 1,
        })
    }
    spans
}
fn terminal3(s: &Substrates, g: &Graph, node: usize) -> ThreePrime {
    let c = g.call(node);
    if c.index + 1 == s.strands[c.strand].sequence.len() {
        s.strands[c.strand].three_prime.clone()
    } else {
        ThreePrime::Unknown
    }
}
#[derive(Default)]
struct TerminalOverrides {
    five: BTreeMap<usize, FivePrime>,
    three: BTreeMap<usize, ThreePrime>,
}
fn effective_five(s: &Substrates, g: &Graph, node: usize, ends: &TerminalOverrides) -> FivePrime {
    ends.five.get(&node).copied().unwrap_or_else(|| {
        let c = g.call(node);
        if c.index == 0 {
            s.strands[c.strand].five_prime
        } else {
            FivePrime::Unknown
        }
    })
}
fn effective_three(s: &Substrates, g: &Graph, node: usize, ends: &TerminalOverrides) -> ThreePrime {
    ends.three
        .get(&node)
        .cloned()
        .unwrap_or_else(|| terminal3(s, g, node))
}
fn output(s: &Substrates, g: &Graph, released: Vec<String>) -> Result<ProductGraph, Error> {
    output_with_terminals(s, g, released, &TerminalOverrides::default())
}
fn output_with_terminals(
    s: &Substrates,
    g: &Graph,
    released: Vec<String>,
    ends: &TerminalOverrides,
) -> Result<ProductGraph, Error> {
    let n = g.pairs.len();
    let mut visited = vec![false; n];
    let mut chains = vec![];
    let mut chain_of = vec![0; n];
    let mut unpaired = vec![];
    let starts = (0..n)
        .filter(|&i| g.prev[i].is_none())
        .chain(0..n)
        .collect::<Vec<_>>();
    for start in starts {
        if visited[start] {
            continue;
        }
        let mut nodes = vec![];
        let mut current = start;
        let circular = g.prev[start].is_some();
        loop {
            if visited[current] {
                if current != start {
                    return Err(bad("backbone graph merges unexpectedly"));
                }
                break;
            }
            visited[current] = true;
            chain_of[current] = chains.len();
            nodes.push(current);
            match g.next[current] {
                Some(next) => current = next,
                None => break,
            }
        }
        let mut sequence = String::with_capacity(nodes.len());
        let mut sources: Vec<SourcePlacement> = vec![];
        for (index, &node) in nodes.iter().enumerate() {
            sequence.push(g.base(s, node) as char);
            let c = g.call(node);
            if let Some(last) = sources.last_mut()
                && last.source_strand == c.strand
                && last.end == c.index
            {
                last.end += 1;
                continue;
            }
            sources.push(SourcePlacement {
                source_strand: c.strand,
                start: c.index,
                end: c.index + 1,
                output_start: index,
            });
        }
        let mut i = 0;
        while i < nodes.len() {
            if g.pairs[nodes[i]].is_some() {
                i += 1;
                continue;
            }
            let first = i;
            let mut displaced = BTreeSet::new();
            while i < nodes.len() && g.pairs[nodes[i]].is_none() {
                if let Some(step) = g.displaced[nodes[i]] {
                    displaced.insert(step);
                }
                i += 1
            }
            unpaired.push(UnpairedRun {
                chain: chains.len(),
                start: first,
                end: i,
                terminal_five: !circular && first == 0,
                terminal_three: !circular && i == nodes.len(),
                displaced_by: displaced.into_iter().collect(),
            });
        }
        let ends = if circular {
            ChainEnds::CovalentlyCircular
        } else {
            ChainEnds::Linear {
                five_prime: effective_five(s, g, start, ends),
                three_prime: effective_three(s, g, *nodes.last().expect("nonempty"), ends),
            }
        };
        chains.push(Chain {
            sequence: Sequence::dna(&sequence).map_err(|e| Error::Invalid(e.to_string()))?,
            ends,
            sources,
        });
    }
    let mut nicks = vec![];
    for node in 0..n {
        if g.next[node].is_some() {
            continue;
        }
        let Some(partner) = g.pairs[node] else {
            continue;
        };
        let Some(previous) = g.prev[partner] else {
            continue;
        };
        let Some(five) = g.pairs[previous] else {
            continue;
        };
        if five == node || g.prev[five].is_some() {
            continue;
        }
        nicks.push(Nick {
            three: g.call(node),
            five: g.call(five),
            three_prime: effective_three(s, g, node, ends),
            five_prime: effective_five(s, g, five, ends),
        });
    }
    // Connected molecular complexes under actual covalent bonds AND declared base pairs.
    let mut parent = (0..chains.len()).collect::<Vec<_>>();
    fn root(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x]
        }
        x
    }
    for (a, b) in g.pairs.iter().enumerate() {
        if let Some(b) = *b {
            let x = root(&mut parent, chain_of[a]);
            let y = root(&mut parent, chain_of[b]);
            parent[x] = y;
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..chains.len() {
        let r = root(&mut parent, i);
        groups.entry(r).or_default().push(i)
    }
    let mut counts = BTreeMap::<usize, (usize, usize)>::new();
    for (node, &chain) in chain_of.iter().enumerate() {
        let r = root(&mut parent, chain);
        let entry = counts.entry(r).or_default();
        if g.pairs[node].is_some() {
            entry.0 += 1
        } else {
            entry.1 += 1
        }
    }
    let mut components = vec![];
    for (root_id, group) in groups {
        let (paired_bases, unpaired_bases) = counts.remove(&root_id).unwrap_or_default();
        let closed = group
            .iter()
            .all(|&i| matches!(chains[i].ends, ChainEnds::CovalentlyCircular));
        components.push(Component {
            chains: group,
            paired_nucleotides: paired_bases,
            unpaired_nucleotides: unpaired_bases,
            covalently_closed: closed,
        });
    }
    Ok(ProductGraph {
        chains,
        source_pairing: pairing_spans(g),
        nicks,
        unpaired,
        components,
        released_activations: released,
    })
}
fn displacement(g: &Graph, removed: &[(usize, usize)]) -> Option<Displacement> {
    let first = removed.last()?.1;
    let continuation = g.prev[first];
    let paired = continuation.and_then(|call| g.pairs[call]);
    let kind = match (continuation, paired) {
        (None, _) => DisplacementKind::FreeDisplacedChain,
        (Some(_), Some(_)) => DisplacementKind::DuplexAnchoredThreePrimeFlap,
        (Some(_), None) => DisplacementKind::UnpairedTailExtension,
    };
    Some(Displacement {
        kind,
        calls_five_to_three: removed
            .iter()
            .rev()
            .map(|&(_, call)| g.call(call))
            .collect(),
        five_prime_continuation: continuation.map(|call| g.call(call)),
        continuation_pair: paired.map(|call| g.call(call)),
    })
}
#[derive(Clone)]
struct ExecutionState {
    graph: Graph,
    evidence: Vec<TransferEvidence>,
    released: Vec<String>,
}
impl ExecutionState {
    fn new(s: &Substrates, n: usize) -> Result<Self, Error> {
        Ok(Self {
            graph: Graph::new(s, n)?,
            evidence: vec![],
            released: vec![],
        })
    }
    fn apply(&mut self, s: &Substrates, steps: &[Transfer]) {
        for (offset, t) in steps.iter().enumerate() {
            let step = self.evidence.len();
            let _ = offset;
            let donor_call = Call {
                strand: t.donor,
                index: s.strands[t.donor].sequence.len() - 1,
            };
            let acceptor_call = Call {
                strand: t.acceptor,
                index: 0,
            };
            let d = self.graph.flat(donor_call);
            let a = self.graph.flat(acceptor_call);
            let (outcome, reason, removed_pairs, added_pairs, released_activation, displacement) =
                match eligibility(s, &self.graph, t) {
                    Ok((remove, add)) => {
                        // Classify the old tail's 5′ continuation before changing any pair.
                        let displacement = displacement(&self.graph, &remove);
                        let rp = remove.iter().map(|&(x, y)| self.graph.pair(x, y)).collect();
                        let ap = add.iter().map(|&(x, y)| self.graph.pair(x, y)).collect();
                        for (x, y) in remove {
                            self.graph.pairs[x] = None;
                            self.graph.pairs[y] = None;
                            self.graph.displaced[y] = Some(step);
                        }
                        for (x, y) in add {
                            self.graph.pairs[x] = Some(y);
                            self.graph.pairs[y] = Some(x)
                        }
                        self.graph.next[d] = Some(a);
                        self.graph.prev[a] = Some(d);
                        let id = activation(&s.strands[t.donor].three_prime)
                            .expect("eligible activation")
                            .to_owned();
                        self.released.push(id.clone());
                        (
                            Outcome::Applied,
                            "one phosphotyrosyl→5′OH transfer applied; no other bond/repair inferred",
                            rp,
                            ap,
                            Some(id),
                            displacement,
                        )
                    }
                    Err((outcome, reason)) => (outcome, reason, vec![], vec![], None, None),
                };
            self.evidence.push(TransferEvidence {
                step,
                transfer: t.clone(),
                outcome,
                reason: reason.into(),
                donor_call,
                acceptor_call,
                removed_pairs,
                added_pairs,
                released_activation,
                displacement,
            });
        }
    }
    fn finish(
        self,
        s: &Substrates,
        steps: &[Transfer],
        l: &Limits,
        n: usize,
    ) -> Result<Report, Error> {
        let batch_revision = batch_revision(s, steps, l)?;
        let product = output(s, &self.graph, self.released)?;
        bounded(
            product.source_pairing.len(),
            l.max_pair_spans,
            Budget::PairSpans,
        )?;
        let report=Report{engine:ENGINE.into(),model:PRIMARY_MODEL.into(),batch_revision,limits:l.clone(),source_revisions:s.strands.iter().map(hash).collect::<Result<_,_>>()?,all_requested_transfers_applied:self.evidence.iter().all(|e|e.outcome==Outcome::Applied),unsupported_or_unknown:self.evidence.iter().any(|e|matches!(e.outcome,Outcome::Unsupported|Outcome::Unknown)),fully_specified_dna:self.graph.fully_concrete(s),steps:self.evidence,product,source_bases:n,work_complete:true,assumptions:vec!["Source strands, pairing and terminal activation are explicitly declared chemical hypotheses; supplier names/reference sequences cannot establish activation.".into(),"Only specified elementary transfers are considered; no kinetic mixture completeness, efficacy, host repair, ligase, polymerase or proprietary kit parity is established.".into(),"All source nucleotides, displaced calls, unreacted chains, existing nicks and retained activation are preserved. Unmodeled chemical modifications remain outside this terminal/backbone model.".into(),"Exact TA and CACC invasion are supported geometry classes; excluded mismatches are unsupported, not proved unable to react.".into(),"CACC invasion classifies a duplex-anchored3′flap, a completely displaced4ntchain, or displacement with unpaired5′continuation. Only anchored duplex continuation matches the standard directional diagram; extended graph branches imply no kit efficiency.".into(),"Unknown/unsupported steps retain their incoming state. Later transfers describe that explicit retained-state branch, not an experimental distribution. Released enzyme recycling/secondary cleavage is not simulated.".into()]};
        audit_sources(s, &report)?;
        report_size(&report, l.max_report_bytes)?;
        Ok(report)
    }
}
pub fn simulate(s: &Substrates, steps: &[Transfer], l: &Limits) -> Result<Report, Error> {
    let n = preflight(s, steps, l)?;
    let mut state = ExecutionState::new(s, n)?;
    state.apply(s, steps);
    state.finish(s, steps, l, n)
}
/// A read-only checkpoint over an exact whole-source scope. No graph is input.
#[derive(Debug, Clone)]
pub struct Checkpoint {
    pub after_transfers: usize,
    pub strand_start: usize,
    pub strand_end: usize,
    pub activation_unmap: BTreeMap<String, String>,
    pub transfer_map: Vec<usize>,
}
/// Fresh native execution plus documentary source-scope snapshots.
pub struct Executed {
    report: Report,
    checkpoints: Vec<ProductGraph>,
    work_cells: usize,
}
impl Executed {
    pub fn reserve_work(&mut self, cells: usize, peak_bytes: usize) -> Result<(), Error> {
        self.work_cells = add(self.work_cells, cells)?;
        bounded(
            self.work_cells,
            self.report.limits.max_work_cells,
            Budget::WorkCells,
        )?;
        bounded(
            peak_bytes,
            self.report.limits.max_work_bytes,
            Budget::WorkBytes,
        )
    }
    pub fn report(&self) -> &Report {
        &self.report
    }
    pub fn checkpoints(&self) -> &[ProductGraph] {
        &self.checkpoints
    }
    pub fn into_report(self) -> Report {
        self.report
    }
}
fn scoped_product(
    s: &Substrates,
    g: &Graph,
    released: &[String],
    c: &Checkpoint,
) -> Result<ProductGraph, Error> {
    let mut p = output(s, g, released.to_vec())?;
    let inside = |i: usize| c.strand_start <= i && i < c.strand_end;
    let mut indices = BTreeMap::new();
    let mut chains = vec![];
    let unmap_three = |t: &mut ThreePrime| {
        if let ThreePrime::Phosphotyrosyl { activation_id, .. } = t {
            if let Some(old) = c.activation_unmap.get(activation_id) {
                *activation_id = old.clone();
            }
        }
    };
    for (i, mut chain) in p.chains.drain(..).enumerate() {
        let count = chain
            .sources
            .iter()
            .filter(|m| inside(m.source_strand))
            .count();
        if count == 0 {
            continue;
        }
        if count != chain.sources.len() {
            return Err(bad("checkpoint cuts an actual covalent chain"));
        }
        for m in &mut chain.sources {
            m.source_strand -= c.strand_start;
        }
        if let ChainEnds::Linear { three_prime, .. } = &mut chain.ends {
            unmap_three(three_prime)
        }
        indices.insert(i, chains.len());
        chains.push(chain);
    }
    let mut pairing = vec![];
    for mut span in p.source_pairing {
        let a = inside(span.strand_a);
        let b = inside(span.strand_b);
        if a != b {
            return Err(bad("checkpoint discards an actual base-pair partner"));
        }
        if a {
            span.strand_a -= c.strand_start;
            span.strand_b -= c.strand_start;
            pairing.push(span)
        }
    }
    let mut nicks = vec![];
    for mut nick in p.nicks {
        let a = inside(nick.three.strand);
        let b = inside(nick.five.strand);
        if a != b {
            return Err(bad("checkpoint cuts an actual nick pairing scope"));
        }
        if a {
            nick.three.strand -= c.strand_start;
            nick.five.strand -= c.strand_start;
            unmap_three(&mut nick.three_prime);
            nicks.push(nick)
        }
    }
    let mut unpaired = vec![];
    for mut u in p.unpaired {
        if let Some(i) = indices.get(&u.chain) {
            u.chain = *i;
            for step in &mut u.displaced_by {
                *step = c
                    .transfer_map
                    .iter()
                    .position(|x| x == step)
                    .ok_or_else(|| bad("checkpoint displacement has no retained native owner"))?;
            }
            unpaired.push(u)
        }
    }
    let mut components = vec![];
    for mut component in p.components {
        let count = component
            .chains
            .iter()
            .filter(|i| indices.contains_key(i))
            .count();
        if count == 0 {
            continue;
        }
        if count != component.chains.len() {
            return Err(bad(
                "checkpoint discards another connected source component",
            ));
        }
        for i in &mut component.chains {
            *i = indices[&*i]
        }
        components.push(component)
    }
    let released_activations = p
        .released_activations
        .into_iter()
        .filter_map(|id| c.activation_unmap.get(&id).cloned())
        .collect();
    Ok(ProductGraph {
        chains,
        source_pairing: pairing,
        nicks,
        unpaired,
        components,
        released_activations,
    })
}
/// Checkpoints observe only a fresh private execution; supplied reports cannot
/// initialize it. Scope extraction retains all calls/components within its native owner.
fn execute_prefix(
    s: &Substrates,
    steps: &[Transfer],
    l: &Limits,
    checkpoints: &[Checkpoint],
) -> Result<(ExecutionState, Vec<ProductGraph>, usize), Error> {
    let n = preflight(s, steps, l)?;
    if checkpoints.len() > 128 {
        return Err(bad("checkpoint count bound"));
    }
    bounded(
        mul(mul(n, 64)?, checkpoints.len().saturating_add(1))?,
        l.max_work_cells,
        Budget::WorkCells,
    )?;
    bounded(
        mul(n, add(256, mul(checkpoints.len(), 96)?)?)?,
        l.max_work_bytes,
        Budget::WorkBytes,
    )?;
    for c in checkpoints {
        if c.after_transfers > steps.len()
            || c.strand_start >= c.strand_end
            || c.strand_end > s.strands.len()
            || c.transfer_map.len() > l.max_steps
            || c.transfer_map.iter().any(|i| *i >= c.after_transfers)
            || c.activation_unmap.len() > l.max_strands
            || c.activation_unmap
                .iter()
                .any(|(a, b)| !text(a, 256) || !text(b, 256))
        {
            return Err(bad("checkpoint source/transcript scope bound"));
        }
    }
    let mut state = ExecutionState::new(s, n)?;
    let mut out = vec![None; checkpoints.len()];
    for at in 0..=steps.len() {
        for (i, c) in checkpoints.iter().enumerate() {
            if c.after_transfers == at {
                out[i] = Some(scoped_product(s, &state.graph, &state.released, c)?);
            }
        }
        if at < steps.len() {
            state.apply(s, &steps[at..at + 1]);
        }
    }
    let checkpoints = out
        .into_iter()
        .map(|p| p.ok_or_else(|| bad("checkpoint was not observed")))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((state, checkpoints, n))
}
pub fn simulate_with_checkpoints(
    s: &Substrates,
    steps: &[Transfer],
    l: &Limits,
    checkpoints: &[Checkpoint],
) -> Result<Executed, Error> {
    let (state, checkpoints, n) = execute_prefix(s, steps, l, checkpoints)?;
    let report = state.finish(s, steps, l, n)?;
    report_size(&(&report, &checkpoints), l.max_report_bytes)?;
    let work_cells = add(
        mul(mul(n, 64)?, checkpoints.len() + 1)?,
        mul(steps.len(), 128)?,
    )?;
    Ok(Executed {
        report,
        checkpoints,
        work_cells,
    })
}
fn report_size<T: Serialize>(v: &T, max: usize) -> Result<(), Error> {
    struct Count {
        used: usize,
        max: usize,
    }
    impl std::io::Write for Count {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.used = self
                .used
                .checked_add(b.len())
                .filter(|n| *n <= self.max)
                .ok_or_else(|| std::io::Error::other("report budget"))?;
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(&mut Count { used: 0, max }, v).map_err(|_| Error::Limit {
        kind: Budget::ReportBytes,
        required: max.saturating_add(1),
        limit: max,
    })
}
/// Independent source-call reconstruction; this audits ancestry, not whether an
/// imported positive chemical report is current. `replay` reruns eligibility.
pub fn audit_sources(s: &Substrates, r: &Report) -> Result<(), Error> {
    preflight(s, &[], &r.limits)?;
    if r.source_revisions.len() != s.strands.len() {
        return Err(bad("source identity count mismatch"));
    }
    for (source, revision) in s.strands.iter().zip(&r.source_revisions) {
        if hash(source)? != *revision {
            return Err(bad("immutable strand identity mismatch"));
        }
    }

    bounded(r.product.chains.len(), s.strands.len(), Budget::Strands)?;
    let mut coverage = vec![false; s.strands.len()];
    let mut total = 0;
    for chain in &r.product.chains {
        let mut cursor = 0;
        bounded(chain.sources.len(), s.strands.len(), Budget::Strands)?;
        for p in &chain.sources {
            if p.source_strand >= s.strands.len()
                || coverage[p.source_strand]
                || p.start != 0
                || p.end != s.strands[p.source_strand].sequence.len()
                || p.output_start != cursor
            {
                return Err(bad(
                    "source placement is not exact once-only whole-strand ancestry",
                ));
            }
            let next = add(cursor, p.end)?;
            if chain.sequence.as_str().as_bytes().get(cursor..next)
                != Some(s.strands[p.source_strand].sequence.as_str().as_bytes())
            {
                return Err(bad("output bases differ from immutable source calls"));
            }
            coverage[p.source_strand] = true;
            cursor = next;
        }
        if cursor != chain.sequence.len() {
            return Err(bad("incomplete chain mapping"));
        }
        total = add(total, cursor)?;
    }
    if coverage.iter().any(|x| !*x) || total != r.source_bases {
        return Err(bad("source nucleotide conservation failed"));
    }
    Ok(())
}
pub fn replay(s: &Substrates, steps: &[Transfer], l: &Limits, r: &Report) -> Result<(), Error> {
    if r.batch_revision != batch_revision(s, steps, l)?
        || r.engine != ENGINE
        || r.model != PRIMARY_MODEL
        || &r.limits != l
    {
        return Err(bad("report source/model/options binding mismatch"));
    }
    audit_sources(s, r)?;
    if &simulate(s, steps, l)? != r {
        return Err(bad(
            "report differs from fresh physical graph reconstruction",
        ));
    }
    Ok(())
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct OrientationCandidates {
    pub engine: String,
    pub candidates: Vec<Report>,
    pub checkpoints: Vec<ProductGraph>,
    pub pending_transfers: Vec<Vec<Transfer>>,
    pub satisfied_prefix_steps: Vec<Vec<usize>>,
    pub orientation_search_complete: bool,
    pub biochemical_mixture_exhaustive: bool,
}
/// Two explicit orientations of ONE vector and ONE insert. All other mixture
/// paths remain outside this enumeration, even if only one candidate applies.
pub fn orientations(
    s: &Substrates,
    prefix: &[Transfer],
    vector: [usize; 2],
    insert: [usize; 2],
    directional: bool,
    l: &Limits,
    checkpoints: &[Checkpoint],
) -> Result<OrientationCandidates, Error> {
    let mut candidates = vec![];
    let (initial, checkpoints, n) = execute_prefix(s, prefix, l, checkpoints)?;
    let mut pending_transfers = vec![];
    let mut satisfied_prefix_steps = vec![];
    if vector
        .iter()
        .chain(insert.iter())
        .copied()
        .collect::<BTreeSet<_>>()
        .len()
        != 4
    {
        return Err(bad("four distinct vector/insert source strands required"));
    }
    bounded(mul(n, 2)?, l.max_output_bases, Budget::OutputBases)?;
    bounded(
        add(
            mul(mul(n, 64)?, checkpoints.len() + 3)?,
            mul(add(prefix.len(), 4)?, 128)?,
        )?,
        l.max_work_cells,
        Budget::WorkCells,
    )?;
    bounded(
        mul(n, add(512, mul(checkpoints.len(), 96)?)?)?,
        l.max_work_bytes,
        Budget::WorkBytes,
    )?;
    for flip in [false, true] {
        let accept = if flip { [insert[1], insert[0]] } else { insert };
        let steps = [
            Transfer {
                donor: vector[0],
                acceptor: accept[0],
                geometry: if directional {
                    Geometry::DirectionalCacc
                } else {
                    Geometry::ThreePrimeTA
                },
            },
            Transfer {
                donor: vector[1],
                acceptor: accept[1],
                geometry: if directional {
                    Geometry::Flush
                } else {
                    Geometry::ThreePrimeTA
                },
            },
        ];
        let mut pending = vec![];
        let mut satisfied = vec![];
        for t in steps {
            let donor = initial.graph.offsets[t.donor] + s.strands[t.donor].sequence.len() - 1;
            let acceptor = initial.graph.offsets[t.acceptor];
            let prior = initial
                .evidence
                .iter()
                .find(|e| e.transfer == t && e.outcome == Outcome::Applied);
            if initial.graph.next[donor] == Some(acceptor)
                && initial.graph.prev[acceptor] == Some(donor)
            {
                let e = prior
                    .ok_or_else(|| bad("existing bond lacks an exact elementary transfer owner"))?;
                satisfied.push(e.step);
            } else {
                pending.push(t);
            }
        }
        let all = prefix.iter().chain(&pending).cloned().collect::<Vec<_>>();
        preflight(s, &all, l)?;
        let mut state = initial.clone();
        state.apply(s, &pending);
        candidates.push(state.finish(s, &all, l, n)?);
        pending_transfers.push(pending);
        satisfied_prefix_steps.push(satisfied);
    }
    let report = OrientationCandidates {
        engine: ENGINE.into(),
        candidates,
        checkpoints,
        pending_transfers,
        satisfied_prefix_steps,
        orientation_search_complete: true,
        biochemical_mixture_exhaustive: false,
    };
    report_size(&report, l.max_report_bytes)?;
    Ok(report)
}

/// Source-bound paired-contour proof, preserving the full physical graph.
pub mod projection;

/// Explicit source-reconstructed downstream DNA ligase nick sealing.
pub mod nick_sealing;

/// Paired-contour proof after explicit source-reconstructed treatment.
pub mod processed_projection;

pub mod processing;

/// Contour proof after extensible accounted chemical actions.
pub mod processing_projection;
