//! Reference-free affine-overlap/read-layout consensus under explicit evidence thresholds.
//! No external aligner, Phred calibration, laboratory identity or empirical uniqueness is inferred.
use crate::{analysis, oligo_formats, traces};
use collagene_model::{Alphabet, Error, Location, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const CONTIG_ENGINE_VERSION: &str = "contig-1";
const MAX_READS: usize = 32;
const MAX_CALLS: usize = 10_000;
const MAX_OUTPUT: usize = 16 * 1024 * 1024;
const NEG: i32 = -1_000_000_000;
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)] #[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))] #[serde(deny_unknown_fields)] $($t)*};}
dto! {pub struct ReadInput{pub id:String,pub molecule:Molecule}}
dto! {pub struct ContigTrimOptions{pub quality_windows:bool,pub raw_score_threshold:u8,pub window_size:usize,pub minimum_passing:usize,pub min_retained_bases:usize,pub allow_unscored:bool}}
impl Default for ContigTrimOptions {
    fn default() -> Self {
        Self {
            quality_windows: false,
            raw_score_threshold: 20,
            window_size: 5,
            minimum_passing: 4,
            min_retained_bases: 20,
            allow_unscored: true,
        }
    }
}
dto! {pub struct OverlapScoring{pub match_score:i32,pub ambiguous_score:i32,pub mismatch_score:i32,pub gap_open:i32,pub gap_extend:i32}}
impl Default for OverlapScoring {
    fn default() -> Self {
        Self {
            match_score: 2,
            ambiguous_score: 0,
            mismatch_score: -3,
            gap_open: -4,
            gap_extend: -1,
        }
    }
}
dto! {pub struct ContigOptions{pub trim:ContigTrimOptions,pub scoring:OverlapScoring,pub min_overlap:usize,pub max_error_bp:u16,pub both_strands:bool,pub min_support:usize,pub require_both_strands:bool,pub max_pair_cells:usize,pub max_total_cells:usize,pub max_edges:usize,pub max_columns:usize}}
impl Default for ContigOptions {
    fn default() -> Self {
        Self {
            trim: Default::default(),
            scoring: Default::default(),
            min_overlap: 20,
            max_error_bp: 1000,
            both_strands: true,
            min_support: 1,
            require_both_strands: false,
            max_pair_cells: 4_000_000,
            max_total_cells: 128_000_000,
            max_edges: 4096,
            max_columns: 100_000,
        }
    }
}
dto! {pub struct ReadReport{pub id:String,pub name:String,pub source_revision:String,pub source_length:usize,pub source_call_length:usize,pub source_reverse:bool,pub retained_start:usize,pub retained_end:usize,pub trim_start:usize,pub trim_end:usize,pub quality_source:String,pub eligible_calls:usize,pub duplicate_group:Option<String>,pub warnings:Vec<String>,pub error:Option<String>}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum OverlapKind {
    Extension,
    Containment,
    Duplicate,
}
dto! {pub struct Overlap{pub id:usize,pub from:usize,pub to:usize,pub from_reverse:bool,pub to_reverse:bool,pub score:i32,pub from_start:usize,pub from_end:usize,pub to_start:usize,pub to_end:usize,pub reference_aligned:String,pub query_aligned:String,pub paired_bases:usize,pub certain_matches:usize,pub error_columns:usize,pub columns:usize,pub optimal_paths_saturated:u8,pub accepted:bool,pub kind:OverlapKind,pub reason:Option<String>,pub selected:bool}}
dto! {pub struct Observation{pub read_index:usize,pub call_index:usize,pub stored_call_index:usize,pub sample_index:Option<u32>,pub stored_sample_index:Option<u32>,pub base:char,pub raw_score:Option<u8>,pub eligible:bool}}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ConsensusState {
    SupportedBase,
    SupportedGap,
    ConflictingCalls,
    GapConflict,
    InsufficientSupport,
    LowQualityOrAmbiguousOnly,
}
dto! {pub struct ConsensusColumn{pub layout_index:usize,pub consensus_position:Option<usize>,pub base:Option<char>,pub state:ConsensusState,pub counts:[usize;4],pub observations:Vec<Observation>,pub gap_reads:Vec<usize>,pub eligible_gap_reads:Vec<usize>}}
dto! {pub struct CallPlacement{pub call_index:usize,pub stored_call_index:usize,pub column:usize}}
dto! {pub struct Placement{pub read_index:usize,pub reverse:bool,pub start:usize,pub end:usize,pub calls:Vec<CallPlacement>}}
dto! {pub struct Contig{pub id:usize,pub name:String,pub sequence:String,pub molecule:Option<Molecule>,pub columns:Vec<ConsensusColumn>,pub placements:Vec<Placement>,pub complete:bool,pub layout_unique:bool,pub warnings:Vec<String>}}
dto! {pub struct UnassembledRead{pub read_index:usize,pub reason:String}}
dto! {pub struct ContigReport{pub engine_version:String,pub batch_revision:String,pub options_revision:String,pub options:ContigOptions,pub reads:Vec<ReadReport>,pub overlaps:Vec<Overlap>,pub contigs:Vec<Contig>,pub unassembled:Vec<UnassembledRead>,pub complete:bool,pub cells_used:usize,pub limits:Vec<String>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("contig assembly: {s}"))
}
struct Counter(usize);
impl Write for Counter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > MAX_OUTPUT.saturating_sub(self.0) {
            return Err(std::io::Error::other("contig size bound"));
        }
        self.0 += b.len();
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded(v: &impl Serialize) -> Result<(), Error> {
    serde_json::to_writer(Counter(0), v).map_err(|_| Error::TooLarge)
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
}
fn location(l: &Location, depth: usize, n: &mut usize) -> Result<(), Error> {
    *n = n.saturating_add(1);
    if depth > 32 || *n > 100_000 {
        return Err(Error::TooLarge);
    }
    match l {
        Location::Complement(x) => location(x, depth + 1, n)?,
        Location::Join(xs) | Location::Order(xs) => {
            for x in xs {
                location(x, depth + 1, n)?
            }
        }
        _ => (),
    }
    Ok(())
}
impl ContigOptions {
    pub fn validate(&self) -> Result<(), Error> {
        let t = &self.trim;
        let s = &self.scoring;
        if t.window_size == 0
            || t.window_size > 1000
            || t.minimum_passing == 0
            || t.minimum_passing > t.window_size
            || t.min_retained_bases == 0
            || t.min_retained_bases > MAX_CALLS
            || self.min_overlap == 0
            || self.min_overlap > MAX_CALLS
            || self.max_error_bp > 10_000
            || self.min_support == 0
            || self.min_support > MAX_READS
            || self.max_pair_cells == 0
            || self.max_pair_cells > 16_000_000
            || self.max_total_cells == 0
            || self.max_total_cells > 512_000_000
            || self.max_edges == 0
            || self.max_edges > 4096
            || self.max_columns == 0
            || self.max_columns > 100_000
            || s.match_score < 1
            || s.match_score > 20
            || s.ambiguous_score > 0
            || s.ambiguous_score < -20
            || s.mismatch_score >= 0
            || s.mismatch_score < -50
            || s.gap_open >= 0
            || s.gap_open < -100
            || s.gap_extend >= 0
            || s.gap_extend < -100
        {
            return Err(bad(
                "invalid trimming, affine score, overlap, support or resource bounds",
            ));
        }
        Ok(())
    }
}
fn preflight(reads: &[ReadInput], o: &ContigOptions) -> Result<(), Error> {
    o.validate()?;
    if reads.is_empty() || reads.len() > MAX_READS {
        return Err(bad("select1..32 read records"));
    }
    let mut ids = BTreeSet::new();
    let mut nodes = 0usize;
    let mut calls = 0usize;
    for r in reads {
        if r.id.is_empty()
            || r.id.len() > 256
            || r.id.chars().any(char::is_control)
            || !ids.insert(&r.id)
        {
            return Err(bad(
                "read role IDs must be distinct nonempty bounded strings",
            ));
        }
        let m = &r.molecule;
        if m.name.len() > 4096
            || m.sequence.len() > MAX_CALLS
            || m.features.len() > 4096
            || m.metadata.len() > 4096
        {
            return Err(Error::TooLarge);
        }
        nodes = nodes.saturating_add(1 + m.features.len() + m.metadata.len());
        if nodes > 100_000 {
            return Err(Error::TooLarge);
        }
        calls = calls.saturating_add(m.sequence.len());
        if calls > 100_000 {
            return Err(Error::TooLarge);
        }
        for f in &m.features {
            if f.qualifiers.len() > 10000 {
                return Err(Error::TooLarge);
            }
            nodes = nodes.saturating_add(f.qualifiers.len());
            if nodes > 100_000 {
                return Err(Error::TooLarge);
            }
            location(&f.location, 0, &mut nodes)?;
        }
        if m.sequence.alphabet() != Alphabet::Dna || m.topology != Topology::Linear {
            return Err(bad(
                "contig inputs must be explicitly linear DNA called-sequence documents",
            ));
        }
    }
    bounded(&reads)?;
    for r in reads {
        r.molecule.validate()?;
    }
    Ok(())
}
pub fn options_revision(o: &ContigOptions) -> Result<String, Error> {
    o.validate()?;
    Ok(
        blake3::hash(&serde_json::to_vec(o).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
pub fn batch_revision(reads: &[ReadInput], o: &ContigOptions) -> Result<String, Error> {
    preflight(reads, o)?;
    let rows = reads
        .iter()
        .map(|r| (&r.id, r.molecule.revision()))
        .collect::<Vec<_>>();
    Ok(blake3::hash(
        &serde_json::to_vec(&(CONTIG_ENGINE_VERSION, rows, options_revision(o)?))
            .map_err(|_| Error::TooLarge)?,
    )
    .to_hex()
    .to_string())
}
#[derive(Clone)]
struct Prepared {
    sequence: Vec<u8>,
    scores: Vec<Option<u8>>,
    eligible: Vec<bool>,
    peaks: Option<Vec<u32>>,
    samples: usize,
    source_reverse: bool,
    retained_start: usize,
    retained_end: usize,
    source_calls: usize,
    trim_start: usize,
    trim_end: usize,
    evidence: String,
}
impl Prepared {
    fn len(&self) -> usize {
        self.trim_end - self.trim_start
    }
    fn call(&self, at: usize, reverse: bool) -> usize {
        if reverse {
            self.trim_end - 1 - at
        } else {
            self.trim_start + at
        }
    }
    fn base(&self, at: usize, reverse: bool) -> u8 {
        let b = self.sequence[self.call(at, reverse)];
        if reverse {
            collagene_model::complement(b, Alphabet::Dna)
        } else {
            b
        }
    }
    fn observation(&self, read_index: usize, at: usize, reverse: bool) -> Observation {
        let call_index = self.call(at, reverse);
        let sample = self.peaks.as_ref().map(|p| p[call_index]);
        Observation {
            read_index,
            call_index,
            stored_call_index: if self.source_reverse {
                self.retained_end - 1 - call_index
            } else {
                self.retained_start + call_index
            },
            sample_index: sample,
            stored_sample_index: sample.map(|s| {
                if self.source_reverse {
                    (self.samples - 1) as u32 - s
                } else {
                    s
                }
            }),
            base: self.base(at, reverse) as char,
            raw_score: self.scores[call_index],
            eligible: self.eligible[call_index],
        }
    }
    fn text(&self, reverse: bool) -> Vec<u8> {
        (0..self.len()).map(|i| self.base(i, reverse)).collect()
    }
}
fn prepare(input: &ReadInput, o: &ContigOptions) -> Result<(Prepared, ReadReport), Error> {
    let m = &input.molecule;
    let n = m.sequence.len();
    let mut scores = vec![None; n];
    let mut peaks = None;
    let mut samples = 0;
    let mut source_reverse = false;
    let mut retained_start = 0;
    let mut retained_end = n;
    let mut source_calls = n;
    let mut warnings = vec![];
    let quality_source;
    let evidence;
    let is_trace = m.metadata.contains_key("trace_original_base64");
    let is_fastq = m.metadata.contains_key(oligo_formats::FASTQ_ORIGINAL_KEY);
    if is_trace && is_fastq {
        return Err(bad("choose a single quality-evidence source per read"));
    }
    if is_trace {
        let t = traces::trace_from_molecule(m)?;
        source_reverse = t.orientation == traces::TraceOrientation::ReverseComplement;
        source_calls = t.sequence.len();
        retained_end = source_calls;
        peaks = Some(t.peak_positions.clone());
        samples = t.channels.len();
        evidence = format!("trace:{}", t.source_hash);
        match &t.qualities {
            Some(traces::QualityEvidence::AbifBasecallerScores(q)) => {
                quality_source = "ABIF native raw basecaller scores".to_string();
                scores = q.iter().copied().map(Some).collect()
            }
            Some(traces::QualityEvidence::ScfAccuracy(q)) => {
                quality_source = "SCF called-base native accuracy".to_string();
                scores = t
                    .sequence
                    .as_str()
                    .bytes()
                    .enumerate()
                    .map(|(i, b)| match b {
                        b'A' => Some(q.a[i]),
                        b'C' => Some(q.c[i]),
                        b'G' => Some(q.g[i]),
                        b'T' => Some(q.t[i]),
                        _ => None,
                    })
                    .collect()
            }
            None => quality_source = "Missing".into(),
        };
        warnings.extend(t.warnings);
        warnings.push(
            "Native trace scores are threshold evidence, not calibrated Phred error probabilities."
                .into(),
        );
    } else if is_fastq {
        let e = oligo_formats::quality_evidence(m)?;
        if !e.current_valid {
            return Err(bad(
                "current calls do not match immutable FASTQ quality mapping",
            ));
        }
        source_reverse = e.reverse;
        retained_start = e.retained_start;
        retained_end = e.retained_end;
        source_calls = e.source_call_length;
        scores = e.qualities.iter().copied().map(Some).collect();
        quality_source = match e.encoding {
            oligo_formats::QualityEncoding::Phred33 => "FASTQ explicit Phred+33",
            oligo_formats::QualityEncoding::Phred64 => "FASTQ explicit Phred+64",
        }
        .into();
        evidence = format!("fastq:{}", e.original_record_blake3);
        warnings.extend(e.warnings);
    } else {
        quality_source = "Missing".into();
        let rc = m.sequence.reverse_complement()?;
        let canonical = m.sequence.as_str().min(rc.as_str());
        evidence = format!("unscored:{}", blake3::hash(canonical.as_bytes()).to_hex());
    }
    let missing = scores.iter().all(Option::is_none);
    if missing && !o.trim.allow_unscored {
        return Err(bad(
            "read lacks native quality; explicitly allow unscored records before assembly",
        ));
    }
    let eligible = scores
        .iter()
        .zip(m.sequence.as_str().bytes())
        .map(|(q, b)| {
            analysis::base_mask(b).count_ones() == 1
                && if missing {
                    o.trim.allow_unscored
                } else {
                    q.is_some_and(|q| q >= o.trim.raw_score_threshold)
                }
        })
        .collect::<Vec<_>>();
    let (mut start, mut end) = (0, n);
    if o.trim.quality_windows && !missing {
        let w = o.trim.window_size.min(n);
        let required = o.trim.minimum_passing.min(w);
        let mut count = eligible[..w].iter().filter(|v| **v).count();
        let mut first = None;
        let mut last = 0;
        for at in 0..=n - w {
            if at > 0 {
                count -= usize::from(eligible[at - 1]);
                count += usize::from(eligible[at + w - 1]);
            }
            if count >= required {
                first.get_or_insert(at);
                last = at + w;
            }
        }
        (start, end) = first.map_or((0, 0), |s| (s, last));
    }
    if end - start < o.trim.min_retained_bases {
        return Err(bad(
            "no interval meets the retained-call and trimming policy",
        ));
    }
    if missing {
        warnings.push("Unscored canonical calls support consensus only under the explicitly declared allow_unscored policy; no confidence probability is invented.".into());
    }
    if warnings.len() > 63 {
        warnings.truncate(63);
        warnings.push("Additional native source warnings omitted from preview; retain the original read evidence.".into());
    }
    let p = Prepared {
        sequence: m.sequence.as_str().as_bytes().to_vec(),
        scores,
        eligible,
        peaks,
        samples,
        source_reverse,
        retained_start,
        retained_end,
        source_calls,
        trim_start: start,
        trim_end: end,
        evidence,
    };
    let r = ReadReport {
        id: input.id.clone(),
        name: m.name.clone(),
        source_revision: m.revision(),
        source_length: n,
        source_call_length: source_calls,
        source_reverse,
        retained_start,
        retained_end,
        trim_start: start,
        trim_end: end,
        quality_source,
        eligible_calls: p.eligible[start..end].iter().filter(|v| **v).count(),
        duplicate_group: None,
        warnings,
        error: None,
    };
    Ok((p, r))
}
#[derive(Clone, Copy)]
struct Value {
    score: i32,
    count: u8,
}
const EMPTY: Value = Value {
    score: NEG,
    count: 0,
};
fn best(values: [Value; 3]) -> (Value, u8) {
    let mut out = EMPTY;
    let mut state = 0;
    for (i, v) in values.into_iter().enumerate() {
        if v.count == 0 {
            continue;
        }
        if v.score > out.score {
            out = v;
            state = 1 << i;
        } else if v.score == out.score {
            out.count = out.count.saturating_add(v.count).min(2);
            state |= 1 << i;
        }
    }
    (out, state)
}
fn add(v: Value, x: i32) -> Value {
    if v.count == 0 {
        EMPTY
    } else {
        Value {
            score: v.score + x,
            count: v.count,
        }
    }
}
fn set(trace: &mut [u16], cell: usize, state: usize, previous: u8) {
    let shift = state * 3;
    trace[cell] = (trace[cell] & !(7 << shift)) | ((previous as u16) << shift);
}
fn sub(a: u8, b: u8, s: &OverlapScoring) -> i32 {
    let x = analysis::base_mask(a);
    let y = analysis::base_mask(b);
    if x & y == 0 {
        s.mismatch_score
    } else if x.count_ones() == 1 && x == y {
        s.match_score
    } else {
        s.ambiguous_score
    }
}
#[derive(Clone, Copy)]
struct Endpoint {
    i: usize,
    j: usize,
    state: usize,
    value: Value,
}
fn eligibility_score_floor(o: &ContigOptions) -> i32 {
    let e = i64::from(o.max_error_bp);
    let denominator = 10_000 - e;
    if denominator == 0 {
        return 1;
    }
    let p = [
        o.scoring.ambiguous_score,
        o.scoring.mismatch_score,
        o.scoring.gap_open,
        o.scoring.gap_extend,
    ]
    .into_iter()
    .min()
    .unwrap() as i64;
    let coefficient = i64::from(o.scoring.match_score) * denominator + p * e;
    if coefficient <= 0 {
        1
    } else {
        (((o.min_overlap as i64) * coefficient + denominator - 1) / denominator).max(1) as i32
    }
}
fn trace_endpoint(
    trace: &[u16],
    a: &[u8],
    b: &[u8],
    o: &ContigOptions,
    from: usize,
    to: usize,
    ar: bool,
    br: bool,
    end: Endpoint,
    remaining: &mut usize,
    pa: &Prepared,
    pb: &Prepared,
) -> Result<(Overlap, bool), ()> {
    let (mut i, mut j, mut state) = (end.i, end.j, end.state);
    let mut aa = vec![];
    let mut bb = vec![];
    while j > 0 && i > 0 {
        if *remaining == 0 {
            return Err(());
        }
        *remaining -= 1;
        let mask = (trace[i * (b.len() + 1) + j] >> (state * 3)) & 7;
        let pred = mask.trailing_zeros();
        match state {
            0 => {
                aa.push(a[i - 1]);
                bb.push(b[j - 1]);
                i -= 1;
                j -= 1
            }
            1 => {
                aa.push(b'-');
                bb.push(b[j - 1]);
                j -= 1
            }
            2 => {
                aa.push(a[i - 1]);
                bb.push(b'-');
                i -= 1
            }
            _ => break,
        }
        if pred > 2 {
            break;
        }
        state = pred as usize;
    }
    aa.reverse();
    bb.reverse();
    let (mut paired, mut certain, mut errors, mut qualified) = (0, 0, 0, 0);
    let (mut ai, mut bi) = (i, j);
    for (x, y) in aa.iter().zip(&bb) {
        if *x != b'-' && *y != b'-' {
            paired += 1;
            let u = analysis::base_mask(*x);
            let v = analysis::base_mask(*y);
            if u == v && u.count_ones() == 1 {
                certain += 1;
                if pa.eligible[pa.call(ai, ar)] && pb.eligible[pb.call(bi, br)] {
                    qualified += 1
                }
            } else {
                errors += 1
            }
        } else {
            errors += 1
        }
        ai += usize::from(*x != b'-');
        bi += usize::from(*y != b'-');
    }
    let columns = aa.len();
    let kind = if i == 0 && end.i == a.len() && j == 0 && end.j == b.len() && a == b {
        OverlapKind::Duplicate
    } else if end.j == b.len() {
        OverlapKind::Containment
    } else {
        OverlapKind::Extension
    };
    let eligible = end.value.score > 0
        && qualified >= o.min_overlap
        && columns > 0
        && (errors as u64) * 10_000 <= (columns as u64) * u64::from(o.max_error_bp)
        && j == 0;
    Ok((
        Overlap {
            id: 0,
            from,
            to,
            from_reverse: ar,
            to_reverse: br,
            score: end.value.score,
            from_start: i,
            from_end: end.i,
            to_start: j,
            to_end: end.j,
            reference_aligned: String::from_utf8(aa).unwrap(),
            query_aligned: String::from_utf8(bb).unwrap(),
            paired_bases: paired,
            certain_matches: certain,
            error_columns: errors,
            columns,
            optimal_paths_saturated: end.value.count,
            accepted: false,
            kind,
            reason: None,
            selected: false,
        },
        eligible,
    ))
}
// Only used when a representative failed: exhaust the optimal-path DAG under
// the shared work budget to prove that no alternate path qualifies. A possible
// qualifying alternate remains unresolved, rather than inventing its placement.
fn alternative_may_qualify(
    trace: &[u16],
    a: &[u8],
    b: &[u8],
    o: &ContigOptions,
    end: Endpoint,
    ar: bool,
    br: bool,
    pa: &Prepared,
    pb: &Prepared,
    remaining: &mut usize,
) -> bool {
    let mut stack = vec![(end.i, end.j, end.state, 0usize, 0usize, 0usize)];
    while let Some((i, j, state, mut qualified, mut errors, mut columns)) = stack.pop() {
        if *remaining == 0 {
            return true;
        }
        *remaining -= 1;
        if i == 0 || j == 0 {
            if j == 0
                && qualified >= o.min_overlap
                && columns > 0
                && (errors as u64) * 10_000 <= (columns as u64) * u64::from(o.max_error_bp)
            {
                return true;
            }
            continue;
        }
        columns += 1;
        let (ni, nj) = match state {
            0 => {
                let x = a[i - 1];
                let y = b[j - 1];
                if x == y && analysis::base_mask(x).count_ones() == 1 {
                    if pa.eligible[pa.call(i - 1, ar)] && pb.eligible[pb.call(j - 1, br)] {
                        qualified += 1
                    }
                } else {
                    errors += 1
                }
                (i - 1, j - 1)
            }
            1 => {
                errors += 1;
                (i, j - 1)
            }
            _ => {
                errors += 1;
                (i - 1, j)
            }
        };
        let mask = (trace[i * (b.len() + 1) + j] >> (state * 3)) & 7;
        for pred in 0..3 {
            if mask & (1 << pred) != 0 {
                if stack.len() >= 100_000 {
                    return true;
                }
                stack.push((ni, nj, pred, qualified, errors, columns));
            }
        }
    }
    false
}
/// Each boundary class is independent: an ineligible containment cannot mask a
/// legitimate extension. Uninspected co-optimal eligibility remains unresolved.
fn overlaps(
    a: &[u8],
    b: &[u8],
    o: &ContigOptions,
    from: usize,
    to: usize,
    ar: bool,
    br: bool,
    pa: &Prepared,
    pb: &Prepared,
    budget: usize,
) -> (Vec<Overlap>, usize, bool) {
    let (n, m) = (a.len(), b.len());
    let mut trace = vec![0u16; (n + 1) * (m + 1)];
    let mut previous = vec![[EMPTY; 3]; m + 1];
    previous[0][0] = Value { score: 0, count: 1 };
    let mut ends = [vec![], vec![]];
    for i in 1..=n {
        let mut row = vec![[EMPTY; 3]; m + 1];
        row[0][0] = Value { score: 0, count: 1 };
        for j in 1..=m {
            let (p, ps) =
                best(previous[j - 1].map(|v| add(v, sub(a[i - 1], b[j - 1], &o.scoring))));
            row[j][0] = p;
            set(&mut trace, i * (m + 1) + j, 0, ps);
            if j > 1 {
                let (q, qs) = best([
                    add(row[j - 1][0], o.scoring.gap_open),
                    add(row[j - 1][1], o.scoring.gap_extend),
                    add(row[j - 1][2], o.scoring.gap_open),
                ]);
                row[j][1] = q;
                set(&mut trace, i * (m + 1) + j, 1, qs);
            }
            let (d, ds) = best([
                add(previous[j][0], o.scoring.gap_open),
                add(previous[j][1], o.scoring.gap_open),
                add(previous[j][2], o.scoring.gap_extend),
            ]);
            row[j][2] = d;
            set(&mut trace, i * (m + 1) + j, 2, ds);
            if i == n && j < m {
                for (state, value) in row[j].iter().copied().enumerate() {
                    if value.count > 0 {
                        ends[0].push(Endpoint { i, j, state, value });
                    }
                }
            }
            if j == m && row[j][0].count > 0 {
                ends[1].push(Endpoint {
                    i,
                    j,
                    state: 0,
                    value: row[j][0],
                });
            }
        }
        previous = row;
    }
    let mut remaining = budget;
    let mut complete = true;
    let mut result = vec![];
    let floor = eligibility_score_floor(o);
    let mut qa = vec![0usize; n + 1];
    let mut qb = vec![0usize; m + 1];
    for i in 0..n {
        qa[i + 1] = qa[i] + usize::from(pa.eligible[pa.call(i, ar)]);
    }
    for j in 0..m {
        qb[j + 1] = qb[j] + usize::from(pb.eligible[pb.call(j, br)]);
    }
    for mut endpoints in ends {
        endpoints.sort_by_key(|e| (std::cmp::Reverse(e.value.score), e.i, e.j, e.state));
        let raw_best = endpoints.first().map(|e| e.value.score);
        let raw_paths = endpoints
            .iter()
            .take_while(|e| Some(e.value.score) == raw_best)
            .fold(0u8, |n, e| n.saturating_add(e.value.count).min(2));
        let mut diagnostic = None;
        let mut winner: Option<Overlap> = None;
        let mut unresolved: Option<i32> = None;
        for end in endpoints {
            if winner.as_ref().is_some_and(|w| end.value.score < w.score) {
                break;
            }
            let possible = end.value.score >= floor && qa[end.i].min(qb[end.j]) >= o.min_overlap;
            if !possible && diagnostic.is_some() {
                if end.value.score < floor {
                    break;
                } else {
                    continue;
                }
            }
            let (mut edge, eligible) = match trace_endpoint(
                &trace,
                a,
                b,
                o,
                from,
                to,
                ar,
                br,
                end,
                &mut remaining,
                pa,
                pb,
            ) {
                Ok(x) => x,
                Err(()) => {
                    complete = false;
                    break;
                }
            };
            if diagnostic.is_none() {
                let mut d = edge.clone();
                d.optimal_paths_saturated = raw_paths;
                diagnostic = Some(d)
            }
            if possible
                && end.value.count > 1
                && !eligible
                && alternative_may_qualify(&trace, a, b, o, end, ar, br, pa, pb, &mut remaining)
            {
                unresolved = Some(unresolved.map_or(end.value.score, |s| s.max(end.value.score)));
            }
            if eligible {
                if let Some(w) = &mut winner {
                    w.optimal_paths_saturated = w
                        .optimal_paths_saturated
                        .saturating_add(edge.optimal_paths_saturated)
                        .min(2);
                } else {
                    edge.accepted = edge.optimal_paths_saturated == 1;
                    winner = Some(edge);
                }
            }
        }
        let mut chosen = winner.or(diagnostic);
        if let Some(e) = &mut chosen {
            if unresolved.is_some_and(|s| s >= e.score) {
                complete = false;
                e.accepted = false;
                e.reason=Some("Uninspected co-optimal paths could change threshold-qualified boundary eligibility; this class remains incomplete.".into());
            } else if e.optimal_paths_saturated > 1 && potential_stats(e, o) {
                e.accepted = false;
                e.reason=Some("Multiple co-optimal qualified paths/endpoints; no arbitrary unique boundary is asserted.".into());
            } else if !e.accepted {
                e.reason=Some("No unique optimal-per-boundary candidate met eligible-anchor and error thresholds.".into());
            }
            if !complete {
                e.accepted = false;
                e.reason=Some("Boundary eligibility/path search was limited or unresolved; no unique join is established.".into());
            }
            result.push(e.clone());
        }
    }
    (result, budget - remaining, complete)
}
fn potential_stats(e: &Overlap, o: &ContigOptions) -> bool {
    e.score > 0
        && e.certain_matches >= o.min_overlap
        && e.columns > 0
        && (e.error_columns as u64) * 10_000 <= (e.columns as u64) * u64::from(o.max_error_bp)
}
fn eligible_overlap(e: &Overlap, prep: &[Option<Prepared>], o: &ContigOptions) -> usize {
    let a = prep[e.from].as_ref().unwrap();
    let b = prep[e.to].as_ref().unwrap();
    let (mut ai, mut bi) = (e.from_start, e.to_start);
    let mut n = 0;
    for (x, y) in e.reference_aligned.bytes().zip(e.query_aligned.bytes()) {
        if x != b'-'
            && y != b'-'
            && x == y
            && analysis::base_mask(x).count_ones() == 1
            && a.eligible[a.call(ai, e.from_reverse)]
            && b.eligible[b.call(bi, e.to_reverse)]
        {
            n += 1;
        }
        ai += usize::from(x != b'-');
        bi += usize::from(y != b'-');
    }
    if n < o.min_overlap { 0 } else { n }
}
#[derive(Clone, Default)]
struct ColumnData {
    observations: Vec<Observation>,
    gaps: BTreeSet<usize>,
    eligible_gaps: BTreeSet<usize>,
}
#[derive(Clone)]
struct Mapping {
    reverse: bool,
    columns: Vec<usize>,
}
#[derive(Clone)]
struct Layout {
    columns: Vec<ColumnData>,
    reads: BTreeMap<usize, Mapping>,
    warnings: Vec<String>,
}
impl Layout {
    fn root(read: usize, p: &Prepared, reverse: bool) -> Self {
        Self {
            columns: (0..p.len())
                .map(|at| ColumnData {
                    observations: vec![p.observation(read, at, reverse)],
                    ..Default::default()
                })
                .collect(),
            reads: BTreeMap::from([(
                read,
                Mapping {
                    reverse,
                    columns: (0..p.len()).collect(),
                },
            )]),
            warnings: vec![],
        }
    }
    fn flank(map: &Mapping, p: &Prepared, at: usize) -> bool {
        let boundary = map.columns.partition_point(|c| *c < at);
        boundary > 0
            && boundary < map.columns.len()
            && p.eligible[p.call(boundary - 1, map.reverse)]
            && p.eligible[p.call(boundary, map.reverse)]
    }
    fn insert(&mut self, at: usize, count: usize, prep: &[Option<Prepared>]) {
        let mut new = vec![ColumnData::default(); count];
        for (read, map) in &self.reads {
            if map.columns.first().is_some_and(|c| *c < at)
                && map.columns.last().is_some_and(|c| *c >= at)
            {
                let eligible = Self::flank(map, prep[*read].as_ref().unwrap(), at);
                for c in &mut new {
                    c.gaps.insert(*read);
                    if eligible {
                        c.eligible_gaps.insert(*read);
                    }
                }
            }
        }
        for map in self.reads.values_mut() {
            for c in &mut map.columns {
                if *c >= at {
                    *c += count
                }
            }
        }
        self.columns.splice(at..at, new);
    }
    fn merge(
        &mut self,
        e: &Overlap,
        prep: &[Option<Prepared>],
        o: &ContigOptions,
    ) -> Result<(), Error> {
        if self.reads.contains_key(&e.to) {
            return Err(bad("read already placed"));
        }
        let a = prep[e.from].as_ref().unwrap();
        let b = prep[e.to].as_ref().unwrap();
        if self
            .reads
            .get(&e.from)
            .is_none_or(|m| m.reverse != e.from_reverse)
        {
            return Err(bad("overlap anchor orientation disagrees with layout"));
        }
        let mut map = vec![usize::MAX; b.len()];
        let (mut ai, mut bi) = (e.from_start, e.to_start);
        let aa = e.reference_aligned.as_bytes();
        let bb = e.query_aligned.as_bytes();
        let mut k = 0;
        while k < aa.len() {
            if aa[k] == b'-' {
                let start = k;
                while k < aa.len() && aa[k] == b'-' {
                    k += 1;
                }
                let count = k - start;
                let am = &self.reads[&e.from].columns;
                let left = if ai > 0 { am[ai - 1] + 1 } else { am[0] };
                let right = if ai < a.len() {
                    am[ai]
                } else {
                    self.columns.len()
                };
                if right > left {
                    if right - left != count {
                        return Err(bad(
                            "insertion topology conflicts with an existing read-column interval",
                        ));
                    }
                    for j in 0..count {
                        let c = &self.columns[left + j];
                        if c.observations.iter().any(|x| {
                            analysis::base_mask(x.base as u8) & analysis::base_mask(bb[start + j])
                                == 0
                        }) {
                            return Err(bad(
                                "insertion calls conflict with existing inserted columns",
                            ));
                        }
                    }
                } else {
                    if self.columns.len() + count > o.max_columns {
                        return Err(Error::TooLarge);
                    }
                    self.insert(left, count, prep);
                }
                for j in 0..count {
                    map[bi] = left + j;
                    self.columns[left + j]
                        .observations
                        .push(b.observation(e.to, bi, e.to_reverse));
                    bi += 1;
                }
                continue;
            }
            let col = self.reads[&e.from].columns[ai];
            if bb[k] != b'-' {
                map[bi] = col;
                self.columns[col]
                    .observations
                    .push(b.observation(e.to, bi, e.to_reverse));
                bi += 1;
            }
            ai += 1;
            k += 1;
        }
        if e.to_start != 0 {
            return Err(bad("query prefix placement is unresolved"));
        }
        for at in bi..b.len() {
            if self.columns.len() >= o.max_columns {
                return Err(Error::TooLarge);
            }
            let col = self.columns.len();
            map[at] = col;
            self.columns.push(ColumnData {
                observations: vec![b.observation(e.to, at, e.to_reverse)],
                ..Default::default()
            });
        }
        if map.iter().any(|c| *c == usize::MAX) || map.windows(2).any(|w| w[0] >= w[1]) {
            return Err(bad("read call mapping is incomplete or nonmonotone"));
        }
        let mapping = Mapping {
            reverse: e.to_reverse,
            columns: map,
        };
        let first = mapping.columns[0];
        let last = *mapping.columns.last().unwrap();
        let covered = mapping.columns.iter().copied().collect::<BTreeSet<_>>();
        for col in first..=last {
            if !covered.contains(&col) {
                self.columns[col].gaps.insert(e.to);
                if Self::flank(&mapping, b, col) {
                    self.columns[col].eligible_gaps.insert(e.to);
                }
            }
        }
        self.reads.insert(e.to, mapping);
        Ok(())
    }
}
fn groups(reads: &[ReadReport], read: usize) -> String {
    reads[read]
        .duplicate_group
        .clone()
        .unwrap_or_else(|| format!("read:{read}"))
}
fn consensus_state(
    c: &ColumnData,
    reads: &[ReadReport],
    placements: &BTreeMap<usize, Mapping>,
    o: &ContigOptions,
) -> ([usize; 4], Option<char>, ConsensusState) {
    let mut bases: [BTreeSet<String>; 4] = Default::default();
    let mut direction: BTreeMap<String, BTreeSet<bool>> = BTreeMap::new();
    for obs in &c.observations {
        if obs.eligible {
            let i = match obs.base {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                _ => continue,
            };
            let g = groups(reads, obs.read_index);
            bases[i].insert(g.clone());
            direction
                .entry(g)
                .or_default()
                .insert(placements[&obs.read_index].reverse ^ reads[obs.read_index].source_reverse);
        }
    }
    let gaps = c
        .eligible_gaps
        .iter()
        .map(|r| groups(reads, *r))
        .collect::<BTreeSet<_>>();
    for r in &c.eligible_gaps {
        direction
            .entry(groups(reads, *r))
            .or_default()
            .insert(placements[r].reverse ^ reads[*r].source_reverse);
    }
    let counts = bases.each_ref().map(BTreeSet::len);
    let kinds = counts.iter().filter(|n| **n > 0).count();
    if kinds > 1 {
        return (counts, Some('N'), ConsensusState::ConflictingCalls);
    }
    if kinds > 0 && !gaps.is_empty() {
        return (counts, Some('N'), ConsensusState::GapConflict);
    }
    let good = if kinds == 1 {
        bases.iter().find(|s| !s.is_empty()).unwrap()
    } else {
        &gaps
    };
    let strands = good
        .iter()
        .filter_map(|g| direction.get(g))
        .filter(|s| s.len() == 1)
        .flat_map(|s| s.iter().copied())
        .collect::<BTreeSet<_>>();
    let support = good.len() >= o.min_support && (!o.require_both_strands || strands.len() == 2);
    if support {
        if kinds == 1 {
            let at = counts.iter().position(|n| *n > 0).unwrap();
            return (
                counts,
                Some(b"ACGT"[at] as char),
                ConsensusState::SupportedBase,
            );
        } else {
            return (counts, None, ConsensusState::SupportedGap);
        }
    }
    if kinds > 0 || !gaps.is_empty() {
        (counts, Some('N'), ConsensusState::InsufficientSupport)
    } else {
        (counts, Some('N'), ConsensusState::LowQualityOrAmbiguousOnly)
    }
}
fn materialize(
    id: usize,
    layout: &Layout,
    reads: &[ReadReport],
    prep: &[Option<Prepared>],
    o: &ContigOptions,
    unique: bool,
    graph_complete: bool,
) -> Contig {
    let mut sequence = String::new();
    let columns = layout
        .columns
        .iter()
        .enumerate()
        .map(|(index, c)| {
            let (counts, base, state) = consensus_state(c, reads, &layout.reads, o);
            let position = base.map(|b| {
                let p = sequence.len();
                sequence.push(b);
                p
            });
            ConsensusColumn {
                layout_index: index,
                consensus_position: position,
                base,
                state,
                counts,
                observations: c.observations.clone(),
                gap_reads: c.gaps.iter().copied().collect(),
                eligible_gap_reads: c.eligible_gaps.iter().copied().collect(),
            }
        })
        .collect::<Vec<_>>();
    let placements = layout
        .reads
        .iter()
        .map(|(read, m)| {
            let p = prep[*read].as_ref().unwrap();
            Placement {
                read_index: *read,
                reverse: m.reverse,
                start: m.columns[0],
                end: m.columns.last().unwrap() + 1,
                calls: m
                    .columns
                    .iter()
                    .enumerate()
                    .map(|(at, col)| {
                        let obs = p.observation(*read, at, m.reverse);
                        CallPlacement {
                            call_index: obs.call_index,
                            stored_call_index: obs.stored_call_index,
                            column: *col,
                        }
                    })
                    .collect(),
            }
        })
        .collect();
    let complete = graph_complete
        && unique
        && columns.iter().all(|c| {
            matches!(
                c.state,
                ConsensusState::SupportedBase | ConsensusState::SupportedGap
            )
        });
    let mut warnings = layout.warnings.clone();
    if layout.reads.len() == 1 {
        warnings
            .push("Singleton read evidence; no overlap-supported extension is established.".into())
    }
    if !unique {
        warnings.push("Branch, repeat, cycle, unresolved placement or resource limit prevents a unique layout under the declared reciprocal-best model.".into())
    }
    if !complete {
        warnings.push("N denotes conflict or insufficient eligible support; gap-supported columns are retained as evidence but omitted from consensus DNA.".into())
    }
    Contig {
        id,
        name: format!("Contig {}", id + 1),
        sequence,
        molecule: None,
        columns,
        placements,
        complete,
        layout_unique: unique,
        warnings,
    }
}
fn potential(e: &Overlap, prep: &[Option<Prepared>], o: &ContigOptions) -> bool {
    e.score > 0
        && e.certain_matches >= o.min_overlap
        && e.columns > 0
        && (e.error_columns as u64) * 10_000 <= (e.columns as u64) * u64::from(o.max_error_bp)
        && eligible_overlap(e, prep, o) > 0
}
fn node(read: usize, reverse: bool) -> usize {
    read * 2 + usize::from(reverse)
}
fn add_unassembled(r: &mut ContigReport, read: usize, reason: &str) {
    if !r.unassembled.iter().any(|x| x.read_index == read) {
        r.unassembled.push(UnassembledRead {
            read_index: read,
            reason: reason.into(),
        });
    }
}
fn choose(values: &[usize], edges: &[Overlap]) -> Option<usize> {
    let best = values.iter().copied().max_by_key(|i| edges[*i].score)?;
    if edges[best].optimal_paths_saturated != 1
        || values
            .iter()
            .filter(|i| edges[**i].score == edges[best].score)
            .count()
            != 1
    {
        None
    } else {
        Some(best)
    }
}
pub fn assemble(inputs: &[ReadInput], o: &ContigOptions) -> Result<ContigReport, Error> {
    let batch = batch_revision(inputs, o)?;
    let mut r=ContigReport{engine_version:CONTIG_ENGINE_VERSION.into(),batch_revision:batch,options_revision:options_revision(o)?,options:o.clone(),reads:vec![],overlaps:vec![],contigs:vec![],unassembled:vec![],complete:true,cells_used:0,limits:vec![],warnings:vec!["Reference-free reciprocal-best affine boundary-overlap model; no external target or reference bases are used.".into(),"Graph/model completion and unique selected layouts do not prove a unique physical molecule or biological origin. Suboptimal alignment paths and all possible assemblies are not exhaustively enumerated.".into(),"Read-record support is not a calibrated error probability or a count of independent physical molecules. Duplicate retained evidence contributes at most one support group per column.".into()]};
    let mut prep = vec![];
    for (i, input) in inputs.iter().enumerate() {
        match prepare(input, o) {
            Ok((p, read)) => {
                prep.push(Some(p));
                r.reads.push(read)
            }
            Err(e) => {
                r.complete = false;
                r.reads.push(ReadReport {
                    id: input.id.clone(),
                    name: input.molecule.name.clone(),
                    source_revision: input.molecule.revision(),
                    source_length: input.molecule.sequence.len(),
                    source_call_length: input.molecule.sequence.len(),
                    source_reverse: false,
                    retained_start: 0,
                    retained_end: input.molecule.sequence.len(),
                    trim_start: 0,
                    trim_end: 0,
                    quality_source: "Unavailable".into(),
                    eligible_calls: 0,
                    duplicate_group: None,
                    warnings: vec![],
                    error: Some(e.to_string()),
                });
                prep.push(None);
                add_unassembled(
                    &mut r,
                    i,
                    "Read quality/call evidence could not be prepared under the declared policy.",
                );
            }
        }
    }
    let mut evidence: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, p) in prep.iter().enumerate() {
        if let Some(p) = p {
            evidence.entry(p.evidence.clone()).or_default().push(i);
        }
    }
    for (key, indices) in evidence {
        if indices.len() > 1 {
            for i in indices {
                r.reads[i].duplicate_group = Some(key.clone());
                r.reads[i].warnings.push("Duplicate retained evidence/canonical unscored sequence group: support is not multiplied by repeating the same evidence.".into());
            }
        }
    }
    let mut overlap_bytes = 0usize;
    'pairs: for a in 0..inputs.len() {
        for b in 0..inputs.len() {
            if a == b {
                continue;
            }
            let (Some(pa), Some(pb)) = (&prep[a], &prep[b]) else {
                continue;
            };
            for ar in [false, true]
                .into_iter()
                .take(if o.both_strands { 2 } else { 1 })
            {
                for br in [false, true]
                    .into_iter()
                    .take(if o.both_strands { 2 } else { 1 })
                {
                    let count = (pa.len() + 1)
                        .checked_mul(pb.len() + 1)
                        .ok_or(Error::TooLarge)?;
                    if count > o.max_pair_cells
                        || count > o.max_total_cells.saturating_sub(r.cells_used)
                    {
                        r.complete = false;
                        r.limits.push("Pair/cumulative affine matrix-cell budget reached; absent edges and unique layouts are not established.".into());
                        break 'pairs;
                    }
                    if r.overlaps.len() >= o.max_edges {
                        r.complete = false;
                        r.limits.push("Overlap evidence count budget reached; additional alternatives may be missing.".into());
                        break 'pairs;
                    }
                    r.cells_used += count;
                    let (edges, work, complete) = overlaps(
                        &pa.text(ar),
                        &pb.text(br),
                        o,
                        a,
                        b,
                        ar,
                        br,
                        pa,
                        pb,
                        o.max_total_cells - r.cells_used,
                    );
                    r.cells_used += work;
                    for mut edge in edges {
                        let bytes = edge.reference_aligned.len() + edge.query_aligned.len() + 1024;
                        if r.overlaps.len() >= o.max_edges
                            || overlap_bytes.saturating_add(bytes) > MAX_OUTPUT / 2
                        {
                            r.complete = false;
                            r.limits.push("Overlap evidence text/count budget reached; additional alternatives may be missing.".into());
                            break 'pairs;
                        }
                        overlap_bytes += bytes;
                        edge.id = r.overlaps.len();
                        r.overlaps.push(edge);
                    }
                    if !complete {
                        r.complete = false;
                        if !r.limits.iter().any(|s| s.starts_with("Boundary traceback")) {
                            r.limits.push("Boundary traceback budget or co-optimal eligibility remained unresolved; absent edges and unique layouts are not established.".into());
                        }
                    }
                }
            }
        }
    }
    let graph_complete = r.complete;
    let n = inputs.len();
    let mut contained = BTreeSet::new();
    for e in &r.overlaps {
        if e.accepted
            && e.kind != OverlapKind::Extension
            && prep[e.from].as_ref().unwrap().len() >= prep[e.to].as_ref().unwrap().len()
            && (prep[e.from].as_ref().unwrap().len() > prep[e.to].as_ref().unwrap().len()
                || e.from < e.to)
            && !r.overlaps.iter().any(|x| {
                x.kind == OverlapKind::Extension
                    && (x.from == e.to || x.to == e.to)
                    && potential(x, &prep, o)
                    && x.score > e.score
            })
        {
            contained.insert(e.to);
        }
    }
    let mut outgoing = vec![vec![]; n * 2];
    let mut incoming = vec![vec![]; n * 2];
    let mut ambiguous = BTreeSet::new();
    for (i, p) in prep.iter().enumerate() {
        if let Some(p) = p {
            let text = p.text(false);
            let mut pi = vec![0usize; text.len()];
            for j in 1..text.len() {
                let mut k = pi[j - 1];
                while k > 0 && text[j] != text[k] {
                    k = pi[k - 1]
                }
                if text[j] == text[k] {
                    k += 1
                }
                pi[j] = k;
            }
            if pi.last().copied().unwrap_or(0) >= o.min_overlap {
                ambiguous.insert(i);
                r.reads[i].warnings.push("A proper repeated prefix/suffix meets the overlap length; repeat-copy placement is not inferred from a full-length best match alone.".into());
            }
        }
    }
    for (i, e) in r.overlaps.iter().enumerate() {
        if !potential(e, &prep, o) {
            continue;
        }
        if e.kind != OverlapKind::Extension {
            if e.optimal_paths_saturated > 1 {
                ambiguous.insert(e.to);
            }
            continue;
        }
        if contained.contains(&e.from) || contained.contains(&e.to) {
            continue;
        }
        outgoing[node(e.from, e.from_reverse)].push(i);
        incoming[node(e.to, e.to_reverse)].push(i);
    }
    let mut next = vec![None; n * 2];
    let mut previous = vec![None; n * 2];
    for v in 0..n * 2 {
        let a = choose(&outgoing[v], &r.overlaps);
        if !outgoing[v].is_empty() && a.is_none() {
            ambiguous.insert(v / 2);
            for e in &outgoing[v] {
                ambiguous.insert(r.overlaps[*e].to);
            }
        }
        if let Some(e) = a {
            let to = node(r.overlaps[e].to, r.overlaps[e].to_reverse);
            if choose(&incoming[to], &r.overlaps) == Some(e) {
                next[v] = Some(e);
                previous[to] = Some(e);
            }
        }
        if !incoming[v].is_empty() && choose(&incoming[v], &r.overlaps).is_none() {
            ambiguous.insert(v / 2);
            for e in &incoming[v] {
                ambiguous.insert(r.overlaps[*e].from);
            }
        }
    }
    if o.both_strands {
        let old = next.clone();
        for (v, entry) in old.iter().enumerate() {
            if let Some(eid) = entry {
                let e = &r.overlaps[*eid];
                let mirror = node(e.to, !e.to_reverse);
                let good = old[mirror].is_some_and(|m| {
                    let x = &r.overlaps[m];
                    x.to == e.from && x.to_reverse == !e.from_reverse && x.score == e.score
                });
                if !good {
                    next[v] = None;
                    let to = node(e.to, e.to_reverse);
                    if previous[to] == Some(*eid) {
                        previous[to] = None;
                    }
                    ambiguous.insert(e.from);
                    ambiguous.insert(e.to);
                }
            }
        }
    }
    let mut used = BTreeSet::new();
    let mut layouts: Vec<(Layout, bool)> = vec![];
    for start in 0..n * 2 {
        let read = start / 2;
        if prep[read].is_none()
            || contained.contains(&read)
            || used.contains(&read)
            || previous[start].is_some()
            || (!o.both_strands && start % 2 == 1)
        {
            continue;
        }
        let mut layout = Layout::root(read, prep[read].as_ref().unwrap(), start % 2 == 1);
        let mut path = BTreeSet::from([read]);
        let mut v = start;
        let mut unique = graph_complete && !ambiguous.contains(&read);
        while let Some(eid) = next[v] {
            let e = r.overlaps[eid].clone();
            if path.contains(&e.to) || used.contains(&e.to) {
                unique = false;
                layout.warnings.push("Selected overlap path revisits an already placed read; no repeat/cycle collapse was forced.".into());
                break;
            }
            let mut candidate = layout.clone();
            match candidate.merge(&e, &prep, o) {
                Ok(()) => {
                    layout = candidate;
                    path.insert(e.to);
                    unique &= !ambiguous.contains(&e.to);
                    r.overlaps[eid].selected = true;
                    v = node(e.to, e.to_reverse);
                }
                Err(error) => {
                    unique = false;
                    r.overlaps[eid].reason = Some(format!(
                        "Pairwise overlap retained, but layout was not forced: {error}"
                    ));
                    layout
                        .warnings
                        .push(format!("A gapped layout join remained unresolved: {error}"));
                    break;
                }
            }
        }
        used.extend(path);
        layouts.push((layout, unique));
    }
    // Remaining noncontained reads belong to cycles, excluded branches or unsupported layouts.
    for i in 0..n {
        if prep[i].is_none() || contained.contains(&i) || used.contains(&i) {
            continue;
        }
        let mut layout = Layout::root(i, prep[i].as_ref().unwrap(), false);
        let cycle = next[node(i, false)].is_some() || next[node(i, true)].is_some();
        let unique = graph_complete && !ambiguous.contains(&i) && !cycle;
        if cycle {
            layout.warnings.push("Cycle or conflicting selected orientation: explicit circular/repeat assembly is not inferred.".into());
        }
        add_unassembled(
            &mut r,
            i,
            if unique {
                "No accepted overlap connects this singleton under the declared model."
            } else {
                "Branch, cycle or alternative placement retained as an unresolved singleton."
            },
        );
        used.insert(i);
        layouts.push((layout, unique));
    }
    let mut pending = contained.clone();
    loop {
        let mut progressed = false;
        for read in pending.clone() {
            let mut choices: Vec<(usize, usize, Layout, Vec<(usize, bool, Vec<usize>)>)> = vec![];
            for (li, (layout, _)) in layouts.iter().enumerate() {
                for e in r
                    .overlaps
                    .iter()
                    .filter(|e| e.accepted && e.kind != OverlapKind::Extension && e.to == read)
                {
                    if layout
                        .reads
                        .get(&e.from)
                        .is_some_and(|m| m.reverse == e.from_reverse)
                    {
                        let mut copy = layout.clone();
                        if copy.merge(e, &prep, o).is_ok() {
                            let signature = copy
                                .reads
                                .iter()
                                .map(|(r, m)| (*r, m.reverse, m.columns.clone()))
                                .collect();
                            choices.push((li, e.id, copy, signature));
                        }
                    }
                }
            }
            if !choices.is_empty() {
                let first = (choices[0].0, choices[0].3.clone());
                if choices.iter().all(|c| c.0 == first.0 && c.3 == first.1) {
                    let (li, eid, layout, _) = choices.remove(0);
                    layouts[li].0 = layout;
                    layouts[li].1 &= !ambiguous.contains(&read);
                    r.overlaps[eid].selected = true;
                    pending.remove(&read);
                    used.insert(read);
                    progressed = true;
                }
            }
        }
        if !progressed {
            break;
        }
    }
    for read in pending {
        let mut layout = Layout::root(read, prep[read].as_ref().unwrap(), false);
        layout.warnings.push("Contained/duplicate read has multiple or unsupported call-column placements; no arbitrary placement is asserted.".into());
        add_unassembled(
            &mut r,
            read,
            "Contained/duplicate placement remained ambiguous or conflicted with existing gaps.",
        );
        used.insert(read);
        layouts.push((layout, false));
    }
    for (layout, unique) in layouts {
        if layout.reads.len() == 1 {
            add_unassembled(
                &mut r,
                *layout.reads.keys().next().unwrap(),
                "Singleton without an unambiguous accepted extension.",
            );
        }
        let id = r.contigs.len();
        r.contigs.push(materialize(
            id,
            &layout,
            &r.reads,
            &prep,
            o,
            unique,
            graph_complete,
        ));
    }
    if !graph_complete {
        for c in &mut r.contigs {
            c.complete = false;
            c.layout_unique = false;
        }
    }
    validate_evidence(&r, &prep)?;
    attach_products(&mut r)?;
    while bounded(&r).is_err() {
        if let Some(c) = r.contigs.pop() {
            r.complete = false;
            if !r.limits.iter().any(|s| s.contains("Report output")) {
                r.limits.push("Report output budget reached; some layouts are omitted and their reads remain explicitly unassembled.".into());
            }
            let removed = c
                .placements
                .iter()
                .map(|p| p.read_index)
                .collect::<BTreeSet<_>>();
            for read in &removed {
                add_unassembled(
                    &mut r,
                    *read,
                    "Layout was omitted at the report-output budget.",
                );
            }
            for e in &mut r.overlaps {
                if removed.contains(&e.from) || removed.contains(&e.to) {
                    e.selected = false;
                }
            }
        } else {
            return Err(Error::TooLarge);
        }
    }
    r.validate()?;
    Ok(r)
}
/// A separately written vote audit: group observations first, then classify support.
fn audit_vote(
    c: &ColumnData,
    reads: &[ReadReport],
    placements: &BTreeMap<usize, Mapping>,
    o: &ContigOptions,
) -> ([usize; 4], Option<char>, ConsensusState) {
    let mut votes: BTreeMap<String, (BTreeSet<char>, bool, BTreeSet<bool>)> = BTreeMap::new();
    for x in &c.observations {
        if x.eligible && b"ACGT".contains(&(x.base as u8)) {
            let v = votes.entry(groups(reads, x.read_index)).or_default();
            v.0.insert(x.base);
            v.2.insert(placements[&x.read_index].reverse ^ reads[x.read_index].source_reverse);
        }
    }
    for read in &c.eligible_gaps {
        let v = votes.entry(groups(reads, *read)).or_default();
        v.1 = true;
        v.2.insert(placements[read].reverse ^ reads[*read].source_reverse);
    }
    let mut counts = [0; 4];
    let mut gaps = 0;
    let mut dirs: BTreeSet<bool> = BTreeSet::new();
    for (bases, gap, orientation) in votes.values() {
        for b in bases {
            counts[b"ACGT".iter().position(|x| *x == *b as u8).unwrap()] += 1;
        }
        gaps += usize::from(*gap);
        if (!bases.is_empty() || *gap) && orientation.len() == 1 {
            dirs.extend(orientation.iter().copied());
        }
    }
    let present = counts.iter().filter(|n| **n > 0).count();
    if present > 1 {
        return (counts, Some('N'), ConsensusState::ConflictingCalls);
    }
    if present == 1 && gaps > 0 {
        return (counts, Some('N'), ConsensusState::GapConflict);
    }
    let support = if present == 1 {
        *counts.iter().max().unwrap()
    } else {
        gaps
    };
    if support >= o.min_support && (!o.require_both_strands || dirs.len() == 2) {
        if present == 1 {
            return (
                counts,
                Some(b"ACGT"[counts.iter().position(|n| *n > 0).unwrap()] as char),
                ConsensusState::SupportedBase,
            );
        }
        return (counts, None, ConsensusState::SupportedGap);
    }
    if support > 0 {
        (counts, Some('N'), ConsensusState::InsufficientSupport)
    } else {
        (counts, Some('N'), ConsensusState::LowQualityOrAmbiguousOnly)
    }
}
fn validate_evidence(r: &ContigReport, prep: &[Option<Prepared>]) -> Result<(), Error> {
    let mut used = BTreeSet::new();
    for c in &r.contigs {
        let mut layout = Layout {
            columns: vec![ColumnData::default(); c.columns.len()],
            reads: BTreeMap::new(),
            warnings: vec![],
        };
        for place in &c.placements {
            if !used.insert(place.read_index) {
                return Err(bad("read occurs in multiple emitted contig layouts"));
            }
            let p = prep
                .get(place.read_index)
                .and_then(Option::as_ref)
                .ok_or_else(|| bad("unprepared read placed"))?;
            if place.calls.len() != p.len() {
                return Err(bad("not every retained source call has a column mapping"));
            }
            let mut map = vec![];
            for (at, call) in place.calls.iter().enumerate() {
                let expected = p.observation(place.read_index, at, place.reverse);
                if expected.stored_call_index >= p.source_calls
                    || call.call_index != expected.call_index
                    || call.stored_call_index != expected.stored_call_index
                    || call.column >= c.columns.len()
                {
                    return Err(bad(
                        "read call mapping differs from the exact source evidence",
                    ));
                }
                layout.columns[call.column].observations.push(expected);
                map.push(call.column);
            }
            if map.windows(2).any(|w| w[0] >= w[1])
                || place.start != map[0]
                || place.end != map.last().unwrap() + 1
            {
                return Err(bad("nonmonotone layout call mapping"));
            }
            layout.reads.insert(
                place.read_index,
                Mapping {
                    reverse: place.reverse,
                    columns: map,
                },
            );
        }
        for (read, map) in &layout.reads {
            let p = prep[*read].as_ref().unwrap();
            let calls = map.columns.iter().copied().collect::<BTreeSet<_>>();
            for at in map.columns[0]..=*map.columns.last().unwrap() {
                if !calls.contains(&at) {
                    layout.columns[at].gaps.insert(*read);
                    if Layout::flank(map, p, at) {
                        layout.columns[at].eligible_gaps.insert(*read);
                    }
                }
            }
        }
        let mut sequence = String::new();
        for (i, got) in c.columns.iter().enumerate() {
            let expected = &layout.columns[i];
            let mut a = got.observations.clone();
            let mut b = expected.observations.clone();
            a.sort_by_key(|x| (x.read_index, x.call_index));
            b.sort_by_key(|x| (x.read_index, x.call_index));
            if a != b
                || got.gap_reads != expected.gaps.iter().copied().collect::<Vec<_>>()
                || got.eligible_gap_reads
                    != expected.eligible_gaps.iter().copied().collect::<Vec<_>>()
            {
                return Err(bad(
                    "column calls/gap evidence does not reconstruct from source placements",
                ));
            }
            let (counts, base, state) = audit_vote(expected, &r.reads, &layout.reads, &r.options);
            let pos = base.map(|base| {
                let at = sequence.len();
                sequence.push(base);
                at
            });
            if got.layout_index != i
                || got.counts != counts
                || got.base != base
                || got.state != state
                || got.consensus_position != pos
            {
                return Err(bad(
                    "consensus column is not supported by independently recounted source calls",
                ));
            }
        }
        if sequence != c.sequence {
            return Err(bad(
                "consensus DNA differs from the audited observed columns",
            ));
        }
    }
    Ok(())
}
fn identity(r: &ContigReport) -> Result<String, Error> {
    let contigs = r
        .contigs
        .iter()
        .map(|c| {
            (
                &c.id,
                &c.name,
                &c.sequence,
                &c.columns,
                &c.placements,
                c.complete,
                c.layout_unique,
                &c.warnings,
            )
        })
        .collect::<Vec<_>>();
    let value = (
        &r.engine_version,
        &r.batch_revision,
        &r.options_revision,
        &r.options,
        &r.reads,
        &r.overlaps,
        contigs,
        &r.unassembled,
        r.complete,
        r.cells_used,
        &r.limits,
        &r.warnings,
    );
    bounded(&value)?;
    Ok(
        blake3::hash(&serde_json::to_vec(&value).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
pub fn report_revision(report: &ContigReport) -> Result<String, Error> {
    report.validate()?;
    identity(report)
}
fn attach_products(r: &mut ContigReport) -> Result<(), Error> {
    let reads = serde_json::to_string(
        &r.reads
            .iter()
            .map(|read| (&read.id, &read.source_revision))
            .collect::<Vec<_>>(),
    )
    .map_err(|_| Error::TooLarge)?;
    let options = serde_json::to_string(&r.options).map_err(|_| Error::TooLarge)?;
    for c in &mut r.contigs {
        if c.sequence.is_empty() {
            continue;
        }
        let mut m = Molecule::dna(&c.name, &c.sequence, Topology::Linear)?;
        for (k, v) in [
            ("engine", r.engine_version.clone()),
            ("batch_revision", r.batch_revision.clone()),
            ("options_revision", r.options_revision.clone()),
            ("id", c.id.to_string()),
            ("options", options.clone()),
            ("reads", reads.clone()),
            ("complete", c.complete.to_string()),
            ("layout_unique", c.layout_unique.to_string()),
        ] {
            m.metadata.insert(format!("contig:{k}"), v);
        }
        c.molecule = Some(m);
    }
    Ok(())
}
impl ContigReport {
    /// Structural historical validation only; a producer's engine/status has no current authority.
    pub fn validate(&self) -> Result<(), Error> {
        let o = &self.options;
        o.validate()?;
        if self.engine_version.is_empty()
            || self.engine_version.len() > 128
            || !digest(&self.batch_revision)
            || self.options_revision != options_revision(o)?
            || self.reads.is_empty()
            || self.reads.len() > MAX_READS
            || self.overlaps.len() > o.max_edges
            || self.contigs.len() > self.reads.len()
            || self.cells_used > o.max_total_cells
            || self.limits.len() > 16
            || self.warnings.len() > 16
            || self.unassembled.len() > self.reads.len()
        {
            return Err(bad("invalid report envelope"));
        }
        let mut ids = BTreeSet::new();
        for read in &self.reads {
            if read.id.is_empty()
                || read.id.len() > 256
                || !ids.insert(&read.id)
                || read.name.len() > 4096
                || !digest(&read.source_revision)
                || read.source_length == 0
                || read.source_length > MAX_CALLS
                || read.source_call_length > 100_000
                || read.retained_start > read.retained_end
                || read.retained_end > read.source_call_length
                || read.retained_end - read.retained_start != read.source_length
                || read.trim_start > read.trim_end
                || read.trim_end > read.source_length
                || read.eligible_calls > read.trim_end - read.trim_start
                || read.warnings.len() > 64
            {
                return Err(bad("invalid read-evidence descriptor"));
            }
        }
        for (i, e) in self.overlaps.iter().enumerate() {
            if e.id != i
                || e.from >= self.reads.len()
                || e.to >= self.reads.len()
                || e.from == e.to
                || e.from_start > e.from_end
                || e.to_start > e.to_end
                || e.from_end > self.reads[e.from].trim_end - self.reads[e.from].trim_start
                || e.to_end > self.reads[e.to].trim_end - self.reads[e.to].trim_start
                || e.columns != e.reference_aligned.len()
                || e.columns != e.query_aligned.len()
                || e.columns > MAX_CALLS * 2
                || e.optimal_paths_saturated == 0
                || e.optimal_paths_saturated > 2
                || e.certain_matches > e.paired_bases
                || e.paired_bases > e.columns
                || e.error_columns > e.columns
                || !o.both_strands && (e.from_reverse || e.to_reverse)
                || e.accepted && e.optimal_paths_saturated != 1
                || e.selected && !e.accepted
            {
                return Err(bad("invalid overlap descriptor"));
            }
            let a = e.reference_aligned.as_bytes();
            let b = e.query_aligned.as_bytes();
            if a.iter().chain(b).any(|v| !b"ACGTRYSWKMBDHVN-".contains(v))
                || a.iter().zip(b).any(|(a, b)| *a == b'-' && *b == b'-')
                || a.iter().filter(|x| **x != b'-').count() != e.from_end - e.from_start
                || b.iter().filter(|x| **x != b'-').count() != e.to_end - e.to_start
            {
                return Err(bad("invalid overlap call strings"));
            }
        }
        let mut used = BTreeSet::new();
        let mut columns = 0usize;
        for (i, c) in self.contigs.iter().enumerate() {
            columns = columns.saturating_add(c.columns.len());
            if c.id != i
                || c.name.is_empty()
                || c.name.len() > 4096
                || c.columns.is_empty()
                || c.columns.len() > o.max_columns
                || c.sequence.len() > o.max_columns
                || c.placements.is_empty()
                || c.placements.len() > self.reads.len()
                || c.warnings.len() > 64
                || c.complete && !c.layout_unique
            {
                return Err(bad("invalid contig layout"));
            }
            let mut owners = BTreeMap::new();
            for p in &c.placements {
                if p.read_index >= self.reads.len()
                    || !used.insert(p.read_index)
                    || !o.both_strands && p.reverse
                {
                    return Err(bad("invalid or duplicated read placement"));
                }
                let read = &self.reads[p.read_index];
                if p.calls.len() != read.trim_end - read.trim_start
                    || p.calls.is_empty()
                    || p.start != p.calls[0].column
                    || p.end != p.calls.last().unwrap().column + 1
                    || p.end > c.columns.len()
                {
                    return Err(bad("invalid placement interval"));
                }
                let mut last = None;
                for (at, call) in p.calls.iter().enumerate() {
                    let expected = if p.reverse {
                        read.trim_end - 1 - at
                    } else {
                        read.trim_start + at
                    };
                    let stored = if read.source_reverse {
                        read.retained_end - 1 - expected
                    } else {
                        read.retained_start + expected
                    };
                    if call.call_index != expected
                        || call.stored_call_index != stored
                        || call.column >= c.columns.len()
                        || last.is_some_and(|x| x >= call.column)
                        || owners
                            .insert((p.read_index, call.call_index), call.column)
                            .is_some()
                    {
                        return Err(bad("invalid source-call coordinate mapping"));
                    }
                    last = Some(call.column);
                }
            }
            let mut sequence = String::new();
            for (j, col) in c.columns.iter().enumerate() {
                if col.layout_index != j
                    || col.observations.is_empty()
                    || col.observations.len() > c.placements.len()
                    || col.gap_reads.len() > c.placements.len()
                    || col.eligible_gap_reads.len() > col.gap_reads.len()
                {
                    return Err(bad("invalid consensus column shape"));
                }
                let mut seen = BTreeSet::new();
                for obs in &col.observations {
                    if !seen.insert(obs.read_index)
                        || owners.get(&(obs.read_index, obs.call_index)) != Some(&j)
                        || obs.read_index >= self.reads.len()
                        || !obs.base.is_ascii()
                        || !b"ACGTRYSWKMBDHVN".contains(&(obs.base as u8))
                    {
                        return Err(bad("invalid column observation source"));
                    }
                    let read = &self.reads[obs.read_index];
                    let stored = if read.source_reverse {
                        read.retained_end - 1 - obs.call_index
                    } else {
                        read.retained_start + obs.call_index
                    };
                    if obs.stored_call_index != stored
                        || obs.eligible && !b"ACGT".contains(&(obs.base as u8))
                    {
                        return Err(bad("invalid observation evidence"));
                    }
                }
                let gaps = col.gap_reads.iter().copied().collect::<BTreeSet<_>>();
                let eligible = col
                    .eligible_gap_reads
                    .iter()
                    .copied()
                    .collect::<BTreeSet<_>>();
                if gaps.len() != col.gap_reads.len()
                    || eligible.len() != col.eligible_gap_reads.len()
                    || !eligible.is_subset(&gaps)
                    || gaps.iter().any(|r| {
                        seen.contains(r)
                            || !c
                                .placements
                                .iter()
                                .any(|p| p.read_index == *r && p.start <= j && p.end > j)
                    })
                {
                    return Err(bad("invalid gap evidence"));
                }
                match col.state {
                    ConsensusState::SupportedBase => {
                        if col
                            .base
                            .is_none_or(|b| !b.is_ascii() || !b"ACGT".contains(&(b as u8)))
                        {
                            return Err(bad("supported base must be concrete"));
                        }
                    }
                    ConsensusState::SupportedGap => {
                        if col.base.is_some() {
                            return Err(bad("supported gap must not invent a base"));
                        }
                    }
                    _ => {
                        if col.base != Some('N') {
                            return Err(bad("uncertain columns must retain explicit N"));
                        }
                    }
                }
                if c.complete
                    && !matches!(
                        col.state,
                        ConsensusState::SupportedBase | ConsensusState::SupportedGap
                    )
                {
                    return Err(bad("complete contig contains unresolved consensus"));
                }
                let expected = col.base.map(|b| {
                    let n = sequence.len();
                    sequence.push(b);
                    n
                });
                if col.consensus_position != expected {
                    return Err(bad("invalid emitted-consensus coordinate"));
                }
            }
            if sequence != c.sequence {
                return Err(bad("column sequence differs from reported consensus"));
            }
            if let Some(m) = &c.molecule {
                if m.name != c.name
                    || m.sequence.alphabet() != Alphabet::Dna
                    || m.sequence.as_str() != c.sequence
                    || m.topology != Topology::Linear
                    || !m.features.is_empty()
                    || m.original_genbank.is_some()
                    || m.metadata.len() > 16
                    || m.metadata.get("contig:engine") != Some(&self.engine_version)
                    || m.metadata.get("contig:batch_revision") != Some(&self.batch_revision)
                    || m.metadata.get("contig:options_revision") != Some(&self.options_revision)
                    || m.metadata.get("contig:id") != Some(&c.id.to_string())
                    || m.metadata.get("contig:complete") != Some(&c.complete.to_string())
                    || m.metadata.get("contig:layout_unique") != Some(&c.layout_unique.to_string())
                    || m.metadata.get("contig:options")
                        != Some(&serde_json::to_string(&self.options).map_err(|_| Error::TooLarge)?)
                    || m.metadata.get("contig:reads")
                        != Some(
                            &serde_json::to_string(
                                &self
                                    .reads
                                    .iter()
                                    .map(|read| (&read.id, &read.source_revision))
                                    .collect::<Vec<_>>(),
                            )
                            .map_err(|_| Error::TooLarge)?,
                        )
                {
                    return Err(bad("generated consensus molecule binding mismatch"));
                }
                m.validate()?;
            }
        }
        if columns > o.max_columns.saturating_mul(self.reads.len()) {
            return Err(Error::TooLarge);
        }
        let mut unplaced = BTreeSet::new();
        for u in &self.unassembled {
            if u.read_index >= self.reads.len()
                || !unplaced.insert(u.read_index)
                || u.reason.is_empty()
                || u.reason.len() > 4096
            {
                return Err(bad("invalid unassembled read descriptor"));
            }
        }
        bounded(self)
    }
    pub fn validate_bindings_for(
        &self,
        inputs: &[ReadInput],
        options: &ContigOptions,
    ) -> Result<(), Error> {
        let batch = batch_revision(inputs, options)?;
        self.validate()?;
        if self.engine_version != CONTIG_ENGINE_VERSION
            || self.batch_revision != batch
            || self.options != *options
            || self.reads.len() != inputs.len()
            || self.reads.iter().zip(inputs).any(|(r, i)| {
                r.id != i.id
                    || r.source_revision != i.molecule.revision()
                    || r.source_length != i.molecule.sequence.len()
            })
        {
            return Err(bad(
                "report does not match current engine/read/options bindings",
            ));
        }
        Ok(())
    }
}
/// Re-run the native model; a serialized report cannot grant consensus authority.
pub fn consensus_documents(
    inputs: &[ReadInput],
    options: &ContigOptions,
    report: &ContigReport,
) -> Result<Vec<Molecule>, Error> {
    report.validate_bindings_for(inputs, options)?;
    let fresh = assemble(inputs, options)?;
    if &fresh != report {
        return Err(bad(
            "saved assembly differs from fresh source reconstruction",
        ));
    }
    Ok(fresh
        .contigs
        .into_iter()
        .filter_map(|c| c.molecule)
        .collect())
}
