//! Bounded synonymous DNA design and an independent sequence/constraint checker.
//! Preferences are declared integer scores, not expression or organism fitness models.
use crate::{
    Alphabet, Error, Location, Molecule, Sequence, Topology, analysis::base_mask, translation,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};

pub const CODON_ENGINE_VERSION: &str = "codon-design-1";
pub const MAX_RESIDUES: usize = 10_000;
const MAX_FRONTIER_BYTES: usize = 16 * 1024 * 1024;
const MAX_SOURCE_BYTES: usize = 4 * 1024 * 1024;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CodonWeight {
    pub codon: String,
    pub weight: u32,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct GcBounds {
    pub minimum_bp: u16,
    pub maximum_bp: u16,
}
impl Default for GcBounds {
    fn default() -> Self {
        Self {
            minimum_bp: 0,
            maximum_bp: 10_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct WindowGc {
    pub window: usize,
    pub bounds: GcBounds,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DesignOptions {
    pub genetic_code: u8,
    pub codon_weights: Vec<CodonWeight>,
    pub forbidden_motifs: Vec<String>,
    pub left_context: String,
    pub right_context: String,
    pub gc: GcBounds,
    pub max_homopolymer: usize,
    pub window_gc: Option<WindowGc>,
    pub beam_width: usize,
    pub max_work: u64,
}
impl Default for DesignOptions {
    fn default() -> Self {
        Self {
            genetic_code: 1,
            codon_weights: vec![],
            forbidden_motifs: vec![],
            left_context: String::new(),
            right_context: String::new(),
            gc: GcBounds::default(),
            max_homopolymer: 12,
            window_gc: None,
            beam_width: 128,
            max_work: 20_000_000,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum DesignStatus {
    CandidateFound,
    Unsatisfiable,
    SearchIncomplete,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CheckReport {
    pub engine_version: String,
    pub source_revision: String,
    pub candidate_revision: String,
    pub options_revision: String,
    pub target_protein: String,
    pub translated_protein: String,
    pub complete: bool,
    pub satisfies: Option<bool>,
    pub gc_count: usize,
    pub coding_length: usize,
    pub score: u64,
    pub violations: Vec<String>,
    pub work_used: u64,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DesignReport {
    pub engine_version: String,
    pub source_revision: String,
    pub options_revision: String,
    pub target_protein: String,
    pub options: DesignOptions,
    pub candidate: Option<Molecule>,
    pub candidate_check: Option<CheckReport>,
    pub score: Option<u64>,
    pub status: DesignStatus,
    pub search_complete: bool,
    pub optimality_proven: bool,
    pub work_used: u64,
    pub pruned_states: usize,
    pub limits: Vec<String>,
    pub warnings: Vec<String>,
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("codon design: {s}"))
}
fn bounds(b: GcBounds) -> Result<(), Error> {
    if b.minimum_bp > b.maximum_bp || b.maximum_bp > 10_000 {
        Err(bad(
            "GC bounds must be ordered integer basis points in0..10000",
        ))
    } else {
        Ok(())
    }
}
fn concrete(s: &str) -> bool {
    s.bytes().all(|b| b"ACGT".contains(&b))
}
fn gc(b: u8) -> bool {
    matches!(b, b'C' | b'G')
}
fn within(count: usize, len: usize, b: GcBounds) -> bool {
    count as u64 * 10_000 >= len as u64 * u64::from(b.minimum_bp)
        && count as u64 * 10_000 <= len as u64 * u64::from(b.maximum_bp)
}
struct CounterWriter(usize, usize);
impl Write for CounterWriter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > self.1.saturating_sub(self.0) {
            return Err(std::io::Error::other("size limit"));
        }
        self.0 += b.len();
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded_json<T: Serialize>(v: &T, max: usize) -> Result<(), Error> {
    serde_json::to_writer(CounterWriter(0, max), v).map_err(|_| Error::TooLarge)
}
impl DesignOptions {
    pub fn validate(&self) -> Result<(), Error> {
        if ![1, 2, 11].contains(&self.genetic_code) {
            return Err(bad("supported genetic codes are1,2,11"));
        }
        bounds(self.gc)?;
        if self.beam_width == 0
            || self.beam_width > 4096
            || self.max_work == 0
            || self.max_work > 100_000_000
            || self.max_homopolymer == 0
            || self.max_homopolymer > 30_000
        {
            return Err(bad(
                "beam1..4096, work1..100000000 and homopolymer1..30000 are required",
            ));
        }
        if self.left_context.len() > 2048
            || self.right_context.len() > 2048
            || !concrete(&self.left_context)
            || !concrete(&self.right_context)
            || self.codon_weights.len() > 64
            || self.forbidden_motifs.len() > 64
        {
            return Err(bad(
                "contexts must be concrete uppercase DNA <=2048 bases; at most64 codon weights/motifs",
            ));
        }
        let mut codons = BTreeSet::new();
        for w in &self.codon_weights {
            if w.codon.len() != 3
                || !concrete(&w.codon)
                || w.weight > 1_000_000
                || !codons.insert(&w.codon)
            {
                return Err(bad(
                    "codon weights require distinct uppercase DNA triplets and integers0..1000000",
                ));
            }
        }
        for m in &self.forbidden_motifs {
            if m.is_empty() || m.len() > 128 || !m.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
            {
                return Err(bad(
                    "forbidden motifs must be uppercase IUPAC DNA of length1..128",
                ));
            }
        }
        if let Some(w) = &self.window_gc {
            bounds(w.bounds)?;
            if w.window == 0 || w.window > 2048 {
                return Err(bad("window length must be1..2048"));
            }
        }
        bounded_json(self, 32 * 1024)
    }
}
pub fn options_revision(o: &DesignOptions) -> Result<String, Error> {
    o.validate()?;
    Ok(
        blake3::hash(&serde_json::to_vec(o).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
fn location_bound(l: &Location, depth: usize, items: &mut usize) -> Result<(), Error> {
    *items += 1;
    if depth > 32 || *items > 100_000 {
        return Err(Error::TooLarge);
    }
    match l {
        Location::Complement(l) => location_bound(l, depth + 1, items)?,
        Location::Join(v) | Location::Order(v) => {
            for l in v {
                location_bound(l, depth + 1, items)?
            }
        }
        _ => (),
    };
    Ok(())
}
fn source_bounds(m: &Molecule) -> Result<(), Error> {
    if m.sequence.len() > MAX_RESIDUES * 3
        || m.features.len() > 4096
        || m.metadata.len() > 4096
        || m.name.len() > 4096
    {
        return Err(Error::TooLarge);
    }
    let mut items = 0usize;
    let mut bytes = m
        .name
        .len()
        .saturating_add(m.sequence.len())
        .saturating_add(m.original_genbank.as_ref().map_or(0, String::len));
    let mut reserve = |n: usize| -> Result<(), Error> {
        items = items.saturating_add(1);
        bytes = bytes.saturating_add(n).saturating_add(64);
        if items > 100_000 || bytes > MAX_SOURCE_BYTES {
            Err(Error::TooLarge)
        } else {
            Ok(())
        }
    };
    reserve(0)?;
    for (k, v) in &m.metadata {
        reserve(k.len().saturating_add(v.len()))?;
    }
    for f in &m.features {
        reserve(f.kind.len())?;
        for (k, v) in &f.qualifiers {
            reserve(k.len().saturating_add(v.as_ref().map_or(0, String::len)))?;
        }
    }
    let mut nodes = 0;
    for f in &m.features {
        location_bound(&f.location, 0, &mut nodes)?;
    }
    bounded_json(m, MAX_SOURCE_BYTES)?;
    m.validate()
}

fn guard_evidence(m: &Molecule, code: u8) -> Result<(), Error> {
    for k in m.metadata.keys() {
        if [
            "protein:modifications",
            "protein:modified",
            "modifications",
            "duplex_model",
            "legacy:unknown_physical_state",
        ]
        .iter()
        .any(|v| k.eq_ignore_ascii_case(v))
        {
            return Err(bad(
                "retained chemical-state evidence is outside the sequence-only coding model",
            ));
        }
    }
    for f in &m.features {
        if [
            "MOD_RES",
            "modified_residue",
            "cross-link",
            "cross_link",
            "bond",
            "LIPID",
            "CARBOHYD",
            "DISULFID",
        ]
        .iter()
        .any(|k| f.kind.eq_ignore_ascii_case(k))
        {
            return Err(bad(
                "modified/cross-linked source is outside the coding model",
            ));
        }
        for (k, v) in &f.qualifiers {
            if [
                "transl_except",
                "exception",
                "ribosomal_slippage",
                "frameshift",
                "pseudo",
                "pseudogene",
                "trans_splicing",
            ]
            .iter()
            .any(|q| k.eq_ignore_ascii_case(q))
            {
                return Err(bad(
                    "retained translation exceptions require an explicit supported model",
                ));
            }
            if k.eq_ignore_ascii_case("transl_table")
                && v.as_ref()
                    .is_none_or(|s| s.trim().parse::<u8>().ok() != Some(code))
            {
                return Err(bad(
                    "retained genetic code conflicts with the selected code",
                ));
            }
            if k.eq_ignore_ascii_case("codon_start") && v.as_ref().is_none_or(|s| s.trim() != "1") {
                return Err(bad("whole-CDS design requires retained codon_start1"));
            }
        }
    }
    Ok(())
}
fn target(m: &Molecule, o: &DesignOptions) -> Result<String, Error> {
    source_bounds(m)?;
    guard_evidence(m, o.genetic_code)?;
    if m.topology != Topology::Linear {
        return Err(bad(
            "select a linear protein or complete linear DNA coding sequence",
        ));
    }
    let p = match m.sequence.alphabet() {
        Alphabet::Protein => m.sequence.as_str().to_string(),
        Alphabet::Dna => {
            if !concrete(m.sequence.as_str()) || m.sequence.len() % 3 != 0 {
                return Err(bad("DNA source must be concrete complete codons in frame0"));
            }
            translation::translate(
                &m.sequence,
                &translation::TranslationOptions {
                    genetic_code: o.genetic_code,
                    ..Default::default()
                },
            )?
            .protein
        }
        Alphabet::Rna => {
            return Err(bad(
                "RNA coding sources are not supported by this DNA design profile",
            ));
        }
    };
    if p == "*" {
        return Err(bad(
            "a required protein must contain at least one actual amino acid",
        ));
    }
    if p.len() > MAX_RESIDUES {
        return Err(Error::TooLarge);
    }
    if p.bytes()
        .enumerate()
        .any(|(i, b)| !b"ACDEFGHIKLMNPQRSTVWY".contains(&b) && !(b == b'*' && i + 1 == p.len()))
    {
        return Err(bad(
            "source must contain canonical residues with at most one terminal stop",
        ));
    }
    Ok(p)
}
/// Preflight a source and derive its content revision without a search/constraint scan.
pub fn source_revision(source: &Molecule, options: &DesignOptions) -> Result<String, Error> {
    options.validate()?;
    target(source, options)?;
    Ok(source.revision())
}

fn warnings(m: &Molecule) -> Vec<String> {
    let mut out=vec!["Ordinary codon decoding only: no initiation recoding, expression, host fitness, mRNA folding or chemical synthesis-success inference.".into(),"Contexts constrain junctions and complete windows but are not included in the returned coding DNA; global GC applies only to coding DNA including any terminal stop codon.".into()];
    if !m.features.is_empty() || !m.metadata.is_empty() || m.original_genbank.is_some() {
        out.push("Source annotations and original evidence remain on the source; their coordinates and chemistry are not copied to designed DNA.".into());
    }
    out
}
fn weight(o: &DesignOptions, c: &str) -> u32 {
    o.codon_weights
        .iter()
        .find(|w| w.codon == c)
        .map_or(1, |w| w.weight)
}
fn codons(code: u8) -> Result<Vec<(String, u8)>, Error> {
    let mut out = vec![];
    for a in b"ACGT" {
        for b in b"ACGT" {
            for c in b"ACGT" {
                let text = String::from_utf8(vec![*a, *b, *c]).unwrap();
                let s = Sequence::dna(&text)?;
                let amino = translation::translate(
                    &s,
                    &translation::TranslationOptions {
                        genetic_code: code,
                        ..Default::default()
                    },
                )?
                .protein
                .as_bytes()[0];
                out.push((text, amino));
            }
        }
    }
    Ok(out)
}
struct Work {
    used: u64,
    limit: u64,
}
impl Work {
    fn take(&mut self, n: u64) -> Result<(), ()> {
        if n > self.limit.saturating_sub(self.used) {
            return Err(());
        }
        self.used += n;
        Ok(())
    }
}
fn patterns(o: &DesignOptions) -> Vec<Vec<u8>> {
    let mut out = vec![];
    for m in &o.forbidden_motifs {
        for text in [
            m.clone(),
            Sequence::dna(m)
                .unwrap()
                .reverse_complement()
                .unwrap()
                .as_str()
                .to_string(),
        ] {
            let p = text.bytes().map(base_mask).collect::<Vec<_>>();
            if !out.contains(&p) {
                out.push(p);
            }
        }
    }
    out
}
fn at(left: &[u8], body: &[u8], right: &[u8], i: usize) -> u8 {
    if i < left.len() {
        left[i]
    } else if i < left.len() + body.len() {
        body[i - left.len()]
    } else {
        right[i - left.len() - body.len()]
    }
}
/// Incremental search predicate: only windows/motifs/runs ending in newly appended bases.
fn append_ok(
    left: &[u8],
    body: &[u8],
    right: &[u8],
    start: usize,
    o: &DesignOptions,
    ps: &[Vec<u8>],
    work: &mut Work,
) -> Result<bool, ()> {
    let n = left.len() + body.len() + right.len();
    for end in start..n {
        work.take(1)?;
        for p in ps {
            if p.len() <= end + 1 {
                let mut hit = true;
                for (j, mask) in p.iter().enumerate() {
                    work.take(1)?;
                    if base_mask(at(left, body, right, end + 1 - p.len() + j)) & mask == 0 {
                        hit = false;
                        break;
                    }
                }
                if hit {
                    return Ok(false);
                }
            }
        }
        let base = at(left, body, right, end);
        let mut run = 1;
        for i in (0..end).rev() {
            work.take(1)?;
            if at(left, body, right, i) != base {
                break;
            }
            run += 1;
            if run > o.max_homopolymer {
                return Ok(false);
            }
        }
        if let Some(w) = &o.window_gc {
            if end + 1 >= w.window {
                let mut count = 0;
                for i in end + 1 - w.window..=end {
                    work.take(1)?;
                    count += usize::from(gc(at(left, body, right, i)));
                }
                if !within(count, w.window, w.bounds) {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}
/// Independent full scan. Does not call the search's incremental constraint predicate.
pub fn check(
    source: &Molecule,
    candidate: &Molecule,
    o: &DesignOptions,
) -> Result<CheckReport, Error> {
    o.validate()?;
    let protein = target(source, o)?;
    source_bounds(candidate)?;
    guard_evidence(candidate, o.genetic_code)?;
    if candidate.sequence.alphabet() != Alphabet::Dna
        || candidate.topology != Topology::Linear
        || !concrete(candidate.sequence.as_str())
    {
        return Err(bad("candidate must be concrete linear DNA"));
    }
    let mut r = CheckReport {
        engine_version: CODON_ENGINE_VERSION.into(),
        source_revision: source.revision(),
        candidate_revision: candidate.revision(),
        options_revision: options_revision(o)?,
        target_protein: protein,
        translated_protein: String::new(),
        complete: true,
        satisfies: None,
        gc_count: 0,
        coding_length: candidate.sequence.len(),
        score: 0,
        violations: vec![],
        work_used: 0,
        warnings: warnings(source),
    };
    let mut work = Work {
        used: 0,
        limit: o.max_work,
    };
    let full = [
        o.left_context.as_bytes(),
        candidate.sequence.as_str().as_bytes(),
        o.right_context.as_bytes(),
    ]
    .concat();
    let scan = (|| -> Result<(), ()> {
        if candidate.sequence.len() % 3 != 0 {
            r.violations
                .push("Candidate has trailing bases outside complete codons.".into());
        }
        work.take(candidate.sequence.len() as u64)?;
        let translation = translation::translate(
            &candidate.sequence,
            &translation::TranslationOptions {
                genetic_code: o.genetic_code,
                ..Default::default()
            },
        )
        .map_err(|_| ())?;
        r.translated_protein = translation.protein;
        if r.translated_protein != r.target_protein {
            r.violations.push(
                "Candidate translation differs from the independently bound required protein."
                    .into(),
            );
        }
        for c in candidate.sequence.as_str().as_bytes().chunks_exact(3) {
            work.take(1)?;
            let w = weight(o, std::str::from_utf8(c).unwrap());
            r.score += u64::from(w);
            if w == 0 && r.violations.len() < 256 {
                r.violations
                    .push("Candidate uses an explicitly disallowed zero-weight codon.".into());
            }
        }
        r.gc_count = candidate
            .sequence
            .as_str()
            .bytes()
            .filter(|b| gc(*b))
            .count();
        if !within(r.gc_count, candidate.sequence.len(), o.gc) {
            r.violations
                .push("Coding GC count is outside the declared global bounds.".into());
        }
        // Separate complete left-to-right motif/run/window traversal.
        for m in patterns(o) {
            for (start, window) in full.windows(m.len()).enumerate() {
                let mut hit = true;
                for (b, mask) in window.iter().zip(&m) {
                    work.take(1)?;
                    if base_mask(*b) & mask == 0 {
                        hit = false;
                        break;
                    }
                }
                if hit && r.violations.len() < 256 {
                    r.violations.push(format!(
                        "Forbidden motif (either strand) at construct offset {start}."
                    ));
                }
            }
        }
        let mut previous = None;
        let mut run = 0;
        for (i, b) in full.iter().enumerate() {
            work.take(1)?;
            run = if previous == Some(*b) { run + 1 } else { 1 };
            previous = Some(*b);
            if run == o.max_homopolymer + 1 && r.violations.len() < 256 {
                r.violations.push(format!(
                    "Homopolymer exceeds its limit at construct offset {i}."
                ));
            }
        }
        if let Some(w) = &o.window_gc {
            let mut count = 0;
            for (i, b) in full.iter().enumerate() {
                work.take(1)?;
                count += usize::from(gc(*b));
                if i >= w.window {
                    count -= usize::from(gc(full[i - w.window]));
                }
                if i + 1 >= w.window
                    && !within(count, w.window, w.bounds)
                    && r.violations.len() < 256
                {
                    r.violations.push(format!(
                        "Window GC outside bounds at construct offset {}.",
                        i + 1 - w.window
                    ));
                }
            }
        }
        Ok(())
    })();
    r.complete = scan.is_ok();
    r.work_used = work.used;
    r.satisfies = r.complete.then_some(r.violations.is_empty());
    if r.violations.len() == 256 {
        r.warnings.push(
            "At most256 violation details are retained; further failures may be omitted.".into(),
        );
    }
    if !r.complete {
        r.warnings
            .push("Checker work budget exhausted; no satisfaction claim is available.".into());
    }
    Ok(r)
}
// Prefixes share compact arena nodes; frontier states retain only constraint lookback.
// Lexicographic ranks preserve full-prefix order without comparing/copying full DNA.
#[derive(Clone, Copy)]
struct Node {
    parent: u32,
    codon: [u8; 3],
}
#[derive(Clone)]
struct State {
    tail: Vec<u8>,
    score: u64,
    gc: usize,
    node: u32,
    lex_rank: usize,
    parent_rank: usize,
    codon: [u8; 3],
}
fn order(a: &State, b: &State) -> std::cmp::Ordering {
    b.score
        .cmp(&a.score)
        .then(a.parent_rank.cmp(&b.parent_rank))
        .then(a.codon.cmp(&b.codon))
}
fn ranking_work(n: usize) -> u64 {
    (n as u64).saturating_mul(usize::BITS.saturating_sub(n.leading_zeros()).max(1) as u64)
}
pub fn design(source: &Molecule, o: &DesignOptions) -> Result<DesignReport, Error> {
    o.validate()?;
    let protein = target(source, o)?;
    let mut r = DesignReport {
        engine_version: CODON_ENGINE_VERSION.into(),
        source_revision: source.revision(),
        options_revision: options_revision(o)?,
        target_protein: protein.clone(),
        options: o.clone(),
        candidate: None,
        candidate_check: None,
        score: None,
        status: DesignStatus::SearchIncomplete,
        search_complete: true,
        optimality_proven: false,
        work_used: 0,
        pruned_states: 0,
        limits: vec![],
        warnings: warnings(source),
    };
    let all = codons(o.genetic_code)?;
    let ps = patterns(o);
    let mut work = Work {
        used: 0,
        limit: o.max_work,
    };
    let lookback = o
        .max_homopolymer
        .max(
            ps.iter()
                .map(|p| p.len().saturating_sub(1))
                .max()
                .unwrap_or(0),
        )
        .max(o.window_gc.as_ref().map_or(0, |w| w.window - 1));
    const MAX_NODES: usize = 1_000_000;
    let mut arena = Vec::<Node>::with_capacity((o.beam_width * protein.len()).min(MAX_NODES));
    let left = o.left_context.as_bytes();
    let tail = left[left.len().saturating_sub(lookback)..].to_vec();
    let mut states = vec![State {
        tail,
        score: 0,
        gc: 0,
        node: u32::MAX,
        lex_rank: 0,
        parent_rank: 0,
        codon: [0; 3],
    }];
    match append_ok(left, b"", b"", 0, o, &ps, &mut work) {
        Ok(true) => (),
        Ok(false) => states.clear(),
        Err(()) => {
            states.clear();
            r.search_complete = false;
            r.limits
                .push("Search work budget exhausted while checking the left context.".into());
        }
    }
    let mut finished = 0;
    for (i, aa) in protein.bytes().enumerate() {
        if states.is_empty() {
            break;
        }
        let new_len = (i + 1) * 3;
        let suffix_len = (left.len() + new_len).min(lookback + 3);
        let capacity = o
            .beam_width
            .min((MAX_FRONTIER_BYTES / 8 / (suffix_len + 128)).max(1));
        if states.len() > capacity {
            states.sort_unstable_by(order);
            r.pruned_states += states.len() - capacity;
            states.truncate(capacity);
            r.search_complete = false;
        }
        let choices = all
            .iter()
            .filter(|(c, a)| *a == aa && weight(o, c) > 0)
            .collect::<Vec<_>>();
        let mut next = Vec::with_capacity((states.len() * choices.len()).min(capacity * 6));
        let mut exhausted = false;
        'states: for s in &states {
            for (c, _) in &choices {
                if work
                    .take(3 + 2 * (s.tail.len() as u64).div_ceil(64))
                    .is_err()
                {
                    exhausted = true;
                    break 'states;
                }
                let count = s.gc + c.bytes().filter(|b| gc(*b)).count();
                let remaining = (protein.len() - i - 1) * 3;
                if count as u64 * 10_000 > protein.len() as u64 * 3 * u64::from(o.gc.maximum_bp)
                    || (count + remaining) as u64 * 10_000
                        < protein.len() as u64 * 3 * u64::from(o.gc.minimum_bp)
                {
                    continue;
                }
                let mut tail = Vec::with_capacity(s.tail.len() + 3);
                tail.extend_from_slice(&s.tail);
                tail.extend_from_slice(c.as_bytes());
                match append_ok(b"", &tail, b"", s.tail.len(), o, &ps, &mut work) {
                    Ok(true) => {
                        if tail.len() > lookback {
                            tail.drain(..tail.len() - lookback);
                        }
                        next.push(State {
                            tail,
                            score: s.score + u64::from(weight(o, c)),
                            gc: count,
                            node: s.node,
                            lex_rank: 0,
                            parent_rank: s.lex_rank,
                            codon: c.as_bytes().try_into().unwrap(),
                        })
                    }
                    Ok(false) => (),
                    Err(()) => {
                        exhausted = true;
                        break 'states;
                    }
                }
            }
        }
        if work.take(ranking_work(next.len())).is_err() {
            exhausted = true;
        }
        next.sort_unstable_by(order);
        if next.len() > capacity {
            r.pruned_states += next.len() - capacity;
            next.truncate(capacity);
            r.search_complete = false;
        }
        if arena.len() + next.len() > MAX_NODES {
            r.search_complete = false;
            r.limits.push(
                "Prefix arena memory limit reached; reduce the beam or target length.".into(),
            );
            states.clear();
            break;
        }
        let mut lex = (0..next.len()).collect::<Vec<_>>();
        lex.sort_unstable_by_key(|i| (next[*i].parent_rank, next[*i].codon));
        if work
            .take(ranking_work(lex.len()) + next.len() as u64)
            .is_err()
        {
            exhausted = true;
        }
        for (rank, index) in lex.into_iter().enumerate() {
            next[index].lex_rank = rank;
        }
        for s in &mut next {
            let node = arena.len() as u32;
            arena.push(Node {
                parent: s.node,
                codon: s.codon,
            });
            s.node = node;
        }
        states = next;
        finished = i + 1;
        if exhausted {
            r.search_complete = false;
            r.limits.push(
                "Search work budget exhausted; feasibility/optimality may remain unresolved."
                    .into(),
            );
            break;
        }
    }
    if r.pruned_states > 0 {
        r.limits.push("Beam or frontier-memory pruning occurred; neither global optimality nor infeasibility is proven.".into());
    }
    let mut valid = vec![];
    if finished == protein.len() {
        for s in states {
            match append_ok(
                b"",
                &s.tail,
                o.right_context.as_bytes(),
                s.tail.len(),
                o,
                &ps,
                &mut work,
            ) {
                Ok(true) => {
                    if within(s.gc, protein.len() * 3, o.gc) {
                        valid.push(s)
                    }
                }
                Ok(false) => (),
                Err(()) => {
                    r.search_complete = false;
                    r.limits.push(
                        "Search work budget exhausted while checking the right context.".into(),
                    );
                    break;
                }
            }
        }
    }
    valid.sort_unstable_by(order);
    if let Some(s) = valid.first() {
        if work.take((protein.len() * 3) as u64).is_err() {
            r.search_complete = false;
            r.limits.push(
                "Search work budget exhausted before reconstructing the final coding DNA.".into(),
            );
        } else {
            let mut path = Vec::with_capacity(protein.len());
            let mut at = s.node;
            while at != u32::MAX {
                let n = arena[at as usize];
                path.push(n.codon);
                at = n.parent;
            }
            let mut sequence = String::with_capacity(protein.len() * 3);
            for c in path.into_iter().rev() {
                sequence.push_str(std::str::from_utf8(&c).unwrap());
            }
            let mut m = Molecule::dna(
                format!(
                    "{} codon design",
                    source.name.chars().take(100).collect::<String>()
                ),
                &sequence,
                Topology::Linear,
            )?;
            m.metadata
                .insert("codon_design:engine".into(), CODON_ENGINE_VERSION.into());
            m.metadata
                .insert("codon_design:source_revision".into(), source.revision());
            m.metadata.insert(
                "codon_design:options_revision".into(),
                r.options_revision.clone(),
            );
            m.metadata
                .insert("codon_design:required_protein".into(), protein);
            m.metadata.insert(
                "codon_design:options".into(),
                serde_json::to_string(o).map_err(|_| Error::TooLarge)?,
            );
            let verified = check(source, &m, o)?;
            if verified.complete && verified.satisfies == Some(true) {
                r.status = DesignStatus::CandidateFound;
                r.score = Some(s.score);
                r.candidate = Some(m);
                r.candidate_check = Some(verified);
                r.optimality_proven = r.search_complete;
            } else {
                r.search_complete = false;
                r.limits.push("Independent checker could not establish all constraints; no candidate is released.".into());
            }
        }
    }
    if r.candidate.is_none() && r.search_complete {
        r.status = DesignStatus::Unsatisfiable;
    }
    r.work_used = work.used;
    Ok(r)
}

fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
}
impl CheckReport {
    pub fn validate(&self) -> Result<(), Error> {
        if self.engine_version != CODON_ENGINE_VERSION
            || ![
                &self.source_revision,
                &self.candidate_revision,
                &self.options_revision,
            ]
            .iter()
            .all(|s| digest(s))
            || self.target_protein.len() > MAX_RESIDUES
            || self.translated_protein.len() > MAX_RESIDUES
            || self.coding_length > MAX_RESIDUES * 3
            || self.gc_count > self.coding_length
            || self.work_used > 100_000_000
            || self.violations.len() > 256
            || self.complete.then_some(self.violations.is_empty()) != self.satisfies
        {
            return Err(bad("invalid checker report"));
        }
        bounded_json(self, 4 * 1024 * 1024)
    }
    /// Envelope/source binding only; does not scan motifs or rerun the checker.
    pub fn validate_bindings_for(
        &self,
        source: &Molecule,
        candidate: &Molecule,
        o: &DesignOptions,
    ) -> Result<(), Error> {
        self.validate()?;
        o.validate()?;
        let required = target(source, o)?;
        source_bounds(candidate)?;
        if self.source_revision != source.revision()
            || self.candidate_revision != candidate.revision()
            || self.options_revision != options_revision(o)?
            || self.target_protein != required
        {
            return Err(bad("checker source/candidate/options binding mismatch"));
        }
        Ok(())
    }
    pub fn validate_for(
        &self,
        source: &Molecule,
        candidate: &Molecule,
        o: &DesignOptions,
    ) -> Result<(), Error> {
        self.validate()?;
        if *self != check(source, candidate, o)? {
            return Err(bad(
                "checker report does not match the bound source/candidate/options",
            ));
        }
        Ok(())
    }
}
impl DesignReport {
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        if self.engine_version != CODON_ENGINE_VERSION
            || !digest(&self.source_revision)
            || self.options_revision != options_revision(&self.options)?
            || self.target_protein.is_empty()
            || self.target_protein.len() > MAX_RESIDUES
            || self.work_used > self.options.max_work
            || (self.search_complete && (self.pruned_states > 0 || !self.limits.is_empty()))
            || self.limits.len() > 10
            || self.optimality_proven != (self.candidate.is_some() && self.search_complete)
            || matches!(self.status, DesignStatus::CandidateFound) != self.candidate.is_some()
            || matches!(self.status, DesignStatus::Unsatisfiable) && !self.search_complete
            || matches!(self.status, DesignStatus::SearchIncomplete) && self.search_complete
        {
            return Err(bad("invalid design report"));
        }
        match (&self.candidate, &self.candidate_check, self.score) {
            (Some(m), Some(c), Some(score)) => {
                source_bounds(m)?;
                c.validate()?;
                if c.source_revision != self.source_revision
                    || c.candidate_revision != m.revision()
                    || c.options_revision != self.options_revision
                    || c.target_protein != self.target_protein
                    || !c.complete
                    || c.satisfies != Some(true)
                    || c.score != score
                {
                    return Err(bad("candidate/check binding mismatch"));
                }
            }
            (None, None, None) => (),
            _ => return Err(bad("incomplete candidate/check binding")),
        };
        bounded_json(self, 8 * 1024 * 1024)
    }
    /// Envelope/source binding only; no search or constraint scan is performed.
    pub fn validate_bindings_for_source(
        &self,
        source: &Molecule,
        o: &DesignOptions,
    ) -> Result<(), Error> {
        self.validate()?;
        o.validate()?;
        let required = target(source, o)?;
        if self.source_revision != source.revision()
            || self.options != *o
            || self.target_protein != required
        {
            return Err(bad("design source/options binding mismatch"));
        }
        if let (Some(m), Some(c)) = (&self.candidate, &self.candidate_check) {
            c.validate_bindings_for(source, m, o)?;
        }
        Ok(())
    }
    pub fn validate_for_source(&self, source: &Molecule, o: &DesignOptions) -> Result<(), Error> {
        self.validate()?;
        o.validate()?;
        let required = target(source, o)?;
        if self.source_revision != source.revision()
            || self.options != *o
            || self.target_protein != required
        {
            return Err(bad("design source/options binding mismatch"));
        }
        if let (Some(m), Some(c)) = (&self.candidate, &self.candidate_check) {
            c.validate_for(source, m, o)?;
        }
        Ok(())
    }
}
