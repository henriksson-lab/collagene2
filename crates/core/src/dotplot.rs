//! Source-bound, bounded window dot plots. No alignment path or biological ancestry is inferred.
use collagene_model::{Alphabet, Error, Location, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, io::Write};

pub const DOTPLOT_ENGINE_VERSION: &str = "dotplot-1";
const MAX_SOURCE_BYTES: usize = 16 * 1024 * 1024;
const MAX_SEQUENCE: usize = 10_000_000;
const MAX_INDEX_BYTES: usize = 96 * 1024 * 1024;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DotPlotOptions {
    pub window: usize,
    pub min_matches: usize,
    pub stride: usize,
    pub both_strands: bool,
    pub wrap_circular: bool,
    pub max_work: u64,
    pub max_hits: usize,
    pub max_index_entries: usize,
}
impl Default for DotPlotOptions {
    fn default() -> Self {
        Self {
            window: 15,
            min_matches: 15,
            stride: 1,
            both_strands: true,
            wrap_circular: true,
            max_work: 20_000_000,
            max_hits: 100_000,
            max_index_entries: 1_000_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DotHit {
    pub id: usize,
    pub reference_start: usize,
    /// Original-query coordinate aligned to the first reference residue; decreases on reverse hits.
    pub query_start: usize,
    pub reverse: bool,
    pub window: usize,
    pub compatible_matches: usize,
    pub certain_matches: usize,
    /// The match threshold is met by definite singleton identities alone.
    pub certain: bool,
    pub reference_spans: Vec<Span>,
    pub query_spans: Vec<Span>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DotPlotReport {
    pub engine_version: String,
    pub reference_revision: String,
    pub query_revision: String,
    pub options_revision: String,
    pub reference_length: usize,
    pub query_length: usize,
    pub alphabet: Alphabet,
    pub reference_circular: bool,
    pub query_circular: bool,
    pub options: DotPlotOptions,
    pub hits: Vec<DotHit>,
    pub strategy: String,
    pub work_used: u64,
    pub examined_windows: u64,
    pub total_windows: u64,
    /// Every grid window was classified and all qualifying hits were retained.
    pub grid_complete: bool,
    /// Every unsampled-domain window was classified and every qualifying hit retained.
    pub complete: bool,
    pub sampled: bool,
    pub limits: Vec<String>,
    pub warnings: Vec<String>,
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("dot plot: {s}"))
}
impl DotPlotOptions {
    pub fn validate(&self) -> Result<(), Error> {
        if self.window == 0
            || self.window > 1024
            || self.min_matches == 0
            || self.min_matches > self.window
            || self.stride == 0
            || self.stride > 1_000_000
            || self.max_work == 0
            || self.max_work > 200_000_000
            || self.max_hits == 0
            || self.max_hits > 100_000
            || self.max_index_entries == 0
            || self.max_index_entries > 1_000_000
        {
            return Err(bad(
                "window1..1024, matches1..window, stride1..1000000, work1..200000000, hits1..100000 and index entries1..1000000 are required",
            ));
        }
        Ok(())
    }
}
struct Counter(usize);
impl Write for Counter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > MAX_SOURCE_BYTES.saturating_sub(self.0) {
            return Err(std::io::Error::other("dot plot evidence limit"));
        }
        self.0 += b.len();
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded<T: Serialize>(v: &T) -> Result<(), Error> {
    serde_json::to_writer(Counter(0), v).map_err(|_| Error::TooLarge)
}
fn location(l: &Location, depth: usize, count: &mut usize) -> Result<(), Error> {
    *count = count.saturating_add(1);
    if depth > 32 || *count > 100_000 {
        return Err(Error::TooLarge);
    }
    match l {
        Location::Complement(v) => location(v, depth + 1, count)?,
        Location::Join(v) | Location::Order(v) => {
            for x in v {
                location(x, depth + 1, count)?
            }
        }
        _ => (),
    }
    Ok(())
}
fn source_bound(m: &Molecule) -> Result<(), Error> {
    if m.sequence.len() > MAX_SEQUENCE
        || m.features.len() > 100_000
        || m.metadata.len() > 4096
        || m.name.len() > 4096
    {
        return Err(Error::TooLarge);
    }
    let mut n = 0;
    for f in &m.features {
        if f.qualifiers.len() > 100_000 {
            return Err(Error::TooLarge);
        }
        location(&f.location, 0, &mut n)?;
    }
    bounded(m)?;
    m.validate()
}
fn preflight(reference: &Molecule, query: &Molecule, o: &DotPlotOptions) -> Result<(), Error> {
    o.validate()?;
    source_bound(reference)?;
    source_bound(query)?;
    if reference.sequence.alphabet() != query.sequence.alphabet() {
        return Err(bad(
            "reference and query must have the same explicitly declared alphabet",
        ));
    }
    if reference.sequence.alphabet() == Alphabet::Protein && o.both_strands {
        return Err(bad(
            "protein dot plots do not have a reverse-complement strand; disable both strands",
        ));
    }
    if o.window > reference.sequence.len() || o.window > query.sequence.len() {
        return Err(bad(
            "window exceeds a source length; repeated circular laps are not modeled",
        ));
    }
    Ok(())
}
pub fn source_revisions(
    reference: &Molecule,
    query: &Molecule,
    o: &DotPlotOptions,
) -> Result<(String, String), Error> {
    preflight(reference, query, o)?;
    Ok((reference.revision(), query.revision()))
}
pub fn options_revision(o: &DotPlotOptions) -> Result<String, Error> {
    o.validate()?;
    Ok(
        blake3::hash(&serde_json::to_vec(o).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
fn mask(b: u8, alphabet: Alphabet) -> u32 {
    if alphabet != Alphabet::Protein {
        return u32::from(crate::analysis::base_mask(b));
    }
    const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWYOU*";
    if let Some(i) = AA.iter().position(|a| *a == b) {
        return 1 << i;
    }
    match b {
        b'B' => mask(b'D', alphabet) | mask(b'N', alphabet),
        b'Z' => mask(b'E', alphabet) | mask(b'Q', alphabet),
        b'J' => mask(b'I', alphabet) | mask(b'L', alphabet),
        b'X' => (1 << 22) - 1,
        _ => 0,
    }
}
fn complement(b: u8, alphabet: Alphabet) -> u8 {
    match b {
        b'A' => {
            if alphabet == Alphabet::Rna {
                b'U'
            } else {
                b'T'
            }
        }
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'S' | b'W' | b'N' => b,
        _ => b,
    }
}
fn query_at(q: &[u8], oriented: usize, reverse: bool, alphabet: Alphabet) -> u8 {
    let i = oriented % q.len();
    if reverse {
        complement(q[q.len() - 1 - i], alphabet)
    } else {
        q[i]
    }
}
fn starts(n: usize, circular: bool, w: usize) -> usize {
    if circular { n } else { n - w + 1 }
}
fn grid_count(n: usize, stride: usize) -> u64 {
    n.div_ceil(stride) as u64
}
fn spans(n: usize, start: usize, length: usize, reverse: bool) -> Vec<Span> {
    if !reverse {
        let first = length.min(n - start);
        let mut out = vec![Span {
            start,
            end: start + first,
        }];
        if first < length {
            out.push(Span {
                start: 0,
                end: length - first,
            });
        }
        out
    } else {
        let first = length.min(start + 1);
        let mut out = vec![Span {
            start: start + 1 - first,
            end: start + 1,
        }];
        if first < length {
            out.push(Span {
                start: n - (length - first),
                end: n,
            });
        }
        out
    }
}
struct Work {
    used: u64,
    limit: u64,
}
impl Work {
    fn take(&mut self, n: u64) -> Result<(), ()> {
        if n > self.limit.saturating_sub(self.used) {
            Err(())
        } else {
            self.used += n;
            Ok(())
        }
    }
}
fn word_hash(mut at: impl FnMut(usize) -> u8, w: usize, work: &mut Work) -> Result<u64, ()> {
    work.take(w as u64)?;
    let mut h = 0xcbf29ce484222325u64;
    for i in 0..w {
        h ^= u64::from(at(i));
        h = h.wrapping_mul(0x100000001b3);
    }
    Ok(h)
}
fn hit(
    reference: &Molecule,
    query: &Molecule,
    o: &DotPlotOptions,
    x: usize,
    y: usize,
    reverse: bool,
    work: &mut Work,
) -> Result<Option<DotHit>, ()> {
    let a = reference.sequence.as_str().as_bytes();
    let b = query.sequence.as_str().as_bytes();
    let alphabet = reference.sequence.alphabet();
    let mut compatible = 0;
    let mut certain = 0;
    for i in 0..o.window {
        work.take(1)?;
        let u = mask(a[(x + i) % a.len()], alphabet);
        let v = mask(query_at(b, y + i, reverse, alphabet), alphabet);
        if u & v != 0 {
            compatible += 1;
            if u == v && u.count_ones() == 1 {
                certain += 1;
            }
        }
        if compatible + o.window - i - 1 < o.min_matches {
            return Ok(None);
        }
    }
    if compatible < o.min_matches {
        return Ok(None);
    }
    let original = if reverse { b.len() - 1 - y } else { y };
    Ok(Some(DotHit {
        id: 0,
        reference_start: x,
        query_start: original,
        reverse,
        window: o.window,
        compatible_matches: compatible,
        certain_matches: certain,
        certain: certain >= o.min_matches,
        reference_spans: spans(a.len(), x, o.window, false),
        query_spans: spans(b.len(), original, o.window, reverse),
    }))
}
fn finish(r: &mut DotPlotReport, work: Work) {
    r.work_used = work.used;
    r.complete = r.grid_complete && !r.sampled;
    r.hits
        .sort_by_key(|h| (h.reference_start, h.query_start, h.reverse));
    for (i, h) in r.hits.iter_mut().enumerate() {
        h.id = i;
    }
}
fn add(r: &mut DotPlotReport, h: DotHit, output_bytes: &mut usize) -> bool {
    let mut count = Counter(0);
    if serde_json::to_writer(&mut count, &h).is_err()
        || output_bytes
            .saturating_add(count.0)
            .saturating_add(32 * 1024 + 6)
            > MAX_SOURCE_BYTES
    {
        r.grid_complete = false;
        r.limits
            .push("Output JSON budget exhausted; further matching windows may be omitted.".into());
        return false;
    }
    *output_bytes = output_bytes.saturating_add(count.0 + 6);
    if r.hits.len() >= r.options.max_hits {
        r.grid_complete = false;
        r.limits
            .push("Hit budget exhausted; further matching windows may be omitted.".into());
        false
    } else {
        r.hits.push(h);
        true
    }
}
pub fn analyze(
    reference: &Molecule,
    query: &Molecule,
    o: &DotPlotOptions,
) -> Result<DotPlotReport, Error> {
    let (reference_revision, query_revision) = source_revisions(reference, query, o)?;
    let a = reference.sequence.as_str().as_bytes();
    let b = query.sequence.as_str().as_bytes();
    let alphabet = reference.sequence.alphabet();
    let rc = o.wrap_circular && reference.topology == Topology::Circular;
    let qc = o.wrap_circular && query.topology == Topology::Circular;
    let nx = starts(a.len(), rc, o.window);
    let ny = starts(b.len(), qc, o.window);
    let orientations = if o.both_strands { 2 } else { 1 };
    let mut r=DotPlotReport{engine_version:DOTPLOT_ENGINE_VERSION.into(),reference_revision,query_revision,options_revision:options_revision(o)?,reference_length:a.len(),query_length:b.len(),alphabet,reference_circular:rc,query_circular:qc,options:o.clone(),hits:vec![],strategy:String::new(),work_used:0,examined_windows:0,total_windows:grid_count(nx,o.stride)*grid_count(ny,o.stride)*orientations,grid_complete:true,complete:false,sampled:o.stride>1,limits:vec![],warnings:vec!["Window similarity is not an alignment path, assembly, coverage proof or experimental verification.".into(),"Compatible matches include unresolved identities; certain matches require identical singleton residues on both inputs.".into()]};
    if r.sampled {
        r.warnings.push("Stride samples reference and oriented-query window starts. Complete refers to the full unsampled domain; grid_complete refers only to the declared sampled grid.".into());
    }
    if !o.wrap_circular
        && (reference.topology == Topology::Circular || query.topology == Topology::Circular)
    {
        r.warnings.push("Circular records are treated as linearized at their stored origins; crossing-origin windows are not included.".into());
    }
    let exact = o.min_matches == o.window
        && a.iter()
            .chain(b)
            .all(|v| mask(*v, alphabet).count_ones() == 1);
    let mut work = Work {
        used: 0,
        limit: o.max_work,
    };
    let mut output_bytes = 0usize;
    if exact {
        r.strategy = "exact-concrete-word-index/full-window-collision-check".into();
        let mut index: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        let mut entries = 0;
        let mut bytes = 0usize;
        for x in (0..nx).step_by(o.stride) {
            if entries >= o.max_index_entries || bytes.saturating_add(104) > MAX_INDEX_BYTES {
                r.grid_complete = false;
                r.limits.push("Reference index entry/memory budget exhausted; no partial index is used to imply missing matches.".into());
                finish(&mut r, work);
                return Ok(r);
            }
            let key = match word_hash(|i| a[(x + i) % a.len()], o.window, &mut work) {
                Ok(k) => k,
                Err(()) => {
                    r.grid_complete = false;
                    r.limits.push(
                        "Work budget exhausted while constructing the reference index.".into(),
                    );
                    finish(&mut r, work);
                    return Ok(r);
                }
            };
            let fresh = !index.contains_key(&key);
            index.entry(key).or_default().push(x);
            bytes = bytes.saturating_add(if fresh { 96 } else { 16 });
            entries += 1;
        }
        'orientations: for orientation in 0..orientations {
            let reverse = orientation == 1;
            for y in (0..ny).step_by(o.stride) {
                let key = match word_hash(
                    |i| query_at(b, y + i, reverse, alphabet),
                    o.window,
                    &mut work,
                ) {
                    Ok(k) => k,
                    Err(()) => {
                        r.grid_complete = false;
                        r.limits
                            .push("Work budget exhausted while scanning query words.".into());
                        break 'orientations;
                    }
                };
                // Different hashes prove nonidentity for singleton words. Matching
                // buckets are individually checked, including during truncated rows.
                let candidates = index.get(&key);
                r.examined_windows +=
                    grid_count(nx, o.stride) - candidates.map_or(0, |xs| xs.len() as u64);
                if let Some(xs) = candidates {
                    for x in xs {
                        match hit(reference, query, o, *x, y, reverse, &mut work) {
                            Ok(Some(h)) => {
                                r.examined_windows += 1;
                                if !add(&mut r, h, &mut output_bytes) {
                                    break 'orientations;
                                }
                            }
                            Ok(None) => r.examined_windows += 1,
                            Err(()) => {
                                r.grid_complete = false;
                                r.limits.push(
                                    "Work budget exhausted while checking indexed candidates."
                                        .into(),
                                );
                                break 'orientations;
                            }
                        }
                    }
                }
            }
        }
    } else {
        r.strategy = "exhaustive-window-grid/compatible-residue-masks".into();
        'orientations: for orientation in 0..orientations {
            let reverse = orientation == 1;
            for x in (0..nx).step_by(o.stride) {
                for y in (0..ny).step_by(o.stride) {
                    match hit(reference, query, o, x, y, reverse, &mut work) {
                        Ok(Some(h)) => {
                            r.examined_windows += 1;
                            if !add(&mut r, h, &mut output_bytes) {
                                break 'orientations;
                            }
                        }
                        Ok(None) => r.examined_windows += 1,
                        Err(()) => {
                            r.grid_complete = false;
                            r.limits.push("Work budget exhausted; unexamined windows may contain additional hits.".into());
                            break 'orientations;
                        }
                    }
                }
            }
        }
    }
    finish(&mut r, work);
    Ok(r)
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
}
impl DotPlotReport {
    pub fn validate(&self) -> Result<(), Error> {
        let o = &self.options;
        o.validate()?;
        if ![
            "exact-concrete-word-index/full-window-collision-check",
            "exhaustive-window-grid/compatible-residue-masks",
        ]
        .contains(&self.strategy.as_str())
            || self
                .limits
                .iter()
                .chain(&self.warnings)
                .any(|s| s.len() > 4096)
        {
            return Err(bad("invalid strategy or report text bounds"));
        }

        if self.engine_version != DOTPLOT_ENGINE_VERSION
            || ![
                &self.reference_revision,
                &self.query_revision,
                &self.options_revision,
            ]
            .iter()
            .all(|v| digest(v))
            || self.options_revision != options_revision(o)?
            || self.reference_length == 0
            || self.query_length == 0
            || self.reference_length > MAX_SEQUENCE
            || self.query_length > MAX_SEQUENCE
            || o.window > self.reference_length
            || o.window > self.query_length
            || self.hits.len() > o.max_hits
            || self.work_used > o.max_work
            || self.examined_windows > self.total_windows
            || self.hits.len() as u64 > self.examined_windows
            || (!o.wrap_circular && (self.reference_circular || self.query_circular))
            || self.sampled != (o.stride > 1)
            || self.complete != (self.grid_complete && !self.sampled)
            || self.limits.len() > 8
            || self.warnings.len() > 8
            || self.grid_complete && !self.limits.is_empty()
            || self.alphabet == Alphabet::Protein && o.both_strands
        {
            return Err(bad("invalid report envelope"));
        }
        let total = grid_count(
            starts(self.reference_length, self.reference_circular, o.window),
            o.stride,
        ) * grid_count(
            starts(self.query_length, self.query_circular, o.window),
            o.stride,
        ) * if o.both_strands { 2 } else { 1 };
        if self.total_windows != total || self.grid_complete && self.examined_windows != total {
            return Err(bad("invalid report window coverage"));
        }
        let mut previous = None;
        for (i, h) in self.hits.iter().enumerate() {
            let oriented_y = if h.reverse {
                self.query_length - 1 - h.query_start.min(self.query_length - 1)
            } else {
                h.query_start
            };
            let key = (h.reference_start, h.query_start, h.reverse);
            if h.id != i
                || h.reference_start
                    >= starts(self.reference_length, self.reference_circular, o.window)
                || oriented_y >= starts(self.query_length, self.query_circular, o.window)
                || h.query_start >= self.query_length
                || h.reference_start % o.stride != 0
                || oriented_y % o.stride != 0
                || h.reverse && !o.both_strands
                || h.window != o.window
                || h.compatible_matches < o.min_matches
                || h.compatible_matches > o.window
                || h.certain_matches > h.compatible_matches
                || h.certain != (h.certain_matches >= o.min_matches)
                || h.reference_spans
                    != spans(self.reference_length, h.reference_start, o.window, false)
                || h.query_spans != spans(self.query_length, h.query_start, o.window, h.reverse)
                || previous.is_some_and(|p| p >= key)
            {
                return Err(bad("invalid hit coordinates/counts"));
            }
            previous = Some(key);
        }
        bounded(self)
    }
    /// Bounded source/options/coordinate binding; no dot-plot search is repeated here.
    pub fn validate_bindings_for(
        &self,
        reference: &Molecule,
        query: &Molecule,
        o: &DotPlotOptions,
    ) -> Result<(), Error> {
        self.validate()?;
        let (a, b) = source_revisions(reference, query, o)?;
        if self.reference_revision != a
            || self.query_revision != b
            || self.options != *o
            || self.alphabet != reference.sequence.alphabet()
            || self.reference_length != reference.sequence.len()
            || self.query_length != query.sequence.len()
            || self.reference_circular
                != (o.wrap_circular && reference.topology == Topology::Circular)
            || self.query_circular != (o.wrap_circular && query.topology == Topology::Circular)
        {
            return Err(bad("source/query/options binding mismatch"));
        }
        Ok(())
    }
}
/// An SVG export of the same original-coordinate windows. Each origin crossing is split.
pub fn export_svg(report: &DotPlotReport) -> Result<String, Error> {
    report.validate()?;
    let mut out = String::from(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" height=\"950\" viewBox=\"0 0 900 950\"><rect width=\"900\" height=\"950\" fill=\"white\"/><g font-family=\"sans-serif\" fill=\"#23483a\"><text x=\"70\" y=\"28\" font-size=\"19\">Source-bound window dot plot</text>",
    );
    out.push_str(&format!("<text x=\"70\" y=\"50\" font-size=\"12\">window {}, minimum {} compatible; {} hits; complete {}; sampled {}</text><text x=\"70\" y=\"69\" font-size=\"10\">reference {}</text><text x=\"70\" y=\"84\" font-size=\"10\">query {}</text><rect x=\"70\" y=\"110\" width=\"780\" height=\"780\" fill=\"#f7faf7\" stroke=\"#9aafa1\"/>",report.options.window,report.options.min_matches,report.hits.len(),report.complete,report.sampled,report.reference_revision,report.query_revision));
    // Carry the exact reproducibility envelope in the standalone figure too.
    // Metadata is producer evidence, never imported verification authority.
    let metadata = serde_json::json!({
        "engine_version": report.engine_version,
        "reference_revision": report.reference_revision,
        "query_revision": report.query_revision,
        "options_revision": report.options_revision,
        "options": report.options,
        "alphabet": report.alphabet,
        "reference_length": report.reference_length,
        "query_length": report.query_length,
        "reference_circular": report.reference_circular,
        "query_circular": report.query_circular,
        "grid_complete": report.grid_complete,
        "complete": report.complete,
        "sampled": report.sampled,
        "retained_hits": report.hits.len(),
        "examined_windows": report.examined_windows,
        "total_windows": report.total_windows,
        "work_used": report.work_used,
        "authority": "Unverified producer report; recompute from bound sources and options."
    });
    let metadata = serde_json::to_string(&metadata)
        .map_err(|_| Error::TooLarge)?
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;");
    out.push_str(&format!(
        "<metadata id=\"collagene-dotplot\">{metadata}</metadata>"
    ));
    for h in &report.hits {
        let mut offset = 0;
        while offset < h.window {
            let x = (h.reference_start + offset) % report.reference_length;
            let y = if h.reverse {
                (h.query_start + report.query_length - offset) % report.query_length
            } else {
                (h.query_start + offset) % report.query_length
            };
            let until_x = report.reference_length - x;
            let until_y = if h.reverse {
                y + 1
            } else {
                report.query_length - y
            };
            let len = (h.window - offset).min(until_x).min(until_y);
            let x0 = 70. + 780. * x as f64 / report.reference_length as f64;
            let y0 = 110.
                + 780. * if h.reverse { (y + 1) as f64 } else { y as f64 }
                    / report.query_length as f64;
            let x1 = 70. + 780. * (x + len) as f64 / report.reference_length as f64;
            let y1 = 110.
                + 780.
                    * if h.reverse {
                        (y + 1 - len) as f64
                    } else {
                        (y + len) as f64
                    }
                    / report.query_length as f64;
            let color = if !h.certain {
                "#a96818"
            } else if h.reverse {
                "#ad426c"
            } else {
                "#267b61"
            };
            out.push_str(&format!("<path d=\"M{x0:.4},{y0:.4} L{x1:.4},{y1:.4}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"1\"/>"));
            offset += len;
        }
        if out.len() > 16 * 1024 * 1024 {
            return Err(Error::TooLarge);
        }
    }
    let unit = if report.alphabet == Alphabet::Protein {
        "residues"
    } else {
        "bases"
    };
    out.push_str(&format!("<text x=\"70\" y=\"916\" font-size=\"12\">Reference: 0–{} {unit}; query vertical: 0–{} {unit}; reverse hits descend in query coordinates.</text><text x=\"70\" y=\"937\" font-size=\"11\">Green: certain forward; pink: certain reverse; amber: ambiguity-dependent. Incomplete plots do not prove absence.</text></g></svg>",report.reference_length,report.query_length));
    if out.len() > 16 * 1024 * 1024 {
        return Err(Error::TooLarge);
    }
    Ok(out)
}
