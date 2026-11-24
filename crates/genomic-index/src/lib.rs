//! BAI/CSI exchange and candidate BGZF chunk lookup. An index is an untrusted
//! claim about another file, not verification of that file or its alignments.
use std::{collections::BTreeSet, fmt};
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Error(&'static str);
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.0)
    }
}
impl std::error::Error for Error {}
type Result<T> = std::result::Result<T, Error>;
fn fail<T>(s: &'static str) -> Result<T> {
    Err(Error(s))
}
fn cap(n: usize, max: usize) -> Result<usize> {
    if n > max {
        fail("index resource limit exceeded")
    } else {
        Ok(n)
    }
}
fn add(a: usize, b: usize, max: usize) -> Result<usize> {
    cap(a.checked_add(b).ok_or(Error("size overflow"))?, max)
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Bai,
    Csi,
}
/// Ordered packed BGZF address. Low16 bits are uncompressed block offset;
/// high48 bits are compressed-file byte offset. It is never a byte count.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct VirtualOffset(u64);
impl VirtualOffset {
    pub fn from_raw(raw: u64) -> Self {
        Self(raw)
    }
    pub fn new(compressed: u64, uncompressed: u16) -> Result<Self> {
        if compressed >= 1u64 << 48 {
            return fail("compressed offset exceeds48bits");
        }
        Ok(Self((compressed << 16) | u64::from(uncompressed)))
    }
    pub fn raw(self) -> u64 {
        self.0
    }
    pub fn packed(self) -> u64 {
        self.0
    }
    pub fn compressed(self) -> u64 {
        self.0 >> 16
    }
    pub fn uncompressed(self) -> u16 {
        self.0 as u16
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Chunk {
    start: VirtualOffset,
    end: VirtualOffset,
}
impl Chunk {
    pub fn start(self) -> VirtualOffset {
        self.start
    }
    pub fn end(self) -> VirtualOffset {
        self.end
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Metadata {
    pub reference_start: VirtualOffset,
    pub reference_end: VirtualOffset,
    pub mapped: u64,
    pub placed_unmapped: u64,
    pub csi_loffset: Option<VirtualOffset>,
}
#[derive(Debug)]
pub struct Bin {
    id: u32,
    loffset: Option<VirtualOffset>,
    chunks: Vec<Chunk>,
}
impl Bin {
    pub fn id(&self) -> u32 {
        self.id
    }
    pub fn loffset(&self) -> Option<VirtualOffset> {
        self.loffset
    }
    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }
}
#[derive(Debug)]
pub struct Reference {
    bins: Vec<Bin>,
    linear: Vec<VirtualOffset>,
    metadata: Option<Metadata>,
}
impl Reference {
    pub fn bins(&self) -> &[Bin] {
        &self.bins
    }
    pub fn linear_offsets(&self) -> &[VirtualOffset] {
        &self.linear
    }
    pub fn metadata(&self) -> Option<Metadata> {
        self.metadata
    }
}
#[derive(Clone, Debug)]
pub struct Limits {
    pub source_bytes: usize,
    pub references: usize,
    pub bins: usize,
    pub chunks: usize,
    pub linear_entries: usize,
    pub aux_bytes: usize,
    pub retained_bytes: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            source_bytes: 32 << 20,
            references: 4096,
            bins: 250_000,
            chunks: 1_000_000,
            linear_entries: 1_000_000,
            aux_bytes: 1 << 20,
            retained_bytes: 96 << 20,
        }
    }
}
impl Limits {
    fn validate(&self) -> Result<()> {
        if self.source_bytes > 64 << 20
            || self.references > 16384
            || self.bins > 1_000_000
            || self.chunks > 2_000_000
            || self.linear_entries > 2_000_000
            || self.aux_bytes > 4 << 20
            || self.retained_bytes > 256 << 20
        {
            return fail("limits exceed hard ceilings");
        }
        Ok(())
    }
}
#[derive(Clone, Debug)]
pub struct QueryLimits {
    pub chunks: usize,
    pub work: usize,
}
impl Default for QueryLimits {
    fn default() -> Self {
        Self {
            chunks: 100_000,
            work: 4_000_000,
        }
    }
}
/// Sealed complete index. The original file and auxiliary bytes remain exact;
/// no Deserialize or mutable access can manufacture validated parser state.
#[derive(Debug)]
pub struct Index {
    format: Format,
    min_shift: u8,
    depth: u8,
    original: Vec<u8>,
    aux: Vec<u8>,
    references: Vec<Reference>,
    unplaced: Option<u64>,
    charge: usize,
}
struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}
impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let e = add(self.at, n, self.bytes.len())?;
        let b = &self.bytes[self.at..e];
        self.at = e;
        Ok(b)
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn count(&mut self, max: usize) -> Result<usize> {
        let n = self.u32()?;
        if n > i32::MAX as u32 {
            return fail("negative or excessive index count");
        }
        cap(n as usize, max)
    }
    fn remaining(&self) -> usize {
        self.bytes.len() - self.at
    }
}
#[derive(Default)]
struct Counts {
    references: usize,
    bins: usize,
    chunks: usize,
    linear: usize,
    aux: usize,
}
struct Header {
    min_shift: u8,
    depth: u8,
    aux: usize,
    references: usize,
}
fn scheme(min: u32, depth: u32) -> Result<(u8, u8)> {
    if depth > 10 || min > 63 || min.checked_add(depth * 3).is_none_or(|n| n > 63) {
        return fail("unsupported CSI shift/depth or coordinate range");
    }
    Ok((min as u8, depth as u8))
}
fn bin_limit(depth: u8) -> u64 {
    ((1u64 << (3 * (u32::from(depth) + 1))) - 1) / 7
}
fn header(c: &mut Cursor<'_>, format: Format, l: &Limits) -> Result<Header> {
    let magic = c.take(4)?;
    if magic
        != match format {
            Format::Bai => b"BAI\x01",
            Format::Csi => b"CSI\x01",
        }
    {
        return fail("explicit index format/magic mismatch");
    }
    let (min_shift, depth, aux) = if format == Format::Csi {
        let min = c.u32()?;
        let depth = c.u32()?;
        let (min, depth) = scheme(min, depth)?;
        let n = c.count(l.aux_bytes)?;
        c.take(n)?;
        (min, depth, n)
    } else {
        (14, 5, 0)
    };
    let references = c.count(l.references)?;
    Ok(Header {
        min_shift,
        depth,
        aux,
        references,
    })
}
fn preflight(bytes: &[u8], format: Format, l: &Limits) -> Result<(Header, Counts, usize)> {
    l.validate()?;
    cap(bytes.len(), l.source_bytes)?;
    let mut c = Cursor { bytes, at: 0 };
    let h = header(&mut c, format, l)?;
    let mut n = Counts {
        references: h.references,
        aux: h.aux,
        ..Default::default()
    };
    let limit = bin_limit(h.depth);
    let metadata = limit + 1;
    for _ in 0..h.references {
        let bins = c.count(l.bins)?;
        n.bins = add(n.bins, bins, l.bins)?;
        for _ in 0..bins {
            let id = u64::from(c.u32()?);
            if id >= limit && id != metadata {
                return fail("reserved/out-of-range index bin");
            }
            if format == Format::Csi {
                c.u64()?;
            }
            let chunks = c.count(l.chunks)?;
            if id == metadata && chunks != 2 {
                return fail("metadata pseudo-bin requires exactly two entries");
            }
            n.chunks = add(n.chunks, chunks, l.chunks)?;
            for j in 0..chunks {
                let start = c.u64()?;
                let end = c.u64()?;
                if (id != metadata || j == 0) && start > end {
                    return fail("reversed chunk/reference offsets");
                }
            }
        }
        if format == Format::Bai {
            let intervals = c.count(l.linear_entries.min(1 << 17))?;
            n.linear = add(n.linear, intervals, l.linear_entries)?;
            c.take(
                intervals
                    .checked_mul(8)
                    .ok_or(Error("linear size overflow"))?,
            )?;
        }
    }
    if !matches!(c.remaining(), 0 | 8) {
        return fail("truncated or trailing index data");
    }
    if c.remaining() == 8 {
        c.u64()?;
    }
    // Include Vec headers, original bytes, auxiliary clone, bin sorting scratch,
    // chunks and linear arrays. Checked before any retained vectors/copies.
    let mut charge = add(bytes.len(), n.aux, l.retained_bytes)?;
    for (count, width) in [
        (n.references, 128usize),
        (n.bins, 128),
        (n.chunks, 32),
        (n.linear, 16),
    ] {
        charge = add(
            charge,
            count
                .checked_mul(width)
                .ok_or(Error("retention size overflow"))?,
            l.retained_bytes,
        )?;
    }
    Ok((h, n, charge))
}
impl Index {
    pub fn parse(bytes: &[u8], format: Format, limits: &Limits) -> Result<Self> {
        let (h, _, charge) = preflight(bytes, format, limits)?;
        let mut c = Cursor { bytes, at: 0 };
        header(&mut c, format, limits)?;
        let aux = if format == Format::Csi {
            bytes[16..16 + h.aux].to_vec()
        } else {
            Vec::new()
        };
        let mut references = Vec::with_capacity(h.references);
        let metadata_id = bin_limit(h.depth) + 1;
        for _ in 0..h.references {
            let count = c.count(limits.bins)?;
            let mut bins = Vec::with_capacity(count);
            let mut metadata = None;
            for _ in 0..count {
                let id = c.u32()?;
                let loffset = if format == Format::Csi {
                    Some(VirtualOffset(c.u64()?))
                } else {
                    None
                };
                let count = c.count(limits.chunks)?;
                if u64::from(id) == metadata_id {
                    if metadata.is_some() {
                        return fail("duplicate metadata pseudo-bin");
                    }
                    metadata = Some(Metadata {
                        reference_start: VirtualOffset(c.u64()?),
                        reference_end: VirtualOffset(c.u64()?),
                        mapped: c.u64()?,
                        placed_unmapped: c.u64()?,
                        csi_loffset: loffset,
                    });
                } else {
                    let mut chunks = Vec::with_capacity(count);
                    for _ in 0..count {
                        chunks.push(Chunk {
                            start: VirtualOffset(c.u64()?),
                            end: VirtualOffset(c.u64()?),
                        });
                    }
                    bins.push(Bin {
                        id,
                        loffset,
                        chunks,
                    });
                }
            }
            bins.sort_unstable_by_key(|b| b.id);
            if bins.windows(2).any(|w| w[0].id == w[1].id) {
                return fail("duplicate real index bin");
            }
            let mut linear = Vec::new();
            if format == Format::Bai {
                let n = c.count(limits.linear_entries)?;
                linear.reserve_exact(n);
                for _ in 0..n {
                    linear.push(VirtualOffset(c.u64()?));
                }
            }
            references.push(Reference {
                bins,
                linear,
                metadata,
            });
        }
        let unplaced = if c.remaining() == 8 {
            Some(c.u64()?)
        } else {
            None
        };
        Ok(Self {
            format,
            min_shift: h.min_shift,
            depth: h.depth,
            original: bytes.to_vec(),
            aux,
            references,
            unplaced,
            charge,
        })
    }
    pub fn format(&self) -> Format {
        self.format
    }
    pub fn min_shift(&self) -> u8 {
        self.min_shift
    }
    pub fn depth(&self) -> u8 {
        self.depth
    }
    pub fn coordinate_limit(&self) -> u64 {
        1u64 << (u32::from(self.min_shift) + 3 * u32::from(self.depth))
    }
    pub fn original_bytes(&self) -> &[u8] {
        &self.original
    }
    pub fn auxiliary_bytes(&self) -> &[u8] {
        &self.aux
    }
    pub fn references(&self) -> &[Reference] {
        &self.references
    }
    pub fn reference_count(&self) -> usize {
        self.references.len()
    }
    pub fn unplaced_count(&self) -> Option<u64> {
        self.unplaced
    }
    pub fn retained_charge(&self) -> usize {
        self.charge
    }
    /// Candidate chunks for zero-based [start,end); alignment records MUST be read
    /// and independently filtered. Metadata/linear/loffset never imply coverage.
    /// Deliberately conservative: no offset pruning based on imported summaries.
    pub fn query(
        &self,
        reference: usize,
        start: u64,
        end: u64,
        limits: &QueryLimits,
    ) -> Result<Vec<Chunk>> {
        if start >= end || end > self.coordinate_limit() {
            return fail("query requires nonempty in-range half-open interval");
        }
        if limits.chunks > 2_000_000 || limits.work > 64_000_000 {
            return fail("query limits exceed hard ceilings");
        }
        let r = self
            .references
            .get(reference)
            .ok_or(Error("unknown ordinal reference"))?;
        let mut count = 0;
        let mut work = 0;
        let mut ranges = [(0usize, 0usize); 11];
        for level in 0..=self.depth {
            let shift = u32::from(self.min_shift) + 3 * u32::from(self.depth - level);
            let base = ((1u64 << (3 * u32::from(level))) - 1) / 7;
            let first = base + (start >> shift);
            let last = base + ((end - 1) >> shift);
            let lo = r.bins.partition_point(|b| u64::from(b.id) < first);
            let hi = r.bins.partition_point(|b| u64::from(b.id) <= last);
            ranges[usize::from(level)] = (lo, hi);
            let search_levels = usize::BITS - r.bins.len().leading_zeros();
            work = add(work, hi - lo + 1 + 2 * search_levels as usize, limits.work)?;
            for b in &r.bins[lo..hi] {
                count = add(count, b.chunks.len(), limits.chunks)?;
            }
        }
        let log = if count < 2 {
            0
        } else {
            usize::BITS - (count - 1).leading_zeros()
        };
        add(
            work,
            count
                .checked_mul(log as usize + 3)
                .ok_or(Error("query work overflow"))?,
            limits.work,
        )?;
        let mut chunks = Vec::with_capacity(count);
        for (lo, hi) in ranges.into_iter().take(usize::from(self.depth) + 1) {
            chunks.extend(
                r.bins[lo..hi]
                    .iter()
                    .flat_map(|b| b.chunks.iter().copied())
                    .filter(|c| c.start < c.end),
            );
        }
        chunks.sort_unstable_by_key(|c| (c.start, c.end));
        let mut out: Vec<Chunk> = Vec::with_capacity(chunks.len());
        for c in chunks {
            if let Some(last) = out.last_mut() {
                if c.start <= last.end {
                    last.end = last.end.max(c.end);
                    continue;
                }
            }
            out.push(c);
        }
        Ok(out)
    }
    /// Caller-declared association, not an authenticated relationship established
    /// by the BAI/CSI format. References are ordinal, never inferred by name.
    pub fn bind<'a>(&'a self, association: &Association) -> Result<BoundIndex<'a>> {
        if association.source_blake3.len() != 64
            || !association
                .source_blake3
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            || association.source_bytes >= 1u64 << 48
            || association.reference_names.len() != self.references.len()
        {
            return fail("invalid declared source association");
        }
        let mut names = BTreeSet::new();
        let mut size = 0;
        for name in &association.reference_names {
            if name.trim().is_empty()
                || name.len() > 1024
                || name.chars().any(char::is_control)
                || !names.insert(name.as_str())
            {
                return fail("reference names must be bounded unique nonblank strings");
            }
            size = add(size, name.len(), 4 << 20)?;
        }
        let valid = |o: VirtualOffset| {
            o.compressed() < association.source_bytes
                || (o.compressed() == association.source_bytes && o.uncompressed() == 0)
        };
        for r in &self.references {
            for b in &r.bins {
                if b.loffset.is_some_and(|o| !valid(o)) {
                    return fail("CSI offset exceeds declared source bytes");
                }
                for c in &b.chunks {
                    if !valid(c.start) || !valid(c.end) {
                        return fail("chunk offset exceeds declared source bytes");
                    }
                }
            }
            if r.linear.iter().any(|o| !valid(*o))
                || r.metadata.is_some_and(|m| {
                    !valid(m.reference_start)
                        || !valid(m.reference_end)
                        || m.csi_loffset.is_some_and(|o| !valid(o))
                })
            {
                return fail("reference offset exceeds declared source bytes");
            }
        }
        Ok(BoundIndex {
            index: self,
            association: association.clone(),
        })
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Association {
    pub source_blake3: String,
    pub source_bytes: u64,
    pub reference_names: Vec<String>,
}
pub struct BoundIndex<'a> {
    index: &'a Index,
    association: Association,
}
impl BoundIndex<'_> {
    pub fn association(&self) -> &Association {
        &self.association
    }
    pub fn query(
        &self,
        expected_source_blake3: &str,
        reference: usize,
        start: u64,
        end: u64,
        limits: &QueryLimits,
    ) -> Result<Vec<Chunk>> {
        if expected_source_blake3 != self.association.source_blake3 {
            return fail("declared source revision changed");
        }
        self.index.query(reference, start, end, limits)
    }
}
