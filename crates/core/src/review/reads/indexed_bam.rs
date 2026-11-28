//! Source-bound BAI/CSI selection over a privately completed BAM collection.
//!
//! The index is an untrusted acceleration hint. This initial retained-collection
//! adapter independently scans reported coordinates to detect every omitted
//! mapped record in the requested interval. It does not establish actual-base
//! agreement, unique mapping, coverage or biological absence.
use crate::review::reads::store::Collection;
use collagene_genomic_index::{Format, Index, Limits, QueryLimits};
use collagene_model::Error;
use std::{
    borrow::Cow,
    io::{Cursor, Read},
};

/// Private live parser state; a serialized query result cannot recreate it.
pub struct IndexedBam {
    source_blake3: String,
    index_blake3: String,
    index: Index,
}
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Selection {
    pub source_blake3: String,
    pub index_blake3: String,
    pub reference_name: String,
    /// Zero-based, half-open reference coordinates.
    pub start: u64,
    pub end: u64,
    pub candidate_chunks: usize,
    pub records_audited: usize,
    pub record_indices: Vec<usize>,
    pub interpretation: String,
}
fn bad(message: impl std::fmt::Display) -> Error {
    Error::Format(format!("indexed BAM: {message}"))
}
impl IndexedBam {
    /// Bind the parsed index to actual immutable BAM bytes and header order.
    /// Equal dictionary counts do not authenticate an index's provenance.
    pub fn bind(
        collection: &Collection,
        bytes: &[u8],
        format: Format,
        limits: &Limits,
    ) -> Result<Self, Error> {
        if collection.format() != crate::review::reads::store::Format::Bam {
            return Err(bad("a privately parsed BAM collection is required"));
        }
        let decoded = index_transport(bytes, format, limits)?;
        let index = Index::parse(&decoded, format, limits).map_err(bad)?;
        if index.reference_count() != collection.header().references().len() {
            return Err(bad(
                "index reference count differs from actual BAM dictionary",
            ));
        }
        Ok(Self {
            source_blake3: collection.source_blake3().into(),
            index_blake3: blake3::hash(bytes).to_hex().to_string(),
            index,
        })
    }
    pub fn source_blake3(&self) -> &str {
        &self.source_blake3
    }
    pub fn index_blake3(&self) -> &str {
        &self.index_blake3
    }

    /// Return only mapped records whose *reported* alignment span intersects
    /// the explicit interval. CIGAR deletions/skips consume reference; insertions
    /// and clips do not. A zero-reference-consuming mapped record occupies one
    /// reference position, matching BAM binning conventions. Placed unmapped
    /// records are deliberately excluded. Full retained scanning is bounded;
    /// resource exhaustion and an index omitting any matching record are errors,
    /// never a partial result or proof of absence.
    pub fn select(
        &self,
        collection: &Collection,
        reference_name: &str,
        start: u64,
        end: u64,
        query_limits: &QueryLimits,
        max_audit_records: usize,
        max_matches: usize,
    ) -> Result<Selection, Error> {
        if collection.source_blake3() != self.source_blake3 {
            return Err(bad("source bytes changed; rebind the index"));
        }
        if max_audit_records == 0
            || max_audit_records > 1_000_000
            || max_matches == 0
            || max_matches > 100_000
            || collection.len() > max_audit_records
        {
            return Err(Error::TooLarge);
        }
        let (reference_index, reference) = collection
            .header()
            .references()
            .iter()
            .enumerate()
            .find(|(_, reference)| reference.name == reference_name)
            .ok_or_else(|| bad("explicit reference name is absent from actual BAM dictionary"))?;
        if start >= end || end > u64::from(reference.length) {
            return Err(bad(
                "interval must be nonempty and inside the actual BAM reference",
            ));
        }
        let query = self
            .index
            .query(reference_index, start, end, query_limits)
            .map_err(bad)?;
        let source_bytes = collection.original_bytes().len() as u64;
        for chunk in &query {
            for raw in [chunk.start().raw(), chunk.end().raw()] {
                let compressed = raw >> 16;
                if compressed > source_bytes || (compressed == source_bytes && raw & 65535 != 0) {
                    return Err(bad("index virtual offset exceeds the actual source bytes"));
                }
            }
        }
        let mut indices = Vec::new();
        for position in 0..collection.len() {
            let record = collection
                .record(position)
                .ok_or_else(|| bad("missing retained record"))?;
            let fields = record.fields();
            if fields.flags & 4 != 0 || fields.reference.as_deref() != Some(reference_name) {
                continue;
            }
            let a = u64::from(
                fields
                    .start
                    .ok_or_else(|| bad("mapped record has no reported start"))?,
            );
            let length = fields
                .cigar
                .iter()
                .filter(|op| op.kind.consumes_reference())
                .try_fold(0u64, |sum, op| sum.checked_add(u64::from(op.len)))
                .ok_or_else(|| bad("reported CIGAR reference span overflow"))?;
            let b = a
                .checked_add(length.max(1))
                .ok_or_else(|| bad("reported alignment span overflow"))?;
            if a >= end || b <= start {
                continue;
            }
            let raw = record
                .virtual_offset()
                .ok_or_else(|| bad("BAM record lacks a physical virtual offset"))?
                .packed()
                .map_err(bad)?;
            let candidate = query.partition_point(|chunk| chunk.end().raw() <= raw);
            if !query
                .get(candidate)
                .is_some_and(|chunk| chunk.start().raw() <= raw && raw < chunk.end().raw())
            {
                return Err(bad(format!(
                    "index omits actual overlapping mapped record {position}; no selection published"
                )));
            }
            if indices.len() == max_matches {
                return Err(Error::TooLarge);
            }
            indices.push(position);
        }
        Ok(Selection {
            source_blake3: self.source_blake3.clone(), index_blake3: self.index_blake3.clone(),
            reference_name: reference_name.into(), start, end,
            candidate_chunks: query.len(), records_audited: collection.len(), record_indices: indices,
            interpretation: "Reported mapped-record spans only; all retained records were independently scanned to detect index omissions in this interval. Index/header association is not authenticated. This is not actual-reference checking, coverage, variant evidence, unique mapping or biological absence. Unmapped records are excluded; skipped/deleted reference spans can intersect the query.".into(),
        })
    }
}

/// CSI is commonly transported in BGZF. Decode strict blocks, CRC, terminal EOF
/// and exact input consumption under separate source/decoded ceilings before
/// handing bytes to the format-only index parser. BAI remains a raw format here.
fn index_transport<'a>(
    bytes: &'a [u8],
    format: Format,
    limits: &Limits,
) -> Result<Cow<'a, [u8]>, Error> {
    if bytes.len() > limits.source_bytes || limits.source_bytes > 64 << 20 {
        return Err(Error::TooLarge);
    }
    if format != Format::Csi || !bytes.starts_with(&[31, 139]) {
        return Ok(Cow::Borrowed(bytes));
    }
    let policy = collagene_read_alignments::Limits {
        max_input_bytes: limits.source_bytes as u64,
        max_uncompressed_bytes: limits.source_bytes as u64,
        max_bgzf_blocks: 65_536,
        ..Default::default()
    };
    let mut reader =
        collagene_read_alignments::BgzfReader::new(Cursor::new(bytes), policy).map_err(bad)?;
    let mut decoded = Vec::new();
    (&mut reader)
        .take(limits.source_bytes as u64 + 1)
        .read_to_end(&mut decoded)
        .map_err(bad)?;
    if decoded.len() > limits.source_bytes {
        return Err(Error::TooLarge);
    }
    if reader.source_blake3().is_none() || reader.source_bytes() != bytes.len() as u64 {
        return Err(bad("CSI BGZF transport did not reach an exact checked EOF"));
    }
    Ok(Cow::Owned(decoded))
}

/// Build a fresh BAI from actual parsed record offsets and reported spans.
/// Actual coordinate ordering is checked; an SO:coordinate header alone is
/// insufficient. Placed unmapped records are indexed too; unplaced records
/// must follow all placed records. No linear/pruning or count pseudo-bins are
/// emitted. This is an index of source assertions, not an alignment verifier.
pub fn build_bai(collection: &Collection, max_output_bytes: usize) -> Result<Vec<u8>, Error> {
    use std::collections::BTreeMap;
    if collection.format() != crate::review::reads::store::Format::Bam {
        return Err(bad("BAI construction requires a completed BAM collection"));
    }
    if max_output_bytes == 0 || max_output_bytes > 32 << 20 || collection.len() > 1_000_000 {
        return Err(Error::TooLarge);
    }
    let references = collection.header().references();
    if references.len() > 4096 {
        return Err(Error::TooLarge);
    }
    let names: BTreeMap<&str, usize> = references
        .iter()
        .enumerate()
        .map(|(i, r)| (r.name.as_str(), i))
        .collect();
    let mut bins: Vec<BTreeMap<u32, Vec<(u64, u64)>>> =
        (0..references.len()).map(|_| BTreeMap::new()).collect();
    let mut encoded = 8usize
        .checked_add(references.len() * 8)
        .ok_or(Error::TooLarge)?;
    if encoded > max_output_bytes {
        return Err(Error::TooLarge);
    }
    let mut charge = references.len() * 128;
    let mut previous = None;
    let mut unplaced_started = false;
    for position in 0..collection.len() {
        let record = collection
            .record(position)
            .ok_or_else(|| bad("missing retained record"))?;
        let fields = record.fields();
        let (Some(name), Some(start)) = (fields.reference.as_deref(), fields.start) else {
            unplaced_started = true;
            continue;
        };
        if unplaced_started {
            return Err(bad(
                "placed record follows unplaced records; coordinate sorting is required",
            ));
        }
        let reference = *names
            .get(name)
            .ok_or_else(|| bad("record reference is absent from actual BAM dictionary"))?;
        if previous.is_some_and(|last| last > (reference, start)) {
            return Err(bad("actual BAM records are not coordinate sorted"));
        }
        previous = Some((reference, start));
        let length = fields
            .cigar
            .iter()
            .filter(|op| op.kind.consumes_reference())
            .try_fold(0u64, |sum, op| sum.checked_add(u64::from(op.len)))
            .ok_or_else(|| bad("reported CIGAR span overflow"))?;
        let a = u64::from(start);
        let b = a
            .checked_add(length.max(1))
            .ok_or_else(|| bad("reported alignment span overflow"))?;
        if b > 1 << 29 || b > u64::from(references[reference].length) {
            return Err(bad(
                "reported span cannot be represented by BAI; review reference bounds or use CSI",
            ));
        }
        let last = b - 1;
        let bin = [(14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)]
            .into_iter()
            .find_map(|(shift, base)| {
                ((a >> shift) == (last >> shift)).then_some(base + (a >> shift) as u32)
            })
            .unwrap_or(0);
        let begin = record
            .virtual_offset()
            .ok_or_else(|| bad("missing actual BAM offset"))?
            .packed()
            .map_err(bad)?;
        let end = if let Some(next) = collection.record(position + 1) {
            next.virtual_offset()
                .ok_or_else(|| bad("missing next actual BAM offset"))?
                .packed()
                .map_err(bad)?
        } else {
            (collection.original_bytes().len() as u64) << 16
        };
        if begin >= end {
            return Err(bad("actual record offsets are not strictly ordered"));
        }
        let map = &mut bins[reference];
        if !map.contains_key(&bin) {
            encoded = encoded
                .checked_add(8)
                .filter(|n| *n <= max_output_bytes)
                .ok_or(Error::TooLarge)?;
            charge = charge
                .checked_add(128)
                .filter(|n| *n <= 128 << 20)
                .ok_or(Error::TooLarge)?;
        }
        let chunks = map.entry(bin).or_default();
        if let Some(chunk) = chunks.last_mut().filter(|chunk| chunk.1 == begin) {
            chunk.1 = end;
        } else {
            encoded = encoded
                .checked_add(16)
                .filter(|n| *n <= max_output_bytes)
                .ok_or(Error::TooLarge)?;
            charge = charge
                .checked_add(64)
                .filter(|n| *n <= 128 << 20)
                .ok_or(Error::TooLarge)?;
            chunks.push((begin, end));
        }
    }
    let mut output = Vec::with_capacity(encoded);
    output.extend_from_slice(b"BAI\x01");
    output.extend_from_slice(&(references.len() as i32).to_le_bytes());
    for bins in bins {
        output.extend_from_slice(&(bins.len() as i32).to_le_bytes());
        for (bin, chunks) in bins {
            output.extend_from_slice(&bin.to_le_bytes());
            output.extend_from_slice(&(chunks.len() as i32).to_le_bytes());
            for (a, b) in chunks {
                output.extend_from_slice(&a.to_le_bytes());
                output.extend_from_slice(&b.to_le_bytes());
            }
        }
        output.extend_from_slice(&0i32.to_le_bytes());
    }
    debug_assert_eq!(output.len(), encoded);
    Ok(output)
}
