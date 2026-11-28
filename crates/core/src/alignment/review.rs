//! Source-bound pairwise review using the existing native affine algorithms.
//! Stored receipts and coordinate maps are reference evidence, not sample proof.
use crate::alignment::{
    self, Alignment, AlignmentMode, AlignmentOptions, ComparisonOptions, MoleculeComparison,
};
use collagene_model::{Error, Molecule, Topology};
use serde::{Deserialize, Serialize};

pub const ENGINE_VERSION: &str = "pairwise-review-1";
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Options {
    Stored(AlignmentOptions),
    Whole(ComparisonOptions),
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine_version: String,
    pub context_revision: String,
    pub reference_revision: String,
    pub query_revision: String,
    pub reference_length: usize,
    pub query_length: usize,
    pub options: Options,
    #[cfg_attr(feature = "json-schema", schemars(with = "serde_json::Value"))]
    pub comparison: MoleculeComparison,
}
fn bad(message: &str) -> Error {
    Error::Format(format!("Pairwise review: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
impl Options {
    pub fn alignment(self) -> AlignmentOptions {
        match self {
            Self::Stored(o) => o,
            Self::Whole(o) => o.alignment,
        }
    }
    pub fn validate(self) -> Result<(), Error> {
        let options = self.alignment();
        let s = options.scoring;
        if !(1..=10_000).contains(&s.match_score)
            || !(-10_000..=0).contains(&s.mismatch_score)
            || s.ambiguous_score < s.mismatch_score
            || s.ambiguous_score > s.match_score
            || !(-10_000..=-1).contains(&s.gap_open)
            || !(-10_000..=-1).contains(&s.gap_extend)
            || options.max_cells == 0
            || options.max_cells > alignment::HARD_MAX_CELLS
        {
            return Err(bad("invalid native scoring or cell bounds"));
        }
        if let Self::Whole(o) = self {
            if o.alignment.mode != AlignmentMode::Global
                || o.max_origins == 0
                || o.max_origins > 64
                || o.max_total_cells == 0
                || o.max_total_cells > 2 * alignment::HARD_MAX_CELLS
            {
                return Err(bad(
                    "whole-molecule search requires global mode and native origin/total-cell bounds",
                ));
            }
        }
        Ok(())
    }
}
pub fn revision(
    reference: &Molecule,
    query: &Molecule,
    options: &Options,
) -> Result<String, Error> {
    collagene_model::bounds::model_preflight(&[reference, query])?;
    reference.validate()?;
    query.validate()?;
    options.validate()?;
    if reference.sequence.alphabet() == collagene_model::Alphabet::Protein
        || reference.sequence.alphabet() != query.sequence.alphabet()
    {
        return Err(bad("matching DNA or RNA alphabets are required"));
    }
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, reference, query, options))
}
pub fn compare(reference: &Molecule, query: &Molecule, options: &Options) -> Result<Report, Error> {
    let context_revision = revision(reference, query, options)?;
    let comparison=match options {
        Options::Stored(o)=>MoleculeComparison {
            alignment:alignment::align(&reference.sequence,&query.sequence,o)?,
            reference_revision:reference.revision(),query_revision:query.revision(),
            reference_origin:0,query_reversed:false,
            reference_topology:reference.topology,query_topology:query.topology,
            evaluated_candidates:1,exhaustive_reference_origins:true,exhaustive_query_orientations:true,optimality_proven:true,
            warnings:vec!["Stored origin/orientation policy: no circular rotation or reverse-complement search was requested. Local/query-in-reference flanks outside the alignment are not verified.".into()],
        },
        Options::Whole(o)=>alignment::compare_molecules(reference,query,o)?,
    };
    let report = Report {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        reference_revision: reference.revision(),
        query_revision: query.revision(),
        reference_length: reference.sequence.len(),
        query_length: query.sequence.len(),
        options: *options,
        comparison,
    };
    report.validate_sources(reference, query)?;
    Ok(report)
}
impl Report {
    /// Structural receipt guard; this does not rerun alignment or prove optimality.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.options.validate()?;
        let c = &self.comparison;
        let a = &c.alignment;
        if self.engine_version != ENGINE_VERSION
            || !digest(expected)
            || self.context_revision != expected
            || !digest(&self.reference_revision)
            || !digest(&self.query_revision)
            || self.reference_revision != c.reference_revision
            || self.query_revision != c.query_revision
            || self.reference_length == 0
            || self.query_length == 0
            || self.reference_length > collagene_model::MAX_SEQUENCE_LEN
            || self.query_length > collagene_model::MAX_SEQUENCE_LEN
            || c.reference_origin >= self.reference_length
            || (c.reference_topology != Topology::Circular && c.reference_origin != 0)
            || c.evaluated_candidates == 0
            || c.evaluated_candidates > 128
            || a.algorithm != "collagene-affine-pairwise-1"
            || a.options != self.options.alignment()
            || !digest(&a.reference_hash)
            || !digest(&a.query_hash)
            || a.reference_start > a.reference_end
            || a.reference_end > self.reference_length
            || a.query_start > a.query_end
            || a.query_end > self.query_length
            || a.reference_aligned.len() != a.query_aligned.len()
            || a.reference_aligned.len() > alignment::MAX_ALIGNMENT_COLUMNS
            || a.operations.len() > alignment::MAX_ALIGNMENT_COLUMNS
            || a.differences.len() > alignment::MAX_ALIGNMENT_COLUMNS
        {
            return Err(bad("invalid comparison receipt identity or bounds"));
        }
        match a.options.mode {
            AlignmentMode::Global
                if a.reference_start != 0
                    || a.reference_end != self.reference_length
                    || a.query_start != 0
                    || a.query_end != self.query_length =>
            {
                return Err(bad("global alignment must cover both complete inputs"));
            }
            AlignmentMode::QueryInReference
                if a.query_start != 0 || a.query_end != self.query_length =>
            {
                return Err(bad(
                    "query-in-reference alignment must cover the complete query",
                ));
            }
            AlignmentMode::Local
                if (a.reference_aligned.is_empty()
                    && (a.reference_start != 0
                        || a.reference_end != 0
                        || a.query_start != 0
                        || a.query_end != 0))
                    || (!a.reference_aligned.is_empty() && a.score <= 0) =>
            {
                return Err(bad(
                    "local alignment must have positive score or canonical empty coordinates",
                ));
            }
            _ => {}
        }
        if c.optimality_proven
            && !(c.exhaustive_reference_origins && c.exhaustive_query_orientations)
            && !(a.options.mode == AlignmentMode::Global
                && a.options.scoring.ambiguous_score >= 0
                && a.reference_aligned == a.query_aligned
                && !a.reference_aligned.contains('-'))
        {
            return Err(bad(
                "non-exhaustive optimality requires the native symbolic identity shortcut",
            ));
        }
        match self.options {
            Options::Stored(_) => {
                if c.reference_origin != 0
                    || c.query_reversed
                    || c.evaluated_candidates != 1
                    || !c.exhaustive_reference_origins
                    || !c.exhaustive_query_orientations
                    || !c.optimality_proven
                {
                    return Err(bad(
                        "stored policy does not perform orientation/origin search",
                    ));
                }
            }
            Options::Whole(o) => {
                if (!o.both_query_strands && c.query_reversed)
                    || (!o.reference_rotations && c.reference_origin != 0)
                {
                    return Err(bad(
                        "returned orientation/origin differs from requested policy",
                    ));
                }
            }
        }
        validate_alignment(a, self.reference_length, self.query_length)?;
        crate::encoding::bounded_json::hash(self)?;
        Ok(())
    }
    /// Bind letters, native coordinate systems and options to actual source models.
    /// Alignment optimality remains a fresh-computation claim, not archive authority.
    pub fn validate_sources(&self, reference: &Molecule, query: &Molecule) -> Result<(), Error> {
        let expected = revision(reference, query, &self.options)?;
        self.validate_binding(&expected)?;
        let c = &self.comparison;
        let a = &c.alignment;
        if self.reference_revision != reference.revision()
            || self.query_revision != query.revision()
            || self.reference_length != reference.sequence.len()
            || self.query_length != query.sequence.len()
            || c.reference_topology != reference.topology
            || c.query_topology != query.topology
        {
            return Err(bad("comparison source models changed"));
        }
        let rb = reference.sequence.as_str().as_bytes();
        let oriented = if c.query_reversed {
            query.sequence.reverse_complement()?
        } else {
            query.sequence.clone()
        };
        let qb = oriented.as_str().as_bytes();
        let mut rh = blake3::Hasher::new();
        rh.update(&rb[c.reference_origin..]);
        rh.update(&rb[..c.reference_origin]);
        if rh.finalize().to_hex().as_str() != a.reference_hash
            || blake3::hash(qb).to_hex().as_str() != a.query_hash
        {
            return Err(bad("oriented source hashes differ from the alignment"));
        }
        let (mut r, mut q) = (a.reference_start, a.query_start);
        for (ra, qa) in a.reference_aligned.bytes().zip(a.query_aligned.bytes()) {
            if ra != b'-' {
                if ra != rb[(r + c.reference_origin) % rb.len()] {
                    return Err(bad("reference columns differ from native source"));
                }
                r += 1;
            }
            if qa != b'-' {
                if qa != qb[q] {
                    return Err(bad("query columns differ from native source"));
                }
                q += 1;
            }
        }
        Ok(())
    }
    /// Original-reference half-open spans; an insertion is the explicit [(x,x)] caret.
    /// Origin-crossing selections retain their two ordered parts.
    pub fn map_reference(&self, start: usize, end: usize) -> Result<Vec<(usize, usize)>, Error> {
        if start > end
            || end > self.reference_length
            || self.reference_length == 0
            || self.comparison.reference_origin >= self.reference_length
        {
            return Err(bad("reference navigation lies outside its source"));
        }
        let origin = self.comparison.reference_origin;
        if origin == 0 {
            return Ok(vec![(start, end)]);
        }
        let n = self.reference_length;
        let s = (origin + start) % n;
        let length = end - start;
        if length == 0 {
            return Ok(vec![(s, s)]);
        }
        let first = length.min(n - s);
        let mut spans = vec![(s, s + first)];
        if first < length {
            spans.push((0, length - first));
        }
        Ok(spans)
    }
    /// Original-query half-open interval or deletion caret, preserving reverse orientation.
    pub fn map_query(&self, start: usize, end: usize) -> Result<(usize, usize), Error> {
        if start > end || end > self.query_length {
            return Err(bad("query navigation lies outside its source"));
        }
        Ok(if self.comparison.query_reversed {
            (self.query_length - end, self.query_length - start)
        } else {
            (start, end)
        })
    }
}

fn validate_alignment(
    a: &Alignment,
    reference_length: usize,
    query_length: usize,
) -> Result<(), Error> {
    use crate::alignment::{CigarKind as K, Difference, DifferenceKind as D};
    let (mut r, mut q) = (a.reference_start, a.query_start);
    let mut offset = 0usize;
    let mut differences = Vec::new();
    let mut cigar = String::new();
    let (mut certain, mut ambiguous, mut mismatch, mut inserted, mut deleted) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    let mut score = 0i64;
    let mut previous = None;
    for run in &a.operations {
        if previous == Some(run.kind) {
            return Err(bad("adjacent CIGAR runs require native canonical merging"));
        }
        previous = Some(run.kind);
        let end = offset
            .checked_add(run.length)
            .filter(|e| *e <= a.reference_aligned.len())
            .ok_or_else(|| bad("operation span exceeds alignment"))?;
        if run.length == 0 {
            return Err(bad("empty alignment operation"));
        }
        let mut rb = String::new();
        let mut qb = String::new();
        let (rs, qs) = (r, q);
        for (ra, qa) in a.reference_aligned.as_bytes()[offset..end]
            .iter()
            .copied()
            .zip(a.query_aligned.as_bytes()[offset..end].iter().copied())
        {
            let rm = crate::sequence::analysis::base_mask(ra);
            let qm = crate::sequence::analysis::base_mask(qa);
            let valid = match run.kind {
                K::Match => {
                    ra != b'-'
                        && qa != b'-'
                        && rm != 0
                        && qm != 0
                        && rm.count_ones() == 1
                        && qm.count_ones() == 1
                        && ra == qa
                }
                K::Ambiguous => {
                    ra != b'-'
                        && qa != b'-'
                        && rm != 0
                        && qm != 0
                        && rm & qm != 0
                        && (rm.count_ones() > 1 || qm.count_ones() > 1)
                }
                K::Mismatch => ra != b'-' && qa != b'-' && rm != 0 && qm != 0 && rm & qm == 0,
                K::Insertion => ra == b'-' && qa != b'-' && qm != 0,
                K::Deletion => qa == b'-' && ra != b'-' && rm != 0,
            };
            if !valid {
                return Err(bad(
                    "alignment columns disagree with their typed operations",
                ));
            }
            if ra != b'-' {
                rb.push(ra as char);
                r += 1;
            }
            if qa != b'-' {
                qb.push(qa as char);
                q += 1;
            }
        }
        let (symbol, kind, unit) = match run.kind {
            K::Match => {
                certain += run.length;
                ('=', None, a.options.scoring.match_score)
            }
            K::Ambiguous => {
                ambiguous += run.length;
                (
                    'M',
                    Some(D::AmbiguousCompatible),
                    a.options.scoring.ambiguous_score,
                )
            }
            K::Mismatch => {
                mismatch += run.length;
                ('X', Some(D::Substitution), a.options.scoring.mismatch_score)
            }
            K::Insertion => {
                inserted += run.length;
                ('I', Some(D::Insertion), 0)
            }
            K::Deletion => {
                deleted += run.length;
                ('D', Some(D::Deletion), 0)
            }
        };
        if matches!(run.kind, K::Insertion | K::Deletion) {
            score += i64::from(a.options.scoring.gap_open)
                + (run.length as i64 - 1) * i64::from(a.options.scoring.gap_extend);
        } else {
            score += run.length as i64 * i64::from(unit);
        }
        use std::fmt::Write;
        write!(cigar, "{}{symbol}", run.length).map_err(|_| bad("CIGAR formatting failed"))?;
        if let Some(kind) = kind {
            differences.push(Difference {
                kind,
                reference_start: rs,
                reference_end: r,
                query_start: qs,
                query_end: q,
                reference_bases: rb,
                query_bases: qb,
            });
        }
        offset = end;
    }
    let paired = certain + ambiguous + mismatch;
    let stats = &a.stats;
    let identity = if paired == 0 {
        None
    } else {
        Some(certain as f64 / paired as f64)
    };
    let compatible = if paired == 0 {
        None
    } else {
        Some((certain + ambiguous) as f64 / paired as f64)
    };
    if offset != a.reference_aligned.len()
        || r != a.reference_end
        || q != a.query_end
        || score != i64::from(a.score)
        || cigar != a.cigar
        || differences != a.differences
        || stats.columns != offset
        || stats.paired_bases != paired
        || stats.certain_matches != certain
        || stats.ambiguous_compatible != ambiguous
        || stats.mismatches != mismatch
        || stats.inserted_query_bases != inserted
        || stats.deleted_reference_bases != deleted
        || stats.reference_span_bases != r - a.reference_start
        || stats.query_span_bases != q - a.query_start
        || stats.reference_coverage != paired as f64 / reference_length as f64
        || stats.query_coverage != paired as f64 / query_length as f64
        || stats.certain_identity != identity
        || stats.compatible_identity != compatible
    {
        return Err(bad(
            "alignment coordinates, score, differences or statistics are inconsistent",
        ));
    }
    Ok(())
}
