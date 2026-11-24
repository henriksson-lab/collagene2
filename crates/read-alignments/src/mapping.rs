use crate::*;
/// One declared CIGAR operation; reference coordinates remain unwrapped, zero-based half-open.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MapRun {
    pub kind: CigarKind,
    pub query_start: u64,
    pub query_end: u64,
    pub reference_start: u64,
    pub reference_end: u64,
    /// Implied original full read, including hard clips. This is not evidence that omitted calls exist.
    pub original_query_span: Option<(u64, u64)>,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoordinateMap {
    pub(crate) reference: Option<String>,
    pub(crate) reference_length: Option<u32>,
    pub(crate) circular: bool,
    pub(crate) dictionary_bound: bool,
    pub(crate) reverse: bool,
    pub(crate) reported_query_length: u64,
    pub(crate) implied_original_length: Option<u64>,
    pub(crate) sequence_present: bool,
    pub(crate) usable_reported_mapping: bool,
    pub(crate) diagnostics: Vec<String>,
    pub(crate) runs: Vec<MapRun>,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryProjection {
    Aligned {
        reference_position: u64,
        unwrapped_position: u64,
        assertion: CigarKind,
        implied_original_query: Option<u64>,
    },
    Inserted {
        reference_boundary: u64,
        implied_original_query: Option<u64>,
    },
    SoftClipped {
        implied_original_query: Option<u64>,
    },
    Unavailable,
    OutsideQuery,
}
impl CoordinateMap {
    pub fn reference(&self) -> Option<&str> {
        self.reference.as_deref()
    }
    pub fn reference_length(&self) -> Option<u32> {
        self.reference_length
    }
    pub fn circular(&self) -> bool {
        self.circular
    }
    pub fn dictionary_bound(&self) -> bool {
        self.dictionary_bound
    }
    pub fn reverse(&self) -> bool {
        self.reverse
    }
    pub fn reported_query_length(&self) -> u64 {
        self.reported_query_length
    }
    pub fn implied_original_length(&self) -> Option<u64> {
        self.implied_original_length
    }
    pub fn sequence_present(&self) -> bool {
        self.sequence_present
    }
    pub fn usable_reported_mapping(&self) -> bool {
        self.usable_reported_mapping
    }
    pub fn diagnostics(&self) -> &[String] {
        &self.diagnostics
    }
    pub fn runs(&self) -> &[MapRun] {
        &self.runs
    }
    pub fn query_position(&self, q: u64) -> QueryProjection {
        if q >= self.reported_query_length {
            return QueryProjection::OutsideQuery;
        }
        if !self.usable_reported_mapping {
            return QueryProjection::Unavailable;
        }
        let i = self.runs.partition_point(|r| r.query_end <= q);
        let Some(r) = self
            .runs
            .get(i)
            .filter(|r| r.query_start <= q && q < r.query_end)
        else {
            return QueryProjection::Unavailable;
        };
        let original = r.original_query_span.map(|(start, end)| {
            if self.reverse {
                end - 1 - (q - r.query_start)
            } else {
                start + (q - r.query_start)
            }
        });
        match r.kind {
            CigarKind::Match | CigarKind::Equal | CigarKind::Difference => {
                let p = r.reference_start + q - r.query_start;
                let reference_position = if self.circular {
                    p % u64::from(self.reference_length.unwrap())
                } else {
                    p
                };
                QueryProjection::Aligned {
                    reference_position,
                    unwrapped_position: p,
                    assertion: r.kind,
                    implied_original_query: original,
                }
            }
            CigarKind::Insertion => QueryProjection::Inserted {
                reference_boundary: r.reference_start,
                implied_original_query: original,
            },
            CigarKind::SoftClip => QueryProjection::SoftClipped {
                implied_original_query: original,
            },
            _ => QueryProjection::Unavailable,
        }
    }
}
pub(crate) fn project(record: &Record, h: &Header) -> Result<CoordinateMap> {
    if record.dictionary != h.dictionary_identity {
        return Err(Error::Invalid("mapping dictionary identity mismatch"));
    }
    let f = &record.fields;
    let reference = f.reference.as_deref().and_then(|name| h.reference(name));
    let circular = reference.is_some_and(|r| r.circular);
    let reverse = f.flags & 16 != 0;
    let mut diagnostics=vec!["CIGAR/flags are reported alignment assertions, not verified reference matches, variant calls or unique coverage. M does not distinguish matches and mismatches.".into()];
    let query = f
        .cigar
        .iter()
        .filter(|o| o.kind.consumes_query())
        .try_fold(0u64, |n, o| n.checked_add(u64::from(o.len)))
        .ok_or(Error::Invalid("query span overflow"))?;
    let left = f
        .cigar
        .first()
        .filter(|o| o.kind == CigarKind::HardClip)
        .map(|o| u64::from(o.len))
        .unwrap_or(0);
    let right = if f.cigar.len() > 1 {
        f.cigar
            .last()
            .filter(|o| o.kind == CigarKind::HardClip)
            .map(|o| u64::from(o.len))
            .unwrap_or(0)
    } else {
        0
    };
    let implied = (!f.cigar.is_empty())
        .then(|| query.checked_add(left).and_then(|n| n.checked_add(right)))
        .flatten();
    if !f.cigar.is_empty() && implied.is_none() {
        return Err(Error::Invalid("original CIGAR span overflow"));
    }
    let reported_query_length = f.sequence.as_ref().map(|s| s.len() as u64).unwrap_or(query);
    if f.sequence.is_none() {
        diagnostics.push("SEQ is absent. Query extents come only from the reported CIGAR; missing sequence/quality calls were not synthesized.".into());
    }
    if left + right > 0 {
        diagnostics.push("Hard-clipped calls are absent from SEQ; original coordinates include only the implied omitted extents. Reverse orientation swaps the original clip ends.".into());
    }
    if f.flags & 256 != 0 {
        diagnostics.push("Secondary alignment retained; no unique-placement inference.".into());
    }
    if f.flags & 2048 != 0 {
        diagnostics.push("Supplementary alignment retained; no independent-read or structural-variant inference.".into());
    }
    let mut usable =
        f.flags & 4 == 0 && f.reference.is_some() && f.start.is_some() && !f.cigar.is_empty();
    if f.flags & 4 != 0 {
        diagnostics.push("Unmapped flag prevents use of RNAME/POS/CIGAR as an alignment, even when those fields are populated.".into());
    }
    if reference.is_none() {
        diagnostics.push("Reference dictionary identity/length is unavailable; coordinates remain unbound reported values.".into());
    }
    let mut runs = Vec::with_capacity(f.cigar.len());
    let mut q = 0;
    let mut r = u64::from(f.start.unwrap_or(0));
    for op in &f.cigar {
        let qe = q + if op.kind.consumes_query() {
            u64::from(op.len)
        } else {
            0
        };
        let re = r
            .checked_add(if op.kind.consumes_reference() {
                u64::from(op.len)
            } else {
                0
            })
            .ok_or(Error::Invalid("reference span overflow"))?;
        let original_query_span = if op.kind.consumes_query() {
            implied.map(|n| {
                if reverse {
                    (n - left - qe, n - left - q)
                } else {
                    (left + q, left + qe)
                }
            })
        } else {
            None
        };
        runs.push(MapRun {
            kind: op.kind,
            query_start: q,
            query_end: qe,
            reference_start: r,
            reference_end: re,
            original_query_span,
        });
        q = qe;
        r = re;
    }
    if let Some(reference) = reference
        && (f.start.is_some_and(|p| p >= reference.length)
            || !circular && r > u64::from(reference.length))
    {
        usable = false;
        diagnostics.push(
            "Reported alignment exceeds reference dictionary bounds; projections are withheld."
                .into(),
        );
    }

    Ok(CoordinateMap {
        reference: f.reference.clone(),
        reference_length: reference.map(|r| r.length),
        circular,
        dictionary_bound: reference.is_some(),
        reverse,
        reported_query_length,
        implied_original_length: implied,
        sequence_present: f.sequence.is_some(),
        usable_reported_mapping: usable,
        diagnostics,
        runs,
    })
}
