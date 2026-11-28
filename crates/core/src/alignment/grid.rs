//! Checked storage and score certificates for the current affine recurrence.
//! Grid storage and certificates serve one current native affine recurrence.
use super::{AlignmentMode, Scoring};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::mem::size_of;

pub(crate) const NEG: i32 = -1_000_000_000;
pub(crate) const MAX_MATRIX_BYTES: usize = 96 * 1024 * 1024;
pub const MAX_ATTEMPTS: usize = 32;
fn bad(message: &str) -> Error {
    Error::Format(message.into())
}
macro_rules! data {
    ($item:item) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        $item
    };
}
data! {
    pub enum GlobalStorage {
        Full,
        Adaptive {
            initial_radius: usize,
            maximum_radius: usize,
            max_attempts: usize,
        },
    }
}
impl Default for GlobalStorage {
    fn default() -> Self {
        Self::Adaptive {
            initial_radius: 16,
            maximum_radius: 4096,
            max_attempts: 10,
        }
    }
}
impl GlobalStorage {
    pub fn validate(self, mode: AlignmentMode) -> Result<(), Error> {
        if let Self::Adaptive {
            initial_radius,
            maximum_radius,
            max_attempts,
        } = self
        {
            if mode != AlignmentMode::Global {
                return Err(bad("Adaptive row bands require global alignment"));
            }
            if initial_radius > maximum_radius
                || maximum_radius > super::MAX_ALIGNMENT_COLUMNS
                || !(1..=MAX_ATTEMPTS).contains(&max_attempts)
            {
                return Err(bad("Invalid adaptive alignment radii or attempt bound"));
            }
        }
        Ok(())
    }
}
data! {
    pub enum Domain {
        Full,
        /// Endpoint-aware diagonal corridor, not a symmetric band around zero.
        Corridor {
            radius: usize,
        },
    }
}
impl Domain {
    fn diagonals(self, n: usize, m: usize) -> Result<Option<(i64, i64)>, Error> {
        let Self::Corridor { radius } = self else {
            return Ok(None);
        };
        let n = i64::try_from(n).map_err(|_| Error::TooLarge)?;
        let m = i64::try_from(m).map_err(|_| Error::TooLarge)?;
        let radius = i64::try_from(radius).map_err(|_| Error::TooLarge)?;
        let delta = n.checked_sub(m).ok_or(Error::TooLarge)?;
        Ok(Some((
            delta.min(0).checked_sub(radius).ok_or(Error::TooLarge)?,
            delta.max(0).checked_add(radius).ok_or(Error::TooLarge)?,
        )))
    }
    pub(crate) fn row(self, i: usize, n: usize, m: usize) -> Result<(usize, usize), Error> {
        if i > n {
            return Err(bad("Row index exceeds alignment reference"));
        }
        let Some((lower, upper)) = self.diagonals(n, m)? else {
            return Ok((0, m));
        };
        let i = i64::try_from(i).map_err(|_| Error::TooLarge)?;
        let end = i64::try_from(m).map_err(|_| Error::TooLarge)?;
        let start = i.checked_sub(upper).ok_or(Error::TooLarge)?.max(0);
        let end = i.checked_sub(lower).ok_or(Error::TooLarge)?.min(end);
        if start > end {
            return Err(bad("Alignment corridor has a disconnected row"));
        }
        Ok((
            usize::try_from(start).map_err(|_| Error::TooLarge)?,
            usize::try_from(end).map_err(|_| Error::TooLarge)?,
        ))
    }
    /// Exact admitted-cell count without allocating a row table. Intermediate
    /// rectangle arithmetic uses u128 so narrow bands also work on wasm32.
    pub fn cell_count(self, n: usize, m: usize) -> Result<usize, Error> {
        let n128 = n as u128;
        let m128 = m as u128;
        let rectangle = (n128 + 1).checked_mul(m128 + 1).ok_or(Error::TooLarge)?;
        let Some((lower, upper)) = self.diagonals(n, m)? else {
            return usize::try_from(rectangle).map_err(|_| Error::TooLarge);
        };
        fn above(n: u128, m: u128, threshold: u128) -> Result<u128, Error> {
            let rows = n.saturating_sub(threshold);
            let rising = rows.min(m + 1);
            let triangle = rising.checked_mul(rising + 1).ok_or(Error::TooLarge)? / 2;
            triangle
                .checked_add((rows - rising).checked_mul(m + 1).ok_or(Error::TooLarge)?)
                .ok_or(Error::TooLarge)
        }
        let upper_outside = above(n128, m128, upper as u128)?;
        let lower_outside = above(m128, n128, u128::from(lower.unsigned_abs()))?;
        let cells = rectangle
            .checked_sub(upper_outside)
            .and_then(|count| count.checked_sub(lower_outside))
            .ok_or(Error::TooLarge)?;
        usize::try_from(cells).map_err(|_| Error::TooLarge)
    }
    /// A conservative bound for every complete global path omitted by this domain.
    /// It concerns the declared scoring model, never biological correctness.
    pub fn off_band_bound(
        self,
        n: usize,
        m: usize,
        scoring: Scoring,
    ) -> Result<OffBandBound, Error> {
        if matches!(self, Self::Full) {
            return Ok(OffBandBound::NoEscapingPath);
        }
        // These are the mathematical premises, not a replacement for the native
        // scoring/options validator. Invalid negative-gap premises cannot prove.
        if scoring.gap_open >= 0 || scoring.gap_extend >= 0 {
            return Err(bad("Off-band proof requires strictly negative gap scores"));
        }
        let Self::Corridor { radius } = self else {
            unreachable!()
        };
        let n = i64::try_from(n).map_err(|_| Error::TooLarge)?;
        let m = i64::try_from(m).map_err(|_| Error::TooLarge)?;
        let radius = i64::try_from(radius).map_err(|_| Error::TooLarge)?;
        let total = n.checked_add(m).ok_or(Error::TooLarge)?;
        let delta = n.checked_sub(m).ok_or(Error::TooLarge)?.unsigned_abs();
        let minimum_gaps = i64::try_from(delta)
            .map_err(|_| Error::TooLarge)?
            .checked_add(
                radius
                    .checked_add(1)
                    .and_then(|v| v.checked_mul(2))
                    .ok_or(Error::TooLarge)?,
            )
            .ok_or(Error::TooLarge)?;
        if minimum_gaps > total {
            return Ok(OffBandBound::NoEscapingPath);
        }
        let maximum_pairs = (total - minimum_gaps) / 2;
        // max(0, ...) makes the bound safe independently of a caller's scoring
        // validation; admitted native match scores are already strictly positive.
        let maximum_substitution = i64::from(
            scoring
                .match_score
                .max(scoring.ambiguous_score)
                .max(scoring.mismatch_score)
                .max(0),
        );
        let gap_unit = i64::from(scoring.gap_open.max(scoring.gap_extend));
        let upper = maximum_pairs
            .checked_mul(maximum_substitution)
            .and_then(|paired| {
                minimum_gaps
                    .checked_mul(gap_unit)
                    .and_then(|gaps| paired.checked_add(gaps))
            })
            .ok_or(Error::TooLarge)?;
        Ok(OffBandBound::Finite(upper))
    }
}
data! {
    pub enum OffBandBound {
        Finite(i64),
        NoEscapingPath,
    }
}
impl OffBandBound {
    /// Strict inequality excludes off-band ties too; it does not enumerate
    /// equally scoring paths still inside the admitted domain.
    pub fn certifies(self, selected_score: i32) -> bool {
        match self {
            Self::NoEscapingPath => true,
            Self::Finite(upper) => i64::from(selected_score) > upper,
        }
    }
}
data! {
    pub struct Attempt {
        pub domain: Domain,
        pub evaluated_cells: usize,
        pub score: i32,
        pub off_band_bound: OffBandBound,
        pub unrestricted_optimum_proven: bool,
    }
}
impl Attempt {
    pub fn completed(
        domain: Domain,
        n: usize,
        m: usize,
        scoring: Scoring,
        score: i32,
        evaluated_cells: usize,
    ) -> Result<Self, Error> {
        if score == NEG || evaluated_cells == 0 {
            return Err(bad(
                "A completed attempt requires a reachable scored path and actual cells",
            ));
        }
        let bound = domain.off_band_bound(n, m, scoring)?;
        Ok(Self {
            domain,
            evaluated_cells,
            score,
            off_band_bound: bound,
            unrestricted_optimum_proven: bound.certifies(score),
        })
    }
}
data! {
    pub enum StopReason {
        Certified,
        ExactSymbolicIdentity,
        RadiusLimit,
        AttemptLimit,
        CellLimit,
        MemoryLimit,
        ComparisonIncumbent,
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SearchEvidence {
    pub storage: GlobalStorage,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Attempt, 32>")]
    pub attempts: Vec<Attempt>,
    pub evaluated_cells: usize,
    pub unrestricted_optimum_proven: bool,
    pub stop_reason: StopReason,
}
impl SearchEvidence {
    /// Documentary arithmetic and policy validation, not a new alignment run.
    /// The source/report owner must additionally bind actual strings and scores.
    pub fn validate(
        &self,
        n: usize,
        m: usize,
        mode: AlignmentMode,
        scoring: Scoring,
        selected_score: i32,
        max_cells: usize,
    ) -> Result<(), Error> {
        self.storage.validate(mode)?;
        if max_cells == 0
            || max_cells > super::HARD_MAX_CELLS
            || selected_score.unsigned_abs() > 500_000_000
        {
            return Err(bad("Invalid alignment evidence score or cell cap"));
        }
        if self.stop_reason == StopReason::ExactSymbolicIdentity {
            if !self.attempts.is_empty()
                || self.evaluated_cells != 0
                || !self.unrestricted_optimum_proven
                || mode != AlignmentMode::Global
                || n != m
                || scoring.ambiguous_score < 0
            {
                return Err(bad("Invalid exact-symbolic-identity evidence"));
            }
            return Ok(());
        }
        let (maximum, mut expected_domain) = match self.storage {
            GlobalStorage::Full => (1, Domain::Full),
            GlobalStorage::Adaptive {
                initial_radius,
                max_attempts,
                ..
            } => (
                max_attempts,
                Domain::Corridor {
                    radius: initial_radius.min(n.min(m)),
                },
            ),
        };
        if self.attempts.is_empty() || self.attempts.len() > maximum {
            return Err(bad("Invalid alignment attempt count"));
        }
        let mut cells = 0usize;
        let mut preceding_score = None;
        for (index, attempt) in self.attempts.iter().enumerate() {
            if attempt.domain != expected_domain
                || attempt.score.unsigned_abs() > 500_000_000
                || attempt.evaluated_cells != attempt.domain.cell_count(n, m)?
                || preceding_score.is_some_and(|score| attempt.score < score)
            {
                return Err(bad("Alignment attempt domain, cell count or score changed"));
            }
            let expected = Attempt::completed(
                attempt.domain,
                n,
                m,
                scoring,
                attempt.score,
                attempt.evaluated_cells,
            )?;
            if attempt != &expected
                || (index + 1 < self.attempts.len() && attempt.unrestricted_optimum_proven)
            {
                return Err(bad("Alignment attempt proof or continuation changed"));
            }
            cells = cells
                .checked_add(attempt.evaluated_cells)
                .filter(|cells| *cells <= max_cells)
                .ok_or(Error::TooLarge)?;
            preceding_score = Some(attempt.score);
            if let (GlobalStorage::Adaptive { maximum_radius, .. }, Domain::Corridor { radius }) =
                (self.storage, expected_domain)
            {
                let next = radius
                    .checked_add(1)
                    .and_then(|one| radius.checked_mul(2).map(|two| one.max(two)))
                    .ok_or(Error::TooLarge)?
                    .min(maximum_radius)
                    .min(n.min(m));
                if index + 1 < self.attempts.len() && next <= radius {
                    return Err(bad("Alignment attempted the same band again"));
                }
                expected_domain = Domain::Corridor { radius: next };
            }
        }
        let last = self.attempts.last().expect("nonempty attempt evidence");
        if cells != self.evaluated_cells
            || selected_score != last.score
            || self.unrestricted_optimum_proven != last.unrestricted_optimum_proven
            || (self.stop_reason == StopReason::Certified) != self.unrestricted_optimum_proven
        {
            return Err(bad(
                "Selected alignment attempt or cumulative proof changed",
            ));
        }
        match (self.stop_reason, self.storage, last.domain) {
            (
                StopReason::RadiusLimit,
                GlobalStorage::Adaptive { maximum_radius, .. },
                Domain::Corridor { radius },
            ) if radius == maximum_radius.min(n.min(m)) => {}
            (StopReason::AttemptLimit, GlobalStorage::Adaptive { max_attempts, .. }, _)
                if self.attempts.len() == max_attempts => {}
            (
                StopReason::Certified
                | StopReason::CellLimit
                | StopReason::MemoryLimit
                | StopReason::ComparisonIncumbent,
                _,
                _,
            ) => {}
            _ => return Err(bad("Alignment stop reason disagrees with admitted policy")),
        }
        Ok(())
    }
}
data! {
    pub enum CandidateDisposition {
        CertifiedOwnOptimum,
        DischargedByIncumbent,
        Limited,
        Unsearched,
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CandidateEvidence {
    pub reference_origin: usize,
    pub query_reversed: bool,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub score: Option<i32>,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub search: Option<SearchEvidence>,
    pub disposition: CandidateDisposition,
}
/// One scientific cell ledger. It is not serialized or reset by a band retry.
pub(crate) struct CellWork {
    limit: usize,
    used: usize,
}
impl CellWork {
    pub fn new(limit: usize) -> Result<Self, Error> {
        if limit == 0 || limit > 2 * super::HARD_MAX_CELLS {
            return Err(Error::TooLarge);
        }
        Ok(Self { limit, used: 0 })
    }
    pub fn used(&self) -> usize {
        self.used
    }
    pub fn remaining(&self) -> usize {
        self.limit - self.used
    }
}
/// One alignment's cap borrows the whole-comparison ledger. Creating another
/// attempt does not create another budget or return already evaluated cells.
pub(crate) struct CellBudget<'a> {
    shared: &'a mut CellWork,
    limit: usize,
    used: usize,
}
impl<'a> CellBudget<'a> {
    pub fn new(shared: &'a mut CellWork, limit: usize) -> Result<Self, Error> {
        Self::resume(shared, limit, 0)
    }
    pub fn resume(shared: &'a mut CellWork, limit: usize, used: usize) -> Result<Self, Error> {
        if limit == 0 || limit > super::HARD_MAX_CELLS || used > limit {
            return Err(Error::TooLarge);
        }
        Ok(Self {
            shared,
            limit,
            used,
        })
    }
    pub fn used(&self) -> usize {
        self.used
    }
    pub fn remaining(&self) -> usize {
        (self.limit - self.used).min(self.shared.remaining())
    }
    pub fn admit(&self, cells: usize) -> Result<(), Error> {
        if cells > self.remaining() {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
    pub fn charge_row(&mut self, cells: usize) -> Result<(), Error> {
        self.admit(cells)?;
        self.used += cells;
        self.shared.used += cells;
        Ok(())
    }
}
/// Layout owns its original-coordinate addressing and preallocation admission.
/// A rectangular layout needs no per-row offset allocation.
pub(crate) struct RowLayout {
    n: usize,
    m: usize,
    domain: Domain,
    offsets: Vec<usize>,
    cells: usize,
    maximum_width: usize,
    bytes: usize,
}
impl RowLayout {
    pub fn new(n: usize, m: usize, domain: Domain, remaining_cells: usize) -> Result<Self, Error> {
        let admitted_cells = domain.cell_count(n, m)?;
        if admitted_cells > remaining_cells {
            return Err(Error::TooLarge);
        }
        let rows = n.checked_add(1).ok_or(Error::TooLarge)?;
        let full_width = m.checked_add(1).ok_or(Error::TooLarge)?;
        if rows > remaining_cells {
            return Err(Error::TooLarge);
        }
        let (offsets, cells, maximum_width) = if matches!(domain, Domain::Full) {
            (
                vec![],
                rows.checked_mul(full_width).ok_or(Error::TooLarge)?,
                full_width,
            )
        } else {
            // Admit metadata before its first allocation. Every corridor row is
            // nonempty, so row count alone is already a lower cell-work bound.
            let offset_count = rows.checked_add(1).ok_or(Error::TooLarge)?;
            if offset_count
                .checked_mul(size_of::<usize>())
                .ok_or(Error::TooLarge)?
                > MAX_MATRIX_BYTES
            {
                return Err(Error::TooLarge);
            }
            let mut offsets = Vec::new();
            offsets
                .try_reserve_exact(offset_count)
                .map_err(|_| Error::TooLarge)?;
            offsets.push(0usize);
            let mut cells = 0usize;
            let mut width = 0usize;
            for i in 0..rows {
                let (start, end) = domain.row(i, n, m)?;
                let length = end
                    .checked_sub(start)
                    .and_then(|v| v.checked_add(1))
                    .ok_or(Error::TooLarge)?;
                cells = cells
                    .checked_add(length)
                    .filter(|value| *value <= remaining_cells)
                    .ok_or(Error::TooLarge)?;
                width = width.max(length);
                offsets.push(cells);
            }
            (offsets, cells, width)
        };
        if cells != admitted_cells || cells > remaining_cells {
            return Err(Error::TooLarge);
        }
        let bytes = maximum_width
            .checked_mul(2 * size_of::<[i32; 3]>())
            .and_then(|scores| cells.checked_add(scores))
            .and_then(|bytes| {
                offsets
                    .len()
                    .checked_mul(size_of::<usize>())
                    .and_then(|offsets| bytes.checked_add(offsets))
            })
            .filter(|bytes| *bytes <= MAX_MATRIX_BYTES)
            .ok_or(Error::TooLarge)?;
        Ok(Self {
            n,
            m,
            domain,
            offsets,
            cells,
            maximum_width,
            bytes,
        })
    }
    pub fn cells(&self) -> usize {
        self.cells
    }
    pub fn bytes(&self) -> usize {
        self.bytes
    }
    pub fn row(&self, i: usize) -> Result<(usize, usize), Error> {
        self.domain.row(i, self.n, self.m)
    }
    pub fn row_width(&self, i: usize) -> Result<usize, Error> {
        let (a, b) = self.row(i)?;
        Ok(b - a + 1)
    }
    fn index(&self, i: usize, j: usize) -> Result<usize, Error> {
        let (start, end) = self.row(i)?;
        if j < start || j > end {
            return Err(bad("Traceback coordinate is outside admitted row domain"));
        }
        let offset = if self.offsets.is_empty() {
            i.checked_mul(self.m + 1).ok_or(Error::TooLarge)?
        } else {
            self.offsets[i]
        };
        offset
            .checked_add(j - start)
            .filter(|value| *value < self.cells)
            .ok_or(Error::TooLarge)
    }
}
pub(crate) struct ScoreRow {
    start: usize,
    length: usize,
    values: Vec<[i32; 3]>,
}
impl ScoreRow {
    fn new(width: usize) -> Result<Self, Error> {
        let mut values = Vec::new();
        values
            .try_reserve_exact(width)
            .map_err(|_| Error::TooLarge)?;
        values.resize(width, [NEG; 3]);
        Ok(Self {
            start: 0,
            length: 0,
            values,
        })
    }
    fn reset(&mut self, start: usize, end: usize) -> Result<(), Error> {
        let length = end
            .checked_sub(start)
            .and_then(|v| v.checked_add(1))
            .ok_or(Error::TooLarge)?;
        if length > self.values.len() {
            return Err(bad("Score row exceeds admitted allocation"));
        }
        self.start = start;
        self.length = length;
        self.values[..length].fill([NEG; 3]);
        Ok(())
    }
    pub fn get(&self, j: usize) -> [i32; 3] {
        j.checked_sub(self.start)
            .filter(|offset| *offset < self.length)
            .map_or([NEG; 3], |offset| self.values[offset])
    }
    pub fn set(&mut self, j: usize, value: [i32; 3]) -> Result<(), Error> {
        let offset = j
            .checked_sub(self.start)
            .filter(|offset| *offset < self.length)
            .ok_or_else(|| bad("Score write is outside admitted row"))?;
        self.values[offset] = value;
        Ok(())
    }
}
/// Storage only: the single native recurrence will evaluate and charge rows.
pub(crate) struct Grid {
    layout: RowLayout,
    trace: Vec<u8>,
    previous: ScoreRow,
    current: ScoreRow,
}
impl Grid {
    pub fn new(layout: RowLayout) -> Result<Self, Error> {
        if layout.bytes() > MAX_MATRIX_BYTES {
            return Err(Error::TooLarge);
        }
        let mut trace = Vec::new();
        trace
            .try_reserve_exact(layout.cells)
            .map_err(|_| Error::TooLarge)?;
        trace.resize(layout.cells, 0x3f);
        let previous = ScoreRow::new(layout.maximum_width)?;
        let current = ScoreRow::new(layout.maximum_width)?;
        Ok(Self {
            layout,
            trace,
            previous,
            current,
        })
    }
    pub fn layout(&self) -> &RowLayout {
        &self.layout
    }
    pub fn begin(&mut self, i: usize) -> Result<(), Error> {
        let (start, end) = self.layout.row(i)?;
        self.current.reset(start, end)
    }
    pub fn rows(&mut self) -> (&ScoreRow, &mut ScoreRow) {
        (&self.previous, &mut self.current)
    }
    pub fn finish_row(&mut self) {
        std::mem::swap(&mut self.previous, &mut self.current);
    }
    pub fn final_row(&self) -> &ScoreRow {
        &self.previous
    }
    pub fn set_predecessor(
        &mut self,
        i: usize,
        j: usize,
        state: u8,
        previous: u8,
    ) -> Result<(), Error> {
        if state > 2 || previous > 3 {
            return Err(bad("Invalid affine predecessor state"));
        }
        let index = self.layout.index(i, j)?;
        let shift = state * 2;
        self.trace[index] = (self.trace[index] & !(3 << shift)) | (previous << shift);
        Ok(())
    }
    pub fn predecessor(&self, i: usize, j: usize, state: u8) -> Result<u8, Error> {
        if state > 2 {
            return Err(bad("Invalid affine traceback state"));
        }
        Ok((self.trace[self.layout.index(i, j)?] >> (state * 2)) & 3)
    }
}

#[cfg(test)]
#[path = "grid_tests.rs"]
mod tests;
