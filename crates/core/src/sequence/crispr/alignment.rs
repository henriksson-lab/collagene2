//! Bounded global guide/target alignments, independent of PAM and record geometry.
//!
//! The supplied DNA is already oriented like the guide (T/U equivalent), not the
//! complementary strand. Every base in both inputs is consumed. Terminal gaps
//! and mixed gap kinds are included; callers must declare any narrower geometry.
//! DNA-bulge = extra target DNA, RNA-bulge = extra guide RNA, following Lin et al.
//! https://pmc.ncbi.nlm.nih.gov/articles/PMC4066799/ . These names describe sequence
//! alignments, not an assertion about physical heteroduplex structure or cleavage.
use crate::sequence::analysis::base_mask;
use collagene_model::{Alphabet, Error, Sequence};
use serde::{Deserialize, Serialize};

pub const MODEL_ID: &str = "global-full-input-paired-dna-rna-1";

const MAX_GUIDE: usize = 64;
const MAX_BULGE_BASES: usize = 4;
const MAX_ALIGNMENTS: usize = 4096;
const MAX_WORK: usize = 500_000_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AlignmentLimits {
    pub max_mismatches: usize,
    pub max_dna_bulge_bases: usize,
    pub max_rna_bulge_bases: usize,
    pub max_total_bulge_bases: usize,
    pub max_alignments: usize,
}
impl AlignmentLimits {
    fn validate(self, guide_len: usize) -> Result<(), Error> {
        if self.max_mismatches > guide_len
            || self.max_dna_bulge_bases > MAX_BULGE_BASES
            || self.max_rna_bulge_bases > MAX_BULGE_BASES
            || self.max_total_bulge_bases > MAX_BULGE_BASES
            || !(1..=MAX_ALIGNMENTS).contains(&self.max_alignments)
        {
            return Err(Error::Unsupported("CRISPR alignment limits".into()));
        }
        Ok(())
    }
}

/// Shared across all candidate intervals in a job. No reset on each window.
#[derive(Debug)]
pub struct AlignmentWork {
    limit: usize,
    used: usize,
    column_limit: usize,
    admitted_columns: usize,
}
impl AlignmentWork {
    pub fn new(limit: usize, column_limit: usize) -> Result<Self, Error> {
        if limit > MAX_WORK || column_limit > 1_000_000 {
            return Err(Error::TooLarge);
        }
        Ok(Self {
            limit,
            used: 0,
            column_limit,
            admitted_columns: 0,
        })
    }
    pub fn used(&self) -> usize {
        self.used
    }
    pub(crate) fn spend(&mut self, units: usize) -> bool {
        if units > self.limit - self.used {
            self.used = self.limit;
            false
        } else {
            self.used += units;
            true
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, PartialOrd, Ord)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum AlignmentColumn {
    /// Indices into the supplied guide and oriented DNA, both zero based.
    Paired {
        guide: usize,
        target: usize,
    },
    DnaBulge {
        target: usize,
    },
    RnaBulge {
        guide: usize,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct GuideAlignment {
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, AlignmentColumn, 132>")]
    pub columns: Vec<AlignmentColumn>,
    pub min_mismatches: usize,
    pub max_mismatches: usize,
    pub dna_bulge_bases: usize,
    pub rna_bulge_bases: usize,
    /// Every independent IUPAC realization satisfies the mismatch threshold for
    /// this particular alignment. Bulge bases are counted separately.
    pub definite: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum AlignmentLimitReached {
    Work,
    Alignments,
    RetainedColumns,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AlignmentEvidence {
    /// All retained possible alignments, in paired / DNA-gap / RNA-gap DFS order.
    /// Repeated bases may yield multiple equally valid gap placements.
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, GuideAlignment, 4096>")]
    pub alignments: Vec<GuideAlignment>,
    /// True only when all paths within the declared bounds have been examined.
    /// An empty incomplete result MUST NOT be interpreted as absence.
    pub complete: bool,
    pub limit_reached: Option<AlignmentLimitReached>,
    /// Work consumed by this interval, not the entire shared job.
    pub work: usize,
}

#[derive(Clone, Copy)]
struct Position {
    guide: usize,
    target: usize,
    minimum: usize,
    maximum: usize,
    dna: usize,
    rna: usize,
}

struct Search<'a> {
    guide: &'a [u8],
    target: &'a [u8],
    limits: AlignmentLimits,
    work: &'a mut AlignmentWork,
    path: Vec<AlignmentColumn>,
    evidence: AlignmentEvidence,
}
impl Search<'_> {
    fn stop(&mut self, reason: AlignmentLimitReached) {
        self.evidence.complete = false;
        self.evidence.limit_reached = Some(reason);
    }

    fn visit(&mut self, p: Position) {
        if !self.evidence.complete {
            return;
        }
        // One unit covers a constant-size state, including pruning and edges.
        if !self.work.spend(1) {
            self.stop(AlignmentLimitReached::Work);
            return;
        }
        if p.minimum > self.limits.max_mismatches {
            return;
        }
        let remaining_guide = self.guide.len() - p.guide;
        let remaining_target = self.target.len() - p.target;
        let needed_dna = remaining_target.saturating_sub(remaining_guide);
        let needed_rna = remaining_guide.saturating_sub(remaining_target);
        if needed_dna > self.limits.max_dna_bulge_bases - p.dna
            || needed_rna > self.limits.max_rna_bulge_bases - p.rna
            || needed_dna + needed_rna > self.limits.max_total_bulge_bases - p.dna - p.rna
        {
            return;
        }
        if remaining_guide == 0 && remaining_target == 0 {
            // Only an additional actual qualifying path makes the cap incomplete.
            if self.evidence.alignments.len() == self.limits.max_alignments {
                self.stop(AlignmentLimitReached::Alignments);
                return;
            }
            if self.path.len()
                > self
                    .work
                    .column_limit
                    .saturating_sub(self.work.admitted_columns)
            {
                self.stop(AlignmentLimitReached::RetainedColumns);
                return;
            }
            if !self.work.spend(self.path.len()) {
                self.stop(AlignmentLimitReached::Work);
                return;
            }
            self.work.admitted_columns += self.path.len();
            self.evidence.alignments.push(GuideAlignment {
                columns: self.path.clone(),
                min_mismatches: p.minimum,
                max_mismatches: p.maximum,
                dna_bulge_bases: p.dna,
                rna_bulge_bases: p.rna,
                definite: p.maximum <= self.limits.max_mismatches,
            });
            return;
        }
        if remaining_guide > 0 && remaining_target > 0 {
            let q = self.guide[p.guide];
            let t = self.target[p.target];
            self.path.push(AlignmentColumn::Paired {
                guide: p.guide,
                target: p.target,
            });
            self.visit(Position {
                guide: p.guide + 1,
                target: p.target + 1,
                minimum: p.minimum + usize::from(q & t == 0),
                // Even N/N can disagree: only equal concrete calls are certain.
                maximum: p.maximum + usize::from(q != t || !q.is_power_of_two()),
                ..p
            });
            self.path.pop();
        }
        if p.dna + p.rna < self.limits.max_total_bulge_bases {
            if remaining_target > 0 && p.dna < self.limits.max_dna_bulge_bases {
                self.path
                    .push(AlignmentColumn::DnaBulge { target: p.target });
                self.visit(Position {
                    target: p.target + 1,
                    dna: p.dna + 1,
                    ..p
                });
                self.path.pop();
            }
            if remaining_guide > 0 && p.rna < self.limits.max_rna_bulge_bases {
                self.path.push(AlignmentColumn::RnaBulge { guide: p.guide });
                self.visit(Position {
                    guide: p.guide + 1,
                    rna: p.rna + 1,
                    ..p
                });
                self.path.pop();
            }
        }
    }
}

/// Align exactly one supplied interval, never a reference-wide search.
///
/// For each paired column, minimum mismatch is one iff the IUPAC sets are
/// disjoint; maximum is zero iff both sets are the same singleton. Realizations
/// are independent at different sequence positions. Gaps are separate from
/// mismatch counts. Every path whose minimum satisfies the threshold is retained
/// until an explicit bound is reached; no scoring or best-path selection occurs.
pub fn align(
    guide: &Sequence,
    oriented_target: &Sequence,
    limits: AlignmentLimits,
    work: &mut AlignmentWork,
) -> Result<AlignmentEvidence, Error> {
    input_shape(guide, oriented_target)?;
    limits.validate(guide.len())?;
    let before = work.used();
    let mut evidence = AlignmentEvidence {
        alignments: Vec::new(),
        complete: true,
        limit_reached: None,
        work: 0,
    };
    // Bound preprocessing as well as exploration; no allocation on exhausted jobs.
    if !work.spend(guide.len() + oriented_target.len()) {
        evidence.complete = false;
        evidence.limit_reached = Some(AlignmentLimitReached::Work);
        evidence.work = work.used() - before;
        return Ok(evidence);
    }
    let guide = guide.as_str().bytes().map(base_mask).collect::<Vec<_>>();
    let target = oriented_target
        .as_str()
        .bytes()
        .map(base_mask)
        .collect::<Vec<_>>();
    let mut search = Search {
        guide: &guide,
        target: &target,
        limits,
        work,
        path: Vec::with_capacity(MAX_GUIDE + MAX_BULGE_BASES),
        evidence,
    };
    search.visit(Position {
        guide: 0,
        target: 0,
        minimum: 0,
        maximum: 0,
        dna: 0,
        rna: 0,
    });
    search.evidence.work = search.work.used() - before;
    Ok(search.evidence)
}

#[cfg(test)]
mod tests;

fn input_shape(guide: &Sequence, oriented_target: &Sequence) -> Result<(), Error> {
    if !matches!(guide.alphabet(), Alphabet::Dna | Alphabet::Rna)
        || oriented_target.alphabet() != Alphabet::Dna
        || !(1..=MAX_GUIDE).contains(&guide.len())
        || !(1..=MAX_GUIDE + MAX_BULGE_BASES).contains(&oriented_target.len())
    {
        return Err(Error::Unsupported(
            "CRISPR alignment requires bounded nucleotide guide and oriented DNA".into(),
        ));
    }
    Ok(())
}
impl GuideAlignment {
    /// Retained path integrity; current genomic/PAM authority requires source replay.
    pub fn validate_for(
        &self,
        guide: &Sequence,
        target: &Sequence,
        limits: AlignmentLimits,
    ) -> Result<(), Error> {
        input_shape(guide, target)?;
        limits.validate(guide.len())?;
        if self.columns.is_empty()
            || self.columns.len() > guide.len() + target.len()
            || self.columns.len() > 132
        {
            return Err(Error::TooLarge);
        }
        let (mut qi, mut ti, mut dna, mut rna, mut minimum, mut maximum) = (0, 0, 0, 0, 0, 0);
        for column in &self.columns {
            match *column {
                AlignmentColumn::Paired {
                    guide: g,
                    target: t,
                } if g == qi && t == ti && qi < guide.len() && ti < target.len() => {
                    let q = base_mask(guide.as_str().as_bytes()[qi]);
                    let t = base_mask(target.as_str().as_bytes()[ti]);
                    minimum += usize::from(q & t == 0);
                    maximum += usize::from(q != t || !q.is_power_of_two());
                    qi += 1;
                    ti += 1;
                }
                AlignmentColumn::DnaBulge { target: t } if t == ti && ti < target.len() => {
                    ti += 1;
                    dna += 1;
                }
                AlignmentColumn::RnaBulge { guide: g } if g == qi && qi < guide.len() => {
                    qi += 1;
                    rna += 1;
                }
                _ => {
                    return Err(Error::Unsupported(
                        "CRISPR alignment path consumes invalid indices".into(),
                    ));
                }
            }
        }
        if qi != guide.len()
            || ti != target.len()
            || dna > limits.max_dna_bulge_bases
            || rna > limits.max_rna_bulge_bases
            || dna + rna > limits.max_total_bulge_bases
            || minimum > limits.max_mismatches
            || (minimum, maximum, dna, rna)
                != (
                    self.min_mismatches,
                    self.max_mismatches,
                    self.dna_bulge_bases,
                    self.rna_bulge_bases,
                )
            || self.definite != (maximum <= limits.max_mismatches)
        {
            return Err(Error::Unsupported(
                "CRISPR retained alignment bounds/counts differ".into(),
            ));
        }
        Ok(())
    }
}
impl AlignmentEvidence {
    pub fn validate_for(
        &self,
        guide: &Sequence,
        target: &Sequence,
        limits: AlignmentLimits,
    ) -> Result<usize, Error> {
        input_shape(guide, target)?;
        limits.validate(guide.len())?;
        if self.alignments.len() > limits.max_alignments
            || self.complete != self.limit_reached.is_none()
            || self.work > MAX_WORK
        {
            return Err(Error::Unsupported(
                "CRISPR alignment evidence completeness/bounds".into(),
            ));
        }
        let mut columns = 0usize;
        let mut paths = std::collections::BTreeSet::new();
        for path in &self.alignments {
            path.validate_for(guide, target, limits)?;
            if !paths.insert(&path.columns) {
                return Err(Error::Unsupported("duplicate CRISPR alignment path".into()));
            }
            columns = columns
                .checked_add(path.columns.len())
                .ok_or(Error::TooLarge)?;
        }
        Ok(columns)
    }
}
