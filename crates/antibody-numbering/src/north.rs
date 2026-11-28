//! North–Lehmann–Dunbrack (2011) structural CDR definitions, conditional on
//! independently established framework correspondence, not displayed numbering.
use crate::{
    Error, Insertion, MAX_SOURCE_RESIDUES, Position, Result,
    alternate_tracks::{Issue, WindowIssue},
    domain_tracks::{
        AssignmentEvidence, CdrRegion, CdrStatus, NumberingScheme, Region, Residue, imgt_region,
    },
    hash,
    numbering_windows::Window,
    v_mapping::Locus,
};
use std::ops::Range;

pub fn catalog_revision() -> String {
    hash(b"antibody-north-region-rules-1\0", &(
        "https://doi.org/10.1016/j.jmb.2010.10.030; Defining complementarity determining regions; Figures1,11",
        "https://www.imgt.org/PDF/DCI/27_55-77_2003.pdf; Table2 pp61-62",
        "IMGT first24..40; heavy second55..66; light second55..69; third105..117; light second includes preceding residue",
        "source intervals from invariant V framework flanks23/41,54/67heavy or54/70light,104/J118; left excluded and right excluded; no interior numbering requirement; missing boundaries unresolved; known intervals and intervening framework source order; no AHo numbering conversion",
    ),65536).expect("fixed North definition rules")
}

/// Membership for an already established IMGT label. This checks syntax, not
/// antibody identity or numbering evidence; a bare label grants no authority.
pub fn region(locus: &Locus, label: &Position) -> Result<Region> {
    imgt_region(label)?;
    Ok(match label.number {
        1..=23 => Region::Fr1,
        24..=40 => Region::Cdr1,
        41..=54 => Region::Fr2,
        n if n <= if *locus == Locus::Heavy { 66 } else { 69 } => Region::Cdr2,
        67..=104 => Region::Fr3,
        105..=117 => Region::Cdr3,
        118..=128 => Region::Fr4,
        _ => return Err(Error::Invalid("unsupported North IMGT position")),
    })
}

/// Runtime source-bound region projection. No serialized claim can construct
/// this owner; the caller supplies freshly reconstructed IMGT correspondence.
pub(crate) struct Projection {
    pub regions: Vec<Option<Region>>,
    pub cdrs: Vec<CdrRegion>,
    pub issues: Vec<Issue>,
}
impl Projection {
    fn anchor(rows: &[Residue], number: u16) -> Option<usize> {
        rows.iter()
            .find(|row| {
                row.assignment.as_ref().is_some_and(|assignment| {
                    assignment.position
                        == Position {
                            number,
                            insertion: Insertion::None,
                        }
                        && assignment.evidence
                            == if number == 118 {
                                AssignmentEvidence::JFr4Reference
                            } else {
                                AssignmentEvidence::InvariantVReference
                            }
                })
            })
            .map(|row| row.source_index)
    }
    fn assign(&mut self, rows: &[Residue], span: Range<usize>, kind: Region) -> Result<()> {
        for (row, region) in rows.iter().zip(&mut self.regions) {
            if span.contains(&row.source_index) {
                if region.as_ref().is_some_and(|old| old != &kind) {
                    return Err(Error::Binding);
                }
                *region = Some(kind.clone());
            }
        }
        Ok(())
    }
    pub(crate) fn build(
        rows: &[Residue],
        locus: &Locus,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        if rows.is_empty() || rows.len() > MAX_SOURCE_RESIDUES {
            return Err(Error::Limit("North source rows"));
        }
        let start = rows[0].source_index;
        let end = start.checked_add(rows.len()).ok_or(Error::Binding)?;
        if end > MAX_SOURCE_RESIDUES {
            return Err(Error::Limit("North source coordinates"));
        }
        let mut regions = Vec::with_capacity(rows.len());
        for (i, row) in rows.iter().enumerate() {
            if cancel() {
                return Err(Error::Cancelled);
            }
            if row.source_index != start + i || row.residue == '*' {
                return Err(Error::Binding);
            }
            regions.push(
                row.assignment
                    .as_ref()
                    .map(|a| region(locus, &a.position))
                    .transpose()?,
            );
        }
        let mut out = Self {
            regions,
            cdrs: Vec::with_capacity(3),
            issues: vec![],
        };
        let mut spans: [Option<Range<usize>>; 3] = [None, None, None];
        for (index, (kind, window, left, right)) in [
            (Region::Cdr1, Window::CysToTrp, 23, 41),
            (
                Region::Cdr2,
                Window::MiddleLoop,
                54,
                if *locus == Locus::Heavy { 67 } else { 70 },
            ),
            (Region::Cdr3, Window::CysToJ, 104, 118),
        ]
        .into_iter()
        .enumerate()
        {
            if cancel() {
                return Err(Error::Cancelled);
            }
            let a = Self::anchor(rows, left);
            let b = Self::anchor(rows, right);
            let status = if let Some((a, b)) = a.zip(b) {
                if a >= b {
                    return Err(Error::Binding);
                }
                if a + 1 == b {
                    out.issues.push(Issue {
                        numbering_scheme: NumberingScheme::Imgt,
                        window,
                        issue: WindowIssue::UnsupportedLength { residues: 0 },
                    });
                    CdrStatus::Unresolved
                } else {
                    let span = a + 1..b;
                    out.assign(rows, span.clone(), kind.clone())?;
                    spans[index] = Some(span.clone());
                    CdrStatus::Mapped {
                        start: span.start,
                        end: span.end,
                    }
                }
            } else {
                let mut positions = vec![];
                if a.is_none() {
                    positions.push(left);
                }
                if b.is_none() {
                    positions.push(right);
                }
                out.issues.push(Issue {
                    numbering_scheme: NumberingScheme::Imgt,
                    window,
                    issue: WindowIssue::MissingInvariantPositions { positions },
                });
                CdrStatus::Unresolved
            };
            out.cdrs.push(CdrRegion { kind, status });
        }
        if let Some(first) = &spans[0] {
            out.assign(rows, start..first.start, Region::Fr1)?;
        }
        for (index, kind) in [(0, Region::Fr2), (1, Region::Fr3)] {
            if let (Some(left), Some(right)) = (&spans[index], &spans[index + 1]) {
                if left.end > right.start {
                    return Err(Error::Binding);
                }
                out.assign(rows, left.end..right.start, kind)?;
            }
        }
        if let Some(last) = &spans[2] {
            out.assign(rows, last.end..end, Region::Fr4)?;
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain_tracks::Assignment;
    fn rows() -> Vec<Residue> {
        // A projection-kernel fixture, not an automatically detected antibody.
        // Six unnumbered interior residues are enclosed by explicit reference
        // flanks; each independent literal two-residue loop must be retained.
        "CASWIAGYCRDWG"
            .chars()
            .enumerate()
            .map(|(source_index, residue)| {
                let number = match source_index {
                    0 => Some(23),
                    3 => Some(41),
                    4 => Some(54),
                    7 => Some(67),
                    8 => Some(104),
                    11 => Some(118),
                    12 => Some(119),
                    _ => None,
                };
                Residue {
                    source_index,
                    residue,
                    region: None,
                    assignment: number.map(|number| Assignment {
                        position: Position {
                            number,
                            insertion: Insertion::None,
                        },
                        evidence: if number >= 118 {
                            AssignmentEvidence::JFr4Reference
                        } else {
                            AssignmentEvidence::InvariantVReference
                        },
                    }),
                }
            })
            .collect()
    }
    #[test]
    fn literal_source_scopes_do_not_require_interior_number_labels() {
        let projection = Projection::build(&rows(), &Locus::Heavy, &mut || false).unwrap();
        assert_eq!(
            projection
                .cdrs
                .iter()
                .map(|r| r.status.clone())
                .collect::<Vec<_>>(),
            [
                CdrStatus::Mapped { start: 1, end: 3 },
                CdrStatus::Mapped { start: 5, end: 7 },
                CdrStatus::Mapped { start: 9, end: 11 },
            ]
        );
        assert_eq!(
            projection.regions,
            [
                Some(Region::Fr1),
                Some(Region::Cdr1),
                Some(Region::Cdr1),
                Some(Region::Fr2),
                Some(Region::Fr2),
                Some(Region::Cdr2),
                Some(Region::Cdr2),
                Some(Region::Fr3),
                Some(Region::Fr3),
                Some(Region::Cdr3),
                Some(Region::Cdr3),
                Some(Region::Fr4),
                Some(Region::Fr4)
            ]
        );
        assert!(projection.issues.is_empty());
    }
    #[test]
    fn missing_or_forged_flank_does_not_create_whole_loop_coverage() {
        for assignment in [
            None,
            Some(Assignment {
                position: Position {
                    number: 54,
                    insertion: Insertion::None,
                },
                evidence: AssignmentEvidence::ConditionalLoopRule,
            }),
        ] {
            let mut rows = rows();
            rows[4].assignment = assignment;
            let projection = Projection::build(&rows, &Locus::Heavy, &mut || false).unwrap();
            assert_eq!(
                projection.cdrs[0].status,
                CdrStatus::Mapped { start: 1, end: 3 }
            );
            assert_eq!(projection.cdrs[1].status, CdrStatus::Unresolved);
            assert_eq!(
                projection.cdrs[2].status,
                CdrStatus::Mapped { start: 9, end: 11 }
            );
            assert_eq!(projection.regions[5], None);
            assert!(
                matches!(&projection.issues[0].issue,WindowIssue::MissingInvariantPositions{positions} if positions==&[54])
            );
        }
    }
    #[test]
    fn interrupted_coordinates_stops_and_cancelled_work_never_form_a_projection() {
        let mut gap = rows();
        gap[2].source_index = 99;
        assert_eq!(
            Projection::build(&gap, &Locus::Heavy, &mut || false).err(),
            Some(Error::Binding)
        );
        let mut stop = rows();
        stop[5].residue = '*';
        assert_eq!(
            Projection::build(&stop, &Locus::Heavy, &mut || false).err(),
            Some(Error::Binding)
        );
        assert_eq!(
            Projection::build(&rows(), &Locus::Heavy, &mut || true).err(),
            Some(Error::Cancelled)
        );
    }
}
