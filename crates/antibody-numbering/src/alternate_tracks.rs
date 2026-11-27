//! Conditional alternate labels from invariant framework correspondences and
//! verified window contexts. CDR definitions are chosen independently.
use super::*;
use crate::{domain_candidates as d, domain_tracks::*, numbering_windows as w, v_mapping::Locus};
data! {pub enum WindowIssue {
    MissingInvariantPositions{positions:Vec<u16>},NonContiguousFlanks,
    UnsupportedLength{residues:usize},NonCanonicalAnchors,
}}
data! {pub struct Issue {pub numbering_scheme:NumberingScheme,pub window:w::Window,pub issue:WindowIssue}}
pub fn catalog_revision() -> String {
    hash(b"antibody-alternate-track-rules-1\0", &(
        w::catalog_revision(),
        "https://www.imgt.org/PDF/DCI/27_55-77_2003.pdf; Table2 framework correspondences",
        "heavy I1..9=K1..9 I11..26=K10..25 I41..55=K36..50 I66..72=K58..64 I74..91=K65..82 I92..94=K82ABC I95..104=K83..92; light I1..26=K1..26 I39..55=K33..49 I66..72=K53..59 I74..80=K60..66 I83..104=K67..88; J heavy103+offset light98+offset; invariant references only; CysTrp requires adjacent23/24/25/26 and39/40/41; no invented framework insertions; IMGT and Kabat definitions via independent source labels",
    ),65536).expect("fixed alternate track rules")
}
fn fixed(locus: &Locus, p: &Position) -> Option<Position> {
    if p.insertion != Insertion::None {
        return None;
    }
    let n = p.number;
    let target = if *locus == Locus::Heavy {
        match n {
            1..=9 => n,
            11..=26 => n - 1,
            41..=55 => n - 5,
            66..=72 => n - 8,
            74..=91 => n - 9,
            92..=94 => {
                return Some(Position {
                    number: 82,
                    insertion: Insertion::Letters(char::from(b'A' + (n - 92) as u8).to_string()),
                });
            }
            95..=104 => n - 12,
            _ => return None,
        }
    } else {
        match n {
            1..=26 => n,
            39..=55 => n - 6,
            66..=72 => n - 13,
            74..=80 => n - 14,
            83..=104 => n - 16,
            _ => return None,
        }
    };
    Some(position(target))
}
fn reference_index(rows: &[Residue], number: u16) -> Option<usize> {
    rows.iter()
        .find(|r| {
            r.assignment.as_ref().is_some_and(|a| {
                a.evidence == AssignmentEvidence::InvariantVReference
                    && a.position == position(number)
            })
        })
        .map(|r| r.source_index)
}
fn put(
    rows: &mut [Residue],
    start: usize,
    index: usize,
    p: Position,
    evidence: AssignmentEvidence,
) -> Result<()> {
    let row = rows
        .get_mut(index.checked_sub(start).ok_or(Error::Binding)?)
        .ok_or(Error::Binding)?;
    if row.assignment.as_ref().is_some_and(|a| a.position != p) {
        return Err(Error::Binding);
    }
    row.assignment = Some(Assignment {
        position: p,
        evidence,
    });
    Ok(())
}
fn window(
    raw: &[Residue],
    rows: &mut [Residue],
    candidate: &d::Candidate,
    scheme: &NumberingScheme,
    kind: w::Window,
    bytes: &[u8],
) -> Result<Option<Issue>> {
    let start = rows.first().ok_or(Error::Binding)?.source_index;
    let required: &[u16] = match kind {
        w::Window::CysToTrp => &[23, 24, 25, 26, 39, 40, 41],
        w::Window::MiddleLoop => &[55, 66],
        w::Window::CysToJ => &[104],
    };
    let positions = required
        .iter()
        .map(|&n| reference_index(raw, n))
        .collect::<Vec<_>>();
    let missing = required
        .iter()
        .zip(&positions)
        .filter_map(|(&n, p)| p.is_none().then_some(n))
        .collect::<Vec<_>>();
    let issue = |reason| {
        Some(Issue {
            numbering_scheme: scheme.clone(),
            window: kind.clone(),
            issue: reason,
        })
    };
    if !missing.is_empty() {
        return Ok(issue(WindowIssue::MissingInvariantPositions {
            positions: missing,
        }));
    }
    let p = positions
        .into_iter()
        .map(Option::unwrap)
        .collect::<Vec<_>>();
    let (left, right) = match kind {
        w::Window::CysToTrp => {
            if p[0] + 1 != p[1]
                || p[1] + 1 != p[2]
                || p[2] + 1 != p[3]
                || p[4] + 1 != p[5]
                || p[5] + 1 != p[6]
            {
                return Ok(issue(WindowIssue::NonContiguousFlanks));
            }
            (p[0], p[6])
        }
        w::Window::MiddleLoop => (p[0], p[1]),
        w::Window::CysToJ => (p[0], candidate.cdr3.end),
    };
    let begin = left.checked_add(1).ok_or(Error::Binding)?;
    let length = right.checked_sub(begin).ok_or(Error::Binding)?;
    let body = bytes.get(left..=right).ok_or(Error::Binding)?;
    let heavy = candidate.locus == Locus::Heavy;
    let canonical = match kind {
        w::Window::CysToTrp => bytes[left] == b'C' && bytes[right] == b'W',
        w::Window::MiddleLoop => true,
        w::Window::CysToJ => bytes[left] == b'C' && bytes[right] == if heavy { b'W' } else { b'F' },
    };
    if !canonical || body.contains(&b'*') {
        return Ok(issue(WindowIssue::NonCanonicalAnchors));
    }
    let convention = match scheme {
        NumberingScheme::Kabat1991 => w::Convention::Kabat1991,
        NumberingScheme::Chothia1997 => w::Convention::Chothia1997,
        NumberingScheme::Imgt => return Err(Error::Binding),
    };
    let labels = match w::layout(&convention, &candidate.locus, &kind, length) {
        Ok(l) => l,
        Err(Error::Invalid(_) | Error::Limit(_)) => {
            return Ok(issue(WindowIssue::UnsupportedLength { residues: length }));
        }
        Err(e) => return Err(e),
    };
    for (i, label) in labels.residues.into_iter().enumerate() {
        put(
            rows,
            start,
            begin + i,
            label,
            AssignmentEvidence::PublishedWindowRule,
        )?;
    }
    Ok(None)
}
fn converted(
    raw: &[Residue],
    candidate: &d::Candidate,
    scheme: &NumberingScheme,
    bytes: &[u8],
) -> Result<(Vec<Residue>, Vec<Issue>)> {
    let mut rows = raw
        .iter()
        .map(|r| Residue {
            source_index: r.source_index,
            residue: r.residue,
            assignment: None,
            region: None,
        })
        .collect::<Vec<_>>();
    let start = rows.first().ok_or(Error::Binding)?.source_index;
    for r in raw {
        if let Some(a) = r
            .assignment
            .as_ref()
            .filter(|a| a.evidence == AssignmentEvidence::InvariantVReference)
        {
            if let Some(p) = fixed(&candidate.locus, &a.position) {
                put(
                    &mut rows,
                    start,
                    r.source_index,
                    p,
                    AssignmentEvidence::PublishedFrameworkCorrespondence,
                )?;
            }
        }
    }
    for (i, r) in candidate.fr4.iter().enumerate() {
        put(
            &mut rows,
            start,
            r.source_index,
            position(
                (if candidate.locus == Locus::Heavy {
                    103
                } else {
                    98
                }) + i as u16,
            ),
            AssignmentEvidence::JFr4Reference,
        )?;
    }
    let mut issues = vec![];
    for kind in [
        w::Window::CysToTrp,
        w::Window::MiddleLoop,
        w::Window::CysToJ,
    ] {
        if let Some(issue) = window(raw, &mut rows, candidate, scheme, kind, bytes)? {
            issues.push(issue);
        }
    }
    let mut labels = std::collections::BTreeSet::new();
    if rows
        .iter()
        .filter_map(|r| r.assignment.as_ref())
        .any(|a| !labels.insert(a.position.to_string()))
    {
        return Err(Error::Binding);
    }
    Ok((rows, issues))
}
fn kabat_cdrs(rows: &[Residue], issues: &[Issue], locus: &Locus) -> Result<Vec<CdrRegion>> {
    let heavy = *locus == Locus::Heavy;
    [
        (
            Region::Cdr1,
            w::Window::CysToTrp,
            if heavy { 31 } else { 24 },
            if heavy { 35 } else { 34 },
        ),
        (
            Region::Cdr2,
            w::Window::MiddleLoop,
            50,
            if heavy { 65 } else { 56 },
        ),
        (
            Region::Cdr3,
            w::Window::CysToJ,
            if heavy { 95 } else { 89 },
            if heavy { 102 } else { 97 },
        ),
    ]
    .into_iter()
    .map(|(kind, window, first, last)| {
        let mut members = vec![];
        for (i, r) in rows.iter().enumerate() {
            if r.assignment
                .as_ref()
                .map(|a| kabat_region(locus, &a.position))
                .transpose()?
                .as_ref()
                == Some(&kind)
            {
                members.push(i);
            }
        }
        let endpoints = [first, last].iter().all(|n| {
            rows.iter().any(|r| {
                r.assignment
                    .as_ref()
                    .is_some_and(|a| a.position == position(*n))
            })
        });
        let mut complete = endpoints && !issues.iter().any(|i| i.window == window);
        if let Some((&a, &b)) = members.first().zip(members.last()) {
            for row in &rows[a..=b] {
                complete &= row
                    .assignment
                    .as_ref()
                    .map(|a| kabat_region(locus, &a.position))
                    .transpose()?
                    .as_ref()
                    == Some(&kind);
            }
            Ok(CdrRegion {
                kind,
                status: if complete {
                    CdrStatus::Mapped {
                        start: rows[a].source_index,
                        end: rows[b].source_index + 1,
                    }
                } else {
                    CdrStatus::Unresolved
                },
            })
        } else {
            Ok(CdrRegion {
                kind,
                status: CdrStatus::Unresolved,
            })
        }
    })
    .collect()
}
fn imgt_cdrs(candidate: &d::Candidate) -> Vec<CdrRegion> {
    let status = |a: &d::LoopAssignment| match a {
        d::LoopAssignment::Numbered(n) => CdrStatus::Mapped {
            start: n.start,
            end: n.end,
        },
        d::LoopAssignment::Unresolved { .. } => CdrStatus::Unresolved,
    };
    vec![
        CdrRegion {
            kind: Region::Cdr1,
            status: status(&candidate.cdr1),
        },
        CdrRegion {
            kind: Region::Cdr2,
            status: status(&candidate.cdr2),
        },
        CdrRegion {
            kind: Region::Cdr3,
            status: CdrStatus::Mapped {
                start: candidate.cdr3.start,
                end: candidate.cdr3.end,
            },
        },
    ]
}
pub(crate) fn assemble(
    raw: &[Residue],
    candidate: &d::Candidate,
    selection: &Selection,
    bytes: &[u8],
    cancel: &mut impl FnMut() -> bool,
) -> Result<(Vec<Residue>, Vec<Issue>, Vec<CdrRegion>)> {
    if cancel() {
        return Err(Error::Cancelled);
    }
    let need_kabat = selection.numbering_scheme == NumberingScheme::Kabat1991
        || selection.region_definition == RegionDefinition::Kabat1991;
    let kabat = if need_kabat {
        Some(converted(
            raw,
            candidate,
            &NumberingScheme::Kabat1991,
            bytes,
        )?)
    } else {
        None
    };
    let chothia = if selection.numbering_scheme == NumberingScheme::Chothia1997 {
        Some(converted(
            raw,
            candidate,
            &NumberingScheme::Chothia1997,
            bytes,
        )?)
    } else {
        None
    };
    let mut rows = match selection.numbering_scheme {
        NumberingScheme::Imgt => raw.to_vec(),
        NumberingScheme::Kabat1991 => kabat.as_ref().ok_or(Error::Binding)?.0.clone(),
        NumberingScheme::Chothia1997 => chothia.as_ref().ok_or(Error::Binding)?.0.clone(),
    };
    for (i, row) in rows.iter_mut().enumerate() {
        if cancel() {
            return Err(Error::Cancelled);
        }
        row.region = match selection.region_definition {
            RegionDefinition::Imgt => raw[i].region.clone(),
            RegionDefinition::Kabat1991 => kabat.as_ref().ok_or(Error::Binding)?.0[i]
                .assignment
                .as_ref()
                .map(|a| kabat_region(&candidate.locus, &a.position))
                .transpose()?,
        };
    }
    let cdrs = match selection.region_definition {
        RegionDefinition::Imgt => imgt_cdrs(candidate),
        RegionDefinition::Kabat1991 => {
            let (rows, issues) = kabat.as_ref().ok_or(Error::Binding)?;
            kabat_cdrs(rows, issues, &candidate.locus)?
        }
    };
    let mut issues = vec![];
    if let Some((_, i)) = kabat {
        issues.extend(i);
    }
    if let Some((_, i)) = chothia {
        issues.extend(i);
    }
    Ok((rows, issues, cdrs))
}
