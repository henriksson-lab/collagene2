//! Declaration-relative geometry and literal pinned identity. These facts never
//! infer an assembly event, annotation function or a changed frame in an ancestor.
use super::*;

fn single(l: &Location, reverse: bool) -> Option<(usize, usize, bool)> {
    match l {
        Location::Range {
            start,
            end,
            before: false,
            after: false,
        } => Some((*start, *end, reverse)),
        Location::Complement(inner) => single(inner, !reverse),
        _ => None,
    }
}
fn phase_at(checks: &Checks<'_>, cds: usize, boundary: usize) -> Option<(usize, char)> {
    let (start, end, reverse) = single(&checks.m.features[cds].location, false)?;
    if !(start..=end).contains(&boundary) {
        return None;
    }
    let frame = checks.frames.get(&cds)?;
    let offset = if reverse {
        end - boundary
    } else {
        boundary - start
    };
    let offset = offset.checked_sub(frame.offset)?;
    // At a codon boundary no single codon crosses the boundary. Internal-stop
    // status of the full CDS is a separate existing translation-owner finding.
    let crossing = if offset % 3 == 0 {
        '-'
    } else {
        *frame.protein.as_bytes().get(offset / 3)? as char
    };
    Some((offset % 3, crossing))
}
fn register(phase: usize, expected: u8, crossing: char) -> (Verdict, String) {
    if phase != usize::from(expected) {
        (
            Verdict::Fail,
            format!(
                "Observed coding offset modulo 3 is {phase}; the explicit declaration requires {expected}."
            ),
        )
    } else if crossing == '*' {
        (
            Verdict::Fail,
            "The boundary crosses a definite stop codon in the existing CDS translation.".into(),
        )
    } else if crossing == 'X' {
        (Verdict::Indeterminate,"The boundary crosses an unresolved codon; an uninterrupted coding junction is not established.".into())
    } else {
        (
            Verdict::Pass,
            format!(
                "Observed coding offset modulo 3 is {phase}, matching the declared register; no definite or unresolved codon crosses this boundary."
            ),
        )
    }
}
pub(super) fn inspect(checks: &mut Checks<'_>) -> Result<(), Error> {
    let m = checks.m;
    let o = checks.o;
    for claim in &o.reference_claims {
        let i = claim.feature_index;
        if !checks.assumed(i) {
            checks.push(
                Rule::ReferenceIntegrity,
                Verdict::Indeterminate,
                &[i],
                "Exact element identity depends on a current declared feature-role/coordinate assumption; the annotation is unreviewed or undeclared.",
            )?;
            continue;
        }
        let Some(reference) = o.references.iter().find(|r| r.id == claim.reference_id) else {
            checks.push(
                Rule::ReferenceIntegrity,
                Verdict::Indeterminate,
                &[i],
                "No pinned reference is available for the declared element; annotation alone cannot establish identity.",
            )?;
            continue;
        };
        if reference.sequence_revision != reference_revision(&reference.sequence)? {
            checks.push(
                Rule::ReferenceIntegrity,
                Verdict::Indeterminate,
                &[i],
                "Declared reference identity does not match its captured sequence; no integrity comparison is admitted.",
            )?;
            continue;
        }
        let extracted = match m.features[i].location.extract(&m.sequence) {
            Ok(s) => s,
            Err(e) => {
                checks.push(
                    Rule::ReferenceIntegrity,
                    Verdict::Indeterminate,
                    &[i],
                    &format!(
                        "Feature does not define an exact extractable reference interval: {e}"
                    ),
                )?;
                continue;
            }
        };
        if reference.sequence.alphabet() != m.sequence.alphabet()
            || !extracted
                .bytes()
                .chain(reference.sequence.as_str().bytes())
                .all(|b| crate::sequence::analysis::base_mask(b).is_power_of_two())
        {
            checks.push(
                Rule::ReferenceIntegrity,
                Verdict::Indeterminate,
                &[i],
                "Reference comparison requires the same nucleotide alphabet and fully concrete calls; IUPAC ambiguity is not silently resolved.",
            )?;
        } else {
            let same = extracted == reference.sequence.as_str();
            checks.push(
                Rule::ReferenceIntegrity,
                if same {Verdict::Pass}else{Verdict::Fail},
                &[i],
                &format!("Literal oriented feature sequence {} pinned reference {} ({}). This checks sequence identity only, not biological role or provenance authenticity.",if same {"matches"}else{"differs from"},reference.id,reference.sequence_revision),
            )?;
        }
    }
    for claim in &o.fusions {
        let indices = [claim.cds_feature, claim.part_feature];
        let geometry = single(&m.features[claim.cds_feature].location, false)
            .zip(single(&m.features[claim.part_feature].location, false));
        if !checks.assumed(claim.cds_feature) || !checks.assumed(claim.part_feature) {
            checks.push(
                Rule::FusionRegister,
                Verdict::Indeterminate,
                &indices,
                "Fusion roles require explicit current annotation assumptions for both CDS and component.",
            )?;
            continue;
        }
        let Some(((cs, ce, cr), (ps, pe, pr))) = geometry else {
            checks.push(
                Rule::FusionRegister,
                Verdict::Indeterminate,
                &indices,
                "Compound or unresolved feature paths need an explicit biological-coordinate mapping; this foundation does not infer fusion register across them.",
            )?;
            continue;
        };
        if cr != pr || ps < cs || pe > ce {
            checks.push(
                Rule::FusionRegister,
                Verdict::Indeterminate,
                &indices,
                "Component is not contained in the declared CDS on the same strand; adjacency alone does not establish a shared coding path.",
            )?;
            continue;
        }
        let phases =
            phase_at(checks, claim.cds_feature, ps).zip(phase_at(checks, claim.cds_feature, pe));
        let Some((a, b)) = phases else {
            checks.push(
                Rule::FusionRegister,
                Verdict::Indeterminate,
                &indices,
                "No resolved complete CDS frame is available for this component.",
            )?;
            continue;
        };
        let (va, ra) = register(a.0, claim.expected_offset_mod3, a.1);
        let (vb, rb) = register(b.0, claim.expected_offset_mod3, b.1);
        let verdict = if va == Verdict::Fail || vb == Verdict::Fail {
            Verdict::Fail
        } else if va == Verdict::Indeterminate || vb == Verdict::Indeterminate {
            Verdict::Indeterminate
        } else {
            Verdict::Pass
        };
        checks.push(
            Rule::FusionRegister,
            verdict,
            &indices,
            &format!("Both component boundaries must share the explicitly declared phase. Start: {ra} End: {rb} This does not establish fusion function, absence of upstream stops or experimental assembly."),
        )?;
    }
    for claim in &o.junctions {
        let i = claim.cds_feature;
        if !checks.assumed(i) {
            checks.push(
                Rule::DeclaredJunctionRegister,
                Verdict::Indeterminate,
                &[i],
                "Declared junction has no current CDS-role/coordinate assumption.",
            )?;
            continue;
        }
        let Some((phase, crossing)) = phase_at(checks, i, claim.boundary) else {
            checks.push(
                Rule::DeclaredJunctionRegister,
                Verdict::Indeterminate,
                &[i],
                "Declared boundary is outside a resolved single-path complete CDS frame; no register is inferred.",
            )?;
            continue;
        };
        let (verdict, reason) = register(phase, claim.expected_offset_mod3, crossing);
        checks.push(
            Rule::DeclaredJunctionRegister,
            verdict,
            &[i],
            &format!("Declared boundary {}: {reason} This current-sequence fact does not show that an assembly/edit introduced a frame change; history identifiers are documentary only. Fresh causal junction review remains owned by assembly::junctions::inspect_plan.",claim.boundary),
        )?;
    }
    Ok(())
}
