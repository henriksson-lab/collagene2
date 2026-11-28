//! Hypothetical CDS streams on the actual path; no annotation/edit ancestry.
use super::projection::Index;
use super::*;
use crate::{alignment::review, sequence::translation};
use collagene_model::{
    Error, Molecule, Sequence, Topology,
    annotations::biological_path::{BiologicalCall, BiologicalPath, PathLimits},
};
struct Coding {
    path: BiologicalPath,
    offset: usize,
    translation: translation::CdsTranslation,
    oriented: String,
}
impl Coding {
    fn load(
        design: &Molecule,
        revision: &str,
        selection: &CdsSelection,
        work: &mut Work,
    ) -> Result<Self, Error> {
        let assumption = &selection.design_cds;
        if assumption.source_revision != revision {
            return Err(Error::Unsupported(
                "CDS source pin differs from selected design event".into(),
            ));
        }
        let feature = design
            .features
            .get(assumption.feature_index)
            .ok_or_else(|| {
                Error::Unsupported("Declared design CDS feature is unavailable".into())
            })?;
        collagene_model::bounds::feature_preflight(&[std::slice::from_ref(feature)])?;
        if crate::construct_qc::annotation_revision(feature)? != assumption.annotation_revision {
            return Err(Error::Unsupported(
                "CDS annotation pin differs from selected design feature".into(),
            ));
        }
        let path = BiologicalPath::from_location(
            &feature.location,
            design.sequence.len(),
            PathLimits {
                max_calls: work.limits.max_coding_calls,
                max_spans: 4096,
            },
        )?;
        if path.len() > work.limits.max_coding_calls - work.summary.coding_calls {
            return Err(Error::TooLarge);
        }
        let mut effective = feature.clone();
        if !effective
            .qualifiers
            .iter()
            .any(|(key, _)| key == "transl_table")
        {
            effective.qualifiers.push((
                "transl_table".into(),
                Some(selection.genetic_code.to_string()),
            ));
        }
        let translated = translation::translate_cds_validated(design, &effective, revision)?;
        if translated.translation.options.genetic_code != selection.genetic_code {
            return Err(Error::Unsupported(
                "Selected genetic code conflicts with the CDS annotation table".into(),
            ));
        }
        if !work.coding(path.len()) {
            return Err(Error::TooLarge);
        }
        let oriented = feature.location.extract(&design.sequence)?;
        if oriented.len() != path.len() {
            return Err(Error::Unsupported(
                "Canonical annotation extraction/path lengths disagree".into(),
            ));
        }
        Ok(Self {
            path,
            offset: usize::from(translated.codon_start - 1),
            translation: translated,
            oriented,
        })
    }
}
#[derive(Clone, Copy)]
struct ProjectedCall {
    query: BiologicalCall,
    design: Option<usize>,
}
fn clipped(message: &str) -> String {
    let mut bytes = 0usize;
    message
        .chars()
        .filter(|c| !c.is_control())
        .take_while(|c| {
            bytes += c.len_utf8();
            bytes <= 2048
        })
        .collect()
}
impl CdsReport {
    pub(super) fn empty(selection: &CdsSelection, code: crate::construct_qc::CodePin) -> Self {
        Self {
            id: selection.id.clone(),
            design_cds: selection.design_cds.clone(),
            code,
            design_translation: None,
            projected_translation: None,
            mapping_complete: false,
            register_complete: false,
            consequences_complete: false,
            register_segments: vec![],
            consequences: vec![],
            findings: vec![],
        }
    }
    fn note(&mut self, rule: Rule, verdict: crate::construct_qc::Verdict, message: &str) {
        if self.findings.len() < 16 {
            self.findings.push(Finding {
                rule,
                verdict,
                message: clipped(message),
            });
        }
    }
    pub(super) fn unavailable(
        selection: &CdsSelection,
        code: crate::construct_qc::CodePin,
        message: &str,
    ) -> Self {
        let mut row = Self::empty(selection, code);
        row.note(
            Rule::CoordinateScope,
            crate::construct_qc::Verdict::Indeterminate,
            message,
        );
        row
    }
}
fn oriented_position(position: usize, report: &review::Report) -> usize {
    (position + report.reference_length - report.comparison.reference_origin)
        % report.reference_length
}
fn query_call(
    position: usize,
    reverse: bool,
    report: &review::Report,
) -> Result<BiologicalCall, Error> {
    let (start, end) = report.map_query(position, position + 1)?;
    if end != start + 1 {
        return Err(contract::bad(
            "Projected call lost its exact source interval",
        ));
    }
    Ok(BiologicalCall {
        position: start,
        reverse: reverse ^ report.comparison.query_reversed,
    })
}
fn base(molecule: &Molecule, call: BiologicalCall) -> u8 {
    let base = molecule.sequence.as_str().as_bytes()[call.position];
    if call.reverse {
        collagene_model::complement(base, molecule.sequence.alphabet())
    } else {
        base
    }
}

fn coordinates(calls: impl IntoIterator<Item = BiologicalCall>) -> Vec<CallSpan> {
    calls
        .into_iter()
        .map(|call| CallSpan {
            start: call.position,
            end: call.position + 1,
            reverse: call.reverse,
        })
        .collect()
}
fn add_call(
    calls: &mut Vec<ProjectedCall>,
    call: ProjectedCall,
    work: &mut Work,
) -> Result<(), Error> {
    if !work.coding(1) {
        return Err(Error::TooLarge);
    }
    calls.try_reserve(1).map_err(|_| Error::TooLarge)?;
    calls.push(call);
    Ok(())
}
fn add_insertions(
    calls: &mut Vec<ProjectedCall>,
    index: &Index,
    boundary: usize,
    reverse: bool,
    report: &review::Report,
    work: &mut Work,
) -> Result<(), Error> {
    if !work.mapping(1) {
        return Err(Error::TooLarge);
    }
    if let Some((start, length)) = index.insertion(boundary) {
        if reverse {
            for position in (start..start + length).rev() {
                add_call(
                    calls,
                    ProjectedCall {
                        query: query_call(position, reverse, report)?,
                        design: None,
                    },
                    work,
                )?;
            }
        } else {
            for position in start..start + length {
                add_call(
                    calls,
                    ProjectedCall {
                        query: query_call(position, reverse, report)?,
                        design: None,
                    },
                    work,
                )?;
            }
        }
    }
    Ok(())
}
fn project(
    coding: &Coding,
    design: &Molecule,
    query: &Molecule,
    report: &review::Report,
    index: &Index,
    axis: &ConsensusAxis,
    work: &mut Work,
) -> Result<Vec<ProjectedCall>, Error> {
    let circular = query.topology == Topology::Circular
        || matches!(axis, ConsensusAxis::DeclaredSingleCircularTurn { .. });
    let mut output = Vec::new();
    let first = coding.path.call_at(0).ok_or(Error::EmptySequence)?;
    let last = coding.path.call_at(coding.path.len() - 1).unwrap();
    let first_position = oriented_position(first.position, report);
    let last_position = oriented_position(last.position, report);
    // Boundary insertions cannot be silently assigned inside or outside a gene.
    for boundary in [
        if first.reverse {
            first_position + 1
        } else {
            first_position
        },
        if last.reverse {
            last_position
        } else {
            last_position + 1
        },
    ] {
        if !work.mapping(1) {
            return Err(Error::TooLarge);
        }
        if index.insertion(boundary).is_some() {
            return Err(Error::Unsupported(
                "Insertion touches a declared CDS boundary; coding membership is unresolved".into(),
            ));
        }
    }
    if index.at(first_position).is_none() || index.at(last_position).is_none() {
        return Err(Error::Unsupported("Selected path deletes a CDS boundary anchor; a replacement gene boundary is not inferred".into()));
    }
    for ordinal in 0..coding.path.len() {
        if !work.mapping(1) {
            return Err(Error::TooLarge);
        }
        let call = coding.path.call_at(ordinal).unwrap();
        let position = oriented_position(call.position, report);
        if let Some(query_position) = index.at(position) {
            add_call(
                &mut output,
                ProjectedCall {
                    query: query_call(query_position, call.reverse, report)?,
                    design: Some(ordinal),
                },
                work,
            )?;
        }
        let Some(next) = coding.path.call_at(ordinal + 1) else {
            continue;
        };
        let next_position = oriented_position(next.position, report);
        // Even a joined path that skips an intron may turn across the chosen
        // coordinate origin. A linear supplied axis does not imply that turn.
        let turns = if call.reverse {
            next_position > position
        } else {
            next_position < position
        };
        if turns && !circular {
            return Err(Error::Unsupported("Selected CDS order turns across the linear consensus coordinate origin; declare a single circular-coordinate turn explicitly".into()));
        }
        let adjacent = call.next(design.sequence.len(), design.topology) == Some(next);
        if !adjacent {
            for boundary in [
                if call.reverse { position } else { position + 1 },
                if next.reverse {
                    next_position + 1
                } else {
                    next_position
                },
            ] {
                if !work.mapping(1) {
                    return Err(Error::TooLarge);
                }
                if index.insertion(boundary).is_some() {
                    return Err(Error::Unsupported(
                        "Insertion at a joined CDS span boundary has unresolved coding membership"
                            .into(),
                    ));
                }
            }
            continue;
        }
        let crossing = if call.reverse {
            position == 0 && next_position + 1 == report.reference_length
        } else {
            position + 1 == report.reference_length && next_position == 0
        };
        if crossing && !circular {
            return Err(Error::Unsupported("Selected CDS requires adjacency across the linear consensus ends; declare a single circular-coordinate turn explicitly".into()));
        }
        let boundaries = if crossing {
            if call.reverse {
                vec![0, report.reference_length]
            } else {
                vec![report.reference_length, 0]
            }
        } else {
            vec![if call.reverse {
                position
            } else {
                next_position
            }]
        };
        for boundary in boundaries {
            if index.insertion(boundary).is_some()
                && (index.at(position).is_none() || index.at(next_position).is_none())
            {
                return Err(Error::Unsupported(
                    "Inserted coding calls lack two paired interior anchors on this selected path"
                        .into(),
                ));
            }
            add_insertions(&mut output, index, boundary, call.reverse, report, work)?;
        }
    }
    Ok(output)
}
fn register(row: &mut CdsReport, calls: &[ProjectedCall], offset: usize, work: &mut Work) {
    row.register_complete = true;
    for (projected, call) in calls.iter().enumerate().skip(offset) {
        if !work.mapping(1) {
            row.register_complete = false;
            break;
        }
        let Some(design) = call.design.filter(|ordinal| *ordinal >= offset) else {
            continue;
        };
        let (design, projected) = (design - offset, projected - offset);
        let delta = ((projected as i64 - design as i64).rem_euclid(3)) as u8;
        if let Some(last) = row.register_segments.last_mut() {
            if last.design_end == design
                && last.projected_end == projected
                && last.offset_mod3 == delta
            {
                last.design_end += 1;
                last.projected_end += 1;
                continue;
            }
        }
        if work.summary.segments == work.limits.max_segments {
            row.register_complete = false;
            break;
        }
        work.summary.segments += 1;
        row.register_segments.push(RegisterSegment {
            design_start: design,
            design_end: design + 1,
            projected_start: projected,
            projected_end: projected + 1,
            offset_mod3: delta,
        });
    }
    if !row.register_complete {
        row.note(
            Rule::Register,
            crate::construct_qc::Verdict::Indeterminate,
            "Conditional register entry budget reached; omitted segments are not absent",
        );
    }
}
fn consequences(
    row: &mut CdsReport,
    coding: &Coding,
    calls: &[ProjectedCall],
    query: &Molecule,
    work: &mut Work,
) {
    let projected = row
        .projected_translation
        .as_ref()
        .expect("native projected translation");
    let mut entries = Vec::new();
    let mut complete = true;
    for (index, amino) in projected.protein.bytes().enumerate() {
        if !work.mapping(1) {
            complete = false;
            break;
        }
        let start = coding.offset + index * 3;
        let triple = &calls[start..start + 3];
        let ordinals = triple
            .iter()
            .map(|call| call.design)
            .collect::<Option<Vec<_>>>();
        let anchor = ordinals
            .as_ref()
            .filter(|ordinals| {
                ordinals[0] >= coding.offset
                    && (ordinals[0] - coding.offset) % 3 == 0
                    && ordinals[1] == ordinals[0] + 1
                    && ordinals[2] == ordinals[0] + 2
            })
            .map(|ordinals| (ordinals[0] - coding.offset) / 3);
        let design_codon = anchor.map(|anchor| {
            coding.oriented[coding.offset + anchor * 3..coding.offset + anchor * 3 + 3].to_owned()
        });
        let design_amino = anchor
            .and_then(|anchor| {
                coding
                    .translation
                    .translation
                    .protein
                    .as_bytes()
                    .get(anchor)
                    .copied()
            })
            .map(char::from);
        let projected_codon =
            String::from_utf8(triple.iter().map(|call| base(query, call.query)).collect())
                .expect("validated ASCII nucleotide calls");
        let internal = amino == b'*' && index + 1 < projected.translated_codons;
        let changed = design_codon.as_deref() != Some(projected_codon.as_str());
        if !changed && !internal {
            continue;
        }
        if work.summary.consequences == work.limits.max_consequences {
            complete = false;
            break;
        }
        let unresolved = amino == b'X' || design_amino == Some('X');
        let kind = if internal {
            ConsequenceKind::DefiniteInternalStop
        } else if unresolved {
            ConsequenceKind::Unresolved
        } else if design_amino == Some('*') && amino != b'*' {
            ConsequenceKind::StopLost
        } else if anchor.is_none() {
            if triple.iter().all(|call| call.design.is_none()) {
                ConsequenceKind::InsertedCodon
            } else {
                ConsequenceKind::Unresolved
            }
        } else if design_amino == Some(char::from(amino)) {
            ConsequenceKind::Synonymous
        } else {
            ConsequenceKind::Missense
        };
        let design_coordinates = coordinates(
            triple
                .iter()
                .filter_map(|call| call.design.and_then(|ordinal| coding.path.call_at(ordinal))),
        );
        let message = if anchor.is_some() {
            "Codon observation is conditional on this selected alignment path, annotation and genetic code; correspondence uniqueness and native edit cause are not established."
        } else {
            "Projected codon does not anchor one complete design-frame codon. Any definite stop is an absolute conditional stream observation, not proof of a newly introduced mutation or native edit cause."
        };
        entries.push(CodonConsequence {
            projected_codon_index: index,
            projected_codon,
            projected_amino_acid: char::from(amino),
            design_codon_index: anchor,
            design_codon,
            design_amino_acid: design_amino,
            kind,
            design_coordinates,
            consensus_coordinates: coordinates(triple.iter().map(|call| call.query)),
            message: message.into(),
        });
        work.summary.consequences += 1;
    }
    row.consequences = entries;
    row.consequences_complete = complete;
    if !complete {
        row.note(Rule::CodonConsequence,crate::construct_qc::Verdict::Indeterminate,"Conditional changed/uncertain codon entry budget reached; unreported consequences are not absent");
    }
}
pub(super) fn assess(
    selection: &CdsSelection,
    code: crate::construct_qc::CodePin,
    request: &Request,
    design: &Molecule,
    query: &Molecule,
    report: &review::Report,
    index: &Index,
    work: &mut Work,
) -> CdsReport {
    use crate::construct_qc::Verdict;
    let mut row = CdsReport::empty(selection, code);
    row.note(Rule::Scope,Verdict::Informational,"Every coding observation is conditional on the actual selected alignment path. Certified affine score does not establish unique correspondence, annotation truth, read confidence or causal edits.");
    let coding = match Coding::load(design, &request.design.revision, selection, work) {
        Ok(coding) => coding,
        Err(error) => {
            row.note(
                Rule::AnnotationBasis,
                Verdict::Indeterminate,
                &format!("Current CDS admission unavailable: {error}"),
            );
            return row;
        }
    };
    row.design_translation = Some(coding.translation.translation.clone());
    if !coding.translation.warnings.is_empty() {
        row.note(
            Rule::TerminalContext,
            Verdict::Informational,
            &coding.translation.warnings.join(" "),
        );
    }
    let calls = match project(
        &coding,
        design,
        query,
        report,
        index,
        &request.consensus_axis,
        work,
    ) {
        Ok(calls) => calls,
        Err(error) => {
            row.note(
                Rule::CoordinateScope,
                Verdict::Indeterminate,
                &format!("Conditional coding projection unavailable: {error}"),
            );
            return row;
        }
    };
    if coding.offset > 0
        && (calls.len() <= coding.offset
            || calls
                .iter()
                .take(coding.offset + 1)
                .enumerate()
                .any(|(ordinal, call)| call.design != Some(ordinal)))
    {
        row.note(
            Rule::CoordinateScope,
            Verdict::Indeterminate,
            "Changed/deleted partial 5′ prefix cannot silently redefine the declared codon_start",
        );
        return row;
    }
    let sequence = match Sequence::dna(
        &String::from_utf8(calls.iter().map(|call| base(query, call.query)).collect())
            .expect("ASCII calls"),
    ) {
        Ok(sequence) => sequence,
        Err(error) => {
            row.note(
                Rule::CoordinateScope,
                Verdict::Indeterminate,
                &format!("Projected coding stream unavailable: {error}"),
            );
            return row;
        }
    };
    match translation::translate(&sequence, &coding.translation.translation.options) {
        Ok(translated) => row.projected_translation = Some(translated),
        Err(error) => {
            row.note(
                Rule::CodonConsequence,
                Verdict::Indeterminate,
                &format!("Native projected translation unavailable: {error}"),
            );
            return row;
        }
    }
    row.mapping_complete = true;
    register(&mut row, &calls, coding.offset, work);
    consequences(&mut row, &coding, &calls, query, work);
    if row
        .projected_translation
        .as_ref()
        .is_some_and(|translation| translation.trailing_bases != 0)
    {
        row.note(Rule::Register,Verdict::Informational,"Projected coding stream has trailing bases under the declared frame; no complete terminal codon is inferred");
    }
    row
}
