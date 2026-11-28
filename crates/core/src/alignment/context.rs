//! Read-only source annotation/CDS context on an already computed alignment.
//! Codons are translated in complete source-feature order, never per page or gap.
use crate::{alignment::review::Report, sequence::translation};
use collagene_model::{Error, Feature, Location, Molecule};
use serde::{Deserialize, Serialize};
use std::io;

pub const PAGE_COLUMNS: usize = 600;
const MAX_CDS_POSITIONS: usize = 1_500_000;
const MAX_TRANSLATION_WORK: usize = 64_000_000;
const CONTRACT: &str = "pairwise-aligned-source-context-2";
const MAX_BYTES: usize = 32 * 1024 * 1024;
fn bad(message: &str) -> Error {
    Error::Format(format!("Pairwise context: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
struct Counter(usize);
impl io::Write for Counter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| io::Error::other("pairwise context byte limit"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Borrowed counts and byte limits before hashing or copying worker inputs.
/// Actual source correspondence is checked by prepare inside the worker.
pub fn preflight(
    report: &Report,
    reference: &Molecule,
    query: &Molecule,
    _translate: bool,
) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[reference, query])?;
    let a = &report.comparison.alignment;
    if a.reference_aligned.len() > crate::alignment::MAX_ALIGNMENT_COLUMNS
        || a.query_aligned.len() > crate::alignment::MAX_ALIGNMENT_COLUMNS
        || a.operations.len() > crate::alignment::MAX_ALIGNMENT_COLUMNS
        || a.differences.len() > crate::alignment::MAX_ALIGNMENT_COLUMNS
        || report.comparison.warnings.len() > 4096
        || a.warnings.len() > 4096
        || report
            .comparison
            .warnings
            .iter()
            .chain(&a.warnings)
            .any(|s| s.len() > 65_536)
    {
        return Err(Error::TooLarge);
    }
    serde_json::to_writer(&mut Counter(0), &(report, reference, query))
        .map_err(|_| Error::TooLarge)?;
    report.validate_binding(&report.context_revision)?;
    if reference.sequence.len() != report.reference_length
        || query.sequence.len() != report.query_length
    {
        return Err(bad("source lengths disagree with comparison"));
    }
    Ok(())
}
pub fn revision(
    report: &Report,
    reference: &Molecule,
    query: &Molecule,
    translate: bool,
) -> Result<String, Error> {
    preflight(report, reference, query, translate)?;
    crate::encoding::bounded_json::hash(&(CONTRACT, report, reference, query, translate))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Side {
    Reference,
    Query,
}
impl Side {
    pub fn label(self) -> &'static str {
        match self {
            Self::Reference => "R",
            Self::Query => "Q",
        }
    }
}
#[derive(Debug)]
struct Codon {
    positions: [usize; 3],
    residue: u8,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
struct Track {
    side: Side,
    index: usize,
    label: String,
    kind: String,
    spans: Vec<(usize, usize, bool)>,
    translation_calls: Vec<(usize, u8)>,
    notes: Vec<String>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Context {
    pub context_revision: String,
    comparison_revision: String,
    reference_length: usize,
    query_length: usize,
    alignment_columns: usize,
    pub translated: bool,
    prefixes: Vec<(usize, usize)>,
    tracks: Vec<Track>,
    pub warnings: Vec<String>,
}
#[derive(Debug, PartialEq, Eq)]
pub struct PageTrack {
    pub side: Side,
    pub feature_index: usize,
    pub label: String,
    pub kind: String,
    /// One symbol per alignment column: >/< source direction, space outside/gap.
    pub annotation: String,
    /// One symbol per column: residue at codon middle, '.' at other source calls.
    pub translation: Option<String>,
    pub notes: Vec<String>,
}
#[derive(Debug, PartialEq, Eq)]
pub struct Page {
    pub start: usize,
    pub end: usize,
    pub reference: String,
    pub query: String,
    /// Original zero-based source positions; None is an alignment gap.
    pub reference_positions: Vec<Option<usize>>,
    pub query_positions: Vec<Option<usize>>,
    pub tracks: Vec<PageTrack>,
}
fn label(feature: &Feature) -> String {
    let text = feature.label();
    let mut result: String = text.chars().take(160).collect();
    if text.chars().count() > 160 {
        result.push_str("… (abbreviated)");
    }
    result
}
fn exact_length(location: &Location) -> Result<usize, String> {
    match location {
        Location::Range {
            start,
            end,
            before: false,
            after: false,
        } => end
            .checked_sub(*start)
            .ok_or_else(|| "Invalid source range.".into()),
        Location::Complement(inner) => exact_length(inner),
        Location::Join(parts) => parts.iter().try_fold(0usize, |sum, part| {
            sum.checked_add(exact_length(part)?)
                .filter(|n| *n <= MAX_CDS_POSITIONS)
                .ok_or_else(|| "CDS mapping exceeds the 1,500,000-position context budget.".into())
        }),
        _ => Err("Fuzzy, ordered or unresolved CDS location has no exact codon mapping.".into()),
    }
}
fn append_positions(location: &Location, reverse: bool, positions: &mut Vec<usize>) {
    match location {
        Location::Range { start, end, .. } => {
            if reverse {
                positions.extend((*start..*end).rev());
            } else {
                positions.extend(*start..*end);
            }
        }
        Location::Complement(inner) => append_positions(inner, !reverse, positions),
        Location::Join(parts) => {
            if reverse {
                for part in parts.iter().rev() {
                    append_positions(part, true, positions);
                }
            } else {
                for part in parts {
                    append_positions(part, false, positions);
                }
            }
        }
        _ => unreachable!("exact_length checked before position extraction"),
    }
}
fn imprecise_location(location: &Location) -> bool {
    match location {
        Location::Range { before, after, .. } => *before || *after,
        Location::Complement(inner) => imprecise_location(inner),
        Location::Join(parts) => parts.iter().any(imprecise_location),
        Location::Order(_) | Location::Unresolved(_) => true,
    }
}
fn add_tracks(
    molecule: &Molecule,
    side: Side,
    source_revision: &str,
    with_translation: bool,
    position_budget: &mut usize,
    work_budget: &mut usize,
    tracks: &mut Vec<Track>,
    warnings: &mut Vec<String>,
) {
    for (index, feature) in molecule.features.iter().enumerate() {
        let mut track = Track {
            side,
            index,
            label: label(feature),
            kind: feature.kind.clone(),
            spans: feature.location.spans(),
            translation_calls: vec![],
            notes: vec![],
        };
        if feature.needs_review {
            track
                .notes
                .push("Source annotation requires review.".into());
        }
        if imprecise_location(&feature.location) {
            track.notes.push("Fuzzy, ordered or unresolved source location: arrows show only retained coordinate spans, not exact biological endpoints/order.".into());
        }
        if track.spans.is_empty() {
            warnings.push(format!("{} feature {} ({}) has no mapped source spans; inspect its retained location in the source document.",side.label(),index+1,track.label));
        }
        if with_translation && feature.kind.eq_ignore_ascii_case("CDS") {
            let translated = (|| -> Result<Vec<Codon>, String> {
                let length = exact_length(&feature.location)?;
                let positions_next = position_budget
                    .checked_add(length)
                    .filter(|n| *n <= MAX_CDS_POSITIONS)
                    .ok_or("Aggregate CDS mapping exceeds 1,500,000 source positions.")?;
                // The source was validated once. Reuse the native validated CDS
                // path instead of hashing a complete source for every feature.
                let work_next = work_budget.checked_add(length)
                    .filter(|n| *n <= MAX_TRANSLATION_WORK)
                    .ok_or("Aggregate CDS translation exceeds the 64,000,000-base context work budget.")?;
                *position_budget = positions_next;
                *work_budget = work_next;
                let cds = translation::translate_cds_validated(molecule, feature, source_revision)
                    .map_err(|e| e.to_string())?;
                let mut positions = Vec::with_capacity(length);
                append_positions(&feature.location, false, &mut positions);
                track.notes.push(format!("Native code {} ({}) · codon_start {}. Full-source biological CDS order; residues are placed at each codon's middle source base. Terminal '*' remains shown as codon context even though retained CDS protein omits a terminal stop.",cds.translation.options.genetic_code,cds.translation.code_name,cds.codon_start));
                track.notes.extend(cds.warnings);
                let frame = usize::from(cds.codon_start - 1);
                Ok(cds
                    .translation
                    .protein
                    .bytes()
                    .enumerate()
                    .map(|(i, residue)| {
                        let first = frame + i * 3;
                        Codon {
                            positions: [
                                positions[first],
                                positions[first + 1],
                                positions[first + 2],
                            ],
                            residue,
                        }
                    })
                    .collect())
            })();
            match translated {
                Ok(codons) => {
                    let mut calls: Vec<_> = codons
                        .iter()
                        .flat_map(|codon| {
                            codon
                                .positions
                                .iter()
                                .enumerate()
                                .map(move |(offset, position)| {
                                    (*position, if offset == 1 { codon.residue } else { b'.' })
                                })
                        })
                        .collect();
                    calls.sort_unstable_by_key(|(position, _)| *position);
                    for (position, glyph) in calls {
                        if let Some((previous, value)) = track.translation_calls.last_mut() {
                            if *previous == position {
                                mark(value, glyph);
                                continue;
                            }
                        }
                        track.translation_calls.push((position, glyph));
                    }
                }
                Err(error) => {
                    let note = format!("CDS translation withheld: {error}");
                    warnings.push(format!(
                        "{} feature {} ({}): {note}",
                        side.label(),
                        index + 1,
                        track.label
                    ));
                    track.notes.push(note);
                }
            }
        }
        tracks.push(track);
    }
}

/// Source validation and translation happen once for a prepared context. Native
/// scientific input bounds apply to these same compared models, not unrelated docs.
pub fn prepare(
    report: &Report,
    reference: &Molecule,
    query: &Molecule,
    with_translation: bool,
) -> Result<Context, Error> {
    let context_revision = revision(report, reference, query, with_translation)?;
    report.validate_sources(reference, query)?;
    let a = &report.comparison.alignment;
    let mut prefixes = vec![];
    let (mut r, mut q) = (a.reference_start, a.query_start);
    for (column, (rb, qb)) in a
        .reference_aligned
        .bytes()
        .zip(a.query_aligned.bytes())
        .enumerate()
    {
        if column % PAGE_COLUMNS == 0 {
            prefixes.push((r, q));
        }
        r += usize::from(rb != b'-');
        q += usize::from(qb != b'-');
    }
    if prefixes.is_empty() {
        prefixes.push((r, q));
    }
    let mut context = Context {
        context_revision,
        comparison_revision: report.context_revision.clone(),
        reference_length: reference.sequence.len(),
        query_length: query.sequence.len(),
        alignment_columns: a.reference_aligned.len(),
        translated: with_translation,
        prefixes,
        tracks: vec![],
        warnings: vec![],
    };
    let (mut positions, mut work) = (0, reference.sequence.len() + query.sequence.len());
    add_tracks(
        reference,
        Side::Reference,
        &report.reference_revision,
        with_translation,
        &mut positions,
        &mut work,
        &mut context.tracks,
        &mut context.warnings,
    );
    add_tracks(
        query,
        Side::Query,
        &report.query_revision,
        with_translation,
        &mut positions,
        &mut work,
        &mut context.tracks,
        &mut context.warnings,
    );
    context.validate()?;
    Ok(context)
}
fn mark(existing: &mut u8, value: u8) {
    if *existing == b' ' || *existing == b'.' {
        *existing = value;
    } else if *existing != value && value != b'.' {
        *existing = b'!';
    }
}
impl Context {
    /// Standalone bounded receipt validation, not a rerun or proof of translation.
    pub fn validate(&self) -> Result<(), Error> {
        if !digest(&self.context_revision)
            || !digest(&self.comparison_revision)
            || self.reference_length == 0
            || self.query_length == 0
            || self.reference_length > collagene_model::MAX_SEQUENCE_LEN
            || self.query_length > collagene_model::MAX_SEQUENCE_LEN
            || self.alignment_columns > crate::alignment::MAX_ALIGNMENT_COLUMNS
            || self.prefixes.len() != self.alignment_columns.div_ceil(PAGE_COLUMNS).max(1)
            || self.tracks.len() > 4096
            || self.warnings.len() > 4096
        {
            return Err(bad("invalid receipt identity or bounds"));
        }
        let (mut previous_r, mut previous_q) = (0, 0);
        for &(r, q) in &self.prefixes {
            if r > self.reference_length
                || q > self.query_length
                || r < previous_r
                || q < previous_q
            {
                return Err(bad("invalid cached page coordinates"));
            }
            (previous_r, previous_q) = (r, q);
        }
        let (mut spans, mut calls, mut text) = (0usize, 0usize, 0usize);
        let (mut previous_reference, mut previous_query) = (None, None);
        for track in &self.tracks {
            if track.index >= 4096
                || track.label.len() > 800
                || track.kind.len() > 256
                || track.notes.len() > 64
            {
                return Err(bad("invalid source feature context"));
            }
            let n = match track.side {
                Side::Reference => self.reference_length,
                Side::Query => self.query_length,
            };
            let previous_index = match track.side {
                Side::Reference => &mut previous_reference,
                Side::Query => &mut previous_query,
            };
            if previous_index.is_some_and(|index| index >= track.index) {
                return Err(bad("source feature indices are not unique and ordered"));
            }
            *previous_index = Some(track.index);
            spans = spans
                .checked_add(track.spans.len())
                .filter(|n| *n <= 16384)
                .ok_or(Error::TooLarge)?;
            calls = calls
                .checked_add(track.translation_calls.len())
                .filter(|n| *n <= MAX_CDS_POSITIONS)
                .ok_or(Error::TooLarge)?;
            if (!self.translated || !track.kind.eq_ignore_ascii_case("CDS"))
                && !track.translation_calls.is_empty()
            {
                return Err(bad("unexpected translation calls"));
            }
            if track.spans.iter().any(|(a, b, _)| a > b || *b > n) {
                return Err(bad("invalid feature coordinate"));
            }
            let mut ordered: Vec<_> = track.spans.iter().map(|(a, b, _)| (*a, *b)).collect();
            ordered.sort_unstable();
            let mut coverage = Vec::<(usize, usize)>::new();
            for (start, end) in ordered {
                if start == end {
                    continue;
                }
                if let Some((_, last_end)) = coverage.last_mut() {
                    if start <= *last_end {
                        *last_end = (*last_end).max(end);
                        continue;
                    }
                }
                coverage.push((start, end));
            }
            let mut part = 0usize;
            let mut previous = None;
            for &(position, glyph) in &track.translation_calls {
                while coverage.get(part).is_some_and(|(_, end)| *end <= position) {
                    part += 1;
                }
                if position >= n
                    || previous.is_some_and(|p| p >= position)
                    || !(glyph.is_ascii_uppercase() || b".*!".contains(&glyph))
                    || !coverage
                        .get(part)
                        .is_some_and(|(start, end)| *start <= position && position < *end)
                {
                    return Err(bad("invalid cached CDS source call"));
                }
                previous = Some(position);
            }
            for value in std::iter::once(&track.label)
                .chain(std::iter::once(&track.kind))
                .chain(&track.notes)
            {
                if value.len() > 2048 {
                    return Err(Error::TooLarge);
                }
                text = text
                    .checked_add(value.len())
                    .filter(|n| *n <= 8 * 1024 * 1024)
                    .ok_or(Error::TooLarge)?;
            }
        }
        for value in &self.warnings {
            if value.len() > 2048 {
                return Err(Error::TooLarge);
            }
            text = text
                .checked_add(value.len())
                .filter(|n| *n <= 8 * 1024 * 1024)
                .ok_or(Error::TooLarge)?;
        }
        serde_json::to_writer(&mut Counter(0), self).map_err(|_| Error::TooLarge)?;
        Ok(())
    }
    pub fn page(&self, report: &Report, page: usize) -> Result<Page, String> {
        if report.context_revision != self.comparison_revision
            || report.reference_length != self.reference_length
            || report.query_length != self.query_length
            || report.comparison.alignment.reference_aligned.len() != self.alignment_columns
            || report.comparison.alignment.query_aligned.len() != self.alignment_columns
        {
            return Err(
                "The comparison context changed; prepare its annotation context again.".into(),
            );
        }
        let (mut r, mut q) = *self
            .prefixes
            .get(page)
            .ok_or("Alignment context page is unavailable.")?;
        let a = &report.comparison.alignment;
        let start = page * PAGE_COLUMNS;
        let end = (start + PAGE_COLUMNS).min(a.reference_aligned.len());
        let mut result = Page {
            start,
            end,
            reference: a.reference_aligned[start..end].into(),
            query: a.query_aligned[start..end].into(),
            reference_positions: vec![],
            query_positions: vec![],
            tracks: vec![],
        };
        for (rb, qb) in result.reference.bytes().zip(result.query.bytes()) {
            result.reference_positions.push(
                (rb != b'-')
                    .then(|| (r + report.comparison.reference_origin) % report.reference_length),
            );
            result.query_positions.push((qb != b'-').then(|| {
                if report.comparison.query_reversed {
                    report.query_length - 1 - q
                } else {
                    q
                }
            }));
            r += usize::from(rb != b'-');
            q += usize::from(qb != b'-');
        }
        for track in &self.tracks {
            let (positions, reversed) = match track.side {
                Side::Reference => (&result.reference_positions, false),
                Side::Query => (&result.query_positions, report.comparison.query_reversed),
            };
            let mut annotation = vec![b' '; end - start];
            for (column, position) in positions.iter().enumerate() {
                if let Some(position) = position {
                    for (a, b, reverse) in &track.spans {
                        if a <= position && position < b {
                            mark(
                                &mut annotation[column],
                                if *reverse ^ reversed { b'<' } else { b'>' },
                            );
                        }
                    }
                }
            }
            if annotation.iter().all(|c| *c == b' ') {
                continue;
            }
            let mut protein = vec![b' '; end - start];
            for (column, position) in positions.iter().enumerate() {
                if let Some(position) = position {
                    if let Ok(index) = track
                        .translation_calls
                        .binary_search_by_key(position, |(p, _)| *p)
                    {
                        protein[column] = track.translation_calls[index].1;
                    }
                }
            }
            let mut notes = track.notes.clone();
            if annotation.contains(&b'!') || protein.contains(&b'!') {
                notes.push("Overlapping/repeated source mappings disagree at '!'; no single residue/direction is implied there.".into());
            }
            result.tracks.push(PageTrack {
                side: track.side,
                feature_index: track.index,
                label: track.label.clone(),
                kind: track.kind.clone(),
                annotation: String::from_utf8(annotation).expect("ASCII track symbols"),
                translation: (!track.translation_calls.is_empty())
                    .then(|| String::from_utf8(protein).expect("native amino-acid symbols")),
                notes,
            });
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::{
        self, AlignmentOptions, ComparisonOptions, MoleculeComparison,
        review::{self as pairwise_review, Options},
    };
    use collagene_model::Topology;
    fn cds(location: Location) -> Feature {
        Feature {
            kind: "CDS".into(),
            location,
            qualifiers: vec![("label".into(), Some("coding feature".into()))],
            needs_review: false,
        }
    }
    fn identity(source: &Molecule) -> Report {
        pairwise_review::compare(
            source,
            source,
            &Options::Stored(AlignmentOptions::default()),
        )
        .unwrap()
    }
    #[test]
    fn full_cds_codon_phase_survives_page_boundary() {
        let mut m = Molecule::dna(
            "source",
            &format!("{}ATGAAATAACC", "C".repeat(598)),
            Topology::Linear,
        )
        .unwrap();
        m.features.push(cds(Location::range(598, 607)));
        let report = identity(&m);
        let context = prepare(&report, &m, &m, true).unwrap();
        let first = context.page(&report, 0).unwrap();
        let second = context.page(&report, 1).unwrap();
        assert_eq!(&first.tracks[0].translation.as_ref().unwrap()[598..], ".M");
        assert_eq!(
            &second.tracks[0].translation.as_ref().unwrap()[..7],
            "..K..*."
        );
        assert_eq!(second.reference_positions[0], Some(600));
        assert_eq!(second.tracks[0].feature_index, 0);
    }
    #[test]
    fn circular_reference_and_reverse_query_share_literal_codon_columns() {
        // Source CDS joins original 9..12 then0..6: ATG AAA TAA.
        // Display is rotated to9. Query's original reverse CDS is complement3..12.
        let mut r = Molecule::dna("circular", "AAATAACCCATG", Topology::Circular).unwrap();
        r.features.push(cds(Location::Join(vec![
            Location::range(9, 12),
            Location::range(0, 6),
        ])));
        let mut q = Molecule::dna("reverse", "GGGTTATTTCAT", Topology::Linear).unwrap();
        q.features
            .push(cds(Location::Complement(Box::new(Location::range(3, 12)))));
        let options = Options::Whole(ComparisonOptions::default());
        let alignment = alignment::align(
            &collagene_model::Sequence::dna("ATGAAATAACCC").unwrap(),
            &q.sequence.reverse_complement().unwrap(),
            &options.alignment(),
        )
        .unwrap();
        let report = Report {
            engine_version: pairwise_review::ENGINE_VERSION.into(),
            context_revision: pairwise_review::revision(&r, &q, &options).unwrap(),
            reference_revision: r.revision(),
            query_revision: q.revision(),
            reference_length: 12,
            query_length: 12,
            options,
            comparison: MoleculeComparison {
                evaluated_cells: alignment.search.evaluated_cells,
                candidate_evidence: vec![alignment::grid::CandidateEvidence {
                    reference_origin: 9,
                    query_reversed: true,
                    score: Some(alignment.score),
                    search: Some(alignment.search.clone()),
                    disposition: if alignment.search.unrestricted_optimum_proven {
                        alignment::grid::CandidateDisposition::CertifiedOwnOptimum
                    } else {
                        alignment::grid::CandidateDisposition::Limited
                    },
                }],
                alignment,
                reference_revision: r.revision(),
                query_revision: q.revision(),
                reference_origin: 9,
                query_reversed: true,
                reference_topology: r.topology,
                query_topology: q.topology,
                evaluated_candidates: 1,
                exhaustive_reference_origins: false,
                exhaustive_query_orientations: false,
                optimality_proven: false,
                warnings: vec!["Declared test placement; no exhaustive search claim.".into()],
            },
        };
        let page = prepare(&report, &r, &q, true)
            .unwrap()
            .page(&report, 0)
            .unwrap();
        assert_eq!(
            page.reference_positions,
            vec![
                Some(9),
                Some(10),
                Some(11),
                Some(0),
                Some(1),
                Some(2),
                Some(3),
                Some(4),
                Some(5),
                Some(6),
                Some(7),
                Some(8)
            ]
        );
        assert_eq!(
            page.query_positions,
            (0..12).rev().map(Some).collect::<Vec<_>>()
        );
        for track in &page.tracks {
            assert_eq!(&track.annotation[..9], ">>>>>>>>>");
            assert_eq!(track.translation.as_deref(), Some(".M..K..*.   "));
        }
    }
    #[test]
    fn unsupported_cds_keeps_annotation_and_explicit_reason_without_protein() {
        let mut m = Molecule::dna("source", "ATGAAATAA", Topology::Linear).unwrap();
        let mut feature = cds(Location::range(0, 9));
        feature
            .qualifiers
            .push(("transl_table".into(), Some("99".into())));
        m.features.push(feature);
        let report = identity(&m);
        let context = prepare(&report, &m, &m, true).unwrap();
        assert!(
            context
                .warnings
                .iter()
                .any(|s| s.contains("genetic code 99"))
        );
        let page = context.page(&report, 0).unwrap();
        assert_eq!(page.tracks[0].annotation, ">>>>>>>>>");
        assert!(page.tracks[0].translation.is_none());
        let mut changed = m.clone();
        changed.name = "different revision".into();
        assert!(prepare(&report, &changed, &m, true).is_err());
        let mut changed_report = report;
        changed_report.context_revision = "0".repeat(64);
        assert!(context.page(&changed_report, 0).is_err());
    }

    #[test]
    fn gap_columns_do_not_shift_the_original_codon_start() {
        // Original CDS skips its first C: ATG AAA TAA. An insertion in the
        // query makes alignment column4 a reference gap; codon centres still
        // belong to original reference positions2,5,8, at columns2,6,9.
        let mut r = Molecule::dna("offset CDS", "CATGAAATAA", Topology::Linear).unwrap();
        let mut feature = cds(Location::range(0, 10));
        feature
            .qualifiers
            .push(("codon_start".into(), Some("2".into())));
        r.features.push(feature);
        let q = Molecule::dna("inserted query", "CATGCAAATAA", Topology::Linear).unwrap();
        let report =
            pairwise_review::compare(&r, &q, &Options::Stored(AlignmentOptions::default()))
                .unwrap();
        assert_eq!(report.comparison.alignment.reference_aligned, "CATG-AAATAA");
        let page = prepare(&report, &r, &q, true)
            .unwrap()
            .page(&report, 0)
            .unwrap();
        assert_eq!(page.reference_positions[4], None);
        assert_eq!(page.reference_positions[6], Some(5));
        assert_eq!(page.tracks[0].translation.as_deref(), Some(" .M. .K..*."));
        assert!(
            page.tracks[0]
                .notes
                .iter()
                .any(|note| note.contains("codon_start 2"))
        );
    }
}
