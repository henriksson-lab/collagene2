//! Conditional antibody analysis of one exact imported protein alignment row.
//! Alignment gaps have no residue index and never receive numbering labels.
use crate::{
    Alphabet, Error, Molecule,
    alignment_formats::AlignmentDocument,
    antibody_numbering::{self as n, domain_candidates as d, domain_tracks as t},
};
use serde::{Deserialize, Serialize};
pub const ENGINE: &str = "alignment-antibody-3";
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub alignment_revision: String,
    pub row_id: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine: String,
    pub source: Source,
    pub context_revision: String,
    pub column_count: usize,
    /// Index in this vector is the zero-based ungapped residue index; values
    /// are zero-based alignment columns. At most 2,000 values, even for long gaps.
    pub residue_columns: Vec<usize>,
    pub analysis: t::Report,
}
fn bad(message: impl ToString) -> Error {
    Error::Format(message.to_string())
}
fn extract(doc: &AlignmentDocument, row_id: &str) -> Result<(Source, Molecule, Vec<usize>), Error> {
    if doc.alphabet != Alphabet::Protein {
        return Err(bad(
            "Antibody alignment analysis requires an explicit protein alphabet",
        ));
    }
    let selected = doc
        .rows
        .iter()
        .find(|r| r.id == row_id)
        .ok_or_else(|| bad("Unknown alignment row"))?;
    if selected
        .aligned
        .bytes()
        .filter(|b| *b != b'-')
        .take(n::MAX_SOURCE_RESIDUES + 1)
        .count()
        > n::MAX_SOURCE_RESIDUES
    {
        return Err(bad(
            "Antibody alignment row exceeds 2,000 ungapped residues",
        ));
    }
    let molecule = doc
        .extract_rows(&[row_id.to_owned()])?
        .pop()
        .ok_or_else(|| bad("Missing alignment row"))?;
    n::source_revision(&molecule).map_err(bad)?;
    let source = Source {
        alignment_revision: molecule
            .metadata
            .get("alignment_current_blake3")
            .cloned()
            .ok_or_else(|| bad("Missing extracted alignment provenance"))?,
        row_id: row_id.into(),
    };
    let row = doc
        .rows
        .iter()
        .find(|r| r.id == row_id)
        .ok_or_else(|| bad("Unknown alignment row"))?;
    let columns = row
        .aligned
        .bytes()
        .enumerate()
        .filter_map(|(i, b)| (b != b'-').then_some(i))
        .collect::<Vec<_>>();
    Ok((source, molecule, columns))
}
pub fn source(doc: &AlignmentDocument, row_id: &str) -> Result<Source, Error> {
    Ok(extract(doc, row_id)?.0)
}
fn resolve(doc: &AlignmentDocument, pin: &Source) -> Result<(Molecule, Vec<usize>), Error> {
    let (actual, m, columns) = extract(doc, &pin.row_id)?;
    if actual != *pin {
        return Err(bad("Alignment or selected row revision changed"));
    }
    Ok((m, columns))
}
fn context(source: &Source, column_count: usize, native_context: &str) -> Result<String, Error> {
    let bytes = serde_json::to_vec(&(ENGINE, source, column_count, native_context)).map_err(bad)?;
    Ok(blake3::hash(&bytes).to_hex().to_string())
}
pub fn request_revision(
    doc: &AlignmentDocument,
    source: &Source,
    options: &d::Options,
    limits: &d::Limits,
    selection: &t::Selection,
) -> Result<String, Error> {
    let (m, _) = resolve(doc, source)?;
    context(
        source,
        doc.column_count(),
        &t::request_revision(&m, options, limits, selection).map_err(bad)?,
    )
}
pub fn analyze(
    doc: &AlignmentDocument,
    source: &Source,
    options: &d::Options,
    limits: &d::Limits,
    selection: &t::Selection,
    cancel: &mut impl FnMut() -> bool,
) -> Result<Report, Error> {
    if cancel() {
        return Err(bad(n::Error::Cancelled));
    }
    let (m, residue_columns) = resolve(doc, source)?;
    let analysis = t::analyze(&m, options, limits, selection, cancel).map_err(bad)?;
    let column_count = doc.column_count();
    Ok(Report {
        engine: ENGINE.into(),
        source: source.clone(),
        context_revision: context(source, column_count, &analysis.context_revision)?,
        column_count,
        residue_columns,
        analysis,
    })
}
impl Report {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.analysis
            .validate_binding(&self.analysis.context_revision)
            .map_err(bad)?;
        if self.engine != ENGINE
            || self.context_revision != expected
            || context(
                &self.source,
                self.column_count,
                &self.analysis.context_revision,
            )? != expected
            || self.source.alignment_revision.len() != 64
            || !self
                .source
                .alignment_revision
                .bytes()
                .all(|b| b.is_ascii_hexdigit())
            || self.source.row_id.is_empty()
            || self.source.row_id.len() > 4096
            || self.column_count == 0
            || self.column_count > crate::alignment_formats::MAX_ALIGNMENT_CELLS
            || self.residue_columns.is_empty()
            || self.residue_columns.len() > n::MAX_SOURCE_RESIDUES
            || self.residue_columns.windows(2).any(|w| w[0] >= w[1])
            || self
                .residue_columns
                .last()
                .is_some_and(|c| *c >= self.column_count)
            || self
                .analysis
                .tracks
                .iter()
                .any(|t| t.source_end > self.residue_columns.len())
        {
            return Err(bad("Alignment antibody report binding/shape changed"));
        }
        Ok(())
    }
    /// Recheck actual alignment provenance and the gap map, without rerunning DP.
    /// Scores/biological claims still require full source-based validate_for.
    pub fn validate_source(
        &self,
        doc: &AlignmentDocument,
        pin: &Source,
        options: &d::Options,
        limits: &d::Limits,
        selection: &t::Selection,
    ) -> Result<(), Error> {
        let (m, columns) = resolve(doc, pin)?;
        self.validate_binding(&context(
            pin,
            doc.column_count(),
            &t::request_revision(&m, options, limits, selection).map_err(bad)?,
        )?)?;
        if self.source != *pin
            || self.column_count != doc.column_count()
            || self.residue_columns != columns
            || self.analysis.source_revision != n::source_revision(&m).map_err(bad)?
            || self.analysis.evidence.options != *options
            || self.analysis.evidence.limits != *limits
        {
            return Err(bad("Alignment antibody source/projection changed"));
        }
        Ok(())
    }
    pub fn validate_for(
        &self,
        doc: &AlignmentDocument,
        pin: &Source,
        options: &d::Options,
        limits: &d::Limits,
        selection: &t::Selection,
    ) -> Result<(), Error> {
        self.validate_source(doc, pin, options, limits, selection)?;
        let (m, _) = resolve(doc, pin)?;
        self.analysis
            .validate_for(&m, options, limits, selection)
            .map_err(bad)
    }
    pub fn at_column(&self, candidate: usize, column: usize) -> Option<&t::Residue> {
        let index = self.residue_columns.binary_search(&column).ok()?;
        let track = self.analysis.tracks.get(candidate)?;
        track
            .residues
            .get(index.checked_sub(track.first_supported_source_index)?)
    }
}
