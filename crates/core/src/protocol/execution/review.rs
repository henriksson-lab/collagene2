//! Alignment, reads, comparison and review-evidence handlers.
use super::context::RequestContext;
use crate::protocol::bindings::alignment_revision;
use crate::protocol::schema::AnalysisResult;
use collagene_model::{Alphabet, Error, Molecule};

pub(super) fn align(
    context: &RequestContext,
    molecule: Molecule,
    query: Molecule,
    options: crate::alignment::AlignmentOptions,
) -> Result<AnalysisResult, Error> {
    let expected = crate::alignment::review::revision(
        &molecule,
        &query,
        &crate::alignment::review::Options::Stored(options),
    )?;
    if expected != *context.revision {
        return Err(Error::Format(
            "Pairwise alignment sources/settings context changed".into(),
        ));
    }
    context.admit(None)?;
    let outcome = (|| -> Result<AnalysisResult, Error> {
        let report = crate::alignment::align(&molecule.sequence, &query.sequence, &options)?;
        report.validate_search(molecule.sequence.len(), query.sequence.len())?;
        if report.options != options {
            return Err(Error::Format(
                "Pairwise alignment returned different options".into(),
            ));
        }
        Ok(AnalysisResult::Alignment(report))
    })();
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn validate_traces(
    context: &RequestContext,
    molecule: Molecule,
    reads: Vec<Molecule>,
    options: crate::review::reads::trace_validation::TraceValidationOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::review::reads::trace_validation::validate_reads(&molecule, &reads, &options)
            .map(AnalysisResult::TraceValidation);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn import_alignment(
    context: &RequestContext,
    bytes: Vec<u8>,
    format: crate::formats::alignment::AlignmentFormat,
    alphabet: Alphabet,
    name: String,
) -> Result<AnalysisResult, Error> {
    {
        let bytes = &bytes;

        let molecule: Option<&Molecule> = {
            if blake3::hash(bytes).to_hex().as_str() != (*context.revision) {
                return Err(Error::Format("alignment input digest mismatch".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::formats::alignment::import_alignment(&bytes, format, alphabet, &name)
            .map(AnalysisResult::AlignmentDocument);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn export_alignment(
    context: &RequestContext,
    document: crate::formats::alignment::AlignmentDocument,
    format: crate::formats::alignment::AlignmentFormat,
) -> Result<AnalysisResult, Error> {
    {
        let document = &document;

        let molecule: Option<&Molecule> = {
            if alignment_revision(document)? != (*context.revision) {
                return Err(Error::Format("alignment document revision mismatch".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::formats::alignment::export_alignment(&document, format)
            .map(AnalysisResult::AlignmentExport)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn extract_alignment_rows(
    context: &RequestContext,
    document: crate::formats::alignment::AlignmentDocument,
    row_ids: Vec<String>,
) -> Result<AnalysisResult, Error> {
    {
        let document = &document;

        let molecule: Option<&Molecule> = {
            if alignment_revision(document)? != (*context.revision) {
                return Err(Error::Format("alignment document revision mismatch".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = document
        .extract_rows(&row_ids)
        .map(AnalysisResult::AlignmentRows);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn save_alignment(
    context: &RequestContext,
    document: crate::formats::alignment::AlignmentDocument,
) -> Result<AnalysisResult, Error> {
    {
        let document = &document;

        let molecule: Option<&Molecule> = {
            if alignment_revision(document)? != (*context.revision) {
                return Err(Error::Format("alignment document revision mismatch".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::formats::alignment::save_alignment(&document).map(AnalysisResult::AlignmentSaved)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn load_alignment(
    context: &RequestContext,
    bytes: Vec<u8>,
) -> Result<AnalysisResult, Error> {
    {
        let bytes = &bytes;

        let molecule: Option<&Molecule> = {
            if blake3::hash(bytes).to_hex().as_str() != (*context.revision) {
                return Err(Error::Format("alignment input digest mismatch".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::formats::alignment::load_alignment(&bytes).map(AnalysisResult::AlignmentDocument)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn dot_plot(
    context: &RequestContext,
    molecule: Molecule,
    query: Molecule,
    options: crate::alignment::dotplot::DotPlotOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let query = &query;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::alignment::dotplot::source_revisions(molecule, query, options)?.0
                != (*context.revision)
            {
                return Err(Error::Format(
                    "dot plot request source revision mismatch".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::alignment::dotplot::analyze(&molecule, &query, &options)
            .map(AnalysisResult::DotPlot);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn assemble_contigs(
    context: &RequestContext,
    reads: Vec<crate::review::contigs::ReadInput>,
    options: crate::review::contigs::ContigOptions,
) -> Result<AnalysisResult, Error> {
    {
        let reads = &reads;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::review::contigs::batch_revision(reads, options)? != (*context.revision) {
                return Err(Error::Format("contig batch source/options changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        { crate::review::contigs::assemble(&reads, &options).map(AnalysisResult::Contigs) };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn capture_read_check(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    reference: crate::material::lineage::EventBinding,
    reads: Vec<crate::material::lineage::EventBinding>,
    options: crate::review::reads::trace_validation::TraceValidationOptions,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let reference = &reference;
        let reads = &reads;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::review::reads::checks::capture_revision(workspace, reference, reads, options)?
                != (*context.revision)
            {
                return Err(Error::Format("Read-check capture context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::review::reads::checks::capture(&workspace, reference, reads, &options)
            .map(AnalysisResult::ReadCheckCaptured);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn replay_read_check(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    attachment_hash: String,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let attachment_hash = &attachment_hash;

        let molecule: Option<&Molecule> = {
            if crate::review::reads::checks::replay_revision(workspace, attachment_hash)?
                != (*context.revision)
            {
                return Err(Error::Format("Read-check replay context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::review::reads::checks::replay(&workspace, &attachment_hash)
            .map(AnalysisResult::ReadCheckReplayed);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn preview_read_review(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    request: crate::review::reads::review::Request,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::review::reads::review::preview_revision(workspace, request)?
                != (*context.revision)
            {
                return Err(Error::Format("Read review preview context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::review::reads::review::preview(&workspace, &request)
            .map(AnalysisResult::ReadReviewPreview)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn commit_read_review(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
    request: crate::review::reads::review::Request,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::review::reads::review::commit_revision(workspace, request)?
                != (*context.revision)
            {
                return Err(Error::Format("Read review commit context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::review::reads::review::commit(&workspace, &request)
            .map(AnalysisResult::ReadReviewCommitted)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn compare_sequences(
    context: &RequestContext,
    reference: Molecule,
    query: Molecule,
    options: crate::alignment::review::Options,
) -> Result<AnalysisResult, Error> {
    {
        let reference = &reference;
        let query = &query;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::alignment::review::revision(reference, query, options)? != (*context.revision)
            {
                return Err(Error::Format(
                    "Comparison sources or options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = (|| {
        let report = crate::alignment::review::compare(&reference, &query, &options)?;
        report.validate_binding(context.revision)?;
        report.validate_sources(&reference, &query)?;
        if report.options != options {
            return Err(Error::Format(
                "Pairwise comparison returned different options".into(),
            ));
        }
        Ok(AnalysisResult::ReviewedComparison(report))
    })();
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn export_pairwise_comparison(
    context: &RequestContext,
    record: crate::alignment::history::Record,
    reference_name: String,
    query_name: String,
    fresh: bool,
    format: crate::alignment::export::Format,
) -> Result<AnalysisResult, Error> {
    {
        let record = &record;
        let reference_name = &reference_name;
        let query_name = &query_name;
        let fresh = &fresh;
        let format = &format;

        let molecule: Option<&Molecule> = {
            if crate::alignment::export::revision(
                record,
                reference_name,
                query_name,
                *fresh,
                *format,
            )? != (*context.revision)
            {
                return Err(Error::Format("Pairwise export context changed.".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::alignment::export::export(&record, &reference_name, &query_name, fresh, format)
            .map(AnalysisResult::ExportedPairwiseComparison);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn search_pairwise_comparison(
    context: &RequestContext,
    report: crate::alignment::review::Report,
    pattern: String,
    options: crate::alignment::search::Options,
) -> Result<AnalysisResult, Error> {
    {
        let report = &report;
        let pattern = &pattern;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::alignment::search::revision(report, pattern, options)? != (*context.revision)
            {
                return Err(Error::Format("Pairwise search context changed.".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::alignment::search::search(&report, &pattern, options)
            .map(AnalysisResult::SearchedPairwiseComparison);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_pairwise_context(
    context: &RequestContext,
    report: crate::alignment::review::Report,
    reference: Molecule,
    query: Molecule,
    translate: bool,
) -> Result<AnalysisResult, Error> {
    {
        let report = &report;
        let reference = &reference;
        let query = &query;
        let translate = &translate;

        let molecule: Option<&Molecule> = {
            if crate::alignment::context::revision(report, reference, query, *translate)?
                != (*context.revision)
            {
                return Err(Error::Format("Pairwise annotation context changed.".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::alignment::context::prepare(&report, &reference, &query, translate)
            .map(AnalysisResult::PreparedPairwiseContext);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn find_similar_sequences(
    context: &RequestContext,
    query: crate::sequence::similarity::Source,
    subjects: Vec<crate::sequence::similarity::Source>,
    options: crate::sequence::similarity::Options,
) -> Result<AnalysisResult, Error> {
    {
        let query = &query;
        let subjects = &subjects;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::sequence::similarity::context_revision(query, subjects, options)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Sequence similarity sources/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::sequence::similarity::search(&query, &subjects, &options, &mut || false)
            .map(AnalysisResult::SequenceSimilarity);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
