//! Document interchange and bounded reference-transfer handlers.
use super::context::RequestContext;
use crate::protocol::schema::AnalysisResult;
use collagene_model::{Error, Molecule};

pub(super) fn save_workspace(
    context: &RequestContext,
    workspace: crate::documents::workspace::Workspace,
) -> Result<AnalysisResult, Error> {
    {
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            workspace.validate()?;
            let head = workspace.head(&workspace.active_document)?;
            if head.revision() != (*context.revision) {
                return Err(Error::Format(
                    "workspace worker input revision mismatch".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::documents::workspace::save_workspace(&workspace).map(AnalysisResult::WorkspaceSaved)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn load_workspace(
    context: &RequestContext,
    bytes: Vec<u8>,
) -> Result<AnalysisResult, Error> {
    {
        let bytes = &bytes;

        let molecule: Option<&Molecule> = {
            if bytes.len() > crate::documents::workspace::MAX_ARCHIVE_BYTES {
                return Err(Error::TooLarge);
            }
            if blake3::hash(bytes).to_hex().as_str() != (*context.revision) {
                return Err(Error::Format(
                    "workspace archive worker digest mismatch".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::documents::workspace::load_workspace(&bytes).map(AnalysisResult::WorkspaceLoaded)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn copy_reference_fragment(
    context: &RequestContext,
    request: crate::sequence::reference_transfer::RequestCopy,
) -> Result<AnalysisResult, Error> {
    {
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::sequence::reference_transfer::copy_revision(request)? != (*context.revision) {
                return Err(Error::Format(
                    "Reference copy source or selection changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::sequence::reference_transfer::copy(&request).map(AnalysisResult::ReferenceFragment)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn paste_reference_fragment(
    context: &RequestContext,
    request: crate::sequence::reference_transfer::RequestPaste,
) -> Result<AnalysisResult, Error> {
    {
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::sequence::reference_transfer::paste_revision(request)? != (*context.revision)
            {
                return Err(Error::Format(
                    "Reference paste sources or selection changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::sequence::reference_transfer::paste(&request).map(AnalysisResult::ReferencePaste)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn translate_reference(
    context: &RequestContext,
    request: crate::sequence::translation_reference::Request,
) -> Result<AnalysisResult, Error> {
    {
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::sequence::translation_reference::inspect_revision(request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Translation source or selection changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::sequence::translation_reference::inspect(&request)
            .map(AnalysisResult::TranslatedReference)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_translated_reference(
    context: &RequestContext,
    request: crate::sequence::translation_reference::Request,
    name: String,
) -> Result<AnalysisResult, Error> {
    {
        let request = &request;
        let name = &name;

        let molecule: Option<&Molecule> = {
            if crate::sequence::translation_reference::output_revision(request, name)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Translation source, selection or output name changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::sequence::translation_reference::prepare(&request, &name)
            .map(AnalysisResult::TranslatedOutputs)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn review_sequence_export(
    context: &RequestContext,
    molecule: Molecule,
    format: crate::formats::sequence::export::Format,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let format = &format;

        let molecule: Option<&Molecule> = {
            if crate::formats::sequence::export::revision(molecule, format)? != (*context.revision)
            {
                return Err(Error::Format(
                    "Sequence export source/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::formats::sequence::export::review(&molecule, &format)
            .map(AnalysisResult::SequenceExportAssessment)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_sequence_export(
    context: &RequestContext,
    molecule: Molecule,
    format: crate::formats::sequence::export::Format,
    expected_context_revision: String,
    allow_omissions: bool,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let format = &format;
        let expected_context_revision = &expected_context_revision;

        let molecule: Option<&Molecule> = {
            if expected_context_revision != &(*context.revision)
                || crate::formats::sequence::export::revision(molecule, format)?
                    != (*context.revision)
            {
                return Err(Error::Format(
                    "Reviewed sequence export context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::formats::sequence::export::prepare(
        &molecule,
        &format,
        &expected_context_revision,
        allow_omissions,
    )
    .map(AnalysisResult::SequenceExportPrepared);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn export_scf_trace(
    context: &RequestContext,
    trace: crate::review::reads::traces::Trace,
    options: crate::formats::scf::Options,
) -> Result<AnalysisResult, Error> {
    {
        let trace = &trace;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::formats::scf::revision(trace, options)? != (*context.revision) {
                return Err(Error::Format(
                    "Trace export source/view/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        { crate::formats::scf::export(&trace, &options).map(AnalysisResult::ScfTraceExported) };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
