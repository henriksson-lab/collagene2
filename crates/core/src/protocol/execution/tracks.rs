//! Track views and explicit rejection of operations requiring persistent sessions.
use super::context::RequestContext;
use crate::protocol::schema::AnalysisResult;
use collagene_model::{Error, Molecule};

pub(super) fn view_tracks(
    context: &RequestContext,
    molecule: Molecule,
    window: crate::tracks::view::TrackWindow,
    options: crate::tracks::view::TrackOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::tracks::view::view_tracks(&molecule, window, &options)
            .map(AnalysisResult::ViewTracks);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_tracks(
    context: &RequestContext,
    molecule: Molecule,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "cached track jobs require an AnalysisSession".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn cached_view_tracks(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "cached track jobs require an AnalysisSession".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn begin_genome_fasta(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn push_genome_fasta(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn finish_genome_fasta(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn begin_genome_native(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn push_genome_native_block(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn finish_genome_native(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn abort_genome_import(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn register_genome(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn genome_interval_guides(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn start_genome_search(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn advance_genome_search(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn page_genome_search(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn cancel_genome_job(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn read_genome_interval(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = None;
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = Err(Error::Unsupported(
        "genome operations require a persistent private worker session".into(),
    ));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
