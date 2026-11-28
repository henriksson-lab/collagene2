//! Primer, PCR and terminal-preparation handlers over native scientific owners.
use super::context::RequestContext;
use crate::protocol::schema::{AnalysisResult, PrimerEvidence};
use crate::sequence::analysis::Enzyme;
use collagene_model::{Error, Molecule};

pub(super) fn restriction(
    context: &RequestContext,
    molecule: Molecule,
    enzyme: Enzyme,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::sequence::analysis::restriction_sites(&molecule, &enzyme)
            .map(AnalysisResult::Restriction)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn pcr(
    context: &RequestContext,
    molecule: Molecule,
    forward: crate::sequence::analysis::Primer,
    reverse: crate::sequence::analysis::Primer,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::sequence::analysis::simulate_pcr(&molecule, &forward, &reverse)
            .map(AnalysisResult::Pcr);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn duplex_digest(
    context: &RequestContext,
    molecule: Molecule,
    enzymes: Vec<Enzyme>,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::chemistry::Duplex::from_molecule(&molecule)
            .and_then(|d| crate::chemistry::digest(&d, &enzymes))
            .map(AnalysisResult::DuplexDigest)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn pcr_documents(
    context: &RequestContext,
    molecule: Molecule,
    forward: crate::sequence::analysis::Primer,
    reverse: crate::sequence::analysis::Primer,
    template_strand: crate::chemistry::Strand,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::verification::pcr::simulate_pcr_documents(
        &molecule,
        &forward,
        &reverse,
        template_strand,
    )
    .map(AnalysisResult::PcrDocuments);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn analyze_primer(
    context: &RequestContext,
    molecule: Molecule,
    primer: crate::sequence::analysis::Primer,
    conditions: crate::cloning::pcr::primers::ThermoConditions,
    binding_options: crate::cloning::pcr::primers::BindingOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::pcr::primers::analyze_primer(&primer, &conditions).and_then(|analysis| {
            crate::cloning::pcr::primers::find_binding_sites(&molecule, &primer, &binding_options)
                .map(|bindings| {
                    AnalysisResult::PrimerEvidence(PrimerEvidence { analysis, bindings })
                })
        });
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_primers(
    context: &RequestContext,
    molecule: Molecule,
    options: crate::cloning::pcr::primers::DesignOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::pcr::primers::design_primers(&molecule, &options)
            .map(AnalysisResult::PrimerDesign)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn golden_gate(
    context: &RequestContext,
    molecule: Molecule,
    templates: Vec<Molecule>,
    options: crate::cloning::golden_gate::GoldenGateOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::golden_gate::simulate(&templates, &options).map(AnalysisResult::GoldenGate);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn mismatch_pcr(
    context: &RequestContext,
    molecule: Molecule,
    forward: crate::sequence::analysis::Primer,
    reverse: crate::sequence::analysis::Primer,
    options: crate::cloning::pcr::mutagenesis::PcrOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::pcr::mutagenesis::simulate_pcr(&molecule, &forward, &reverse, &options)
            .map(AnalysisResult::MismatchPcr)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inverse_pcr(
    context: &RequestContext,
    molecule: Molecule,
    target: Molecule,
    origin: Option<crate::material::run::MaterialOrigin>,
    options: crate::cloning::pcr::mutagenesis::InversePcrOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let target = &target;
        let origin = &origin;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::cloning::pcr::mutagenesis::inverse_pcr_revision(
                molecule,
                target,
                origin.as_ref(),
                options,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Inverse PCR source, origin, independent target or edit settings changed"
                        .into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::pcr::mutagenesis::design_inverse_pcr(
            &molecule,
            &target,
            origin.as_ref(),
            &options,
        )
        .map(AnalysisResult::InversePcr);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn mismatch_pcr_preview(
    context: &RequestContext,
    molecule: Molecule,
    forward: crate::sequence::analysis::Primer,
    reverse: crate::sequence::analysis::Primer,
    options: crate::cloning::pcr::mutagenesis::PcrOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::verification::pcr::simulate_mutagenesis_preview(
            &molecule, &forward, &reverse, &options,
        )
        .map(AnalysisResult::MismatchPcrPreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn repair_preview(
    context: &RequestContext,
    molecule: Molecule,
    enzymes: Vec<String>,
    fragment_index: usize,
    method: crate::chemistry::EndRepairMethod,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::verification::terminal::simulate_end_repair(
            &molecule,
            &enzymes,
            fragment_index,
            method,
        )
        .map(AnalysisResult::RepairPreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn anneal_preview(
    context: &RequestContext,
    molecule: Molecule,
    bottom: Molecule,
    options: crate::chemistry::AnnealOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::verification::terminal::simulate_oligo_annealing(&molecule, &bottom, &options)
            .map(AnalysisResult::AnnealPreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn annealed_document(
    context: &RequestContext,
    molecule: Molecule,
    bottom: Molecule,
    preview: crate::verification::terminal::AnnealPreview,
    candidate_index: usize,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::verification::terminal::annealed_document(
        &molecule,
        &bottom,
        &preview,
        candidate_index,
    )
    .map(AnalysisResult::PhysicalPreviewDocument);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn tail_preview(
    context: &RequestContext,
    molecule: Molecule,
    recipe: crate::verification::terminal::TailSubstrate,
    base: crate::cloning::ligation::tailing::TailingBase,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::verification::terminal::simulate_tailing(&molecule, &recipe, base)
            .map(AnalysisResult::TailPreview);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn assess_methylation(
    context: &RequestContext,
    molecule: Molecule,
    options: crate::chemistry::methylation::MethylationOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::chemistry::methylation::batch_revision(molecule, options)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "methylation source/options/profile changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::chemistry::methylation::assess(&molecule, &options).map(AnalysisResult::Methylation)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_pcr_outputs(
    context: &RequestContext,
    request: crate::cloning::pcr::outputs::Request,
) -> Result<AnalysisResult, Error> {
    {
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::pcr::outputs::revision(request)? != (*context.revision) {
                return Err(Error::Format("PCR output selection changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        { crate::cloning::pcr::outputs::prepare(&request).map(AnalysisResult::PcrOutputs) };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn prepare_mismatch_pcr_outputs(
    context: &RequestContext,
    request: crate::cloning::pcr::mismatch_outputs::Request,
) -> Result<AnalysisResult, Error> {
    {
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::pcr::mismatch_outputs::revision(request)? != (*context.revision) {
                return Err(Error::Format(
                    "Mismatch-PCR output selection changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::pcr::mismatch_outputs::prepare(&request)
            .map(AnalysisResult::MismatchPcrOutputs)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn score_ligation_fidelity(
    context: &RequestContext,
    dataset: crate::cloning::ligation::fidelity::Dataset,
    ends: Vec<String>,
) -> Result<AnalysisResult, Error> {
    {
        let dataset = &dataset;
        let ends = &ends;

        let molecule: Option<&Molecule> = {
            if crate::cloning::ligation::fidelity::revision(dataset, ends)? != (*context.revision) {
                return Err(Error::Format(
                    "Ligation fidelity data or ends changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::ligation::fidelity::score(&dataset, &ends)
            .map(AnalysisResult::LigationFidelity)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn list_ligation_profiles(context: &RequestContext) -> Result<AnalysisResult, Error> {
    {
        let molecule: Option<&Molecule> = {
            if crate::cloning::ligation::profiles::catalog_revision()? != (*context.revision) {
                return Err(Error::Format(
                    "Published ligation catalog context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        { crate::cloning::ligation::profiles::list().map(AnalysisResult::LigationProfiles) };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn score_published_ligation(
    context: &RequestContext,
    profile: crate::cloning::ligation::profiles::ProfileRef,
    ends: Vec<String>,
    physical_competitors: Vec<String>,
) -> Result<AnalysisResult, Error> {
    {
        let profile = &profile;
        let ends = &ends;
        let physical_competitors = &physical_competitors;

        let molecule: Option<&Molecule> = {
            if crate::cloning::ligation::profiles::score_revision_with_competitors(
                profile,
                ends,
                physical_competitors,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Published ligation profile/end context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::ligation::profiles::score_with_competitors(
            &profile,
            &ends,
            &physical_competitors,
        )
        .map(AnalysisResult::PublishedLigationScore);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn score_ligation_competition(
    context: &RequestContext,
    dataset: crate::cloning::ligation::fidelity::Dataset,
    ends: Vec<String>,
    physical_competitors: Vec<String>,
) -> Result<AnalysisResult, Error> {
    {
        let dataset = &dataset;
        let ends = &ends;
        let physical_competitors = &physical_competitors;

        let molecule: Option<&Molecule> = {
            if crate::cloning::ligation::fidelity::competition_revision(
                dataset,
                ends,
                physical_competitors,
            )? != (*context.revision)
            {
                return Err(Error::Format("Ligation competition context changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::ligation::fidelity::score_competition(
            &dataset,
            &ends,
            &physical_competitors,
        )
        .map(AnalysisResult::LigationCompetition);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_pcr_digest(
    context: &RequestContext,
    molecule: Molecule,
    request: crate::cloning::pcr::digest::Request,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let request = &request;

        let molecule: Option<&Molecule> = {
            if crate::cloning::pcr::digest::inspection_revision(molecule, request)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "PCR digest original source, primers or catalog options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::cloning::pcr::digest::inspect(&molecule, &request)
            .map(AnalysisResult::PcrDigestInspection)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_material_pcr_digest(
    context: &RequestContext,
    molecule: Molecule,
    origin: crate::material::run::MaterialOrigin,
    request: crate::cloning::pcr::digest::Request,
    workspace: crate::documents::workspace::Workspace,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let origin = &origin;
        let request = &request;
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::pcr::digest::inspection_revision_from_material(
                molecule, origin, request, workspace,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Saved-run PCR digest origin/source/primer/catalog context changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::pcr::digest::inspect_from_material(
            &molecule, &origin, &request, &workspace,
        )
        .map(AnalysisResult::MaterialPcrDigestInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_mismatch_pcr_digest(
    context: &RequestContext,
    molecule: Molecule,
    origin: Option<crate::material::run::MaterialOrigin>,
    request: crate::cloning::pcr::mismatch_digest::Request,
    workspace: Option<crate::documents::workspace::Workspace>,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let origin = &origin;
        let request = &request;
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::pcr::mismatch_digest::inspection_revision(
                molecule,
                origin.as_ref(),
                request,
                workspace.as_ref(),
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Mismatch PCR digest original source/origin/policy/catalog context changed"
                        .into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::cloning::pcr::mismatch_digest::inspect(
        &molecule,
        origin.as_ref(),
        &request,
        workspace.as_ref(),
    )
    .map(AnalysisResult::MismatchPcrDigestInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_overlap_extension_primers(
    context: &RequestContext,
    source: Molecule,
    request: crate::cloning::overlap::extension::primers::Request,
    workspace: Option<crate::documents::workspace::Workspace>,
) -> Result<AnalysisResult, Error> {
    {
        let source = &source;
        let request = &request;
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::overlap::extension::primers::design_revision(
                source,
                request,
                workspace.as_ref(),
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Overlap-extension primer source, independent edits or options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::overlap::extension::primers::design(&source, &request, workspace.as_ref())
            .map(AnalysisResult::OverlapExtensionPrimersDesigned);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn inspect_overlap_extension_digest(
    context: &RequestContext,
    source: Molecule,
    origin: Option<crate::material::run::MaterialOrigin>,
    request: crate::cloning::overlap::extension::digest::Request,
    workspace: Option<crate::documents::workspace::Workspace>,
) -> Result<AnalysisResult, Error> {
    {
        let source = &source;
        let origin = &origin;
        let request = &request;
        let workspace = &workspace;

        let molecule: Option<&Molecule> = {
            if crate::cloning::overlap::extension::digest::inspection_revision(
                source,
                origin.as_ref(),
                request,
                workspace.as_ref(),
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Overlap-extension digest source, origin, oligos or catalog/options changed"
                        .into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::cloning::overlap::extension::digest::inspect(
            &source,
            origin.as_ref(),
            &request,
            workspace.as_ref(),
        )
        .map(AnalysisResult::OverlapExtensionDigestInspection);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
