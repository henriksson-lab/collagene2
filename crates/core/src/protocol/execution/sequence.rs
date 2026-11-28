//! Sequence, annotations, translation and biological-analysis handlers.
use super::context::RequestContext;
use crate::protocol::schema::{AnalysisResult, AnnotationExportReport};
use collagene_model::{Error, Molecule};

pub(super) fn search(
    context: &RequestContext,
    molecule: Molecule,
    query: String,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::sequence::analysis::find_matches(&molecule, &query, true).map(AnalysisResult::Search)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn translate(
    context: &RequestContext,
    molecule: Molecule,
    options: crate::sequence::translation::TranslationOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::sequence::translation::translate(&molecule.sequence, &options)
            .map(AnalysisResult::Translation)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn translate_cds(
    context: &RequestContext,
    molecule: Molecule,
    feature_index: usize,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = molecule
        .features
        .get(feature_index)
        .ok_or_else(|| Error::Coordinates("CDS annotation no longer exists".into()))
        .and_then(|f| crate::sequence::translation::translate_cds(&molecule, f))
        .map(AnalysisResult::CdsTranslation);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn find_orfs(
    context: &RequestContext,
    molecule: Molecule,
    options: crate::sequence::translation::OrfOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        { crate::sequence::translation::find_orfs(&molecule, &options).map(AnalysisResult::Orfs) };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn scan_feature_library(
    context: &RequestContext,
    molecule: Molecule,
    library: crate::sequence::annotations::library::FeatureLibrary,
    options: crate::sequence::annotations::library::ScanOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::sequence::annotations::library::scan(&molecule, &library, &options)
            .map(AnalysisResult::FeatureLibraryScan);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn apply_feature_library(
    context: &RequestContext,
    molecule: Molecule,
    library: crate::sequence::annotations::library::FeatureLibrary,
    report: crate::sequence::annotations::library::ScanReport,
    selected: Vec<String>,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::sequence::annotations::library::apply_selected(
            &molecule, &library, &report, &selected,
        )
        .map(AnalysisResult::AnnotatedDocument);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn import_annotations(
    context: &RequestContext,
    molecule: Molecule,
    text: String,
    format: crate::formats::annotations::AnnotationFormat,
    seqid: String,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::formats::annotations::import_annotations(&molecule, &text, format, &seqid)
            .map(AnalysisResult::AnnotationImport);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn apply_annotations(
    context: &RequestContext,
    molecule: Molecule,
    report: crate::formats::annotations::AnnotationImportReport,
    selected: Vec<usize>,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::formats::annotations::apply_annotations(&molecule, &report, &selected)
            .map(AnalysisResult::AnnotatedDocument);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn export_annotations(
    context: &RequestContext,
    molecule: Molecule,
    format: crate::formats::annotations::AnnotationFormat,
    seqid: String,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::formats::annotations::export_annotations(&molecule, format, &seqid).map(|r| {
            AnalysisResult::AnnotationExport(AnnotationExportReport {
                text: r.text,
                warnings: r.warnings,
            })
        });
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn search_document(
    context: &RequestContext,
    molecule: Molecule,
    options: crate::sequence::search::SearchOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::sequence::search::search_document(&molecule, &options)
            .map(AnalysisResult::DocumentSearch)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn protein_properties(
    context: &RequestContext,
    molecule: Molecule,
    options: crate::protein::properties::ProteinOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;

        let molecule: Option<&Molecule> = Some(molecule);
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::protein::properties::analyze_protein(&molecule, &options)
            .map(AnalysisResult::ProteinProperties)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_codons(
    context: &RequestContext,
    molecule: Molecule,
    options: crate::sequence::codon_design::DesignOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let options = &options;

        let molecule: Option<&Molecule> = {
            // Bound source evidence and derive the required protein before hashing.
            // These jobs must not take the generic molecule-revision path first.
            if crate::sequence::codon_design::source_revision(molecule, options)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "codon request source revision mismatch".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::sequence::codon_design::design(&molecule, &options).map(AnalysisResult::CodonDesign)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn check_codon_design(
    context: &RequestContext,
    source: Molecule,
    candidate: Molecule,
    options: crate::sequence::codon_design::DesignOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &source;
        let options = &options;

        let molecule: Option<&Molecule> = {
            // Bound source evidence and derive the required protein before hashing.
            // These jobs must not take the generic molecule-revision path first.
            if crate::sequence::codon_design::source_revision(molecule, options)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "codon request source revision mismatch".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::sequence::codon_design::check(&source, &candidate, &options)
            .map(AnalysisResult::CodonCheck);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_crispr_guides(
    context: &RequestContext,
    molecule: Molecule,
    options: crate::sequence::crispr::GuideOptions,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::sequence::crispr::guide_revision(molecule, options)? != (*context.revision) {
                return Err(Error::Format("CRISPR guide source/options changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::sequence::crispr::design_guides(&molecule, &options)
            .map(AnalysisResult::CrisprGuides)
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn search_crispr_off_targets(
    context: &RequestContext,
    query: crate::sequence::crispr::GuideQuery,
    references: Vec<crate::sequence::crispr::ReferenceInput>,
    options: crate::sequence::crispr::OffTargetOptions,
) -> Result<AnalysisResult, Error> {
    {
        let query = &query;
        let references = &references;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::sequence::crispr::off_target_revision(query, references, options)?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "CRISPR off-target query/references/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::sequence::crispr::search_off_targets(&query, &references, &options)
            .map(AnalysisResult::CrisprOffTargets);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn design_edit_template(
    context: &RequestContext,
    source: Molecule,
    target: Molecule,
    options: crate::sequence::template::EditTemplateOptions,
    guide: Option<crate::sequence::template::GuideBinding>,
) -> Result<AnalysisResult, Error> {
    {
        let source = &source;
        let target = &target;
        let options = &options;
        let guide = &guide;

        let molecule: Option<&Molecule> = {
            if crate::sequence::template::batch_revision(source, target, options, guide.as_ref())?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "edit-template source/target/options/guide changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::sequence::template::design(&source, &target, &options, guide.as_ref())
            .map(AnalysisResult::EditTemplate);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn check_edit_template(
    context: &RequestContext,
    source: Molecule,
    target: Molecule,
    options: crate::sequence::template::EditTemplateOptions,
    guide: Option<crate::sequence::template::GuideBinding>,
    donor: String,
) -> Result<AnalysisResult, Error> {
    {
        let source = &source;
        let target = &target;
        let options = &options;
        let guide = &guide;
        let donor = &donor;

        let molecule: Option<&Molecule> = {
            if crate::sequence::template::check_batch_revision(
                source,
                target,
                options,
                guide.as_ref(),
                donor,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "edit-template candidate/source/target/options/guide changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::sequence::template::check(&source, &target, &options, guide.as_ref(), &donor)
            .map(AnalysisResult::EditTemplateCheck);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn analyze_rna_chemistry(
    context: &RequestContext,
    molecule: Molecule,
    declaration: crate::rna_chemistry::Declaration,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let declaration = &declaration;

        let molecule: Option<&Molecule> = {
            if crate::rna_chemistry::request_revision(molecule, declaration)
                .map_err(|e| Error::Format(e.to_string()))?
                != (*context.revision)
            {
                return Err(Error::Format("RNA source/declaration changed".into()));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> =
        crate::rna_chemistry::analyze(&molecule, &declaration)
            .map(AnalysisResult::RnaChemistry)
            .map_err(|e| Error::Format(e.to_string()));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn screen_protein_liabilities(
    context: &RequestContext,
    molecule: Molecule,
    options: crate::protein_liabilities::Options,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let options = &options;

        let molecule: Option<&Molecule> = {
            if crate::protein_liabilities::request_revision(molecule, options)
                .map_err(|e| Error::Format(e.to_string()))?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Protein liability source/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = {
        crate::protein_liabilities::screen(&molecule, &options, &mut || false)
            .map(AnalysisResult::ProteinLiabilities)
            .map_err(|e| Error::Format(e.to_string()))
    };
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn analyze_antibody(
    context: &RequestContext,
    molecule: Molecule,
    options: crate::antibody_numbering::domain_candidates::Options,
    limits: crate::antibody_numbering::domain_candidates::Limits,
    selection: crate::antibody_numbering::domain_tracks::Selection,
) -> Result<AnalysisResult, Error> {
    {
        let molecule = &molecule;
        let options = &options;
        let limits = &limits;
        let selection = &selection;

        let molecule: Option<&Molecule> = {
            if crate::antibody_numbering::domain_tracks::request_revision(
                molecule, options, limits, selection,
            )
            .map_err(|e| Error::Format(e.to_string()))?
                != (*context.revision)
            {
                return Err(Error::Format(
                    "Antibody source/options/limits changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::antibody_numbering::domain_tracks::analyze(
        &molecule,
        &options,
        &limits,
        &selection,
        &mut || false,
    )
    .map(AnalysisResult::Antibody)
    .map_err(|e| Error::Format(e.to_string()));
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}

pub(super) fn analyze_alignment_antibody(
    context: &RequestContext,
    alignment: crate::formats::alignment::AlignmentDocument,
    source: crate::protein::antibody_alignment::Source,
    options: crate::antibody_numbering::domain_candidates::Options,
    limits: crate::antibody_numbering::domain_candidates::Limits,
    selection: crate::antibody_numbering::domain_tracks::Selection,
) -> Result<AnalysisResult, Error> {
    {
        let alignment = &alignment;
        let source = &source;
        let options = &options;
        let limits = &limits;
        let selection = &selection;

        let molecule: Option<&Molecule> = {
            if crate::protein::antibody_alignment::request_revision(
                alignment, source, options, limits, selection,
            )? != (*context.revision)
            {
                return Err(Error::Format(
                    "Alignment antibody source/options changed".into(),
                ));
            }
            None
        };
        context.admit(molecule)?;
    }
    let outcome: Result<AnalysisResult, Error> = crate::protein::antibody_alignment::analyze(
        &alignment,
        &source,
        &options,
        &limits,
        &selection,
        &mut || false,
    )
    .map(AnalysisResult::AlignmentAntibody);
    Ok(outcome.unwrap_or_else(|error| AnalysisResult::Error(error.to_string())))
}
