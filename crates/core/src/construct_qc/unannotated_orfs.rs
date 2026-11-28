//! Bounded candidates from the canonical ORF owner, compared to declared CDS paths.
//! An unannotated candidate is informational, never evidence of unintended expression.
use super::*;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ScanOptions {
    pub start_policy: translation::StartPolicy,
    pub incomplete: translation::IncompletePolicy,
    /// Explicit user threshold; no universal biological minimum is inferred.
    pub minimum_amino_acids: usize,
    pub max_candidates: usize,
    pub both_strands: bool,
}
impl ScanOptions {
    pub(super) fn validate(&self) -> Result<(), Error> {
        if !(1..=MAX_BASES / 3).contains(&self.minimum_amino_acids)
            || !(1..=128).contains(&self.max_candidates)
        {
            return Err(bad(
                "ORF QC requires a declared positive length threshold and 1–128 candidates",
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum AnnotationRelation {
    CoveredByDeclaredCds,
    NoDeclaredCdsCoverage,
    Unresolved,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Boundary {
    StopCodon,
    AmbiguousStop,
    LinearBoundary,
    OneTurnLimit,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Candidate {
    pub location: Location,
    pub reverse: bool,
    pub amino_acids: usize,
    pub nucleotide_length: usize,
    pub crosses_origin: bool,
    pub has_start_codon: bool,
    pub has_stop_codon: bool,
    pub partial_5prime: bool,
    pub partial_3prime: bool,
    pub boundary: Boundary,
    pub ambiguous_codons: usize,
    pub annotation_relation: AnnotationRelation,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub covering_cds: Option<usize>,
    pub reason: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ScanEvidence {
    pub options: ScanOptions,
    pub genetic_code: u8,
    /// Complete only for the declared canonical start/partial policies. This is
    /// not exhaustive coding-potential detection; possible ambiguous starts are
    /// outside the definite-start model even when this field is true.
    pub enumeration_complete: bool,
    pub ambiguity_present: bool,
    pub context_dependent_stops: bool,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Candidate, 128>")]
    pub candidates: Vec<Candidate>,
    pub reason: String,
}

fn annotation_relation(
    checks: &Checks<'_>,
    orf: &translation::Orf,
) -> (AnnotationRelation, Option<usize>, String) {
    let mut unresolved = false;
    let candidate_interval = geometry::single(&orf.location, false);
    for (index, feature) in checks.m.features.iter().enumerate() {
        if !feature.kind.eq_ignore_ascii_case("CDS") {
            continue;
        }
        // Equal exact biological paths also support origin-spanning annotations.
        // A current complete translator frame is required, not just a feature name.
        if checks.assumed(index)
            && checks
                .frames
                .get(&index)
                .is_some_and(|f| f.offset == 0 && f.genetic_code == checks.o.genetic_code)
            && feature.location == orf.location
        {
            return (AnnotationRelation::CoveredByDeclaredCds, Some(index), "Candidate has the same biological path as a current complete declared CDS in the selected code; annotation truth and expression remain unverified.".into());
        }
        let Some((cs, ce, cr)) = geometry::single(&feature.location, false) else {
            // Unsupported paths may overlap; absence of a simple mapping cannot
            // be interpreted as lack of an annotation covering this candidate.
            unresolved = true;
            continue;
        };
        let Some((os, oe, orf_reverse)) = candidate_interval else {
            unresolved = true;
            continue;
        };
        if cr != orf_reverse || os >= ce || cs >= oe {
            continue;
        }
        let Some(frame) = checks.frames.get(&index).filter(|_| checks.assumed(index)) else {
            unresolved = true;
            continue;
        };
        if frame.genetic_code != checks.o.genetic_code {
            unresolved = true;
            continue;
        }
        if os >= cs && oe <= ce {
            let offset = if cr { ce - oe } else { os - cs };
            if offset.checked_sub(frame.offset).is_some_and(|n| n % 3 == 0) {
                return (AnnotationRelation::CoveredByDeclaredCds, Some(index), "Candidate is wholly contained on the same declared CDS strand and codon register. Nested accepted starts are retained; this relation does not establish biological use of either start.".into());
            }
        }
    }
    if unresolved {
        (AnnotationRelation::Unresolved, None, "Coverage is unresolved because a potentially relevant CDS is unreviewed, lacks a current complete declaration, uses another genetic code or has an unsupported biological path. No claim of annotation absence is made.".into())
    } else {
        (AnnotationRelation::NoDeclaredCdsCoverage, None, "No current declared CDS covers this candidate's full interval on the same strand and codon register. This is a documentary annotation relation, not an unintended-gene or expression finding.".into())
    }
}

pub(super) fn inspect(checks: &mut Checks<'_>) -> Result<Option<ScanEvidence>, Error> {
    let Some(options) = checks.o.orf_scan.clone() else {
        return Ok(None);
    };
    // The canonical owner performs linear successor-boundary scans, one strand
    // at a time; source<=2M, results<=128, total materialized amino acids<=8M.
    // No translated proteins are copied into this QC report.
    let report = translation::find_orfs_validated(
        checks.m,
        &translation::OrfOptions {
            genetic_code: checks.o.genetic_code,
            start_policy: options.start_policy,
            incomplete: options.incomplete,
            min_amino_acids: options.minimum_amino_acids,
            max_orfs: options.max_candidates,
            both_strands: options.both_strands,
        },
        &checks.source_revision,
    )?;
    let ambiguity_present = checks
        .m
        .sequence
        .as_str()
        .bytes()
        .any(|b| !crate::sequence::analysis::base_mask(b).is_power_of_two());
    let context_dependent_stops = !translation::CodeDescriptor::for_id(checks.o.genetic_code)?
        .context_dependent_stop_codons
        .is_empty();
    let reason = format!(
        "Canonical ORF candidate enumeration {} under the declared start, partial, strand and minimum-length policies. Ambiguous possible starts are not enumerated; possible stops terminate uncertain candidates. No candidate or empty result establishes expression, absence of coding potential or biological unintendedness. Annotation coverage is a separate declaration-relative fact.",
        if report.truncated {
            "was limited and is incomplete"
        } else {
            "completed"
        },
    );
    checks.push_literal_finding(
        Rule::UnannotatedOrf,
        if report.truncated || ambiguity_present || context_dependent_stops {
            Verdict::Indeterminate
        } else {
            Verdict::Informational
        },
        &[],
        &reason,
    )?;
    let mut candidates = Vec::with_capacity(report.orfs.len());
    for orf in report.orfs {
        let (annotation_relation, covering_cds, reason) = annotation_relation(checks, &orf);
        let candidate = Candidate {
            location: orf.location,
            reverse: orf.reverse,
            amino_acids: orf.protein.len(),
            nucleotide_length: orf.nucleotide_length,
            crosses_origin: orf.crosses_origin,
            has_start_codon: orf.has_start_codon,
            has_stop_codon: orf.has_stop_codon,
            partial_5prime: orf.partial_5prime,
            partial_3prime: orf.partial_3prime,
            boundary: match orf.termination {
                translation::Termination::StopCodon => Boundary::StopCodon,
                translation::Termination::AmbiguousStop => Boundary::AmbiguousStop,
                translation::Termination::LinearBoundary => Boundary::LinearBoundary,
                translation::Termination::OneTurnLimit => Boundary::OneTurnLimit,
            },
            ambiguous_codons: orf.ambiguous_codons,
            annotation_relation,
            covering_cds,
            reason,
        };
        serde_json::to_writer(&mut checks.budget, &candidate).map_err(|_| Error::TooLarge)?;
        checks.push_literal_finding(
            Rule::UnannotatedOrf,
            if candidate.annotation_relation == AnnotationRelation::Unresolved {
                Verdict::Indeterminate
            } else {
                Verdict::Informational
            },
            std::slice::from_ref(&candidate.location),
            &candidate.reason,
        )?;
        candidates.push(candidate);
    }
    Ok(Some(ScanEvidence {
        options,
        genetic_code: checks.o.genetic_code,
        enumeration_complete: !report.truncated,
        ambiguity_present,
        context_dependent_stops,
        candidates,
        reason,
    }))
}

pub(super) fn validate_reports(report: &Report) -> Result<(), Error> {
    let Some(scan) = &report.orf_scan else {
        return Ok(());
    };
    scan.options.validate()?;
    translation::CodeDescriptor::for_id(scan.genetic_code)?;
    if scan.candidates.len() > scan.options.max_candidates || !text(&scan.reason, 4096) {
        return Err(bad("invalid ORF scan evidence"));
    }
    for candidate in &scan.candidates {
        if !interval_evidence_shape(&candidate.location) {
            return Err(bad("unsupported ORF evidence coordinate path"));
        }
        candidate.location.validate(report.source_length)?;
        if candidate.amino_acids > MAX_BASES / 3
            || candidate.nucleotide_length > report.source_length
            || candidate.ambiguous_codons > candidate.amino_acids + 1
            || !text(&candidate.reason, 4096)
            || candidate
                .covering_cds
                .is_some_and(|i| i >= report.source_feature_count)
            || (candidate.annotation_relation == AnnotationRelation::CoveredByDeclaredCds)
                != candidate.covering_cds.is_some()
        {
            return Err(bad("invalid ORF candidate evidence"));
        }
    }
    Ok(())
}
