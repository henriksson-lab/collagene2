//! Literal conformity to a captured, declared initiation-context pattern.
//! The host scope is an assumption; no initiation probability is calculated.
use super::*;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct PatternProfile {
    pub id: String,
    pub revision: String,
    pub host_context: String,
    pub genetic_code: u8,
    /// IUPAC pattern in biological reading direction, including literal ATG/AUG.
    pub pattern: Sequence,
    /// Zero-based index of the A in the initiating ATG/AUG, not a genomic offset.
    pub start_offset: usize,
    /// Documentary provenance, never authority to infer biological function.
    pub provenance: Citation,
}
impl PatternProfile {
    pub fn compute_revision(&self) -> Result<String, Error> {
        self.validate_shape()?;
        hash(&(
            ENGINE_VERSION,
            &self.id,
            &self.host_context,
            self.genetic_code,
            &self.pattern,
            self.start_offset,
            &self.provenance,
        ))
    }
    fn validate_shape(&self) -> Result<(), Error> {
        if !text(&self.id, 128)
            || !text(&self.host_context, 1024)
            || !matches!(self.pattern.alphabet(), Alphabet::Dna | Alphabet::Rna)
            || !(3..=128).contains(&self.pattern.len())
            || self.start_offset > self.pattern.len() - 3
            || !text(&self.provenance.source, 2048)
            || !self.provenance.source.starts_with("https://")
            || !text(&self.provenance.locator, 1024)
        {
            return Err(bad("invalid bounded initiation-pattern profile"));
        }
        translation::CodeDescriptor::for_id(self.genetic_code)?;
        let masks = self.pattern.as_str().as_bytes()[self.start_offset..self.start_offset + 3]
            .iter()
            .map(|b| crate::sequence::analysis::base_mask(*b))
            .collect::<Vec<_>>();
        if masks != [1, 8, 4] {
            return Err(bad(
                "initiation profile must explicitly contain ATG/AUG at its start offset",
            ));
        }
        Ok(())
    }
    pub fn validate(&self) -> Result<(), Error> {
        self.validate_shape()?;
        if !hex(&self.revision) || self.revision != self.compute_revision()? {
            return Err(bad("initiation-pattern content differs from its exact pin"));
        }
        Ok(())
    }
}

/// A documentary comparison to the experimentally studied ACCATGG sequence.
/// This is one declared mammalian AUG-context pattern, not a universal requirement
/// or a prediction that conformity is necessary/sufficient for expression.
pub fn kozak_1986_pattern() -> Result<PatternProfile, Error> {
    let mut profile = PatternProfile {
        id: "kozak-1986-accatgg".into(),
        revision: String::new(),
        host_context: "Declared mammalian AUG-initiation context; experimental sequence comparison only".into(),
        genetic_code: 1,
        pattern: Sequence::dna("ACCATGG")?,
        start_offset: 3,
        provenance: Citation {
            source: "https://doi.org/10.1016/0092-8674(86)90762-2".into(),
            locator: "Kozak (1986), Cell 44:283–292; cloned preproinsulin single-base substitution study. Literal ACCATGG conformity only; no strength score or universal host requirement.".into(),
        },
    };
    profile.revision = profile.compute_revision()?;
    Ok(profile)
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ContextClaim {
    pub cds_feature: usize,
    pub profile: PatternProfile,
    /// Explicit applicability declaration, not proof of transcript or host context.
    pub rationale: String,
}
impl ContextClaim {
    pub(super) fn validate(&self, m: &Molecule) -> Result<(), Error> {
        if self.cds_feature >= m.features.len() || !text(&self.rationale, 1024) {
            return Err(bad("invalid initiation-context declaration"));
        }
        self.profile.validate()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ContextEvidence {
    pub cds_feature: usize,
    pub profile_id: String,
    pub profile_revision: String,
    pub host_context: String,
    pub provenance: Citation,
    pub verdict: Verdict,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub observed: Option<String>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Location, 1>")]
    pub coordinates: Vec<Location>,
    pub reason: String,
}

fn window(
    m: &Molecule,
    start: usize,
    end: usize,
    reverse: bool,
    profile: &PatternProfile,
) -> Option<Location> {
    let n = m.sequence.len();
    let length = profile.pattern.len();
    if length > n {
        return None;
    }
    // n<=2M and pattern<=128 are admitted before these signed operations.
    let low = if reverse {
        end as isize + profile.start_offset as isize - length as isize
    } else {
        start as isize - profile.start_offset as isize
    };
    if m.topology == Topology::Linear && (low < 0 || low + length as isize > n as isize) {
        return None;
    }
    let low = low.rem_euclid(n as isize) as usize;
    let location = if low + length <= n {
        Location::range(low, low + length)
    } else {
        Location::Join(vec![
            Location::range(low, n),
            Location::range(0, low + length - n),
        ])
    };
    Some(if reverse {
        Location::Complement(Box::new(location))
    } else {
        location
    })
}

fn evaluate(checks: &Checks<'_>, claim: &ContextClaim) -> Result<ContextEvidence, Error> {
    let profile = &claim.profile;
    let mut result = ContextEvidence {
        cds_feature: claim.cds_feature,
        profile_id: profile.id.clone(),
        profile_revision: profile.revision.clone(),
        host_context: profile.host_context.clone(),
        provenance: profile.provenance.clone(),
        verdict: Verdict::Indeterminate,
        observed: None,
        coordinates: vec![],
        reason: "A current complete, explicitly assumed CDS on an exact single interval is required. Compound/fuzzy/partial paths and unreviewed annotations do not establish the initiation boundary or mature transcript context.".into(),
    };
    let feature = &checks.m.features[claim.cds_feature];
    if !checks.assumed(claim.cds_feature) || !feature.kind.eq_ignore_ascii_case("CDS") {
        return Ok(result);
    }
    let Some(frame) = checks.frames.get(&claim.cds_feature) else {
        return Ok(result);
    };
    let Some((start, end, reverse)) = geometry::single(&feature.location, false) else {
        return Ok(result);
    };
    if end - start < 3 || frame.offset != 0 || frame.genetic_code != profile.genetic_code {
        result.reason = "The complete CDS initiation boundary/genetic code does not match the pinned profile applicability declaration.".into();
        return Ok(result);
    }
    let Some(location) = window(checks.m, start, end, reverse, profile) else {
        result.reason = "The full declared context window is unavailable at this linear boundary or would reuse a circular source base; missing sequence is unresolved.".into();
        return Ok(result);
    };
    let observed = location.extract(&checks.m.sequence)?;
    let bytes = observed.as_bytes();
    let start_masks = bytes[profile.start_offset..profile.start_offset + 3]
        .iter()
        .map(|b| crate::sequence::analysis::base_mask(*b))
        .collect::<Vec<_>>();
    result.coordinates.push(location);
    result.observed = Some(observed.clone());
    if start_masks != [1, 8, 4] {
        result.reason = "Observed start is not a definite ATG/AUG. An AUG-specific context pattern is inapplicable or unresolved, including recognized alternative initiators; no nonconformity verdict is inferred.".into();
        return Ok(result);
    }
    let mut possible = false;
    let mut disjoint = false;
    for (actual, expected) in bytes.iter().zip(profile.pattern.as_str().bytes()) {
        let actual = crate::sequence::analysis::base_mask(*actual);
        let expected = crate::sequence::analysis::base_mask(expected);
        disjoint |= actual & expected == 0;
        possible |= actual & !expected != 0;
    }
    result.verdict = if disjoint {
        Verdict::Fail
    } else if possible {
        Verdict::Indeterminate
    } else {
        Verdict::Pass
    };
    result.reason = format!(
        "Observed oriented context {observed} {} the exact declared pattern {}. This is sequence-pattern conformity under the declared host/transcript assumption, not an initiation-strength, transcript-existence or expression finding.",
        if disjoint {
            "differs from"
        } else if possible {
            "only possibly conforms to"
        } else {
            "conforms to"
        },
        profile.pattern.as_str(),
    );
    Ok(result)
}

pub(super) fn inspect(checks: &mut Checks<'_>) -> Result<Vec<ContextEvidence>, Error> {
    let options = checks.o;
    let mut results = Vec::with_capacity(options.initiation_contexts.len());
    for claim in &options.initiation_contexts {
        let result = evaluate(checks, claim)?;
        // Charge both report evidence and the additional primary citation before copies.
        serde_json::to_writer(&mut checks.budget, &(&result, &result.provenance))
            .map_err(|_| Error::TooLarge)?;
        checks.push(
            Rule::InitiationContext,
            result.verdict.clone(),
            &[claim.cds_feature],
            &result.reason,
        )?;
        if let Some(finding) = checks.findings.last_mut() {
            finding.citations.push(result.provenance.clone());
        }
        results.push(result);
    }
    Ok(results)
}
pub(super) fn validate_reports(report: &Report) -> Result<(), Error> {
    if report.initiation_contexts.len() > 64 {
        return Err(Error::TooLarge);
    }
    for evidence in &report.initiation_contexts {
        if evidence.cds_feature >= report.source_feature_count
            || !text(&evidence.profile_id, 128)
            || !hex(&evidence.profile_revision)
            || !text(&evidence.host_context, 1024)
            || !text(&evidence.reason, 4096)
            || !text(&evidence.provenance.source, 2048)
            || !evidence.provenance.source.starts_with("https://")
            || !text(&evidence.provenance.locator, 1024)
            || evidence
                .observed
                .as_ref()
                .is_some_and(|s| s.len() > 128 || Sequence::dna(&s.replace('U', "T")).is_err())
            || evidence.coordinates.len() > 1
            || evidence.observed.is_some() != !evidence.coordinates.is_empty()
            || (evidence.observed.is_none() && evidence.verdict != Verdict::Indeterminate)
        {
            return Err(bad("invalid initiation-pattern evidence"));
        }
        for coordinate in &evidence.coordinates {
            if !interval_evidence_shape(coordinate) {
                return Err(bad("unsupported context evidence coordinate path"));
            }
            coordinate.validate(report.source_length)?;
        }
    }
    Ok(())
}
