//! Source-bound reference translation and explicit, reviewable output preparation.
//! All amino-acid calls come from the existing native translation engine.
use crate::sequence::translation::{self, CdsTranslation, TranslationOptions, TranslationReport};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const ENGINE_VERSION: &str = "translation-reference-1";
const MAX_BASES: usize = collagene_model::MAX_SEQUENCE_LEN;
macro_rules! dto { ($($item:tt)*) => {
    #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
    #[cfg_attr(feature="json-schema", derive(schemars::JsonSchema))]
    #[serde(deny_unknown_fields)] $($item)*
}; }
dto! { pub struct Request { pub source: Molecule, pub selection: Selection } }
dto! { pub enum Selection {
    Region { location: Location, options: TranslationOptions },
    Cds { feature_index: usize },
} }
dto! { pub struct Report {
    pub context_revision: String,
    pub source_revision: String,
    pub selection: Selection,
    pub effective_location: Location,
    pub translation: TranslationReport,
    /// Literal region translation retains every '*'. CDS translation omits one
    /// definite terminal stop, following the existing native CDS semantics.
    pub protein: String,
    #[cfg_attr(feature="json-schema", schemars(with="Option<serde_json::Value>"))]
    pub cds: Option<CdsTranslation>,
    pub warnings: Vec<String>,
} }
dto! { pub struct Outputs {
    pub context_revision: String,
    pub report: Report,
    pub protein: Molecule,
    /// Region-only CDS proposal. Review is required before CDS-directed use.
    pub annotation: Option<Feature>,
} }

fn bad(message: &str) -> Error {
    Error::Unsupported(format!("reference translation: {message}"))
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn name_valid(name: &str) -> bool {
    !name.trim().is_empty() && name.len() <= 1024 && !name.chars().any(char::is_control)
}
fn options_valid(options: &TranslationOptions) -> Result<(), Error> {
    if options.frame > 2 || !matches!(options.genetic_code, 1 | 2 | 11) {
        return Err(bad("unsupported genetic code or frame"));
    }
    Ok(())
}
/// Validate borrowed geometry before any extraction, clone or recursive walk.
/// Exact repeated regions are permitted, with their total extracted length bounded.
fn location_length(location: &Location, length: usize) -> Result<usize, Error> {
    let mut stack = vec![(location, 1usize)];
    let mut nodes = 0usize;
    let mut bases = 0usize;
    while let Some((location, depth)) = stack.pop() {
        nodes += 1;
        if depth > 32 || nodes > 16_384 {
            return Err(Error::TooLarge);
        }
        match location {
            Location::Range {
                start,
                end,
                before: false,
                after: false,
            } if start < end && *end <= length => {
                bases = bases
                    .checked_add(end - start)
                    .filter(|n| *n <= MAX_BASES)
                    .ok_or(Error::TooLarge)?;
            }
            Location::Complement(inner) => stack.push((inner, depth + 1)),
            Location::Join(parts) if !parts.is_empty() => {
                if parts.len() > 16_384usize.saturating_sub(nodes + stack.len()) {
                    return Err(Error::TooLarge);
                }
                stack.extend(parts.iter().map(|part| (part, depth + 1)));
            }
            _ => {
                return Err(bad(
                    "translation selection requires exact nonempty ranges, joins or complements",
                ));
            }
        }
    }
    Ok(bases)
}
fn effective(location: &Location, reverse: bool) -> Location {
    match (location, reverse) {
        (Location::Complement(inner), true) => inner.as_ref().clone(),
        (_, true) => Location::Complement(Box::new(location.clone())),
        (_, false) => location.clone(),
    }
}
fn request_preflight(request: &Request) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[&request.source])?;
    request.source.validate()?;
    if request.source.sequence.alphabet() == Alphabet::Protein {
        return Err(bad("choose a DNA or RNA source"));
    }
    match &request.selection {
        Selection::Region { location, options } => {
            options_valid(options)?;
            location_length(location, request.source.sequence.len())?;
            // Reversing an already-maximal location adds one level; check the
            // effective tree after its bounded clone before extraction.
        }
        Selection::Cds { feature_index } => {
            let feature = request
                .source
                .features
                .get(*feature_index)
                .ok_or_else(|| bad("CDS feature is unavailable"))?;
            location_length(&feature.location, request.source.sequence.len())?;
        }
    }
    Ok(())
}
pub fn inspect_revision(request: &Request) -> Result<String, Error> {
    request_preflight(request)?;
    crate::cloning::assembly::recipe::hash(&(
        ENGINE_VERSION,
        request.source.revision(),
        &request.selection,
    ))
}
pub fn inspect(request: &Request) -> Result<Report, Error> {
    request_preflight(request)?;
    let source_revision = request.source.revision();
    let context_revision = crate::cloning::assembly::recipe::hash(&(
        ENGINE_VERSION,
        &source_revision,
        &request.selection,
    ))?;
    let (effective_location, translation, protein, cds, warnings) = match &request.selection {
        Selection::Region { location, options } => {
            let effective = effective(location, options.reverse);
            location_length(&effective, request.source.sequence.len())?;
            effective.validate(request.source.sequence.len())?;
            let bases = location.extract(&request.source.sequence)?;
            let sequence = Sequence::new(request.source.sequence.alphabet(), &bases)?;
            let translated = translation::translate(&sequence, options)?;
            let protein = translated.protein.clone();
            let warnings = translated.warnings.clone();
            (effective, translated, protein, None, warnings)
        }
        Selection::Cds { feature_index } => {
            let result = translation::translate_cds_validated(
                &request.source,
                &request.source.features[*feature_index],
                &source_revision,
            )?;
            let location = result.location.clone();
            let translation = result.translation.clone();
            let protein = result.protein.clone();
            let warnings = result.warnings.clone();
            (location, translation, protein, Some(result), warnings)
        }
    };
    let report = Report {
        context_revision,
        source_revision,
        selection: request.selection.clone(),
        effective_location,
        translation,
        protein,
        cds,
        warnings,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
pub fn output_revision(request: &Request, name: &str) -> Result<String, Error> {
    if !name_valid(name) {
        return Err(bad(
            "output name must be 1–1024 bytes without control characters",
        ));
    }
    crate::cloning::assembly::recipe::hash(&(
        ENGINE_VERSION,
        "outputs",
        inspect_revision(request)?,
        name,
    ))
}

fn annotation(report: &Report, name: &str) -> Option<Feature> {
    let Selection::Region { options, .. } = &report.selection else {
        return None;
    };
    let mut qualifiers = vec![
        ("label".into(), Some(name.into())),
        (
            "transl_table".into(),
            Some(options.genetic_code.to_string()),
        ),
        ("codon_start".into(), Some((options.frame + 1).to_string())),
        (
            "note".into(),
            Some(format!(
                "Collagene selected-region translation proposal; original initiation policy {:?}. Review required: future CDS-directed initiation and terminal-stop semantics can differ from this literal frame translation. No biological completeness or expression is asserted. Source revision {}; context {}.",
                options.initiation, report.source_revision, report.context_revision
            )),
        ),
    ];
    qualifiers.extend(report.warnings.iter().map(|warning| {
        (
            "note".into(),
            Some(format!("Translation warning: {warning}")),
        )
    }));
    Some(Feature {
        kind: "CDS".into(),
        location: report.effective_location.clone(),
        qualifiers,
        needs_review: true,
    })
}

/// Fresh native reconstruction from actual source calls; no submitted report or
/// stored protein qualifier can authorize the output. Source documents are unchanged.
pub fn prepare(request: &Request, name: &str) -> Result<Outputs, Error> {
    if !name_valid(name) {
        return Err(bad(
            "output name must be 1–1024 bytes without control characters",
        ));
    }
    let report = inspect(request)?;
    if report.protein.is_empty() {
        return Err(bad("translation contains no protein residues to open"));
    }
    let context_revision = crate::cloning::assembly::recipe::hash(&(
        ENGINE_VERSION,
        "outputs",
        &report.context_revision,
        name,
    ))?;
    let metadata = BTreeMap::from([
        ("translation:engine".into(), ENGINE_VERSION.into()),
        (
            "translation:source_revision".into(),
            report.source_revision.clone(),
        ),
        (
            "translation:context_revision".into(),
            report.context_revision.clone(),
        ),
        (
            "translation:output_revision".into(),
            context_revision.clone(),
        ),
        (
            "translation:policy".into(),
            match report.selection {
                Selection::Region { .. } => {
                    "Literal selected-region frame translation; all stops retained"
                }
                Selection::Cds { .. } => {
                    "Native CDS translation; one definite terminal stop omitted"
                }
            }
            .into(),
        ),
    ]);
    let protein = Molecule {
        name: name.into(),
        sequence: Sequence::new(Alphabet::Protein, &report.protein)?,
        topology: Topology::Linear,
        features: Vec::new(),
        metadata,
        original_genbank: None,
    };
    let proposed = annotation(&report, name);
    let outputs = Outputs {
        context_revision,
        report,
        protein,
        annotation: proposed,
    };
    outputs.validate_binding(&outputs.context_revision)?;
    Ok(outputs)
}

fn warnings_valid(warnings: &[String]) -> bool {
    warnings.len() <= 64 && warnings.iter().all(|s| s.len() <= 4096)
}
fn translation_valid(report: &TranslationReport) -> Result<(), Error> {
    options_valid(&report.options)?;
    let remaining = report
        .nucleotide_length
        .saturating_sub(usize::from(report.options.frame));
    if report.nucleotide_length == 0
        || report.nucleotide_length > MAX_BASES
        || report.protein.len() != remaining / 3
        || report.translated_codons != remaining / 3
        || report.trailing_bases != remaining % 3
        || report.code_name.len() > 128
        || !report
            .protein
            .bytes()
            .all(|b| b.is_ascii_uppercase() || b == b'*')
        || report.ambiguous_codons > report.translated_codons
        || report.unresolved_codons > report.ambiguous_codons
        || report.unresolved_codons != report.protein.bytes().filter(|&b| b == b'X').count()
        || report.possible_stop_codons > report.ambiguous_codons
        || report.initiator_methionine_applied
            && (report.options.initiation == translation::InitiationPolicy::Ordinary
                || !report.protein.starts_with('M'))
        || !warnings_valid(&report.warnings)
    {
        return Err(bad("invalid bounded translation report"));
    }
    Ok(())
}
impl Report {
    /// Read-only bounded consistency check. A matching context is not proof of
    /// protein correctness or experimental expression; native computation supplies calls.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !hex(&self.context_revision)
            || !hex(&self.source_revision)
            || self.context_revision != expected
            || !warnings_valid(&self.warnings)
            || self.protein.len() > MAX_BASES / 3
        {
            return Err(bad("translation report context/bounds mismatch"));
        }
        translation_valid(&self.translation)?;
        let selected_length = location_length(&self.effective_location, MAX_BASES)?;
        if selected_length != self.translation.nucleotide_length {
            return Err(bad("effective location differs from translated base count"));
        }
        if let Selection::Region { location, options } = &self.selection {
            location_length(location, MAX_BASES)?;
            options_valid(options)?;
        }
        if self.context_revision
            != crate::cloning::assembly::recipe::hash(&(
                ENGINE_VERSION,
                &self.source_revision,
                &self.selection,
            ))?
        {
            return Err(bad("translation selection binding mismatch"));
        }
        match (&self.selection, &self.cds) {
            (Selection::Region { location, options }, None) => {
                location_length(location, MAX_BASES)?;
                if self.effective_location != effective(location, options.reverse)
                    || options != &self.translation.options
                    || self.protein != self.translation.protein
                    || self.warnings != self.translation.warnings
                {
                    return Err(bad("region translation fields disagree"));
                }
            }
            (Selection::Cds { .. }, Some(cds)) => {
                location_length(&cds.location, MAX_BASES)?;
                if cds.feature_label.len() > 1_048_576
                    || !warnings_valid(&cds.warnings)
                    || cds.source_revision != self.source_revision
                    || cds.location != self.effective_location
                    || cds.translation != self.translation
                    || cds.protein != self.protein
                    || cds.warnings != self.warnings
                    || !(1..=3).contains(&cds.codon_start)
                    || cds.codon_start != self.translation.options.frame + 1
                    || cds.has_terminal_stop != self.translation.protein.ends_with('*')
                {
                    return Err(bad("CDS translation fields disagree"));
                }
                let protein = self
                    .translation
                    .protein
                    .strip_suffix('*')
                    .unwrap_or(&self.translation.protein);
                if self.protein != protein
                    || cds.internal_stops != protein.bytes().filter(|b| *b == b'*').count()
                {
                    return Err(bad("CDS stop semantics disagree"));
                }
                if cds.complete_under_model
                    && (cds.codon_start != 1
                        || !cds.has_terminal_stop
                        || cds.internal_stops != 0
                        || !self.translation.initiator_methionine_applied
                        || self.translation.trailing_bases != 0
                        || self.translation.possible_stop_codons != 0
                        || self.translation.unresolved_codons != 0
                        || cds.stored_translation_matches == Some(false))
                {
                    return Err(bad("CDS completeness contradicts retained codon findings"));
                }
            }
            _ => return Err(bad("translation selection/result kind mismatch")),
        }
        Ok(())
    }
    /// Recompute against actual source calls to detect even coordinated report
    /// tampering. Ordinary private-worker presentation can use validate_binding;
    /// prepare always reconstructs outputs independently of any submitted report.
    pub fn validate_sources(&self, request: &Request) -> Result<(), Error> {
        let context = inspect_revision(request)?;
        self.validate_binding(&context)?;
        if &inspect(request)? != self {
            return Err(bad("translation differs from native source computation"));
        }
        Ok(())
    }
}
impl Outputs {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        collagene_model::bounds::model_preflight(&[&self.protein])?;
        self.protein.validate()?;
        self.report
            .validate_binding(&self.report.context_revision)?;
        if !hex(&self.context_revision)
            || self.context_revision != expected
            || !name_valid(&self.protein.name)
            || self.protein.sequence.alphabet() != Alphabet::Protein
            || self.protein.sequence.is_empty()
            || self.protein.sequence.as_str() != self.report.protein
            || self.protein.topology != Topology::Linear
            || !self.protein.features.is_empty()
            || self.protein.original_genbank.is_some()
            || self.protein.metadata.len() != 5
            || self.context_revision
                != crate::cloning::assembly::recipe::hash(&(
                    ENGINE_VERSION,
                    "outputs",
                    &self.report.context_revision,
                    &self.protein.name,
                ))?
        {
            return Err(bad("translation output fields or context disagree"));
        }
        for (key, value) in [
            ("translation:engine", ENGINE_VERSION),
            (
                "translation:source_revision",
                self.report.source_revision.as_str(),
            ),
            (
                "translation:context_revision",
                self.report.context_revision.as_str(),
            ),
            (
                "translation:output_revision",
                self.context_revision.as_str(),
            ),
            (
                "translation:policy",
                match self.report.selection {
                    Selection::Region { .. } => {
                        "Literal selected-region frame translation; all stops retained"
                    }
                    Selection::Cds { .. } => {
                        "Native CDS translation; one definite terminal stop omitted"
                    }
                },
            ),
        ] {
            if self.protein.metadata.get(key).map(String::as_str) != Some(value) {
                return Err(bad("protein provenance metadata differs"));
            }
        }
        if let Some(feature) = &self.annotation {
            collagene_model::bounds::feature_preflight(&[std::slice::from_ref(feature)])?;
        }
        if self.annotation != annotation(&self.report, &self.protein.name) {
            return Err(bad(
                "reviewable annotation differs from translation proposal",
            ));
        }
        Ok(())
    }
    pub fn validate_sources(&self, request: &Request, name: &str) -> Result<(), Error> {
        self.validate_binding(&output_revision(request, name)?)?;
        if &prepare(request, name)? != self {
            return Err(bad("outputs differ from native source reconstruction"));
        }
        Ok(())
    }
}
