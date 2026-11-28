//! Audited genetic-code translation and bounded ORF candidates, independent of external engines.
//! Coordinates are zero-based half-open; proteins use '*' for stop and 'X' for unresolved ambiguity.
use crate::analysis::base_mask;
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};

const STANDARD: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const MITOCHONDRIAL: &[u8; 64] =
    b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
const MAX_ORFS: usize = 10_000;
const MAX_ORF_AMINO_ACIDS: usize = 8_000_000;
const NO_BOUNDARY: u32 = u32::MAX;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct CodeDescriptor {
    pub id: u8,
    pub name: String,
    pub start_codons: Vec<String>,
}
#[derive(Clone, Copy)]
struct Code {
    id: u8,
    name: &'static str,
    amino_acids: &'static [u8; 64],
    starts: &'static [usize],
}
fn code(id: u8) -> Result<Code, Error> {
    Ok(match id {
        1 => Code {
            id,
            name: "Standard",
            amino_acids: STANDARD,
            starts: &[3, 19, 35],
        },
        2 => Code {
            id,
            name: "Vertebrate mitochondrial",
            amino_acids: MITOCHONDRIAL,
            starts: &[32, 33, 34, 35, 51],
        },
        11 => Code {
            id,
            name: "Bacterial, archaeal and plant plastid",
            amino_acids: STANDARD,
            starts: &[3, 19, 32, 33, 34, 35, 51],
        },
        _ => {
            return Err(Error::Unsupported(format!(
                "genetic code {id} is not in the audited code set (1, 2, 11)"
            )));
        }
    })
}
fn index_codon(index: usize) -> String {
    let bases = b"TCAG";
    [bases[index / 16], bases[index / 4 % 4], bases[index % 4]]
        .into_iter()
        .map(char::from)
        .collect()
}
pub fn genetic_codes() -> Vec<CodeDescriptor> {
    [1, 2, 11]
        .into_iter()
        .map(|id| {
            let c = code(id).unwrap();
            CodeDescriptor {
                id: c.id,
                name: c.name.into(),
                start_codons: c.starts.iter().map(|i| index_codon(*i)).collect(),
            }
        })
        .collect()
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum InitiationPolicy {
    /// Read all codons as internal codons, without biological initiation inference.
    #[default]
    Ordinary,
    /// Recode the first codon as M only when every expansion is a recognized initiator.
    RecognizedStartAsMethionine,
    /// Require a definite recognized first initiator and recode it as M.
    RequireRecognizedStart,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct TranslationOptions {
    pub genetic_code: u8,
    /// Bases omitted at the 5′ end of the chosen orientation, 0/1/2.
    pub frame: u8,
    pub reverse: bool,
    pub initiation: InitiationPolicy,
}
impl Default for TranslationOptions {
    fn default() -> Self {
        Self {
            genetic_code: 1,
            frame: 0,
            reverse: false,
            initiation: InitiationPolicy::Ordinary,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct TranslationReport {
    pub options: TranslationOptions,
    pub code_name: String,
    /// Includes all stop codons; translation does not silently stop at the first one.
    pub protein: String,
    pub nucleotide_length: usize,
    pub translated_codons: usize,
    pub trailing_bases: usize,
    pub ambiguous_codons: usize,
    pub unresolved_codons: usize,
    /// Codons that could be stop or sense; no definite stop is inferred for them.
    pub possible_stop_codons: usize,
    pub initiator_methionine_applied: bool,
    pub warnings: Vec<String>,
}
#[derive(Clone, Copy)]
struct Outcome {
    amino: u8,
    start_all: bool,
    start_any: bool,
    stop_all: bool,
    stop_any: bool,
    ambiguous: bool,
}
struct Decoder {
    code: Code,
    cache: Vec<Option<Outcome>>,
}
impl Decoder {
    fn new(code: Code) -> Self {
        Self {
            code,
            cache: vec![None; 4096],
        }
    }
    fn codon(&mut self, codon: [u8; 3]) -> Outcome {
        let masks = codon.map(base_mask);
        let key = ((masks[0] as usize) << 8) | ((masks[1] as usize) << 4) | masks[2] as usize;
        if let Some(outcome) = self.cache[key] {
            return outcome;
        }
        // All inputs come from a validated nucleotide Sequence, so none of these masks is empty.
        let mut amino = None;
        let mut unique = true;
        let mut starts = 0;
        let mut stops = 0;
        let mut count = 0;
        let bits = [8, 2, 1, 4]; // TCAG indexing of the NCBI code table.
        for a in 0..4 {
            for b in 0..4 {
                for c in 0..4 {
                    if masks[0] & bits[a] == 0 || masks[1] & bits[b] == 0 || masks[2] & bits[c] == 0
                    {
                        continue;
                    }
                    let index = a * 16 + b * 4 + c;
                    let aa = self.code.amino_acids[index];
                    if amino.is_some_and(|previous| previous != aa) {
                        unique = false;
                    }
                    amino.get_or_insert(aa);
                    count += 1;
                    starts += usize::from(self.code.starts.contains(&index));
                    stops += usize::from(aa == b'*');
                }
            }
        }
        let outcome = Outcome {
            amino: if unique { amino.unwrap() } else { b'X' },
            start_all: starts == count,
            start_any: starts > 0,
            stop_all: stops == count,
            stop_any: stops > 0,
            ambiguous: masks.iter().any(|m| m.count_ones() > 1),
        };
        self.cache[key] = Some(outcome);
        outcome
    }
}
fn nucleotide(sequence: &Sequence) -> Result<(), Error> {
    if sequence.alphabet() == Alphabet::Protein {
        Err(Error::Unsupported(
            "nucleotide translation requires DNA or RNA".into(),
        ))
    } else {
        Ok(())
    }
}
fn triple(bytes: &[u8], position: usize, circular: bool) -> [u8; 3] {
    if circular {
        [
            bytes[position % bytes.len()],
            bytes[(position + 1) % bytes.len()],
            bytes[(position + 2) % bytes.len()],
        ]
    } else {
        [bytes[position], bytes[position + 1], bytes[position + 2]]
    }
}

pub fn translate(
    sequence: &Sequence,
    options: &TranslationOptions,
) -> Result<TranslationReport, Error> {
    nucleotide(sequence)?;
    let code = code(options.genetic_code)?;
    if options.frame > 2 {
        return Err(Error::Coordinates(
            "translation frame must be 0, 1 or 2".into(),
        ));
    }
    let oriented = if options.reverse {
        sequence.reverse_complement()?
    } else {
        sequence.clone()
    };
    let bytes = oriented.as_str().as_bytes();
    let start = (options.frame as usize).min(bytes.len());
    let codons = (bytes.len() - start) / 3;
    let trailing = (bytes.len() - start) % 3;
    if codons == 0 && options.initiation == InitiationPolicy::RequireRecognizedStart {
        return Err(Error::Format(
            "a complete first codon is required for initiation".into(),
        ));
    }
    let mut decoder = Decoder::new(code);
    let mut protein = String::with_capacity(codons);
    let (mut ambiguous, mut unresolved, mut possible_stop) = (0, 0, 0);
    let mut recoded = false;
    let mut warnings = vec![];
    for i in 0..codons {
        let outcome = decoder.codon(triple(bytes, start + 3 * i, false));
        let mut aa = outcome.amino;
        if i == 0 && options.initiation != InitiationPolicy::Ordinary {
            if options.initiation == InitiationPolicy::RequireRecognizedStart && !outcome.start_all
            {
                return Err(if outcome.start_any {
                    Error::Unsupported(
                        "first ambiguous codon is not definitely a recognized initiator".into(),
                    )
                } else {
                    Error::Format(
                        "first codon is not a recognized initiator in this genetic code".into(),
                    )
                });
            }
            if outcome.start_all {
                aa = b'M';
                recoded = true;
            } else if outcome.start_any {
                warnings.push("First codon could be an initiator, but ambiguity prevents definite initiation recoding.".into());
            }
        }
        ambiguous += usize::from(outcome.ambiguous);
        unresolved += usize::from(aa == b'X');
        possible_stop += usize::from(outcome.stop_any && !outcome.stop_all);
        protein.push(aa as char);
    }
    if trailing > 0 {
        warnings.push(format!("{trailing} trailing nucleotide(s) do not form a complete codon and were not translated."));
    }
    if codons == 0 {
        warnings.push("This frame contains no complete codon.".into());
    }
    if possible_stop > 0 {
        warnings.push(format!(
            "{possible_stop} ambiguous codon(s) could encode a stop; their amino acid is X."
        ));
    }
    Ok(TranslationReport {
        options: options.clone(),
        code_name: code.name.into(),
        protein,
        nucleotide_length: sequence.len(),
        translated_codons: codons,
        trailing_bases: trailing,
        ambiguous_codons: ambiguous,
        unresolved_codons: unresolved,
        possible_stop_codons: possible_stop,
        initiator_methionine_applied: recoded,
        warnings,
    })
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct CdsTranslation {
    pub source_revision: String,
    pub feature_label: String,
    pub location: Location,
    /// INSDC one-based offset of the first complete codon (1/2/3).
    pub codon_start: u8,
    pub translation: TranslationReport,
    /// CDS protein excludes a definite terminal stop, but preserves internal '*'.
    pub protein: String,
    pub has_terminal_stop: bool,
    pub internal_stops: usize,
    pub complete_under_model: bool,
    pub stored_translation_matches: Option<bool>,
    pub warnings: Vec<String>,
}
fn qualifier<'a>(feature: &'a Feature, key: &str) -> Result<Option<&'a str>, Error> {
    let values: Vec<_> = feature
        .qualifiers
        .iter()
        .filter(|(k, _)| k == key)
        .collect();
    if values.len() > 1 {
        return Err(Error::Format(format!("duplicate CDS /{key} qualifiers")));
    }
    values
        .first()
        .map(|(_, v)| {
            v.as_deref()
                .filter(|v| !v.trim().is_empty())
                .ok_or_else(|| Error::Format(format!("CDS /{key} needs a value")))
        })
        .transpose()
}
fn scalar(feature: &Feature, key: &str, default: u8) -> Result<u8, Error> {
    let Some(value) = qualifier(feature, key)? else {
        return Ok(default);
    };
    let value = value.trim();
    if !value.bytes().all(|b| b.is_ascii_digit()) {
        return Err(Error::Format(format!(
            "CDS /{key} must be a positive integer"
        )));
    }
    value
        .parse::<u8>()
        .ok()
        .filter(|n| *n > 0)
        .ok_or_else(|| Error::Format(format!("CDS /{key} is outside the supported integer range")))
}
// Bound extraction before allocation, and traverse complement joins in biological order.
// A repeated-range CDS may be much larger than its source even when every range is valid.
fn cds_sequence(sequence: &Sequence, location: &Location) -> Result<Sequence, Error> {
    fn length(location: &Location) -> Result<usize, Error> {
        let n = match location {
            Location::Range {
                start,
                end,
                before: false,
                after: false,
            } => end - start,
            Location::Complement(inner) => return length(inner),
            Location::Join(parts) => parts.iter().try_fold(0usize, |n, p| {
                n.checked_add(length(p)?)
                    .filter(|sum| *sum <= collagene_model::MAX_SEQUENCE_LEN)
                    .ok_or(Error::TooLarge)
            })?,
            _ => {
                return Err(Error::Unsupported(
                    "CDS location does not define an exact extractable sequence".into(),
                ));
            }
        };
        if n > collagene_model::MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        Ok(n)
    }
    fn append(sequence: &Sequence, location: &Location, reverse: bool, text: &mut String) {
        match location {
            Location::Range { start, end, .. } => {
                let part = &sequence.as_str()[*start..*end];
                if reverse {
                    text.extend(
                        part.bytes()
                            .rev()
                            .map(|b| collagene_model::complement(b, sequence.alphabet()) as char),
                    );
                } else {
                    text.push_str(part);
                }
            }
            Location::Complement(inner) => append(sequence, inner, !reverse, text),
            Location::Join(parts) => {
                if reverse {
                    for part in parts.iter().rev() {
                        append(sequence, part, true, text);
                    }
                } else {
                    for part in parts {
                        append(sequence, part, false, text);
                    }
                }
            }
            _ => unreachable!("validated exact CDS location"),
        }
    }
    let mut text = String::with_capacity(length(location)?);
    append(sequence, location, false, &mut text);
    Sequence::new(sequence.alphabet(), &text)
}

pub fn translate_cds(molecule: &Molecule, feature: &Feature) -> Result<CdsTranslation, Error> {
    molecule.validate()?;
    translate_cds_validated(molecule, feature, &molecule.revision())
}

/// Internal browser path: caller owns a validated immutable molecule and its revision.
pub(crate) fn translate_cds_validated(
    molecule: &Molecule,
    feature: &Feature,
    source_revision: &str,
) -> Result<CdsTranslation, Error> {
    nucleotide(&molecule.sequence)?;
    if !feature.kind.eq_ignore_ascii_case("CDS") {
        return Err(Error::Format(
            "CDS translation requires a CDS feature".into(),
        ));
    }
    feature.location.validate(molecule.sequence.len())?;
    if feature.needs_review {
        return Err(Error::Unsupported(
            "CDS annotation requires review before feature-directed translation".into(),
        ));
    }
    for (key, _) in &feature.qualifiers {
        if matches!(
            key.as_str(),
            "transl_except"
                | "exception"
                | "ribosomal_slippage"
                | "frameshift"
                | "pseudo"
                | "pseudogene"
                | "trans_splicing"
        ) {
            return Err(Error::Unsupported(format!(
                "CDS /{key} requires biology outside the simple codon model"
            )));
        }
    }
    let codon_start = scalar(feature, "codon_start", 1)?;
    if !(1..=3).contains(&codon_start) {
        return Err(Error::Format("CDS /codon_start must be 1, 2 or 3".into()));
    }
    let genetic_code = scalar(feature, "transl_table", 1)?;
    code(genetic_code)?;
    let sequence = cds_sequence(&molecule.sequence, &feature.location)?;
    let partial = feature.qualifiers.iter().any(|(key, _)| key == "partial");
    let options = TranslationOptions {
        genetic_code,
        frame: codon_start - 1,
        reverse: false,
        initiation: if codon_start == 1 && !partial {
            InitiationPolicy::RecognizedStartAsMethionine
        } else {
            InitiationPolicy::Ordinary
        },
    };
    let translation = translate(&sequence, &options)?;
    let has_terminal_stop = translation.protein.ends_with('*');
    let protein = if has_terminal_stop {
        translation.protein[..translation.protein.len() - 1].to_owned()
    } else {
        translation.protein.clone()
    };
    let internal_stops = protein.bytes().filter(|b| *b == b'*').count();
    let mut warnings = translation.warnings.clone();
    if codon_start > 1 {
        warnings.push(
            "codon_start skips an incomplete 5′ codon; no initiator methionine is inferred.".into(),
        );
    }
    if partial {
        warnings.push("The feature is marked partial; a complete CDS is not asserted.".into());
    }
    if !has_terminal_stop {
        warnings.push("CDS has no definite terminal stop in the selected code.".into());
    }
    if internal_stops > 0 {
        warnings.push(format!("CDS has {internal_stops} internal stop codon(s)."));
    }
    if codon_start == 1 && !partial && !translation.initiator_methionine_applied {
        warnings.push("CDS does not begin with a definite recognized initiator.".into());
    }
    let stored_translation_matches = qualifier(feature, "translation")?.map(|s| {
        s.chars()
            .filter(|c| !c.is_ascii_whitespace())
            .eq(protein.chars())
    });
    if stored_translation_matches == Some(false) {
        warnings.push(
            "Stored /translation differs from the independently computed CDS protein.".into(),
        );
    }
    if translation.unresolved_codons > 0 {
        warnings.push("CDS protein contains unresolved amino acids; complete translation evidence is not asserted.".into());
    }
    let complete_under_model = !partial
        && codon_start == 1
        && translation.initiator_methionine_applied
        && has_terminal_stop
        && internal_stops == 0
        && translation.trailing_bases == 0
        && translation.possible_stop_codons == 0
        && translation.unresolved_codons == 0
        && stored_translation_matches != Some(false);
    Ok(CdsTranslation {
        source_revision: source_revision.into(),
        feature_label: feature.label().into(),
        location: feature.location.clone(),
        codon_start,
        translation,
        protein,
        has_terminal_stop,
        internal_stops,
        complete_under_model,
        stored_translation_matches,
        warnings,
    })
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum StartPolicy {
    #[default]
    AtgOnly,
    CodeStarts,
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum IncompletePolicy {
    #[default]
    CompleteOnly,
    /// Includes 5′/3′ linear-edge candidates; circular candidates still need a definite stop.
    IncludeLinearPartials,
    /// Additionally includes accepted circular starts that reach one turn or an uncertain stop.
    IncludeAllPartials,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct OrfOptions {
    pub genetic_code: u8,
    pub start_policy: StartPolicy,
    pub incomplete: IncompletePolicy,
    pub min_amino_acids: usize,
    pub max_orfs: usize,
    pub both_strands: bool,
}
impl Default for OrfOptions {
    fn default() -> Self {
        Self {
            genetic_code: 1,
            start_policy: StartPolicy::AtgOnly,
            incomplete: IncompletePolicy::CompleteOnly,
            min_amino_acids: 20,
            max_orfs: 1000,
            both_strands: true,
        }
    }
}
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum Termination {
    StopCodon,
    AmbiguousStop,
    LinearBoundary,
    OneTurnLimit,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Orf {
    /// Covers coding nucleotides plus a definite terminal stop codon when present.
    pub location: Location,
    pub reverse: bool,
    /// Offset class of the start in the selected strand orientation, not a global circular frame.
    pub frame: u8,
    pub oriented_start: usize,
    pub nucleotide_length: usize,
    /// Excludes a definite terminal '*'. Accepted initiating codons are represented as M.
    pub protein: String,
    pub crosses_origin: bool,
    pub has_start_codon: bool,
    pub has_stop_codon: bool,
    pub partial_5prime: bool,
    pub partial_3prime: bool,
    pub termination: Termination,
    pub ambiguous_codons: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct OrfReport {
    pub source_revision: String,
    pub options: OrfOptions,
    pub orfs: Vec<Orf>,
    pub truncated: bool,
    pub warnings: Vec<String>,
}
fn accepted_start(codon: [u8; 3], outcome: Outcome, policy: StartPolicy) -> bool {
    match policy {
        StartPolicy::AtgOnly => codon.map(base_mask) == [1, 8, 4],
        StartPolicy::CodeStarts => outcome.start_all,
    }
}
fn location(start: usize, length: usize, n: usize, reverse: bool) -> Location {
    let parts = if start + length <= n {
        vec![(start, start + length)]
    } else {
        vec![(start, n), (0, start + length - n)]
    };
    let parts: Vec<_> = if reverse {
        parts
            .into_iter()
            .rev()
            .map(|(s, e)| Location::range(n - e, n - s))
            .collect()
    } else {
        parts
            .into_iter()
            .map(|(s, e)| Location::range(s, e))
            .collect()
    };
    let inner = if parts.len() == 1 {
        parts.into_iter().next().unwrap()
    } else {
        Location::Join(parts)
    };
    if reverse {
        Location::Complement(Box::new(inner))
    } else {
        inner
    }
}

/// O(n) boundary scan plus bounded output materialization; every codon stays within one source turn.
pub fn find_orfs(molecule: &Molecule, options: &OrfOptions) -> Result<OrfReport, Error> {
    molecule.validate()?;
    find_orfs_validated(molecule, options, &molecule.revision())
}

/// Internal browser path: caller owns a validated immutable molecule and its revision.
pub(crate) fn find_orfs_validated(
    molecule: &Molecule,
    options: &OrfOptions,
    source_revision: &str,
) -> Result<OrfReport, Error> {
    nucleotide(&molecule.sequence)?;
    let code = code(options.genetic_code)?;
    if options.max_orfs == 0
        || options.max_orfs > MAX_ORFS
        || options.min_amino_acids > collagene_model::MAX_SEQUENCE_LEN / 3
    {
        return Err(Error::Format(
            "ORF limits require 1–10,000 results and a feasible minimum length".into(),
        ));
    }
    let mut report=OrfReport{source_revision:source_revision.into(),options:options.clone(),orfs:vec![],truncated:false,warnings:vec!["ORFs are codon-model candidates, not evidence of expression. Ambiguous possible-stop codons form uncertain boundaries. Circular candidates cover at most one source turn.".into()]};
    let n = molecule.sequence.len();
    if n < 3 {
        report
            .warnings
            .push("Molecule is shorter than one complete codon.".into());
        return Ok(report);
    }
    let circular = molecule.topology == Topology::Circular;
    let mut total_amino = 0usize;
    for reverse in [false, true]
        .into_iter()
        .take(if options.both_strands { 2 } else { 1 })
    {
        let sequence = if reverse {
            molecule.sequence.reverse_complement()?
        } else {
            molecule.sequence.clone()
        };
        let bytes = sequence.as_str().as_bytes();
        let mut decoder = Decoder::new(code);
        let mut next = vec![NO_BOUNDARY; n];
        let mut latest = [NO_BOUNDARY; 3];
        // Rolling per-frame successors avoid allocating a doubled circular genome.
        let upper = if circular { 2 * n + 2 } else { n };
        for p in (0..upper).rev() {
            if !circular && p + 3 > n {
                continue;
            }
            let outcome = decoder.codon(triple(bytes, p, circular));
            if outcome.stop_any {
                latest[p % 3] = p as u32;
            }
            if p < n {
                next[p] = latest[p % 3];
            }
        }
        for p in 0..n {
            if !circular && p + 3 > n {
                continue;
            }
            let codon = triple(bytes, p, circular);
            let first = decoder.codon(codon);
            let start = accepted_start(codon, first, options.start_policy);
            let partial_prefix = !circular
                && p < 3
                && options.incomplete != IncompletePolicy::CompleteOnly
                && !start;
            if (!start && !partial_prefix) || first.stop_any {
                continue;
            }
            let max_codons = if circular { n / 3 } else { (n - p) / 3 };
            let limit = p + max_codons * 3;
            let boundary = next[p];
            let boundary = if boundary != NO_BOUNDARY && (boundary as usize) + 3 <= limit {
                Some(boundary as usize)
            } else {
                None
            };
            let (length, has_stop, termination) = if let Some(stop) = boundary {
                let ending = decoder.codon(triple(bytes, stop, circular));
                if ending.stop_all {
                    (stop - p + 3, true, Termination::StopCodon)
                } else {
                    (stop - p, false, Termination::AmbiguousStop)
                }
            } else {
                (
                    max_codons * 3,
                    false,
                    if circular {
                        Termination::OneTurnLimit
                    } else {
                        Termination::LinearBoundary
                    },
                )
            };
            if !has_stop
                && (options.incomplete == IncompletePolicy::CompleteOnly
                    || circular && options.incomplete != IncompletePolicy::IncludeAllPartials)
            {
                continue;
            }
            let amino = length / 3 - usize::from(has_stop);
            if amino < options.min_amino_acids || amino == 0 {
                continue;
            }
            if report.orfs.len() >= options.max_orfs
                || total_amino.saturating_add(amino) > MAX_ORF_AMINO_ACIDS
            {
                report.truncated = true;
                report.warnings.push("ORF result/materialization limit reached; refine minimum length or search a smaller region. Search results are not exhaustive.".into());
                return Ok(report);
            }
            let mut protein = String::with_capacity(amino);
            let mut ambiguous = 0;
            for i in 0..amino {
                let outcome = decoder.codon(triple(bytes, p + i * 3, circular));
                ambiguous += usize::from(outcome.ambiguous);
                protein.push(if i == 0 && start { b'M' } else { outcome.amino } as char);
            }
            total_amino += amino;
            report.orfs.push(Orf {
                location: location(p, length, n, reverse),
                reverse,
                frame: (p % 3) as u8,
                oriented_start: p,
                nucleotide_length: length,
                protein,
                crosses_origin: p + length > n,
                has_start_codon: start,
                has_stop_codon: has_stop,
                partial_5prime: !start,
                partial_3prime: !has_stop,
                termination,
                ambiguous_codons: ambiguous,
            });
        }
    }
    Ok(report)
}
