//! Finite declared base-editor reference states, not biological outcome prediction.
//! Guide/PAM geometry and codon calls remain owned by the canonical CRISPR and
//! translation engines. No activity/efficiency score or physical repair is inferred.
use super::{GuideCandidate, GuideOptions};
use crate::sequence::translation;
use collagene_model::{Alphabet, Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub const ENGINE_VERSION: &str = "base-editing-1";
pub const MAX_REPORT_BYTES: usize = 2 * 1024 * 1024;
pub const MAX_REQUEST_BYTES: usize = 64 * 1024;
const RULE_SPEC: &str = "collagene://base-editing-rules/v1";
pub const STATEMENT: &str = "Reference states enumerate declared within-window substitutions, including the unedited reference. Finite combinatorial completeness does not establish biological outcome completeness, activity, efficiency, repair chemistry, isolation or expression. Outside-window edits, other substitutions, indels and off-target outcomes are not enumerated.";

macro_rules! dto {
    ($($item:tt)*) => {
        #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        $($item)*
    };
}
mod consequences;
pub use consequences::{CdsEvidence, CdsRequest, CodonConsequence, ConsequenceKind};
#[cfg(test)]
mod tests;

dto! {pub struct Citation {pub source: String, pub locator: String}}
dto! {pub enum Conversion {CytosineTransition, AdenineTransition}}
impl Conversion {
    fn bases(&self) -> (u8, u8) {
        match self {
            Self::CytosineTransition => (b'C', b'T'),
            Self::AdenineTransition => (b'A', b'G'),
        }
    }
}
dto! {
    pub struct EditorProfile {
        pub id: String,
        pub edition: String,
        pub revision: String,
        /// Exact compatible guide/PAM geometry pin; nominal cut marks do not
        /// assert that this editor physically cleaves DNA.
        pub guide_profile_revision: String,
        pub conversion: Conversion,
        /// Inclusive, one-based positions in biological guide orientation.
        /// Position 1 is the guide end distal to a three-prime PAM.
        pub window_first: usize,
        pub window_last: usize,
        pub custom: bool,
        pub rationale: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Citation, 8>")]
        pub citations: Vec<Citation>,
    }
}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("base editing: {message}"))
}
fn text(value: &str, maximum: usize) -> bool {
    !value.trim().is_empty() && value.len() <= maximum && !value.chars().any(char::is_control)
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn bounded(value: &impl Serialize, maximum: usize) -> Result<(), Error> {
    crate::encoding::bounded_json::preflight(value, maximum)
}
fn hash(value: &impl Serialize) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(value)
}
fn citation_valid(citation: &Citation) -> bool {
    text(&citation.source, 2048)
        && citation.source.starts_with("https://")
        && text(&citation.locator, 512)
}
impl EditorProfile {
    fn shape(&self) -> Result<(), Error> {
        if !text(&self.id, 128)
            || !text(&self.edition, 64)
            || !hex(&self.guide_profile_revision)
            || !text(&self.rationale, 2048)
            || self.window_first == 0
            || self.window_first > self.window_last
            || self.window_last > 64
            || self.window_last - self.window_first + 1 > 12
            || self.citations.is_empty()
            || self.citations.len() > 8
            || self
                .citations
                .iter()
                .any(|citation| !citation_valid(citation))
        {
            return Err(bad("invalid bounded declared editor profile"));
        }
        Ok(())
    }
    pub fn compute_revision(&self) -> Result<String, Error> {
        self.shape()?;
        hash(&(
            ENGINE_VERSION,
            &self.id,
            &self.edition,
            &self.guide_profile_revision,
            &self.conversion,
            self.window_first,
            self.window_last,
            self.custom,
            &self.rationale,
            &self.citations,
        ))
    }
    pub fn validate(&self) -> Result<(), Error> {
        self.shape()?;
        if !hex(&self.revision) || self.revision != self.compute_revision()? {
            return Err(bad("editor profile content differs from its pin"));
        }
        if (!self.custom && *self != be1_2016()?)
            || (self.custom && self.id == "komor-2016-be1-c4-c8")
        {
            return Err(bad(
                "unknown/altered builtin profile; custom declarations require a distinct identity",
            ));
        }
        Ok(())
    }
}
/// The primary study describes a typical BE1 window at PAM-distal positions4–8.
/// This nominal finite C→T reference-state model does not model its C→U chemical
/// intermediate, repair, rates or universal/outside-window exclusion.
pub fn be1_2016() -> Result<EditorProfile, Error> {
    let mut profile = EditorProfile {
        id: "komor-2016-be1-c4-c8".into(), edition: "1".into(), revision: String::new(),
        guide_profile_revision: super::profile_revision(&super::spcas9())?,
        conversion: Conversion::CytosineTransition, window_first: 4, window_last: 8,
        custom: false,
        rationale: "Declared typical BE1 window, not an absolute boundary of biological activity. Enumerate conditional resolved C:G→T:A reference substitutions; do not simulate C→U chemistry, repair or experimental yield.".into(),
        citations: vec![Citation {
            source: "https://doi.org/10.1038/nature17946".into(),
            locator: "Komor et al. (2016), Nature 533:420–424, Results and Fig.1–2: BE1 window typically protospacer positions4–8 counted from PAM-distal position1; C→U intermediate and resolved C→T/G→A substitutions.".into(),
        }],
    };
    profile.revision = profile.compute_revision()?;
    Ok(profile)
}
/// The current admitted builtins are owned here, not duplicated in transports.
pub fn catalogue() -> Result<Vec<EditorProfile>, Error> {
    Ok(vec![be1_2016()?])
}
pub fn catalogue_revision() -> Result<String, Error> {
    hash(&(ENGINE_VERSION, "profiles", catalogue()?))
}
dto! {
    pub struct Request {
        pub guide_options: GuideOptions,
        pub guide: GuideCandidate,
        pub profile: EditorProfile,
        /// One-based guide positions explicitly intended to change. Empty means
        /// no intent declaration; such outcomes do not claim intent satisfaction.
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 12>")]
        pub intended_positions: Vec<usize>,
        pub max_outcomes: usize,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub cds: Option<CdsRequest>,
    }
}
impl Request {
    pub fn validate(&self, molecule: &Molecule) -> Result<(), Error> {
        bounded(self, MAX_REQUEST_BYTES)?;
        self.profile.validate()?;
        if molecule.sequence.alphabet() != Alphabet::Dna
            || !(1..=4096).contains(&self.max_outcomes)
            || self.intended_positions.len() > 12
            || self.profile.window_last > self.guide_options.profile.guide_length
            || self.guide_options.profile.pam_side != super::PamSide::ThreePrime
            || self.profile.guide_profile_revision
                != super::profile_revision(&self.guide_options.profile)?
        {
            return Err(bad(
                "DNA, compatible pinned three-prime-PAM geometry and bounded outcomes are required",
            ));
        }
        let mut intended = BTreeSet::new();
        if self.intended_positions.iter().any(|position| {
            *position == 0
                || *position > self.guide_options.profile.guide_length
                || !intended.insert(*position)
        }) {
            return Err(bad(
                "intended positions must be distinct existing one-based guide positions",
            ));
        }
        super::validate_candidate(molecule, &self.guide_options, &self.guide)?;
        if let Some(cds) = &self.cds {
            cds.validate(molecule)?;
        }
        Ok(())
    }
}
dto! {pub enum Rule {GuideBinding, WindowDeclaration, IntendedEdits, FiniteEnumeration, CodonConsequences}}
dto! {pub enum Verdict {Pass, Fail, Indeterminate, Informational}}
impl Rule {
    fn locator(&self) -> String {
        format!(
            "{self:?}: declared finite reference-state predicate; no biological probability or exhaustiveness claim"
        )
    }
}

dto! {
    pub struct Finding {
        pub rule: Rule,
        pub verdict: Verdict,
        pub message: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Citation, 10>")]
        pub citations: Vec<Citation>,
    }
}
fn finding(
    rule: Rule,
    verdict: Verdict,
    message: impl Into<String>,
    profile: &EditorProfile,
) -> Finding {
    let mut citations = vec![Citation {
        source: RULE_SPEC.into(),
        locator: rule.locator(),
    }];
    citations.extend(profile.citations.iter().cloned());
    Finding {
        rule,
        verdict,
        message: message.into(),
        citations,
    }
}
dto! {
    pub struct EditableSite {
        pub guide_position: usize,
        pub reference_position: usize,
        pub oriented_from: char,
        pub oriented_to: char,
        pub reference_from: char,
        pub reference_to: char,
        pub intended: bool,
    }
}
dto! {
    pub struct ReferenceState {
        /// Changed guide positions, not source-forward coordinates or site indices.
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 12>")]
        pub edited_positions: Vec<usize>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 12>")]
        pub bystander_positions: Vec<usize>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub satisfies_intent: Option<bool>,
        pub reference_guide: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, CodonConsequence, 12>")]
        pub codon_consequences: Vec<CodonConsequence>,
    }
}
dto! {
    pub struct Report {
        pub engine_version: String,
        pub context_revision: String,
        pub source_revision: String,
        pub source_length: usize,
        pub source_topology: Topology,
        pub request: Request,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, EditableSite, 12>")]
        pub sites: Vec<EditableSite>,
        pub finite_state_count: usize,
        pub finite_enumeration_complete: bool,
        /// Always false: this is not an experimental outcome enumeration.
        pub biological_outcomes_exhaustive: bool,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, ReferenceState, 4096>")]
        pub states: Vec<ReferenceState>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub cds: Option<CdsEvidence>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Finding, 8>")]
        pub findings: Vec<Finding>,
        pub statement: String,
        pub report_revision: String,
    }
}
pub fn revision(molecule: &Molecule, request: &Request) -> Result<String, Error> {
    request.validate(molecule)?;
    hash(&(
        ENGINE_VERSION,
        super::ENGINE_VERSION,
        molecule.revision(),
        request,
    ))
}
pub fn assess(molecule: &Molecule, request: &Request) -> Result<Report, Error> {
    let context_revision = revision(molecule, request)?;
    let (from, to) = request.profile.conversion.bases();
    let mut sites = Vec::new();
    for guide_position in request.profile.window_first..=request.profile.window_last {
        if request.guide.guide_sequence.as_bytes()[guide_position - 1] != from {
            continue;
        }
        let reference_position = request.guide.guide_positions[guide_position - 1];
        let reference_from = if request.guide.reverse {
            collagene_model::complement(from, Alphabet::Dna)
        } else {
            from
        };
        let reference_to = if request.guide.reverse {
            collagene_model::complement(to, Alphabet::Dna)
        } else {
            to
        };
        sites.push(EditableSite {
            guide_position,
            reference_position,
            oriented_from: from as char,
            oriented_to: to as char,
            reference_from: reference_from as char,
            reference_to: reference_to as char,
            intended: request.intended_positions.contains(&guide_position),
        });
    }
    let (cds, codons) = consequences::prepare(molecule, request.cds.as_ref(), &sites)?;
    let finite_state_count = 1usize << sites.len();
    let mut report = Report {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: molecule.revision(),
        source_length: molecule.sequence.len(),
        source_topology: molecule.topology,
        request: request.clone(),
        sites,
        finite_state_count,
        finite_enumeration_complete: true,
        biological_outcomes_exhaustive: false,
        states: vec![],
        cds,
        findings: vec![],
        statement: STATEMENT.into(),
        report_revision: String::new(),
    };
    let eligible_intent = request.intended_positions.iter().all(|position| {
        report
            .sites
            .iter()
            .any(|site| site.guide_position == *position)
    });
    report.findings.push(finding(Rule::GuideBinding, Verdict::Pass,
        "The exact selected guide/PAM was freshly reconstructed from the current source and compatible pinned geometry; this is reference binding geometry, not activity.", &request.profile));
    report.findings.push(finding(Rule::WindowDeclaration, Verdict::Informational,
        "Window positions and conversions are declared model assumptions. The source citation does not make them universal biological boundaries.", &request.profile));
    report.findings.push(finding(Rule::IntendedEdits,
        if request.intended_positions.is_empty() {Verdict::Indeterminate} else if eligible_intent {Verdict::Pass} else {Verdict::Fail},
        if request.intended_positions.is_empty() {"No intended edit was declared; states do not claim an intended path."}
        else if eligible_intent {"Every requested position has the required concrete convertible base inside this declared window. This establishes finite-model reachability only."}
        else {"At least one intended position is not convertible inside this declared window. This model cannot produce that requested path; activity outside the model is not excluded."}, &request.profile));
    let mut bytes = crate::verification::evidence::evidence_bytes(&report, MAX_REPORT_BYTES)?;
    for mask in 0..finite_state_count {
        if report.states.len() >= request.max_outcomes {
            report.finite_enumeration_complete = false;
            break;
        }
        let mut reference_guide = request.guide.guide_sequence.as_bytes().to_vec();
        let selected = report
            .sites
            .iter()
            .enumerate()
            .filter(|(index, _)| mask & (1usize << index) != 0)
            .map(|(_, site)| site)
            .collect::<Vec<_>>();
        for site in &selected {
            reference_guide[site.guide_position - 1] = site.oriented_to as u8;
        }
        let state = ReferenceState {
            edited_positions: selected.iter().map(|site| site.guide_position).collect(),
            bystander_positions: if request.intended_positions.is_empty() {
                vec![]
            } else {
                selected
                    .iter()
                    .filter(|site| !site.intended)
                    .map(|site| site.guide_position)
                    .collect()
            },
            satisfies_intent: (!request.intended_positions.is_empty()).then(|| {
                request
                    .intended_positions
                    .iter()
                    .all(|position| selected.iter().any(|site| site.guide_position == *position))
            }),
            reference_guide: String::from_utf8(reference_guide)
                .map_err(|_| bad("non-ASCII reference guide"))?,
            codon_consequences: codons
                .as_ref()
                .map(|codons| codons.assess(&selected))
                .transpose()?
                .unwrap_or_default(),
        };
        let state_bytes = crate::verification::evidence::evidence_bytes(&state, MAX_REPORT_BYTES)?;
        // Reserve summary/codon-scope findings and final hashes before retaining data.
        if bytes.saturating_add(state_bytes + 1) > MAX_REPORT_BYTES - 128 * 1024 {
            report.finite_enumeration_complete = false;
            break;
        }
        bytes += state_bytes + 1;
        report.states.push(state);
    }
    report.findings.push(finding(Rule::FiniteEnumeration,
        if report.finite_enumeration_complete {Verdict::Informational} else {Verdict::Indeterminate},
        format!("Retained {} of {} exact subsets, including the unedited state when present. {} This finite combinatorial set is not an exhaustive list of biological outcomes.", report.states.len(), finite_state_count,
            if report.finite_enumeration_complete {"Declared finite enumeration completed."} else {"The explicit outcome/report bound stopped enumeration; omitted subsets are unresolved."}), &request.profile));
    if let Some(cds) = &report.cds {
        let unresolved = report.states.iter().any(|state| {
            state
                .codon_consequences
                .iter()
                .any(|codon| codon.kind == ConsequenceKind::Unresolved)
        });
        let mut check = finding(
            Rule::CodonConsequences,
            if cds.available && !unresolved {
                Verdict::Informational
            } else {
                Verdict::Indeterminate
            },
            if unresolved {
                format!(
                    "{} At least one retained ordinary codon call is unresolved.",
                    cds.reason
                )
            } else {
                cds.reason.clone()
            },
            &request.profile,
        );
        check.citations.push(Citation {
            source: "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi".into(),
            locator: "Pinned NCBI genetic-code table, evaluated by the canonical native translation owner; ordinary codon calls only.".into(),
        });
        report.findings.push(check);
    }
    bounded(&report, MAX_REPORT_BYTES)?;
    report.report_revision = report.content_revision()?;
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl Report {
    fn content_revision(&self) -> Result<String, Error> {
        hash(&(
            (
                &self.engine_version,
                &self.context_revision,
                &self.source_revision,
                self.source_length,
                self.source_topology,
            ),
            (
                &self.request,
                &self.sites,
                self.finite_state_count,
                self.finite_enumeration_complete,
                self.biological_outcomes_exhaustive,
            ),
            (&self.states, &self.cds, &self.findings, &self.statement),
        ))
    }
    /// Documentary integrity only. Fresh source/guide replay is required for use.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        bounded(self, MAX_REPORT_BYTES)?;
        bounded(&self.request, MAX_REQUEST_BYTES)?;
        self.request.profile.validate()?;
        self.request.guide_options.validate()?;
        self.request.guide_options.profile.current()?;
        if self.engine_version != ENGINE_VERSION
            || self.statement != STATEMENT
            || !hex(&self.context_revision)
            || self.context_revision != expected
            || !hex(&self.source_revision)
            || self.source_length == 0
            || self.source_length > collagene_model::MAX_SEQUENCE_LEN
            || self.sites.len() > 12
            || self.states.len() > self.request.max_outcomes
            || !(1..=4096).contains(&self.request.max_outcomes)
            || self.finite_state_count != 1usize << self.sites.len()
            || self.states.len() > self.finite_state_count
            || self.finite_enumeration_complete != (self.states.len() == self.finite_state_count)
            || self.biological_outcomes_exhaustive
            || self.findings.len() > 8
            || self.report_revision != self.content_revision()?
        {
            return Err(bad(
                "report bindings, current version, scope or content differ",
            ));
        }
        let mut positions = BTreeSet::new();
        for site in &self.sites {
            let (from, to) = self.request.profile.conversion.bases();
            let reference_from = if self.request.guide.reverse {
                collagene_model::complement(from, Alphabet::Dna)
            } else {
                from
            };
            let reference_to = if self.request.guide.reverse {
                collagene_model::complement(to, Alphabet::Dna)
            } else {
                to
            };
            if site.guide_position < self.request.profile.window_first
                || site.guide_position > self.request.profile.window_last
                || site.guide_position > self.request.guide_options.profile.guide_length
                || site.reference_position >= self.source_length
                || !positions.insert(site.guide_position)
                || site.oriented_from != from as char
                || site.oriented_to != to as char
                || site.reference_from != reference_from as char
                || site.reference_to != reference_to as char
                || self
                    .request
                    .guide
                    .guide_positions
                    .get(site.guide_position - 1)
                    != Some(&site.reference_position)
                || site.intended
                    != self
                        .request
                        .intended_positions
                        .contains(&site.guide_position)
            {
                return Err(bad("invalid editable-site coordinates"));
            }
        }
        for (mask, state) in self.states.iter().enumerate() {
            let expected_positions = self
                .sites
                .iter()
                .enumerate()
                .filter(|(index, _)| mask & (1usize << index) != 0)
                .map(|(_, site)| site.guide_position)
                .collect::<Vec<_>>();
            let expected_bystanders = if self.request.intended_positions.is_empty() {
                vec![]
            } else {
                expected_positions
                    .iter()
                    .filter(|position| !self.request.intended_positions.contains(position))
                    .copied()
                    .collect()
            };
            let satisfies_intent = (!self.request.intended_positions.is_empty()).then(|| {
                self.request
                    .intended_positions
                    .iter()
                    .all(|position| expected_positions.contains(position))
            });
            let mut expected_guide = self.request.guide.guide_sequence.as_bytes().to_vec();
            for site in &self.sites {
                if expected_positions.contains(&site.guide_position) {
                    let Some(base) = expected_guide.get_mut(site.guide_position - 1) else {
                        return Err(bad("missing declared guide call"));
                    };
                    *base = site.oriented_to as u8;
                }
            }
            if state.edited_positions != expected_positions
                || state.bystander_positions != expected_bystanders
                || state.satisfies_intent != satisfies_intent
                || state.codon_consequences.len() > 12
                || state.reference_guide.as_bytes() != expected_guide.as_slice()
                || !state
                    .reference_guide
                    .bytes()
                    .all(|base| b"ACGT".contains(&base))
            {
                return Err(bad("invalid finite reference-state evidence"));
            }
        }
        for finding in &self.findings {
            if !text(&finding.message, 4096)
                || finding.citations.is_empty()
                || finding.citations.len() > 10
                || !finding.citations.iter().any(|citation| {
                    citation.source == RULE_SPEC && citation.locator == finding.rule.locator()
                })
                || finding.citations.iter().any(|citation| {
                    !(citation.source == RULE_SPEC && text(&citation.locator, 512))
                        && !citation_valid(citation)
                })
            {
                return Err(bad("invalid cited editor-rule finding"));
            }
        }
        consequences::validate_report(self)?;
        Ok(())
    }
    pub fn validate_for(&self, molecule: &Molecule, request: &Request) -> Result<(), Error> {
        self.validate_binding(&revision(molecule, request)?)?;
        if *self != assess(molecule, request)? {
            return Err(bad("report differs from fresh declared-model replay"));
        }
        Ok(())
    }
}
