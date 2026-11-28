//! Named verifier predicates and documentary evidence, never imported authority.
use crate::planning::plan::Source;
use collagene_model::{Error, Location};
use serde::{Deserialize, Serialize};

pub const MAX_FINDINGS: usize = 8192;
pub const MAX_FINDING_BYTES: usize = 4 * 1024 * 1024;
pub const RULE_SPEC: &str = "collagene://verification-rules/v1";

macro_rules! rules {
    ($($name:ident),+ $(,)?) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        pub enum RuleId { $($name),+ }
        impl RuleId {
            pub const ALL: &'static [Self] = &[$(Self::$name),+];
            pub fn id(self) -> &'static str {
                match self { $(Self::$name => stringify!($name)),+ }
            }
        }
    };
}
rules! {
    Admission, PlanShape, TargetBinding, InputIdentity, InputShape, InputRevision,
    SourceDependency, SourceChemistry, AnnotationProjection, PcrPrimer, PcrBinding,
    PcrSelection, MismatchPcrPolicy, OverlapGeometry, HomologyPolicy,
    HomologyReconstruction, RestrictionProfile, RestrictionCleavage,
    ReactionConditions, DigestSelection, EndRepair, EndTreatment, TailAddition,
    OligoAnnealing, GoldenGatePolicy, GoldenGateAssembly, LigationGeometry,
    LigationChemistry, LigationAlternatives, RecombinationReference, TopoDeclaration,
    TopoReconstruction, TopoSelection, TopoTargetPolicy, SavedMaterial,
    MaterialOrigin, ProductIntegrity, ProductBudget, EvidenceBudget, CheckBudget,
    ReferenceOnlyAuthority, TargetAgreement, StepCompletion, DependencyBlocked,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum CheckVerdict {
    Pass,
    Fail,
    Indeterminate,
    Informational,
    Blocked,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum CheckScope {
    Plan,
    Input { input_id: String },
    Step { step_index: usize },
    Target { step_index: usize },
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CheckCitation {
    pub source: String,
    pub locator: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CheckFinding {
    pub rule: RuleId,
    pub verdict: CheckVerdict,
    pub scope: CheckScope,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, 1000>")]
    pub subjects: Vec<Source>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Location, 64>")]
    pub coordinates: Vec<Location>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, CheckCitation, 8>")]
    pub citations: Vec<CheckCitation>,
    pub message: String,
}
impl CheckFinding {
    /// The named predicate is supplied by its guard, not inferred from prose.
    /// Every predicate carries its versioned internal specification locator;
    /// scientific/profile references can be appended at the actual domain site.
    pub fn new(
        rule: RuleId,
        verdict: CheckVerdict,
        scope: CheckScope,
        message: impl Into<String>,
    ) -> Self {
        Self {
            rule,
            verdict,
            scope,
            subjects: vec![],
            coordinates: vec![],
            citations: vec![CheckCitation {
                source: RULE_SPEC.into(),
                locator: rule.id().into(),
            }],
            message: message.into(),
        }
    }
    pub fn with_subjects(mut self, subjects: Vec<Source>) -> Self {
        self.subjects = subjects;
        self
    }
    pub fn with_coordinates(mut self, coordinates: Vec<Location>) -> Self {
        self.coordinates = coordinates;
        self
    }
    pub fn with_citations(mut self, citations: Vec<CheckCitation>) -> Self {
        self.citations.extend(citations);
        self
    }
    /// Plain-text presentation; scientific records retain the complete typed finding.
    pub fn summary(&self) -> String {
        format!(
            "{:?} · {:?} · {:?}: {}",
            self.rule, self.verdict, self.scope, self.message
        )
    }
    pub fn validate(&self) -> Result<(), Error> {
        let bounded_text = |value: &str, limit: usize| {
            value.len() <= limit && !value.chars().any(char::is_control)
        };
        if self.message.trim().is_empty()
            || self.message.len() > 4096
            || self.subjects.len() > 1000
            || self.coordinates.len() > 64
            || self.citations.is_empty()
            || self.citations.len() > 8
        {
            return Err(Error::TooLarge);
        }
        if matches!(&self.scope, CheckScope::Input{input_id} if !bounded_text(input_id,128))
            || matches!(&self.scope, CheckScope::Step{step_index}|CheckScope::Target{step_index} if *step_index>=1000)
            || self
                .subjects
                .iter()
                .any(|source| matches!(source,Source::Input(id) if !bounded_text(id,128)))
            || !self
                .citations
                .iter()
                .any(|citation| citation.source == RULE_SPEC && citation.locator == self.rule.id())
            || self.citations.iter().any(|citation| {
                !(citation.source == RULE_SPEC || citation.source.starts_with("https://"))
                    || !bounded_text(&citation.source, 2048)
                    || citation.locator.trim().is_empty()
                    || !bounded_text(&citation.locator, 512)
            })
        {
            return Err(Error::Format(
                "invalid named verification check evidence".into(),
            ));
        }
        // This checks encoding/resource shape only. A coordinate's applicability
        // to a particular source is established at the scientific guard.
        for location in &self.coordinates {
            collagene_model::bounds::ModelBudget::default().location(location)?;
            crate::encoding::bounded_json::preflight(location, 16 * 1024)?;
        }
        Ok(())
    }
}

impl RuleId {
    pub fn predicate(self) -> &'static str {
        match self {
            Self::Admission => {
                "The bounded request and nested source representation can be admitted before hashing or reconstruction."
            }
            Self::PlanShape => {
                "The plan has an admissible nonempty step graph and bounded operation fields."
            }
            Self::TargetBinding => {
                "The separately supplied intended target matches its explicitly bound revision."
            }
            Self::InputIdentity => {
                "Each plan-local input identifier is unambiguous and admissible."
            }
            Self::InputShape => {
                "The referenced source molecule satisfies its native sequence and annotation invariants."
            }
            Self::InputRevision => {
                "The source snapshot equals the revision explicitly bound by the input."
            }
            Self::SourceDependency => {
                "Every operation source names an existing input or a strictly earlier plan step."
            }
            Self::SourceChemistry => {
                "The operation has the required freshly reconstructed strand, pairing and terminal state."
            }
            Self::AnnotationProjection => {
                "Annotations are projected from actual source coordinates with unresolved mappings retained."
            }
            Self::PcrPrimer => {
                "Each supplied primer satisfies the explicit DNA, domain, tail and orientation admission rules."
            }
            Self::PcrBinding => {
                "The canonical declared PCR model finds supported binding geometry on the exact source."
            }
            Self::PcrSelection => {
                "The selected PCR geometry exists and has the authority required for the next operation."
            }
            Self::MismatchPcrPolicy => {
                "Intentional mismatches satisfy the separately declared mismatch and three-prime policy."
            }
            Self::OverlapGeometry => {
                "Declared overlap lengths, order, orientation and actual paired sequences agree."
            }
            Self::HomologyPolicy => {
                "The selected homology method and its explicitly declared recommendations apply."
            }
            Self::HomologyReconstruction => {
                "Actual source strands, declared processing and overlaps reconstruct the reported product."
            }
            Self::RestrictionProfile => {
                "The exact enzyme geometry or catalogue selection matches the current admitted profile."
            }
            Self::RestrictionCleavage => {
                "Fresh positional and cleavage evidence supports the reported restriction products under the declared model."
            }
            Self::ReactionConditions => {
                "Documentary formulation, buffer and temperature evidence matches the explicitly selected conditions."
            }
            Self::DigestSelection => {
                "The selected digest fragment exists in the actual bounded reconstruction."
            }
            Self::EndRepair => {
                "Only the explicitly requested accessible-end repair is applied to actual strand state."
            }
            Self::EndTreatment => {
                "Only the requested terminal chemistry transform is applied; no other repair is inferred."
            }
            Self::TailAddition => {
                "The declared nucleotide addition applies to the modeled actual endpoints."
            }
            Self::OligoAnnealing => {
                "The declared oligos have supported antiparallel pairing under the admitted exact model."
            }
            Self::GoldenGatePolicy => {
                "The explicit one-copy pool, recut and enumeration limits are admitted."
            }
            Self::GoldenGateAssembly => {
                "Fresh source cleavage and bounded assembly enumeration support the selected candidate."
            }
            Self::LigationGeometry => {
                "The intended fragment order and orientations form the declared molecular geometry."
            }
            Self::LigationChemistry => {
                "Actual strand end groups and pairing permit the intended ligation under the explicit model."
            }
            Self::LigationAlternatives => {
                "Retained alternative one-copy ligation paths and enumeration limits remain explicit."
            }
            Self::RecombinationReference => {
                "The supplied finite reference-site profile and selected reference product agree."
            }
            Self::TopoDeclaration => {
                "Declared TOPO source roles, strands, pairing and transcript satisfy native admission."
            }
            Self::TopoReconstruction => {
                "The actual declared TOPO graph and ordered actions reconstruct the retained evidence."
            }
            Self::TopoSelection => {
                "The selected TOPO component or projection exists and has the stated representability."
            }
            Self::TopoTargetPolicy => {
                "The explicitly allowed TOPO target comparison applies to the designated final step."
            }
            Self::SavedMaterial => {
                "Saved material is freshly replayed from its exact closed source history."
            }
            Self::MaterialOrigin => {
                "The selected native origin and source event are bound to the reconstructed material."
            }
            Self::ProductIntegrity => {
                "The staged product satisfies native molecule, annotation and provenance invariants."
            }
            Self::ProductBudget => {
                "Accumulated reconstructed products and actual material remain inside the shared bounds."
            }
            Self::EvidenceBudget => {
                "Accumulated retained scientific evidence remains inside the shared byte bounds."
            }
            Self::CheckBudget => {
                "The report retains explicit completion or limit outcomes within its shared check budget."
            }
            Self::ReferenceOnlyAuthority => {
                "Reference-only ancestry is not promoted into physical material authority."
            }
            Self::TargetAgreement => {
                "The designated final reconstructed product matches the independent target under the declared comparison policy."
            }
            Self::StepCompletion => {
                "This original plan step has a staged reconstructed product; this alone is not experimental success."
            }
            Self::DependencyBlocked => {
                "A direct required input or producer is unavailable, so this dependent step was not reconstructed."
            }
        }
    }
}
