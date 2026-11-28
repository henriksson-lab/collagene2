//! Declared research-vector graphs on exact reference and annotation premises.
//! No role inference, biological absence, production-system or safety verdict.
use super::{AnnotationAssumption, ReferenceClaim, ReferenceElement, Verdict};
use collagene_model::{Alphabet, Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub const ENGINE_VERSION: &str = "vector-class-qc-1";
pub const MAX_REQUEST_BYTES: usize = 2 * 1024 * 1024;
pub const MAX_REPORT_BYTES: usize = 4 * 1024 * 1024;
const MAX_BASES: usize = 2_000_000;
pub const STATEMENT: &str = super::STATEMENT;

macro_rules! data {
    ($($item:tt)*) => {
        #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
        #[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
        #[serde(deny_unknown_fields)]
        $($item)*
    };
}
mod assessment;
pub mod profiles;
#[cfg(test)]
mod tests;
pub use assessment::{BindingEvidence, BoundaryMeasurement, Finding, Report, Rule, Span, assess};
pub use profiles::{Profile, ProfileRef, catalogue, catalogue_revision};

fn bad(message: &str) -> Error {
    Error::Unsupported(format!("declared vector-class QC: {message}"))
}
fn text(value: &str, limit: usize) -> bool {
    !value.trim().is_empty() && value.len() <= limit && !value.chars().any(char::is_control)
}
fn bounded(value: &impl Serialize, maximum: usize) -> Result<(), Error> {
    crate::encoding::bounded_json::preflight(value, maximum)
}
fn hash(value: &impl Serialize) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(value)
}

data! {
    pub struct RoleBinding {
        pub id: String,
        pub role_id: String,
        pub feature_index: usize,
        pub reference_id: String,
    }
}
data! {
    pub struct BoundarySelection {
        /// Outer edges of these explicit feature bindings delimit one traversal.
        /// reverse is the declared graph axis, not an inferred transcription path.
        pub start_binding: String,
        pub end_binding: String,
        pub reverse: bool,
    }
}
data! {
    pub struct DeclaredSpacing {
        pub id: String,
        pub upstream_binding: String,
        pub downstream_binding: String,
        pub minimum_gap: usize,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub maximum_gap: Option<usize>,
        pub rationale: String,
    }
}
data! {
    pub struct Request {
        pub source_revision: String,
        pub profile: ProfileRef,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, AnnotationAssumption, 64>")]
        pub annotations: Vec<AnnotationAssumption>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, ReferenceElement, 64>")]
        pub references: Vec<ReferenceElement>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, RoleBinding, 64>")]
        pub bindings: Vec<RoleBinding>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub boundaries: Option<BoundarySelection>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, DeclaredSpacing, 32>")]
        pub spacing: Vec<DeclaredSpacing>,
        /// Optional exhaustive literal-reference requirements, including maximum0
        /// for a declared forbidden literal. No supplied sequence has an inferred role.
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, super::required_parts::RequiredPart, 64>")]
        pub literal_parts: Vec<super::required_parts::RequiredPart>,
        pub literal_part_limits: super::required_parts::RequiredPartLimits,
    }
}
impl Request {
    fn options(&self) -> super::Options {
        super::Options {
            annotations: self.annotations.clone(),
            references: self.references.clone(),
            reference_claims: self
                .bindings
                .iter()
                .map(|binding| ReferenceClaim {
                    feature_index: binding.feature_index,
                    reference_id: binding.reference_id.clone(),
                })
                .collect(),
            required_parts: self.literal_parts.clone(),
            required_part_limits: self.literal_part_limits.clone(),
            ..Default::default()
        }
    }
    fn shape(&self) -> Result<Profile, Error> {
        bounded(self, MAX_REQUEST_BYTES)?;
        let profile = profiles::resolve(&self.profile)?;
        if self.annotations.len() > 64
            || self.references.len() > 64
            || self.bindings.len() > 64
            || self.spacing.len() > 32
            || self.literal_parts.len() > 64
            || self.source_revision.len() != 64
            || !self
                .source_revision
                .bytes()
                .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
        {
            return Err(Error::TooLarge);
        }
        let mut ids = BTreeSet::new();
        let mut role_features = BTreeSet::new();
        for binding in &self.bindings {
            if !text(&binding.id, 128)
                || !ids.insert(binding.id.as_str())
                || !role_features.insert((binding.role_id.as_str(), binding.feature_index))
                || !text(&binding.reference_id, 128)
                || !profile
                    .definition
                    .roles
                    .iter()
                    .any(|role| role.id == binding.role_id)
            {
                return Err(bad("invalid or duplicate declared role binding"));
            }
        }
        if self.boundaries.as_ref().is_some_and(|boundaries| {
            boundaries.start_binding == boundaries.end_binding
                || !ids.contains(boundaries.start_binding.as_str())
                || !ids.contains(boundaries.end_binding.as_str())
        }) {
            return Err(bad("boundaries require two distinct existing bindings"));
        }
        let mut spacing_ids = BTreeSet::new();
        for spacing in &self.spacing {
            if !text(&spacing.id, 128)
                || !spacing_ids.insert(spacing.id.as_str())
                || !ids.contains(spacing.upstream_binding.as_str())
                || !ids.contains(spacing.downstream_binding.as_str())
                || spacing.upstream_binding == spacing.downstream_binding
                || spacing.minimum_gap > MAX_BASES
                || spacing
                    .maximum_gap
                    .is_some_and(|maximum| maximum < spacing.minimum_gap || maximum > MAX_BASES)
                || !text(&spacing.rationale, 1024)
            {
                return Err(bad("invalid explicitly declared spacing constraint"));
            }
        }
        Ok(profile)
    }
    pub fn validate(&self, source: &Molecule) -> Result<(), Error> {
        // Admit recursive source annotation trees before hashing or cloning.
        collagene_model::bounds::model_preflight(&[source])?;
        self.shape()?;
        if source.sequence.alphabet() != Alphabet::Dna {
            return Err(bad("declared DNA vector required"));
        }
        self.options().validate(source)?;
        if self.source_revision != source.revision() {
            return Err(Error::Format("vector-class source revision differs".into()));
        }
        Ok(())
    }
}
pub fn revision(source: &Molecule, request: &Request) -> Result<String, Error> {
    request.validate(source)?;
    hash(&(
        ENGINE_VERSION,
        super::ENGINE_VERSION,
        crate::sequence::search::SEARCH_ENGINE_VERSION,
        source.revision(),
        request,
    ))
}
