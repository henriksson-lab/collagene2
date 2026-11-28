//! Current installed standard membership, exact source intervals and report binding.
//! The independent grammar crate checks part format; it does not simulate assembly.
use crate::{cloning::restriction::catalog::builtin_catalog, encoding::bounded_json};
use collagene_assembly_standards::{
    assessment::{Assessment, Limits, Part},
    grammar::Successors,
    profile::{BUILTIN_IDS, CatalogPin, Definition, Profile},
};
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};

pub const ENGINE_VERSION: &str = "assembly-standard-review-1";
pub const MAX_REPORT_BYTES: usize = 2 * 1024 * 1024;
const MAX_REQUEST_BYTES: usize = 16 * 1024 * 1024;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SourcePart {
    pub id: String,
    pub molecule: Molecule,
    pub start: usize,
    pub length: usize,
    pub reverse: bool,
    pub slot: String,
    pub body_start: usize,
    pub body_end: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub profile_id: String,
    pub profile_identity: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, SourcePart, 128>")]
    pub parts: Vec<SourcePart>,
    pub limits: Limits,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProfileSummary {
    pub profile_id: String,
    pub profile_identity: String,
    pub name: String,
    pub scope: String,
    pub definition: Definition,
    pub catalog_pin: CatalogPin,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SourceBinding {
    pub id: String,
    pub revision: String,
    pub start: usize,
    pub length: usize,
    pub reverse: bool,
    pub source_length: usize,
    pub topology: Topology,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine_version: String,
    pub context_revision: String,
    pub profile: ProfileSummary,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, SourceBinding, 128>")]
    pub sources: Vec<SourceBinding>,
    pub assessment: Assessment,
    pub successors: Successors,
    pub report_revision: String,
}

fn bad(s: &str) -> Error {
    Error::Unsupported(format!("assembly standards: {s}"))
}
fn lib_error(e: collagene_assembly_standards::Error) -> Error {
    match e {
        collagene_assembly_standards::Error::Limit(_) => Error::TooLarge,
        _ => bad(&e.to_string()),
    }
}
fn summary(profile: &Profile) -> ProfileSummary {
    ProfileSummary {
        profile_id: profile.definition().id.clone(),
        profile_identity: profile.identity().into(),
        name: profile.definition().name.clone(),
        scope: profile.definition().scope.clone(),
        definition: profile.definition().clone(),
        catalog_pin: profile.catalog_pin().clone(),
    }
}
pub fn catalogue() -> Result<Vec<ProfileSummary>, Error> {
    let c = builtin_catalog()?;
    BUILTIN_IDS
        .iter()
        .map(|id| {
            Profile::builtin(id, &c)
                .map(|profile| summary(&profile))
                .map_err(lib_error)
        })
        .collect()
}
impl SourceBinding {
    pub fn validate(&self) -> Result<(), Error> {
        if self.source_length == 0
            || self.length == 0
            || self.length > self.source_length
            || self.start >= self.source_length
            || self.source_length > collagene_model::MAX_SEQUENCE_LEN
            || (self.topology == Topology::Linear && self.length > self.source_length - self.start)
            || self.id.is_empty()
            || self.id.len() > 128
            || self.id.chars().any(char::is_control)
            || self.revision.len() != 64
            || !self
                .revision
                .bytes()
                .all(|b| b.is_ascii_digit() || matches!(b, b'a'..=b'f'))
        {
            return Err(bad("invalid source interval binding"));
        }
        Ok(())
    }
    /// Convert one oriented part interval to original ascending source spans.
    pub fn map_interval(&self, start: usize, end: usize) -> Result<Vec<(usize, usize)>, Error> {
        self.validate()?;
        if start >= end || end > self.length {
            return Err(bad("finding outside declared part"));
        }
        let offset = if self.reverse {
            self.length - end
        } else {
            start
        };
        let at = (self.start + offset) % self.source_length;
        let length = end - start;
        Ok(if length <= self.source_length - at {
            vec![(at, at + length)]
        } else {
            vec![
                (at, self.source_length),
                (0, length - (self.source_length - at)),
            ]
        })
    }
}
impl SourcePart {
    fn binding(&self) -> SourceBinding {
        SourceBinding {
            id: self.id.clone(),
            revision: self.molecule.revision(),
            start: self.start,
            length: self.length,
            reverse: self.reverse,
            source_length: self.molecule.sequence.len(),
            topology: self.molecule.topology,
        }
    }
    fn validate(&self) -> Result<(), Error> {
        if self.id.is_empty() || self.id.len() > 128 || self.id.chars().any(char::is_control) {
            return Err(bad("invalid part identifier"));
        }
        if self.molecule.sequence.alphabet() != Alphabet::Dna || self.length > 2_000_000 {
            return Err(bad("bounded DNA part required"));
        }
        self.binding().validate()?;
        if self.body_start >= self.body_end
            || self.body_end > self.length
            || self.slot.is_empty()
            || self.slot.len() > 96
            || self.slot.chars().any(char::is_control)
        {
            return Err(bad("invalid body or slot"));
        }
        Ok(())
    }
    fn part(&self) -> Result<Part, Error> {
        self.validate()?;
        let n = self.molecule.sequence.len();
        let source = self.molecule.sequence.as_str().as_bytes();
        let bytes = (0..self.length)
            .map(|i| source[(self.start + i) % n])
            .collect::<Vec<_>>();
        let mut sequence =
            Sequence::dna(std::str::from_utf8(&bytes).map_err(|_| bad("invalid DNA"))?)?;
        if self.reverse {
            sequence = sequence.reverse_complement()?;
        }
        Ok(Part {
            id: self.id.clone(),
            slot: self.slot.clone(),
            sequence,
            body_start: self.body_start,
            body_end: self.body_end,
        })
    }
}
impl Request {
    pub fn validate(&self) -> Result<(), Error> {
        if self.parts.is_empty()
            || self.parts.len() > 128
            || self.profile_id.len() > 128
            || self.profile_identity.len() != 64
        {
            return Err(bad("invalid profile or part count"));
        }
        collagene_model::bounds::model_preflight(
            &self.parts.iter().map(|p| &p.molecule).collect::<Vec<_>>(),
        )?;
        bounded_json::preflight(self, MAX_REQUEST_BYTES)?;
        self.limits.validate().map_err(lib_error)?;
        let c = builtin_catalog()?;
        let profile = Profile::builtin(&self.profile_id, &c).map_err(lib_error)?;
        if profile.identity() != self.profile_identity {
            return Err(bad("current profile identity changed"));
        }
        let mut bases = 0usize;
        let mut ids = std::collections::BTreeSet::new();
        for p in &self.parts {
            p.molecule.validate()?;
            p.validate()?;
            bases = bases.checked_add(p.length).ok_or(Error::TooLarge)?;
            if bases > self.limits.max_bases || !ids.insert(&p.id) {
                return Err(bad("duplicate part or aggregate base bound"));
            }
        }
        if self.parts.len() > self.limits.max_parts {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
pub fn revision(request: &Request) -> Result<String, Error> {
    request.validate()?;
    bounded_json::hash(&(ENGINE_VERSION, request))
}
pub fn assess(request: &Request) -> Result<Report, Error> {
    let context_revision = revision(request)?;
    let c = builtin_catalog()?;
    let profile = Profile::builtin(&request.profile_id, &c).map_err(lib_error)?;
    let parts = request
        .parts
        .iter()
        .map(SourcePart::part)
        .collect::<Result<Vec<_>, _>>()?;
    let assessment = profile.assess(&parts, &request.limits).map_err(lib_error)?;
    let successors = profile
        .legal_successors(&parts.iter().map(|p| p.slot.clone()).collect::<Vec<_>>())
        .map_err(lib_error)?;
    let mut report = Report {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        profile: summary(&profile),
        sources: request.parts.iter().map(SourcePart::binding).collect(),
        assessment,
        successors,
        report_revision: String::new(),
    };
    report.report_revision = report.content_revision()?;
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl Report {
    fn content_revision(&self) -> Result<String, Error> {
        bounded_json::hash(&(
            &self.engine_version,
            &self.context_revision,
            &self.profile,
            &self.sources,
            &self.assessment,
            &self.successors,
        ))
    }
    /// Current transport integrity only; not authentication of a producer's result.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        bounded_json::preflight(self, MAX_REPORT_BYTES)?;
        if self.engine_version != ENGINE_VERSION
            || self.context_revision != expected
            || expected.len() != 64
            || self.assessment.engine != collagene_assembly_standards::ENGINE_VERSION
            || self.profile.definition.schema_version != 1
            || self.report_revision != self.content_revision()?
            || self.sources.len() != self.assessment.parts.len()
            || self.profile.profile_identity != self.assessment.profile_identity
            || self.sources.len() > 128
        {
            return Err(bad("report context or identity differs"));
        }
        for (source, part) in self.sources.iter().zip(&self.assessment.parts) {
            source.validate()?;
            if source.id != part.id || source.length != part.sequence_length {
                return Err(bad("report source correspondence differs"));
            }
        }
        Ok(())
    }
    pub fn validate_for(&self, request: &Request) -> Result<(), Error> {
        self.validate_binding(&revision(request)?)?;
        if bounded_json::hash(self)? != bounded_json::hash(&assess(request)?)? {
            return Err(bad("fresh assessment differs"));
        }
        Ok(())
    }
}
