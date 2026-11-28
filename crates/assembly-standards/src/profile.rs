//! Immutable installed grammar and its exact catalog binding.
use crate::assessment::{Assessment, Examiner, Limits, Part};
use crate::grammar::{Grammar, Successors};
use crate::{Error, MAX_PROFILE_BYTES, Result, bounded_json, mask, require, text};
use collagene_enzyme_catalog::{Catalog, Geometry};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub const RFC10_ID: &str = "biobrick-rfc10-part-sequence-2007";
pub const RFC1000_L0_ID: &str = "igem-rfc1000-level0-documented-body-2026";
pub const BUILTIN_IDS: &[&str] = &[RFC10_ID, RFC1000_L0_ID];

fn builtin_bytes(id: &str) -> Result<&'static [u8]> {
    match id {
        RFC10_ID => Ok(include_bytes!("../data/rfc10-2007/profile.json")),
        RFC1000_L0_ID => Ok(include_bytes!("../data/rfc1000-l0-2026/profile.json")),
        _ => Err(Error::UnsupportedProfile),
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ProfileKind {
    Builtin,
    Custom,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub id: String,
    pub author: String,
    pub url: String,
    pub revision: String,
    pub reviewed_on: String,
    pub locator: String,
    pub source_license: String,
    pub upstream_sha256: Option<String>,
    pub snapshot_note: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CatalogPin {
    pub edition: String,
    pub content_hash: String,
}
impl CatalogPin {
    pub fn current(catalog: &Catalog) -> Self {
        Self {
            edition: catalog.edition().edition.clone(),
            content_hash: catalog.content_hash().into(),
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct EnzymeBinding {
    pub id: String,
    pub geometry: Geometry,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Requirement {
    Required,
    Advisory,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AvoidedSite {
    pub id: String,
    pub name: String,
    pub recognition: String,
    pub requirement: Requirement,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 16>")]
    pub source_ids: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Definition {
    pub schema_version: u16,
    pub id: String,
    pub revision: String,
    pub name: String,
    pub scope: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 16>")]
    pub limitations: Vec<String>,
    pub curation_license: String,
    pub expected_catalog_edition: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, 16>")]
    pub sources: Vec<Source>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, EnzymeBinding, 32>")]
    pub enzymes: Vec<EnzymeBinding>,
    pub grammar: Grammar,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, AvoidedSite, 64>")]
    pub avoided_sites: Vec<AvoidedSite>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 64>")]
    pub default_avoided_sites: Vec<String>,
}
impl Definition {
    pub fn validate(&self, catalog: &Catalog) -> Result<()> {
        // Native DTO proposals receive the same bounded serialization before
        // normalization, hashing or any installation clone.
        bounded_json(self, MAX_PROFILE_BYTES)?;
        if self.schema_version != 1 {
            return Err(Error::UnsupportedVersion);
        }
        require(
            text(&self.id, 96)
                && text(&self.revision, 128)
                && text(&self.name, 256)
                && text(&self.scope, 1024)
                && text(&self.curation_license, 128)
                && self.limitations.len() <= 16
                && self.limitations.iter().all(|s| text(s, 2048))
                && !self.sources.is_empty()
                && self.sources.len() <= 16
                && !self.enzymes.is_empty()
                && self.enzymes.len() <= 32
                && self.avoided_sites.len() <= 64
                && self.default_avoided_sites.len() <= 64,
            "profile inventory/metadata",
        )?;
        if self.expected_catalog_edition != catalog.edition().edition {
            return Err(Error::CatalogMismatch);
        }
        self.grammar.validate()?;
        let mut sources = BTreeSet::new();
        for source in &self.sources {
            require(
                text(&source.id, 96)
                    && sources.insert(source.id.as_str())
                    && text(&source.author, 256)
                    && text(&source.url, 2048)
                    && source.url.starts_with("https://")
                    && text(&source.revision, 256)
                    && text(&source.reviewed_on, 32)
                    && text(&source.locator, 1024)
                    && text(&source.source_license, 128)
                    && text(&source.snapshot_note, 1024)
                    && source.upstream_sha256.as_ref().is_none_or(|s| hex(s)),
                "source provenance",
            )?;
        }
        let mut enzymes = BTreeSet::new();
        for binding in &self.enzymes {
            require(
                text(&binding.id, 96) && enzymes.insert(binding.id.as_str()),
                "enzyme identity",
            )?;
            let installed = catalog
                .edition()
                .enzymes
                .iter()
                .find(|e| e.id == binding.id)
                .ok_or(Error::CatalogMismatch)?;
            if installed.geometry != binding.geometry {
                return Err(Error::CatalogMismatch);
            }
        }
        let mut sites = BTreeSet::new();
        for site in &self.avoided_sites {
            require(
                text(&site.id, 96)
                    && sites.insert(site.id.as_str())
                    && text(&site.name, 128)
                    && !site.recognition.is_empty()
                    && site.recognition.len() <= 32
                    && site.recognition.bytes().all(|b| mask(b) != 0)
                    && !site.source_ids.is_empty()
                    && site.source_ids.len() <= 16
                    && site.source_ids.iter().all(|s| sources.contains(s.as_str())),
                "avoided-site rule",
            )?;
        }
        let refs_valid = |refs: &[String]| {
            let mut unique = BTreeSet::new();
            refs.iter()
                .all(|id| sites.contains(id.as_str()) && unique.insert(id.as_str()))
        };
        require(
            refs_valid(&self.default_avoided_sites),
            "default avoided sites",
        )?;
        for slot in &self.grammar.slots {
            require(
                enzymes.contains(slot.left_enzyme.as_str())
                    && enzymes.contains(slot.right_enzyme.as_str())
                    && refs_valid(&slot.extra_avoided_sites)
                    && !slot.source_ids.is_empty()
                    && slot.source_ids.iter().all(|s| sources.contains(s.as_str())),
                "slot evidence/enzyme/site references",
            )?;
        }
        Ok(())
    }
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

/// A validated installed definition. Neither imported source citations nor a
/// checksum establish upstream authenticity or experimental authority.
pub struct Profile {
    definition: Definition,
    catalog_pin: CatalogPin,
    identity: String,
    kind: ProfileKind,
}
impl Profile {
    pub fn new(definition: Definition, pin: CatalogPin, catalog: &Catalog) -> Result<Self> {
        definition.validate(catalog)?;
        if pin != CatalogPin::current(catalog) {
            return Err(Error::CatalogMismatch);
        }
        let kind = if BUILTIN_IDS.contains(&definition.id.as_str()) {
            let canonical: Definition = serde_json::from_slice(builtin_bytes(&definition.id)?)
                .map_err(|_| Error::Invalid("embedded standard definition"))?;
            if definition != canonical {
                return Err(Error::Invalid(
                    "reserved built-in profile identity/facts altered",
                ));
            }
            ProfileKind::Builtin
        } else {
            ProfileKind::Custom
        };
        let bytes = bounded_json(&(kind, &definition, &pin), MAX_PROFILE_BYTES)?;
        let mut h = blake3::Hasher::new();
        h.update(b"assembly-standard-profile-1\0");
        h.update(&bytes);
        Ok(Self {
            definition,
            catalog_pin: pin,
            identity: h.finalize().to_hex().to_string(),
            kind,
        })
    }
    /// Import a current definition, with an explicit separately supplied exact
    /// catalog pin. There is no old-format reader or implicit stale rebinding.
    pub fn from_json(bytes: &[u8], pin: CatalogPin, catalog: &Catalog) -> Result<Self> {
        if bytes.len() > MAX_PROFILE_BYTES {
            return Err(Error::Limit("profile JSON bytes"));
        }
        let definition =
            serde_json::from_slice(bytes).map_err(|_| Error::Invalid("profile JSON"))?;
        Self::new(definition, pin, catalog)
    }
    pub fn builtin(id: &str, catalog: &Catalog) -> Result<Self> {
        Self::from_json(builtin_bytes(id)?, CatalogPin::current(catalog), catalog)
    }
    pub fn definition(&self) -> &Definition {
        &self.definition
    }
    pub fn catalog_pin(&self) -> &CatalogPin {
        &self.catalog_pin
    }
    pub fn identity(&self) -> &str {
        &self.identity
    }
    pub fn kind(&self) -> ProfileKind {
        self.kind
    }
    pub fn validate(&self, catalog: &Catalog) -> Result<()> {
        if self.catalog_pin != CatalogPin::current(catalog) {
            return Err(Error::CatalogMismatch);
        }
        self.definition.validate(catalog)
    }
    pub fn legal_successors(&self, prefix: &[String]) -> Result<Successors> {
        self.definition.grammar.legal_successors(prefix)
    }
    pub fn assess(&self, parts: &[Part], limits: &Limits) -> Result<Assessment> {
        Examiner::assess(self, parts, limits)
    }
}
