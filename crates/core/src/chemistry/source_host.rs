//! Exact retained source-host provenance, independent of molecular chemistry.
//! Current membership is required; an old hash, strain label or caller-provided
//! source URL cannot confer current documentary authority or modification calls.
use collagene_enzyme_catalog::{Catalog, Source, hosts::HostPreset};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, sync::OnceLock};

pub const MAX_BYTES: usize = 64 * 1024;
const MAX_SOURCES: usize = 8;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SourceHost {
    pub catalog_edition: String,
    pub catalog_content_hash: String,
    pub preset: HostPreset,
    /// Exactly the cited source records, ordered by source ID. No extra sources.
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, 8>")]
    pub sources: Vec<Source>,
}
fn bad(message: &str) -> Error {
    Error::Format(format!("source host: {message}"))
}
fn current() -> Result<&'static Catalog, Error> {
    static CATALOG: OnceLock<Result<Catalog, String>> = OnceLock::new();
    CATALOG
        .get_or_init(|| Catalog::builtin().map_err(|error| error.to_string()))
        .as_ref()
        .map_err(|message| bad(message))
}
impl SourceHost {
    /// Capture the current immutable preset, catalogue identity and complete
    /// cited provenance together. This does not declare molecular modification.
    pub fn from_current(id: &str) -> Result<Self, Error> {
        if id.is_empty() || id.len() > 96 || id.chars().any(char::is_control) {
            return Err(bad("invalid exact preset identity"));
        }
        let catalog = current()?;
        let preset = catalog.host(id).map_err(|error| bad(&error.to_string()))?;
        let ids = Self::cited_ids(preset)?;
        let sources = ids
            .into_iter()
            .map(|id| {
                catalog
                    .edition()
                    .sources
                    .iter()
                    .find(|source| source.id == id)
                    .ok_or_else(|| bad("cited source is absent from current catalogue"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        // Borrowed tuple admission precedes deep copies. Reserve more than the
        // named transport fields' JSON overhead relative to this tuple.
        crate::encoding::bounded_json::preflight(
            &(
                catalog.edition().edition.as_str(),
                catalog.content_hash(),
                preset,
                &sources,
            ),
            MAX_BYTES - 256,
        )?;
        Ok(Self {
            catalog_edition: catalog.edition().edition.clone(),
            catalog_content_hash: catalog.content_hash().to_owned(),
            preset: preset.clone(),
            sources: sources.into_iter().cloned().collect(),
        })
    }
    fn cited_ids(preset: &HostPreset) -> Result<BTreeSet<&str>, Error> {
        preset.validate().map_err(|error| bad(&error.to_string()))?;
        let ids = preset
            .citations
            .iter()
            .map(|citation| citation.source.as_str())
            .collect::<BTreeSet<_>>();
        if ids.is_empty() || ids.len() > MAX_SOURCES {
            return Err(Error::TooLarge);
        }
        Ok(ids)
    }
    /// Validate the entire current pin and exact ordered cited source membership.
    /// This rejects historical editions rather than rebinding retained provenance.
    pub fn validate(&self) -> Result<(), Error> {
        if self.sources.is_empty() || self.sources.len() > MAX_SOURCES {
            return Err(Error::TooLarge);
        }
        crate::encoding::bounded_json::preflight(self, MAX_BYTES)?;
        let catalog = current()?;
        if self.catalog_edition != catalog.edition().edition
            || self.catalog_content_hash != catalog.content_hash()
        {
            return Err(bad(
                "catalogue edition/content pin differs from current catalogue",
            ));
        }
        catalog
            .validate_host(&self.preset)
            .map_err(|error| bad(&error.to_string()))?;
        let ids = Self::cited_ids(&self.preset)?;
        if ids.len() != self.sources.len() {
            return Err(bad("missing, duplicate or additional cited source records"));
        }
        for (id, source) in ids.into_iter().zip(&self.sources) {
            if source.id != id
                || catalog
                    .edition()
                    .sources
                    .iter()
                    .find(|candidate| candidate.id == id)
                    != Some(source)
            {
                return Err(bad("cited source record or canonical source order differs"));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_enzyme_catalog::hosts::MethylaseExpectation as M;

    fn fixture() -> SourceHost {
        SourceHost::from_current("neb-dam-dcm-c2925").unwrap()
    }

    #[test]
    fn capture_retains_exact_supplier_preset_pin_and_complete_canonical_sources() {
        let host = fixture();
        host.validate().unwrap();
        assert_eq!(host.catalog_edition, "2026-11-28.2");
        assert_eq!(host.catalog_content_hash, current().unwrap().content_hash());
        assert_eq!(
            (host.preset.dam, host.preset.dcm, host.preset.cpg),
            (M::Deficient, M::Deficient, M::Unknown)
        );
        assert_eq!(
            host.sources
                .iter()
                .map(|source| source.id.as_str())
                .collect::<Vec<_>>(),
            ["neb-dam-dcm-cells", "neb-host-selection"]
        );
        assert_eq!(
            host.sources[0].url,
            "https://www.neb.com/en/products/c2925-dam-dcm-competent-e-coli"
        );
        assert!(
            host.sources.iter().all(
                |source| source.retrieved_on == "2026-11-28" && source.facts_blake3.len() == 64
            )
        );
        assert!(std::ptr::eq(current().unwrap(), current().unwrap()));
        assert!(SourceHost::from_current("DH5alpha").is_err());
    }
    #[test]
    fn current_validation_rejects_changed_pin_preset_and_every_provenance_dimension() {
        let host = fixture();
        let mut altered = host.clone();
        altered.catalog_edition = "2026-11-28.1".into();
        assert!(altered.validate().is_err());
        let mut altered = host.clone();
        altered.catalog_content_hash = "0".repeat(64);
        assert!(altered.validate().is_err());
        let mut altered = host.clone();
        altered.preset.dam = M::Active;
        assert!(altered.validate().is_err());
        let mut altered = host.clone();
        altered.preset.citations[0].locator = "different evidence".into();
        assert!(altered.validate().is_err());
        for dimension in 0..7 {
            let mut altered = host.clone();
            let source = &mut altered.sources[0];
            match dimension {
                0 => source.publisher = "different publisher".into(),
                1 => source.url = "https://example.invalid/claimed-source".into(),
                2 => source.retrieved_on = "2026-11-27".into(),
                3 => source.edition = "different source version".into(),
                4 => source.snapshot_note = "different capture scope".into(),
                5 => source.facts_blake3 = "0".repeat(64),
                _ => {
                    source.snapshot = Some(collagene_enzyme_catalog::Snapshot {
                        sha256: "0".repeat(64),
                        bytes: 1,
                    })
                }
            }
            assert!(
                altered.validate().is_err(),
                "changed source dimension {dimension}"
            );
        }
    }
    #[test]
    fn missing_extra_duplicate_and_noncanonical_source_records_are_refused() {
        let host = fixture();
        let mut missing = host.clone();
        missing.sources.pop();
        assert!(missing.validate().is_err());
        let mut duplicate = host.clone();
        duplicate.sources[1] = duplicate.sources[0].clone();
        assert!(duplicate.validate().is_err());
        let mut extra = host.clone();
        extra.sources.push(
            current()
                .unwrap()
                .edition()
                .sources
                .iter()
                .find(|source| source.id == "agilent-jm110")
                .unwrap()
                .clone(),
        );
        assert!(extra.validate().is_err());
        let mut reordered = host.clone();
        reordered.sources.reverse();
        assert!(reordered.validate().is_err());
    }
    #[test]
    fn current_json_and_postcard_preserve_every_retained_source_without_fallback() {
        let host = fixture();
        let json = serde_json::to_value(&host).unwrap();
        let decoded: SourceHost = serde_json::from_value(json.clone()).unwrap();
        assert_eq!(decoded, host);
        decoded.validate().unwrap();
        let bytes = postcard::to_allocvec(&host).unwrap();
        let decoded: SourceHost = postcard::from_bytes(&bytes).unwrap();
        assert_eq!(decoded, host);
        decoded.validate().unwrap();
        let mut missing = json.clone();
        missing.as_object_mut().unwrap().remove("sources");
        assert!(serde_json::from_value::<SourceHost>(missing).is_err());
        let mut missing = json.clone();
        missing
            .as_object_mut()
            .unwrap()
            .remove("catalog_content_hash");
        assert!(serde_json::from_value::<SourceHost>(missing).is_err());
        let mut extra = json;
        extra["verified"] = true.into();
        assert!(serde_json::from_value::<SourceHost>(extra).is_err());
    }
    #[test]
    fn source_count_and_bytes_are_bounded_before_current_membership_comparison() {
        let host = fixture();
        let mut excess = host.clone();
        excess.sources = vec![host.sources[0].clone(); 9];
        assert!(matches!(excess.validate(), Err(Error::TooLarge)));
        let json = serde_json::to_value(excess).unwrap();
        assert!(serde_json::from_value::<SourceHost>(json).is_err());
        let mut oversized = host;
        oversized.sources[0].url = "x".repeat(MAX_BYTES + 1);
        assert!(matches!(oversized.validate(), Err(Error::TooLarge)));
    }
}
