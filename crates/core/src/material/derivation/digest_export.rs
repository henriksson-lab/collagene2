//! Detached digest references reconstructed from the same tracked-material
//! execution used for native derivations. Exported summaries never restore it.
use super::{BoundSource, Derivation, ModificationSummary, SourceRole};
use crate::{chemistry::DuplexSummary, cloning::restriction::setups::Recipe};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};

pub const ENGINE_VERSION: &str = "material-digest-export-1";
const LOSS: &str = "Detached sequence export does not preserve reusable physical material, exact native history, both strand sequences and offsets, terminal chemistry, internal nicks or positional modification evidence. Any included comments and summaries are documentary only. Retain the native workspace and replay its exact origin for subsequent cloning. Complete cleavage and experimental fragment isolation remain assumptions; no intended target is assessed.";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub recipe: Recipe,
    pub fragment_index: usize,
    pub name: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Product {
    pub engine_version: String,
    pub context_revision: String,
    pub inspection_revision: String,
    pub molecule: Molecule,
    pub duplex: DuplexSummary,
    pub modifications: ModificationSummary,
    pub incomplete: bool,
    pub findings: Vec<String>,
    pub findings_total: usize,
    pub findings_truncated: bool,
}

impl Request {
    pub fn revision(&self, source: &Molecule) -> Result<String, Error> {
        self.recipe.validate()?;
        super::name(&self.name)?;
        if self.fragment_index >= self.recipe.options.max_fragments {
            return Err(super::bad(
                "export fragment exceeds the selected digest limit",
            ));
        }
        let inspection = self.recipe.inspection_revision(source)?;
        crate::encoding::bounded_json::hash(&(
            ENGINE_VERSION,
            inspection,
            self.fragment_index,
            &self.name,
        ))
    }

    pub fn materialize(&self, source: &Molecule) -> Result<Product, Error> {
        let context_revision = self.revision(source)?;
        let request = super::Request {
            sources: vec![BoundSource {
                role: SourceRole::Template,
                event: self.recipe.source.clone(),
                origin: self.recipe.origin.clone(),
            }],
            derivation: Derivation::CatalogDigest {
                selection: self.recipe.selection.clone(),
                options: self.recipe.options.clone(),
                fragment_index: self.fragment_index,
            },
            name: self.name.clone(),
            save_primers: false,
        };
        request.validate_shallow()?;
        let fresh = super::fresh_template(source, &request, &context_revision)?;
        let (mut molecule, tracked, incomplete, mut findings, total, truncated) =
            fresh.into_parts();
        molecule
            .metadata
            .insert("source:material_digest_export:losses".into(), LOSS.into());
        findings.insert(0, LOSS.into());
        let (findings, findings_total, findings_truncated) = super::limitations(
            findings,
            total.checked_add(1).ok_or(Error::TooLarge)?,
            truncated,
        );
        let product = Product {
            engine_version: ENGINE_VERSION.into(),
            context_revision,
            inspection_revision: self.recipe.inspection_revision(source)?,
            molecule,
            duplex: tracked.duplex.summary()?,
            modifications: tracked.summary()?,
            incomplete,
            findings,
            findings_total,
            findings_truncated,
        };
        product.validate_bindings_for(self, source)?;
        crate::protocol::save_document(&product.molecule)?;
        Ok(product)
    }
}

impl Product {
    /// Documentary transport checks only. This cannot restore tracked strands.
    pub fn validate(&self) -> Result<(), Error> {
        self.molecule.validate()?;
        super::summary_shape(&self.modifications)?;
        if self.engine_version != ENGINE_VERSION
            || !super::hex(&self.context_revision)
            || !super::hex(&self.inspection_revision)
            || self.molecule.sequence.len() != self.duplex.top_length
            || self.molecule.topology != self.duplex.topology
            || self.duplex.top_length != self.modifications.top_length
            || self.duplex.bottom_length != self.modifications.bottom_length
            || self.findings.len() != self.findings_total.min(128)
            || self.findings.iter().any(|s| s.chars().count() > 2048)
            || (self.findings_total > self.findings.len() && !self.findings_truncated)
            || self
                .molecule
                .metadata
                .get("source:material_digest_export:losses")
                .map(String::as_str)
                != Some(LOSS)
        {
            return Err(super::bad(
                "detached digest export context or summary differs",
            ));
        }
        super::encoded(self, super::MAX_CAPSULE_BYTES)?;
        Ok(())
    }

    pub fn validate_bindings_for(&self, request: &Request, source: &Molecule) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != request.revision(source)?
            || self.inspection_revision != request.recipe.inspection_revision(source)?
            || self.molecule.name != request.name
        {
            return Err(super::bad(
                "detached digest export source/origin/options changed",
            ));
        }
        Ok(())
    }
}
