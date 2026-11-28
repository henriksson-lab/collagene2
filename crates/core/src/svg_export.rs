//! Bounded portable reference-map export, independent of viewport rendering.
use crate::{restriction_inspection::InspectionSite, svg};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub source: Molecule,
    pub source_revision: String,
    pub restriction_sites: Option<Vec<InspectionSite>>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Rendered {
    pub source_revision: String,
    pub svg: String,
    pub warnings: Vec<String>,
}
impl Request {
    pub fn validate(&self) -> Result<(), Error> {
        collagene_model::bounds::model_preflight(&[&self.source])?;
        self.source.validate()?;
        if self.source_revision != self.source.revision() {
            return Err(Error::Format("Map export source revision changed".into()));
        }
        if self
            .restriction_sites
            .as_ref()
            .is_some_and(|s| s.len() > svg::MAX_RESTRICTION_SITES)
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
impl Rendered {
    /// Documentary identity and output limits, not a grant of biological authority.
    pub fn validate_binding(&self, revision: &str) -> Result<(), Error> {
        if self.source_revision != revision
            || revision.len() != 64
            || !revision
                .bytes()
                .all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase())
            || self.svg.len() > svg::MAX_SVG_BYTES
            || !self.svg.starts_with("<svg ")
            || !self.svg.ends_with("</svg>")
            || self.warnings.len() > 32
            || self.warnings.iter().any(|s| s.len() > 4096)
        {
            return Err(Error::Format(
                "Map export identity or output bounds differ".into(),
            ));
        }
        Ok(())
    }
}
pub fn render(request: &Request) -> Result<Rendered, Error> {
    request.validate()?;
    let scene = svg::layout_map(&request.source)?;
    let output = match &request.restriction_sites {
        Some(sites) => svg::render_map_with_restriction_sites(
            &request.source,
            &request.source_revision,
            sites,
        )?,
        None => svg::render_map(&request.source)?,
    };
    let mut warnings = scene.warnings;
    warnings.push("Portable SVG uses the whole document in its canonical orientation and layout. Viewport rotation, zoom, selection shading, hidden-feature filters, search/PCR overlays and selected-site highlighting are not reproduced. This is a reference map, not evidence of experimental cleavage or material isolation.".into());
    if request
        .restriction_sites
        .as_ref()
        .is_some_and(|s| s.len() > svg::MAX_RESTRICTION_LABELS)
    {
        warnings.push(format!("All supplied restriction markers are retained; at most {} restriction labels are printed. Full coordinates remain in the SVG evidence.", svg::MAX_RESTRICTION_LABELS));
    }
    let result = Rendered {
        source_revision: request.source_revision.clone(),
        svg: output,
        warnings,
    };
    result.validate_binding(&request.source_revision)?;
    Ok(result)
}
