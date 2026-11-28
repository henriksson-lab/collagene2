//! Bounded portable reference-map SVG/PNG export, independent of viewport rendering.
use crate::{
    cloning::restriction::inspection::InspectionSite,
    output::{
        raster::{self, Format, Image},
        svg,
    },
};
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub source: Molecule,
    pub source_revision: String,
    pub restriction_sites: Option<Vec<InspectionSite>>,
    pub format: Format,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Rendered {
    pub source_revision: String,
    pub image: Image,
    pub warnings: Vec<String>,
}
impl Request {
    pub fn validate(&self) -> Result<(), Error> {
        self.format.validate()?;
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
            || self.warnings.len() > 32
            || self.warnings.iter().any(|s| s.len() > 4096)
        {
            return Err(Error::Format(
                "Map export identity or output bounds differ".into(),
            ));
        }
        if matches!(&self.image, Image::Svg(svg) if svg.len() > super::svg::MAX_SVG_BYTES) {
            return Err(Error::TooLarge);
        }
        self.image.format().validate()?;
        self.image.validate_binding(revision)
    }
}
impl Request {
    pub fn render(&self) -> Result<Rendered, Error> {
        self.validate()?;
        let scene = svg::layout_map(&self.source)?;
        let output = match &self.restriction_sites {
            Some(sites) => {
                svg::render_map_with_restriction_sites(&self.source, &self.source_revision, sites)?
            }
            None => svg::render_map(&self.source)?,
        };
        let mut warnings = scene.warnings;
        warnings.push("The portable map uses the whole document in its canonical orientation and layout. Viewport rotation, zoom, selection shading, hidden-feature filters, search/PCR overlays and selected-site highlighting are not reproduced. This is a reference map, not evidence of experimental cleavage or material isolation.".into());
        if self
            .restriction_sites
            .as_ref()
            .is_some_and(|s| s.len() > svg::MAX_RESTRICTION_LABELS)
        {
            warnings.push(format!("All supplied restriction markers are retained; at most {} restriction labels are printed. Full coordinates remain in the SVG output or original source evidence; PNG is pixels only.", svg::MAX_RESTRICTION_LABELS));
        }
        let image = match self.format {
            Format::Svg => Image::Svg(output),
            Format::Png(options) => {
                let (png, raster_warnings) =
                    raster::render_generated_svg(&output, options, &self.source_revision)?;
                warnings.extend(raster_warnings);
                Image::Png(png)
            }
        };
        let result = Rendered {
            source_revision: self.source_revision.clone(),
            image,
            warnings,
        };
        result.validate_binding(&self.source_revision)?;
        Ok(result)
    }
}
