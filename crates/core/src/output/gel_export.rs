//! Static gel images retain an exact report/options binding, never physical sample authority.
use super::{
    gel,
    raster::{self, Format, Image},
};
use collagene_model::Error;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Request {
    pub report: gel::GelReport,
    pub report_revision: String,
    pub format: Format,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Rendered {
    pub report_revision: String,
    pub context_revision: String,
    pub image: Image,
    pub warnings: Vec<String>,
}
impl Request {
    pub fn new(report: gel::GelReport, format: Format) -> Result<Self, Error> {
        format.validate()?;
        let report_revision = report.revision()?;
        Ok(Self {
            report,
            report_revision,
            format,
        })
    }
    pub fn revision(&self) -> Result<String, Error> {
        self.format.validate()?;
        // No replay/material claim is inferred from a report supplied for documentary rendering.
        let report_revision = self.report.revision()?;
        if report_revision != self.report_revision {
            return Err(Error::Format("Gel export report revision differs".into()));
        }
        let bytes =
            serde_json::to_vec(&("collagene-gel-image-1", &self.report_revision, self.format))
                .map_err(|e| Error::Format(e.to_string()))?;
        Ok(blake3::hash(&bytes).to_hex().to_string())
    }
    pub fn render(&self) -> Result<Rendered, Error> {
        let context_revision = self.revision()?;
        let svg = gel::render_svg(&self.report)?;
        let (image, warnings) = match self.format {
            Format::Svg => (Image::Svg(svg), vec!["Predicted gel under the exact supplied report/calibration; no observed-gel comparison or experimental sample authority.".into()]),
            Format::Png(options) => {
                let (png, mut warnings) = raster::render_generated_svg(&svg, options, &context_revision)?;
                warnings.push("Predicted gel under the exact supplied report/calibration; no observed-gel comparison or experimental sample authority.".into());
                (Image::Png(png), warnings)
            }
        };
        let result = Rendered {
            report_revision: self.report_revision.clone(),
            context_revision,
            image,
            warnings,
        };
        result.validate_binding(&result.context_revision)?;
        Ok(result)
    }
}
impl Rendered {
    pub fn validate_binding(&self, revision: &str) -> Result<(), Error> {
        if self.context_revision != revision
            || self.report_revision.len() != 64
            || !self
                .report_revision
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            || self.warnings.len() > 32
            || self.warnings.iter().any(|w| w.len() > 4096)
        {
            return Err(Error::Format(
                "Gel image identity or warning bounds differ".into(),
            ));
        }
        let expected = serde_json::to_vec(&(
            "collagene-gel-image-1",
            &self.report_revision,
            self.image.format(),
        ))
        .map_err(|e| Error::Format(e.to_string()))?;
        if blake3::hash(&expected).to_hex().as_str() != revision {
            return Err(Error::Format(
                "Gel image format/context binding differs".into(),
            ));
        }
        self.image.format().validate()?;
        self.image.validate_binding(revision)
    }
}
