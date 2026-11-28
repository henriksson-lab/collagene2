//! Exact source/view capture for portable maps. Exporting never changes biology.
use collagene_core::{
    cloning::restriction::inspection::InspectionSite,
    material::lineage::EventBinding,
    output::{
        raster::{Format, Image},
        svg,
    },
};
use collagene_model::{Molecule, Topology};
use serde::Serialize;

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct Selection {
    /// Zero-based reference boundaries; wrap crosses the circular origin.
    pub start: usize,
    pub end: usize,
    pub wrap: bool,
}
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ViewSettings {
    pub selection: Option<Selection>,
    pub rotation_radians: f64,
    pub zoom: f64,
    pub selected_feature: Option<usize>,
    pub highlighted_spans: Vec<(usize, usize)>,
    /// App-owned presentation settings key (label/display modes, not authority).
    pub display_context: String,
    /// Sorted, unique original feature indices hidden by the viewport.
    pub hidden_features: Vec<usize>,
}
#[derive(Clone, Debug, Serialize)]
pub struct RestrictionLayer {
    pub source_revision: String,
    /// Pin the scan/filter settings even if two settings display identical sites.
    pub display_context: String,
    pub sites: Vec<InspectionSite>,
    pub selected_site: Option<usize>,
}
#[derive(Clone)]
pub struct Snapshot {
    pub source: Molecule,
    pub event: EventBinding,
    pub epoch: u64,
    pub view: ViewSettings,
    /// None means restriction markers are hidden, not a request to rescan DNA.
    pub restrictions: Option<RestrictionLayer>,
}
#[derive(Clone)]
pub struct Capture {
    pub snapshot: Snapshot,
    pub context_revision: String,
    pub format: Format,
}
#[derive(Clone)]
pub struct RenderRequest {
    pub source: Molecule,
    pub source_revision: String,
    pub context_revision: String,
    pub restriction_sites: Option<Vec<InspectionSite>>,
    pub format: Format,
}
#[derive(Clone)]
pub struct Rendered {
    pub source_revision: String,
    pub context_revision: String,
    pub image: Image,
    pub warnings: Vec<String>,
}
#[derive(Clone)]
pub struct Prepared {
    pub capture: Capture,
    pub rendered: Rendered,
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
}
fn identifier(s: &str) -> bool {
    !s.is_empty() && s.len() <= 128 && !s.chars().any(|c| c.is_control() || c.is_whitespace())
}
impl Capture {
    pub fn new(snapshot: Snapshot, format: Format) -> Result<Self, String> {
        format.validate().map_err(|e| e.to_string())?;
        collagene_model::bounds::model_preflight(&[&snapshot.source]).map_err(|e| e.to_string())?;
        snapshot.source.validate().map_err(|e| e.to_string())?;
        let n = snapshot.source.sequence.len();
        if !identifier(&snapshot.event.document_id)
            || !snapshot
                .event
                .event_id
                .strip_prefix("event-")
                .is_some_and(hex)
            || snapshot.event.revision != snapshot.source.revision()
            || !snapshot.view.rotation_radians.is_finite()
            || !snapshot.view.zoom.is_finite()
            || snapshot.view.zoom <= 0.
            || snapshot.view.zoom > 1_000_000.
            || snapshot.view.display_context.len() > 4096
            || snapshot
                .view
                .selected_feature
                .is_some_and(|i| i >= snapshot.source.features.len())
            || snapshot.view.hidden_features.len() > snapshot.source.features.len()
            || snapshot
                .view
                .hidden_features
                .iter()
                .any(|i| *i >= snapshot.source.features.len())
            || snapshot
                .view
                .hidden_features
                .windows(2)
                .any(|w| w[0] >= w[1])
        {
            return Err("Map export source event or view settings are invalid.".into());
        }
        if snapshot.view.highlighted_spans.len() > 4096
            || snapshot
                .view
                .highlighted_spans
                .iter()
                .any(|(a, b)| a >= b || *b > n)
        {
            return Err("Map export highlighted ranges exceed source bounds.".into());
        }
        if let Some(range) = &snapshot.view.selection {
            if range.start > n
                || range.end > n
                || (range.wrap
                    && (snapshot.source.topology != Topology::Circular || range.start <= range.end))
                || (!range.wrap && range.start > range.end)
            {
                return Err("Map export selection is outside the exact source.".into());
            }
        }
        if let Some(layer) = &snapshot.restrictions {
            if layer.source_revision != snapshot.event.revision
                || layer.display_context.len() > 4096
                || layer.sites.len() > svg::MAX_RESTRICTION_SITES
                || layer
                    .selected_site
                    .is_some_and(|i| !layer.sites.iter().any(|s| s.index == i))
            {
                return Err(
                    "Displayed restriction sites belong to a stale or oversized map.".into(),
                );
            }
            // Bind bounded coordinate evidence, not arbitrary unbounded foreign text.
            for site in &layer.sites {
                if site.enzyme.len() > 256 || site.spans.len() > 8 || site.diagnostics.len() > 8 {
                    return Err(
                        "Displayed restriction-site annotation exceeds export bounds.".into(),
                    );
                }
            }
        }
        let bytes = serde_json::to_vec(&(
            "collagene-map-export-v2",
            format,
            &snapshot.event,
            snapshot.epoch,
            &snapshot.view,
            &snapshot.restrictions,
        ))
        .map_err(|e| e.to_string())?;
        let context_revision = blake3::hash(&bytes).to_hex().to_string();
        Ok(Self {
            snapshot,
            context_revision,
            format,
        })
    }
    pub fn matches(&self, current: Snapshot) -> Result<(), String> {
        if Self::new(current, self.format)?.context_revision != self.context_revision {
            return Err(
                "Document, event, selection or map settings changed; prepare a new image.".into(),
            );
        }
        Ok(())
    }
    pub fn matches_capture(&self, current: &Self) -> Result<(), String> {
        if self.context_revision != current.context_revision {
            return Err(
                "Document, event, view or export format changed; prepare a new image.".into(),
            );
        }
        Ok(())
    }
    pub fn request(&self) -> RenderRequest {
        RenderRequest {
            source: self.snapshot.source.clone(),
            source_revision: self.snapshot.event.revision.clone(),
            context_revision: self.context_revision.clone(),
            restriction_sites: self.snapshot.restrictions.as_ref().map(|s| s.sites.clone()),
            format: self.format,
        }
    }
    pub fn accept(self, rendered: Rendered) -> Result<Prepared, String> {
        if rendered.context_revision != self.context_revision
            || rendered.source_revision != self.snapshot.event.revision
            || rendered.image.format() != self.format
            || rendered.warnings.len() > 128
            || rendered.warnings.iter().any(|w| w.len() > 4096)
        {
            return Err(
                "Image response does not match the captured source/view or its bounds.".into(),
            );
        }
        rendered
            .image
            .validate_binding(&self.snapshot.event.revision)
            .map_err(|e| e.to_string())?;
        Ok(Prepared {
            capture: self,
            rendered,
        })
    }
}
impl Prepared {
    pub fn filename(&self) -> String {
        let mut stem: String = self
            .capture
            .snapshot
            .source
            .name
            .chars()
            .map(|c| {
                if c.is_control()
                    || matches!(c, '/' | '\\' | ':' | '<' | '>' | '"' | '|' | '?' | '*')
                {
                    '_'
                } else {
                    c
                }
            })
            .collect();
        while stem.len() > 180 {
            stem.pop();
        }
        let stem = stem.trim().trim_matches('.').trim();
        if stem.is_empty() {
            format!(
                "map-{}.{}",
                &self.capture.snapshot.event.revision[..12],
                self.capture.format.extension()
            )
        } else {
            format!("map-{stem}.{}", self.capture.format.extension())
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn snapshot() -> Snapshot {
        let source = Molecule::dna("literal plasmid", "ACGTTGCA", Topology::Circular).unwrap();
        Snapshot {
            event: EventBinding {
                document_id: "document-1".into(),
                event_id: format!("event-{}", "1".repeat(64)),
                revision: source.revision(),
            },
            source,
            epoch: 7,
            view: ViewSettings {
                selection: None,
                rotation_radians: 0.25,
                zoom: 1.,
                selected_feature: None,
                highlighted_spans: vec![],
                display_context: "reference map".into(),
                hidden_features: vec![],
            },
            restrictions: None,
        }
    }
    #[test]
    fn equal_sequence_other_document_event_epoch_or_view_invalidates_capture() {
        let initial = snapshot();
        let bound = Capture::new(initial.clone(), Format::Svg).unwrap();
        bound.matches(initial.clone()).unwrap();
        for changed in 0..8 {
            let mut next = initial.clone();
            match changed {
                0 => next.event.document_id = "document-2".into(),
                1 => next.event.event_id = format!("event-{}", "2".repeat(64)),
                2 => next.epoch += 1,
                3 => next.view.rotation_radians += 0.1,
                4 => next.view.zoom = 2.,
                5 => {
                    next.view.selection = Some(Selection {
                        start: 1,
                        end: 3,
                        wrap: false,
                    })
                }
                6 => next.view.highlighted_spans = vec![(1, 3)],
                _ => next.view.display_context = "different label settings".into(),
            }
            assert!(bound.matches(next).is_err());
        }
        let mut stale = initial;
        stale.source.name = "changed name".into();
        assert!(Capture::new(stale, Format::Svg).is_err());
    }
    #[test]
    fn renderer_reply_and_restriction_filter_context_cannot_cross_captures() {
        let mut initial = snapshot();
        initial.restrictions = Some(RestrictionLayer {
            source_revision: initial.event.revision.clone(),
            display_context: "single cutters".into(),
            sites: vec![],
            selected_site: None,
        });
        let first = Capture::new(initial.clone(), Format::Svg).unwrap();
        let request = first.request();
        let rendered = collagene_core::output::map_export::Request {
            source: request.source,
            source_revision: request.source_revision,
            restriction_sites: request.restriction_sites,
            format: request.format,
        }
        .render()
        .unwrap();
        let output = Rendered {
            source_revision: rendered.source_revision,
            context_revision: request.context_revision,
            image: rendered.image,
            warnings: rendered.warnings,
        };
        let mut other = initial.clone();
        other.restrictions.as_mut().unwrap().display_context = "all enzymes".into();
        assert!(first.matches(other.clone()).is_err());
        assert!(
            Capture::new(other, Format::Svg)
                .unwrap()
                .accept(output.clone())
                .is_err()
        );
        let prepared = first.accept(output).unwrap();
        assert_eq!(prepared.filename(), "map-literal plasmid.svg");
        let mut invalid = initial;
        invalid.view.rotation_radians = f64::NAN;
        assert!(Capture::new(invalid, Format::Svg).is_err());
    }
    #[test]
    fn image_format_and_pixel_width_are_part_of_the_captured_map_context() {
        let original = snapshot();
        let svg = Capture::new(original.clone(), Format::Svg).unwrap();
        let png = Capture::new(
            original.clone(),
            Format::Png(collagene_core::output::raster::PngOptions { width: 400 }),
        )
        .unwrap();
        let wider = Capture::new(
            original.clone(),
            Format::Png(collagene_core::output::raster::PngOptions { width: 800 }),
        )
        .unwrap();
        assert!(svg.matches_capture(&png).is_err());
        assert!(png.matches_capture(&wider).is_err());
        assert!(
            Capture::new(
                original,
                Format::Png(collagene_core::output::raster::PngOptions { width: 0 })
            )
            .is_err()
        );
        let request = png.request();
        let rendered = collagene_core::output::map_export::Request {
            source: request.source,
            source_revision: request.source_revision,
            restriction_sites: request.restriction_sites,
            format: request.format,
        }
        .render()
        .unwrap();
        let returned = Rendered {
            source_revision: rendered.source_revision,
            context_revision: request.context_revision,
            image: rendered.image,
            warnings: rendered.warnings,
        };
        assert!(svg.accept(returned.clone()).is_err());
        let prepared = png.accept(returned).unwrap();
        assert_eq!(prepared.filename(), "map-literal plasmid.png");
        assert_eq!(prepared.capture.format.mime_type(), "image/png");
    }
}
