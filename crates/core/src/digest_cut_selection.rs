//! Choose one actual catalog-digested fragment by explicit cut instances/free ends.
//! This adapter never merges reference pieces into a purported physical product.
use crate::{
    assembly_recipe,
    catalog_digest::{
        self, CatalogDigestProduct, CatalogRestrictionInspection, CatalogSelection,
        InspectionOptions, catalog,
    },
    restriction_inspection::{FragmentSummary, RestrictionInspection},
};
use collagene_model::{Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "digest-cut-selection-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)] pub enum Boundary {SourceStart,SourceEnd,Cut{site_index:usize}}}
dto! {#[derive(PartialEq,Eq)] pub struct Choice {pub first:Boundary,pub second:Boundary,pub complement:bool}}
dto! {pub struct Report {pub context_revision:String,pub choice:Choice,pub inspection:CatalogRestrictionInspection,pub product:CatalogDigestProduct}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("digest cut selection: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn choice_shape(choice: &Choice, topology: Topology, max_sites: usize) -> Result<(), Error> {
    for b in [&choice.first, &choice.second] {
        match b {
            Boundary::Cut { site_index } if *site_index >= max_sites || *site_index >= 100_000 => {
                return Err(Error::TooLarge);
            }
            Boundary::SourceStart | Boundary::SourceEnd if topology != Topology::Linear => {
                return Err(bad(
                    "source start/end boundaries apply only to an actually linear source",
                ));
            }
            _ => {}
        }
    }
    if topology == Topology::Linear && choice.complement {
        return Err(bad(
            "linear complement would be disconnected pieces; select a single actual fragment explicitly",
        ));
    }
    Ok(())
}
fn context(batch: &str, choice: &Choice) -> Result<String, Error> {
    if !hex(batch) {
        return Err(bad("invalid expected catalog inspection batch revision"));
    }
    assembly_recipe::hash(&(ENGINE_VERSION, batch, choice))
}
#[allow(clippy::too_many_arguments)]
pub fn selection_revision(
    catalog: &catalog::Catalog,
    selection: &CatalogSelection,
    source_role: &str,
    source: &Molecule,
    options: &InspectionOptions,
    expected_batch_revision: &str,
    choice: &Choice,
) -> Result<String, Error> {
    collagene_model::bounds::model_preflight(&[source])?;
    catalog_digest::source_preflight(catalog, selection, source_role, source, options)?;
    choice_shape(choice, source.topology, options.max_sites)?;
    if !hex(expected_batch_revision) {
        return Err(bad("invalid expected catalog inspection batch revision"));
    }
    let current = catalog_digest::batch_revision(catalog, selection, source_role, source, options)?;
    if current != expected_batch_revision {
        return Err(bad(
            "source, catalog selection or inspection options changed; inspect again before choosing cuts",
        ));
    }
    context(&current, choice)
}
#[derive(Clone, Copy)]
struct Endpoint {
    top: usize,
    bottom: i64,
}
fn endpoint(
    boundary: &Boundary,
    r: &RestrictionInspection,
    source: Option<&crate::chemistry::DuplexSummary>,
) -> Result<Endpoint, Error> {
    match boundary {
        Boundary::SourceStart => Ok(Endpoint {
            top: 0,
            bottom: source.map_or(0, |s| s.bottom_start),
        }),
        Boundary::SourceEnd => Ok(Endpoint {
            top: r.source_length,
            bottom: match source {
                Some(s) => s
                    .bottom_start
                    .checked_add(i64::try_from(s.bottom_length).map_err(|_| Error::TooLarge)?)
                    .ok_or(Error::TooLarge)?,
                None => i64::try_from(r.source_length).map_err(|_| Error::TooLarge)?,
            },
        }),
        Boundary::Cut { site_index } => {
            let s = r
                .sites
                .get(*site_index)
                .filter(|s| s.index == *site_index)
                .ok_or_else(|| bad("selected cut instance is absent from this exact inspection"))?;
            if !s.certain || !s.diagnostics.is_empty() {
                return Err(bad(
                    "selected site is possible, outside substrate or has ambiguous strand geometry; no physical boundary is assigned",
                ));
            }
            let top = s
                .top_cut
                .ok_or_else(|| bad("selected site has no determined top-strand cut"))?;
            s.bottom_cut
                .ok_or_else(|| bad("selected site has no determined bottom-strand cut"))?;
            let stagger = s
                .bottom_boundary
                .checked_sub(s.top_boundary)
                .ok_or(Error::TooLarge)?;
            let bottom = i64::try_from(top)
                .map_err(|_| Error::TooLarge)?
                .checked_add(stagger)
                .ok_or(Error::TooLarge)?;
            Ok(Endpoint { top, bottom })
        }
    }
}
fn selected<'a>(
    inspection: &'a CatalogRestrictionInspection,
    choice: &Choice,
) -> Result<&'a FragmentSummary, Error> {
    selected_inner(inspection, choice, None)
}
fn selected_inner<'a>(
    inspection: &'a CatalogRestrictionInspection,
    choice: &Choice,
    source: Option<&crate::chemistry::DuplexSummary>,
) -> Result<&'a FragmentSummary, Error> {
    selected_raw_inner(&inspection.inspection, choice, source)
}
pub(crate) fn selected_fragment_raw_with_source<'a>(
    inspection: &'a crate::restriction_inspection::RestrictionInspection,
    choice: &Choice,
    source: &crate::chemistry::DuplexSummary,
) -> Result<&'a FragmentSummary, Error> {
    inspection.validate_with_source(source)?;
    source_geometry_preflight(inspection, source)?;
    selected_raw_inner(inspection, choice, Some(source))
}
fn selected_raw_inner<'a>(
    r: &'a crate::restriction_inspection::RestrictionInspection,
    choice: &Choice,
    source: Option<&crate::chemistry::DuplexSummary>,
) -> Result<&'a FragmentSummary, Error> {
    choice_shape(choice, r.topology, r.options.max_sites)?;
    let (mut first, mut second) = (
        endpoint(&choice.first, r, source)?,
        endpoint(&choice.second, r, source)?,
    );
    if !r.scan_complete || !r.products_determined || !r.limits.is_empty() {
        return Err(bad(
            "fresh recognition/product determination is incomplete; inspect diagnostics and bounds before selecting physical cuts",
        ));
    }
    let n = r.source_length;
    if n == 0 {
        return Err(bad("empty source"));
    }
    let (length, bottom_end) = if r.topology == Topology::Circular {
        if choice.complement {
            std::mem::swap(&mut first, &mut second);
        }
        if first.top >= n || second.top >= n {
            return Err(bad("circular cut lies outside one source turn"));
        }
        let wrap = second.top <= first.top;
        let end = second
            .top
            .checked_add(if wrap { n } else { 0 })
            .ok_or(Error::TooLarge)?;
        let bottom = second
            .bottom
            .checked_add(if wrap {
                i64::try_from(n).map_err(|_| Error::TooLarge)?
            } else {
                0
            })
            .ok_or(Error::TooLarge)?;
        (end - first.top, bottom)
    } else {
        if first.top > second.top {
            std::mem::swap(&mut first, &mut second);
        }
        if second.top > n || first.top == second.top {
            return Err(bad(
                "linear endpoints must delimit a nonempty source interval",
            ));
        }
        (second.top - first.top, second.bottom)
    };
    let bottom_length = bottom_end
        .checked_sub(first.bottom)
        .filter(|n| *n > 0)
        .and_then(|n| usize::try_from(n).ok())
        .ok_or_else(|| bad("selected boundaries do not delimit a determined bottom strand"))?;
    let mut matches = r.fragments.iter().filter(|f| {
        f.duplex.topology == Topology::Linear
            && f.top_start == first.top
            && f.duplex.top_length == length
            && f.bottom_start == first.bottom
            && f.duplex.bottom_length == bottom_length
    });
    let f=matches.next().ok_or_else(||bad("additional internal cuts split the requested path, or these endpoints do not identify one determined linear physical fragment; no sliced/merged reference substitute is produced"))?;
    if matches.next().is_some() {
        return Err(bad(
            "chosen endpoints match multiple physical fragment records; selection is ambiguous",
        ));
    }
    Ok(f)
}
/// Bounded coordinate consistency only; no bases, chemistry or material are
/// reconstructed from this documentary summary.
pub(crate) fn source_geometry_preflight(
    r: &RestrictionInspection,
    s: &crate::chemistry::DuplexSummary,
) -> Result<(), Error> {
    use crate::chemistry::{Polarity, Strand};
    if s.model.is_empty()
        || s.model.len() > 256
        || s.top_length != r.source_length
        || s.topology != r.topology
        || s.bottom_length == 0
        || s.bottom_length > collagene_model::MAX_SEQUENCE_LEN
        || s.nicks.len() > 10_000
    {
        return Err(bad("source duplex geometry differs from inspection"));
    }
    let n = i64::try_from(s.top_length).map_err(|_| Error::TooLarge)?;
    let end = s
        .bottom_start
        .checked_add(i64::try_from(s.bottom_length).map_err(|_| Error::TooLarge)?)
        .ok_or(Error::TooLarge)?;
    if s.bottom_start >= n
        || end <= 0
        || s.covalently_closed_in_model != (s.topology == Topology::Circular && s.nicks.is_empty())
    {
        return Err(bad("invalid source paired-strand contour"));
    }
    if s.topology == Topology::Circular {
        if s.bottom_start != 0 || s.bottom_length != s.top_length || s.ends.is_some() {
            return Err(bad("invalid circular source geometry"));
        }
    } else {
        let (left, right) = s
            .ends
            .as_ref()
            .ok_or_else(|| bad("linear source free-end geometry is unavailable"))?;
        let lp = if s.bottom_start > 0 {
            Polarity::FivePrime
        } else if s.bottom_start < 0 {
            Polarity::ThreePrime
        } else {
            Polarity::Blunt
        };
        let rp = if end > n {
            Polarity::FivePrime
        } else if end < n {
            Polarity::ThreePrime
        } else {
            Polarity::Blunt
        };
        for (end, polarity, length) in [
            (left, lp, s.bottom_start.unsigned_abs()),
            (right, rp, (end - n).unsigned_abs()),
        ] {
            if end.polarity != polarity
                || end.overhang.len() as u64 != length
                || end.overhang.len() > collagene_model::MAX_SEQUENCE_LEN
                || !end.overhang.bytes().all(|b| b"ACGT".contains(&b))
            {
                return Err(bad("source overhang geometry changed"));
            }
        }
    }
    let mut seen = std::collections::BTreeSet::new();
    for nick in &s.nicks {
        let (low, high) = match nick.strand {
            Strand::Top => (0, n),
            Strand::Bottom => (s.bottom_start, end),
        };
        let valid = if s.topology == Topology::Circular {
            nick.position >= 0 && nick.position < n
        } else {
            nick.position > low && nick.position < high
        };
        if !valid || !seen.insert((nick.strand == Strand::Bottom, nick.position)) {
            return Err(bad("invalid source nick coordinate"));
        }
    }
    Ok(())
}
/// Select using the actual source strand boundaries. This is a readonly geometry
/// lookup, never restoration of physical state from an imported summary.
pub fn selected_fragment_with_source<'a>(
    inspection: &'a CatalogRestrictionInspection,
    choice: &Choice,
    source: &crate::chemistry::DuplexSummary,
) -> Result<&'a FragmentSummary, Error> {
    inspection.validate()?;
    source_geometry_preflight(&inspection.inspection, source)?;
    selected_inner(inspection, choice, Some(source))
}
/// Readonly selection from an already bound inspection. No cleavage or product
/// materialization runs here, and imported evidence never acquires freshness.
pub fn selected_fragment<'a>(
    inspection: &'a CatalogRestrictionInspection,
    choice: &Choice,
) -> Result<&'a FragmentSummary, Error> {
    inspection.validate()?;
    selected(inspection, choice)
}
#[allow(clippy::too_many_arguments)]
pub fn preview(
    catalog: &catalog::Catalog,
    selection: &CatalogSelection,
    source_role: &str,
    source: &Molecule,
    options: &InspectionOptions,
    expected_batch_revision: &str,
    choice: &Choice,
) -> Result<Report, Error> {
    let context_revision = selection_revision(
        catalog,
        selection,
        source_role,
        source,
        options,
        expected_batch_revision,
        choice,
    )?;
    let inspection = catalog_digest::inspect(catalog, selection, source_role, source, options)?;
    if inspection.batch_revision != expected_batch_revision {
        return Err(bad(
            "fresh full inspection does not match the selected source/options/catalog batch",
        ));
    }
    let index = selected(&inspection, choice)?.index;
    let product = catalog_digest::product(catalog, selection, source_role, source, options, index)?;
    let report = Report {
        context_revision,
        choice: choice.clone(),
        inspection,
        product,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl Report {
    /// Readonly receipt checks. No digest rerun or restoration of physical state.
    pub fn validate(&self) -> Result<(), Error> {
        collagene_model::bounds::model_preflight(&[&self.product.product.molecule])?;
        self.inspection.validate()?;
        self.product.validate()?;
        if self.inspection.engine_version != catalog_digest::ENGINE_VERSION
            || self.product.engine_version != catalog_digest::ENGINE_VERSION
            || self.inspection.selection != self.product.selection
            || self.context_revision != context(&self.inspection.batch_revision, &self.choice)?
        {
            return Err(bad("catalog/choice context binding changed"));
        }
        let r = &self.inspection.inspection;
        let p = &self.product.product;
        if p.inspection_revision != r.batch_revision
            || p.source_role != r.source_role
            || p.source_revision != r.source_revision
            || p.options != r.options
            || p.geometry_revision != r.geometry_revision
            || p.options_revision != r.options_revision
        {
            return Err(bad(
                "product does not belong to the selected source and full inspection",
            ));
        }
        let actual = selected(&self.inspection, &self.choice)?;
        if actual.index != p.fragment_index
            || assembly_recipe::hash(actual)? != assembly_recipe::hash(&p.summary)?
        {
            return Err(bad(
                "returned product differs from the actual selected two-strand fragment",
            ));
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad(
                "receipt belongs to a different source/catalog/options/cut choice",
            ));
        }
        Ok(())
    }
}
