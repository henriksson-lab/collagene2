//! Independently implemented reference-only gel predictions and diagnostic digest ranking.
//! No observed images, empirical topology inference, external software or provider datasets.
use crate::{
    Alphabet, Error, Molecule, Topology,
    analysis::{self, Primer},
    chemistry::{self, Duplex},
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::fmt::Write;

const MAX_LANES: usize = 32;
const MAX_FRAGMENTS: usize = 1024;
const MAX_REPORT_BYTES: usize = 16 * 1024 * 1024;
const MODEL: &str = "Piecewise linear migration versus log10(mean strand length), calibration-bound linear dsDNA approximation v1";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CalibrationPoint {
    pub size_bp: f64,
    /// Fraction of the illustrated separation length measured from the well, 0..1.
    pub migration: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Calibration {
    pub label: String,
    pub points: Vec<CalibrationPoint>,
    /// Bands within this normalized migration span co-migrate under this declared model.
    pub resolution: f64,
}
impl Default for Calibration {
    fn default() -> Self {
        Self {
            label: "Illustrative log-size ladder; not an empirical gel calibration".into(),
            points: vec![
                CalibrationPoint {
                    size_bp: 100.,
                    migration: 0.9,
                },
                CalibrationPoint {
                    size_bp: 1000.,
                    migration: 0.5,
                },
                CalibrationPoint {
                    size_bp: 10000.,
                    migration: 0.1,
                },
            ],
            resolution: 0.02,
        }
    }
}
impl Calibration {
    pub fn validate(&self) -> Result<(), Error> {
        if self.label.len() > 1024
            || self.label.chars().any(char::is_control)
            || self.points.len() < 2
            || self.points.len() > 64
            || !self.resolution.is_finite()
            || !(0.001..=0.2).contains(&self.resolution)
        {
            return Err(Error::Unsupported(
                "gel calibration needs 2–64 ordered points and finite resolution 0.001–0.2".into(),
            ));
        }
        for p in &self.points {
            if !p.size_bp.is_finite()
                || !(1.0..=crate::model::MAX_SEQUENCE_LEN as f64).contains(&p.size_bp)
                || !p.migration.is_finite()
                || !(0.0..=1.0).contains(&p.migration)
            {
                return Err(Error::Unsupported(
                    "invalid gel calibration size or normalized migration".into(),
                ));
            }
        }
        if self
            .points
            .windows(2)
            .any(|p| p[0].size_bp >= p[1].size_bp || p[0].migration <= p[1].migration)
        {
            return Err(Error::Unsupported(
                "gel ladder sizes must increase while migration strictly decreases".into(),
            ));
        }
        Ok(())
    }
    pub fn position(&self, size_bp: f64) -> Result<Option<f64>, Error> {
        self.validate()?;
        if !size_bp.is_finite() || size_bp <= 0. {
            return Err(Error::Unsupported(
                "fragment size must be positive and finite".into(),
            ));
        }
        if size_bp < self.points[0].size_bp || size_bp > self.points.last().unwrap().size_bp {
            return Ok(None);
        }
        let upper = self.points.partition_point(|p| p.size_bp < size_bp);
        if self.points[upper].size_bp == size_bp {
            return Ok(Some(self.points[upper].migration));
        }
        let (a, b) = (&self.points[upper - 1], &self.points[upper]);
        // Relative log1p avoids cancellation for very close, but distinct, ladder sizes.
        let fraction = ((size_bp - a.size_bp) / a.size_bp).ln_1p()
            / ((b.size_bp - a.size_bp) / a.size_bp).ln_1p();
        Ok(Some(a.migration + (b.migration - a.migration) * fraction))
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum Preparation {
    Uncut,
    Digest {
        enzymes: Vec<String>,
    },
    Pcr {
        forward: Primer,
        reverse: Primer,
    },
    /// Full exact-event/catalog recipe. Saved reports never supply fresh bands.
    CatalogDigest {
        recipe: crate::digest_setups::Recipe,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct LaneInput {
    pub label: String,
    pub molecule: Molecule,
    pub preparation: Preparation,
    pub mass_ng: Option<f64>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MobilityState {
    LinearDuplex,
    CircularConformationUnknown,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Fragment {
    pub top_nt: usize,
    pub bottom_nt: usize,
    /// Mass-equivalent dsDNA base-pair count. Overhang-bearing fragments may be half-integral.
    pub size_bp: f64,
    pub state: MobilityState,
    pub migration: Option<f64>,
    pub amount_ng: Option<f64>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Band {
    pub fragment_indexes: Vec<usize>,
    pub migration: f64,
    pub min_size_bp: f64,
    pub max_size_bp: f64,
    pub migration_start: f64,
    pub migration_end: f64,
    pub amount_ng: Option<f64>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LaneReport {
    pub label: String,
    pub source_revision: String,
    pub preparation: String,
    pub fragments: Vec<Fragment>,
    pub bands: Vec<Band>,
    pub unpositioned_fragments: Vec<usize>,
    pub complete_under_model: bool,
    pub warnings: Vec<String>,
    pub error: Option<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GelReport {
    pub model: String,
    pub calibration: Calibration,
    pub lanes: Vec<LaneReport>,
    pub assumptions: Vec<String>,
}
fn error(message: &str) -> Error {
    Error::Unsupported(message.into())
}
fn enzyme_set(names: &[String]) -> Result<Vec<analysis::Enzyme>, Error> {
    if names.is_empty() || names.len() > 16 {
        return Err(error(
            "gel digest needs 1–16 distinct embedded enzyme names",
        ));
    }
    let catalog = analysis::starter_enzymes();
    let mut seen = BTreeSet::new();
    names
        .iter()
        .map(|name| {
            if !seen.insert(name) {
                return Err(error("duplicate digest enzyme"));
            }
            catalog
                .iter()
                .find(|e| e.name == *name)
                .cloned()
                .ok_or_else(|| error("unknown embedded digest enzyme"))
        })
        .collect()
}
fn groups(fragments: &[Fragment], resolution: f64) -> Vec<Band> {
    let mut positions: Vec<_> = fragments
        .iter()
        .enumerate()
        .filter_map(|(i, f)| f.migration.map(|p| (p, i)))
        .collect();
    positions.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut bands: Vec<Band> = vec![];
    for (position, index) in positions {
        let fragment = &fragments[index];
        if let Some(last) = bands
            .last_mut()
            .filter(|b| position - b.migration_start <= resolution + 1e-12)
        {
            last.fragment_indexes.push(index);
            last.migration_end = position;
            last.migration = (last.migration_start + position) / 2.;
            last.min_size_bp = last.min_size_bp.min(fragment.size_bp);
            last.max_size_bp = last.max_size_bp.max(fragment.size_bp);
            last.amount_ng = last.amount_ng.zip(fragment.amount_ng).map(|(a, b)| a + b);
        } else {
            bands.push(Band {
                fragment_indexes: vec![index],
                migration: position,
                migration_start: position,
                migration_end: position,
                min_size_bp: fragment.size_bp,
                max_size_bp: fragment.size_bp,
                amount_ng: fragment.amount_ng,
            });
        }
    }
    bands
}
fn lane(
    label: &str,
    molecule: &Molecule,
    preparation: &Preparation,
    mass_ng: Option<f64>,
    calibration: &Calibration,
) -> Result<LaneReport, Error> {
    molecule.validate()?;
    if label.len() > 256
        || label.chars().any(char::is_control)
        || mass_ng.is_some_and(|m| !m.is_finite() || m <= 0. || m > 1e9)
    {
        return Err(error("invalid gel lane label or loaded mass"));
    }
    let mut report = LaneReport {
        label: label.into(),
        source_revision: molecule.revision(),
        preparation: match preparation {
            Preparation::Uncut => "Uncut reference",
            Preparation::Digest { .. } => "Complete restriction digest",
            Preparation::Pcr { .. } => "Possible exact-PCR products",
            Preparation::CatalogDigest { .. } => "Catalog-pinned complete restriction digest",
        }
        .into(),
        fragments: vec![],
        bands: vec![],
        unpositioned_fragments: vec![],
        complete_under_model: true,
        warnings: vec![],
        error: None,
    };
    let result = (|| -> Result<(), Error> {
        let mut lengths: Vec<(usize, usize, Topology)> = vec![];
        match preparation {
            Preparation::Uncut => {
                let duplex = Duplex::from_molecule(molecule)?;
                lengths.push((duplex.top.len(), duplex.bottom.len(), duplex.topology));
            }
            Preparation::Digest { enzymes } => {
                let duplex = Duplex::from_molecule(molecule)?;
                let fragments = chemistry::digest(&duplex, &enzyme_set(enzymes)?)?;
                if fragments.len() > MAX_FRAGMENTS {
                    return Err(Error::TooLarge);
                }
                for f in fragments {
                    lengths.push((f.duplex.top.len(), f.duplex.bottom.len(), f.duplex.topology));
                }
                report.preparation = format!("Complete digest: {}", enzymes.join(" + "));
                report.warnings.push("Assumes complete cleavage; methylation, star activity and buffer/temperature effects are not evaluated.".into());
            }
            Preparation::Pcr { forward, reverse } => {
                if forward.name.len() > 256 || reverse.name.len() > 256 {
                    return Err(error("gel PCR primer names exceed 256 bytes"));
                }
                if chemistry::has_unmapped_physical_evidence(molecule) {
                    return Err(error(
                        "PCR of unmodeled physical-intermediate evidence is not a gel reference",
                    ));
                }
                let products = analysis::simulate_pcr(molecule, forward, reverse)?;
                if products.len() > MAX_FRAGMENTS {
                    return Err(Error::TooLarge);
                }
                report.complete_under_model = products.len() <= 1;
                for p in products {
                    lengths.push((p.sequence.len(), p.sequence.len(), Topology::Linear));
                }
                report.warnings.push("These are possible exact-PCR products, not predicted yields or a claim that every alternative accumulates. Primer roles do not exclude competing orientations or single-oligo products.".into());
                if lengths.len() != 1 && mass_ng.is_some() {
                    report.warnings.push("Loaded mass cannot be allocated across PCR alternatives; band masses remain unknown.".into());
                }
            }
            Preparation::CatalogDigest { recipe } => {
                recipe.validate()?;
                if recipe.source.revision != report.source_revision {
                    return Err(error(
                        "Pinned digest source revision differs from the actual lane molecule; resolve its exact retained source event.",
                    ));
                }
                let catalog = crate::catalog_digest::builtin_catalog()?;
                recipe.selection.validate_for(&catalog)?;
                let fresh = crate::catalog_digest::inspect(
                    &catalog,
                    &recipe.selection,
                    &crate::restriction_history::source_role(&recipe.source)?,
                    molecule,
                    &recipe.options,
                )?;
                report.warnings.push("Assumes complete ideal cleavage; methylation remains unassessed. Supplier conditions, partial digestion, star activity and experimental cleavage are not evaluated.".into());
                report
                    .warnings
                    .extend(fresh.inspection.diagnostics.iter().cloned());
                if !fresh.inspection.scan_complete || !fresh.inspection.products_determined {
                    return Err(error(
                        "Catalog digest products are withheld: the selected inspection limits or unresolved recognition/strand/physical evidence do not establish a complete fragment set.",
                    ));
                }
                if fresh.inspection.fragments.len() > MAX_FRAGMENTS {
                    return Err(error(
                        "Catalog digest exceeds the gel limit of 1024 fragments; no fragments are omitted and the selected inspection limits were not changed.",
                    ));
                }
                report.preparation = format!(
                    "Catalog digest: {} · edition {}",
                    recipe
                        .selection
                        .enzymes
                        .iter()
                        .map(|choice| choice.enzyme.name.as_str())
                        .collect::<Vec<_>>()
                        .join(" + "),
                    recipe.selection.catalog_edition
                );
                for fragment in fresh.inspection.fragments {
                    let duplex = fragment.duplex;
                    lengths.push((duplex.top_length, duplex.bottom_length, duplex.topology));
                }
            }
        }
        let total_nt = lengths
            .iter()
            .try_fold(0usize, |n, (a, b, _)| {
                n.checked_add(*a).and_then(|n| n.checked_add(*b))
            })
            .ok_or(Error::TooLarge)?;
        let amounts_known = !matches!(preparation, Preparation::Pcr { .. }) || lengths.len() == 1;
        for (a, b, topology) in lengths {
            let size_bp = (a + b) as f64 / 2.;
            let state = if topology == Topology::Circular {
                MobilityState::CircularConformationUnknown
            } else {
                MobilityState::LinearDuplex
            };
            let migration = if state == MobilityState::LinearDuplex {
                calibration.position(size_bp)?
            } else {
                None
            };
            if a != b {
                report.warnings.push("Staggered-end fragment uses mean strand length as an approximate linear-dsDNA mass surrogate; overhang-specific mobility is not modeled.".into());
            }
            if state == MobilityState::CircularConformationUnknown {
                report.warnings.push("Uncut circular conformation is unknown: nicked/open-circular, relaxed and supercoiled species do not receive a linear-DNA migration position.".into());
            }
            if migration.is_none() && state == MobilityState::LinearDuplex {
                report.warnings.push("A linear fragment is outside the supplied ladder interval; no migration extrapolation is asserted.".into());
            }
            let amount_ng = if amounts_known {
                mass_ng.map(|m| m * (a + b) as f64 / total_nt as f64)
            } else {
                None
            };
            report.fragments.push(Fragment {
                top_nt: a,
                bottom_nt: b,
                size_bp,
                state,
                migration,
                amount_ng,
            });
        }
        Ok(())
    })();
    if let Err(e) = result {
        report.complete_under_model = false;
        report.error = Some(e.to_string());
        report.fragments.clear();
    }
    report.warnings.sort();
    report.warnings.dedup();
    report.unpositioned_fragments = report
        .fragments
        .iter()
        .enumerate()
        .filter_map(|(i, f)| f.migration.is_none().then_some(i))
        .collect();
    report.bands = groups(&report.fragments, calibration.resolution);
    Ok(report)
}
fn assumptions() -> Vec<String> {
    vec![MODEL.into(),"Positions are normalized illustrative/calibration coordinates, not a gel percentage, voltage, run time or kinetic prediction. No model extrapolation outside the supplied ladder interval.".into(),"Band grouping uses a maximum migration span (complete-link threshold), not transitive single-link chaining. Every modeled fragment is retained in the lane table.".into(),"Mass allocation for complete digests assumes equimolar copies from each source molecule and fluorescence proportional to DNA mass; no detection threshold or exposure model is assigned. Unknown amounts use uniform illustrative strokes.".into(),"Circular conformation, incompletely cleaved mixtures, protein/RNA and single-strand mobilities are not inferred. This tool does not compare observed gels or integrate OpenGel.".into()]
}
fn bound_to<T: Serialize + ?Sized>(value: &T, limit: usize) -> Result<(), Error> {
    struct Counter {
        count: usize,
        limit: usize,
    }
    impl std::io::Write for Counter {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.count = self
                .count
                .checked_add(b.len())
                .filter(|n| *n <= self.limit)
                .ok_or_else(|| std::io::Error::other("gel report bound"))?;
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(&mut Counter { count: 0, limit }, value).map_err(|_| Error::TooLarge)
}
fn bound<T: Serialize + ?Sized>(value: &T) -> Result<(), Error> {
    bound_to(value, MAX_REPORT_BYTES)
}
pub fn predict_gel(inputs: &[LaneInput], calibration: &Calibration) -> Result<GelReport, Error> {
    calibration.validate()?;
    if inputs.is_empty() || inputs.len() > MAX_LANES {
        return Err(error("gel needs 1–32 source lanes"));
    }
    for input in inputs {
        if let Preparation::CatalogDigest { recipe } = &input.preparation {
            crate::model::bounds::model_preflight(&[&input.molecule])?;
            recipe.validate()?;
        }
        input.molecule.validate()?;
    }
    bound_to(inputs, crate::formats::MAX_IMPORT_BYTES)?;
    let bases = inputs
        .iter()
        .try_fold(0usize, |n, i| n.checked_add(i.molecule.sequence.len()))
        .ok_or(Error::TooLarge)?;
    if bases > crate::formats::MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    let lanes = inputs
        .iter()
        .map(|input| {
            lane(
                &input.label,
                &input.molecule,
                &input.preparation,
                input.mass_ng,
                calibration,
            )
        })
        .collect::<Result<_, _>>()?;
    let result = GelReport {
        model: MODEL.into(),
        calibration: calibration.clone(),
        lanes,
        assumptions: assumptions(),
    };
    bound(&result)?;
    Ok(result)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct DiagnosticOptions {
    pub enzymes: Vec<String>,
    pub include_double_digests: bool,
    pub max_results: usize,
    pub max_evaluations: usize,
    pub max_total_bases: usize,
    pub calibration: Calibration,
}
impl Default for DiagnosticOptions {
    fn default() -> Self {
        Self {
            enzymes: analysis::starter_enzymes()
                .into_iter()
                .map(|e| e.name)
                .collect(),
            include_double_digests: true,
            max_results: 10,
            max_evaluations: 128,
            max_total_bases: 50_000_000,
            calibration: Calibration::default(),
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DigestCandidate {
    pub enzymes: Vec<String>,
    pub profiles: Vec<LaneReport>,
    pub distinguished_pairs: Vec<(usize, usize)>,
    pub indistinguishable_pairs: Vec<(usize, usize)>,
    pub unknown_pairs: Vec<(usize, usize)>,
    pub minimum_diagnostic_separation: Option<f64>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiagnosticReport {
    pub source_revisions: Vec<String>,
    pub options: DiagnosticOptions,
    pub candidates: Vec<DigestCandidate>,
    pub evaluated_digests: usize,
    pub work_bases: usize,
    pub exhaustive_search: bool,
    pub warnings: Vec<String>,
}
fn distinguish(a: &LaneReport, b: &LaneReport, resolution: f64) -> Option<(bool, f64)> {
    if !a.complete_under_model
        || !b.complete_under_model
        || a.error.is_some()
        || b.error.is_some()
        || !a.unpositioned_fragments.is_empty()
        || !b.unpositioned_fragments.is_empty()
    {
        return None;
    }
    if a.bands.len() != b.bands.len() {
        return Some((true, 1.));
    }
    let separation = a
        .bands
        .iter()
        .zip(&b.bands)
        .map(|(a, b)| (a.migration - b.migration).abs())
        .fold(0., f64::max);
    Some((separation > resolution + 1e-12, separation))
}
pub fn diagnostic_digests(
    molecules: &[Molecule],
    options: &DiagnosticOptions,
) -> Result<DiagnosticReport, Error> {
    options.calibration.validate()?;
    if !(2..=16).contains(&molecules.len())
        || options.max_results == 0
        || options.max_results > 20
        || options.max_evaluations == 0
        || options.max_evaluations > 256
        || options.max_total_bases == 0
        || options.max_total_bases > 200_000_000
    {
        return Err(error("invalid diagnostic candidate/output/work limits"));
    }
    enzyme_set(&options.enzymes)?;
    for m in molecules {
        m.validate()?;
        if m.sequence.alphabet() != Alphabet::Dna {
            return Err(error("diagnostic digest requires DNA documents"));
        }
    }
    bound_to(molecules, crate::formats::MAX_IMPORT_BYTES)?;
    let total_bases = molecules
        .iter()
        .try_fold(0usize, |n, m| n.checked_add(m.sequence.len()))
        .ok_or(Error::TooLarge)?;
    let mut names = options.enzymes.clone();
    names.sort();
    let mut sets: Vec<Vec<String>> = names.iter().map(|n| vec![n.clone()]).collect();
    if options.include_double_digests {
        for i in 0..names.len() {
            for j in i + 1..names.len() {
                sets.push(vec![names[i].clone(), names[j].clone()]);
            }
        }
    }
    let mut report=DiagnosticReport{source_revisions:molecules.iter().map(Molecule::revision).collect(),options:options.clone(),candidates:vec![],evaluated_digests:0,work_bases:0,exhaustive_search:true,warnings:vec!["Discrimination concerns modeled band positions/counts, not experimental sensitivity, staining intensity or sequence uniqueness. Unknown/unpositioned profiles never count as negative or positive identification evidence.".into()]};
    for enzymes in sets {
        let work = total_bases
            .checked_mul(enzymes.len())
            .ok_or(Error::TooLarge)?;
        if report.evaluated_digests >= options.max_evaluations
            || report
                .work_bases
                .checked_add(work)
                .is_none_or(|n| n > options.max_total_bases)
        {
            report.exhaustive_search = false;
            break;
        }
        report.work_bases += work;
        report.evaluated_digests += 1;
        let preparation = Preparation::Digest {
            enzymes: enzymes.clone(),
        };
        let profiles = molecules
            .iter()
            .map(|m| {
                lane(
                    &m.name.chars().take(60).collect::<String>(),
                    m,
                    &preparation,
                    None,
                    &options.calibration,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut candidate = DigestCandidate {
            enzymes,
            profiles,
            distinguished_pairs: vec![],
            indistinguishable_pairs: vec![],
            unknown_pairs: vec![],
            minimum_diagnostic_separation: None,
        };
        for a in 0..molecules.len() {
            for b in a + 1..molecules.len() {
                match distinguish(
                    &candidate.profiles[a],
                    &candidate.profiles[b],
                    options.calibration.resolution,
                ) {
                    Some((true, separation)) => {
                        candidate.distinguished_pairs.push((a, b));
                        candidate.minimum_diagnostic_separation = Some(
                            candidate
                                .minimum_diagnostic_separation
                                .map_or(separation, |s| s.min(separation)),
                        );
                    }
                    Some((false, _)) => candidate.indistinguishable_pairs.push((a, b)),
                    None => candidate.unknown_pairs.push((a, b)),
                }
            }
        }
        report.candidates.push(candidate);
        report.candidates.sort_by(|a, b| {
            b.distinguished_pairs
                .len()
                .cmp(&a.distinguished_pairs.len())
                .then(a.unknown_pairs.len().cmp(&b.unknown_pairs.len()))
                .then_with(|| {
                    b.minimum_diagnostic_separation
                        .unwrap_or(0.)
                        .total_cmp(&a.minimum_diagnostic_separation.unwrap_or(0.))
                })
                .then(a.enzymes.len().cmp(&b.enzymes.len()))
                .then(a.enzymes.cmp(&b.enzymes))
        });
        report.candidates.truncate(options.max_results);
    }
    if !report.exhaustive_search {
        report.warnings.push("Search budget exhausted; retained candidates are best among evaluated enzyme sets, not proven globally optimal.".into());
    }
    bound(&report)?;
    Ok(report)
}
fn esc(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}
/// Vector diagram generated from a bounded report. Every unknown/out-of-range fragment is
/// called out in text; no arbitrary circular band is drawn as if it were calibrated linear DNA.
pub fn render_svg(report: &GelReport) -> Result<String, Error> {
    report.calibration.validate()?;
    bound(report)?;
    if report.lanes.is_empty() || report.lanes.len() > MAX_LANES {
        return Err(Error::TooLarge);
    }
    let width = (160 + report.lanes.len() * 90).max(920);
    let height = 680;
    let mut out = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\" role=\"img\" aria-label=\"Predicted linear DNA gel with ladder and unknown fragments listed\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/><text x=\"18\" y=\"28\" font-family=\"sans-serif\" font-size=\"16\">Collagene · predicted gel</text><text x=\"18\" y=\"49\" font-family=\"sans-serif\" font-size=\"10\">Illustrative/calibration model; no observed-gel comparison</text><rect x=\"50\" y=\"80\" width=\"{}\" height=\"450\" fill=\"#132a27\"/><g fill=\"white\" font-family=\"sans-serif\" font-size=\"9\">",
        width - 65
    );
    let metadata=serde_json::to_string(&serde_json::json!({"model":report.model,"calibration":report.calibration,"source_revisions":report.lanes.iter().map(|l|l.source_revision.as_str()).collect::<Vec<_>>(),"assumptions":report.assumptions})).map_err(|e|Error::Format(e.to_string()))?;
    let _ = write!(
        out,
        "<metadata id=\"collagene-gel-evidence\">{}</metadata>",
        esc(&metadata)
    );
    let _ = write!(out, "<text x=\"70\" y=\"101\">Ladder</text>");
    for point in &report.calibration.points {
        let y = 110. + point.migration * 400.;
        let _ = write!(
            out,
            "<path d=\"M70 {y:.3}h40\" stroke=\"#b5eed6\" stroke-width=\"3\"/><text x=\"5\" y=\"{:.3}\" fill=\"#173b30\">{:.0} bp</text>",
            y + 3.,
            point.size_bp
        );
    }
    for (i, lane) in report.lanes.iter().enumerate() {
        if lane.label.len() > 256
            || lane.fragments.len() > MAX_FRAGMENTS
            || lane.bands.len() > MAX_FRAGMENTS
        {
            return Err(Error::TooLarge);
        }
        if lane.fragments.iter().any(|f| {
            !f.size_bp.is_finite()
                || f.size_bp <= 0.
                || f.migration
                    .is_some_and(|p| !p.is_finite() || !(0.0..=1.0).contains(&p))
                || f.amount_ng.is_some_and(|m| !m.is_finite() || m < 0.)
        }) {
            return Err(error("invalid serialized gel fragment"));
        }
        let x = 150 + i * 90;
        let title = esc(&lane.label);
        let short = esc(&lane.label.chars().take(12).collect::<String>());
        let _ = write!(
            out,
            "<g data-source-revision=\"{}\"><title>{title} · {}</title><text x=\"{x}\" y=\"101\">{short}</text><rect x=\"{x}\" y=\"80\" width=\"45\" height=\"8\" fill=\"#e9f5e9\"/>",
            esc(&lane.source_revision),
            esc(&lane.preparation)
        );
        for band in &lane.bands {
            if !band.migration.is_finite()
                || !(0.0..=1.0).contains(&band.migration)
                || band.amount_ng.is_some_and(|v| !v.is_finite() || v < 0.)
                || band
                    .fragment_indexes
                    .iter()
                    .any(|i| *i >= lane.fragments.len())
                || !band.min_size_bp.is_finite()
                || !band.max_size_bp.is_finite()
                || band.min_size_bp <= 0.
                || band.min_size_bp > band.max_size_bp
            {
                return Err(error("invalid serialized gel band"));
            }
            let y = 110. + band.migration * 400.;
            let _ = write!(
                out,
                "<path d=\"M{x} {y:.3}h45\" stroke=\"#b5eed6\" stroke-width=\"3\"><title>{:.2}–{:.2} bp; {} source fragments; mass {}</title></path>",
                band.min_size_bp,
                band.max_size_bp,
                band.fragment_indexes.len(),
                band.amount_ng
                    .map(|m| format!("{m:.3} ng"))
                    .unwrap_or_else(|| "unknown".into())
            );
        }
        let unknown = lane.unpositioned_fragments.len();
        let _ = write!(
            out,
            "<text x=\"{x}\" y=\"552\" fill=\"#203d32\">{} bands</text><text x=\"{x}\" y=\"568\" fill=\"#8a3e28\">{unknown} unplaced</text></g>",
            lane.bands.len()
        );
        if lane.error.is_some() {
            let _ = write!(
                out,
                "<text x=\"{x}\" y=\"585\" fill=\"#a22\">unsupported</text>"
            );
        }
    }
    out.push_str("</g><g font-family=\"sans-serif\" font-size=\"11\" fill=\"#203d32\"><text x=\"18\" y=\"615\">Unplaced: unknown circular conformation or size outside ladder. Inspect the lane table/report.</text><text x=\"18\" y=\"635\">Co-migration preserves member fragments. Uniform strokes do not predict staining intensity.</text><text x=\"18\" y=\"655\">Linear dsDNA, complete cleavage and explicit resolution assumptions only.</text></g></svg>");
    if out.len() > MAX_REPORT_BYTES {
        return Err(Error::TooLarge);
    }
    Ok(out)
}
