//! Revision-bound, bounded reference browser tracks. No viewport operation clones the source.
use crate::translation::{self, InitiationPolicy, TranslationOptions};
use collagene_model::{Alphabet, Error, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

pub const TRACK_ENGINE_VERSION: &str = "browser-tracks-1";
pub const MAX_TRACK_ITEMS: usize = 1000;
pub const MAX_TRACK_RESIDUES: usize = 10000;
pub const MAX_TRACK_WINDOW: usize = 100000;
pub const MAX_TRACK_WORK: usize = 1000000;
pub const MAX_TRACK_PARTS: usize = 10000;
pub const MAX_TRACK_WARNINGS: usize = 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct TrackWindow {
    pub start: usize,
    pub end: usize,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum TrackOrientation {
    Forward,
    Reverse,
    Both,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum KnownCdsPolicy {
    Off,
    Validated,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct TrackOptions {
    pub features: bool,
    pub primers: bool,
    pub enzyme_recognition: bool,
    pub enzyme_cuts: bool,
    pub translations: bool,
    pub orfs: bool,
    pub enzyme_names: Vec<String>,
    pub genetic_code: u8,
    pub frames: Vec<u8>,
    pub orientation: TrackOrientation,
    pub known_cds: KnownCdsPolicy,
    pub max_items: usize,
    pub max_residues: usize,
    pub max_window_bases: usize,
    pub max_work_units: usize,
    pub min_orf_amino_acids: usize,
}
impl Default for TrackOptions {
    fn default() -> Self {
        Self {
            features: true,
            primers: true,
            enzyme_recognition: false,
            enzyme_cuts: false,
            translations: false,
            orfs: false,
            enzyme_names: vec!["EcoRI".into(), "BamHI".into(), "BsaI".into()],
            genetic_code: 1,
            frames: vec![0, 1, 2],
            orientation: TrackOrientation::Both,
            known_cds: KnownCdsPolicy::Off,
            max_items: 1000,
            max_residues: 4000,
            max_window_bases: 10000,
            max_work_units: 100000,
            min_orf_amino_acids: 20,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum TrackKind {
    Feature,
    Primer,
    Enzyme,
    Translation,
    Cds,
    Orf,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ReferencePart {
    pub start: usize,
    pub end: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ResidueMapping {
    /// Zero-based residue index in the full oriented frame/CDS/ORF, never the viewport.
    pub residue_index: usize,
    pub amino_acid: String,
    /// Codon bases in biological order; reverse orientation parts have ascending endpoints.
    pub codon_parts: Vec<ReferencePart>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CutPositions {
    /// Reference strand boundaries, not recognition starts. Missing means outside a linear molecule.
    pub top: Option<usize>,
    pub bottom: Option<usize>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct TrackItem {
    pub id: String,
    pub kind: TrackKind,
    pub label: String,
    pub parts: Vec<ReferencePart>,
    /// None for fuzzy, unordered, or mixed-direction annotation evidence.
    pub reverse: Option<bool>,
    pub source_feature: Option<usize>,
    pub enzyme: Option<String>,
    pub certain: bool,
    pub cuts: Option<CutPositions>,
    pub frame: Option<u8>,
    /// Visible residues only; each residue retains its full-source index and codon coordinates.
    pub protein: Option<String>,
    pub residues: Vec<ResidueMapping>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct TrackWork {
    pub used_units: usize,
    pub limit_units: usize,
    pub emitted_residues: usize,
    /// Preparation cost, performed once per immutable source, outside viewport budget.
    pub source_bases: usize,
    pub indexed_spans: usize,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct TrackReport {
    pub source_revision: String,
    pub engine: String,
    pub version: String,
    pub window: TrackWindow,
    pub options: TrackOptions,
    /// False when any configured scientific track was omitted or capped. Disabled tracks do not count.
    pub complete: bool,
    pub work: TrackWork,
    pub warnings: Vec<String>,
    pub items: Vec<TrackItem>,
}
#[derive(Debug)]
struct Span {
    start: usize,
    end: usize,
    feature: usize,
}
#[derive(Default)]
struct IndexState {
    reverse: Option<bool>,
    mixed: bool,
    length: usize,
}
fn index_location(
    location: &Location,
    reverse: bool,
    feature: usize,
    spans: &mut Vec<Span>,
    state: &mut IndexState,
) -> Result<(), Error> {
    match location {
        Location::Range { start, end, .. } => {
            if spans.len() >= crate::tracks::MAX_INDEX_SPANS {
                return Err(Error::TooLarge);
            }
            state.mixed |= state.reverse.is_some_and(|r| r != reverse);
            state.reverse.get_or_insert(reverse);
            state.length = state.length.saturating_add(end - start);
            spans.push(Span {
                start: *start,
                end: *end,
                feature,
            });
        }
        Location::Complement(inner) => index_location(inner, !reverse, feature, spans, state)?,
        Location::Join(parts) | Location::Order(parts) => {
            for part in parts {
                index_location(part, reverse, feature, spans, state)?;
            }
        }
        Location::Unresolved(_) => {}
    }
    Ok(())
}
/// Owns an immutable shared document, so cached revision/index cannot become stale.
/// Construction validates and hashes the complete source once. Clone the Arc, not the chromosome.
#[derive(Debug)]
pub struct PreparedTracks {
    molecule: Arc<Molecule>,
    revision: String,
    spans: Vec<Span>,
    prefix_end: Vec<usize>,
    feature_meta: Vec<(Option<bool>, usize)>,
    feature_labels: Vec<String>,
    qualifier_work: Vec<usize>,
    unplaced_features: usize,
    unplaced_primers: usize,
    unplaced_cds: usize,
}
impl PreparedTracks {
    pub fn new(molecule: Arc<Molecule>) -> Result<Self, Error> {
        molecule.validate()?;
        let revision = molecule.revision();
        let mut spans = Vec::new();
        let (mut unplaced_features, mut unplaced_primers, mut unplaced_cds) = (0, 0, 0);
        let mut feature_meta = vec![(None, 0); molecule.features.len()];
        let mut feature_labels = vec![String::new(); molecule.features.len()];
        let mut qualifier_work = vec![0; molecule.features.len()];
        for (feature, f) in molecule
            .features
            .iter()
            .enumerate()
            .filter(|(_, f)| !f.kind.eq_ignore_ascii_case("source"))
        {
            let mut state = IndexState::default();
            index_location(&f.location, false, feature, &mut spans, &mut state)?;
            if state.length == 0 {
                if f.kind.eq_ignore_ascii_case("primer_bind") {
                    unplaced_primers += 1;
                } else {
                    unplaced_features += 1;
                }
                if f.kind.eq_ignore_ascii_case("CDS") {
                    unplaced_cds += 1;
                }
            }
            feature_meta[feature] = (
                if !state.mixed && crate::tracks::terminal_span(&f.location).is_some() {
                    state.reverse
                } else {
                    None
                },
                state.length,
            );
            feature_labels[feature] = f.label().chars().take(512).collect();
            qualifier_work[feature] = f.qualifiers.iter().fold(0usize, |sum, (key, value)| {
                sum.saturating_add(key.len())
                    .saturating_add(value.as_ref().map_or(0, String::len))
                    .saturating_add(1)
            });
        }
        spans.sort_by_key(|s| (s.start, s.end, s.feature));
        let mut end = 0;
        let prefix_end = spans
            .iter()
            .map(|s| {
                end = end.max(s.end);
                end
            })
            .collect();
        Ok(Self {
            molecule,
            revision,
            spans,
            prefix_end,
            feature_meta,
            feature_labels,
            qualifier_work,
            unplaced_features,
            unplaced_primers,
            unplaced_cds,
        })
    }
    pub fn source_revision(&self) -> &str {
        &self.revision
    }
    pub fn report(
        &self,
        window: TrackWindow,
        options: &TrackOptions,
    ) -> Result<TrackReport, Error> {
        self.report_for_revision(&self.revision, window, options)
    }
    /// Explicit revision check for asynchronous callers; a stale request is an error, never a partial report.
    pub fn report_for_revision(
        &self,
        revision: &str,
        window: TrackWindow,
        options: &TrackOptions,
    ) -> Result<TrackReport, Error> {
        if revision != self.revision {
            return Err(Error::Format("stale browser-track source revision".into()));
        }
        validate_options(&self.molecule, window, options)?;
        let mut out = TrackReport {
            source_revision: self.revision.clone(),
            engine: "collagene-core".into(),
            version: TRACK_ENGINE_VERSION.into(),
            window,
            options: options.clone(),
            complete: true,
            work: TrackWork {
                used_units: 0,
                limit_units: options.max_work_units,
                emitted_residues: 0,
                source_bases: self.molecule.sequence.len(),
                indexed_spans: self.spans.len(),
            },
            warnings: vec![],
            items: vec![],
        };
        if (options.features && self.unplaced_features > 0)
            || (options.primers && self.unplaced_primers > 0)
            || (options.known_cds == KnownCdsPolicy::Validated && self.unplaced_cds > 0)
        {
            warn(
                &mut out,
                "Some enabled annotations have unresolved locations and cannot be placed in any reference viewport.",
            );
        }
        self.annotations(&mut out);
        if options.enzyme_recognition || options.enzyme_cuts {
            self.enzymes(&mut out)?;
        }
        if options.translations {
            self.frames(&mut out)?;
        }
        if options.known_cds == KnownCdsPolicy::Validated {
            self.cds(&mut out)?;
        }
        if options.orfs {
            self.orfs(&mut out)?;
        }
        if out.warnings.len() > MAX_TRACK_WARNINGS {
            out.warnings.truncate(MAX_TRACK_WARNINGS - 1);
            out.warnings
                .push("Additional track warnings were omitted.".into());
            out.complete = false;
        }
        Ok(out)
    }
}
/// One-shot convenience; repeatedly scrolling callers should retain PreparedTracks instead.
pub fn view_tracks(
    molecule: &Molecule,
    window: TrackWindow,
    options: &TrackOptions,
) -> Result<TrackReport, Error> {
    PreparedTracks::new(Arc::new(molecule.clone()))?.report(window, options)
}
fn validate_options(m: &Molecule, w: TrackWindow, o: &TrackOptions) -> Result<(), Error> {
    let n = m.sequence.len();
    if w.start >= n
        || w.end <= w.start
        || w.end - w.start > n
        || (m.topology == Topology::Linear && w.end > n)
        || w.end > n.saturating_mul(2)
    {
        return Err(Error::Coordinates(
            "track window must be a nonempty reference interval of at most one turn".into(),
        ));
    }
    if o.max_items == 0
        || o.max_items > MAX_TRACK_ITEMS
        || o.max_residues == 0
        || o.max_residues > MAX_TRACK_RESIDUES
        || o.max_window_bases == 0
        || o.max_window_bases > MAX_TRACK_WINDOW
        || w.end - w.start > o.max_window_bases
        || o.max_work_units == 0
        || o.max_work_units > MAX_TRACK_WORK
        || o.frames.len() > 3
        || o.frames.iter().any(|f| *f > 2)
        || o.frames.iter().collect::<BTreeSet<_>>().len() != o.frames.len()
        || o.enzyme_names.len() > 12
        || o.enzyme_names.iter().collect::<BTreeSet<_>>().len() != o.enzyme_names.len()
        || o.min_orf_amino_acids > collagene_model::MAX_SEQUENCE_LEN / 3
    {
        return Err(Error::Format(
            "invalid browser-track configuration or budget".into(),
        ));
    }
    if !translation::genetic_codes()
        .iter()
        .any(|c| c.id == o.genetic_code)
    {
        return Err(Error::Unsupported("unsupported track genetic code".into()));
    }
    let enzymes = crate::analysis::starter_enzymes();
    if o.enzyme_names
        .iter()
        .any(|name| !enzymes.iter().any(|e| &e.name == name))
    {
        return Err(Error::Unsupported(
            "track enzyme is absent from the factual starter catalog".into(),
        ));
    }
    if m.sequence.alphabet() == Alphabet::Protein
        && (o.translations || o.orfs || o.known_cds != KnownCdsPolicy::Off)
    {
        return Err(Error::Unsupported(
            "nucleotide translation tracks require DNA or RNA".into(),
        ));
    }
    if m.sequence.alphabet() != Alphabet::Dna && (o.enzyme_recognition || o.enzyme_cuts) {
        return Err(Error::Unsupported("restriction tracks require DNA".into()));
    }
    Ok(())
}
fn warn(out: &mut TrackReport, message: &str) {
    out.complete = false;
    if out.warnings.len() >= MAX_TRACK_WARNINGS {
        if let Some(last) = out.warnings.last_mut() {
            *last = "Additional track warnings were omitted.".into();
        }
        return;
    }
    if !out.warnings.iter().any(|s| s == message) {
        out.warnings.push(message.into());
    }
}
fn spend(out: &mut TrackReport, units: usize) -> bool {
    if units > out.work.limit_units - out.work.used_units {
        warn(
            out,
            "Track work budget reached; omitted results are not evidence of absence.",
        );
        false
    } else {
        out.work.used_units += units;
        true
    }
}
fn emit(out: &mut TrackReport, item: TrackItem) -> bool {
    if out.items.iter().map(|t| t.parts.len()).sum::<usize>() + item.parts.len() > MAX_TRACK_PARTS {
        out.work.emitted_residues -= item.residues.len();
        warn(
            out,
            "Track reference-part limit reached; additional results are omitted.",
        );
        return false;
    }
    if out.items.len() >= out.options.max_items {
        out.work.emitted_residues -= item.residues.len();
        warn(
            out,
            "Track item limit reached; additional results are omitted.",
        );
        false
    } else {
        out.items.push(item);
        true
    }
}
fn item(
    id: String,
    kind: TrackKind,
    label: String,
    parts: Vec<ReferencePart>,
    reverse: Option<bool>,
) -> TrackItem {
    TrackItem {
        id,
        kind,
        label,
        parts,
        reverse,
        source_feature: None,
        enzyme: None,
        certain: true,
        cuts: None,
        frame: None,
        protein: None,
        residues: vec![],
    }
}
fn window_parts(w: TrackWindow, n: usize) -> Vec<ReferencePart> {
    if w.end <= n {
        vec![ReferencePart {
            start: w.start,
            end: w.end,
        }]
    } else {
        vec![
            ReferencePart {
                start: w.start,
                end: n,
            },
            ReferencePart {
                start: 0,
                end: w.end - n,
            },
        ]
    }
}
fn overlaps(parts: &[ReferencePart], w: TrackWindow, n: usize) -> bool {
    window_parts(w, n)
        .iter()
        .any(|v| parts.iter().any(|p| p.start < v.end && p.end > v.start))
}
fn visible_boundary(pos: usize, w: TrackWindow, n: usize) -> bool {
    window_parts(w, n)
        .iter()
        .any(|p| pos >= p.start && pos <= p.end)
}
fn location_parts(l: &Location) -> Vec<ReferencePart> {
    l.spans()
        .into_iter()
        .map(|(start, end, _)| ReferencePart { start, end })
        .collect()
}
fn direction(l: &Location) -> Option<bool> {
    crate::tracks::terminal_span(l)?;
    let spans = l.spans();
    let reverse = spans.first()?.2;
    spans.iter().all(|s| s.2 == reverse).then_some(reverse)
}
fn allowed(o: TrackOrientation, reverse: bool) -> bool {
    o == TrackOrientation::Both || (o == TrackOrientation::Reverse) == reverse
}
fn positions_parts(positions: &[usize], reverse: bool) -> Vec<ReferencePart> {
    let mut parts: Vec<ReferencePart> = vec![];
    for &p in positions {
        if let Some(last) = parts.last_mut() {
            if !reverse && last.end == p {
                last.end = p + 1;
                continue;
            }
            if reverse && last.start == p + 1 {
                last.start = p;
                continue;
            }
        }
        parts.push(ReferencePart {
            start: p,
            end: p + 1,
        });
    }
    parts
}
impl PreparedTracks {
    fn visible_features(&self, out: &mut TrackReport) -> BTreeMap<usize, Vec<ReferencePart>> {
        let mut selected: BTreeMap<usize, Vec<ReferencePart>> = BTreeMap::new();
        let mut count = 0;
        let mut seen = BTreeSet::new();
        for part in window_parts(out.window, self.molecule.sequence.len()) {
            let first = self.prefix_end.partition_point(|&e| e <= part.start);
            for s in self.spans[first..]
                .iter()
                .take_while(|s| s.start < part.end)
            {
                if !spend(out, 1) {
                    return selected;
                }
                if s.end > part.start {
                    if !seen.insert((s.feature, s.start, s.end)) {
                        continue;
                    }
                    if count >= MAX_TRACK_PARTS {
                        warn(
                            out,
                            "Track reference-part limit reached; additional results are omitted.",
                        );
                        return selected;
                    }
                    let parts = selected.entry(s.feature).or_default();
                    let span = ReferencePart {
                        start: s.start,
                        end: s.end,
                    };
                    parts.push(span);
                    count += 1;
                }
            }
        }
        selected
    }
    fn annotations(&self, out: &mut TrackReport) {
        if !out.options.features && !out.options.primers {
            return;
        }
        for (index, visible_parts) in self.visible_features(out) {
            let f = &self.molecule.features[index];
            let primer = f.kind.eq_ignore_ascii_case("primer_bind");
            if (primer && !out.options.primers) || (!primer && !out.options.features) {
                continue;
            }
            let mut t = item(
                format!("feature:{index}"),
                if primer {
                    TrackKind::Primer
                } else {
                    TrackKind::Feature
                },
                self.feature_labels[index].clone(),
                visible_parts,
                self.feature_meta[index].0,
            );
            t.source_feature = Some(index);
            t.certain = !f.needs_review && t.reverse.is_some();
            if !emit(out, t) {
                break;
            }
        }
    }
    fn enzymes(&self, out: &mut TrackReport) -> Result<(), Error> {
        let m = &self.molecule;
        let n = m.sequence.len();
        let circular = m.topology == Topology::Circular;
        let source = m.sequence.as_str().as_bytes();
        for enzyme in crate::analysis::starter_enzymes()
            .into_iter()
            .filter(|e| out.options.enzyme_names.contains(&e.name))
            .collect::<Vec<_>>()
        {
            let pattern = Sequence::dna(&enzyme.recognition)?;
            let rc = pattern.reverse_complement()?;
            let len = pattern.len();
            if len > n {
                continue;
            }
            for (p, reverse) in [(pattern.as_str(), false), (rc.as_str(), true)]
                .into_iter()
                .take(if pattern == rc { 1 } else { 2 })
            {
                if pattern != rc && !allowed(out.options.orientation, reverse) {
                    continue;
                }
                // Include upstream recognition context and Type IIS sites whose cuts fall in the window.
                let halo = len
                    + enzyme
                        .top_cut
                        .unsigned_abs()
                        .max(enzyme.bottom_cut.unsigned_abs()) as usize;
                let mut starts = BTreeSet::new();
                for part in window_parts(out.window, n) {
                    for raw in part.start as i64 - halo as i64..part.end as i64 + halo as i64 {
                        if circular {
                            starts.insert(raw.rem_euclid(n as i64) as usize);
                        } else if raw >= 0 && raw as usize + len <= n {
                            starts.insert(raw as usize);
                        }
                    }
                }
                for start in starts {
                    if !spend(out, len) {
                        return Ok(());
                    }
                    let mut certain = true;
                    let mut matched = true;
                    for (j, b) in p.bytes().enumerate() {
                        let a = crate::analysis::base_mask(source[(start + j) % n]);
                        let b = crate::analysis::base_mask(b);
                        if a & b == 0 {
                            matched = false;
                            break;
                        }
                        certain &= a & !b == 0;
                    }
                    if !matched {
                        continue;
                    }
                    let parts = positions_parts(
                        &(0..len).map(|j| (start + j) % n).collect::<Vec<_>>(),
                        false,
                    );
                    let (top, bottom) = if reverse {
                        (
                            start as i64 + len as i64 - i64::from(enzyme.bottom_cut),
                            start as i64 + len as i64 - i64::from(enzyme.top_cut),
                        )
                    } else {
                        (
                            start as i64 + i64::from(enzyme.top_cut),
                            start as i64 + i64::from(enzyme.bottom_cut),
                        )
                    };
                    let cut = |p: i64| {
                        if circular {
                            Some(p.rem_euclid(n as i64) as usize)
                        } else {
                            (p >= 0 && p <= n as i64).then_some(p as usize)
                        }
                    };
                    let cuts = CutPositions {
                        top: cut(top),
                        bottom: cut(bottom),
                    };
                    let show = (out.options.enzyme_recognition && overlaps(&parts, out.window, n))
                        || (out.options.enzyme_cuts
                            && [cuts.top, cuts.bottom]
                                .into_iter()
                                .flatten()
                                .any(|p| visible_boundary(p, out.window, n)));
                    if !show {
                        continue;
                    }
                    let mut t = item(
                        format!("enzyme:{}:{start}:{reverse}", enzyme.name),
                        TrackKind::Enzyme,
                        enzyme.name.clone(),
                        parts,
                        Some(reverse),
                    );
                    t.enzyme = Some(enzyme.name.clone());
                    t.certain = certain;
                    if out.options.enzyme_cuts {
                        if cuts.top.is_none() || cuts.bottom.is_none() {
                            out.warnings.push(format!(
                                "{} at {} lacks sequence for both modeled cuts.",
                                enzyme.name, start
                            ));
                        }
                        t.cuts = Some(cuts);
                    }
                    if !emit(out, t) {
                        return Ok(());
                    }
                }
            }
        }
        Ok(())
    }
    fn frames(&self, out: &mut TrackReport) -> Result<(), Error> {
        let n = self.molecule.sequence.len();
        if n < 3 {
            return Ok(());
        }
        let circular = self.molecule.topology == Topology::Circular;
        if circular && n % 3 != 0 {
            out.warnings.push("Circular frame display is anchored to the reference origin; this source length does not preserve phase after a full turn.".into());
        }
        let bytes = self.molecule.sequence.as_str().as_bytes();
        for reverse in [false, true] {
            if !allowed(out.options.orientation, reverse) {
                continue;
            }
            for frame in out.options.frames.clone() {
                let mut starts = BTreeSet::new();
                for p in out.window.start..out.window.end {
                    if !spend(out, 1) {
                        return Ok(());
                    }
                    let x = if reverse { n - 1 - p % n } else { p % n };
                    let f = usize::from(frame);
                    if x >= f {
                        let start = f + (x - f) / 3 * 3;
                        if circular || start + 3 <= n {
                            starts.insert(start);
                        }
                    }
                    // A final circular codon can cover reference-origin bases even when frame == 0.
                    if circular {
                        let last = f + (n - 1 - f) / 3 * 3;
                        if last + 3 > n && x + n < last + 3 {
                            starts.insert(last);
                        }
                    }
                }
                let mut codons = String::new();
                let mut mappings = vec![];
                for start in starts {
                    if out.work.emitted_residues >= out.options.max_residues {
                        warn(
                            out,
                            "Track residue limit reached; additional residues are omitted.",
                        );
                        break;
                    }
                    if !spend(out, 3) {
                        break;
                    }
                    let positions: Vec<_> = (start..start + 3)
                        .map(|x| if reverse { n - 1 - x % n } else { x % n })
                        .collect();
                    for &p in &positions {
                        codons.push(if reverse {
                            collagene_model::complement(bytes[p], self.molecule.sequence.alphabet())
                        } else {
                            bytes[p]
                        } as char);
                    }
                    mappings.push(ResidueMapping {
                        residue_index: (start - usize::from(frame)) / 3,
                        amino_acid: String::new(),
                        codon_parts: positions_parts(&positions, reverse),
                    });
                    out.work.emitted_residues += 1;
                }
                if mappings.is_empty() {
                    continue;
                }
                let decoded = translation::translate(
                    &Sequence::new(self.molecule.sequence.alphabet(), &codons)?,
                    &TranslationOptions {
                        genetic_code: out.options.genetic_code,
                        frame: 0,
                        reverse: false,
                        initiation: InitiationPolicy::Ordinary,
                    },
                )?;
                out.warnings.extend(decoded.warnings.iter().map(|w| {
                    format!(
                        "Frame {}{}: {w}",
                        if reverse { "-" } else { "+" },
                        frame + 1
                    )
                }));
                let certain = decoded.unresolved_codons == 0 && decoded.possible_stop_codons == 0;
                let protein = decoded.protein;
                for (mapping, aa) in mappings.iter_mut().zip(protein.chars()) {
                    mapping.amino_acid = aa.to_string();
                }
                let mut t = item(
                    format!("frame:{reverse}:{frame}"),
                    TrackKind::Translation,
                    format!("{}{}", if reverse { "−" } else { "+" }, frame + 1),
                    window_parts(out.window, n),
                    Some(reverse),
                );
                t.frame = Some(frame);
                t.certain = certain;
                t.protein = Some(protein);
                t.residues = mappings;
                if !emit(out, t) {
                    return Ok(());
                }
            }
        }
        Ok(())
    }
}

/// Exact ordered feature coordinates in biological 5′→3′ order.
fn coding_positions(l: &Location) -> Option<(Vec<usize>, bool)> {
    let reverse = direction(l)?;
    fn walk(l: &Location, reverse: bool, positions: &mut Vec<usize>) {
        match l {
            Location::Range { start, end, .. } => {
                if reverse {
                    positions.extend((*start..*end).rev());
                } else {
                    positions.extend(*start..*end);
                }
            }
            Location::Complement(inner) => walk(inner, !reverse, positions),
            Location::Join(parts) => {
                if reverse {
                    for p in parts.iter().rev() {
                        walk(p, reverse, positions);
                    }
                } else {
                    for p in parts {
                        walk(p, reverse, positions);
                    }
                }
            }
            _ => {}
        }
    }
    let mut positions = vec![];
    walk(l, false, &mut positions);
    Some((positions, reverse))
}
impl PreparedTracks {
    fn cds(&self, out: &mut TrackReport) -> Result<(), Error> {
        let n = self.molecule.sequence.len();
        for (index, visible_parts) in self.visible_features(out) {
            let f = &self.molecule.features[index];
            if !f.kind.eq_ignore_ascii_case("CDS") {
                continue;
            }
            let Some(reverse) = self.feature_meta[index].0 else {
                warn(
                    out,
                    "CDS translation omitted: fuzzy, unordered, or mixed-direction location has no supported coding frame.",
                );
                continue;
            };
            if !allowed(out.options.orientation, reverse) {
                continue;
            }
            let size = self.feature_meta[index].1;
            // Extraction, coordinate materialization, and translation are all bounded before allocation.
            if !spend(
                out,
                size.saturating_mul(4)
                    .saturating_add(self.qualifier_work[index])
                    .saturating_add(1),
            ) {
                continue;
            }
            let translated =
                match translation::translate_cds_validated(&self.molecule, f, &self.revision) {
                    Ok(t) => t,
                    Err(error) => {
                        warn(
                            out,
                            &format!("CDS feature {index} translation omitted: {error}"),
                        );
                        continue;
                    }
                };
            let (positions, _) = coding_positions(&f.location).expect("exact direction checked");
            let mut residues = vec![];
            let mut protein = String::new();
            for (residue_index, aa) in translated.protein.chars().enumerate() {
                let start = usize::from(translated.codon_start - 1) + 3 * residue_index;
                let Some(codon) = positions.get(start..start + 3) else {
                    break;
                };
                let parts = positions_parts(codon, reverse);
                if !overlaps(&parts, out.window, n) {
                    continue;
                }
                if out.work.emitted_residues >= out.options.max_residues {
                    warn(
                        out,
                        "Track residue limit reached; additional residues are omitted.",
                    );
                    break;
                }
                residues.push(ResidueMapping {
                    residue_index,
                    amino_acid: aa.to_string(),
                    codon_parts: parts,
                });
                protein.push(aa);
                out.work.emitted_residues += 1;
            }
            out.warnings.extend(
                translated
                    .warnings
                    .iter()
                    .map(|w| format!("CDS feature {index}: {w}")),
            );
            let mut t = item(
                format!("cds:{index}"),
                TrackKind::Cds,
                self.feature_labels[index].clone(),
                visible_parts,
                Some(reverse),
            );
            t.source_feature = Some(index);
            t.frame = Some(translated.codon_start - 1);
            t.protein = Some(protein);
            t.residues = residues;
            t.certain = translated.complete_under_model;
            if !emit(out, t) {
                break;
            }
        }
        Ok(())
    }
    fn orfs(&self, out: &mut TrackReport) -> Result<(), Error> {
        let n = self.molecule.sequence.len();
        // Whole-source ORFs retain real starts/stops outside the window. Never call a cropped
        // window a complete ORF. Precharge conservative scan + worst-case output materialization.
        let remaining = out.work.limit_units - out.work.used_units;
        let scan = n.saturating_mul(8).saturating_add(16);
        if remaining < scan.saturating_add(n) {
            warn(
                out,
                "ORF track omitted: whole-source scan exceeds remaining work budget.",
            );
            return Ok(());
        }
        let max_orfs = ((remaining - scan) / n).min(out.options.max_items).max(1);
        if !spend(out, scan + max_orfs * n) {
            return Ok(());
        }
        let options = translation::OrfOptions {
            genetic_code: out.options.genetic_code,
            min_amino_acids: out.options.min_orf_amino_acids,
            max_orfs,
            both_strands: true,
            ..Default::default()
        };
        let report = translation::find_orfs_validated(&self.molecule, &options, &self.revision)?;
        if report.truncated {
            warn(
                out,
                "ORF prediction cap reached before viewport filtering; additional candidates may overlap the window.",
            );
        }
        out.warnings.extend(report.warnings);
        for orf in report.orfs {
            if !allowed(out.options.orientation, orf.reverse)
                || !out.options.frames.contains(&orf.frame)
            {
                continue;
            }
            let parts = location_parts(&orf.location);
            if !overlaps(&parts, out.window, n) {
                continue;
            }
            let mut residues = vec![];
            let mut protein = String::new();
            for (residue_index, aa) in orf.protein.chars().enumerate() {
                let positions: Vec<_> = (0..3)
                    .map(|j| {
                        let x = (orf.oriented_start + 3 * residue_index + j) % n;
                        if orf.reverse { n - 1 - x } else { x }
                    })
                    .collect();
                let codon_parts = positions_parts(&positions, orf.reverse);
                if !overlaps(&codon_parts, out.window, n) {
                    continue;
                }
                if out.work.emitted_residues >= out.options.max_residues {
                    warn(
                        out,
                        "Track residue limit reached; additional residues are omitted.",
                    );
                    break;
                }
                residues.push(ResidueMapping {
                    residue_index,
                    amino_acid: aa.to_string(),
                    codon_parts,
                });
                protein.push(aa);
                out.work.emitted_residues += 1;
            }
            let mut t = item(
                format!(
                    "orf:{}:{}:{}",
                    orf.reverse, orf.oriented_start, orf.nucleotide_length
                ),
                TrackKind::Orf,
                "Predicted ORF".into(),
                parts,
                Some(orf.reverse),
            );
            t.frame = Some(orf.frame);
            t.certain = orf.has_start_codon
                && orf.has_stop_codon
                && !orf.partial_5prime
                && !orf.partial_3prime
                && orf.ambiguous_codons == 0;
            t.protein = Some(protein);
            t.residues = residues;
            if !emit(out, t) {
                break;
            }
        }
        Ok(())
    }
}

impl TrackReport {
    /// Structural trust boundary for decoded worker/MCP reports. It cannot verify biology
    /// without the source; consumers must additionally match source_revision to their document.
    pub fn validate(&self) -> Result<(), Error> {
        let bad = || Error::Format("invalid bounded browser-track report".into());
        let n = self.work.source_bases;
        if self.engine != "collagene-core"
            || self.version != TRACK_ENGINE_VERSION
            || self.source_revision.len() != 64
            || !self.source_revision.bytes().all(|b| b.is_ascii_hexdigit())
            || n == 0
            || n > collagene_model::MAX_SEQUENCE_LEN
            || self.window.start >= n
            || self.window.end <= self.window.start
            || self.window.end - self.window.start > n
            || self.window.end > 2 * n
            || self.window.end - self.window.start > self.options.max_window_bases
            || self.options.max_window_bases == 0
            || self.options.max_window_bases > MAX_TRACK_WINDOW
            || self.options.max_items == 0
            || self.options.max_items > MAX_TRACK_ITEMS
            || self.options.max_residues == 0
            || self.options.max_residues > MAX_TRACK_RESIDUES
            || self.options.max_work_units == 0
            || self.options.max_work_units > MAX_TRACK_WORK
            || self.work.limit_units != self.options.max_work_units
            || self.work.used_units > self.work.limit_units
            || self.work.indexed_spans > crate::tracks::MAX_INDEX_SPANS
            || self.items.len() > self.options.max_items
            || self.work.emitted_residues > self.options.max_residues
            || self.options.frames.len() > 3
            || self.options.frames.iter().any(|f| *f > 2)
            || self.options.frames.iter().collect::<BTreeSet<_>>().len()
                != self.options.frames.len()
            || self.options.enzyme_names.len() > 12
            || self
                .options
                .enzyme_names
                .iter()
                .collect::<BTreeSet<_>>()
                .len()
                != self.options.enzyme_names.len()
            || !translation::genetic_codes()
                .iter()
                .any(|c| c.id == self.options.genetic_code)
            || self.options.min_orf_amino_acids > collagene_model::MAX_SEQUENCE_LEN / 3
        {
            return Err(bad());
        }
        let catalog = crate::analysis::starter_enzymes();
        if self
            .options
            .enzyme_names
            .iter()
            .any(|s| !catalog.iter().any(|e| &e.name == s))
        {
            return Err(bad());
        }
        if self.items.iter().map(|t| t.parts.len()).sum::<usize>() > MAX_TRACK_PARTS {
            return Err(bad());
        }
        if self.warnings.len() > MAX_TRACK_WARNINGS || self.warnings.iter().any(|w| w.len() > 4096)
        {
            return Err(bad());
        }
        let mut ids = BTreeSet::new();
        let mut residues = 0usize;
        let valid_parts = |parts: &[ReferencePart]| {
            !parts.is_empty() && parts.iter().all(|p| p.start < p.end && p.end <= n)
        };
        for item in &self.items {
            if item.id.is_empty()
                || item.id.len() > 128
                || item.label.len() > 2048
                || !ids.insert(&item.id)
                || !valid_parts(&item.parts)
                || item.frame.is_some_and(|f| f > 2)
                || item
                    .cuts
                    .as_ref()
                    .is_some_and(|c| [c.top, c.bottom].into_iter().flatten().any(|p| p > n))
                || item
                    .source_feature
                    .is_some_and(|f| f >= crate::tracks::MAX_INDEX_SPANS)
            {
                return Err(bad());
            }
            if item.kind == TrackKind::Enzyme {
                if item
                    .enzyme
                    .as_ref()
                    .is_none_or(|name| !self.options.enzyme_names.contains(name))
                    || item.cuts.is_some() != self.options.enzyme_cuts
                    || (!self.options.enzyme_recognition && !self.options.enzyme_cuts)
                {
                    return Err(bad());
                }
            } else if item.enzyme.is_some() || item.cuts.is_some() {
                return Err(bad());
            }
            if (item.kind == TrackKind::Feature && !self.options.features)
                || (item.kind == TrackKind::Primer && !self.options.primers)
                || (item.kind == TrackKind::Translation && !self.options.translations)
                || (item.kind == TrackKind::Cds && self.options.known_cds == KnownCdsPolicy::Off)
                || (item.kind == TrackKind::Orf && !self.options.orfs)
            {
                return Err(bad());
            }
            let protein_kind = matches!(
                item.kind,
                TrackKind::Translation | TrackKind::Cds | TrackKind::Orf
            );
            if protein_kind != (item.protein.is_some())
                || (protein_kind && (item.frame.is_none() || item.reverse.is_none()))
                || (!protein_kind && !item.residues.is_empty())
            {
                return Err(bad());
            }
            if let Some(protein) = &item.protein {
                if !protein.is_ascii()
                    || protein.len() != item.residues.len()
                    || !protein
                        .bytes()
                        .all(|b| b"ACDEFGHIKLMNPQRSTVWYBXZJUO*".contains(&b))
                {
                    return Err(bad());
                }
                let mut seen = BTreeSet::new();
                for (mapping, aa) in item.residues.iter().zip(protein.bytes()) {
                    if mapping.amino_acid.as_bytes() != [aa]
                        || !seen.insert(mapping.residue_index)
                        || !valid_parts(&mapping.codon_parts)
                        || mapping
                            .codon_parts
                            .iter()
                            .map(|p| p.end - p.start)
                            .sum::<usize>()
                            != 3
                        || !overlaps(&mapping.codon_parts, self.window, n)
                    {
                        return Err(bad());
                    }
                }
            }
            residues = residues.checked_add(item.residues.len()).ok_or_else(bad)?;
        }
        if residues != self.work.emitted_residues {
            return Err(bad());
        }
        Ok(())
    }
}
