//! Exact oriented overlap checks and explicitly advisory published method guidance.
//! No reagent, kinetics, proprietary chemistry or reaction-mixture uniqueness model.
use crate::{Alphabet, Error, Location, Molecule, Topology, model::complement};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

pub const ENGINE_VERSION: &str = "overlap-profiles-1";
pub const MODEL: &str = "exact-oriented-end-homology/vendor-guidance-2026-11-24";
const MAX_FRAGMENTS: usize = 1000;
const MAX_SOURCE_BASES: usize = 16 * 1024 * 1024;
const MAX_SOURCE_BYTES: usize = 32 * 1024 * 1024;
const MAX_OVERLAP: usize = 10_000;
const MAX_HOMOLOGY_BASES: usize = 1_000_000;
const MAX_COMPARISONS: usize = 50_000_000;
const MAX_OCCURRENCES: usize = 10_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum OverlapMethod {
    Gibson,
    HiFi,
    InFusion,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProfileOptions {
    pub method: OverlapMethod,
    pub scan_repeated_overlaps: bool,
    pub max_comparisons: usize,
    pub max_occurrences: usize,
}
impl Default for ProfileOptions {
    fn default() -> Self {
        Self {
            method: OverlapMethod::HiFi,
            scan_repeated_overlaps: true,
            max_comparisons: 5_000_000,
            max_occurrences: 1000,
        }
    }
}
impl ProfileOptions {
    pub fn validate(&self) -> Result<(), Error> {
        if self.max_comparisons == 0
            || self.max_comparisons > MAX_COMPARISONS
            || self.max_occurrences == 0
            || self.max_occurrences > MAX_OCCURRENCES
        {
            return Err(Error::Unsupported(
                "overlap profile limits must be comparisons1..50,000,000 and occurrences1..10,000"
                    .into(),
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ReferenceSpan {
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Occurrence {
    pub fragment: usize,
    pub start: usize,
    pub end: usize,
    /// Strand in the ORIGINAL input reference. Palindromic words are canonically forward.
    pub reverse: bool,
    pub assigned_junction: bool,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct JunctionReport {
    pub index: usize,
    pub upstream: usize,
    pub downstream: usize,
    pub overlap_bases: usize,
    pub upstream_span: Option<ReferenceSpan>,
    pub downstream_span: Option<ReferenceSpan>,
    /// Ordered 5′→3′ strings in the declared assembly orientation, not raw reference slices.
    pub homology: String,
    pub downstream_homology: String,
    pub exact_match: bool,
    pub gc_fraction: Option<f64>,
    pub palindromic: bool,
    pub recommended_min: Option<usize>,
    pub recommended_max: Option<usize>,
    /// Only a LENGTH guideline; it is not full kit compatibility or empirical success.
    pub within_recommended: Option<bool>,
    pub occurrences: Vec<Occurrence>,
    pub repeat_search_complete: bool,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ProfileReport {
    pub source_revisions: Vec<String>,
    pub fragment_lengths: Vec<usize>,
    pub orientations: Vec<bool>,
    pub overlaps: Vec<usize>,
    pub circular: bool,
    pub options: ProfileOptions,
    pub engine_version: String,
    pub model: String,
    pub geometry_valid: bool,
    pub recommendation_compliant: Option<bool>,
    /// Requested repeat-site scan completed. Never means exhaustive alternative assemblies.
    pub complete: bool,
    pub comparisons: usize,
    pub source_bases: usize,
    pub junctions: Vec<JunctionReport>,
    pub findings: Vec<String>,
    pub warnings: Vec<String>,
    pub exhaustive_alternative_search: bool,
}
fn bad(s: impl Into<String>) -> Error {
    Error::Format(s.into())
}
fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| complement(b, Alphabet::Dna) as char)
        .collect()
}
fn base(m: &Molecule, reverse: bool, i: usize) -> u8 {
    if reverse {
        complement(
            m.sequence.as_str().as_bytes()[m.sequence.len() - 1 - i],
            Alphabet::Dna,
        )
    } else {
        m.sequence.as_str().as_bytes()[i]
    }
}
fn oriented_text(m: &Molecule, reverse: bool, start: usize, end: usize) -> String {
    (start..end).map(|i| base(m, reverse, i) as char).collect()
}
fn span(n: usize, reverse: bool, start: usize, end: usize) -> ReferenceSpan {
    if reverse {
        ReferenceSpan {
            start: n - end,
            end: n - start,
            reverse,
        }
    } else {
        ReferenceSpan {
            start,
            end,
            reverse,
        }
    }
}
fn guideline(method: OverlapMethod, fragments: usize) -> (Option<usize>, Option<usize>) {
    match (method, fragments) {
        (OverlapMethod::HiFi, 2..=3) => (Some(15), Some(20)),
        (OverlapMethod::HiFi, 4..=6) => (Some(20), Some(30)),
        (OverlapMethod::Gibson, 2..=3) => (Some(15), Some(25)),
        (OverlapMethod::Gibson, 4..=6) => (Some(20), Some(80)),
        (OverlapMethod::InFusion, 2) => (Some(15), None),
        (OverlapMethod::InFusion, 3..=1000) => (Some(20), None),
        _ => (None, None),
    }
}
fn within(length: usize, min: Option<usize>, max: Option<usize>) -> Option<bool> {
    min.map(|min| length >= min && max.is_none_or(|max| length <= max))
}
fn gc(s: &str) -> Option<f64> {
    (!s.is_empty())
        .then(|| s.bytes().filter(|b| matches!(b, b'G' | b'C')).count() as f64 / s.len() as f64)
}
fn is_hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|c| c.is_ascii_digit() || matches!(c, b'a'..=b'f'))
}
fn source_preflight(molecules: &[Molecule]) -> Result<usize, Error> {
    fn reserve(bytes: &mut usize, items: &mut usize, size: usize) -> Result<(), Error> {
        *bytes = bytes.saturating_add(size).saturating_add(64);
        *items = items.saturating_add(1);
        if *bytes > MAX_SOURCE_BYTES || *items > 500_000 {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
    fn location(
        l: &Location,
        depth: usize,
        bytes: &mut usize,
        items: &mut usize,
    ) -> Result<(), Error> {
        if depth > 64 {
            return Err(Error::TooLarge);
        }
        reserve(bytes, items, 96)?;
        match l {
            Location::Complement(l) => location(l, depth + 1, bytes, items)?,
            Location::Join(parts) | Location::Order(parts) => {
                for p in parts {
                    location(p, depth + 1, bytes, items)?
                }
            }
            Location::Unresolved(s) => reserve(bytes, items, s.len())?,
            Location::Range { .. } => {}
        }
        Ok(())
    }
    let (mut bytes, mut items, mut bases) = (0usize, 0usize, 0usize);
    for m in molecules {
        bases = bases.saturating_add(m.sequence.len());
        if bases > MAX_SOURCE_BASES {
            return Err(Error::TooLarge);
        }
        reserve(
            &mut bytes,
            &mut items,
            m.name.len().saturating_add(m.sequence.len()),
        )?;
        for (k, v) in &m.metadata {
            reserve(&mut bytes, &mut items, k.len().saturating_add(v.len()))?;
        }
        if let Some(s) = &m.original_genbank {
            reserve(&mut bytes, &mut items, s.len())?;
        }
        for f in &m.features {
            reserve(&mut bytes, &mut items, f.kind.len())?;
            location(&f.location, 0, &mut bytes, &mut items)?;
            for (k, v) in &f.qualifiers {
                reserve(
                    &mut bytes,
                    &mut items,
                    k.len().saturating_add(v.as_ref().map_or(0, String::len)),
                )?;
            }
        }
    }
    Ok(bases)
}
/// All-forward convenience. The oriented entry point retains ORIGINAL revisions and coordinates.
pub fn check_fragments(
    molecules: &[Molecule],
    overlaps: &[usize],
    circular: bool,
    options: &ProfileOptions,
) -> Result<ProfileReport, Error> {
    if molecules.len() > MAX_FRAGMENTS {
        return Err(Error::TooLarge);
    }
    check_oriented_fragments(
        molecules,
        &vec![false; molecules.len()],
        overlaps,
        circular,
        options,
    )
}
pub fn check_oriented_fragments(
    molecules: &[Molecule],
    reverse: &[bool],
    overlaps: &[usize],
    circular: bool,
    options: &ProfileOptions,
) -> Result<ProfileReport, Error> {
    options.validate()?;
    if molecules.is_empty()
        || molecules.len() > MAX_FRAGMENTS
        || reverse.len() != molecules.len()
        || overlaps.len() != molecules.len() - 1 + usize::from(circular)
    {
        return Err(bad(
            "profile requires1..1000 ordered sources, one orientation per source and one overlap per junction",
        ));
    }
    let total_homology = overlaps
        .iter()
        .try_fold(0usize, |sum, n| sum.checked_add(*n).ok_or(Error::TooLarge))?;
    if overlaps.iter().any(|n| *n > MAX_OVERLAP) || total_homology > MAX_HOMOLOGY_BASES {
        return Err(Error::TooLarge);
    }
    if total_homology > options.max_comparisons {
        return Err(Error::Unsupported(
            "comparison budget is too small to establish the specified exact junction geometry"
                .into(),
        ));
    }
    let source_bases = source_preflight(molecules)?;
    for m in molecules {
        m.validate()?;
        if m.topology != Topology::Linear
            || m.sequence.alphabet() != Alphabet::Dna
            || m.sequence.as_str().bytes().any(|b| !b"ACGT".contains(&b))
        {
            return Err(Error::Unsupported("profile inputs must be linear unambiguous DNA; no ambiguity or topology is repaired by assumption".into()));
        }
    }
    let n = molecules.len();
    let (min, max) = guideline(options.method, n);
    let mut report=ProfileReport{source_revisions:molecules.iter().map(Molecule::revision).collect(),fragment_lengths:molecules.iter().map(|m|m.sequence.len()).collect(),orientations:reverse.to_vec(),overlaps:overlaps.to_vec(),circular,options:options.clone(),engine_version:ENGINE_VERSION.into(),model:MODEL.into(),geometry_valid:true,recommendation_compliant:min.map(|_|true),complete:true,comparisons:0,source_bases,junctions:vec![],findings:vec![],warnings:vec!["This report checks an explicit ordered exact-homology candidate. Published length guidance is advisory and does not establish reagent compatibility, kinetics, yield, purity, transformation outcome or empirical success.".into(),"Sequence profiles do not establish physical end chemistry, accessibility, methylation or covalent completion. The originating physical graph and assembly reconstruction remain separate checks.".into(),"A complete motif scan is not an exhaustive alternative assembly search: copy number, order/orientation combinations, shorter microhomologies, mismatches and reaction products are not enumerated.".into()],exhaustive_alternative_search:false};
    if min.is_none() {
        report.warnings.push("The selected fragment count is outside the versioned published length-guidance table; no recommendation status is inferred.".into());
    }
    if options.method == OverlapMethod::InFusion {
        report.warnings.push("In-Fusion length screening uses minimum15nt for two fragments or20nt for multiple inserts; no upper cutoff or guarantee for longer homologies is inferred from that guidance.".into());
    }
    for (i, &k) in overlaps.iter().enumerate() {
        let a = i;
        let b = (i + 1) % n;
        let an = molecules[a].sequence.len();
        let bn = molecules[b].sequence.len();
        let bounded = k > 0 && k < an && k < bn;
        let (upstream_span, downstream_span, homology, downstream) = if bounded {
            (
                Some(span(an, reverse[a], an - k, an)),
                Some(span(bn, reverse[b], 0, k)),
                oriented_text(&molecules[a], reverse[a], an - k, an),
                oriented_text(&molecules[b], reverse[b], 0, k),
            )
        } else {
            (None, None, String::new(), String::new())
        };
        let exact = bounded
            && homology
                .bytes()
                .zip(downstream.bytes())
                .fold(true, |equal, (x, y)| {
                    report.comparisons += 1;
                    equal && x == y
                });
        if !exact {
            report.geometry_valid = false;
            report.findings.push(format!("Junction {} has zero/exhaustive/out-of-range homology or nonidentical oriented terminal bases; exact assembly is not established.",i+1));
        }
        let recommended = within(k, min, max);
        if recommended == Some(false) {
            report.recommendation_compliant = Some(false);
            report.findings.push(format!("Junction {} overlap length {k} is outside the versioned {:?} length guidance; this is an advisory finding, not a sequence incompatibility.",i+1,options.method));
        }
        let palindromic = !homology.is_empty() && rc(&homology) == homology;
        report.junctions.push(JunctionReport {
            index: i,
            upstream: a,
            downstream: b,
            overlap_bases: k,
            upstream_span,
            downstream_span,
            gc_fraction: gc(&homology),
            palindromic,
            homology,
            downstream_homology: downstream,
            exact_match: exact,
            recommended_min: min,
            recommended_max: max,
            within_recommended: recommended,
            occurrences: vec![],
            repeat_search_complete: !options.scan_repeated_overlaps,
        });
    }
    // Same explicit unsupported state as the existing exact ordered reconstruction.
    if report.geometry_valid {
        for (i, m) in molecules.iter().enumerate() {
            let incoming = if i > 0 {
                overlaps[i - 1]
            } else if circular {
                overlaps[n - 1]
            } else {
                0
            };
            let outgoing = if i < n - 1 {
                overlaps[i]
            } else if circular {
                overlaps[n - 1]
            } else {
                0
            };
            if incoming.saturating_add(outgoing) >= m.sequence.len() {
                return Err(Error::Unsupported(format!(
                    "fragment {} has intersecting/exhaustive overlap arms; a more detailed assembly model is required",
                    i + 1
                )));
            }
        }
        let length = source_bases
            .checked_sub(total_homology)
            .ok_or(Error::TooLarge)?;
        if length == 0 || length > crate::model::MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
    }
    if options.scan_repeated_overlaps {
        let mut stopped = false;
        let mut occurrence_count = 0;
        for j in &mut report.junctions {
            if stopped {
                j.repeat_search_complete = false;
                continue;
            }
            if j.homology.is_empty() {
                j.repeat_search_complete = true;
                continue;
            }
            let pattern = j.homology.as_bytes();
            let reversed = rc(&j.homology);
            let patterns = if j.palindromic {
                vec![(pattern, false)]
            } else {
                vec![(pattern, false), (reversed.as_bytes(), true)]
            };
            'scan: for (fragment, m) in molecules.iter().enumerate() {
                if m.sequence.len() < pattern.len() {
                    continue;
                }
                for (p, motif_reverse) in &patterns {
                    for start in 0..=m.sequence.len() - p.len() {
                        let mut matched = true;
                        for (offset, &letter) in p.iter().enumerate() {
                            if report.comparisons == options.max_comparisons {
                                stopped = true;
                                break 'scan;
                            }
                            report.comparisons += 1;
                            if base(m, reverse[fragment], start + offset) != letter {
                                matched = false;
                                break;
                            }
                        }
                        if matched {
                            if occurrence_count == options.max_occurrences {
                                stopped = true;
                                break 'scan;
                            }
                            let s =
                                span(m.sequence.len(), reverse[fragment], start, start + p.len());
                            let strand = if j.palindromic {
                                false
                            } else {
                                reverse[fragment] ^ *motif_reverse
                            };
                            let assigned = [
                                (j.upstream, &j.upstream_span),
                                (j.downstream, &j.downstream_span),
                            ]
                            .iter()
                            .any(|(f, span)| {
                                *f == fragment
                                    && span.as_ref().is_some_and(|s2| {
                                        s.start == s2.start
                                            && s.end == s2.end
                                            && (j.palindromic || strand == s2.reverse)
                                    })
                            });
                            j.occurrences.push(Occurrence {
                                fragment,
                                start: s.start,
                                end: s.end,
                                reverse: strand,
                                assigned_junction: assigned,
                            });
                            occurrence_count += 1;
                        }
                    }
                }
            }
            j.repeat_search_complete = !stopped;
        }
        report.complete = !stopped;
        if stopped {
            report.warnings.push("Repeated-overlap search stopped at its comparison or occurrence limit; omitted matches do not establish uniqueness or absence.".into());
        }
        for j in &report.junctions {
            let extras = j
                .occurrences
                .iter()
                .filter(|o| !o.assigned_junction)
                .count();
            if extras > 0 {
                report.findings.push(format!("Junction {} homology has {extras} additional concrete occurrence(s) among the supplied fragments/strands. These expose possible competing sequence contexts, not enumerated alternative products.",j.index+1));
            }
        }
    } else {
        for j in &mut report.junctions {
            j.repeat_search_complete = false;
        }
        report.warnings.push(
            "Repeated-overlap search was disabled; no uniqueness inference is available.".into(),
        );
    }
    report.validate()?;
    Ok(report)
}
impl ProfileReport {
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        let n = self.source_revisions.len();
        if n == 0
            || n > MAX_FRAGMENTS
            || self.fragment_lengths.len() != n
            || self.orientations.len() != n
            || self.source_revisions.iter().any(|s| !is_hex(s))
            || self
                .fragment_lengths
                .iter()
                .any(|n| *n == 0 || *n > crate::model::MAX_SEQUENCE_LEN)
            || self
                .fragment_lengths
                .iter()
                .try_fold(0usize, |s, n| s.checked_add(*n))
                != Some(self.source_bases)
            || self.source_bases > MAX_SOURCE_BASES
            || self.overlaps.len() != n - 1 + usize::from(self.circular)
            || self.junctions.len() != self.overlaps.len()
            || self.overlaps.iter().any(|k| *k > MAX_OVERLAP)
            || self
                .overlaps
                .iter()
                .try_fold(0usize, |s, k| s.checked_add(*k))
                .is_none_or(|k| k > MAX_HOMOLOGY_BASES)
            || self.comparisons > self.options.max_comparisons
            || self.engine_version != ENGINE_VERSION
            || self.model != MODEL
            || self.exhaustive_alternative_search
            || self.findings.len() > 4000
            || self.warnings.len() > 16
            || self
                .findings
                .iter()
                .chain(&self.warnings)
                .any(|s| s.len() > 2048)
        {
            return Err(bad("invalid overlap-profile report identity or bounds"));
        }
        let (min, max) = guideline(self.options.method, n);
        let mut geometry = true;
        let mut recommended = min.map(|_| true);
        let mut count = 0usize;
        let mut scan_complete = true;
        let mut charged = 0usize;
        for (i, j) in self.junctions.iter().enumerate() {
            let a = i;
            let b = (i + 1) % n;
            let k = self.overlaps[i];
            let an = self.fragment_lengths[a];
            let bn = self.fragment_lengths[b];
            let bounded = k > 0 && k < an && k < bn;
            let expected_up = bounded.then(|| span(an, self.orientations[a], an - k, an));
            let expected_down = bounded.then(|| span(bn, self.orientations[b], 0, k));
            let expected_length = if bounded { k } else { 0 };
            let exact = bounded && j.homology == j.downstream_homology;
            if j.index != i
                || j.upstream != a
                || j.downstream != b
                || j.overlap_bases != k
                || j.upstream_span != expected_up
                || j.downstream_span != expected_down
                || j.homology.len() != expected_length
                || j.downstream_homology.len() != expected_length
                || j.homology
                    .bytes()
                    .chain(j.downstream_homology.bytes())
                    .any(|b| !b"ACGT".contains(&b))
                || j.exact_match != exact
                || j.gc_fraction != gc(&j.homology)
                || j.palindromic != (!j.homology.is_empty() && rc(&j.homology) == j.homology)
                || j.recommended_min != min
                || j.recommended_max != max
                || j.within_recommended != within(k, min, max)
            {
                return Err(bad("invalid overlap-profile junction evidence"));
            }
            geometry &= exact;
            if j.within_recommended == Some(false) {
                recommended = Some(false);
            }
            charged += expected_length;
            if !self.options.scan_repeated_overlaps {
                if j.repeat_search_complete || !j.occurrences.is_empty() {
                    return Err(bad("disabled repeat scan has claimed evidence"));
                }
            } else {
                scan_complete &= j.repeat_search_complete;
            }
            if j.occurrences.len() > self.options.max_occurrences.saturating_sub(count) {
                return Err(bad("repeated-overlap occurrence budget exceeded"));
            }
            let mut seen = BTreeSet::new();
            for o in &j.occurrences {
                count += 1;
                if o.fragment >= n
                    || o.start >= o.end
                    || o.end > self.fragment_lengths[o.fragment]
                    || o.end - o.start != k
                    || (!seen.insert((o.fragment, o.start, o.end, o.reverse)))
                    || (j.palindromic && o.reverse)
                {
                    return Err(bad("invalid repeated-overlap occurrence"));
                }
                let assigned =
                    [(a, &j.upstream_span), (b, &j.downstream_span)]
                        .iter()
                        .any(|(f, s)| {
                            *f == o.fragment
                                && s.as_ref().is_some_and(|s| {
                                    s.start == o.start
                                        && s.end == o.end
                                        && (j.palindromic || s.reverse == o.reverse)
                                })
                        });
                if o.assigned_junction != assigned {
                    return Err(bad("invalid assigned junction occurrence"));
                }
            }
        }
        if self.geometry_valid != geometry
            || self.recommendation_compliant != recommended
            || self.comparisons < charged
            || count > self.options.max_occurrences
            || self.complete != (scan_complete || !self.options.scan_repeated_overlaps)
        {
            return Err(bad(
                "invalid overlap-profile completeness or derived status",
            ));
        }
        if geometry {
            for (i, &length) in self.fragment_lengths.iter().enumerate() {
                let incoming = if i > 0 {
                    self.overlaps[i - 1]
                } else if self.circular {
                    self.overlaps[n - 1]
                } else {
                    0
                };
                let outgoing = if i < n - 1 {
                    self.overlaps[i]
                } else if self.circular {
                    self.overlaps[n - 1]
                } else {
                    0
                };
                if incoming.saturating_add(outgoing) >= length {
                    return Err(bad("unmodeled intersecting overlap arms"));
                }
            }
            if self.source_bases - self.overlaps.iter().sum::<usize>()
                > crate::model::MAX_SEQUENCE_LEN
            {
                return Err(Error::TooLarge);
            }
        }
        Ok(())
    }
}
