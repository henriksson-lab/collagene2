//! Bounded recognition inspection and lazy, assumed-complete duplex products.
//! Recognition evidence is separate from substrate geometry and experimental cleavage.
use crate::{
    Alphabet, Error, Location, Molecule, Sequence, Topology,
    analysis::{self, Enzyme},
    chemistry,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const ENGINE_VERSION: &str = "restriction-inspection-1";
const INSPECTION_ASSUMPTIONS: &str = "Recognition and reference top-strand lengths do not establish experimental cleavage. Products assume complete cleavage of intact, unmodified-in-effect duplex DNA; methylation, reagent conditions and reaction success are not assessed.";
const SCAN_ASSUMPTIONS: &str = "Recognition evidence does not establish experimental cleavage. This scan computes no physical fragments or reference top-strand product lengths; methylation, reagent conditions and reaction success are not assessed.";
const MAX_SOURCE_BYTES: usize = 32 * 1024 * 1024;
const MAX_REPORT_BYTES: usize = 16 * 1024 * 1024;
macro_rules! dto {($($x:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
macro_rules! enumeration {($($x:tt)*)=>{#[derive(Debug,Clone,Copy,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]$($x)*};}
dto! {#[derive(PartialEq,Eq)] pub struct InspectionOptions {pub max_comparisons:u64,pub max_sites:usize,pub max_fragments:usize,pub max_report_bytes:usize}}
impl Default for InspectionOptions {
    fn default() -> Self {
        Self {
            max_comparisons: 50_000_000,
            max_sites: 10_000,
            max_fragments: 10_000,
            max_report_bytes: 4 * 1024 * 1024,
        }
    }
}
dto! {#[derive(PartialEq,Eq)] pub struct EnzymeSpec {pub name:String,pub recognition:String,pub top_cut:i32,pub bottom_cut:i32}}
dto! {#[derive(PartialEq,Eq)] pub struct ReferenceSpan {pub start:usize,pub end:usize}}
enumeration! {pub enum SiteDiagnostic {PossibleRecognition,CutsOutsideSubstrate,AmbiguousStrandGeometry}}
enumeration! {pub enum InspectionLimit {Comparisons,Sites,Fragments,ReportBytes}}
dto! {#[derive(PartialEq,Eq)] pub struct InspectionSite {pub index:usize,pub enzyme_index:usize,pub enzyme:String,pub recognition_start:usize,pub recognition_length:usize,pub reverse:bool,pub certain:bool,pub spans:Vec<ReferenceSpan>,pub top_boundary:i64,pub bottom_boundary:i64,pub top_cut:Option<usize>,pub bottom_cut:Option<usize>,pub diagnostics:Vec<SiteDiagnostic>}}
dto! {pub struct FragmentSummary {pub index:usize,pub top_start:usize,pub bottom_start:i64,pub duplex:chemistry::DuplexSummary}}
dto! {pub struct RestrictionInspection {pub engine_version:String,pub source_role:String,pub source_revision:String,pub options_revision:String,pub geometry_revision:String,pub batch_revision:String,pub source_length:usize,pub topology:Topology,pub enzymes:Vec<EnzymeSpec>,pub options:InspectionOptions,pub sites:Vec<InspectionSite>,pub scan_complete:bool,pub products_determined:bool,pub comparisons:u64,pub limits:Vec<InspectionLimit>,pub ref_top_lengths:Option<Vec<usize>>,pub fragments:Vec<FragmentSummary>,pub diagnostics:Vec<String>}}
pub const RECIPE_KEY: &str = "restriction_inspection:recipe";
pub const MAX_RECIPE_BYTES: usize = 32 * 1024;
// Unverified producer recipe. It records a reproducible request, never imported physical authority.
dto! {pub struct UnverifiedDigestRecipe {pub schema_version:u16,pub engine_version:String,pub source_role:String,pub source_revision:String,pub geometry_revision:String,pub options_revision:String,pub inspection_revision:String,pub batch_revision:String,pub enzymes:Vec<EnzymeSpec>,pub options:InspectionOptions,pub fragment_index:usize}}
dto! {pub struct DigestProduct {pub engine_version:String,pub batch_revision:String,pub inspection_revision:String,pub source_role:String,pub source_revision:String,pub geometry_revision:String,pub options_revision:String,pub fragment_index:usize,pub molecule:Molecule,pub summary:FragmentSummary,pub enzymes:Vec<EnzymeSpec>,pub options:InspectionOptions,pub projection_findings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("restriction inspection: {s}"))
}
struct Bound {
    used: usize,
    max: usize,
}
impl Write for Bound {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if b.len() > self.max.saturating_sub(self.used) {
            return Err(std::io::Error::other("bounded scientific evidence"));
        }
        self.used += b.len();
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn size<T: Serialize>(v: &T, max: usize) -> Result<usize, Error> {
    let mut w = Bound { used: 0, max };
    serde_json::to_writer(&mut w, v).map_err(|_| Error::TooLarge)?;
    Ok(w.used)
}
fn hash<T: Serialize>(v: &T) -> Result<String, Error> {
    Ok(
        blake3::hash(&serde_json::to_vec(v).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'a'..=b'f'))
}
fn specs(enzymes: &[Enzyme]) -> Vec<EnzymeSpec> {
    enzymes
        .iter()
        .map(|e| EnzymeSpec {
            name: e.name.clone(),
            recognition: e.recognition.clone(),
            top_cut: e.top_cut,
            bottom_cut: e.bottom_cut,
        })
        .collect()
}
fn enzymes(specs: &[EnzymeSpec]) -> Vec<Enzyme> {
    specs
        .iter()
        .map(|e| Enzyme {
            name: e.name.clone(),
            recognition: e.recognition.clone(),
            top_cut: e.top_cut,
            bottom_cut: e.bottom_cut,
        })
        .collect()
}
fn role_valid(role: &str) -> bool {
    !role.is_empty() && role.len() <= 256 && role.bytes().all(|b| b.is_ascii_graphic() || b == b' ')
}
impl InspectionOptions {
    pub fn validate(&self) -> Result<(), Error> {
        if self.max_comparisons == 0
            || self.max_comparisons > 1_000_000_000
            || self.max_sites == 0
            || self.max_sites > 10_000
            || self.max_fragments == 0
            || self.max_fragments > 10_000
            || self.max_report_bytes < 4096
            || self.max_report_bytes > MAX_REPORT_BYTES
        {
            return Err(bad("invalid bounds"));
        }
        Ok(())
    }
}
fn geometry_valid(enzymes: &[Enzyme]) -> Result<(), Error> {
    if enzymes.is_empty() || enzymes.len() > 64 {
        return Err(bad("select 1–64 enzymes"));
    }
    let mut names = BTreeSet::new();
    for e in enzymes {
        if e.name.is_empty()
            || e.name.len() > 64
            || !e.name.bytes().all(|b| b.is_ascii_graphic() || b == b' ')
            || !names.insert(e.name.to_ascii_lowercase())
            || e.recognition.is_empty()
            || e.recognition.len() > 64
        {
            return Err(bad("invalid or duplicate enzyme definition"));
        }
        Sequence::dna(&e.recognition)?;
    }
    Ok(())
}
fn loc(l: &Location, depth: usize, nodes: &mut usize) -> Result<(), Error> {
    *nodes = nodes.saturating_add(1);
    if depth > 32 || *nodes > 100_000 {
        return Err(Error::TooLarge);
    }
    match l {
        Location::Complement(x) => loc(x, depth + 1, nodes)?,
        Location::Join(xs) | Location::Order(xs) => {
            for x in xs {
                loc(x, depth + 1, nodes)?;
            }
        }
        _ => (),
    }
    Ok(())
}
/// Bounds rich sources before validation, serialization allocation or source hashing.
pub fn source_preflight(
    role: &str,
    m: &Molecule,
    es: &[Enzyme],
    o: &InspectionOptions,
) -> Result<(), Error> {
    o.validate()?;
    geometry_valid(es)?;
    if !role_valid(role) {
        return Err(bad("invalid source role"));
    }
    if m.sequence.len() > crate::model::MAX_SEQUENCE_LEN
        || m.features.len() > 100_000
        || m.metadata.len() > 100_000
    {
        return Err(Error::TooLarge);
    }
    let mut nodes = m.features.len().saturating_add(m.metadata.len());
    for f in &m.features {
        nodes = nodes.saturating_add(f.qualifiers.len());
        loc(&f.location, 0, &mut nodes)?;
    }
    if nodes > 100_000 {
        return Err(Error::TooLarge);
    }
    size(m, MAX_SOURCE_BYTES)?;
    m.validate()?;
    if m.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("requires DNA"));
    }
    Ok(())
}
pub fn batch_revision(
    role: &str,
    m: &Molecule,
    es: &[Enzyme],
    o: &InspectionOptions,
) -> Result<String, Error> {
    source_preflight(role, m, es, o)?;
    hash(&(ENGINE_VERSION, role, m.revision(), specs(es), o))
}
pub fn product_revision(
    role: &str,
    m: &Molecule,
    es: &[Enzyme],
    o: &InspectionOptions,
    index: usize,
) -> Result<String, Error> {
    if index >= o.max_fragments {
        return Err(bad("fragment index exceeds bound"));
    }
    hash(&(
        "restriction-product-1",
        batch_revision(role, m, es, o)?,
        index,
    ))
}
fn spans(start: usize, len: usize, n: usize) -> Vec<ReferenceSpan> {
    let first = len.min(n - start);
    let mut out = vec![ReferenceSpan {
        start,
        end: start + first,
    }];
    if first < len {
        out.push(ReferenceSpan {
            start: 0,
            end: len - first,
        });
    }
    out
}
fn boundary(raw: i64, n: usize, topology: Topology) -> Option<usize> {
    if topology == Topology::Circular {
        Some(raw.rem_euclid(n as i64) as usize)
    } else if (0..=n as i64).contains(&raw) {
        Some(raw as usize)
    } else {
        None
    }
}
fn site(
    e: &Enzyme,
    ei: usize,
    index: usize,
    start: usize,
    reverse: bool,
    certain: bool,
    n: usize,
    t: Topology,
) -> InspectionSite {
    let len = e.recognition.len();
    let (a, b) = if reverse {
        (
            start as i64 + len as i64 - i64::from(e.bottom_cut),
            start as i64 + len as i64 - i64::from(e.top_cut),
        )
    } else {
        (
            start as i64 + i64::from(e.top_cut),
            start as i64 + i64::from(e.bottom_cut),
        )
    };
    let mut top_cut = boundary(a, n, t);
    let mut bottom_cut = boundary(b, n, t);
    let mut diagnostics = vec![];
    if !certain {
        diagnostics.push(SiteDiagnostic::PossibleRecognition);
    }
    if top_cut.is_none() || bottom_cut.is_none() {
        diagnostics.push(SiteDiagnostic::CutsOutsideSubstrate);
    }
    // A palindromic footprint cannot select one of two different strand-cut
    // orientations. Retain the declared raw boundaries, not a purported unique cut.
    let palindrome = Sequence::dna(&e.recognition)
        .and_then(|p| Ok(p.reverse_complement()? == p))
        .unwrap_or(false);
    if palindrome && i64::from(e.top_cut) + i64::from(e.bottom_cut) != len as i64 {
        diagnostics.push(SiteDiagnostic::AmbiguousStrandGeometry);
        top_cut = None;
        bottom_cut = None;
    }
    InspectionSite {
        index,
        enzyme_index: ei,
        enzyme: e.name.clone(),
        recognition_start: start,
        recognition_length: len,
        reverse,
        certain,
        spans: spans(start, len, n),
        top_boundary: a,
        bottom_boundary: b,
        top_cut,
        bottom_cut,
        diagnostics,
    }
}
fn summary(index: usize, f: &chemistry::DigestFragment) -> Result<FragmentSummary, Error> {
    Ok(FragmentSummary {
        index,
        top_start: f.top_start,
        bottom_start: f.bottom_start,
        duplex: f.duplex.summary()?,
    })
}
fn push_limit(r: &mut RestrictionInspection, l: InspectionLimit) {
    if !r.limits.contains(&l) {
        r.limits.push(l);
    }
}
/// Recognition evidence only, with explicit bounded partial results. This never
/// reconstructs a duplex, computes fragment lengths, or manufactures products.
/// The same source/geometry/options bindings are used by the full inspection.
pub fn scan(
    role: &str,
    m: &Molecule,
    es: &[Enzyme],
    o: &InspectionOptions,
) -> Result<RestrictionInspection, Error> {
    source_preflight(role, m, es, o)?;
    let mut r = RestrictionInspection {
        engine_version: ENGINE_VERSION.into(),
        source_role: role.into(),
        source_revision: m.revision(),
        options_revision: hash(o)?,
        geometry_revision: hash(&specs(es))?,
        batch_revision: batch_revision(role, m, es, o)?,
        source_length: m.sequence.len(),
        topology: m.topology,
        enzymes: specs(es),
        options: o.clone(),
        sites: vec![],
        scan_complete: true,
        products_determined: false,
        comparisons: 0,
        limits: vec![],
        ref_top_lengths: None,
        fragments: vec![],
        diagnostics: vec![INSPECTION_ASSUMPTIONS.into()],
    };
    let n = m.sequence.len();
    let bases = m.sequence.as_str().as_bytes();
    let mut bytes = size(&r, o.max_report_bytes)?; // fixed overhead must fit before scanning
    'scan: for (ei, e) in es.iter().enumerate() {
        let p = Sequence::dna(&e.recognition)?;
        if p.len() > n {
            continue;
        }
        let rc = p.reverse_complement()?;
        let mut patterns = vec![(p.as_str(), false)];
        if rc != p {
            patterns.push((rc.as_str(), true));
        }
        for (pattern, reverse) in patterns {
            let end = if m.topology == Topology::Circular {
                n
            } else {
                n - p.len() + 1
            };
            for start in 0..end {
                let (mut certain, mut matched) = (true, true);
                for (j, b) in pattern.bytes().enumerate() {
                    if r.comparisons == o.max_comparisons {
                        push_limit(&mut r, InspectionLimit::Comparisons);
                        r.scan_complete = false;
                        break 'scan;
                    }
                    r.comparisons += 1;
                    let a = analysis::base_mask(bases[(start + j) % n]);
                    let b = analysis::base_mask(b);
                    if a & b == 0 {
                        matched = false;
                        break;
                    }
                    if a & !b != 0 {
                        certain = false;
                    }
                }
                if matched {
                    if r.sites.len() == o.max_sites {
                        push_limit(&mut r, InspectionLimit::Sites);
                        r.scan_complete = false;
                        break 'scan;
                    }
                    let s = site(e, ei, r.sites.len(), start, reverse, certain, n, m.topology);
                    let additional = size(&s, o.max_report_bytes)?.saturating_add(1);
                    if additional
                        > o.max_report_bytes
                            .saturating_sub(bytes)
                            .saturating_sub(1024)
                    {
                        push_limit(&mut r, InspectionLimit::ReportBytes);
                        r.scan_complete = false;
                        break 'scan;
                    }
                    bytes += additional;
                    r.sites.push(s);
                }
            }
        }
    }
    // Keep the historical recognition-loop reserve unchanged for inspect(),
    // then describe this narrower API without asserting physical preparation.
    r.diagnostics = vec![SCAN_ASSUMPTIONS.into()];
    size(&r, o.max_report_bytes)?;
    r.validate()?;
    Ok(r)
}

/// Recognition inspection plus lazily summarized assumed-complete duplex products.
/// Existing report bindings and physical-product diagnostics remain unchanged.
pub fn inspect(
    role: &str,
    m: &Molecule,
    es: &[Enzyme],
    o: &InspectionOptions,
) -> Result<RestrictionInspection, Error> {
    inspect_with_digest(role, m, es, o, || {
        chemistry::Duplex::from_molecule(m).and_then(|d| chemistry::digest(&d, es))
    })
}
/// Private fresh-native composition only. Callers must reconstruct this substrate,
/// never deserialize a producer's physical summary as input authority.
pub(crate) fn inspect_fresh_duplex(
    role: &str,
    m: &Molecule,
    es: &[Enzyme],
    o: &InspectionOptions,
    duplex: &chemistry::Duplex,
) -> Result<RestrictionInspection, Error> {
    duplex.validate()?;
    if duplex.top != m.sequence || duplex.topology != m.topology {
        return Err(bad("fresh duplex/reference binding differs"));
    }
    inspect_with_digest(role, m, es, o, || chemistry::digest(duplex, es))
}
fn inspect_with_digest(
    role: &str,
    m: &Molecule,
    es: &[Enzyme],
    o: &InspectionOptions,
    digest: impl FnOnce() -> Result<Vec<chemistry::DigestFragment>, Error>,
) -> Result<RestrictionInspection, Error> {
    let mut r = scan(role, m, es, o)?;
    r.diagnostics = vec![INSPECTION_ASSUMPTIONS.into()];
    if r.scan_complete
        && r.sites
            .iter()
            .all(|s| s.certain && s.top_cut.is_some() && s.bottom_cut.is_some())
    {
        let cuts: Vec<_> = r
            .sites
            .iter()
            .map(|s| analysis::CutSite {
                enzyme: s.enzyme.clone(),
                recognition_start: s.recognition_start,
                reverse: s.reverse,
                top_cut: s.top_cut.unwrap(),
                bottom_cut: s.bottom_cut.unwrap(),
                certain: s.certain,
            })
            .collect();
        r.ref_top_lengths = Some(analysis::digest_lengths(m, &cuts)?);
        match digest() {
            Ok(fs) => {
                if fs.len() > o.max_fragments {
                    push_limit(&mut r, InspectionLimit::Fragments);
                    r.diagnostics.push("Complete physical fragment set exceeds the selected fragment bound; no product set is asserted.".into());
                } else {
                    r.fragments = fs
                        .iter()
                        .enumerate()
                        .map(|(i, f)| summary(i, f))
                        .collect::<Result<_, _>>()?;
                    r.products_determined = true;
                }
            }
            Err(e) => r
                .diagnostics
                .push(format!("Physical products withheld: {e}")),
        }
    } else {
        r.diagnostics.push("Physical products withheld: recognition scan is partial or contains possible sites/unavailable cuts.".into());
    }
    if size(&r, o.max_report_bytes).is_err() {
        r.fragments.clear();
        r.ref_top_lengths = None;
        r.products_determined = false;
        push_limit(&mut r, InspectionLimit::ReportBytes);
        r.diagnostics =
            vec!["Report output bound reached; product and length summaries withheld.".into()];
        while size(&r, o.max_report_bytes).is_err() {
            if r.sites.pop().is_none() {
                return Err(Error::TooLarge);
            }
            r.scan_complete = false;
        }
    }
    r.validate()?;
    Ok(r)
}
/// Recompute the current science and project source annotations into one lazily requested product.
pub fn product(
    role: &str,
    m: &Molecule,
    es: &[Enzyme],
    o: &InspectionOptions,
    index: usize,
) -> Result<DigestProduct, Error> {
    let revision = product_revision(role, m, es, o, index)?;
    let report = inspect(role, m, es, o)?;
    if !report.products_determined {
        return Err(bad(
            "physical product set is not determined under this inspection",
        ));
    }
    let fs = chemistry::digest(&chemistry::Duplex::from_molecule(m)?, es)?;
    let f = fs
        .get(index)
        .ok_or_else(|| bad("fragment index does not exist"))?;
    let base: String = m
        .name
        .chars()
        .take(48)
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.') {
                c
            } else {
                '_'
            }
        })
        .collect();
    let mut molecule = f
        .duplex
        .molecule(format!("{base}_digest_fragment_{}", index + 1))?;
    let mut findings = vec![];
    crate::verification::project_digest_annotations(m, f, &mut molecule, &mut findings)?;
    molecule.metadata.insert(
        "restriction_inspection:engine".into(),
        ENGINE_VERSION.into(),
    );
    molecule.metadata.insert(
        "restriction_inspection:batch_revision".into(),
        revision.clone(),
    );
    molecule.metadata.insert(
        "restriction_inspection:source_revision".into(),
        report.source_revision.clone(),
    );
    molecule.metadata.insert("restriction_inspection:assumption".into(),"Assumed complete cleavage; producer recipe is descriptive and cannot restore physical authority.".into());
    let recipe = UnverifiedDigestRecipe {
        schema_version: 1,
        engine_version: ENGINE_VERSION.into(),
        source_role: role.into(),
        source_revision: report.source_revision.clone(),
        geometry_revision: report.geometry_revision.clone(),
        options_revision: report.options_revision.clone(),
        inspection_revision: report.batch_revision.clone(),
        batch_revision: revision.clone(),
        enzymes: specs(es),
        options: o.clone(),
        fragment_index: index,
    };
    recipe.validate()?;
    molecule.metadata.insert(
        RECIPE_KEY.into(),
        serde_json::to_string(&recipe).map_err(|_| Error::TooLarge)?,
    );
    let p = DigestProduct {
        engine_version: ENGINE_VERSION.into(),
        batch_revision: revision,
        inspection_revision: report.batch_revision,
        source_role: role.into(),
        source_revision: report.source_revision,
        geometry_revision: report.geometry_revision,
        options_revision: report.options_revision,
        fragment_index: index,
        molecule,
        summary: summary(index, f)?,
        enzymes: specs(es),
        options: o.clone(),
        projection_findings: findings,
    };
    p.validate()?;
    size(&p, MAX_SOURCE_BYTES)?;
    Ok(p)
}
fn duplex_same(a: &chemistry::DuplexSummary, b: &chemistry::DuplexSummary) -> bool {
    a.model == b.model
        && a.top_length == b.top_length
        && a.bottom_length == b.bottom_length
        && a.bottom_start == b.bottom_start
        && a.topology == b.topology
        && a.ends == b.ends
        && a.nicks == b.nicks
        && a.covalently_closed_in_model == b.covalently_closed_in_model
}
/// Validate borrowed descriptor shapes before any conversion, clone or hash.
fn specs_valid(es: &[EnzymeSpec]) -> Result<(), Error> {
    if es.is_empty()
        || es.len() > 64
        || es.iter().any(|e| {
            e.name.is_empty()
                || e.name.len() > 64
                || e.recognition.is_empty()
                || e.recognition.len() > 64
        })
    {
        return Err(bad("invalid enzyme descriptor bounds"));
    }
    geometry_valid(&enzymes(es))
}
fn version(s: &str) -> bool {
    !s.is_empty() && s.len() <= 128 && s.bytes().all(|b| b.is_ascii_graphic())
}
fn summary_valid(s: &FragmentSummary, n: usize) -> Result<(), Error> {
    let d = &s.duplex;
    if s.top_start >= n
        || d.top_length == 0
        || d.top_length > n
        || d.bottom_length == 0
        || d.bottom_length > n
        || d.bottom_start.unsigned_abs() > n as u64
        || s.bottom_start.unsigned_abs() > 2 * n as u64
        || d.nicks.len() > 10_000
        || d.model.len() > 256
        || d.ends
            .as_ref()
            .is_some_and(|(a, b)| a.overhang.len() > n || b.overhang.len() > n)
    {
        return Err(bad("invalid fragment summary"));
    }
    if (d.topology == Topology::Circular) != d.ends.is_none() {
        return Err(bad("fragment topology/end mismatch"));
    }
    Ok(())
}
impl RestrictionInspection {
    /// Bounded structural validation only. Historical reports are evidence, never science authority.
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        specs_valid(&self.enzymes)?;
        if !version(&self.engine_version)
            || !role_valid(&self.source_role)
            || self.source_length == 0
            || self.source_length > crate::model::MAX_SEQUENCE_LEN
            || ![
                &self.source_revision,
                &self.options_revision,
                &self.geometry_revision,
                &self.batch_revision,
            ]
            .iter()
            .all(|s| hex(s))
            || self.options_revision != hash(&self.options)?
            || self.geometry_revision != hash(&self.enzymes)?
            || self.batch_revision
                != hash(&(
                    &self.engine_version,
                    &self.source_role,
                    &self.source_revision,
                    &self.enzymes,
                    &self.options,
                ))?
            || self.comparisons > self.options.max_comparisons
            || self.sites.len() > self.options.max_sites
            || self.fragments.len() > self.options.max_fragments
            || self.limits.len() > 4
            || self.diagnostics.len() > 16
            || self.diagnostics.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("invalid report shape/bindings"));
        }
        let mut limits = BTreeSet::new();
        if self.limits.iter().any(|l| !limits.insert(*l as u8)) {
            return Err(bad("duplicate limits"));
        }
        for (i, s) in self.sites.iter().enumerate() {
            let e = self
                .enzymes
                .get(s.enzyme_index)
                .ok_or_else(|| bad("site enzyme index"))?;
            if s.index != i
                || s.spans.is_empty()
                || s.spans.len() > 2
                || s.diagnostics.len() > 3
                || s.enzyme != e.name
                || s.recognition_start >= self.source_length
                || s.recognition_length != e.recognition.len()
                || s.recognition_length > self.source_length
                || self.topology == Topology::Linear
                    && s.recognition_start + s.recognition_length > self.source_length
            {
                return Err(bad("invalid recognition evidence"));
            }
            let expected = site(
                &enzymes(std::slice::from_ref(e))[0],
                s.enzyme_index,
                i,
                s.recognition_start,
                s.reverse,
                s.certain,
                self.source_length,
                self.topology,
            );
            if s != &expected {
                return Err(bad("site coordinates/diagnostics disagree with geometry"));
            }
        }
        if !self.scan_complete && self.limits.is_empty()
            || self.scan_complete
                && self
                    .limits
                    .iter()
                    .any(|l| matches!(l, InspectionLimit::Comparisons | InspectionLimit::Sites))
        {
            return Err(bad("scan completeness disagrees with limits"));
        }
        let eligible = self.scan_complete
            && self
                .sites
                .iter()
                .all(|s| s.certain && s.top_cut.is_some() && s.bottom_cut.is_some());
        if !eligible && (self.ref_top_lengths.is_some() || self.products_determined) {
            return Err(bad(
                "partial/possible evidence cannot establish products or lengths",
            ));
        }
        if let Some(lengths) = &self.ref_top_lengths {
            if lengths.is_empty()
                || lengths.len() > self.options.max_sites + 1
                || lengths.iter().any(|l| *l == 0 || *l > self.source_length)
                || lengths.iter().sum::<usize>() != self.source_length
            {
                return Err(bad("invalid reference lengths"));
            }
        }
        if self.products_determined {
            if self.fragments.is_empty()
                || !self.limits.is_empty()
                || self.ref_top_lengths.is_none()
            {
                return Err(bad("invalid product completeness"));
            }
            let mut total = 0usize;
            let mut bottom_total = 0usize;
            for (i, s) in self.fragments.iter().enumerate() {
                if s.index != i {
                    return Err(bad("fragment order"));
                }
                summary_valid(s, self.source_length)?;
                total = total.saturating_add(s.duplex.top_length);
                bottom_total = bottom_total.saturating_add(s.duplex.bottom_length);
            }
            if total != self.source_length || bottom_total != self.source_length {
                return Err(bad("fragment strand accounting"));
            }
        } else if !self.fragments.is_empty() {
            return Err(bad("undetermined product summaries"));
        }
        size(self, self.options.max_report_bytes)?;
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        role: &str,
        m: &Molecule,
        es: &[Enzyme],
        o: &InspectionOptions,
    ) -> Result<(), Error> {
        source_preflight(role, m, es, o)?;
        self.validate()?;
        if self.engine_version != ENGINE_VERSION
            || self.source_role != role
            || self.source_revision != m.revision()
            || self.source_length != m.sequence.len()
            || self.topology != m.topology
            || self.enzymes != specs(es)
            || self.options != *o
            || self.batch_revision != batch_revision(role, m, es, o)?
        {
            return Err(bad("stale report/source/options/geometry binding"));
        }
        Ok(())
    }
}
impl DigestProduct {
    /// Shape and producer binding checks only; serialized recipes do not restore physical authority.
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        specs_valid(&self.enzymes)?;
        if !version(&self.engine_version)
            || !role_valid(&self.source_role)
            || self.fragment_index >= self.options.max_fragments
            || self.summary.index != self.fragment_index
            || ![
                &self.source_revision,
                &self.options_revision,
                &self.geometry_revision,
                &self.batch_revision,
                &self.inspection_revision,
            ]
            .iter()
            .all(|s| hex(s))
            || self.options_revision != hash(&self.options)?
            || self.geometry_revision != hash(&self.enzymes)?
            || self.inspection_revision
                != hash(&(
                    &self.engine_version,
                    &self.source_role,
                    &self.source_revision,
                    &self.enzymes,
                    &self.options,
                ))?
            || self.batch_revision
                != hash(&(
                    "restriction-product-1",
                    &self.inspection_revision,
                    self.fragment_index,
                ))?
            || self.projection_findings.len() > 101
            || self.projection_findings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad("invalid product shape/binding"));
        }
        let m = &self.molecule;
        let d = &self.summary.duplex;
        if m.sequence.len() != d.top_length
            || m.topology != d.topology
            || m.sequence.alphabet() != Alphabet::Dna
            || m.original_genbank.is_some()
            || m.metadata.contains_key("embl.original_record")
        {
            return Err(bad("derived product/evidence contradiction"));
        }
        source_preflight(&self.source_role, m, &enzymes(&self.enzymes), &self.options)?;
        summary_valid(&self.summary, crate::model::MAX_SEQUENCE_LEN)?;
        let metadata = |key: &str, value: &str| m.metadata.get(key).is_some_and(|v| v == value);
        if !metadata("restriction_inspection:engine", &self.engine_version)
            || !metadata(
                "restriction_inspection:batch_revision",
                &self.batch_revision,
            )
            || !metadata(
                "restriction_inspection:source_revision",
                &self.source_revision,
            )
        {
            return Err(bad("product metadata binding"));
        }
        let recipe =
            unverified_recipe(m)?.ok_or_else(|| bad("missing unverified producer recipe"))?;
        if recipe.engine_version != self.engine_version
            || recipe.source_role != self.source_role
            || recipe.source_revision != self.source_revision
            || recipe.geometry_revision != self.geometry_revision
            || recipe.options_revision != self.options_revision
            || recipe.inspection_revision != self.inspection_revision
            || recipe.batch_revision != self.batch_revision
            || recipe.enzymes != self.enzymes
            || recipe.options != self.options
            || recipe.fragment_index != self.fragment_index
        {
            return Err(bad("product recipe binding contradiction"));
        }
        let raw = m
            .metadata
            .get("duplex_model")
            .ok_or_else(|| bad("missing duplex summary"))?;
        let summary: chemistry::DuplexSummary =
            serde_json::from_str(raw).map_err(|_| bad("invalid duplex summary"))?;
        if !duplex_same(&summary, d) {
            return Err(bad("product summary contradiction"));
        }
        size(self, MAX_SOURCE_BYTES)?;
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        role: &str,
        m: &Molecule,
        es: &[Enzyme],
        o: &InspectionOptions,
        index: usize,
    ) -> Result<(), Error> {
        source_preflight(role, m, es, o)?;
        self.validate()?;
        summary_valid(&self.summary, m.sequence.len())?;
        if self.engine_version != ENGINE_VERSION
            || self.source_role != role
            || self.source_revision != m.revision()
            || self.fragment_index != index
            || self.enzymes != specs(es)
            || self.options != *o
            || self.batch_revision != product_revision(role, m, es, o, index)?
        {
            return Err(bad("stale product/source/options/geometry binding"));
        }
        Ok(())
    }
}

impl UnverifiedDigestRecipe {
    /// Shape/hash binding only; importing this value does not establish cleavage or chemistry.
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        specs_valid(&self.enzymes)?;
        if self.schema_version != 1
            || !version(&self.engine_version)
            || !role_valid(&self.source_role)
            || self.fragment_index >= self.options.max_fragments
            || ![
                &self.source_revision,
                &self.geometry_revision,
                &self.options_revision,
                &self.inspection_revision,
                &self.batch_revision,
            ]
            .iter()
            .all(|s| hex(s))
            || self.geometry_revision != hash(&self.enzymes)?
            || self.options_revision != hash(&self.options)?
            || self.inspection_revision
                != hash(&(
                    &self.engine_version,
                    &self.source_role,
                    &self.source_revision,
                    &self.enzymes,
                    &self.options,
                ))?
            || self.batch_revision
                != hash(&(
                    "restriction-product-1",
                    &self.inspection_revision,
                    self.fragment_index,
                ))?
        {
            return Err(bad("invalid unverified producer recipe"));
        }
        size(self, MAX_RECIPE_BYTES)?;
        Ok(())
    }
}
/// Read bounded native producer data only. Call `product` with the exact source for fresh science.
pub fn unverified_recipe(m: &Molecule) -> Result<Option<UnverifiedDigestRecipe>, Error> {
    let Some(raw) = m.metadata.get(RECIPE_KEY) else {
        return Ok(None);
    };
    if raw.len() > MAX_RECIPE_BYTES {
        return Err(Error::TooLarge);
    }
    let recipe: UnverifiedDigestRecipe =
        serde_json::from_str(raw).map_err(|_| bad("invalid producer recipe JSON"))?;
    recipe.validate()?;
    Ok(Some(recipe))
}
