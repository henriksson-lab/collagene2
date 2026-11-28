//! Bounded one-copy ligation paths from actual native tracked substrates.
//! Intended geometry is never selected by a target or ranked by probability.
use crate::{
    chemistry::Duplex,
    modified_duplex::{self, ModificationSummary, TrackedDuplex},
};
use collagene_model::{Error, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ENGINE_VERSION: &str = "ligation-outcomes-1";
pub const MAX_FRAGMENTS: usize = 1000;
pub const MAX_ENUMERATION_FRAGMENTS: usize = 8;
const MAX_SOURCE_BASES: usize = 2_000_000;
const MAX_BYTES: usize = 32 * 1024 * 1024;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)]pub struct Placement{pub key:String,pub reverse:bool}}
dto! {pub struct Options{pub max_enumeration_fragments:usize,pub max_path_checks:usize,pub max_products:usize,pub max_total_product_bases:usize,pub max_report_bytes:usize,pub max_work:usize}}
impl Default for Options {
    fn default() -> Self {
        Self {
            max_enumeration_fragments: MAX_ENUMERATION_FRAGMENTS,
            max_path_checks: 4096,
            max_products: 100,
            max_total_product_bases: 8_000_000,
            max_report_bytes: MAX_BYTES,
            max_work: 20_000_000,
        }
    }
}
dto! {pub struct Request{pub intended:Vec<Placement>,pub circular:bool,pub allow_nicked:bool,pub options:Options}}
/// Never deserialize a report to construct these operands. The caller retains
/// exact source events and supplies freshly reconstructed native material.
pub struct Operand<'a> {
    pub key: &'a str,
    pub material: &'a TrackedDuplex,
}
dto! {pub struct Source{pub key:String,pub state_revision:String,pub top_length:usize,pub bottom_length:usize}}
dto! {pub struct Product{pub path:Vec<Placement>,pub sequence:Sequence,pub topology:Topology,pub duplex:Duplex,pub modifications:ModificationSummary}}
dto! {#[derive(Copy,PartialEq,Eq)]pub enum Limit{PathChecks,Products,ProductBases,ReportBytes,Work,NativeMaterial,EnumerationFragments}}
dto! {pub struct UnresolvedPath{pub path:Vec<Placement>,pub reason:String}}
dto! {pub struct Report{pub engine_version:String,pub scientific_engine:String,pub material_engine:String,pub context_revision:String,pub request:Request,pub sources:Vec<Source>,pub intended:Product,pub alternatives:Vec<Product>,pub paths_checked:usize,pub incompatible_paths:usize,pub unresolved_paths:usize,pub unresolved_examples:Vec<UnresolvedPath>,pub unresolved_examples_truncated:bool,pub withheld_products:usize,pub work:usize,pub retained_product_bases:usize,pub enumeration_complete:bool,pub exhaustive_within_scope:bool,pub limits:Vec<Limit>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("Ligation outcomes: {s}"))
}
fn count(value: &impl Serialize, limit: usize) -> Result<usize, Error> {
    struct Count {
        n: usize,
        limit: usize,
    }
    impl Write for Count {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.n = self
                .n
                .checked_add(b.len())
                .filter(|n| *n <= self.limit)
                .ok_or_else(|| std::io::Error::other("bounded ligation outcomes"))?;
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut c = Count { n: 0, limit };
    serde_json::to_writer(&mut c, value).map_err(|_| Error::TooLarge)?;
    Ok(c.n)
}
impl Options {
    pub fn validate(&self) -> Result<(), Error> {
        if self.max_enumeration_fragments == 0
            || self.max_enumeration_fragments > MAX_ENUMERATION_FRAGMENTS
            || self.max_path_checks == 0
            || self.max_path_checks > 100_000
            || self.max_products == 0
            || self.max_products > 1000
            || self.max_total_product_bases == 0
            || self.max_total_product_bases > 32_000_000
            || !(4096..=MAX_BYTES).contains(&self.max_report_bytes)
            || self.max_work == 0
            || self.max_work > 100_000_000
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
fn shape(request: &Request) -> Result<(), Error> {
    request.options.validate()?;
    if request.intended.is_empty() || request.intended.len() > MAX_FRAGMENTS {
        return Err(Error::TooLarge);
    }
    let mut keys = BTreeSet::new();
    for p in &request.intended {
        if p.key.is_empty()
            || p.key.len() > 128
            || p.key.chars().any(char::is_control)
            || !keys.insert(&p.key)
        {
            return Err(bad(
                "intended path must use each unique source key exactly once",
            ));
        }
    }
    Ok(())
}
fn sources<'a>(
    operands: &'a [Operand<'a>],
    request: &Request,
) -> Result<(Vec<Source>, Vec<&'a TrackedDuplex>, usize), Error> {
    shape(request)?;
    if operands.len() != request.intended.len() {
        return Err(bad(
            "intended path must include every declared operand exactly once",
        ));
    }
    let mut map = BTreeMap::new();
    let mut bases = 0usize;
    let mut source_bytes = 0usize;
    for operand in operands {
        if operand.key.is_empty()
            || operand.key.len() > 128
            || operand.key.chars().any(char::is_control)
            || map.insert(operand.key, operand.material).is_some()
        {
            return Err(bad("duplicate or invalid source key"));
        }
        bases = bases
            .checked_add(operand.material.duplex.top.len())
            .and_then(|n| n.checked_add(operand.material.duplex.bottom.len()))
            .filter(|n| *n <= MAX_SOURCE_BASES)
            .ok_or(Error::TooLarge)?;
        source_bytes = source_bytes
            .checked_add(count(
                operand.material,
                MAX_BYTES.saturating_sub(source_bytes),
            )?)
            .ok_or(Error::TooLarge)?;
    }
    if request
        .intended
        .iter()
        .any(|p| !map.contains_key(p.key.as_str()))
    {
        return Err(bad("unknown intended source key"));
    }
    let work = bases
        .checked_mul(4)
        .filter(|n| *n <= request.options.max_work)
        .ok_or(Error::TooLarge)?;
    let mut descriptors = vec![];
    let mut materials = vec![];
    for (key, material) in map {
        material.validate()?;
        if material.duplex.topology != Topology::Linear {
            return Err(Error::Format(
                "Ligation outcomes: operands must be linear duplex fragments".into(),
            ));
        }
        let summary = modified_duplex::summary(material)?;
        descriptors.push(Source {
            key: key.into(),
            state_revision: summary.state_revision,
            top_length: material.duplex.top.len(),
            bottom_length: material.duplex.bottom.len(),
        });
        materials.push(material);
    }
    Ok((descriptors, materials, work))
}
/// Full native material identities plus literal path/bounds; no target input.
pub fn revision(operands: &[Operand<'_>], request: &Request) -> Result<String, Error> {
    let (sources, _, _) = sources(operands, request)?;
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        crate::verification::ENGINE_VERSION,
        modified_duplex::ENGINE_VERSION,
        sources,
        request,
    ))
}
fn context(
    sources: &[Source],
    request: &Request,
    science: &str,
    material: &str,
) -> Result<String, Error> {
    crate::assembly_recipe::hash(&(ENGINE_VERSION, science, material, sources, request))
}
struct Search<'a> {
    keys: Vec<String>,
    materials: Vec<[TrackedDuplex; 2]>,
    request: &'a Request,
    intended: Vec<(usize, bool)>,
    report: Report,
    stopped: bool,
}
impl Search<'_> {
    fn limit(&mut self, limit: Limit) {
        if !self.report.limits.contains(&limit) {
            self.report.limits.push(limit);
        }
        self.stopped = true;
    }
    fn path(&self, p: &[(usize, bool)]) -> Vec<Placement> {
        p.iter()
            .map(|&(i, reverse)| Placement {
                key: self.keys[i].clone(),
                reverse,
            })
            .collect()
    }
    fn footprint(&self, p: &[(usize, bool)]) -> Option<usize> {
        p.iter().try_fold(0usize, |n, &(i, r)| {
            let d = &self.materials[i][usize::from(r)].duplex;
            n.checked_add(d.top.len().checked_mul(2)?)
                .and_then(|n| n.checked_add(d.bottom.len()))
        })
    }
    fn compute(&mut self, p: &[(usize, bool)]) -> Result<Product, Error> {
        let source_bases = p
            .iter()
            .map(|&(i, r)| {
                let d = &self.materials[i][usize::from(r)].duplex;
                d.top.len() + d.bottom.len()
            })
            .sum::<usize>();
        let work = source_bases.checked_mul(6).ok_or(Error::TooLarge)?;
        self.report.work = self
            .report
            .work
            .checked_add(work)
            .filter(|n| *n <= self.request.options.max_work)
            .ok_or(Error::TooLarge)?;
        let materials = p
            .iter()
            .map(|&(i, r)| self.materials[i][usize::from(r)].clone())
            .collect::<Vec<_>>();
        let actual =
            modified_duplex::ligate(&materials, self.request.circular, self.request.allow_nicked)?;
        let modifications = modified_duplex::summary(&actual)?;
        Ok(Product {
            path: self.path(p),
            sequence: actual.duplex.top.clone(),
            topology: actual.duplex.topology,
            duplex: actual.duplex,
            modifications,
        })
    }
    fn leaf(&mut self, path: &[(usize, bool)]) {
        if path == self.intended.as_slice() {
            return;
        }
        if self.report.paths_checked >= self.request.options.max_path_checks {
            self.limit(Limit::PathChecks);
            return;
        }
        if self.report.alternatives.len() + 1 >= self.request.options.max_products {
            self.limit(Limit::Products);
            return;
        }
        let Some(footprint) = self.footprint(path) else {
            self.limit(Limit::ProductBases);
            return;
        };
        if footprint
            > self
                .request
                .options
                .max_total_product_bases
                .saturating_sub(self.report.retained_product_bases)
        {
            self.limit(Limit::ProductBases);
            return;
        }
        let work = path
            .iter()
            .map(|&(i, r)| {
                let d = &self.materials[i][usize::from(r)].duplex;
                (d.top.len() + d.bottom.len()) * 6
            })
            .sum::<usize>();
        if work
            > self
                .request
                .options
                .max_work
                .saturating_sub(self.report.work)
        {
            self.limit(Limit::Work);
            return;
        }
        self.report.paths_checked += 1;
        match self.compute(path) {
            Ok(product) => {
                self.report.alternatives.push(product);
                if count(&self.report, self.request.options.max_report_bytes - 512).is_err() {
                    self.report.alternatives.pop();
                    self.report.withheld_products += 1;
                    self.limit(Limit::ReportBytes);
                } else {
                    self.report.retained_product_bases += footprint;
                }
            }
            Err(Error::Format(_))
            | Err(Error::Coordinates(_))
            | Err(Error::InvalidBase { .. })
            | Err(Error::EmptySequence) => self.report.incompatible_paths += 1,
            Err(Error::TooLarge) => {
                self.report.unresolved_paths += 1;
                self.limit(Limit::NativeMaterial);
            }
            Err(error) => {
                self.report.unresolved_paths += 1;
                let example = UnresolvedPath {
                    path: self.path(path),
                    reason: error.to_string().chars().take(256).collect(),
                };
                if self.report.unresolved_examples.len() < 8 {
                    self.report.unresolved_examples.push(example);
                    if count(&self.report, self.request.options.max_report_bytes - 512).is_err() {
                        self.report.unresolved_examples.pop();
                        self.report.unresolved_examples_truncated = true;
                    }
                } else {
                    self.report.unresolved_examples_truncated = true;
                }
            }
        }
    }
    fn enumerate(&mut self, path: &mut Vec<(usize, bool)>, used: &mut [bool]) {
        if self.stopped {
            return;
        }
        if path.len() == self.keys.len() {
            self.leaf(path);
            return;
        }
        for i in 0..self.keys.len() {
            if used[i] {
                continue;
            }
            used[i] = true;
            for reverse in [false, true] {
                path.push((i, reverse));
                self.enumerate(path, used);
                path.pop();
                if self.stopped {
                    break;
                }
            }
            used[i] = false;
            if self.stopped {
                break;
            }
        }
    }
}
/// Inspect one explicitly declared one-copy-per-fragment scope. Every compatible
/// returned alternative remains extractable; capped searches never invalidate the
/// intended compatible path or imply isolation/mixture uniqueness.
pub fn analyze(operands: &[Operand<'_>], request: &Request) -> Result<Report, Error> {
    analyze_fresh(operands, request).map(|fresh| fresh.report)
}
/// Only this module can construct fresh output state; serialized evidence never
/// reconstructs it. The verifier consumes the actual native intended value once.
pub(crate) struct FreshOutcomes {
    report: Report,
    intended: TrackedDuplex,
}
impl FreshOutcomes {
    pub(crate) fn into_parts(self) -> (Report, TrackedDuplex) {
        (self.report, self.intended)
    }
}
pub(crate) fn analyze_fresh(
    operands: &[Operand<'_>],
    request: &Request,
) -> Result<FreshOutcomes, Error> {
    let (sources, materials, work) = sources(operands, request)?;
    let context_revision = context(
        &sources,
        request,
        crate::verification::ENGINE_VERSION,
        modified_duplex::ENGINE_VERSION,
    )?;
    let keys = sources.iter().map(|s| s.key.clone()).collect::<Vec<_>>();
    let intended = request
        .intended
        .iter()
        .map(|p| (keys.iter().position(|k| k == &p.key).unwrap(), p.reverse))
        .collect::<Vec<_>>();
    // The intended product is evaluated first, independently of enumeration and
    // regardless of whether an alternative agrees with any future chosen target.
    let work = work
        .checked_add(
            intended
                .iter()
                .map(|&(i, _)| {
                    let d = &materials[i].duplex;
                    (d.top.len() + d.bottom.len()) * 6
                })
                .sum::<usize>(),
        )
        .filter(|n| *n <= request.options.max_work)
        .ok_or(Error::TooLarge)?;
    let expected_bases = intended
        .iter()
        .try_fold(0usize, |n, &(i, r)| {
            let d = &materials[i].duplex;
            let (top, bottom) = if r {
                (d.bottom.len(), d.top.len())
            } else {
                (d.top.len(), d.bottom.len())
            };
            n.checked_add(top.checked_mul(2)?)
                .and_then(|n| n.checked_add(bottom))
        })
        .filter(|n| *n <= request.options.max_total_product_bases)
        .ok_or(Error::TooLarge)?;
    let intended_materials = intended
        .iter()
        .map(|&(i, reverse)| {
            if reverse {
                materials[i].reverse()
            } else {
                Ok(materials[i].clone())
            }
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let actual =
        modified_duplex::ligate(&intended_materials, request.circular, request.allow_nicked)?;
    let footprint = actual
        .duplex
        .top
        .len()
        .checked_mul(2)
        .and_then(|n| n.checked_add(actual.duplex.bottom.len()))
        .filter(|n| *n <= request.options.max_total_product_bases)
        .ok_or(Error::TooLarge)?;
    if footprint != expected_bases {
        return Err(bad(
            "native ligation did not conserve the declared stored strand lengths",
        ));
    }
    let modifications = modified_duplex::summary(&actual)?;
    let product = Product {
        path: request.intended.clone(),
        sequence: actual.duplex.top.clone(),
        topology: actual.duplex.topology,
        duplex: actual.duplex.clone(),
        modifications,
    };
    let mut report=Report{engine_version:ENGINE_VERSION.into(),scientific_engine:crate::verification::ENGINE_VERSION.into(),material_engine:modified_duplex::ENGINE_VERSION.into(),context_revision,request:request.clone(),sources,intended:product,alternatives:vec![],paths_checked:1,incompatible_paths:0,unresolved_paths:0,unresolved_examples:vec![],unresolved_examples_truncated:false,withheld_products:0,work,retained_product_bases:footprint,enumeration_complete:false,exhaustive_within_scope:false,limits:vec![],warnings:vec!["One copy of every declared fragment only. Up to32 operands retain an intended path; alternative enumeration is not attempted above the explicit request fragment limit (at most8). Other copy numbers, vector-only/self-ligation, multimers, partial assemblies, concentrations, kinetics, yields and experimental isolation are outside this scope; no probabilities are assigned.".into(),"Circular paths fix the first intended fragment and its intended orientation as the coordinate anchor. Linear paths retain distinct declared orders/orientations, which may include equivalent whole-product coordinate views or identical DNA with different source paths.".into(),"Alternative existence does not invalidate or replace the intended compatible path. Every product is conditional on the native ligation model and declared nick permission; reported strand/modification summaries are documentary, not imported material authority.".into()]};
    count(&report, request.options.max_report_bytes - 512)?;
    if materials.len() > request.options.max_enumeration_fragments {
        report.limits.push(Limit::EnumerationFragments);
        report.validate_binding(&report.context_revision)?;
        return Ok(FreshOutcomes {
            report,
            intended: actual,
        });
    }
    let mut orientations = Vec::with_capacity(materials.len());
    for material in materials {
        orientations.push([material.clone(), material.reverse()?]);
    }
    let mut search = Search {
        keys,
        materials: orientations,
        request,
        intended,
        report,
        stopped: false,
    };
    let mut used = vec![false; search.keys.len()];
    let mut path = vec![];
    if request.circular {
        let anchor = search.intended[0];
        used[anchor.0] = true;
        path.push(anchor);
    }
    search.enumerate(&mut path, &mut used);
    search.report.enumeration_complete = !search.stopped;
    search.report.exhaustive_within_scope = !search.stopped && search.report.unresolved_paths == 0;
    search
        .report
        .validate_binding(&search.report.context_revision)?;
    Ok(FreshOutcomes {
        report: search.report,
        intended: actual,
    })
}

/// Every ligation retains its native intended product and explicit alternative-search status.
pub fn operation(
    fragments: Vec<crate::verification::OrientedSource>,
    circular: bool,
    allow_nicked: bool,
) -> crate::verification::Operation {
    crate::verification::Operation::Ligate {
        fragments,
        circular,
        allow_nicked,
        options: Options::default(),
    }
}

fn hex(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|b| b.is_ascii_hexdigit())
}
impl Report {
    /// Complete retained top sequences with stable product IDs. Full keyed paths,
    /// opposite strands and modification evidence remain in the typed JSON.
    pub fn fasta(&self) -> Result<String, Error> {
        self.validate_binding(&self.context_revision)?;
        let mut text = String::new();
        for (index, product) in std::iter::once(&self.intended)
            .chain(&self.alternatives)
            .enumerate()
        {
            let id = if index == 0 {
                "intended".into()
            } else {
                format!("alternative_{index}")
            };
            text.push_str(&format!(
                ">{id} topology={} material={}\n",
                if product.topology == Topology::Circular {
                    "circular"
                } else {
                    "linear"
                },
                product.modifications.state_revision
            ));
            for chunk in product.sequence.as_str().as_bytes().chunks(80) {
                text.push_str(std::str::from_utf8(chunk).map_err(|_| bad("non-ASCII native DNA"))?);
                text.push('\n');
            }
            if text.len() > MAX_BYTES {
                return Err(Error::TooLarge);
            }
        }
        Ok(text)
    }
    /// Documentary/context checks only. A consumer must replay the path from
    /// its actual source material before using it as physical state.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        shape(&self.request)?;
        let accounted_paths = self
            .alternatives
            .len()
            .checked_add(1)
            .and_then(|n| n.checked_add(self.incompatible_paths))
            .and_then(|n| n.checked_add(self.unresolved_paths))
            .and_then(|n| n.checked_add(self.withheld_products));
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || self.scientific_engine.is_empty()
            || self.scientific_engine.len() > 128
            || self.material_engine.is_empty()
            || self.material_engine.len() > 128
            || self.context_revision
                != context(
                    &self.sources,
                    &self.request,
                    &self.scientific_engine,
                    &self.material_engine,
                )?
            || self.sources.len() != self.request.intended.len()
            || self.intended.path != self.request.intended
            || self.alternatives.len() + 1 > self.request.options.max_products
            || self.paths_checked > self.request.options.max_path_checks
            || self.work > self.request.options.max_work
            || self.enumeration_complete != self.limits.is_empty()
            || self.exhaustive_within_scope
                != (self.enumeration_complete && self.unresolved_paths == 0)
            || accounted_paths != Some(self.paths_checked)
            || self.unresolved_examples.len() > 8
            || self.unresolved_examples.len() > self.unresolved_paths
            || self
                .unresolved_examples
                .iter()
                .any(|e| e.reason.chars().count() > 256)
            || self.warnings.len() > 8
            || self.warnings.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("report context, scope or bounded counters changed"));
        }
        let skipped = self.sources.len() > self.request.options.max_enumeration_fragments;
        if skipped != self.limits.contains(&Limit::EnumerationFragments)
            || (skipped && (self.paths_checked != 1 || !self.alternatives.is_empty()))
        {
            return Err(bad("enumeration fragment scope changed"));
        }
        let mut sources = BTreeMap::new();
        for s in &self.sources {
            if !hex(&s.state_revision)
                || s.top_length == 0
                || s.bottom_length == 0
                || s.top_length > 1_000_000
                || s.bottom_length > 1_000_000
                || sources.insert(s.key.as_str(), s).is_some()
            {
                return Err(bad("source identity changed"));
            }
        }
        let mut paths = BTreeSet::new();
        let mut bases = 0usize;
        for product in std::iter::once(&self.intended).chain(&self.alternatives) {
            product.duplex.validate()?;
            if product.path.len() != sources.len()
                || product.sequence != product.duplex.top
                || product.topology != product.duplex.topology
                || (product.topology == Topology::Circular) != self.request.circular
                || !hex(&product.modifications.state_revision)
                || product.modifications.top_length != product.duplex.top.len()
                || product.modifications.bottom_length != product.duplex.bottom.len()
                || product.modifications.generation.len() > 1000
                || (self.request.circular && product.path.first() != self.request.intended.first())
            {
                return Err(bad("product path or strand geometry changed"));
            }
            let mut seen = BTreeSet::new();
            let mut top = 0usize;
            let mut bottom = 0usize;
            for p in &product.path {
                let s = sources
                    .get(p.key.as_str())
                    .ok_or_else(|| bad("unknown product source"))?;
                if !seen.insert(p.key.as_str()) {
                    return Err(bad("product duplicates a source"));
                }
                top = top
                    .checked_add(if p.reverse {
                        s.bottom_length
                    } else {
                        s.top_length
                    })
                    .ok_or(Error::TooLarge)?;
                bottom = bottom
                    .checked_add(if p.reverse {
                        s.top_length
                    } else {
                        s.bottom_length
                    })
                    .ok_or(Error::TooLarge)?;
            }
            if top != product.duplex.top.len()
                || bottom != product.duplex.bottom.len()
                || !paths.insert(
                    product
                        .path
                        .iter()
                        .map(|p| (&p.key, p.reverse))
                        .collect::<Vec<_>>(),
                )
            {
                return Err(bad(
                    "product does not conserve one copy of each actual source strand",
                ));
            }
            bases = bases
                .checked_add(
                    top.checked_mul(2)
                        .and_then(|n| n.checked_add(bottom))
                        .ok_or(Error::TooLarge)?,
                )
                .ok_or(Error::TooLarge)?;
        }
        if bases != self.retained_product_bases
            || bases > self.request.options.max_total_product_bases
            || self.limits.len() > 7
            || self
                .limits
                .iter()
                .enumerate()
                .any(|(i, l)| self.limits[..i].contains(l))
        {
            return Err(bad("retained material or limit accounting changed"));
        }
        count(self, self.request.options.max_report_bytes)?;
        Ok(())
    }
}

#[cfg(test)]
mod native_tests {
    use super::*;
    #[test]
    fn private_intended_material_is_retained_independently_of_documentary_product() {
        let a = TrackedDuplex::exact_pcr(&Sequence::dna("AGTC").unwrap())
            .unwrap()
            .treat(crate::chemistry::EndTreatment::Phosphorylate)
            .unwrap();
        let b = TrackedDuplex::exact_pcr(&Sequence::dna("GACT").unwrap())
            .unwrap()
            .treat(crate::chemistry::EndTreatment::Phosphorylate)
            .unwrap();
        let mut options = Options::default();
        options.max_path_checks = 1;
        let request = Request {
            intended: vec![
                Placement {
                    key: "a".into(),
                    reverse: false,
                },
                Placement {
                    key: "b".into(),
                    reverse: false,
                },
            ],
            circular: true,
            allow_nicked: false,
            options,
        };
        let (mut report, actual) = analyze_fresh(
            &[
                Operand {
                    key: "a",
                    material: &a,
                },
                Operand {
                    key: "b",
                    material: &b,
                },
            ],
            &request,
        )
        .unwrap()
        .into_parts();
        assert_eq!(actual.duplex.top.as_str(), "AGTCGACT");
        assert_eq!(
            modified_duplex::summary(&actual).unwrap().state_revision,
            report.intended.modifications.state_revision
        );
        report.intended.sequence = Sequence::dna("AAAAAAAA").unwrap();
        assert_eq!(actual.duplex.top.as_str(), "AGTCGACT");
    }
}
