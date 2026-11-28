//! Source-bound structural checks under explicit annotation-role assumptions.
//! No expression, packaging, titre, annotation-truth or experimental QC verdict.
use crate::sequence::translation;
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};

pub mod causal;
mod cds;
#[cfg(test)]
mod context_orf_tests;
pub mod element_chain;
#[cfg(test)]
mod element_tests;
mod geometry;
pub mod initiation_context;
pub mod required_parts;
#[cfg(test)]
mod tests;
pub mod unannotated_orfs;
pub mod vector_class;

pub const ENGINE_VERSION: &str = "construct-qc-3";
pub const STATEMENT: &str = "Construct QC checks declared structural requirements against the current sequence and annotations. Passing checks do not establish that a construct will express, package, produce titre or behave as intended; failing or indeterminate checks do not establish that it will not.";
const MAX_FEATURES: usize = 256;
const MAX_BASES: usize = 2_000_000;
const MAX_EXTRACTED_BASES: usize = 1_000_000;
pub const MAX_REPORT_BYTES: usize = 2 * 1024 * 1024;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum Verdict {
    Pass,
    Fail,
    Indeterminate,
    Informational,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum Rule {
    AnnotationBasis,
    CdsFrame,
    CdsStart,
    CdsTerminalStop,
    CdsInternalStop,
    StoredTranslation,
    FusionRegister,
    DeclaredJunctionRegister,
    ReferenceIntegrity,
    ElementChain,
    UpstreamElement,
    DownstreamElement,
    RequiredPart,
    InitiationContext,
    UnannotatedOrf,
    Scope,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Citation {
    pub source: String,
    pub locator: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    pub source_url: String,
    pub accession: String,
    pub retrieved_on: String,
    pub license: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct AnnotationAssumption {
    pub feature_index: usize,
    pub source_revision: String,
    pub annotation_revision: String,
    /// Captured declaration of role/coordinates, never proof of annotation truth.
    pub rationale: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ReferenceElement {
    pub id: String,
    pub sequence: Sequence,
    pub sequence_revision: String,
    pub provenance: Provenance,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ReferenceClaim {
    pub feature_index: usize,
    pub reference_id: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct FusionClaim {
    pub cds_feature: usize,
    pub part_feature: usize,
    pub expected_offset_mod3: u8,
    pub rationale: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct JunctionClaim {
    pub cds_feature: usize,
    pub boundary: usize,
    pub expected_offset_mod3: u8,
    pub rationale: String,
    /// Retained identifiers only; not replay authority or a causal edit proof.
    pub history_document_id: Option<String>,
    pub history_event_id: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Options {
    pub schema_version: u16,
    pub genetic_code: u8,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, AnnotationAssumption, 256>"
    )]
    pub annotations: Vec<AnnotationAssumption>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, ReferenceElement, 256>")]
    pub references: Vec<ReferenceElement>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, ReferenceClaim, 256>")]
    pub reference_claims: Vec<ReferenceClaim>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, FusionClaim, 256>")]
    pub fusions: Vec<FusionClaim>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, JunctionClaim, 256>")]
    pub junctions: Vec<JunctionClaim>,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, element_chain::ElementChain, 64>"
    )]
    pub element_chains: Vec<element_chain::ElementChain>,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, required_parts::RequiredPart, 64>"
    )]
    pub required_parts: Vec<required_parts::RequiredPart>,
    pub required_part_limits: required_parts::RequiredPartLimits,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, initiation_context::ContextClaim, 64>"
    )]
    pub initiation_contexts: Vec<initiation_context::ContextClaim>,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub orf_scan: Option<unannotated_orfs::ScanOptions>,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            schema_version: 3,
            genetic_code: 1,
            annotations: vec![],
            references: vec![],
            reference_claims: vec![],
            fusions: vec![],
            junctions: vec![],
            element_chains: vec![],
            required_parts: vec![],
            required_part_limits: Default::default(),
            initiation_contexts: vec![],
            orf_scan: None,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Finding {
    pub rule: Rule,
    pub verdict: Verdict,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 3>")]
    pub feature_indices: Vec<usize>,
    /// Coordinates may come from literal reference matches without feature IDs.
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Location, 16>")]
    pub coordinates: Vec<Location>,
    pub reason: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Citation, 8>")]
    pub citations: Vec<Citation>,
    /// True means a declared feature role is an assumption, not a verified identity.
    pub under_declared_annotation_assumptions: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CodePin {
    pub descriptor: translation::CodeDescriptor,
    /// All 64 codons in TCAG enumeration order, translated by the existing owner.
    pub ordinary_amino_acids: String,
    pub revision: String,
}
impl CodePin {
    /// Capture the current canonical descriptor and all 64 ordinary codons.
    /// This is table reproducibility evidence, not an independent biology model.
    pub fn for_id(id: u8) -> Result<Self, Error> {
        let descriptor = translation::CodeDescriptor::for_id(id)?;
        let mut codons = String::with_capacity(192);
        for a in b"TCAG" {
            for b in b"TCAG" {
                for c in b"TCAG" {
                    codons.extend([*a as char, *b as char, *c as char]);
                }
            }
        }
        let sequence = Sequence::dna(&codons)?;
        let translated = translation::translate(
            &sequence,
            &translation::TranslationOptions {
                genetic_code: id,
                ..Default::default()
            },
        )?;
        let revision = hash(&(
            crate::verification::ENGINE_VERSION,
            &descriptor,
            &translated,
        ))?;
        Ok(Self {
            descriptor,
            ordinary_amino_acids: translated.protein,
            revision,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ReferencePin {
    pub id: String,
    pub sequence_revision: String,
    pub provenance: Provenance,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine_version: String,
    pub context_revision: String,
    pub source_revision: String,
    pub source_length: usize,
    pub source_feature_count: usize,
    pub source_topology: Topology,
    pub annotations_revision: String,
    pub options_revision: String,
    pub references_revision: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, CodePin, 32>")]
    pub codes: Vec<CodePin>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, ReferencePin, 256>")]
    pub references: Vec<ReferencePin>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Finding, 4096>")]
    pub findings: Vec<Finding>,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, element_chain::ElementChainReport, 64>"
    )]
    pub element_chains: Vec<element_chain::ElementChainReport>,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, required_parts::RequiredPartReport, 64>"
    )]
    pub required_parts: Vec<required_parts::RequiredPartReport>,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, required_parts::ReferenceSearch, 64>"
    )]
    pub required_part_searches: Vec<required_parts::ReferenceSearch>,
    #[serde(
        deserialize_with = "collagene_model::decode::bounded_vec::<_, initiation_context::ContextEvidence, 64>"
    )]
    pub initiation_contexts: Vec<initiation_context::ContextEvidence>,
    #[serde(deserialize_with = "crate::encoding::required_option")]
    #[cfg_attr(feature = "json-schema", schemars(required))]
    pub orf_scan: Option<unannotated_orfs::ScanEvidence>,
    pub required_part_search_work: usize,
    pub required_part_search_engine: String,
    pub statement: String,
    pub report_revision: String,
}

fn bad(s: &str) -> Error {
    Error::Unsupported(format!("construct QC: {s}"))
}
fn text(s: &str, n: usize) -> bool {
    !s.trim().is_empty() && s.len() <= n && !s.chars().any(char::is_control)
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
struct Budget(usize);
impl Write for Budget {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_sub(b.len())
            .ok_or_else(|| std::io::Error::other("construct QC bytes"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded(v: &impl Serialize, n: usize) -> Result<(), Error> {
    serde_json::to_writer(Budget(n), v).map_err(|_| Error::TooLarge)
}
fn hash(v: &impl Serialize) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(v)
}
pub fn annotation_revision(feature: &Feature) -> Result<String, Error> {
    collagene_model::bounds::feature_preflight(&[std::slice::from_ref(feature)])?;
    bounded(feature, 256 * 1024)?;
    hash(feature)
}
pub fn reference_revision(sequence: &Sequence) -> Result<String, Error> {
    if sequence.len() > MAX_EXTRACTED_BASES {
        return Err(Error::TooLarge);
    }
    hash(sequence)
}
fn location_bases(l: &Location) -> Result<usize, Error> {
    match l {
        Location::Range { start, end, .. } => end
            .checked_sub(*start)
            .filter(|n| *n <= MAX_EXTRACTED_BASES)
            .ok_or(Error::TooLarge),
        Location::Complement(inner) => location_bases(inner),
        Location::Join(parts) | Location::Order(parts) => parts.iter().try_fold(0usize, |n, p| {
            n.checked_add(location_bases(p)?)
                .filter(|n| *n <= MAX_EXTRACTED_BASES)
                .ok_or(Error::TooLarge)
        }),
        Location::Unresolved(_) => Ok(0),
    }
}
/// Canonical contiguous interval, possibly wrapping once and complemented.
/// Reject richer restored evidence before recursive coordinate validation.
fn interval_evidence_shape(location: &Location) -> bool {
    let inner = match location {
        Location::Complement(inner) => inner.as_ref(),
        location => location,
    };
    let exact_range = |location: &Location| {
        matches!(
            location,
            Location::Range {
                before: false,
                after: false,
                ..
            }
        )
    };
    match inner {
        Location::Join(parts) => parts.len() == 2 && parts.iter().all(exact_range),
        location => exact_range(location),
    }
}
fn add_bases(total: &mut usize, n: usize) -> Result<(), Error> {
    *total = total
        .checked_add(n)
        .filter(|n| *n <= MAX_EXTRACTED_BASES)
        .ok_or(Error::TooLarge)?;
    Ok(())
}
impl Options {
    pub fn validate(&self, m: &Molecule) -> Result<(), Error> {
        if m.sequence.len() > MAX_BASES
            || m.features.len() > MAX_FEATURES
            || [
                self.annotations.len(),
                self.references.len(),
                self.reference_claims.len(),
                self.fusions.len(),
                self.junctions.len(),
            ]
            .iter()
            .any(|n| *n > 256)
        {
            return Err(Error::TooLarge);
        }
        if self.element_chains.len() > 64
            || self.required_parts.len() > 64
            || self.initiation_contexts.len() > 64
        {
            return Err(Error::TooLarge);
        }
        collagene_model::bounds::model_preflight(&[m])?;
        // Streaming byte admission precedes qualifier copies, hashes and extraction.
        bounded(m, 8 * 1024 * 1024)?;
        bounded(self, 2 * 1024 * 1024)?;
        m.validate()?;
        if !matches!(m.sequence.alphabet(), Alphabet::Dna | Alphabet::Rna)
            || self.schema_version != 3
        {
            return Err(bad("current nucleotide source/options required"));
        }
        translation::CodeDescriptor::for_id(self.genetic_code)?;
        self.required_part_limits.validate()?;
        if let Some(scan) = &self.orf_scan {
            scan.validate()?;
        }
        for claim in &self.initiation_contexts {
            claim.validate(m)?;
        }
        let mut nodes = 0usize;
        let mut extracted = 0usize;
        for f in &m.features {
            let mut stack = vec![&f.location];
            while let Some(l) = stack.pop() {
                nodes += 1;
                if nodes > 4096 {
                    return Err(Error::TooLarge);
                }
                match l {
                    Location::Complement(i) => stack.push(i),
                    Location::Join(p) | Location::Order(p) => {
                        if p.len() > 4096usize.saturating_sub(nodes + stack.len()) {
                            return Err(Error::TooLarge);
                        }
                        stack.extend(p)
                    }
                    _ => (),
                }
            }
            if f.kind.eq_ignore_ascii_case("CDS") {
                add_bases(&mut extracted, location_bases(&f.location)?)?;
            }
        }
        let mut reference_bases = 0;
        let mut ids = BTreeSet::new();
        for r in &self.references {
            add_bases(&mut reference_bases, r.sequence.len())?;
            if !text(&r.id, 128)
                || !ids.insert(&r.id)
                || !hex(&r.sequence_revision)
                || !matches!(r.sequence.alphabet(), Alphabet::Dna | Alphabet::Rna)
                || !r.provenance.source_url.starts_with("https://")
                || !text(&r.provenance.source_url, 2048)
                || !text(&r.provenance.accession, 256)
                || !text(&r.provenance.retrieved_on, 32)
                || !text(&r.provenance.license, 256)
            {
                return Err(bad("invalid exact reference pin/provenance"));
            }
        }
        let mut assumed = BTreeSet::new();
        for a in &self.annotations {
            if a.feature_index >= m.features.len()
                || !assumed.insert(a.feature_index)
                || !hex(&a.source_revision)
                || !hex(&a.annotation_revision)
                || !text(&a.rationale, 1024)
            {
                return Err(bad("stale/invalid annotation-role assumption"));
            }
        }
        for r in &self.reference_claims {
            if r.feature_index >= m.features.len() || !text(&r.reference_id, 128) {
                return Err(bad(
                    "reference claim must identify existing feature/pinned reference",
                ));
            }
            add_bases(
                &mut extracted,
                location_bases(&m.features[r.feature_index].location)?,
            )?;
        }
        for f in &self.fusions {
            if f.cds_feature >= m.features.len()
                || f.part_feature >= m.features.len()
                || f.expected_offset_mod3 > 2
                || !text(&f.rationale, 1024)
            {
                return Err(bad("invalid fusion declaration"));
            }
        }
        for j in &self.junctions {
            if j.cds_feature >= m.features.len()
                || j.boundary > m.sequence.len()
                || j.expected_offset_mod3 > 2
                || !text(&j.rationale, 1024)
                || j.history_document_id
                    .as_ref()
                    .is_some_and(|s| !text(s, 128))
                || j.history_event_id.as_ref().is_some_and(|s| !text(s, 128))
                || j.history_document_id.is_some() != j.history_event_id.is_some()
            {
                return Err(bad("invalid declared junction/history identifiers"));
            }
        }
        let mut chains = BTreeSet::new();
        for chain in &self.element_chains {
            chain.validate(m)?;
            if !chains.insert(&chain.id) {
                return Err(bad("duplicate element-chain ID"));
            }
            for e in [&chain.upstream, &chain.downstream].into_iter().flatten() {
                add_bases(
                    &mut extracted,
                    location_bases(&m.features[e.feature_index].location)?,
                )?;
            }
        }
        let mut requirements = BTreeSet::new();
        for part in &self.required_parts {
            part.validate()?;
            if !requirements.insert(&part.id) {
                return Err(bad("duplicate required-part ID"));
            }
        }
        Ok(())
    }
}

fn code_pins(m: &Molecule, o: &Options) -> Result<Vec<CodePin>, Error> {
    let mut codes = BTreeSet::from([o.genetic_code]);
    for f in &m.features {
        for (k, v) in &f.qualifiers {
            if k == "transl_table" {
                if let Some(id) = v.as_ref().and_then(|s| s.trim().parse::<u8>().ok()) {
                    if translation::CodeDescriptor::for_id(id).is_ok() {
                        codes.insert(id);
                    }
                }
            }
        }
    }
    codes
        .into_iter()
        .map(CodePin::for_id)
        .collect::<Result<Vec<_>, _>>()
}
pub fn revision(m: &Molecule, o: &Options) -> Result<String, Error> {
    o.validate(m)?;
    hash(&(
        ENGINE_VERSION,
        crate::verification::ENGINE_VERSION,
        crate::sequence::search::SEARCH_ENGINE_VERSION,
        m,
        o,
        code_pins(m, o)?,
    ))
}
struct Checks<'a> {
    m: &'a Molecule,
    o: &'a Options,
    source_revision: String,
    context_revision: String,
    findings: Vec<Finding>,
    budget: Budget,
    frames: BTreeMap<usize, cds::CodingFrame>,
    assumptions: BTreeSet<usize>,
    reference_facts: BTreeMap<String, ReferenceFacts>,
}
struct ReferenceFacts {
    pin_matches: bool,
    concrete: bool,
}
impl<'a> Checks<'a> {
    fn new(m: &'a Molecule, o: &'a Options) -> Result<Self, Error> {
        let context_revision = revision(m, o)?;
        let source_revision = m.revision();
        let assumptions = o
            .annotations
            .iter()
            .filter(|a| {
                a.source_revision == source_revision
                    && !m.features[a.feature_index].needs_review
                    && annotation_revision(&m.features[a.feature_index])
                        .is_ok_and(|r| r == a.annotation_revision)
            })
            .map(|a| a.feature_index)
            .collect();
        // Each captured reference is hashed/classified once per job. Reusing one
        // long reference across many tiny feature claims must not amplify that work.
        let reference_facts = o
            .references
            .iter()
            .map(|r| {
                Ok((
                    r.id.clone(),
                    ReferenceFacts {
                        pin_matches: r.sequence_revision == reference_revision(&r.sequence)?,
                        concrete: r
                            .sequence
                            .as_str()
                            .bytes()
                            .all(|b| crate::sequence::analysis::base_mask(b).is_power_of_two()),
                    },
                ))
            })
            .collect::<Result<BTreeMap<_, _>, Error>>()?;
        Ok(Self {
            m,
            o,
            source_revision,
            context_revision,
            findings: vec![],
            budget: Budget(MAX_REPORT_BYTES / 2),
            frames: BTreeMap::new(),
            assumptions,
            reference_facts,
        })
    }
}
impl Checks<'_> {
    fn assumed(&self, index: usize) -> bool {
        self.assumptions.contains(&index)
    }
    fn push(
        &mut self,
        rule: Rule,
        verdict: Verdict,
        indices: &[usize],
        reason: &str,
    ) -> Result<(), Error> {
        if self.findings.len() >= 4096 {
            return Err(Error::TooLarge);
        }
        let locations = indices
            .iter()
            .map(|&i| &self.m.features[i].location)
            .collect::<Vec<_>>();
        let citations = vec![
            Citation {
                source: "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi".into(),
                locator:
                    "Declared genetic-code codon/initiator semantics; not an expression prediction"
                        .into(),
            },
            Citation {
                source: format!("collagene:{ENGINE_VERSION}"),
                locator: format!("{rule:?}: exact declared structural predicate"),
            },
        ];
        serde_json::to_writer(
            &mut self.budget,
            &(&rule, &verdict, indices, &locations, reason, &citations),
        )
        .map_err(|_| Error::TooLarge)?;
        let under_declared_annotation_assumptions =
            !indices.is_empty() && indices.iter().all(|&i| self.assumed(i));
        self.findings.push(Finding {
            rule,
            verdict,
            feature_indices: indices.to_vec(),
            coordinates: locations.into_iter().cloned().collect(),
            reason: reason.into(),
            citations,
            under_declared_annotation_assumptions,
        });
        Ok(())
    }
    fn push_literal_finding(
        &mut self,
        rule: Rule,
        verdict: Verdict,
        coordinates: &[Location],
        reason: &str,
    ) -> Result<(), Error> {
        if coordinates.len() > 16 || self.findings.len() >= 4096 {
            return Err(Error::TooLarge);
        }
        let citations = vec![Citation {
            source: format!("collagene:{ENGINE_VERSION}"),
            locator: format!(
                "{rule:?}: bounded current-source structural evidence under declared options; not functional identity"
            ),
        }];
        serde_json::to_writer(
            &mut self.budget,
            &(&rule, &verdict, coordinates, reason, &citations),
        )
        .map_err(|_| Error::TooLarge)?;
        self.findings.push(Finding {
            rule,
            verdict,
            feature_indices: vec![],
            coordinates: coordinates.to_vec(),
            reason: reason.into(),
            citations,
            under_declared_annotation_assumptions: false,
        });
        Ok(())
    }
}
pub fn assess(m: &Molecule, o: &Options) -> Result<Report, Error> {
    let mut checks = Checks::new(m, o)?;
    let context_revision = checks.context_revision.clone();
    let source_revision = checks.source_revision.clone();
    cds::inspect(&mut checks)?;
    geometry::inspect(&mut checks)?;
    let element_chains = element_chain::inspect(&mut checks)?;
    let required_parts = required_parts::inspect(&mut checks)?;
    let initiation_contexts = initiation_context::inspect(&mut checks)?;
    let orf_scan = unannotated_orfs::inspect(&mut checks)?;
    checks.push(
        Rule::Scope,
        Verdict::Informational,
        &[],
        "Checks cover declared CDS, reference identity, explicit fusion/junction register, declared flanking-element geometry, required literal-reference intervals, requested declared initiation-pattern conformity and requested bounded ORF candidates. History-derived causal junction checks are not evaluated. No declaration, motif match or ORF candidate verifies biological role/function or expression.",
    )?;
    let codes = code_pins(m, o)?;
    let references = o
        .references
        .iter()
        .map(|r| ReferencePin {
            id: r.id.clone(),
            sequence_revision: r.sequence_revision.clone(),
            provenance: r.provenance.clone(),
        })
        .collect();
    let mut report = Report {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision,
        source_length: m.sequence.len(),
        source_feature_count: m.features.len(),
        source_topology: m.topology,
        annotations_revision: hash(&m.features)?,
        options_revision: hash(o)?,
        references_revision: hash(&o.references)?,
        codes,
        references,
        findings: checks.findings,
        element_chains,
        required_parts: required_parts.reports,
        required_part_searches: required_parts.searches,
        initiation_contexts,
        orf_scan,
        required_part_search_work: required_parts.work,
        required_part_search_engine: crate::sequence::search::SEARCH_ENGINE_VERSION.into(),
        statement: STATEMENT.into(),
        report_revision: String::new(),
    };
    bounded(&report, MAX_REPORT_BYTES)?;
    report.report_revision = report.content_revision()?;
    report.validate()?;
    Ok(report)
}
impl Report {
    fn content_revision(&self) -> Result<String, Error> {
        hash(&(
            (
                &self.engine_version,
                &self.context_revision,
                &self.source_revision,
                self.source_length,
                self.source_feature_count,
                self.source_topology,
            ),
            (
                &self.annotations_revision,
                &self.options_revision,
                &self.references_revision,
                &self.codes,
                &self.references,
            ),
            (
                &self.element_chains,
                &self.required_parts,
                &self.required_part_searches,
                self.required_part_search_work,
                &self.required_part_search_engine,
            ),
            (
                &self.findings,
                &self.statement,
                &self.initiation_contexts,
                &self.orf_scan,
            ),
        ))
    }
    /// Shape/content integrity only; a restored report has no current authority.
    pub fn validate(&self) -> Result<(), Error> {
        bounded(self, MAX_REPORT_BYTES)?;
        if self.engine_version != ENGINE_VERSION
            || self.statement != STATEMENT
            || self.findings.len() > 4096
            || self.references.len() > 256
            || self.codes.len() > 32
            || self.codes.is_empty()
            || self.element_chains.len() > 64
            || self.required_parts.len() > 64
            || self.required_part_searches.len() > 64
            || self.required_part_search_work > 50_000_000
            || self.required_part_search_engine != crate::sequence::search::SEARCH_ENGINE_VERSION
            || !(1..=MAX_BASES).contains(&self.source_length)
            || self.source_feature_count > MAX_FEATURES
            || [
                &self.context_revision,
                &self.source_revision,
                &self.annotations_revision,
                &self.options_revision,
                &self.references_revision,
            ]
            .iter()
            .any(|s| !hex(s))
            || self.report_revision != self.content_revision()?
        {
            return Err(bad("report version/shape/content differs"));
        }
        for f in &self.findings {
            if f.feature_indices.len() > 3
                || f.coordinates.len() > 16
                || !text(&f.reason, 4096)
                || f.citations.is_empty()
                || f.citations.len() > 8
                || f.feature_indices
                    .iter()
                    .any(|&i| i >= self.source_feature_count)
                || (f.under_declared_annotation_assumptions && f.feature_indices.is_empty())
            {
                return Err(bad("invalid finding coordinates/citations"));
            }
            for location in &f.coordinates {
                location.validate(self.source_length)?;
            }
        }
        initiation_context::validate_reports(self)?;
        unannotated_orfs::validate_reports(self)?;
        required_parts::validate_reports(self)?;
        element_chain::validate_reports(self)?;
        let mut ids = BTreeSet::new();
        for c in &self.codes {
            if c.descriptor.id == 0
                || !ids.insert(c.descriptor.id)
                || c.ordinary_amino_acids.len() != 64
                || !c
                    .ordinary_amino_acids
                    .bytes()
                    .all(|b| b"ACDEFGHIKLMNPQRSTVWYX*".contains(&b))
                || !hex(&c.revision)
            {
                return Err(bad("invalid captured genetic-code pin"));
            }
        }
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad("report belongs to different source/options"));
        }
        Ok(())
    }
    pub fn validate_for(&self, m: &Molecule, o: &Options) -> Result<(), Error> {
        self.validate_binding(&revision(m, o)?)?;
        if *self != assess(m, o)? {
            return Err(bad("report differs from freshly computed structural facts"));
        }
        Ok(())
    }
}
