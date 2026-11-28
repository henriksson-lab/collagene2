//! Reference-only adapter for the independently bounded reciprocal-exchange engine.
//! Exact sequence maps support annotation projection, never inherited physical state.
use crate::sequence::annotations::{self, MappingSegment};
use collagene_model::{Error, Feature, Location, Molecule, Topology};
pub use collagene_recombination as reference;
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const ENGINE_VERSION: &str = "recombination-reference-adapter-2";
pub const RECIPE_KEY: &str = "recombination:unverified_recipe";
pub const REFERENCE_ONLY_KEY: &str = "recombination:reference_only";
pub const CHEMISTRY_UNMAPPED_KEY: &str = "source:legacy:chemistry_unmapped";
macro_rules! dto{($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
dto! {pub struct ProfileSnapshot{pub profile:reference::Profile,pub revision:String}}
impl ProfileSnapshot {
    pub fn new(profile: reference::Profile) -> Result<Self, Error> {
        let revision = profile.revision()?;
        Ok(Self { profile, revision })
    }
    pub fn manufacturer_context() -> Result<Self, Error> {
        Self::new(reference::manufacturer_context_profile()?)
    }
    pub fn validate(&self) -> Result<(), Error> {
        if !hex(&self.revision) || self.profile.revision()? != self.revision {
            return Err(bad("profile snapshot content/revision mismatch"));
        }
        Ok(())
    }
}
dto! {pub struct AnnotationLimits{pub max_work:usize,pub max_copy_bytes:usize,pub max_features:usize,pub max_report_bytes:usize}}
impl Default for AnnotationLimits {
    fn default() -> Self {
        Self {
            max_work: 4_000_000,
            max_copy_bytes: 16 * 1024 * 1024,
            max_features: 16384,
            max_report_bytes: 32 * 1024 * 1024,
        }
    }
}
impl AnnotationLimits {
    pub fn validate(&self) -> Result<(), Error> {
        if self.max_work == 0
            || self.max_work > 16_000_000
            || self.max_copy_bytes == 0
            || self.max_copy_bytes > 64 * 1024 * 1024
            || self.max_features == 0
            || self.max_features > 32768
            || self.max_report_bytes < 4096
            || self.max_report_bytes > 32 * 1024 * 1024
        {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
}
dto! {pub struct Request{pub profile:ProfileSnapshot,pub reaction:reference::Reaction,pub limits:reference::Limits,pub annotation_limits:AnnotationLimits}}
impl Request {
    pub fn new(profile: ProfileSnapshot, reaction: reference::Reaction) -> Self {
        Self {
            profile,
            reaction,
            limits: reference::Limits::default(),
            annotation_limits: AnnotationLimits::default(),
        }
    }
    pub fn validate(&self) -> Result<(), Error> {
        self.annotation_limits.validate()?;
        self.profile.validate()
    }
}
dto! {pub struct Selection{pub batch_revision:String,pub candidate_id:usize,pub component_id:usize}}
dto! {pub struct ComponentIdentity{pub component_id:usize,pub reference_revision:String,pub length:usize,pub topology:Topology}}
dto! {pub struct SourceIdentity{pub role:String,pub revision:String}}
dto! {pub struct SourceAnnotationEvidence{pub input:usize,pub source:SourceIdentity,pub segments:Vec<MappingSegment>,pub co_sourced_segments:Vec<MappingSegment>,pub mapped_features:usize,pub reviewed_features:usize,pub dropped_features:usize,pub unmapped_features:Vec<Feature>,pub warnings:Vec<String>}}
dto! {pub struct UnverifiedRecipe{pub schema_version:u16,pub engine_version:String,pub request:Request,pub selection:Selection,pub source_identities:Vec<SourceIdentity>,pub reference_batch_revision:String,pub reference_product_revision:String,pub reciprocal_components:Vec<ComponentIdentity>}}
dto! {pub struct ProductDocument{pub engine_version:String,pub selection:Selection,pub reference_batch_revision:String,pub reference_product_revision:String,pub reciprocal_components:Vec<ComponentIdentity>,pub molecule:Molecule,pub annotation_evidence:Vec<SourceAnnotationEvidence>,pub warnings:Vec<String>}}
dto! {pub struct RecombinationReport{pub engine_version:String,pub batch_revision:String,pub reference:reference::Report,pub documents:Vec<ProductDocument>,pub annotation_work:usize,pub annotation_copy_bytes:usize,pub annotation_features:usize,pub warnings:Vec<String>}}
dto! {pub struct TargetAssessment{pub engine_version:String,pub batch_revision:String,pub reference:reference::TargetCheck}}
dto! {pub struct CheckedProduct{pub engine_version:String,pub batch_revision:String,pub product:ProductDocument,pub check:TargetAssessment}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("reference recombination: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|b| b.is_ascii_hexdigit())
}
struct Bound {
    used: usize,
    max: usize,
}
impl Write for Bound {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.used = self
            .used
            .checked_add(b.len())
            .filter(|n| *n <= self.max)
            .ok_or_else(|| std::io::Error::other("bounded recombination evidence"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn size<T: Serialize + ?Sized>(v: &T, max: usize) -> Result<usize, Error> {
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
pub fn source_preflight(molecules: [&Molecule; 2], request: &Request) -> Result<(), Error> {
    request.validate()?;
    reference::source_preflight(molecules, &request.profile.profile, &request.limits, None)
}
pub fn source_target_preflight(
    molecules: [&Molecule; 2],
    request: &Request,
    target: &Molecule,
) -> Result<(), Error> {
    request.validate()?;
    reference::source_preflight(
        molecules,
        &request.profile.profile,
        &request.limits,
        Some(target),
    )
}
pub fn preflight(inputs: &[reference::Input], request: &Request) -> Result<(), Error> {
    request.validate()?;
    reference::preflight(inputs, &request.profile.profile, &request.limits)
}
pub fn batch_revision(inputs: &[reference::Input], request: &Request) -> Result<String, Error> {
    preflight(inputs, request)?;
    let revision = reference::batch_revision(
        inputs,
        &request.profile.profile,
        request.reaction,
        &request.limits,
    )?;
    hash(&(ENGINE_VERSION, revision, request))
}
#[derive(Default)]
pub(crate) struct Work {
    pub(crate) units: usize,
    pub(crate) bytes: usize,
    pub(crate) features: usize,
}
fn nodes(location: &Location) -> usize {
    match location {
        Location::Range { .. } | Location::Unresolved(_) => 1,
        Location::Complement(p) => 1 + nodes(p),
        Location::Join(p) | Location::Order(p) => 1 + p.iter().map(nodes).sum::<usize>(),
    }
}
pub(crate) fn charge(
    source: &Molecule,
    maps: usize,
    limits: &AnnotationLimits,
    work: &mut Work,
) -> Result<(), Error> {
    // The engine preflight bounded location depth/nodes/strings before this recursion.
    let n = source
        .features
        .iter()
        .map(|f| nodes(&f.location) + f.qualifiers.len() + 1)
        .sum::<usize>();
    let units = n.checked_mul((maps + 1) * 32).ok_or(Error::TooLarge)?;
    let bytes = size(&source.features, 16 * 1024 * 1024)?
        .checked_mul(2)
        .and_then(|v| v.checked_add(n.checked_mul((maps + 1) * 128)?))
        .and_then(|v| v.checked_add(source.features.len() * 1024))
        .ok_or(Error::TooLarge)?;
    work.units = work
        .units
        .checked_add(units)
        .filter(|n| *n <= limits.max_work)
        .ok_or(Error::TooLarge)?;
    work.bytes = work
        .bytes
        .checked_add(bytes)
        .filter(|n| *n <= limits.max_copy_bytes)
        .ok_or(Error::TooLarge)?;
    work.features = work
        .features
        .checked_add(source.features.len())
        .filter(|n| *n <= limits.max_features)
        .ok_or(Error::TooLarge)?;
    Ok(())
}
fn destination(m: &MappingSegment, position: usize) -> Option<usize> {
    (m.source_start <= position && position < m.source_end).then(|| {
        m.destination_start
            + if m.reverse {
                m.source_end - 1 - position
            } else {
                position - m.source_start
            }
    })
}
fn normalize(mut maps: Vec<MappingSegment>) -> Result<Vec<MappingSegment>, Error> {
    maps.sort_by_key(|m| m.source_start);
    let mut out: Vec<MappingSegment> = vec![];
    for m in maps {
        if let Some(p) = out.last_mut() {
            if p.source_end > m.source_start {
                return Err(bad("overlapping source annotation maps"));
            }
            if p.source_end == m.source_start
                && p.reverse == m.reverse
                && if p.reverse {
                    m.destination_start + (m.source_end - m.source_start) == p.destination_start
                } else {
                    p.destination_start + (p.source_end - p.source_start) == m.destination_start
                }
            {
                p.source_end = m.source_end;
                p.destination_start = p.destination_start.min(m.destination_start);
                continue;
            }
        }
        out.push(m)
    }
    Ok(out)
}
pub(crate) fn annotation_maps(
    product: &reference::Product,
    input: usize,
) -> Result<(Vec<MappingSegment>, Vec<MappingSegment>), Error> {
    let mut maps = product
        .segments
        .iter()
        .filter(|m| m.input == input)
        .map(|m| MappingSegment {
            source_start: m.source_start,
            source_end: m.source_end,
            destination_start: m.output_start,
            reverse: m.reverse,
        })
        .collect::<Vec<_>>();
    let mut co = vec![];
    for junction in &product.junctions {
        for parent in junction.parents.iter().filter(|p| p.input == input) {
            for (&source, &dest) in parent
                .overlap_positions
                .iter()
                .zip(&junction.product_overlap_positions)
            {
                if let Some(existing) = maps.iter().find_map(|m| destination(m, source)) {
                    if existing != dest {
                        return Err(bad(
                            "co-sourced overlap maps one source position to distinct product positions",
                        ));
                    }
                    continue;
                }
                let map = MappingSegment {
                    source_start: source,
                    source_end: source + 1,
                    destination_start: dest,
                    reverse: parent.reverse != junction.reverse,
                };
                maps.push(map);
                co.push(map);
            }
        }
    }
    Ok((normalize(maps)?, normalize(co)?))
}
pub(crate) fn mark_co_source(features: &mut [Feature], co: &[MappingSegment]) {
    for feature in features {
        if feature.location.spans().iter().any(|(a, b, _)| {
            co.iter().any(|m| {
                *a < m.destination_start + (m.source_end - m.source_start)
                    && m.destination_start < *b
            })
        }) {
            feature.needs_review = true;
            feature.qualifiers.retain(|(k, _)| {
                !matches!(k.as_str(), "translation" | "codon_start" | "transl_table")
            });
            feature.qualifiers.push(("note".into(),Some("Collagene: annotation includes co-sourced crossover overlap; this maps reference bases, not physical-strand ancestry or function.".into())));
        }
    }
}
fn documents(
    inputs: &[reference::Input],
    request: &Request,
    report: &reference::Report,
    batch: &str,
) -> Result<(Vec<ProductDocument>, Work), Error> {
    let identities = inputs
        .iter()
        .map(|i| SourceIdentity {
            role: i.role.clone(),
            revision: i.molecule.revision(),
        })
        .collect::<Vec<_>>();
    let mut out = vec![];
    let mut work = Work::default();
    for candidate in &report.candidates {
        let components = candidate
            .products
            .iter()
            .enumerate()
            .map(|(component_id, p)| ComponentIdentity {
                component_id,
                reference_revision: p.molecule.revision(),
                length: p.molecule.sequence.len(),
                topology: p.molecule.topology,
            })
            .collect::<Vec<_>>();
        for (component_id, product) in candidate.products.iter().enumerate() {
            let selection = Selection {
                batch_revision: batch.into(),
                candidate_id: candidate.id,
                component_id,
            };
            let mut molecule = product.molecule.clone();
            let mut evidence = vec![];
            let mut warnings=vec!["Reference-only recombination component; duplex chemistry, nicks, methylation, reaction efficacy and selection are not established.".into()];
            for (input, source) in inputs.iter().enumerate() {
                let mapping = annotation_maps(product, input);
                let map_count = mapping
                    .as_ref()
                    .map_or(product.segments.len() + 14, |(m, _)| m.len());
                charge(
                    &source.molecule,
                    map_count,
                    &request.annotation_limits,
                    &mut work,
                )?;
                let (segments, co_sourced_segments, projection) = match mapping {
                    Ok((segments, co)) => {
                        let projection = annotations::project_features_with_topology(
                            &source.molecule,
                            &segments,
                            molecule.sequence.len(),
                            molecule.topology,
                        );
                        (segments, co, projection)
                    }
                    Err(error) => (vec![], vec![], Err(error)),
                };
                let (
                    mapped_features,
                    reviewed_features,
                    dropped_features,
                    unmapped_features,
                    findings,
                ) = match projection {
                    Ok(mut projection) => {
                        mark_co_source(&mut projection.features, &co_sourced_segments);
                        let mapped = projection.features.len();
                        let reviewed = projection
                            .features
                            .iter()
                            .filter(|f| f.needs_review)
                            .count();
                        molecule.features.extend(projection.features);
                        (
                            mapped,
                            reviewed,
                            projection.dropped_features,
                            projection.unmapped_features,
                            projection.warnings,
                        )
                    }
                    Err(error) => {
                        let mut original = source.molecule.features.clone();
                        for f in &mut original {
                            f.needs_review = true;
                        }
                        let findings = vec![format!(
                            "Source '{}' annotations retained as reviewed unmapped evidence; projection was unsupported: {error}",
                            source.role
                        )];
                        (0, 0, 0, original, findings)
                    }
                };
                evidence.push(SourceAnnotationEvidence {
                    input,
                    source: identities[input].clone(),
                    segments,
                    co_sourced_segments,
                    mapped_features,
                    reviewed_features,
                    dropped_features,
                    unmapped_features,
                    warnings: findings,
                });
                if crate::chemistry::has_unmapped_physical_evidence(&source.molecule)
                    || source
                        .molecule
                        .metadata
                        .contains_key("tracked_modifications")
                {
                    warnings.push(format!("Source '{}' physical/modification evidence was not propagated; preserve its immutable source document.",source.role));
                }
            }
            molecule
                .metadata
                .insert(REFERENCE_ONLY_KEY.into(), "true".into());
            molecule
                .metadata
                .insert(CHEMISTRY_UNMAPPED_KEY.into(), "true".into());
            let reference_product_revision = product.molecule.revision();
            let producer_bytes =
                size(&(request, &selection, &identities, &components), 256 * 1024)?;
            work.bytes = work
                .bytes
                .checked_add(producer_bytes * 2)
                .filter(|n| *n <= request.annotation_limits.max_copy_bytes)
                .ok_or(Error::TooLarge)?;
            let recipe = UnverifiedRecipe {
                schema_version: 1,
                engine_version: ENGINE_VERSION.into(),
                request: request.clone(),
                selection: selection.clone(),
                source_identities: identities.clone(),
                reference_batch_revision: report.batch_revision.clone(),
                reference_product_revision: reference_product_revision.clone(),
                reciprocal_components: components.clone(),
            };
            molecule.metadata.insert(
                RECIPE_KEY.into(),
                serde_json::to_string(&recipe).map_err(|_| Error::TooLarge)?,
            );
            molecule.validate()?;
            out.push(ProductDocument {
                engine_version: ENGINE_VERSION.into(),
                selection,
                reference_batch_revision: report.batch_revision.clone(),
                reference_product_revision,
                reciprocal_components: components.clone(),
                molecule,
                annotation_evidence: evidence,
                warnings,
            });
        }
    }
    Ok((out, work))
}
pub fn simulate(
    inputs: &[reference::Input],
    request: &Request,
) -> Result<RecombinationReport, Error> {
    let batch_revision = batch_revision(inputs, request)?;
    let reference = reference::recombine(
        inputs,
        &request.profile.profile,
        request.reaction,
        &request.limits,
    )?;
    let (documents, work) = documents(inputs, request, &reference, &batch_revision)?;
    let report=RecombinationReport{engine_version:ENGINE_VERSION.into(),batch_revision,reference,documents,annotation_work:work.units,annotation_copy_bytes:work.bytes,annotation_features:work.features,warnings:vec!["All reciprocal components and source/profile uncertainty remain in the reference report. Producer recipes are descriptive, never imported reaction/chemistry authority. Unsupported annotation projection requires the retained immutable sources.".into()]};
    report.validate()?;
    size(&report, request.annotation_limits.max_report_bytes)?;
    Ok(report)
}
pub fn product(
    inputs: &[reference::Input],
    request: &Request,
    selection: &Selection,
) -> Result<ProductDocument, Error> {
    if !hex(&selection.batch_revision)
        || selection.batch_revision != batch_revision(inputs, request)?
    {
        return Err(bad(
            "selection belongs to a different ordered source/profile/options request",
        ));
    }
    simulate(inputs, request)?
        .documents
        .into_iter()
        .find(|p| &p.selection == selection)
        .ok_or_else(|| bad("selected candidate/component is not a reconstructed reference product"))
}
pub fn check_revision(
    inputs: &[reference::Input],
    request: &Request,
    target: &Molecule,
    policy: &reference::TargetPolicy,
) -> Result<String, Error> {
    if inputs.len() != 2 {
        return Err(Error::TooLarge);
    }
    source_target_preflight([&inputs[0].molecule, &inputs[1].molecule], request, target)?;
    let reference = reference::check_batch_revision(
        inputs,
        &request.profile.profile,
        request.reaction,
        &request.limits,
        target,
        *policy,
    )?;
    hash(&(ENGINE_VERSION, batch_revision(inputs, request)?, reference))
}
pub fn check_target(
    inputs: &[reference::Input],
    request: &Request,
    target: &Molecule,
    policy: &reference::TargetPolicy,
) -> Result<TargetAssessment, Error> {
    let batch_revision = check_revision(inputs, request, target, policy)?;
    let reference = reference::check_target(
        inputs,
        &request.profile.profile,
        request.reaction,
        &request.limits,
        target,
        *policy,
    )?;
    Ok(TargetAssessment {
        engine_version: ENGINE_VERSION.into(),
        batch_revision,
        reference,
    })
}
pub fn checked_product(
    inputs: &[reference::Input],
    request: &Request,
    selection: &Selection,
    target: &Molecule,
    policy: &reference::TargetPolicy,
) -> Result<CheckedProduct, Error> {
    let batch_revision = checked_product_revision(inputs, request, selection, target, policy)?;
    let check = check_target(inputs, request, target, policy)?;
    if !check.reference.matches.iter().any(|m| {
        m.candidate == selection.candidate_id
            && m.product == selection.component_id
            && m.fully_specified_dna
    }) {
        return Err(bad(
            "independent target does not match the specifically selected fully specified reference component",
        ));
    }
    let product = product(inputs, request, selection)?;
    Ok(CheckedProduct {
        engine_version: ENGINE_VERSION.into(),
        batch_revision,
        product,
        check,
    })
}

pub fn borrowed_preflight(
    sources: [(&str, &Molecule); 2],
    request: &Request,
    target: Option<&Molecule>,
) -> Result<(), Error> {
    if sources
        .iter()
        .any(|(role, _)| role.len() > 256 || role.trim().is_empty())
        || sources[0].0 == sources[1].0
    {
        return Err(bad("two distinct bounded ordered source roles required"));
    }
    match target {
        Some(t) => source_target_preflight([sources[0].1, sources[1].1], request, t),
        None => source_preflight([sources[0].1, sources[1].1], request),
    }
}
pub fn checked_product_revision(
    inputs: &[reference::Input],
    request: &Request,
    selection: &Selection,
    target: &Molecule,
    policy: &reference::TargetPolicy,
) -> Result<String, Error> {
    if !hex(&selection.batch_revision)
        || selection.batch_revision != batch_revision(inputs, request)?
        || selection.candidate_id >= 256
        || selection.component_id >= 4
    {
        return Err(bad("invalid component selection binding"));
    }
    hash(&(
        "checked-recombination-component-1",
        check_revision(inputs, request, target, policy)?,
        selection,
    ))
}
fn text(s: &str, n: usize) -> bool {
    s.len() <= n && !s.trim().is_empty()
}
fn findings(values: &[String]) -> Result<(), Error> {
    if values.len() > 1024 || values.iter().any(|s| s.len() > 4096) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn identities(components: &[ComponentIdentity]) -> Result<(), Error> {
    if components.is_empty()
        || components.len() > 4
        || components.iter().enumerate().any(|(n, c)| {
            c.component_id != n
                || !hex(&c.reference_revision)
                || c.length == 0
                || c.length > 2_000_000
        })
    {
        return Err(bad("invalid reciprocal component identities"));
    }
    Ok(())
}
impl ProductDocument {
    /// Historical bounded shape, not authority to materialize any stored/imported claim.
    pub fn validate(&self) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.selection.batch_revision)
            || self.selection.candidate_id >= 256
            || self.selection.component_id >= 4
            || !hex(&self.reference_batch_revision)
            || !hex(&self.reference_product_revision)
            || self.annotation_evidence.len() != 2
        {
            return Err(bad("invalid reference product identity"));
        }
        identities(&self.reciprocal_components)?;
        findings(&self.warnings)?;
        let Some(component) = self.reciprocal_components.get(self.selection.component_id) else {
            return Err(bad("missing selected component"));
        };
        if self.molecule.sequence.alphabet() != collagene_model::Alphabet::Dna
            || component.reference_revision != self.reference_product_revision
            || component.length != self.molecule.sequence.len()
            || component.topology != self.molecule.topology
        {
            return Err(bad("component content identity mismatch"));
        }
        if self.molecule.metadata.len() != 3
            || self.molecule.original_genbank.is_some()
            || self
                .molecule
                .metadata
                .get(REFERENCE_ONLY_KEY)
                .map(String::as_str)
                != Some("true")
            || self
                .molecule
                .metadata
                .get(CHEMISTRY_UNMAPPED_KEY)
                .map(String::as_str)
                != Some("true")
            || self
                .molecule
                .metadata
                .get(RECIPE_KEY)
                .is_none_or(|s| s.len() > 256 * 1024)
        {
            return Err(bad("reference-only producer metadata missing or unbounded"));
        }
        if self.molecule.features.len() > 4096
            || self.annotation_evidence.iter().any(|e| {
                e.unmapped_features.len() > 4096
                    || e.segments.len() > 32
                    || e.co_sourced_segments.len() > 14
            })
        {
            return Err(Error::TooLarge);
        }
        reference::model_preflight(&[&self.molecule])?;
        // Bind actual DNA/topology to the deterministic annotation-free engine
        // identity, rather than accepting only a persisted claimed raw hash.
        let canonical = Molecule::dna(
            format!("reciprocal product {}", self.selection.component_id + 1),
            self.molecule.sequence.as_str(),
            self.molecule.topology,
        )?;
        if canonical.revision() != self.reference_product_revision {
            return Err(bad(
                "actual document DNA/topology differs from raw reference product identity",
            ));
        }
        reference::feature_preflight(
            &self
                .annotation_evidence
                .iter()
                .map(|e| e.unmapped_features.as_slice())
                .collect::<Vec<_>>(),
        )?;
        self.molecule.validate()?;
        for (input, e) in self.annotation_evidence.iter().enumerate() {
            if e.input != input
                || !text(&e.source.role, 256)
                || !hex(&e.source.revision)
                || e.mapped_features > 4096
                || e.reviewed_features > e.mapped_features
                || e.dropped_features > 4096
                || e.unmapped_features.iter().any(|f| !f.needs_review)
            {
                return Err(bad("invalid source annotation evidence"));
            }
            findings(&e.warnings)?;
            for map in e.segments.iter().chain(&e.co_sourced_segments) {
                if map.source_start >= map.source_end
                    || map.source_end > 2_000_000
                    || map
                        .destination_start
                        .checked_add(map.source_end - map.source_start)
                        .is_none_or(|end| end > self.molecule.sequence.len())
                {
                    return Err(bad("invalid annotation map shape"));
                }
            }
        }
        let mapped = self
            .annotation_evidence
            .iter()
            .map(|e| e.mapped_features)
            .sum::<usize>();
        let reviewed = self
            .annotation_evidence
            .iter()
            .map(|e| e.reviewed_features)
            .sum::<usize>();
        if mapped != self.molecule.features.len()
            || reviewed
                != self
                    .molecule
                    .features
                    .iter()
                    .filter(|f| f.needs_review)
                    .count()
        {
            return Err(bad(
                "annotation counters disagree with actual projected features",
            ));
        }
        let recipe: UnverifiedRecipe = serde_json::from_str(&self.molecule.metadata[RECIPE_KEY])
            .map_err(|_| bad("invalid unverified producer recipe"))?;
        recipe.request.validate()?;
        if recipe.schema_version != 1
            || recipe.engine_version != ENGINE_VERSION
            || recipe.selection != self.selection
            || recipe.reference_batch_revision != self.reference_batch_revision
            || recipe.reference_product_revision != self.reference_product_revision
            || recipe.reciprocal_components != self.reciprocal_components
            || recipe.source_identities
                != self
                    .annotation_evidence
                    .iter()
                    .map(|e| e.source.clone())
                    .collect::<Vec<_>>()
        {
            return Err(bad("producer recipe/evidence binding mismatch"));
        }
        size(self, 32 * 1024 * 1024)?;
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        inputs: &[reference::Input],
        request: &Request,
        selection: &Selection,
    ) -> Result<(), Error> {
        self.validate()?;
        let recipe: UnverifiedRecipe = serde_json::from_str(&self.molecule.metadata[RECIPE_KEY])
            .map_err(|_| bad("invalid bounded recipe"))?;
        if &recipe.request != request {
            return Err(bad("producer request differs from current request"));
        }
        if &self.selection != selection
            || selection.batch_revision != batch_revision(inputs, request)?
            || self.reference_batch_revision
                != reference::batch_revision(
                    inputs,
                    &request.profile.profile,
                    request.reaction,
                    &request.limits,
                )?
            || self.annotation_evidence.iter().zip(inputs).any(|(e, i)| {
                e.source.role != i.role
                    || e.source.revision != i.molecule.revision()
                    || e.mapped_features + e.dropped_features + e.unmapped_features.len()
                        != i.molecule.features.len()
            })
        {
            return Err(bad("reference product source/request binding mismatch"));
        }
        Ok(())
    }
}
impl RecombinationReport {
    pub fn validate(&self) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.batch_revision)
            || self.documents.len() > 1024
            || self.annotation_work > 16_000_000
            || self.annotation_copy_bytes > 64 * 1024 * 1024
            || self.annotation_features > 32768
        {
            return Err(Error::TooLarge);
        }
        // Bound all duplicated document payloads before recursive guards, recipe
        // decoding, hashes or per-document serialization. Native callers can forge
        // these rows independently of the small raw reference report.
        let mut bases = 0usize;
        let mut features = 0usize;
        for d in &self.documents {
            bases = bases
                .checked_add(d.molecule.sequence.len())
                .filter(|n| *n <= self.reference.limits.max_product_bases.min(8_000_000))
                .ok_or(Error::TooLarge)?;
            features = features
                .checked_add(d.molecule.features.len())
                .filter(|n| *n <= 32768)
                .ok_or(Error::TooLarge)?;
            if d.molecule.name.len() > 1024
                || d.molecule.metadata.len() != 3
                || d.molecule.original_genbank.is_some()
                || d.annotation_evidence.len() != 2
            {
                return Err(Error::TooLarge);
            }
            for (k, v) in &d.molecule.metadata {
                if k.len() > 128 || v.len() > 256 * 1024 {
                    return Err(Error::TooLarge);
                }
            }
            findings(&d.warnings)?;
            for e in &d.annotation_evidence {
                features = features
                    .checked_add(e.unmapped_features.len())
                    .filter(|n| *n <= 32768)
                    .ok_or(Error::TooLarge)?;
                if e.unmapped_features.len() > 4096
                    || e.segments.len() > 32
                    || e.co_sourced_segments.len() > 14
                    || e.source.role.len() > 256
                    || e.source.revision.len() != 64
                {
                    return Err(Error::TooLarge);
                }
                findings(&e.warnings)?;
            }
        }
        let mut borrowed_nodes = 0usize;
        for d in &self.documents {
            reference::model_preflight(&[&d.molecule])?;
            reference::feature_preflight(
                &d.annotation_evidence
                    .iter()
                    .map(|e| e.unmapped_features.as_slice())
                    .collect::<Vec<_>>(),
            )?;
            borrowed_nodes = borrowed_nodes
                .checked_add(
                    d.molecule
                        .features
                        .iter()
                        .chain(
                            d.annotation_evidence
                                .iter()
                                .flat_map(|e| e.unmapped_features.iter()),
                        )
                        .map(|f| nodes(&f.location) + f.qualifiers.len())
                        .sum::<usize>(),
                )
                .filter(|n| *n <= 1_000_000)
                .ok_or(Error::TooLarge)?;
        }
        // The library raw report has its own pre-serialization tree/byte guards.
        self.reference.validate()?;
        findings(&self.warnings)?;
        size(self, 32 * 1024 * 1024)?;
        if self.documents.len()
            != self
                .reference
                .candidates
                .iter()
                .map(|c| c.products.len())
                .sum::<usize>()
        {
            return Err(bad("not every reciprocal component has a document"));
        }
        let mut count = 0;
        let mut tree_nodes = 0;
        let mut index = 0;
        for c in &self.reference.candidates {
            for (component_id, p) in c.products.iter().enumerate() {
                let d = &self.documents[index];
                index += 1;
                d.validate()?;
                count += d.molecule.features.len()
                    + d.annotation_evidence
                        .iter()
                        .map(|e| e.unmapped_features.len())
                        .sum::<usize>();
                if count > 32768 {
                    return Err(Error::TooLarge);
                }
                tree_nodes += d
                    .molecule
                    .features
                    .iter()
                    .chain(
                        d.annotation_evidence
                            .iter()
                            .flat_map(|e| e.unmapped_features.iter()),
                    )
                    .map(|f| nodes(&f.location) + f.qualifiers.len())
                    .sum::<usize>();
                if tree_nodes > 1_000_000 {
                    return Err(Error::TooLarge);
                }
                if d.selection.batch_revision != self.batch_revision
                    || d.selection.candidate_id != c.id
                    || d.selection.component_id != component_id
                    || d.reference_batch_revision != self.reference.batch_revision
                    || d.reference_product_revision != p.molecule.revision()
                    || d.molecule.sequence != p.molecule.sequence
                    || d.molecule.topology != p.molecule.topology
                {
                    return Err(bad(
                        "annotated document differs from reference reconstruction identity",
                    ));
                }
            }
        }
        size(self, 32 * 1024 * 1024)?;
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        inputs: &[reference::Input],
        request: &Request,
    ) -> Result<(), Error> {
        self.validate()?;
        self.reference.validate_bindings_for(
            inputs,
            &request.profile.profile,
            request.reaction,
            &request.limits,
        )?;
        if self.batch_revision != batch_revision(inputs, request)? {
            return Err(bad("adapter report batch mismatch"));
        }
        if self.annotation_work > request.annotation_limits.max_work
            || self.annotation_copy_bytes > request.annotation_limits.max_copy_bytes
            || self.annotation_features > request.annotation_limits.max_features
        {
            return Err(bad("adapter report exceeds the bound annotation policy"));
        }
        // Current source identities are hashed ONCE, never once per reciprocal
        // document. Self-consistent imported evidence is not current provenance.
        let expected_sources = inputs
            .iter()
            .map(|i| SourceIdentity {
                role: i.role.clone(),
                revision: i.molecule.revision(),
            })
            .collect::<Vec<_>>();
        let expected_feature_counts = inputs
            .iter()
            .map(|i| i.molecule.features.len())
            .collect::<Vec<_>>();
        let expected_charged_features = expected_feature_counts
            .iter()
            .sum::<usize>()
            .checked_mul(self.documents.len())
            .ok_or(Error::TooLarge)?;
        if self.annotation_features != expected_charged_features {
            return Err(bad(
                "reported annotation feature budget differs from current source/document count",
            ));
        }
        for document in &self.documents {
            if document
                .annotation_evidence
                .iter()
                .zip(&expected_feature_counts)
                .any(|(e, &count)| {
                    e.mapped_features + e.dropped_features + e.unmapped_features.len() != count
                })
            {
                return Err(bad(
                    "document annotation coverage counters differ from current source features",
                ));
            }
            if document
                .annotation_evidence
                .iter()
                .map(|e| &e.source)
                .ne(expected_sources.iter())
            {
                return Err(bad(
                    "document source identities differ from current ordered sources",
                ));
            }
            // validate() already bounded and structurally checked this JSON before
            // any native tree recursion, and capped the whole aggregate report.
            let recipe: UnverifiedRecipe =
                serde_json::from_str(&document.molecule.metadata[RECIPE_KEY])
                    .map_err(|_| bad("invalid bounded document recipe"))?;
            if &recipe.request != request || recipe.source_identities != expected_sources {
                return Err(bad(
                    "document producer request/source pin differs from current request",
                ));
            }
        }
        Ok(())
    }
}
impl TargetAssessment {
    pub fn validate(&self) -> Result<(), Error> {
        let r = &self.reference;
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.batch_revision)
            || !hex(&r.batch_revision)
            || !hex(&r.target_revision)
            || r.matches.len() > 2048
        {
            return Err(bad("invalid target assessment shape"));
        }
        findings(&r.warnings)?;
        findings(&r.diagnostics)?;
        let unresolved = !r.search_complete
            || r.possible_sites
            || r.unsupported_geometry
            || r.ambiguous_source_dna
            || r.ambiguous_target_dna;
        let applicability = if unresolved {
            reference::TargetApplicability::Unresolved
        } else {
            reference::TargetApplicability::DeclaredFullSiteProfileOnly
        };
        if r.applicability != applicability
            || r.no_match_proven_under_model != (r.matches.is_empty() && !unresolved)
            || (r.unique_under_declared_profile
                && (!r.search_complete || r.possible_sites || r.unsupported_geometry))
            || (r.ambiguous_target_dna && r.matches.iter().any(|m| m.fully_specified_dna))
        {
            return Err(bad(
                "contradictory target completeness/uncertainty/applicability fields",
            ));
        }
        if r.matches
            .iter()
            .any(|m| m.candidate >= 256 || m.product >= 4 || m.target_origin > 2_000_000)
        {
            return Err(bad("invalid target match shape"));
        }
        size(self, 4 * 1024 * 1024)?;
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        inputs: &[reference::Input],
        request: &Request,
        target: &Molecule,
        policy: &reference::TargetPolicy,
    ) -> Result<(), Error> {
        self.validate()?;
        let expected_revision = check_revision(inputs, request, target, policy)?;
        // Actual ambiguity can be checked once without site search or DP replay.
        let ambiguous_source = inputs.iter().any(|i| {
            i.molecule
                .sequence
                .as_str()
                .bytes()
                .any(|b| !b"ACGT".contains(&b))
        });
        let ambiguous_target = target
            .sequence
            .as_str()
            .bytes()
            .any(|b| !b"ACGT".contains(&b));
        if self.reference.ambiguous_source_dna != ambiguous_source
            || self.reference.ambiguous_target_dna != ambiguous_target
        {
            return Err(bad(
                "target ambiguity flags differ from actual source/target DNA",
            ));
        }
        if self.batch_revision != expected_revision
            || self.reference.batch_revision
                != reference::check_batch_revision(
                    inputs,
                    &request.profile.profile,
                    request.reaction,
                    &request.limits,
                    target,
                    *policy,
                )?
            || self.reference.target_revision != target.revision()
            || self.reference.policy != *policy
        {
            return Err(bad("target assessment request/source mismatch"));
        }
        Ok(())
    }
}
impl CheckedProduct {
    pub fn validate(&self) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION || !hex(&self.batch_revision) {
            return Err(bad("invalid checked component shape"));
        }
        self.product.validate()?;
        self.check.validate()?;
        if self.batch_revision
            != hash(&(
                "checked-recombination-component-1",
                &self.check.batch_revision,
                &self.product.selection,
            ))?
            || !self.check.reference.matches.iter().any(|m| {
                m.candidate == self.product.selection.candidate_id
                    && m.product == self.product.selection.component_id
                    && m.fully_specified_dna
            })
        {
            return Err(bad("checked result lacks a matching selected component"));
        }
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        inputs: &[reference::Input],
        request: &Request,
        selection: &Selection,
        target: &Molecule,
        policy: &reference::TargetPolicy,
    ) -> Result<(), Error> {
        self.validate()?;
        self.product
            .validate_bindings_for(inputs, request, selection)?;
        self.check
            .validate_bindings_for(inputs, request, target, policy)?;
        if !literal_target_match(&self.product.molecule, target, policy) {
            return Err(bad(
                "actual selected DNA does not match the independent target under the declared policy",
            ));
        }
        if self.batch_revision
            != checked_product_revision(inputs, request, selection, target, policy)?
        {
            return Err(bad("checked component request mismatch"));
        }
        Ok(())
    }
}

// Bounded literal comparison only; no recombination scan or biological replay.
fn literal_target_match(
    product: &Molecule,
    target: &Molecule,
    policy: &reference::TargetPolicy,
) -> bool {
    if product.topology != target.topology
        || product.sequence.len() != target.sequence.len()
        || product
            .sequence
            .as_str()
            .bytes()
            .any(|b| !b"ACGT".contains(&b))
        || target
            .sequence
            .as_str()
            .bytes()
            .any(|b| !b"ACGT".contains(&b))
    {
        return false;
    }
    let source = product.sequence.as_str();
    let matches = |query: &str| {
        if product.topology == Topology::Circular {
            let mut doubled = String::with_capacity(source.len() * 2);
            doubled.push_str(source);
            doubled.push_str(source);
            doubled[..doubled.len() - 1].contains(query)
        } else {
            source == query
        }
    };
    if matches(target.sequence.as_str()) {
        return true;
    }
    policy.allow_reverse_complement
        && target
            .sequence
            .reverse_complement()
            .is_ok_and(|q| matches(q.as_str()))
}

pub mod multisite;
