//! Source-bound multi-parent reciprocal reference exchange and conservative
//! annotation projection. No imported recipe establishes chemical strand state.
use crate::cloning::recombination::{self as adapter, reference};
use crate::cloning::recombination::{
    AnnotationLimits, ComponentIdentity, Selection, SourceAnnotationEvidence, SourceIdentity,
};
use collagene_model::annotations;
use collagene_model::{Alphabet, Error, Molecule, Topology};
pub use reference::multisite;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};

pub const ENGINE_VERSION: &str = "multisite-reference-adapter-2";
pub const RECIPE_KEY: &str = "multisite_recombination:unverified_recipe";
pub const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
const MAX_RECIPE_BYTES: usize = 256 * 1024;
const MAX_INPUTS: usize = 5;
const MAX_DOCUMENTS: usize = 256 * 32;
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request{pub reaction:multisite::Request,pub annotation_limits:AnnotationLimits}}
dto! {pub struct Batch{pub engine_version:String,pub batch_revision:String,pub request:Request,pub sources:Vec<SourceIdentity>,pub reference:multisite::Report,pub documents:Vec<ProductDocument>,pub warnings:Vec<String>}}
dto! {pub struct ProductDocument{pub engine_version:String,pub selection:Selection,pub reference_product_revision:String,pub reciprocal_components:Vec<ComponentIdentity>,pub molecule:Molecule,pub annotation_evidence:Vec<SourceAnnotationEvidence>,pub warnings:Vec<String>}}
dto! {pub struct Assessment{pub target_revision:String,pub target_matches:bool,pub complete:bool,pub findings:Vec<String>}}
dto! {pub struct CheckedProduct{pub product:ProductDocument,pub assessment:Assessment}}
dto! {pub struct UnverifiedRecipe{pub schema_version:u16,pub engine_version:String,pub request:Request,pub selection:Selection,pub source_identities:Vec<SourceIdentity>,pub reference_batch_revision:String,pub reference_product_revision:String,pub reciprocal_components:Vec<ComponentIdentity>}}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("multisite reference recombination: {message}"))
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
struct Counter {
    bytes: usize,
    max: usize,
}
impl Write for Counter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|n| *n <= self.max)
            .ok_or_else(|| std::io::Error::other("bounded multisite evidence"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn size(value: &impl Serialize, max: usize) -> Result<usize, Error> {
    let mut out = Counter { bytes: 0, max };
    serde_json::to_writer(&mut out, value).map_err(|_| Error::TooLarge)?;
    Ok(out.bytes)
}
struct HashWriter {
    hash: blake3::Hasher,
    bytes: usize,
}
impl Write for HashWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|n| *n <= MAX_REPORT_BYTES)
            .ok_or_else(|| std::io::Error::other("bounded multisite hash"))?;
        self.hash.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn hash(value: &impl Serialize) -> Result<String, Error> {
    let mut out = HashWriter {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(&mut out, value).map_err(|_| Error::TooLarge)?;
    Ok(out.hash.finalize().to_hex().to_string())
}
fn texts(values: &[String]) -> Result<(), Error> {
    if values.len() > 1024 || values.iter().any(|s| s.len() > 4096) {
        Err(Error::TooLarge)
    } else {
        Ok(())
    }
}
fn sources(inputs: &[reference::Input]) -> Vec<SourceIdentity> {
    inputs
        .iter()
        .map(|input| SourceIdentity {
            role: input.role.clone(),
            revision: input.molecule.revision(),
        })
        .collect()
}
fn request_shape(request: &Request) -> Result<(), Error> {
    request.annotation_limits.validate()?;
    request.reaction.profile.validate()?;
    let l = &request.reaction.limits;
    if l.max_input_bases == 0
        || l.max_input_bases > 2_000_000
        || l.max_comparisons == 0
        || l.max_comparisons > 100_000_000
        || l.max_sites == 0
        || l.max_sites > 1024
        || l.max_combinations == 0
        || l.max_combinations > 65536
        || l.max_candidates == 0
        || l.max_candidates > 256
        || l.max_components_per_candidate == 0
        || l.max_components_per_candidate > 32
        || l.max_product_bases == 0
        || l.max_product_bases > 8_000_000
        || l.max_report_bytes < 1024
        || l.max_report_bytes > MAX_REPORT_BYTES
        || request.reaction.specificities.is_empty()
        || request.reaction.specificities.len() > 16
    {
        return Err(Error::TooLarge);
    }
    let mut ids = BTreeSet::new();
    for id in &request.reaction.specificities {
        if !(1..=16).contains(id) || !ids.insert(*id) {
            return Err(bad("specificity IDs must be distinct in 1–16"));
        }
    }
    size(request, 2 * 1024 * 1024)?;
    Ok(())
}
/// Borrowed declaration bounds before any source cloning, hashing or simulation.
pub fn request_preflight(request: &Request) -> Result<(), Error> {
    request_shape(request)
}

fn binding(raw: &str, request: &Request, sources: &[SourceIdentity]) -> Result<String, Error> {
    hash(&(ENGINE_VERSION, raw, request, sources))
}
/// Borrowed aggregate source guards before rich input cloning or model hashes.
pub fn source_preflight(molecules: &[&Molecule], request: &Request) -> Result<(), Error> {
    request_shape(request)?;
    if !(2..=MAX_INPUTS).contains(&molecules.len()) {
        return Err(Error::TooLarge);
    }
    let mut bases = 0usize;
    for molecule in molecules {
        if molecule.sequence.alphabet() != Alphabet::Dna || molecule.sequence.is_empty() {
            return Err(bad("nonempty DNA parents are required"));
        }
        bases = bases
            .checked_add(molecule.sequence.len())
            .filter(|n| *n <= request.reaction.limits.max_input_bases)
            .ok_or(Error::TooLarge)?;
    }
    reference::model_preflight(molecules)?;
    for molecule in molecules {
        molecule.validate()?;
    }
    Ok(())
}
pub fn preflight(inputs: &[reference::Input], request: &Request) -> Result<(), Error> {
    if inputs.is_empty() || inputs.len() > MAX_INPUTS {
        return Err(Error::TooLarge);
    }
    source_preflight(
        &inputs
            .iter()
            .map(|input| &input.molecule)
            .collect::<Vec<_>>(),
        request,
    )?;
    multisite::preflight(inputs, &request.reaction)
}
pub fn revision(inputs: &[reference::Input], request: &Request) -> Result<String, Error> {
    preflight(inputs, request)?;
    binding(
        &multisite::revision(inputs, &request.reaction)?,
        request,
        &sources(inputs),
    )
}
fn selection_shape(selection: &Selection) -> Result<(), Error> {
    if !hex(&selection.batch_revision)
        || selection.candidate_id >= 256
        || selection.component_id >= 32
    {
        return Err(bad("selected candidate/component shape is invalid"));
    }
    Ok(())
}
fn selected_binding(selection: &Selection) -> Result<String, Error> {
    hash(&(
        ENGINE_VERSION,
        "selected-reference",
        &selection.batch_revision,
        selection,
    ))
}
pub fn product_revision(
    inputs: &[reference::Input],
    request: &Request,
    selection: &Selection,
) -> Result<String, Error> {
    selection_shape(selection)?;
    if selection.batch_revision != revision(inputs, request)? {
        return Err(bad(
            "selected component belongs to another source/profile/options batch",
        ));
    }
    selected_binding(selection)
}
pub fn check_revision(
    inputs: &[reference::Input],
    request: &Request,
    selection: &Selection,
    target: &Molecule,
) -> Result<String, Error> {
    preflight(inputs, request)?;
    reference::model_preflight(&[target])?;
    target.validate()?;
    if target.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("independent target must be DNA"));
    }
    hash(&(
        ENGINE_VERSION,
        "independent-target",
        product_revision(inputs, request, selection)?,
        target.revision(),
    ))
}

fn components(products: &[reference::Product]) -> Vec<ComponentIdentity> {
    products
        .iter()
        .enumerate()
        .map(|(component_id, p)| ComponentIdentity {
            component_id,
            reference_revision: p.molecule.revision(),
            length: p.molecule.sequence.len(),
            topology: p.molecule.topology,
        })
        .collect()
}
fn reserve(work: &mut adapter::Work, bytes: usize, limits: &AnnotationLimits) -> Result<(), Error> {
    work.bytes = work
        .bytes
        .checked_add(bytes)
        .filter(|n| *n <= limits.max_copy_bytes)
        .ok_or(Error::TooLarge)?;
    Ok(())
}
fn document(
    inputs: &[reference::Input],
    request: &Request,
    raw_batch: &str,
    batch: &str,
    candidate_id: usize,
    component_id: usize,
    raw: &reference::Product,
    components: &[ComponentIdentity],
    identities: &[SourceIdentity],
    work: &mut adapter::Work,
) -> Result<ProductDocument, Error> {
    // Reserve every source's mapping/copy budget before cloning any DNA/features.
    let mut mappings = Vec::with_capacity(inputs.len());
    for source in inputs.iter().enumerate() {
        let mapping = adapter::annotation_maps(raw, source.0);
        let count = mapping.as_ref().map_or(
            raw.segments
                .len()
                .saturating_add(raw.junctions.len().saturating_mul(14)),
            |(maps, _)| maps.len(),
        );
        adapter::charge(&source.1.molecule, count, &request.annotation_limits, work)?;
        mappings.push(mapping);
    }
    reserve(
        work,
        size(&raw.molecule, MAX_REPORT_BYTES)?
            .checked_mul(2)
            .ok_or(Error::TooLarge)?,
        &request.annotation_limits,
    )?;
    let selection = Selection {
        batch_revision: batch.into(),
        candidate_id,
        component_id,
    };
    let reference_product_revision = raw.molecule.revision();
    let recipe_bytes = size(
        &(request, &selection, identities, components),
        MAX_RECIPE_BYTES,
    )?;
    reserve(
        work,
        recipe_bytes
            .checked_add(1024)
            .and_then(|n| n.checked_mul(2))
            .ok_or(Error::TooLarge)?,
        &request.annotation_limits,
    )?;
    let mut molecule = Molecule::dna(
        raw.molecule.name.clone(),
        raw.molecule.sequence.as_str(),
        raw.molecule.topology,
    )?;
    let mut evidence = Vec::with_capacity(inputs.len());
    let mut warnings=vec!["Reference-only multi-parent exchange. Chemical strand ancestry, nicks, methylation, host repair, reagent efficacy and a selected experimental winner are not established.".into()];
    for (input, (source, mapping)) in inputs.iter().zip(mappings).enumerate() {
        let (segments, co, projection) = match mapping {
            Ok((segments, co)) => {
                let projected = annotations::project_features_with_topology(
                    &source.molecule,
                    &segments,
                    molecule.sequence.len(),
                    molecule.topology,
                );
                (segments, co, projected)
            }
            Err(error) => (vec![], vec![], Err(error)),
        };
        let (mapped_features, reviewed_features, dropped_features, unmapped_features, findings) =
            match projection {
                Ok(mut projection) => {
                    adapter::mark_co_source(&mut projection.features, &co);
                    let count = projection.features.len();
                    let reviewed = projection
                        .features
                        .iter()
                        .filter(|f| f.needs_review)
                        .count();
                    molecule.features.extend(projection.features);
                    (
                        count,
                        reviewed,
                        projection.dropped_features,
                        projection.unmapped_features,
                        projection.warnings,
                    )
                }
                Err(error) => {
                    let mut features = source.molecule.features.clone();
                    for feature in &mut features {
                        feature.needs_review = true;
                    }
                    (
                        0,
                        0,
                        0,
                        features,
                        vec![format!(
                            "Source '{}' annotations retained as reviewed unmapped evidence: {error}",
                            source.role
                        )],
                    )
                }
            };
        evidence.push(SourceAnnotationEvidence {
            input,
            source: identities[input].clone(),
            segments,
            co_sourced_segments: co,
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
            warnings.push(format!("Source '{}' physical/modification claims were not propagated. Retain its original immutable document.",source.role));
        }
    }
    molecule
        .metadata
        .insert(adapter::REFERENCE_ONLY_KEY.into(), "true".into());
    molecule
        .metadata
        .insert(adapter::CHEMISTRY_UNMAPPED_KEY.into(), "true".into());
    let recipe = UnverifiedRecipe {
        schema_version: 1,
        engine_version: ENGINE_VERSION.into(),
        request: request.clone(),
        selection: selection.clone(),
        source_identities: identities.to_vec(),
        reference_batch_revision: raw_batch.into(),
        reference_product_revision: reference_product_revision.clone(),
        reciprocal_components: components.to_vec(),
    };
    let metadata = serde_json::to_string(&recipe).map_err(|_| Error::TooLarge)?;
    if metadata.len() > MAX_RECIPE_BYTES {
        return Err(Error::TooLarge);
    }
    molecule.metadata.insert(RECIPE_KEY.into(), metadata);
    molecule.validate()?;
    Ok(ProductDocument {
        engine_version: ENGINE_VERSION.into(),
        selection,
        reference_product_revision,
        reciprocal_components: components.to_vec(),
        molecule,
        annotation_evidence: evidence,
        warnings,
    })
}
pub fn simulate(inputs: &[reference::Input], request: &Request) -> Result<Batch, Error> {
    preflight(inputs, request)?;
    let raw = multisite::recombine(inputs, &request.reaction)?;
    let identities = sources(inputs);
    let batch_revision = binding(&raw.batch_revision, request, &identities)?;
    let mut documents = Vec::new();
    let mut work = adapter::Work::default();
    for candidate in &raw.candidates {
        let components = components(&candidate.products);
        for (component, product) in candidate.products.iter().enumerate() {
            if documents.len() >= MAX_DOCUMENTS {
                return Err(Error::TooLarge);
            }
            documents.push(document(
                inputs,
                request,
                &raw.batch_revision,
                &batch_revision,
                candidate.id,
                component,
                product,
                &components,
                &identities,
                &mut work,
            )?);
        }
    }
    let out=Batch{engine_version:ENGINE_VERSION.into(),batch_revision:batch_revision.clone(),request:request.clone(),sources:identities,reference:raw,documents,warnings:vec!["All declared reciprocal reference components and uncertainty are retained. Source-coordinate features are correspondence, not chemical inheritance; stored producer recipes never restore current authority.".into()]};
    out.validate_binding(&batch_revision)?;
    Ok(out)
}
pub fn product(
    inputs: &[reference::Input],
    request: &Request,
    selection: &Selection,
) -> Result<ProductDocument, Error> {
    preflight(inputs, request)?;
    if selection.batch_revision != revision(inputs, request)? {
        return Err(bad(
            "selection belongs to a different ordered source/profile/options request",
        ));
    }
    simulate(inputs, request)?.selected(selection)
}
fn target_matches(product: &Molecule, target: &Molecule) -> bool {
    product.topology == target.topology
        && product.sequence.len() == target.sequence.len()
        && product
            .sequence
            .as_str()
            .bytes()
            .all(|b| b"ACGT".contains(&b))
        && target
            .sequence
            .as_str()
            .bytes()
            .all(|b| b"ACGT".contains(&b))
        && (product.sequence == target.sequence
            || product.topology == Topology::Circular
                && crate::sequence::analysis::rotation_offset(
                    product.sequence.as_str(),
                    target.sequence.as_str(),
                )
                .is_some())
}
pub fn check_target(
    inputs: &[reference::Input],
    request: &Request,
    selection: &Selection,
    target: &Molecule,
) -> Result<CheckedProduct, Error> {
    preflight(inputs, request)?;
    reference::model_preflight(&[target])?;
    target.validate()?;
    if target.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("independent target must be DNA"));
    }
    let batch = simulate(inputs, request)?;
    let product = batch.selected(selection)?;
    let complete = batch.reference.complete
        && batch.reference.unsupported.is_empty()
        && batch.reference.sites.iter().all(|s| s.certain);
    let target_matches = target_matches(&product.molecule, target);
    let out=CheckedProduct{product,assessment:Assessment{target_revision:target.revision(),target_matches,complete,findings:vec![if target_matches{"Selected reference DNA/topology agrees with the independent target under oriented circular-rotation comparison.".into()}else{"Selected reference DNA/topology differs from the independent target.".into()},"Target comparison does not establish chemical strand state, repaired DNA, reaction efficacy, mixture uniqueness or experimental success.".into(),if complete{"Declared reference search completed within its limits.".into()}else{"Reference search remains incomplete or uncertain; target agreement cannot close its search/unsupported-state limits.".into()}]}};
    out.validate_binding(&selection.batch_revision)?;
    Ok(out)
}
fn source_shape(sources: &[SourceIdentity]) -> Result<(), Error> {
    if sources.is_empty() || sources.len() > MAX_INPUTS {
        return Err(Error::TooLarge);
    }
    let mut roles = BTreeSet::new();
    for source in sources {
        if source.role.trim().is_empty()
            || source.role.len() > 256
            || !roles.insert(&source.role)
            || !hex(&source.revision)
        {
            return Err(bad("invalid source identities"));
        }
    }
    Ok(())
}
fn components_shape(components: &[ComponentIdentity]) -> Result<(), Error> {
    if components.is_empty() || components.len() > 32 {
        return Err(Error::TooLarge);
    }
    for (i, c) in components.iter().enumerate() {
        if c.component_id != i
            || !hex(&c.reference_revision)
            || c.length == 0
            || c.length > 2_000_000
        {
            return Err(bad("invalid reciprocal component identity"));
        }
    }
    Ok(())
}
fn document_shape(doc: &ProductDocument) -> Result<(), Error> {
    if doc.molecule.name.len() > 1024
        || doc.molecule.sequence.len() > 2_000_000
        || doc.molecule.features.len() > 32768
        || doc.molecule.metadata.len() != 3
        || doc.molecule.original_genbank.is_some()
        || doc.annotation_evidence.is_empty()
        || doc.annotation_evidence.len() > MAX_INPUTS
    {
        return Err(Error::TooLarge);
    }
    components_shape(&doc.reciprocal_components)?;
    texts(&doc.warnings)?;
    let mut unmapped = 0usize;
    for e in &doc.annotation_evidence {
        unmapped = unmapped
            .checked_add(e.unmapped_features.len())
            .filter(|n| *n <= 32768)
            .ok_or(Error::TooLarge)?;
        if e.segments.len() > 2048
            || e.co_sourced_segments.len() > 256
            || e.mapped_features > 32768
            || e.reviewed_features > e.mapped_features
            || e.dropped_features > 32768
        {
            return Err(Error::TooLarge);
        }
        texts(&e.warnings)?;
    }
    reference::model_preflight(&[&doc.molecule])?;
    reference::feature_preflight(
        &doc.annotation_evidence
            .iter()
            .map(|e| e.unmapped_features.as_slice())
            .collect::<Vec<_>>(),
    )?;
    for (key, value) in &doc.molecule.metadata {
        if key.len() > 128 || value.len() > MAX_RECIPE_BYTES {
            return Err(Error::TooLarge);
        }
    }
    Ok(())
}
impl ProductDocument {
    /// Selected reply context, without source reconstruction.
    pub fn context_revision(&self) -> Result<String, Error> {
        selection_shape(&self.selection)?;
        selected_binding(&self.selection)
    }
    /// Bounded documentary consistency, never fresh recombination or chemistry.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        document_shape(self)?;
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.selection.batch_revision != expected
            || !hex(&self.reference_product_revision)
        {
            return Err(bad("product context differs"));
        }
        if self.molecule.sequence.alphabet() != Alphabet::Dna
            || self
                .molecule
                .metadata
                .get(adapter::REFERENCE_ONLY_KEY)
                .map(String::as_str)
                != Some("true")
            || self
                .molecule
                .metadata
                .get(adapter::CHEMISTRY_UNMAPPED_KEY)
                .map(String::as_str)
                != Some("true")
        {
            return Err(bad("reference-only markers are absent"));
        }
        let Some(component) = self.reciprocal_components.get(self.selection.component_id) else {
            return Err(bad("selected reciprocal component is absent"));
        };
        if component.reference_revision != self.reference_product_revision
            || component.length != self.molecule.sequence.len()
            || component.topology != self.molecule.topology
        {
            return Err(bad("selected component content differs"));
        }
        let canonical = Molecule::dna(
            self.molecule.name.clone(),
            self.molecule.sequence.as_str(),
            self.molecule.topology,
        )?;
        if canonical.revision() != self.reference_product_revision {
            return Err(bad(
                "reference document DNA/topology/name differs from raw product identity",
            ));
        }
        let Some(metadata) = self.molecule.metadata.get(RECIPE_KEY) else {
            return Err(bad("unverified producer recipe absent"));
        };
        let recipe: UnverifiedRecipe = serde_json::from_str(metadata)
            .map_err(|_| bad("invalid unverified producer recipe"))?;
        request_shape(&recipe.request)?;
        source_shape(&recipe.source_identities)?;
        if recipe.schema_version != 1
            || recipe.engine_version != ENGINE_VERSION
            || recipe.selection != self.selection
            || recipe.reference_product_revision != self.reference_product_revision
            || recipe.reciprocal_components != self.reciprocal_components
            || !hex(&recipe.reference_batch_revision)
            || binding(
                &recipe.reference_batch_revision,
                &recipe.request,
                &recipe.source_identities,
            )? != expected
            || recipe.source_identities.len() != self.annotation_evidence.len()
        {
            return Err(bad("unverified recipe/source/component binding differs"));
        }
        let mut mapped = 0usize;
        let mut reviewed = 0usize;
        for (input, e) in self.annotation_evidence.iter().enumerate() {
            if e.input != input
                || e.source != recipe.source_identities[input]
                || e.unmapped_features.iter().any(|f| !f.needs_review)
            {
                return Err(bad(
                    "annotation source identity or retained-review state differs",
                ));
            }
            mapped = mapped
                .checked_add(e.mapped_features)
                .ok_or(Error::TooLarge)?;
            reviewed = reviewed
                .checked_add(e.reviewed_features)
                .ok_or(Error::TooLarge)?;
            for map in e.segments.iter().chain(&e.co_sourced_segments) {
                if map.source_start >= map.source_end
                    || map.source_end > 2_000_000
                    || map
                        .destination_start
                        .checked_add(map.source_end - map.source_start)
                        .is_none_or(|end| end > self.molecule.sequence.len())
                {
                    return Err(bad("annotation map is outside bounded strands"));
                }
            }
        }
        if mapped != self.molecule.features.len()
            || reviewed
                != self
                    .molecule
                    .features
                    .iter()
                    .filter(|f| f.needs_review)
                    .count()
        {
            return Err(bad("annotation counts differ from document"));
        }
        self.molecule.validate()?;
        size(self, recipe.request.annotation_limits.max_report_bytes)?;
        Ok(())
    }
}
impl Batch {
    /// Readonly source/request/raw-product bindings; this never reruns a scan.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        request_shape(&self.request)?;
        source_shape(&self.sources)?;
        texts(&self.warnings)?;
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.batch_revision != expected
            || self.documents.len() > MAX_DOCUMENTS
        {
            return Err(bad("batch context differs"));
        }
        self.reference
            .validate_binding(&self.reference.batch_revision)?;
        if self.reference.profile_revision != self.request.reaction.profile.revision()?
            || self.reference.reaction != self.request.reaction.reaction
            || self.reference.input_revisions
                != self
                    .sources
                    .iter()
                    .map(|s| s.revision.clone())
                    .collect::<Vec<_>>()
            || binding(&self.reference.batch_revision, &self.request, &self.sources)? != expected
        {
            return Err(bad("raw request/source binding differs"));
        }
        let expected_count = self
            .reference
            .candidates
            .iter()
            .try_fold(0usize, |sum, c| {
                sum.checked_add(c.products.len()).ok_or(Error::TooLarge)
            })?;
        if expected_count != self.documents.len() {
            return Err(bad("reciprocal product documents missing or duplicated"));
        }
        let mut seen = BTreeSet::new();
        let mut bases = 0usize;
        let mut features = 0usize;
        let mut unmapped = 0usize;
        let mut source_features = 0usize;
        for document in &self.documents {
            document_shape(document)?;
            bases = bases
                .checked_add(document.molecule.sequence.len())
                .filter(|n| *n <= self.request.reaction.limits.max_product_bases)
                .ok_or(Error::TooLarge)?;
            features = features
                .checked_add(document.molecule.features.len())
                .filter(|n| *n <= self.request.annotation_limits.max_features)
                .ok_or(Error::TooLarge)?;
            for e in &document.annotation_evidence {
                source_features = source_features
                    .checked_add(e.mapped_features)
                    .and_then(|n| n.checked_add(e.dropped_features))
                    .and_then(|n| n.checked_add(e.unmapped_features.len()))
                    .filter(|n| *n <= self.request.annotation_limits.max_features)
                    .ok_or(Error::TooLarge)?;
                unmapped = unmapped
                    .checked_add(e.unmapped_features.len())
                    .filter(|n| *n <= self.request.annotation_limits.max_features)
                    .ok_or(Error::TooLarge)?;
            }
            if !seen.insert((
                document.selection.candidate_id,
                document.selection.component_id,
            )) {
                return Err(bad("duplicate candidate/component document"));
            }
        }
        size(self, self.request.annotation_limits.max_report_bytes)?;
        for document in &self.documents {
            document.validate_binding(expected)?;
            let candidate = self
                .reference
                .candidates
                .iter()
                .find(|c| c.id == document.selection.candidate_id)
                .ok_or_else(|| bad("document candidate absent from raw report"))?;
            let raw = candidate
                .products
                .get(document.selection.component_id)
                .ok_or_else(|| bad("document component absent from raw report"))?;
            if document.reference_product_revision != raw.molecule.revision()
                || document.molecule.sequence != raw.molecule.sequence
                || document.molecule.topology != raw.molecule.topology
                || document.molecule.name != raw.molecule.name
                || document.reciprocal_components != components(&candidate.products)
                || document.annotation_evidence.len() != self.sources.len()
                || document
                    .annotation_evidence
                    .iter()
                    .zip(&self.sources)
                    .any(|(e, s)| e.source != *s)
            {
                return Err(bad(
                    "annotated document differs from actual raw reference product/source",
                ));
            }
        }
        Ok(())
    }
    /// Select an existing document from this validated batch, without a new scan.
    pub fn selected(&self, selection: &Selection) -> Result<ProductDocument, Error> {
        if selection.batch_revision != self.batch_revision {
            return Err(bad(
                "selection is from another source/profile/options batch",
            ));
        }
        self.validate_binding(&selection.batch_revision)?;
        self.documents
            .iter()
            .find(|d| d.selection == *selection)
            .cloned()
            .ok_or_else(|| bad("selected candidate/component was not reconstructed"))
    }
}
impl CheckedProduct {
    /// Selected-source and independent-target reply context; no source replay.
    pub fn context_revision(&self) -> Result<String, Error> {
        selection_shape(&self.product.selection)?;
        if !hex(&self.assessment.target_revision) {
            return Err(bad("target revision is invalid"));
        }
        hash(&(
            ENGINE_VERSION,
            "independent-target",
            selected_binding(&self.product.selection)?,
            &self.assessment.target_revision,
        ))
    }

    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.product.validate_binding(expected)?;
        if !hex(&self.assessment.target_revision) {
            return Err(bad("target revision is invalid"));
        }
        texts(&self.assessment.findings)?;
        size(self, MAX_REPORT_BYTES)?;
        Ok(())
    }
}
