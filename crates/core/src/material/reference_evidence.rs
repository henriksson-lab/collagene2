//! Immutable historical multi-parent reference graphs, separately bound to exact run events.
//! Reading this attachment never restores fresh scientific/material authority.
use crate::{
    cloning::assembly::recipe as assembly_recipe,
    documents::workspace::{Workspace, WorkspaceAttachment},
    material::lineage::{
        self as run_lineage, EventBinding, InputEventBinding, ProductBinding, RunRecord,
    },
    planning::plan::{Operation, Source},
    verification::evidence::ReferenceStepReport,
};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const KIND: &str = "reference-design-run-evidence-v1";
pub const MAX_BYTES: usize = 32 * 1024 * 1024;
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Evidence {
    pub schema_version: u32,
    pub run_id: String,
    pub plan_revision: String,
    pub engine_version: String,
    pub target: EventBinding,
    pub inputs: Vec<InputEventBinding>,
    pub products: Vec<ProductBinding>,
    pub steps: Vec<ReferenceStepReport>,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("Historical reference design evidence: {s}"))
}
struct Bounded(Vec<u8>);
impl Write for Bounded {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > MAX_BYTES)
        {
            return Err(std::io::Error::other("reference run evidence exceeds32MiB"));
        }
        self.0
            .try_reserve(b.len())
            .map_err(|_| std::io::Error::other("reference evidence allocation"))?;
        self.0.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encode(e: &Evidence) -> Result<Vec<u8>, Error> {
    let mut out = Bounded(vec![]);
    serde_json::to_writer(&mut out, e).map_err(|_| Error::TooLarge)?;
    Ok(out.0)
}
fn decode(bytes: &[u8]) -> Result<Evidence, Error> {
    if bytes.len() > MAX_BYTES {
        return Err(Error::TooLarge);
    };
    serde_json::from_slice(bytes).map_err(|e| bad(&e.to_string()))
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn texts(values: &[String]) -> Result<(), Error> {
    if values.len() > 1024 || values.iter().any(|v| v.len() > 4096) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn hit_shape(
    hit: &crate::cloning::recombination::reference::SiteHit,
    inputs: usize,
) -> Result<(), Error> {
    if hit.input >= inputs
        || hit.variant_id.is_empty()
        || hit.variant_id.len() > 256
        || hit.length > 4096
        || hit.start >= 2_000_000
        || hit.switch_boundary > 2_000_000
        || hit.overlap_positions.len() != 7
        || hit.overlap_positions.iter().any(|p| *p >= 2_000_000)
    {
        return Err(bad("historical site shape exceeds documentary bounds"));
    }
    Ok(())
}
/// Structural/documentary limits only. In particular, no current profile scan,
/// family classification, graph reconstruction or engine-label equality runs.
fn batch_shape(batch: &crate::cloning::recombination::multisite::Batch) -> Result<(), Error> {
    if batch.engine_version.is_empty()
        || batch.engine_version.len() > 128
        || !hex(&batch.batch_revision)
        || batch.sources.len() < 2
        || batch.sources.len() > 5
        || batch.documents.len() > 8192
        || batch.reference.engine.is_empty()
        || batch.reference.engine.len() > 128
        || !hex(&batch.reference.batch_revision)
        || !hex(&batch.reference.profile_revision)
        || batch.reference.input_revisions.len() != batch.sources.len()
        || batch.reference.sites.len() > 1024
        || batch.reference.candidates.len() > 256
    {
        return Err(Error::TooLarge);
    }
    texts(&batch.warnings)?;
    texts(&batch.reference.warnings)?;
    texts(&batch.reference.limits)?;
    texts(&batch.reference.unsupported)?;
    let profile = &batch.request.reaction.profile;
    if profile.id.is_empty()
        || profile.id.len() > 256
        || profile.edition.is_empty()
        || profile.edition.len() > 256
        || profile.sources.len() > 64
        || profile.variants.len() > 64
        || profile.limitations.len() > 64
        || batch.request.reaction.specificities.len() > 16
    {
        return Err(Error::TooLarge);
    }
    texts(&profile.limitations)?;
    for source in &profile.sources {
        if source.id.len() > 256
            || source.url.len() > 4096
            || source.map_url.len() > 4096
            || source.vector_name.len() > 1024
            || source.raw_sha256.len() != 64
            || source.map_sha256.len() != 64
        {
            return Err(Error::TooLarge);
        }
    }
    for variant in &profile.variants {
        if variant.id.len() > 256 || variant.source_id.len() > 256 || variant.sequence.len() > 4096
        {
            return Err(Error::TooLarge);
        }
    }
    let count = batch.sources.len();
    for (identity, revision) in batch.sources.iter().zip(&batch.reference.input_revisions) {
        if identity.role.is_empty()
            || identity.role.len() > 256
            || !hex(&identity.revision)
            || identity.revision != *revision
        {
            return Err(bad("historical source identities differ"));
        }
    }
    for hit in &batch.reference.sites {
        hit_shape(hit, count)?;
    }
    let mut bases = 0usize;
    let mut items = 0usize;
    for candidate in &batch.reference.candidates {
        if candidate.sites.len() > 32 || candidate.pairs.len() > 16 || candidate.products.len() > 32
        {
            return Err(Error::TooLarge);
        }
        for hit in &candidate.sites {
            hit_shape(hit, count)?;
        }
        for pair in &candidate.pairs {
            if pair.iter().any(|index| *index >= candidate.sites.len()) {
                return Err(bad("historical pair index outside retained sites"));
            }
        }
        for product in &candidate.products {
            if product.molecule.name.len() > 1024
                || product.molecule.sequence.len() > 2_000_000
                || !product.molecule.features.is_empty()
                || !product.molecule.metadata.is_empty()
                || product.molecule.original_genbank.is_some()
                || product.segments.len() > 74
                || product.junctions.len() > 32
            {
                return Err(Error::TooLarge);
            }
            bases = bases
                .checked_add(product.molecule.sequence.len())
                .filter(|n| *n <= 8_000_000)
                .ok_or(Error::TooLarge)?;
            for map in &product.segments {
                if map.input >= count
                    || map.source_start >= map.source_end
                    || map.source_end > 2_000_000
                    || map.output_start > map.output_end
                    || map.output_end > product.molecule.sequence.len()
                {
                    return Err(bad("historical source segment outside documentary bounds"));
                }
            }
            for junction in &product.junctions {
                if junction.output_boundary > product.molecule.sequence.len()
                    || junction.product_overlap_positions.len() != 7
                    || junction
                        .product_overlap_positions
                        .iter()
                        .any(|p| *p >= product.molecule.sequence.len())
                    || junction.matched_output_variants.len() > 64
                    || junction
                        .matched_output_variants
                        .iter()
                        .any(|id| id.len() > 256)
                {
                    return Err(Error::TooLarge);
                }
                for hit in &junction.parents {
                    hit_shape(hit, count)?;
                }
            }
            items = items
                .checked_add(product.segments.len() + product.junctions.len())
                .filter(|n| *n <= 100_000)
                .ok_or(Error::TooLarge)?;
        }
    }
    let mut annotated_bases = 0usize;
    let mut features = 0usize;
    let mut unmapped = 0usize;
    // Count every nested array before any model recursion, rich JSON hashing or clone.
    for document in &batch.documents {
        if document.engine_version.is_empty()
            || document.engine_version.len() > 128
            || document.molecule.name.len() > 1024
            || document.molecule.sequence.len() > 2_000_000
            || document.molecule.metadata.len() > 16
            || document.molecule.original_genbank.is_some()
            || document.annotation_evidence.len() != count
            || document.reciprocal_components.len() > 32
        {
            return Err(Error::TooLarge);
        }
        annotated_bases = annotated_bases
            .checked_add(document.molecule.sequence.len())
            .filter(|n| *n <= 8_000_000)
            .ok_or(Error::TooLarge)?;
        features = features
            .checked_add(document.molecule.features.len())
            .filter(|n| *n <= 32768)
            .ok_or(Error::TooLarge)?;
        if document.selection.batch_revision != batch.batch_revision
            || !hex(&document.reference_product_revision)
            || document.reciprocal_components.iter().any(|component| {
                !hex(&component.reference_revision) || component.length > 2_000_000
            })
        {
            return Err(bad("historical reference product identity is invalid"));
        }
        texts(&document.warnings)?;
        for (key, value) in &document.molecule.metadata {
            if key.len() > 128 || value.len() > 256 * 1024 {
                return Err(Error::TooLarge);
            }
        }
        for evidence in &document.annotation_evidence {
            if evidence.input >= count
                || evidence.source.role.len() > 256
                || !hex(&evidence.source.revision)
                || evidence.segments.len() > 2048
                || evidence.co_sourced_segments.len() > 256
                || evidence.mapped_features > 32768
                || evidence.reviewed_features > evidence.mapped_features
                || evidence.dropped_features > 32768
            {
                return Err(Error::TooLarge);
            }
            unmapped = unmapped
                .checked_add(evidence.unmapped_features.len())
                .filter(|n| *n <= 32768)
                .ok_or(Error::TooLarge)?;
            texts(&evidence.warnings)?;
            for map in evidence
                .segments
                .iter()
                .chain(&evidence.co_sourced_segments)
            {
                if map.source_start >= map.source_end
                    || map.source_end > 2_000_000
                    || map
                        .destination_start
                        .checked_add(map.source_end - map.source_start)
                        .is_none_or(|end| end > document.molecule.sequence.len())
                {
                    return Err(bad("historical annotation map outside documentary bounds"));
                }
            }
        }
    }
    // Iterative shared model guards reject deep/fanout/qualifier payloads before
    // serialization. Each document's duplicated evidence retains its own budget.
    for document in &batch.documents {
        crate::cloning::recombination::reference::model_preflight(&[&document.molecule])?;
        crate::cloning::recombination::reference::feature_preflight(
            &document
                .annotation_evidence
                .iter()
                .map(|e| e.unmapped_features.as_slice())
                .collect::<Vec<_>>(),
        )?;
        document.molecule.validate()?;
    }
    Ok(())
}
struct Counter(usize);
impl Write for Counter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| std::io::Error::other("reference evidence exceeds32MiB"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn validate(e: &Evidence, r: &RunRecord) -> Result<(), Error> {
    if e.schema_version != 1
        || e.run_id != r.run_id
        || e.plan_revision != r.plan_revision
        || e.engine_version != r.assessment.engine_version
        || e.steps.is_empty()
        || e.steps.len() > 64
        || e.inputs.len() > 1000
        || e.products.len() > 1000
    {
        return Err(bad("run identity, engine or exact event bounds disagree"));
    }
    for step in &e.steps {
        let crate::verification::evidence::ReferenceStepEvidence::MultiSiteRecombination {
            batch,
            ..
        } = &step.evidence;
        batch_shape(batch)?;
    }
    // All nested model/location shapes have been bounded above. Count the whole
    // attachment before rich hashes or selected-document cloning.
    serde_json::to_writer(Counter(0), e).map_err(|_| Error::TooLarge)?;
    if assembly_recipe::hash(&(&e.target, &e.inputs, &e.products))?
        != assembly_recipe::hash(&(&r.target, &r.inputs, &r.products))?
    {
        return Err(bad("exact source/target/product events differ"));
    }
    let expected = r
        .plan
        .steps
        .iter()
        .enumerate()
        .filter_map(|(index, step)| {
            matches!(step, Operation::MultiSiteRecombineReference { .. }).then_some(index)
        })
        .collect::<Vec<_>>();
    if expected.len() != e.steps.len() {
        return Err(bad("reference step evidence is missing or duplicated"));
    }
    for (index, step) in expected.into_iter().zip(&e.steps) {
        let crate::verification::evidence::ReferenceStepEvidence::MultiSiteRecombination {
            selection,
            batch,
        } = &step.evidence;
        let Operation::MultiSiteRecombineReference {
            sources,
            request,
            selection: declared,
        } = &r.plan.steps[index]
        else {
            unreachable!()
        };
        if step.step_index != index
            || selection.batch_revision != declared.batch_revision
            || selection.candidate_id != declared.candidate_id
            || selection.component_id != declared.component_id
            || batch.batch_revision != selection.batch_revision
            || assembly_recipe::hash(request)? != assembly_recipe::hash(&batch.request)?
            || sources.len() != batch.sources.len()
        {
            return Err(bad("reference step declaration or selection changed"));
        }
        for (source, identity) in sources.iter().zip(&batch.sources) {
            let revision = match &source.source {
                Source::Input(id) => r
                    .plan
                    .inputs
                    .iter()
                    .find(|input| &input.id == id)
                    .map(|input| input.molecule.revision())
                    .ok_or_else(|| bad("missing original source input"))?,
                Source::Step(prior) if *prior < index => r
                    .products
                    .iter()
                    .find(|product| product.step_index == *prior)
                    .map(|product| product.revision.clone())
                    .ok_or_else(|| bad("missing exact upstream product event"))?,
                Source::Step(_) => return Err(bad("cyclic or forward source reference")),
            };
            if source.role != identity.role || revision != identity.revision {
                return Err(bad("source role/revision differs from exact run source"));
            }
        }
        let document = batch
            .documents
            .iter()
            .find(|document| {
                document.selection.batch_revision == selection.batch_revision
                    && document.selection.candidate_id == selection.candidate_id
                    && document.selection.component_id == selection.component_id
            })
            .ok_or_else(|| bad("selected reference document absent from history"))?;
        let product = r
            .products
            .iter()
            .find(|product| product.step_index == index)
            .ok_or_else(|| bad("selected reference step has no retained event"))?;
        let mut recorded = document.molecule.clone();
        recorded
            .metadata
            .insert("plan_revision".into(), r.plan_revision.clone());
        recorded
            .metadata
            .insert("engine_version".into(), r.assessment.engine_version.clone());
        if recorded.revision() != product.revision {
            return Err(bad(
                "selected reference document differs from retained product event",
            ));
        }
    }
    Ok(())
}

pub(crate) fn attachment(
    r: &RunRecord,
    steps: Vec<ReferenceStepReport>,
) -> Result<Option<WorkspaceAttachment>, Error> {
    if steps.is_empty() {
        return Ok(None);
    }
    let e = Evidence {
        schema_version: 1,
        run_id: r.run_id.clone(),
        plan_revision: r.plan_revision.clone(),
        engine_version: r.assessment.engine_version.clone(),
        target: r.target.clone(),
        inputs: r.inputs.clone(),
        products: r.products.clone(),
        steps,
    };
    validate(&e, r)?;
    Ok(Some(WorkspaceAttachment {
        kind: KIND.into(),
        name: format!("{}.reference-design-graphs.json", r.run_id),
        media_type: "application/json".into(),
        bytes: encode(&e)?,
    }))
}
/// Structural archive validation only. No enzyme/PCR/contour calculation runs.
pub(crate) fn validate_workspace(w: &Workspace) -> Result<(), Error> {
    let mut records = BTreeMap::new();
    for a in w
        .attachments
        .values()
        .filter(|a| a.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND)
    {
        let r = run_lineage::load_run_record(&a.bytes)?;
        records.insert(r.run_id.clone(), r);
    }
    let mut found = BTreeSet::new();
    for a in w.attachments.values().filter(|a| a.kind == KIND) {
        if a.media_type != "application/json" {
            return Err(bad("attachment media type changed"));
        }
        let e = decode(&a.bytes)?;
        if !found.insert(e.run_id.clone()) {
            return Err(bad("duplicate reference graph attachment"));
        }
        let r = records
            .get(&e.run_id)
            .ok_or_else(|| bad("graph attachment has no owning run"))?;
        validate(&e, r)?;
    }
    for r in records.values() {
        if r.plan
            .steps
            .iter()
            .any(|s| matches!(s, Operation::MultiSiteRecombineReference { .. }))
            && !found.contains(&r.run_id)
        {
            return Err(bad(
                "reference run is missing its historical graph attachment",
            ));
        }
    }
    Ok(())
}
/// Historical document only; callers must replay the literal plan for fresh state.
pub fn load(w: &Workspace, run_id: &str) -> Result<Option<Evidence>, Error> {
    w.validate()?;
    find(w, run_id)
}
pub(crate) fn find(w: &Workspace, run_id: &str) -> Result<Option<Evidence>, Error> {
    for a in w.attachments.values().filter(|a| a.kind == KIND) {
        let e = decode(&a.bytes)?;
        if e.run_id == run_id {
            return Ok(Some(e));
        }
    }
    Ok(None)
}
