//! Target-free digest inspection of freshly reconstructed overlap-extension PCR.
//! Imported summaries never restore physical state; saved origins are replayed.
use crate::{
    Error, Molecule, Topology,
    annotations::{self, MappingSegment},
    assembly_recipe,
    catalog_digest::{self, CatalogSelection},
    chemistry,
    modified_duplex::TrackedDuplex,
    overlap_extension::{self, Geometry},
    overlap_extension_pcr,
    pcr_digest_inspection::GeometrySummary,
    restriction_inspection::{self, InspectionOptions, RestrictionInspection},
    run_lineage::EventBinding,
    run_material::{Capsule, Selection},
    verification::{self, PcrDocument, Status},
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const ENGINE_VERSION: &str = "overlap-extension-digest-inspection-1";
const MAX_GEOMETRIES: usize = 10_000;
const MAX_BYTES: usize = 32 << 20;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request {pub source_role:String,pub preparation:overlap_extension_pcr::Request,pub selection:CatalogSelection,pub options:InspectionOptions}}
dto! {pub struct PcrStage {pub stage:String,pub template_revision:String,pub geometries:Vec<GeometrySummary>,pub product_revision:Option<String>,pub product_length:Option<usize>}}
dto! {pub struct OriginSummary {pub selection:Selection,pub target:EventBinding,pub plan_revision:String,pub engine_version:String,pub status:Status,pub material_revision:String}}
dto! {pub struct Preview {pub engine_version:String,pub context_revision:String,pub source_revision:String,pub origin_revision:Option<String>,pub origin_status:Option<Status>,pub request:Request,pub pcr_stages:Vec<PcrStage>,pub extension:Option<Geometry>,pub origin_evidence:Option<OriginSummary>,pub amplicon:Option<Molecule>,pub inspection:Option<RestrictionInspection>,pub incomplete:bool,pub findings:Vec<String>}}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("overlap-extension digest inspection: {message}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
struct Count(usize);
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| std::io::Error::other("bounded extension digest evidence"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn size(value: &impl Serialize) -> Result<(), Error> {
    serde_json::to_writer(Count(0), value).map_err(|_| Error::TooLarge)
}
fn shape(request: &Request) -> Result<(), Error> {
    if request.source_role.trim().is_empty()
        || request.source_role.len() > 128
        || request.source_role.chars().any(char::is_control)
    {
        return Err(bad("invalid bounded source role"));
    }
    request.preparation.validate_shape()?;
    request.selection.validate()?;
    request.options.validate()
}
pub fn source_preflight(
    source: &Molecule,
    origin: Option<&Capsule>,
    request: &Request,
) -> Result<(), Error> {
    shape(request)?;
    crate::model::bounds::model_preflight(&[source])?;
    size(source)?;
    if let Some(origin) = origin {
        origin.validate_documentary()?;
        if source.revision() != origin.selection.event.revision {
            return Err(bad(
                "saved source differs from the exact pinned origin event",
            ));
        }
    } else if chemistry::has_unmapped_physical_evidence(source)
        || source
            .metadata
            .get(crate::recombination::REFERENCE_ONLY_KEY)
            .is_some()
    {
        return Err(bad(
            "imported physical/reference-only material requires reconstruction of its originating graph",
        ));
    }
    // The borrowed guards ran before these bounded request-field copies. There
    // is no intended target and no copied/sliced original declared as new stock.
    let row = assembly_recipe::FragmentSpec {
        generation: None,
        role: request.source_role.clone(),
        input_id: request.source_role.clone(),
        preparation: assembly_recipe::Preparation::OverlapExtensionPcrCatalogDigest {
            request: request.preparation.clone(),
            selection: request.selection.clone(),
            fragment_index: 0,
            repair: None,
        },
        reverse: false,
    };
    assembly_recipe::preparation_rows_preflight(
        &[(&request.source_role, source)],
        &request.source_role,
        &[row],
    )?;
    let catalog = catalog_digest::builtin_catalog()?;
    catalog_digest::source_preflight(
        &catalog,
        &request.selection,
        &request.source_role,
        source,
        &request.options,
    )
}
fn context(
    source_revision: &str,
    origin_revision: &Option<String>,
    request: &Request,
) -> Result<String, Error> {
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        restriction_inspection::ENGINE_VERSION,
        source_revision,
        origin_revision,
        request,
    ))
}
pub fn inspection_revision(
    source: &Molecule,
    origin: Option<&Capsule>,
    request: &Request,
) -> Result<String, Error> {
    source_preflight(source, origin, request)?;
    context(
        &source.revision(),
        &origin.map(assembly_recipe::hash).transpose()?,
        request,
    )
}
fn pcr_stage(
    stage: &str,
    template: &Molecule,
    result: Result<Vec<PcrDocument>, Error>,
    findings: &mut Vec<String>,
) -> Result<(PcrStage, Option<PcrDocument>), Error> {
    let documents = match result {
        Ok(documents) => documents,
        Err(error) => {
            findings.push(format!("{stage} PCR could not complete: {error}"));
            return Ok((
                PcrStage {
                    stage: stage.into(),
                    template_revision: template.revision(),
                    geometries: vec![],
                    product_revision: None,
                    product_length: None,
                },
                None,
            ));
        }
    };
    if documents.len() > MAX_GEOMETRIES {
        return Err(Error::TooLarge);
    }
    let geometries = documents
        .iter()
        .enumerate()
        .map(|(index, p)| GeometrySummary {
            index,
            forward_start: p.geometry.forward_start,
            reverse_start: p.geometry.reverse_start,
            product_length: p.molecule.sequence.len(),
        })
        .collect();
    let product = if documents.len() == 1 {
        documents.into_iter().next()
    } else {
        findings.push(format!("{stage} PCR retained{} geometries; downstream material is withheld rather than selecting a sequence matching an intended product.",documents.len()));
        None
    };
    let summary = PcrStage {
        stage: stage.into(),
        template_revision: template.revision(),
        geometries,
        product_revision: product.as_ref().map(|p| p.molecule.revision()),
        product_length: product.as_ref().map(|p| p.molecule.sequence.len()),
    };
    if let Some(p) = &product {
        findings.extend(p.findings.iter().take(100).cloned());
    }
    Ok((summary, product))
}
fn project(
    source: &Molecule,
    offset: usize,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<(), Error> {
    if source.features.is_empty() {
        return Ok(());
    }
    let projected = annotations::project_features_with_topology(
        source,
        &[MappingSegment {
            source_start: 0,
            source_end: source.sequence.len(),
            destination_start: offset,
            reverse: false,
        }],
        product.sequence.len(),
        product.topology,
    )?;
    product.features.extend(projected.features);
    if !projected.unmapped_features.is_empty() {
        product.metadata.insert(
            format!("annotation_mapping:unmapped:{}", source.revision()),
            serde_json::to_string(&projected.unmapped_features).map_err(|_| Error::TooLarge)?,
        );
    }
    findings.extend(
        projected
            .warnings
            .iter()
            .take(100)
            .map(|s| s.chars().take(2048).collect()),
    );
    let mut seen = BTreeSet::new();
    let mut features = Vec::with_capacity(product.features.len());
    for feature in product.features.drain(..) {
        if seen.insert(assembly_recipe::hash(&feature)?) {
            features.push(feature);
        }
    }
    product.features = features;
    Ok(())
}
pub fn inspect(
    source: &Molecule,
    origin: Option<&Capsule>,
    request: &Request,
) -> Result<Preview, Error> {
    let context_revision = inspection_revision(source, origin, request)?;
    let mut findings=vec!["Fresh ideal exact-PCR branches, selected reciprocal extension and outer PCR, followed by the pinned assumed-complete digest. There is no intended target input. Geometry uniqueness is conditional on this model and is not experimental isolation or reaction success.".into()];
    let material = origin
        .map(|origin| crate::run_material_capsule::replay(origin, source))
        .transpose()?;
    let origin_evidence = if let Some(material) = &material {
        let report = material.evidence.decode_report()?;
        Some(OriginSummary {
            selection: material.evidence.origin.selection.clone(),
            target: material.evidence.origin.target_event.clone(),
            plan_revision: report.plan_revision,
            engine_version: report.engine_version,
            status: report.status,
            material_revision: material.evidence.material_revision.clone(),
        })
    } else {
        None
    };
    let origin_status = origin_evidence.as_ref().map(|o| o.status);
    let mut incomplete = material.as_ref().is_some_and(|m| m.incomplete);
    if origin_status.is_some() {
        findings.push("Source chemistry was freshly reconstructed against the complete unchanged original saved scope. Origin target/event identity is retained; its incomplete findings remain applicable.".into());
    }
    let initial = |pair: &overlap_extension_pcr::PrimerPair| {
        if let Some(material) = &material {
            verification::simulate_pcr_documents_from_fresh(
                source,
                &pair.forward,
                &pair.reverse,
                &material.tracked,
            )
        } else {
            verification::simulate_pcr_documents(source, &pair.forward, &pair.reverse)
        }
    };
    let (left_stage, left) = pcr_stage(
        "left",
        source,
        initial(&request.preparation.left),
        &mut findings,
    )?;
    let (right_stage, right) = pcr_stage(
        "right",
        source,
        initial(&request.preparation.right),
        &mut findings,
    )?;
    let mut pcr_stages = vec![left_stage, right_stage];
    let mut extension = None;
    let mut amplicon = None;
    let mut inspection = None;
    if let (Some(left), Some(right)) = (left, right) {
        let a = TrackedDuplex::exact_pcr(&left.molecule.sequence)?;
        let b = TrackedDuplex::exact_pcr(&right.molecule.sequence)?;
        match a.overlap_extend(&b, &request.preparation.extension) {
            Ok((extended, geometry)) => {
                incomplete |=
                    !geometry.terminal_pairings.complete || !geometry.terminal_pairings.unique;
                findings.extend(geometry.warnings.iter().take(100).cloned());
                let mut reference = extended
                    .duplex
                    .molecule("Fresh reciprocal overlap-extension substrate")?;
                project(&left.molecule, 0, &mut reference, &mut findings)?;
                project(
                    &right.molecule,
                    a.duplex.top.len() - request.preparation.extension.overlap,
                    &mut reference,
                    &mut findings,
                )?;
                let outer = &request.preparation.outer;
                let (outer_stage, product) = pcr_stage(
                    "outer",
                    &reference,
                    verification::simulate_pcr_documents_from_fresh(
                        &reference,
                        &outer.forward,
                        &outer.reverse,
                        &extended,
                    ),
                    &mut findings,
                )?;
                pcr_stages.push(outer_stage);
                if let Some(product) = product {
                    if incomplete {
                        findings.push("The selected actual digest is conditional: originating/terminal-pair uncertainty remains incomplete. Selecting a fragment never resolves alternative pairings or upgrades whole-plan status.".into());
                    }
                    let physical = TrackedDuplex::exact_pcr(&product.molecule.sequence)?;
                    let catalog = catalog_digest::builtin_catalog()?;
                    let enzymes = request.selection.resolved_enzymes(&catalog)?;
                    let result = restriction_inspection::inspect_fresh_duplex(
                        &request.source_role,
                        &product.molecule,
                        &enzymes,
                        &request.options,
                        &physical.duplex,
                    )?;
                    incomplete |= !result.scan_complete || !result.products_determined;
                    inspection = Some(result);
                    amplicon = Some(product.molecule);
                } else {
                    incomplete = true;
                }
                extension = Some(geometry);
            }
            Err(error) => {
                incomplete = true;
                findings.push(format!("Reciprocal extension unavailable: {error}"));
            }
        }
    } else {
        incomplete = true;
    }
    let preview = Preview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: source.revision(),
        origin_revision: origin.map(assembly_recipe::hash).transpose()?,
        origin_status,
        request: request.clone(),
        pcr_stages,
        extension,
        origin_evidence,
        amplicon,
        inspection,
        incomplete,
        findings,
    };
    preview.validate()?;
    Ok(preview)
}
impl Preview {
    /// Bounded documentary correspondence only. No PCR, source-origin replay,
    /// polymerase extension or digest is performed by this validator.
    pub fn validate(&self) -> Result<(), Error> {
        shape(&self.request)?;
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.context_revision)
            || !hex(&self.source_revision)
            || self.context_revision
                != context(&self.source_revision, &self.origin_revision, &self.request)?
            || self.pcr_stages.len() < 2
            || self.pcr_stages.len() > 3
            || self.findings.len() > 1000
            || self.findings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad(
                "invalid preview context, source or bounded stage evidence",
            ));
        }
        for (i, p) in self.pcr_stages.iter().enumerate() {
            if p.stage != ["left", "right", "outer"][i]
                || !hex(&p.template_revision)
                || i < 2 && p.template_revision != self.source_revision
                || p.geometries.len() > MAX_GEOMETRIES
            {
                return Err(bad("invalid native PCR stage binding"));
            }
            for (index, g) in p.geometries.iter().enumerate() {
                if g.index != index
                    || g.product_length == 0
                    || g.product_length > crate::model::MAX_SEQUENCE_LEN
                    || g.forward_start > crate::model::MAX_SEQUENCE_LEN
                    || g.reverse_start > crate::model::MAX_SEQUENCE_LEN
                {
                    return Err(bad("invalid native PCR geometry summary"));
                }
            }
            match (&p.product_revision, p.product_length) {
                (Some(revision), Some(length))
                    if p.geometries.len() == 1
                        && hex(revision)
                        && length == p.geometries[0].product_length => {}
                (None, None) if p.geometries.len() != 1 => {}
                _ => {
                    return Err(bad(
                        "only a native singleton PCR may expose a selected stage product",
                    ));
                }
            }
        }
        size(self)?;
        match (
            &self.origin_revision,
            self.origin_status,
            &self.origin_evidence,
        ) {
            (None, None, None) => {}
            (Some(revision), Some(status), Some(e))
                if hex(revision)
                    && status != Status::Invalid
                    && status == e.status
                    && hex(&e.plan_revision)
                    && hex(&e.material_revision)
                    && hex(&e.target.revision)
                    && e.selection.event.revision == self.source_revision
                    && e.engine_version == verification::ENGINE_VERSION
                    && e.selection.run_id.len() == 64
                    && hex(&e.selection.run_id)
                    && e.selection.step_index < 1000
                    && [
                        &e.selection.event.document_id,
                        &e.selection.event.event_id,
                        &e.target.document_id,
                        &e.target.event_id,
                    ]
                    .iter()
                    .all(|s| !s.is_empty() && s.len() <= 256) => {}
            _ => return Err(bad("saved-origin documentary binding changed")),
        }
        if let Some(g) = &self.extension {
            if g.engine_version != overlap_extension::ENGINE_VERSION
                || g.options != self.request.preparation.extension
                || g.homology.len() != g.options.overlap
                || g.mappings.len() > 8
                || g.fill.mappings.len() > 8
                || g.source_geometry_revisions.iter().any(|s| !hex(s))
                || g.terminal_pairings.overlaps.len() > g.options.max_candidates
                || g.terminal_pairings.comparisons > g.options.max_comparisons
                || g.terminal_pairings.unique
                    && (!g.terminal_pairings.complete
                        || g.terminal_pairings.first_alternative.is_some())
                || g.fill.source != g.annealed
                || g.fill.method != chemistry::EndRepairMethod::FillFivePrimeOverhangs
                || self.pcr_stages[0].product_length != Some(g.annealed.top.len())
                || self.pcr_stages[1].product_length != Some(g.annealed.bottom.len())
                || self.pcr_stages.len() != 3
            {
                return Err(bad(
                    "extension evidence differs from actual singleton branch lengths/options",
                ));
            }
        } else if self.pcr_stages.len() != 2 || self.amplicon.is_some() || self.inspection.is_some()
        {
            return Err(bad(
                "downstream products require retained native extension evidence",
            ));
        }
        if let Some(m) = &self.amplicon {
            crate::model::bounds::model_preflight(&[m])?;
            m.validate()?;
            let outer = &self.pcr_stages[2];
            let amplicon_revision = m.revision();
            if m.topology != Topology::Linear
                || outer.product_length != Some(m.sequence.len())
                || outer.product_revision.as_deref() != Some(amplicon_revision.as_str())
                || m.metadata.get("pcr:template_revision") != Some(&outer.template_revision)
                || m.metadata.get("pcr:engine_version").map(String::as_str)
                    != Some(verification::ENGINE_VERSION)
                || m.metadata.get("pcr:primers")
                    != Some(
                        &serde_json::to_string(&(
                            &self.request.preparation.outer.forward,
                            &self.request.preparation.outer.reverse,
                        ))
                        .map_err(|_| Error::TooLarge)?,
                    )
            {
                return Err(bad("outer amplicon differs from its native PCR stage"));
            }
        } else if self.inspection.is_some()
            || self
                .pcr_stages
                .get(2)
                .is_some_and(|p| p.product_revision.is_some())
        {
            return Err(bad("digest fragments require a native singleton outer PCR"));
        }
        let pipeline_available = self.extension.is_some() && self.amplicon.is_some();
        let inherited_incomplete = self.origin_status == Some(Status::Incomplete)
            || self
                .extension
                .as_ref()
                .is_none_or(|g| !g.terminal_pairings.complete || !g.terminal_pairings.unique);
        if let Some(report) = &self.inspection {
            if !pipeline_available {
                return Err(bad(
                    "digest evidence requires actual singleton PCRs and the explicit selected extension",
                ));
            }
            let catalog = catalog_digest::builtin_catalog()?;
            let enzymes = self.request.selection.resolved_enzymes(&catalog)?;
            report.validate_bindings_for(
                &self.request.source_role,
                self.amplicon
                    .as_ref()
                    .ok_or_else(|| bad("missing amplicon"))?,
                &enzymes,
                &self.request.options,
            )?;
        }
        let incomplete = inherited_incomplete
            || !pipeline_available
            || self
                .inspection
                .as_ref()
                .is_none_or(|r| !r.scan_complete || !r.products_determined);
        if self.incomplete != incomplete {
            return Err(bad("inherited or capped incompleteness was changed"));
        }
        size(self)
    }
    pub fn validate_binding(&self, revision: &str) -> Result<(), Error> {
        if self.context_revision != revision {
            return Err(bad("stale source/preparation/catalog/origin context"));
        }
        self.validate()
    }
    pub fn selectable(&self) -> bool {
        self.amplicon.is_some()
            && self.inspection.as_ref().is_some_and(|r| {
                r.scan_complete && r.products_determined && !r.fragments.is_empty()
            })
    }
}
