//! Source-backed processed TOPO material. A fresh paired-contour proof can bridge
//! only chemistry represented faithfully by the ordinary two-strand model.
use crate::{
    Error, Molecule, Topology, assembly_recipe, chemistry, modified_duplex as tracked, topo,
    topo_annotations, topo_nick_sealing as treatment, topo_worker,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use topo::physical::{self, processed_projection as processed, projection};
pub const ENGINE_VERSION: &str = "topo-processed-source-product-1";
pub const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request {pub treatment:treatment::Request,pub selection:processed::Selection,pub options:projection::Options}}
dto! {pub struct Report {
    pub engine_version:String,pub context_revision:String,pub preparation_revision:String,
    pub request:Request,pub sources:Vec<topo::SourceIdentity>,pub preparation:topo::PreparationEvidence,
    pub projection:processed::ProcessedReport,pub duplex:Option<chemistry::Duplex>,
    pub tracked:Option<tracked::TrackedDuplex>,pub product:Option<Molecule>,
    pub annotations:Option<topo_annotations::AnnotationEvidence>,pub bridge_notes:Vec<String>,
}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("Processed TOPO source product: {s}"))
}
fn physical_error(e: impl std::fmt::Display) -> Error {
    bad(&e.to_string())
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn shape(r: &Request) -> Result<(), Error> {
    let original = &r.treatment.preparation;
    physical::nick_sealing::request_preflight(&r.treatment.sealing, &original.engine_limits)
        .map_err(physical_error)?;
    if original.recipe.schema_version != 1
        || original.recipe.materials.is_empty()
        || original.recipe.materials.len() > 16
        || original
            .transfers
            .len()
            .checked_add(r.treatment.sealing.nicks.len())
            .is_none_or(|n| n > 16 || n > original.engine_limits.max_steps)
    {
        return Err(Error::TooLarge);
    }
    let l = &r.options.limits;
    if !digest(&r.selection.sealing_context_revision)
        || r.selection.component >= 64
        || r.selection.reference_seed.strand >= 32
        || r.selection.reference_seed.index >= 100_000
        || l.max_pairs == 0
        || l.max_pairs > 500_000
        || l.max_work_cells == 0
        || l.max_work_cells > 128_000_000
        || l.max_work_bytes == 0
        || l.max_work_bytes > 128 * 1024 * 1024
        || l.max_mapping_segments == 0
        || l.max_mapping_segments > 16384
        || l.max_report_bytes == 0
        || l.max_report_bytes > 32 * 1024 * 1024
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
/// Borrowed guards before source hashes, preparation/PCR, source cloning or proof.
pub fn preflight(inputs: &[topo::Input], r: &Request) -> Result<(), Error> {
    shape(r)?;
    treatment::preflight(inputs, &r.treatment)
}
fn binding(
    preparation_revision: &str,
    sources: &[topo::SourceIdentity],
    r: &Request,
) -> Result<String, Error> {
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        processed::ENGINE,
        preparation_revision,
        sources,
        r,
    ))
}
pub fn revision(inputs: &[topo::Input], r: &Request) -> Result<String, Error> {
    preflight(inputs, r)?;
    binding(
        &topo::batch_revision(inputs, &r.treatment.preparation)?,
        &topo::identities(inputs),
        r,
    )
}
fn five(value: physical::FivePrime) -> Result<chemistry::FivePrime, String> {
    Ok(match value {
        physical::FivePrime::Phosphate=>chemistry::FivePrime::Phosphate,
        physical::FivePrime::Hydroxyl=>chemistry::FivePrime::Hydroxyl,
        physical::FivePrime::Unknown=>chemistry::FivePrime::Unknown,
        physical::FivePrime::Blocked=>return Err("Blocked 5′ terminal chemistry is not represented by the ordinary duplex model; no OH/P substitution was made.".into()),
    })
}
fn three(value: &physical::ThreePrime) -> Result<chemistry::ThreePrime, String> {
    Ok(match value {
        physical::ThreePrime::Hydroxyl=>chemistry::ThreePrime::Hydroxyl,
        physical::ThreePrime::Unknown=>chemistry::ThreePrime::Unknown,
        physical::ThreePrime::Blocked=>chemistry::ThreePrime::Blocked,
        physical::ThreePrime::Phosphate=>return Err("3′ phosphate chemistry cannot be collapsed to generic blocked/OH chemistry in the ordinary duplex model.".into()),
        physical::ThreePrime::Phosphotyrosyl{..}=>return Err("A retained phosphotyrosyl adduct cannot be represented faithfully by the ordinary duplex model.".into()),
    })
}
fn terminus(e: &projection::EndChemistry) -> Result<chemistry::Terminus, String> {
    Ok(chemistry::Terminus {
        five_prime: five(e.five_prime)?,
        three_prime: three(&e.three_prime)?,
    })
}
fn geometry(
    view: &projection::DuplexView,
) -> Result<
    (
        Topology,
        chemistry::Terminus,
        chemistry::Terminus,
        Vec<chemistry::Nick>,
    ),
    String,
> {
    if view.top.len() != view.bottom.len() || view.top.is_empty() {
        return Err("Only a complete aligned paired contour is bridgeable.".into());
    }
    let topology = match view.topology {
        projection::ContourTopology::Linear => Topology::Linear,
        projection::ContourTopology::Circular => Topology::Circular,
    };
    let (left, right) = match (topology, &view.ends) {
        (Topology::Linear, Some(ends)) => (terminus(&ends.left)?, terminus(&ends.right)?),
        (Topology::Circular, None) => (chemistry::Terminus::UNKNOWN, chemistry::Terminus::UNKNOWN),
        _ => return Err("Contour topology and exposed terminal evidence disagree.".into()),
    };
    let mut nicks = Vec::with_capacity(view.nicks.len());
    for n in &view.nicks {
        if n.three_prime != physical::ThreePrime::Hydroxyl
            || n.five_prime != physical::FivePrime::Phosphate
        {
            return Err("An internal nick is bridgeable only with explicit 3′OH/5′P chemistry. Unknown, blocked, 5′OH, phosphate or adduct states remain in the retained graph.".into());
        }
        nicks.push(chemistry::Nick {
            strand: match n.strand {
                projection::ReferenceStrand::Top => chemistry::Strand::Top,
                projection::ReferenceStrand::Bottom => chemistry::Strand::Bottom,
            },
            position: i64::try_from(n.boundary)
                .map_err(|_| "Nick coordinate exceeds the duplex model".to_string())?,
        });
    }
    Ok((topology, left, right, nicks))
}
fn bridge(
    view: &projection::DuplexView,
) -> Result<(chemistry::Duplex, tracked::TrackedDuplex), String> {
    let (topology, left, right, nicks) = geometry(view)?;
    let duplex = chemistry::Duplex {
        top: view.top.clone(),
        bottom: view.bottom.clone(),
        bottom_start: 0,
        topology,
        left,
        right,
        nicks,
    };
    duplex.validate().map_err(|e| e.to_string())?;
    let tracked =
        tracked::TrackedDuplex::from_unknown_duplex(&duplex).map_err(|e| e.to_string())?;
    Ok((duplex, tracked))
}
/// Fresh original preparation, one transfer/sealing reconstruction, then contour
/// proof. No imported physical graph, reference summary or target selects DNA.
pub fn inspect(inputs: &[topo::Input], r: &Request) -> Result<Report, Error> {
    preflight(inputs, r)?;
    let original = &r.treatment.preparation;
    let preparation_revision = topo::batch_revision(inputs, original)?;
    let sources = topo::identities(inputs);
    let context_revision = binding(&preparation_revision, &sources, r)?;
    let preparation = topo::prepare(inputs, original)?;
    let projection = processed::project(
        &preparation.substrates,
        &original.transfers,
        &original.engine_limits,
        &r.treatment.sealing,
        &r.selection,
        &r.options,
    )
    .map_err(physical_error)?;
    let (mut duplex, mut tracked, mut product, mut annotations) = (None, None, None, None);
    let mut bridge_notes=vec!["This source-backed graph bridge preserves unknown internal modification calls. Reference features describe correspondence, not inherited atoms, methylation or experimentally established material.".into()];
    bridge_notes.push("The chosen component/orientation is an explicit graph view, not physical isolation or a reaction-mixture winner. Every other component, unpaired call and released activation remains in the retained processing evidence.".into());
    match &projection.outcome {
        projection::Outcome::Withheld(reason)=>bridge_notes.push(format!("Paired-contour proof withheld: {reason:?}; no ordinary duplex material or sequence document was created.")),
        projection::Outcome::Proven(view)=>match bridge(view){
            Err(reason)=>bridge_notes.push(reason),
            Ok((d,t))=>{
                let mut m=d.molecule("Processed TOPO product")?;
                m.metadata.insert("duplex_model".into(),ENGINE_VERSION.into());
                m.metadata.insert("topo:processed_context".into(),context_revision.clone());
                m.metadata.insert("topo:physical_evidence_required".into(),"true".into());
                match topo_annotations::project(inputs,&preparation,view,&m){
                    Ok(evidence)=>{m.features=evidence.features.clone();annotations=Some(evidence);},
                    Err(error)=>bridge_notes.push(format!("Source annotations withheld: {error}. Actual duplex geometry remains retained; no substitute features were invented.")),
                }
                m.validate()?;duplex=Some(d);tracked=Some(t);product=Some(m);
                bridge_notes.push("The standalone product document is reference-only: it retains an unmapped-physical-evidence marker and cannot establish ordinary material by import. Circular placeholder termini are unused; actual end and nick states remain in the graph.".into());
            }
        }
    }
    let report = Report {
        engine_version: ENGINE_VERSION.into(),
        context_revision: context_revision.clone(),
        preparation_revision,
        request: r.clone(),
        sources,
        preparation,
        projection,
        duplex,
        tracked,
        product,
        annotations,
        bridge_notes,
    };
    report.validate_binding(&context_revision)?;
    Ok(report)
}
impl Report {
    /// Bounded readonly binding and representability, not chemistry replay.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        shape(&self.request)?;
        if self.engine_version != ENGINE_VERSION
            || !digest(expected)
            || expected != self.context_revision
            || !digest(&self.preparation_revision)
            || self.sources.is_empty()
            || self.sources.len() > 32
            || self.bridge_notes.len() > 16
            || self.bridge_notes.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("report source/request bounds changed"));
        }
        bounded(&self.request, 2 * 1024 * 1024)?;
        let mut roles = BTreeSet::new();
        for s in &self.sources {
            if s.role.trim().is_empty()
                || s.role.len() > 256
                || !digest(&s.revision)
                || !roles.insert(s.role.as_str())
            {
                return Err(bad("source identity changed"));
            }
        }
        let original = &self.request.treatment.preparation;
        let original_revision = assembly_recipe::hash(&(
            topo::ENGINE_VERSION,
            physical::ENGINE,
            physical::PRIMARY_MODEL,
            projection::ENGINE,
            &self.sources,
            original,
        ))?;
        if original_revision != self.preparation_revision
            || binding(&original_revision, &self.sources, &self.request)? != expected
        {
            return Err(bad("source/preparation/selection context changed"));
        }
        topo_worker::preparation(&self.preparation)?;
        processed::validate_report(
            &self.preparation.substrates,
            &original.transfers,
            &original.engine_limits,
            &self.request.treatment.sealing,
            &self.request.selection,
            &self.request.options,
            &self.projection,
        )
        .map_err(physical_error)?;
        match (
            &self.projection.outcome,
            &self.duplex,
            &self.tracked,
            &self.product,
        ) {
            (projection::Outcome::Proven(view), Some(d), Some(t), Some(m)) => {
                let (topology, left, right, nicks) = geometry(view).map_err(|e| bad(&e))?;
                if d.top != view.top
                    || d.bottom != view.bottom
                    || d.topology != topology
                    || d.bottom_start != 0
                    || d.left != left
                    || d.right != right
                    || d.nicks != nicks
                    || t.duplex != *d
                {
                    return Err(bad("bridged actual duplex geometry changed"));
                }
                if !t.generation.is_empty()
                    || [&t.top_modifications, &t.bottom_modifications]
                        .iter()
                        .any(|s| {
                            s.n6a_default != tracked::CallState::Unknown
                                || s.c5c_default != tracked::CallState::Unknown
                                || !s.spans.is_empty()
                        })
                {
                    return Err(bad("unknown modification calls were promoted"));
                }
                t.validate()?;
                collagene_model::bounds::model_preflight(&[m])?;
                m.validate()?;
                if m.sequence != d.top
                    || m.topology != d.topology
                    || m.metadata.get("duplex_model").map(String::as_str) != Some(ENGINE_VERSION)
                    || m.metadata.get("topo:processed_context") != Some(&self.context_revision)
                    || !chemistry::has_unmapped_physical_evidence(m)
                {
                    return Err(bad("standalone product lost physical-evidence restriction"));
                }
                if let Some(evidence) = &self.annotations {
                    topo_worker::annotation_sources_shape(
                        evidence,
                        &self.sources,
                        m.sequence.len(),
                    )?;
                    if assembly_recipe::hash(&evidence.features)?
                        != assembly_recipe::hash(&m.features)?
                    {
                        return Err(bad("product features differ from annotation evidence"));
                    }
                } else if !m.features.is_empty() {
                    return Err(bad(
                        "product has annotations without correspondence evidence",
                    ));
                }
            }
            (projection::Outcome::Proven(view), None, None, None) => {
                if geometry(view).is_ok() {
                    return Err(bad("representable contour bridge is unexpectedly absent"));
                }
                if self.annotations.is_some() {
                    return Err(bad("annotations without a bridgeable destination"));
                }
            }
            (projection::Outcome::Withheld(_), None, None, None) => {
                if self.annotations.is_some() {
                    return Err(bad("annotations without proven contour"));
                }
            }
            _ => return Err(bad("partial or invented bridge material")),
        }
        bounded(self, MAX_REPORT_BYTES)?;
        Ok(())
    }
}
struct Counter {
    bytes: usize,
    limit: usize,
}
impl std::io::Write for Counter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(b.len())
            .filter(|n| *n <= self.limit)
            .ok_or_else(|| std::io::Error::other("bounded processed TOPO evidence"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded(value: &impl Serialize, limit: usize) -> Result<(), Error> {
    serde_json::to_writer(Counter { bytes: 0, limit }, value).map_err(|_| Error::TooLarge)
}
pub fn replay(inputs: &[topo::Input], r: &Request, report: &Report) -> Result<(), Error> {
    report.validate_binding(&revision(inputs, r)?)?;
    if assembly_recipe::hash(&inspect(inputs, r)?)? != assembly_recipe::hash(report)? {
        return Err(bad("report differs from fresh source-backed processing"));
    }
    Ok(())
}
