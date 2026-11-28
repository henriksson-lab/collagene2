//! Source-backed processed TOPO material. A fresh paired-contour proof can bridge
//! only chemistry represented faithfully by the ordinary two-strand model.
use crate::cloning::topo::physical::{self, projection};
use crate::{
    chemistry::{self, tracked_duplex as tracked},
    cloning::topo::{self, annotations as topo_annotations, worker as topo_worker},
};
use collagene_model::{Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const ENGINE_VERSION: &str = "topo-source-actions-product-1";
pub const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request {pub preparation:topo::Request,pub selection:projection::Selection,pub options:projection::Options}}
dto! {pub struct Report {
    pub engine_version:String,pub context_revision:String,pub preparation_revision:String,
    pub request:Request,pub input_sources:Vec<topo::SourceIdentity>,pub sources:Vec<topo::SourceIdentity>,pub preparation:topo::PreparationEvidence,
    pub projection:projection::ProjectionReport,pub duplex:Option<chemistry::Duplex>,
    pub tracked:Option<tracked::TrackedDuplex>,pub product:Option<Molecule>,
    pub annotations:Option<topo_annotations::AnnotationEvidence>,pub bridge_notes:Vec<String>,
}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("TOPO source actions product: {s}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
/// Resource-only admission before plan hashes or source cloning. Selection pins,
/// immutable calls and source-role relationships are branch-local semantic checks.
pub fn resource_preflight(r: &Request) -> Result<(), Error> {
    physical::steps_resource_preflight(&r.preparation.steps, &r.preparation.engine_limits)
        .map_err(|error| Error::Unsupported(format!("TOPO transcript resources: {error}")))?;
    r.options.limits.validate().map_err(topo::projection_error)
}
fn shape(r: &Request) -> Result<(), Error> {
    resource_preflight(r)?;
    physical::steps_preflight(&r.preparation.steps, &r.preparation.engine_limits)
        .map_err(|error| Error::Unsupported(format!("TOPO transcript: {error}")))?;
    if !digest(&r.selection.source_batch_revision)
        || r.selection.component >= 64
        || r.selection.reference_seed.strand >= 32
        || r.selection.reference_seed.index >= 100_000
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
/// Borrowed request-only bounds; no source hashing, preparation or simulation.
pub fn request_preflight(r: &Request) -> Result<(), Error> {
    shape(r)
}
/// Borrowed guards before source hashes, preparation/PCR, source cloning or proof.
pub fn preflight(inputs: &[topo::Input], r: &Request) -> Result<(), Error> {
    shape(r)?;
    topo::source_preflight(inputs, &r.preparation)
}
fn binding(
    preparation_revision: &str,
    sources: &[topo::SourceIdentity],
    r: &Request,
) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        projection::ENGINE,
        preparation_revision,
        sources,
        r,
    ))
}
pub fn revision(inputs: &[topo::Input], r: &Request) -> Result<String, Error> {
    preflight(inputs, r)?;
    binding(
        &topo::batch_revision(inputs, &r.preparation)?,
        &topo::documentary_sources(inputs, &r.preparation)?,
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
pub(crate) fn geometry(
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
pub(crate) fn bridge(
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

/// Fresh original preparation, one transfer/sealing reconstruction, then contour
/// proof. No imported physical graph, reference summary or target selects DNA.
pub fn inspect(inputs: &[topo::Input], r: &Request) -> Result<Report, Error> {
    preflight(inputs, r)?;
    let preparation_revision = topo::batch_revision(inputs, &r.preparation)?;
    inspect_preflighted(inputs, r, preparation_revision, None, None).map(|(report, _)| report)
}
/// Fresh plan-only PCR synthesis route. Actual tracked Step templates establish
/// eligibility, not inherited atoms/modifications or imported receipt authority.
pub(crate) fn inspect_reconstructed(
    inputs: &[topo::Input],
    r: &Request,
    templates: &[(String, &tracked::TrackedDuplex)],
) -> Result<Report, Error> {
    inspect_reconstructed_native(inputs, r, templates, None, None).map(|(report, _)| report)
}
pub(crate) fn inspect_reconstructed_with_projections(
    inputs: &[topo::Input],
    r: &Request,
    templates: &[(String, &tracked::TrackedDuplex)],
    capture: Option<crate::verification::ancestry::TraceBudget>,
    replay: &mut crate::verification::replay_budget::ReplayBudget<'_>,
) -> Result<
    (
        Report,
        Option<crate::verification::projection::GraphProjections>,
    ),
    Error,
> {
    inspect_reconstructed_native(inputs, r, templates, capture, Some(replay))
}
fn inspect_reconstructed_native(
    inputs: &[topo::Input],
    r: &Request,
    templates: &[(String, &tracked::TrackedDuplex)],
    capture: Option<crate::verification::ancestry::TraceBudget>,
    replay: Option<&mut crate::verification::replay_budget::ReplayBudget<'_>>,
) -> Result<
    (
        Report,
        Option<crate::verification::projection::GraphProjections>,
    ),
    Error,
> {
    shape(r)?;
    let preparation_revision =
        topo::batch_revision_reconstructed(inputs, &r.preparation, templates)?;
    inspect_preflighted(inputs, r, preparation_revision, capture, replay)
}
fn inspect_preflighted(
    inputs: &[topo::Input],
    r: &Request,
    preparation_revision: String,
    capture: Option<crate::verification::ancestry::TraceBudget>,
    replay: Option<&mut crate::verification::replay_budget::ReplayBudget<'_>>,
) -> Result<
    (
        Report,
        Option<crate::verification::projection::GraphProjections>,
    ),
    Error,
> {
    let original = &r.preparation;
    let prepared = topo::prepare_with_ancestry(inputs, original, capture.clone(), replay)?;
    let sources = topo::identities(&prepared.inputs);
    let context_revision = binding(&preparation_revision, &sources, r)?;
    let fresh = topo::execute(&prepared, &original.engine_limits)?;
    let projection = projection::project_fresh(
        &prepared.evidence.substrates,
        &prepared.effective_steps,
        &original.engine_limits,
        &r.selection,
        &r.options,
        fresh,
    )
    .map_err(topo::projection_error)?;
    let preparation = prepared.evidence;
    let mut receipts = None;
    let (mut duplex, mut tracked, mut product, mut annotations) = (None, None, None, None);
    let mut bridge_notes=vec!["This source-backed graph bridge preserves unknown internal modification calls. Reference features describe correspondence, not inherited atoms, methylation or experimentally established material.".into()];
    bridge_notes.push("Any requested ideal flap cleavage describes only explicit graph geometry and declared new terminal chemistry. It is not a named nuclease profile, host-repair model, enzyme efficiency or biochemical mixture guarantee.".into());
    bridge_notes.push("The chosen component/orientation is an explicit graph view, not physical isolation or a reaction-mixture winner. Every other component, unpaired call and released activation remains in the retained processing evidence.".into());
    match &projection.outcome {
        projection::Outcome::Withheld(reason)=>bridge_notes.push(format!("Paired-contour proof withheld: {reason:?}; no ordinary duplex material or sequence document was created.")),
        projection::Outcome::Proven(view)=>match bridge(view){
            Err(reason)=>bridge_notes.push(reason),
            Ok((d,t))=>{
                let mut m=d.molecule("TOPO action product")?;
                m.metadata.insert("duplex_model".into(),ENGINE_VERSION.into());
                m.metadata.insert("topo:actions_context".into(),context_revision.clone());
                m.metadata.insert("topo:physical_evidence_required".into(),"true".into());
                match topo_annotations::project_native(&prepared.inputs,&preparation,view,&m,capture.is_some()){
                    Ok((evidence,captured))=>{receipts=captured;m.features=evidence.features.clone();annotations=Some(evidence);},
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
        input_sources: topo::identities(inputs),
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
    Ok((
        report,
        topo::graph_projections(
            prepared.inputs,
            prepared.source_scopes,
            prepared.source_origin_indices,
            prepared.native_origins,
            receipts,
        ),
    ))
}
impl Report {
    /// Bounded readonly binding and representability, not chemistry replay.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        shape(&self.request)?;
        if self.engine_version != ENGINE_VERSION
            || !digest(expected)
            || expected != self.context_revision
            || !digest(&self.preparation_revision)
            || self.input_sources.len() > 64
            || self.sources.is_empty()
            || self.sources.len() > 64
            || self.bridge_notes.len() > 16
            || self.bridge_notes.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("report source/request bounds changed"));
        }
        bounded(&self.request, MAX_REPORT_BYTES - 65536)?;
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
        let mut input_roles = BTreeSet::new();
        for source in &self.input_sources {
            if source.role.trim().is_empty()
                || source.role.len() > 256
                || !digest(&source.revision)
                || !input_roles.insert(&source.role)
            {
                return Err(bad("root input source identity changed"));
            }
        }
        let expanded =
            topo::expanded_source_identities(&self.input_sources, &self.request.preparation)?;
        if crate::encoding::bounded_json::hash(&expanded)?
            != crate::encoding::bounded_json::hash(&self.sources)?
        {
            return Err(bad("root and expanded native source identities differ"));
        }
        let original = &self.request.preparation;
        let original_revision = crate::encoding::bounded_json::hash(&(
            topo::ENGINE_VERSION,
            physical::ENGINE,
            physical::PRIMARY_MODEL,
            projection::ENGINE,
            &self.input_sources,
            original,
        ))?;
        if original_revision != self.preparation_revision
            || binding(&original_revision, &self.sources, &self.request)? != expected
        {
            return Err(bad("source/preparation/selection context changed"));
        }
        topo_worker::preparation(&self.preparation)?;
        projection::validate_report(
            &self.preparation.substrates,
            &topo::effective_steps(&self.preparation, original)?,
            &original.engine_limits,
            &self.request.selection,
            &self.request.options,
            &self.projection,
        )
        .map_err(topo::projection_error)?;
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
                    || m.metadata.get("topo:actions_context") != Some(&self.context_revision)
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
                    if crate::encoding::bounded_json::hash(&evidence.features)?
                        != crate::encoding::bounded_json::hash(&m.features)?
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
    if crate::encoding::bounded_json::hash(&inspect(inputs, r)?)?
        != crate::encoding::bounded_json::hash(report)?
    {
        return Err(bad("report differs from fresh source-backed processing"));
    }
    Ok(())
}
// Discovery is separate from an explicit selected component/seed. Its graph is
// retained for review but never imported as input to the projection job.
