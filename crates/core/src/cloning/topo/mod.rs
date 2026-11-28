//! Source-backed preparation declarations and fresh physical TOPO reconstruction.
//! No sequence-record metadata establishes activation, end chemistry or kit state.
use crate::{chemistry, cloning::ligation::tailing, sequence::analysis};
use collagene_model::{Alphabet, Error, Molecule, Sequence, Topology};
pub use collagene_topo as physical;
use physical::projection;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const ENGINE_VERSION: &str = "topo-source-adapter-3";
pub const ANNOTATION_ENGINE_VERSION: &str = "topo-annotated-reference-2";
pub const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
pub const RECIPE_KEY: &str = "topo:unverified_generation";
macro_rules! dto { ($($t:tt)*) => { #[derive(Debug,Clone,Serialize,Deserialize)] #[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))] #[serde(deny_unknown_fields)] $($t)* }; }
dto! {pub struct Input { pub role:String, pub molecule:Molecule }}
dto! {pub struct SourceIdentity { pub role:String,pub revision:String }}
// Forward reference interval, optionally read in reverse-complement sense. One
// circular origin crossing is permitted; length never exceeds one source turn.
dto! {pub struct SourceInterval {pub source_role:String,pub start:usize,pub length:usize,pub reverse:bool}}
dto! {pub struct PreparedStrand {pub id:String,pub segments:Vec<SourceInterval>,pub five_prime:physical::FivePrime,pub three_prime:physical::ThreePrime}}
// A declared preparation hypothesis, not inferred manufacturer state or a
// physical excision accounting assertion. Unrepresented reference DNA is not waste.
dto! {pub enum MaterialRecipe {
    Prepared {declaration_id:String,assumption:String,strands:Vec<PreparedStrand>,pairing:Vec<physical::PairSpan>},
    ExactPcr {material_id:String,template_role:String,forward:analysis::Primer,reverse:analysis::Primer,tail:Option<tailing::TailingBase>},
    SavedGraph {material_id:String,origin:crate::cloning::topo::graph_history::Capsule},
}}
dto! {pub struct GenerationRecipe {pub schema_version:u16,pub materials:Vec<MaterialRecipe>}}
dto! {pub struct Request {pub recipe:GenerationRecipe,pub steps:Vec<physical::Step>,pub engine_limits:physical::Limits,pub projection_options:projection::Options}}
dto! {pub enum Origin {
    ReferenceCorrespondence(SourceInterval),
    /// Full freshly computed PCR strand before optional tailing. This identifies
    /// an ideal polymerized copy; it asserts neither inherited atoms nor methylation.
    ExactPcrCopy {material:usize,template_role:String,template_start:usize,template_length:usize,forward_tail_length:usize,reverse_tail_length:usize,reverse:bool},
    SynthesizedTail {material:usize,base:tailing::TailingBase},
}}
dto! {pub struct PreparedMap {pub strand:usize,pub strand_start:usize,pub length:usize,pub origin:Origin}}
dto! {pub struct ReferenceGap {pub source_role:String,pub start:usize,pub end:usize}}
dto! {pub struct PreparationEvidence {pub substrates:physical::Substrates,pub maps:Vec<PreparedMap>,pub not_represented:Vec<ReferenceGap>,pub assumptions:Vec<String>,pub retained_steps:Vec<physical::Step>,pub layout:Vec<MaterialLayout>,pub graph_origins:Vec<GraphOriginEvidence>}}
dto! {pub struct TopoReport {pub engine_version:String,pub batch_revision:String,pub sources:Vec<SourceIdentity>,pub preparation:PreparationEvidence,pub physical:physical::Report,pub incomplete:bool,pub findings:Vec<String>}}
dto! {pub struct Selection {pub batch_revision:String,pub component:usize,pub reference_seed:physical::Call}}
dto! {pub struct Product {pub engine_version:String,pub batch_revision:String,pub selection:Selection,pub sources:Vec<SourceIdentity>,pub preparation:PreparationEvidence,pub projection:projection::ProjectionReport,pub annotations:Option<crate::cloning::topo::annotations::AnnotationEvidence>,pub reference:Option<Molecule>}}
dto! {#[derive(Default)] pub struct AcceptancePolicy {pub allow_nicks:bool,pub allow_unknown_terminal_chemistry:bool,pub allow_unapplied_steps:bool,pub allow_other_components:bool,pub allow_circular_rotation:bool,pub allow_reverse_complement:bool}}
dto! {pub enum TargetOutcome {AcceptedUnderDeclaredGraph,Rejected{reasons:Vec<String>},Withheld{reason:String}}}
dto! {pub struct TargetAssessment {pub engine_version:String,pub batch_revision:String,pub target_revision:String,pub policy:AcceptancePolicy,pub product:Product,pub outcome:TargetOutcome}}
// Whole-graph evidence and a selected annotated contour answer distinct current
// questions. There is one Product document, with documentary annotations only.
dto! {pub struct MaterialLayout{pub owner:Option<String>,pub material:usize,pub source_strands:Vec<usize>,pub source_roles:Vec<(String,String)>}}
dto! {pub struct GraphOriginEvidence{pub selection:crate::cloning::topo::graph_history::Selection,pub strand_start:usize,pub strand_end:usize,pub after_steps:usize,pub activation_unmap:Vec<(String,String)>,pub step_map:Vec<usize>,pub incomplete:bool}}
enum PreparedSourceScope {
    CurrentInput {
        index: usize,
    },
    RetainedGraph {
        role: String,
        scopes: Vec<crate::verification::projection::TopoScope>,
    },
}
pub(crate) struct PreparedContext {
    pub inputs: Vec<Input>,
    source_scopes: Vec<PreparedSourceScope>,
    source_origin_indices: Vec<Option<usize>>,
    native_origins: Vec<crate::material::run::RetainedReplay>,
    pub evidence: PreparationEvidence,
    pub effective_steps: Vec<physical::Step>,
    pub incomplete: bool,
    pub origin_source_incomplete: Vec<bool>,
    pub findings: Vec<String>,
}

fn bad(s: &str) -> Error {
    Error::Unsupported(format!("TOPO source adapter: {s}"))
}
fn invalid(s: &str) -> Error {
    Error::Format(format!("TOPO source adapter: {s}"))
}
fn text(s: &str, max: usize) -> Result<(), Error> {
    if s.is_empty() || s.len() > max || s.trim().is_empty() {
        Err(Error::TooLarge)
    } else {
        Ok(())
    }
}
struct Counter {
    n: usize,
    max: usize,
}
impl Write for Counter {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.n = self
            .n
            .checked_add(b.len())
            .filter(|n| *n <= self.max)
            .ok_or_else(|| std::io::Error::other("bounded TOPO evidence"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded<T: Serialize>(v: &T, max: usize) -> Result<(), Error> {
    serde_json::to_writer(Counter { n: 0, max }, v).map_err(|_| Error::TooLarge)
}
struct Hasher(blake3::Hasher);
impl Write for Hasher {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0.update(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn hash<T: Serialize>(v: &T) -> Result<String, Error> {
    let mut h = Hasher(blake3::Hasher::new());
    serde_json::to_writer(&mut h, v).map_err(|_| Error::TooLarge)?;
    Ok(h.0.finalize().to_hex().to_string())
}
fn engine_error(error: physical::Error) -> Error {
    match error {
        physical::Error::Invalid(reason) => invalid(&reason),
        physical::Error::Limit { .. } => Error::TooLarge,
    }
}
fn projection_error(error: projection::Error) -> Error {
    match error {
        projection::Error::Engine(error) => engine_error(error),
        projection::Error::InvalidSelection(reason) => invalid(reason),
        projection::Error::InvalidLimits => invalid("invalid paired-contour limits"),
        projection::Error::ReportBytes { .. } => Error::TooLarge,
    }
}
fn source<'a>(inputs: &'a [Input], role: &str) -> Result<&'a Molecule, Error> {
    inputs
        .iter()
        .find(|x| x.role == role)
        .map(|x| &x.molecule)
        .ok_or_else(|| invalid("unknown source role"))
}
fn interval(inputs: &[Input], p: &SourceInterval) -> Result<(), Error> {
    text(&p.source_role, 256)?;
    let m = source(inputs, &p.source_role)?;
    let n = m.sequence.len();
    if p.start >= n
        || p.length == 0
        || p.length > n
        || m.topology == Topology::Linear && p.start.checked_add(p.length).is_none_or(|e| e > n)
    {
        return Err(invalid("invalid source interval"));
    }
    Ok(())
}
pub fn source_preflight(inputs: &[Input], r: &Request) -> Result<(), Error> {
    source_preflight_inner(inputs, r, &BTreeSet::new())
}
fn source_preflight_inner(
    inputs: &[Input],
    r: &Request,
    reconstructed_templates: &BTreeSet<String>,
) -> Result<(), Error> {
    physical::steps_preflight(&r.steps, &r.engine_limits).map_err(engine_error)?;
    let roots = r.recipe.materials.iter().filter_map(|m| {
        if let MaterialRecipe::SavedGraph { origin, .. } = m {
            Some(crate::material::capsule::NativeOriginRef::TopoGraph(origin))
        } else {
            None
        }
    });
    crate::material::capsule::native_origins_preflight(roots)?;
    bounded(r, MAX_REPORT_BYTES - 65536)?;
    if r.recipe
        .materials
        .iter()
        .any(|m| matches!(m, MaterialRecipe::SavedGraph { .. }))
    {
        let expanded = expand(inputs, r, None)?;
        let mut roles = reconstructed_templates.clone();
        roles.extend(
            expanded
                .bindings
                .iter()
                .filter(|b| b.origin.is_some())
                .map(|b| b.role.clone()),
        );
        source_preflight_flat(&expanded.inputs, &expanded.request, &roles)
    } else {
        source_preflight_flat(inputs, r, reconstructed_templates)
    }
}
fn source_preflight_flat(
    inputs: &[Input],
    r: &Request,
    reconstructed_templates: &BTreeSet<String>,
) -> Result<(), Error> {
    physical::steps_preflight(&r.steps, &r.engine_limits).map_err(engine_error)?;
    if inputs.len() > 64
        || r.recipe.schema_version != 1
        || r.recipe.materials.is_empty()
        || r.recipe.materials.len() > 64
        || r.steps.len() > 16
    {
        return Err(Error::TooLarge);
    }
    let refs: Vec<_> = inputs.iter().map(|x| &x.molecule).collect();
    collagene_model::bounds::model_preflight(&refs)?;
    let mut roles = BTreeSet::new();
    let mut bases = 0usize;
    for i in inputs {
        text(&i.role, 256)?;
        if !roles.insert(i.role.as_str()) || i.molecule.sequence.alphabet() != Alphabet::Dna {
            return Err(invalid("unique DNA source roles required"));
        }
        bases = bases
            .checked_add(i.molecule.sequence.len())
            .filter(|n| *n <= 100_000)
            .ok_or(Error::TooLarge)?;
    }
    let mut strands = 0usize;
    let mut spans = 0usize;
    let mut parts = 0usize;
    let mut output = 0usize;
    let mut ids = BTreeSet::new();
    for m in &r.recipe.materials {
        match m {
            MaterialRecipe::SavedGraph { .. } => return Err(bad("unexpanded saved graph")),
            MaterialRecipe::Prepared {
                declaration_id,
                assumption,
                strands: ss,
                pairing,
            } => {
                text(declaration_id, 256)?;
                text(assumption, 4096)?;
                strands += ss.len();
                spans += pairing.len();
                if strands > 32 || spans > 1024 {
                    return Err(Error::TooLarge);
                }
                for s in ss {
                    text(&s.id, 256)?;
                    if s.segments.len() > 1024 {
                        return Err(Error::TooLarge);
                    }
                    if !ids.insert(s.id.as_str()) || s.segments.is_empty() {
                        return Err(invalid("unique strand ids and nonempty segments required"));
                    }
                    if let physical::ThreePrime::Phosphotyrosyl { activation_id, .. } =
                        &s.three_prime
                    {
                        text(activation_id, 256)?;
                    }
                    for p in &s.segments {
                        parts += 1;
                        if parts > 2048 {
                            return Err(Error::TooLarge);
                        }
                        interval(inputs, p)?;
                        let original = source(inputs, &p.source_role)?;
                        if chemistry::has_unmapped_physical_evidence(original)
                            || original.metadata.contains_key("tracked_modifications")
                            || original.metadata.keys().any(|k| k.starts_with("topo:"))
                        {
                            // A genuine freshly reconstructed Step can be a PCR
                            // template, but activation from its tracked strands
                            // is not implemented by the Prepared-reference path.
                            if reconstructed_templates.contains(&p.source_role) {
                                return Err(bad(
                                    "Prepared activation hypotheses require original references; activation from actual tracked Step material is not modeled",
                                ));
                            }
                            return Err(invalid(
                                "Prepared activation hypotheses require original references; a saved/marked physical document must use its exact native graph instead of redeclaring stock",
                            ));
                        }
                        output = output
                            .checked_add(p.length)
                            .filter(|n| *n <= 100_000)
                            .ok_or(Error::TooLarge)?;
                    }
                }
                for p in pairing {
                    if p.strand_a >= ss.len() || p.strand_b >= ss.len() {
                        return Err(invalid("prepared pairing indices are material-local"));
                    }
                }
            }
            MaterialRecipe::ExactPcr {
                material_id,
                template_role,
                forward,
                reverse,
                ..
            } => {
                text(material_id, 256)?;
                text(template_role, 256)?;
                let template = source(inputs, template_role)?;
                if chemistry::has_unmapped_physical_evidence(template)
                    && !reconstructed_templates.contains(template_role)
                {
                    return Err(bad(
                        "PCR requires reconstructible original template, not imported physical evidence",
                    ));
                }
                strands += 2;
                if strands > 32 {
                    return Err(Error::TooLarge);
                }
                for p in [forward, reverse] {
                    text(&p.name, 256)?;
                    if p.annealing.len() > 1000
                        || p.tail.len() > 1000
                        || p.annealing.alphabet() != Alphabet::Dna
                    {
                        return Err(Error::TooLarge);
                    }
                }
            }
        }
    }
    // Scalar limits are checked without constructing source DNA or trusting a report.
    let l = &r.engine_limits;
    if l.max_source_bases == 0
        || l.max_source_bases > 500_000
        || l.max_strands == 0
        || l.max_strands > 64
        || l.max_pair_spans == 0
        || l.max_pair_spans > 4096
        || l.max_steps == 0
        || l.max_steps > 64
        || l.max_work_cells == 0
        || l.max_work_cells > 128_000_000
        || l.max_work_bytes == 0
        || l.max_work_bytes > 128 * 1024 * 1024
        || l.max_output_bases == 0
        || l.max_output_bases > 2_000_000
        || l.max_report_bytes == 0
        || l.max_report_bytes > 32 * 1024 * 1024
    {
        return Err(Error::TooLarge);
    }
    let p = &r.projection_options.limits;
    if p.max_pairs == 0
        || p.max_pairs > 500_000
        || p.max_work_cells == 0
        || p.max_work_cells > 128_000_000
        || p.max_work_bytes == 0
        || p.max_work_bytes > 128 * 1024 * 1024
        || p.max_mapping_segments == 0
        || p.max_mapping_segments > 16384
        || p.max_report_bytes == 0
        || p.max_report_bytes > 32 * 1024 * 1024
    {
        return Err(Error::TooLarge);
    }
    bounded(r, 1024 * 1024)?;
    // The existing exact-PCR search is O(template × primer). It is run once
    // to bound candidate materialization and again by the geometry enumerator.
    // Charge both orientations of both passes across the entire recipe before
    // any search, source hashing, or preparation allocation can begin. This is
    // a separate preparation-search allowance under the requested work cap;
    // the physical engine still independently reserves its graph work.
    let mut search_work = 0usize;
    for material in &r.recipe.materials {
        if let MaterialRecipe::ExactPcr {
            template_role,
            forward,
            reverse,
            ..
        } = material
        {
            let template = source(inputs, template_role)?;
            let cells = forward
                .annealing
                .len()
                .checked_add(reverse.annealing.len())
                .and_then(|n| n.checked_mul(template.sequence.len()))
                .and_then(|n| n.checked_mul(4))
                .ok_or(Error::TooLarge)?;
            search_work = search_work
                .checked_add(cells)
                .filter(|n| *n <= r.engine_limits.max_work_cells)
                .ok_or(Error::TooLarge)?;
        }
    }
    for i in inputs {
        i.molecule.validate()?;
    }
    Ok(())
}
pub(crate) fn identities(inputs: &[Input]) -> Vec<SourceIdentity> {
    inputs
        .iter()
        .map(|x| SourceIdentity {
            role: x.role.clone(),
            revision: x.molecule.revision(),
        })
        .collect()
}
pub(crate) fn documentary_batch_revision(inputs: &[Input], r: &Request) -> Result<String, Error> {
    batch_binding(inputs, r)
}
pub(crate) fn documentary_sources(
    inputs: &[Input],
    r: &Request,
) -> Result<Vec<SourceIdentity>, Error> {
    expanded_source_identities(&identities(inputs), r)
}
/// Bounded documentary role/revision expansion; no source DNA allocation, PCR or
/// graph replay. Native ownership is separately checked by the shared walker.
pub(crate) fn expanded_source_identities(
    sources: &[SourceIdentity],
    r: &Request,
) -> Result<Vec<SourceIdentity>, Error> {
    let roots = r.recipe.materials.iter().filter_map(|m| {
        if let MaterialRecipe::SavedGraph { origin, .. } = m {
            Some(crate::material::capsule::NativeOriginRef::TopoGraph(origin))
        } else {
            None
        }
    });
    crate::material::capsule::native_origins_preflight(roots)?;
    fn visit(
        sources: &[SourceIdentity],
        r: &Request,
        owned: bool,
        counter: &mut usize,
        used: &mut BTreeSet<String>,
        out: &mut Vec<SourceIdentity>,
    ) -> Result<(), Error> {
        if *counter > 128 {
            return Err(Error::TooLarge);
        }
        let node = *counter;
        *counter += 1;
        for (i, s) in sources.iter().enumerate() {
            out.push(SourceIdentity {
                role: if owned {
                    resolved_role(node, i, used)?
                } else {
                    s.role.clone()
                },
                revision: s.revision.clone(),
            });
        }
        if out.len() > 64 {
            return Err(Error::TooLarge);
        }
        for m in &r.recipe.materials {
            if let MaterialRecipe::SavedGraph { origin, .. } = m {
                let record = origin.validate_shallow()?;
                let sources = identities(&origin.source_snapshots);
                visit(&sources, &record.request.request, true, counter, used, out)?;
            }
        }
        Ok(())
    }
    let mut out = vec![];
    let mut used = sources
        .iter()
        .map(|s| s.role.clone())
        .collect::<BTreeSet<_>>();
    visit(sources, r, false, &mut 0, &mut used, &mut out)?;
    Ok(out)
}
pub fn batch_revision(inputs: &[Input], r: &Request) -> Result<String, Error> {
    source_preflight(inputs, r)?;
    batch_binding(inputs, r)
}
/// Private synthesis-only route. The caller supplies actual freshly reconstructed
/// Step material, never a deserialized report or a metadata-derived duplex.
pub(crate) fn source_preflight_reconstructed(
    inputs: &[Input],
    r: &Request,
    templates: &[(String, &crate::chemistry::tracked_duplex::TrackedDuplex)],
) -> Result<(), Error> {
    if inputs.len() > 32 || templates.len() > 32 || templates.len() > inputs.len() {
        return Err(Error::TooLarge);
    }
    let mut roles = BTreeSet::new();
    for (role, _) in templates {
        text(role, 256)?;
        if !roles.insert(role.clone()) || inputs.iter().filter(|i| i.role == *role).count() != 1 {
            return Err(invalid(
                "reconstructed PCR templates require unique actual source roles",
            ));
        }
    }
    // This validates all input models, recipe counts/segments and search budgets
    // before recursive tracked validation or source/recipe hashing.
    source_preflight_inner(inputs, r, &roles)?;
    bounded(&templates, MAX_REPORT_BYTES)?;
    let mut template_bases = 0usize;
    for (role, tracked) in templates {
        let mut used = false;
        for material in &r.recipe.materials {
            match material {
                MaterialRecipe::SavedGraph { .. } => {}
                MaterialRecipe::Prepared { strands, .. } => {
                    if strands
                        .iter()
                        .any(|s| s.segments.iter().any(|p| p.source_role == *role))
                    {
                        return Err(bad(
                            "tracked Step material cannot supply reference-only Prepared intervals; only fresh exact-PCR synthesis is supported",
                        ));
                    }
                }
                MaterialRecipe::ExactPcr { template_role, .. } => {
                    used |= template_role == role;
                }
            }
        }
        if !used {
            return Err(invalid(
                "reconstructed template role is not used by any exact-PCR material",
            ));
        }
        template_bases = template_bases
            .checked_add(tracked.duplex.top.len())
            .and_then(|n| n.checked_add(tracked.duplex.bottom.len()))
            .filter(|n| *n <= 200_000)
            .ok_or(Error::TooLarge)?;
        tracked.validate()?;
        let molecule = source(inputs, role)?;
        if tracked.duplex.top != molecule.sequence || tracked.duplex.topology != molecule.topology {
            return Err(invalid(
                "actual tracked PCR template differs from its exact Step sequence/topology",
            ));
        }
        if tracked.duplex.bottom_start != 0
            || tracked.duplex.bottom.len() != tracked.duplex.top.len()
            || !tracked.duplex.nicks.is_empty()
        {
            return Err(bad(
                "exact PCR of nicked or staggered tracked templates is not modeled; actual complete aligned strands are required",
            ));
        }
        if molecule.metadata.keys().any(|key| {
            key.starts_with("variant:")
                || key.starts_with("recombination:")
                || key.starts_with("source:recombination:")
                || key == crate::cloning::recombination::CHEMISTRY_UNMAPPED_KEY
        }) {
            return Err(invalid(
                "reference-only, variant or unmapped reference ancestry cannot establish a physical PCR template",
            ));
        }
    }
    Ok(())
}
pub(crate) fn batch_revision_reconstructed(
    inputs: &[Input],
    r: &Request,
    templates: &[(String, &crate::chemistry::tracked_duplex::TrackedDuplex)],
) -> Result<String, Error> {
    source_preflight_reconstructed(inputs, r, templates)?;
    batch_binding(inputs, r)
}
fn batch_binding(inputs: &[Input], r: &Request) -> Result<String, Error> {
    hash(&(
        ENGINE_VERSION,
        physical::ENGINE,
        physical::PRIMARY_MODEL,
        projection::ENGINE,
        identities(inputs),
        r,
    ))
}
fn segment_dna(inputs: &[Input], p: &SourceInterval) -> Result<Sequence, Error> {
    let m = source(inputs, &p.source_role)?;
    let n = m.sequence.len();
    let dna: String = (0..p.length)
        .map(|j| m.sequence.as_str().as_bytes()[(p.start + j) % n] as char)
        .collect();
    let s = Sequence::dna(&dna)?;
    if p.reverse {
        s.reverse_complement()
    } else {
        Ok(s)
    }
}
fn mark(coverage: &mut [bool], start: usize, length: usize) {
    let n = coverage.len();
    for k in 0..length {
        coverage[(start + k) % n] = true;
    }
}
struct Expanded {
    inputs: Vec<Input>,
    source_scopes: Vec<PreparedSourceScope>,
    source_origin_keys: Vec<Option<String>>,
    request: Request,
    bindings: Vec<crate::cloning::topo::graph_history::SourceBinding>,
    origins: Vec<GraphOriginEvidence>,
    layout: Vec<MaterialLayout>,
}
fn resolved_role(node: usize, index: usize, used: &mut BTreeSet<String>) -> Result<String, Error> {
    let base = format!("graph{node}-source{index}");
    for suffix in 0..=128 {
        let name = if suffix == 0 {
            base.clone()
        } else {
            format!("{base}-{suffix}")
        };
        if used.insert(name.clone()) {
            return Ok(name);
        }
    }
    Err(Error::TooLarge)
}
fn expand(
    inputs: &[Input],
    r: &Request,
    root_bindings: Option<&[crate::cloning::topo::graph_history::SourceBinding]>,
) -> Result<Expanded, Error> {
    let roots = r.recipe.materials.iter().filter_map(|m| {
        if let MaterialRecipe::SavedGraph { origin, .. } = m {
            Some(crate::material::capsule::NativeOriginRef::TopoGraph(origin))
        } else {
            None
        }
    });
    crate::material::capsule::native_origins_preflight(roots)?;
    let mut out = Expanded {
        inputs: vec![],
        source_scopes: vec![],
        source_origin_keys: vec![],
        request: Request {
            recipe: GenerationRecipe {
                schema_version: 1,
                materials: vec![],
            },
            steps: vec![],
            engine_limits: r.engine_limits.clone(),
            projection_options: r.projection_options.clone(),
        },
        bindings: vec![],
        origins: vec![],
        layout: vec![],
    };
    let mut seen = BTreeSet::new();
    let mut used = inputs
        .iter()
        .map(|i| i.role.clone())
        .collect::<BTreeSet<_>>();
    let mut counter = 0;
    expand_node(
        inputs,
        r,
        None,
        root_bindings,
        &[],
        &mut seen,
        &mut used,
        &mut counter,
        &mut out,
    )?;
    Ok(out)
}
fn native_activation_identity(
    strand: &PreparedStrand,
    inputs: &[Input],
    bindings: &[crate::cloning::topo::graph_history::SourceBinding],
) -> Result<String, Error> {
    let enzyme = match &strand.three_prime {
        physical::ThreePrime::Phosphotyrosyl { enzyme, .. } => *enzyme,
        _ => return Err(bad("activation identity requires actual declared adduct")),
    };
    let mut h = blake3::Hasher::new();
    h.update(b"collagene-native-topo-activation-call-v1\0");
    h.update(hash(&(strand.five_prime, enzyme))?.as_bytes());
    let mut calls = 0usize;
    for segment in &strand.segments {
        let binding = bindings
            .iter()
            .find(|b| b.role == segment.source_role)
            .ok_or_else(|| bad("activated graph source lacks its exact original event"))?;
        let molecule = source(inputs, &segment.source_role)?;
        let n = molecule.sequence.len();
        interval(inputs, segment)?;
        let event = hash(&binding.event)?;
        for j in 0..segment.length {
            calls = calls
                .checked_add(1)
                .filter(|n| *n <= 100_000)
                .ok_or(Error::TooLarge)?;
            let position = (segment.start
                + if segment.reverse {
                    segment.length - 1 - j
                } else {
                    j
                })
                % n;
            h.update(event.as_bytes());
            h.update(&(position as u64).to_le_bytes());
            h.update(&[u8::from(segment.reverse)]);
        }
    }
    Ok(format!("activation-{}", h.finalize().to_hex()))
}
fn expand_node(
    inputs: &[Input],
    r: &Request,
    owner: Option<&str>,
    bindings: Option<&[crate::cloning::topo::graph_history::SourceBinding]>,
    scopes: &[crate::verification::projection::TopoScope],
    seen: &mut BTreeSet<String>,
    used: &mut BTreeSet<String>,
    counter: &mut usize,
    out: &mut Expanded,
) -> Result<Vec<usize>, Error> {
    if *counter > 128 {
        return Err(Error::TooLarge);
    }
    let node_step_base = out.request.steps.len();
    let node = *counter;
    *counter += 1;
    let mut roles = std::collections::BTreeMap::new();
    for (i, input) in inputs.iter().enumerate() {
        let role = if owner.is_some() {
            resolved_role(node, i, used)?
        } else {
            input.role.clone()
        };
        roles.insert(input.role.clone(), role.clone());
        let binding = bindings
            .and_then(|bindings| bindings.iter().find(|binding| binding.role == input.role));
        out.source_origin_keys.push(
            binding
                .and_then(|binding| binding.origin.as_ref().map(|origin| (binding, origin)))
                .map(|(binding, origin)| hash(&(&binding.event, origin.selection())))
                .transpose()?,
        );
        out.source_scopes.push(if scopes.is_empty() {
            PreparedSourceScope::CurrentInput { index: i }
        } else {
            PreparedSourceScope::RetainedGraph {
                role: input.role.clone(),
                scopes: scopes.to_vec(),
            }
        });
        out.inputs.push(Input {
            role: role.clone(),
            molecule: input.molecule.clone(),
        });
        if let Some(binding) = binding {
            let mut b = binding.clone();
            b.role = role;
            out.bindings.push(b)
        }
    }
    let mut strand_map = vec![];
    for (material, m) in r.recipe.materials.iter().enumerate() {
        let start = flat_strands(&out.request.recipe.materials);
        match m {
            MaterialRecipe::Prepared {
                declaration_id,
                assumption,
                strands,
                pairing,
            } => {
                let originals = strands;
                let mut strands = strands.clone();
                for (i, s) in strands.iter_mut().enumerate() {
                    s.id = if owner.is_some() {
                        format!("graph{node}-strand{material}-{i}")
                    } else {
                        s.id.clone()
                    };
                    for segment in &mut s.segments {
                        segment.source_role = roles
                            .get(&segment.source_role)
                            .ok_or_else(|| bad("missing original graph source role"))?
                            .clone();
                    }
                    if bindings.is_some()
                        && matches!(s.three_prime, physical::ThreePrime::Phosphotyrosyl { .. })
                    {
                        let id = native_activation_identity(
                            &originals[i],
                            inputs,
                            bindings
                                .ok_or_else(|| bad("native activation source bindings absent"))?,
                        )?;
                        if let physical::ThreePrime::Phosphotyrosyl { activation_id, .. } =
                            &mut s.three_prime
                        {
                            *activation_id = id;
                        }
                    }
                }
                let n = strands.len();
                out.request.recipe.materials.push(MaterialRecipe::Prepared {
                    declaration_id: if owner.is_some() {
                        format!("graph{node}-material{material}")
                    } else {
                        declaration_id.clone()
                    },
                    assumption: assumption.clone(),
                    strands,
                    pairing: pairing.clone(),
                });
                strand_map.extend(start..start + n);
            }
            MaterialRecipe::ExactPcr {
                material_id,
                template_role,
                forward,
                reverse,
                tail,
            } => {
                out.request.recipe.materials.push(MaterialRecipe::ExactPcr {
                    material_id: if owner.is_some() {
                        format!("graph{node}-material{material}")
                    } else {
                        material_id.clone()
                    },
                    template_role: roles
                        .get(template_role)
                        .ok_or_else(|| bad("missing graph PCR source"))?
                        .clone(),
                    forward: forward.clone(),
                    reverse: reverse.clone(),
                    tail: *tail,
                });
                strand_map.extend(start..start + 2);
            }
            MaterialRecipe::SavedGraph { origin, .. } => {
                if !seen.insert(origin.selection.record_hash.clone()) {
                    return Err(bad(
                        "the same saved graph or ancestor is supplied twice; no extra stock instance is inferred",
                    ));
                }
                let record = origin.validate_shallow()?;
                if record.engine_version != crate::cloning::topo::graph_history::ENGINE_VERSION
                    || record.scientific_engine != crate::verification::ENGINE_VERSION
                {
                    return Err(bad(
                        "historical graph model cannot provide current material",
                    ));
                }
                let steps_before = out.request.steps.len();
                let activation_namespace = origin.selection.record_hash.clone();
                let mut child_scopes = scopes.to_vec();
                child_scopes.push(crate::verification::projection::TopoScope {
                    occurrence: u32::try_from(*counter).map_err(|_| Error::TooLarge)?,
                    binding: std::sync::Arc::new(
                        crate::verification::ancestry::OriginBinding::TopoGraph {
                            selection: origin.selection.clone(),
                            context_revision: record.context_revision.clone(),
                            sources: record
                                .request
                                .bindings
                                .iter()
                                .map(|binding| crate::verification::ancestry::TopoRootBinding {
                                    role: binding.role.clone(),
                                    event: Some(binding.event.clone()),
                                })
                                .collect(),
                        },
                    ),
                });
                let indices = expand_node(
                    &origin.source_snapshots,
                    &record.request.request,
                    Some(&activation_namespace),
                    Some(&record.request.bindings),
                    &child_scopes,
                    seen,
                    used,
                    counter,
                    out,
                )?;
                let end = flat_strands(&out.request.recipe.materials);
                let mut activation_unmap = vec![];
                // Own primitive declarations acquire native-owner namespaces;
                // inherited child namespaces are already canonical in this record.
                let mut local = 0;
                for old in &record.request.request.recipe.materials {
                    match old {
                        MaterialRecipe::Prepared { strands, .. } => {
                            for st in strands {
                                if let physical::ThreePrime::Phosphotyrosyl { .. } = &st.three_prime
                                {
                                    let global = indices[local];
                                    let current =
                                        flat_activation(&out.request.recipe.materials, global)
                                            .ok_or_else(|| bad("lost original activated strand"))?;
                                    activation_unmap.push((current.clone(), current));
                                }
                                local += 1;
                            }
                        }
                        MaterialRecipe::ExactPcr { .. } => local += 2,
                        MaterialRecipe::SavedGraph { origin, .. } => {
                            let count = origin.strand_count()?;
                            for k in 0..count {
                                if let Some(current) = flat_activation(
                                    &out.request.recipe.materials,
                                    indices[local + k],
                                ) {
                                    activation_unmap.push((current.clone(), current));
                                }
                            }
                            local += count;
                        }
                    }
                }
                // Include nested transfer prefixes, all in the same physical job.
                let step_map = (steps_before..out.request.steps.len()).collect::<Vec<_>>();
                out.origins.push(GraphOriginEvidence {
                    selection: origin.selection.clone(),
                    strand_start: start,
                    strand_end: end,
                    after_steps: out.request.steps.len(),
                    activation_unmap,
                    step_map,
                    incomplete: record.incomplete,
                });
                strand_map.extend(indices);
            }
        }
        let end = flat_strands(&out.request.recipe.materials);
        out.layout.push(MaterialLayout {
            owner: owner.map(str::to_owned),
            material,
            source_strands: (start..end).collect(),
            source_roles: roles.iter().map(|(a, b)| (a.clone(), b.clone())).collect(),
        });
    }
    for step in &r.steps {
        out.request
            .steps
            .push(remap_step(step, &strand_map, node_step_base)?);
    }
    if out.inputs.len() > 64
        || out.request.recipe.materials.len() > 64
        || out.request.steps.len() > 64
    {
        return Err(Error::TooLarge);
    }
    Ok(strand_map)
}
fn remap_step(
    step: &physical::Step,
    strands: &[usize],
    step_base: usize,
) -> Result<physical::Step, Error> {
    let call = |c: physical::Call| -> Result<physical::Call, Error> {
        Ok(physical::Call {
            strand: *strands
                .get(c.strand)
                .ok_or_else(|| bad("call outside native graph source namespace"))?,
            index: c.index,
        })
    };
    Ok(match step {
        physical::Step::Transfer(t) => physical::Step::Transfer(physical::Transfer {
            donor: call(t.donor)?,
            acceptor: call(t.acceptor)?,
            geometry: t.geometry,
        }),
        physical::Step::Process(physical::processing::Action::CleaveThreePrimeFlap {
            bond,
            displacement_step,
            profile,
            assumption,
        }) => physical::Step::Process(physical::processing::Action::CleaveThreePrimeFlap {
            bond: physical::processing::Bond {
                upstream: call(bond.upstream)?,
                downstream: call(bond.downstream)?,
            },
            displacement_step: step_base
                .checked_add(*displacement_step)
                .ok_or(Error::TooLarge)?,
            profile: *profile,
            assumption: assumption.clone(),
        }),
        physical::Step::Process(physical::processing::Action::SealNick {
            selection,
            ligase_assumption,
        }) => physical::Step::Process(physical::processing::Action::SealNick {
            selection: physical::nick_sealing::NickSelection {
                three: call(selection.three)?,
                five: call(selection.five)?,
            },
            ligase_assumption: ligase_assumption.clone(),
        }),
    })
}
fn flat_strands(materials: &[MaterialRecipe]) -> usize {
    materials
        .iter()
        .map(|m| match m {
            MaterialRecipe::Prepared { strands, .. } => strands.len(),
            MaterialRecipe::ExactPcr { .. } => 2,
            MaterialRecipe::SavedGraph { .. } => unreachable!("expanded primitive recipe"),
        })
        .sum()
}
fn flat_activation(materials: &[MaterialRecipe], index: usize) -> Option<String> {
    let mut offset = 0;
    for m in materials {
        match m {
            MaterialRecipe::Prepared { strands, .. } => {
                if index < offset + strands.len() {
                    return match &strands[index - offset].three_prime {
                        physical::ThreePrime::Phosphotyrosyl { activation_id, .. } => {
                            Some(activation_id.clone())
                        }
                        _ => None,
                    };
                }
                offset += strands.len();
            }
            MaterialRecipe::ExactPcr { .. } => offset += 2,
            MaterialRecipe::SavedGraph { .. } => return None,
        }
    }
    None
}
fn prepare_declared(
    inputs: &[Input],
    r: &Request,
    mut budget: Option<&mut analysis::ExactSearchBudget>,
) -> Result<PreparationEvidence, Error> {
    let mut substrates = physical::Substrates {
        strands: vec![],
        pairing: vec![],
    };
    let mut maps = vec![];
    let mut assumptions=vec!["Prepared-material declarations are hypotheses supplied by the caller; no supplier state, physical ancestry, methylation or manufacturing byproducts are inferred from reference DNA.".into()];
    let mut coverage: Vec<Vec<bool>> = inputs
        .iter()
        .map(|i| vec![false; i.molecule.sequence.len()])
        .collect();
    for (material, m) in r.recipe.materials.iter().enumerate() {
        match m {
            MaterialRecipe::SavedGraph { .. } => return Err(bad("unexpanded saved graph")),
            MaterialRecipe::Prepared {
                declaration_id,
                assumption,
                strands,
                pairing,
            } => {
                assumptions.push(format!("{declaration_id}: {assumption}"));
                let offset = substrates.strands.len();
                for s in strands {
                    let strand = substrates.strands.len();
                    let mut dna = String::new();
                    for p in &s.segments {
                        let seq = segment_dna(inputs, p)?;
                        maps.push(PreparedMap {
                            strand,
                            strand_start: dna.len(),
                            length: p.length,
                            origin: Origin::ReferenceCorrespondence(p.clone()),
                        });
                        dna.push_str(seq.as_str());
                        let i = inputs.iter().position(|x| x.role == p.source_role).unwrap();
                        mark(&mut coverage[i], p.start, p.length);
                    }
                    substrates.strands.push(physical::Strand {
                        role: declaration_id.clone(),
                        id: s.id.clone(),
                        sequence: Sequence::dna(&dna)?,
                        five_prime: s.five_prime,
                        three_prime: s.three_prime.clone(),
                    });
                }
                for p in pairing {
                    let mut p = p.clone();
                    p.strand_a += offset;
                    p.strand_b += offset;
                    substrates.pairing.push(p);
                }
            }
            MaterialRecipe::ExactPcr {
                material_id,
                template_role,
                forward,
                reverse,
                tail,
            } => {
                let template = source(inputs, template_role)?;
                // Existing PCR enumerator has a separate materialization bound. Guard its
                // conservative pair-count × maximum product length before invoking it here.
                let mut f = 0usize;
                let mut b = 0usize;
                for primer in [forward, reverse] {
                    let canonical = analysis::Primer::new(
                        &primer.name,
                        primer.annealing.as_str(),
                        &primer.tail,
                    )?;
                    f += match budget.as_deref_mut() {
                        Some(budget) => analysis::find_matches_with_budget(
                            template,
                            canonical.annealing.as_str(),
                            false,
                            budget,
                        )?,
                        None => {
                            analysis::find_matches(template, canonical.annealing.as_str(), false)?
                        }
                    }
                    .len();
                    let reverse = canonical.annealing.reverse_complement()?;
                    b += match budget.as_deref_mut() {
                        Some(budget) => analysis::find_matches_with_budget(
                            template,
                            reverse.as_str(),
                            false,
                            budget,
                        )?,
                        None => analysis::find_matches(template, reverse.as_str(), false)?,
                    }
                    .len();
                }
                let upper = f
                    .checked_mul(b)
                    .and_then(|x| x.checked_mul(template.sequence.len() + 4000))
                    .filter(|x| *x <= 1_000_000)
                    .ok_or(Error::TooLarge)?;
                let _ = upper;
                let possibilities = match budget.as_deref_mut() {
                    Some(budget) => analysis::simulate_pcr_geometry_view_with_budget(
                        analysis::TemplateView::reference(template),
                        forward,
                        reverse,
                        budget,
                    )?,
                    None => analysis::simulate_pcr_geometry(template, forward, reverse)?,
                };
                if possibilities.len() != 1 {
                    return Err(bad(
                        "PCR preparation requires exactly one exact current geometry",
                    ));
                }
                let g = &possibilities[0];
                let prepared_bases: usize =
                    substrates.strands.iter().map(|s| s.sequence.len()).sum();
                let new_bases = g
                    .product
                    .sequence
                    .len()
                    .checked_add(usize::from(tail.is_some()))
                    .and_then(|n| n.checked_mul(2))
                    .ok_or(Error::TooLarge)?;
                if prepared_bases
                    .checked_add(new_bases)
                    .is_none_or(|n| n > 100_000)
                {
                    return Err(Error::TooLarge);
                }
                let original = chemistry::Duplex::exact_pcr(&g.product.sequence)?;
                let d = if let Some(base) = tail {
                    tailing::tail_blunt_duplex(&original, *base)?.duplex
                } else {
                    original
                };
                let strand = substrates.strands.len();
                for (reverse, dna) in [(false, &d.top), (true, &d.bottom)] {
                    let index = substrates.strands.len();
                    maps.push(PreparedMap {
                        strand: index,
                        strand_start: 0,
                        length: g.product.sequence.len(),
                        origin: Origin::ExactPcrCopy {
                            material,
                            template_role: template_role.clone(),
                            template_start: g.product.forward_start,
                            template_length: g.template_length,
                            forward_tail_length: g.forward_tail_length,
                            reverse_tail_length: g.reverse_tail_length,
                            reverse,
                        },
                    });
                    if let Some(base) = tail {
                        maps.push(PreparedMap {
                            strand: index,
                            strand_start: g.product.sequence.len(),
                            length: 1,
                            origin: Origin::SynthesizedTail {
                                material,
                                base: *base,
                            },
                        });
                    }
                    substrates.strands.push(physical::Strand {
                        role: material_id.clone(),
                        id: format!("{material_id}/{}", if reverse { "bottom" } else { "top" }),
                        sequence: dna.clone(),
                        five_prime: physical::FivePrime::Hydroxyl,
                        three_prime: physical::ThreePrime::Hydroxyl,
                    });
                }
                // Both new 3′ residues are unpaired. The fresh tailing transform uses
                // bottom_start=-1, placing the original complementary portions at zero.
                substrates.pairing.push(physical::PairSpan {
                    strand_a: strand,
                    start_a: 0,
                    strand_b: strand + 1,
                    start_b: 0,
                    length: g.product.sequence.len(),
                });
                let i = inputs
                    .iter()
                    .position(|x| x.role == *template_role)
                    .unwrap();
                mark(&mut coverage[i], g.product.forward_start, g.template_length);
                assumptions.push(format!("{material_id}: unique exact PCR geometry; ideal new 5′OH/3′OH strands{}; no inherited chemical states",if tail.is_some(){" plus explicit ideal single-nucleotide 3′ tailing"}else{""}));
            }
        }
    }
    if substrates
        .strands
        .iter()
        .map(|s| s.sequence.len())
        .sum::<usize>()
        > 100_000
        || substrates.pairing.len() > 1024
    {
        return Err(Error::TooLarge);
    }
    let mut not_represented = vec![];
    for (i, c) in coverage.iter().enumerate() {
        let mut start = 0;
        while start < c.len() {
            if c[start] {
                start += 1;
                continue;
            }
            let mut end = start + 1;
            while end < c.len() && !c[end] {
                end += 1
            }
            if not_represented.len() >= 4096 {
                return Err(Error::TooLarge);
            }
            not_represented.push(ReferenceGap {
                source_role: inputs[i].role.clone(),
                start,
                end,
            });
            start = end;
        }
    }
    physical::preflight(&substrates, &r.steps, &r.engine_limits).map_err(engine_error)?;
    Ok(PreparationEvidence {
        substrates,
        maps,
        not_represented,
        assumptions,
        retained_steps: vec![],
        layout: vec![],
        graph_origins: vec![],
    })
}
pub(crate) fn prepare(inputs: &[Input], r: &Request) -> Result<PreparedContext, Error> {
    prepare_with_ancestry(inputs, r, None, None)
}
fn prepare_with_ancestry(
    inputs: &[Input],
    r: &Request,
    budget: Option<crate::verification::ancestry::TraceBudget>,
    replay: Option<&mut crate::verification::replay_budget::ReplayBudget<'_>>,
) -> Result<PreparedContext, Error> {
    prepare_bound(inputs, r, &[], budget, replay)
}
fn prepare_bound(
    inputs: &[Input],
    r: &Request,
    root_bindings: &[crate::cloning::topo::graph_history::SourceBinding],
    trace_budget: Option<crate::verification::ancestry::TraceBudget>,
    mut replay: Option<&mut crate::verification::replay_budget::ReplayBudget<'_>>,
) -> Result<PreparedContext, Error> {
    let expanded = expand(
        inputs,
        r,
        if root_bindings.is_empty() {
            None
        } else {
            Some(root_bindings)
        },
    )?;
    let roots = expanded
        .bindings
        .iter()
        .filter_map(|b| b.origin.as_ref())
        .map(|o| match o {
            crate::material::run::MaterialOrigin::Run(c) => {
                crate::material::capsule::NativeOriginRef::Run(c)
            }
            crate::material::run::MaterialOrigin::Derived(c) => {
                crate::material::capsule::NativeOriginRef::Derived(c)
            }
        });
    crate::material::capsule::native_origins_preflight(roots)?;
    let mut cache = std::collections::BTreeMap::new();
    let mut incomplete = false;
    let mut conditional_roles = BTreeSet::new();
    let mut findings = vec![];
    for b in &expanded.bindings {
        if let Some(origin) = &b.origin {
            let input = source(&expanded.inputs, &b.role)?;
            let key = hash(&(&b.event, origin.selection()))?;
            if !cache.contains_key(&key) {
                let fresh = if let Some(replay) = replay.as_deref_mut() {
                    origin.replay_scoped(input, replay, trace_budget.clone())?
                } else {
                    (origin.replay(input)?, None)
                };
                cache.insert(key.clone(), fresh);
            }
            let fresh = &cache[&key].0;
            incomplete |= fresh.incomplete();
            if fresh.incomplete() {
                conditional_roles.insert(b.role.clone());
            }
            if fresh.findings_total() > 0 {
                findings.push(format!("Source {}: showing at most one of {} inherited native findings; complete originating closure remains retained{}.",b.role,fresh.findings_total(),if fresh.findings_truncated(){"; original excerpt was bounded"}else{""}));
                findings.extend(
                    fresh
                        .findings()
                        .iter()
                        .take(1)
                        .map(|s| s.chars().take(2048).collect::<String>()),
                );
            }
        }
    }
    let mut templates = vec![];
    for b in &expanded.bindings {
        if let Some(origin) = &b.origin {
            let key = hash(&(&b.event, origin.selection()))?;
            templates.push((b.role.clone(), cache[&key].0.tracked()));
        }
    }
    // Root Step callers already proved their actual PCR-template eligibility in
    // source_preflight_reconstructed. Nested origins are freshly replayed here.
    for (role, tracked) in &templates {
        tracked.validate()?;
        let m = source(&expanded.inputs, role)?;
        if tracked.duplex.top != m.sequence
            || tracked.duplex.topology != m.topology
            || tracked.duplex.bottom_start != 0
            || tracked.duplex.bottom != tracked.duplex.top.reverse_complement()?
            || !tracked.duplex.nicks.is_empty()
        {
            return Err(bad(
                "saved graph PCR requires an actual intact aligned template; staggered/nicked models remain unsupported",
            ));
        }
    }
    let budget = match replay.as_deref_mut() {
        Some(replay) if replay.is_qc_job() => Some(&mut replay.pcr),
        _ => None,
    };
    let mut evidence = prepare_declared(&expanded.inputs, &expanded.request, budget)?;
    let prefix = expanded
        .request
        .steps
        .len()
        .checked_sub(r.steps.len())
        .ok_or(Error::TooLarge)?;
    evidence.retained_steps = expanded.request.steps[..prefix].to_vec();
    evidence.layout = expanded.layout;
    evidence.graph_origins = expanded.origins;
    let origin_source_incomplete=evidence.graph_origins.iter().map(|origin|evidence.maps.iter().any(|map|{
        origin.strand_start<=map.strand&&map.strand<origin.strand_end&&matches!(&map.origin,Origin::ExactPcrCopy{template_role,..}if conditional_roles.contains(template_role))
    })).collect();
    // These flags are compared with actual scoped execution before any result
    // is published; the retained record is never authority for completeness.
    incomplete |= evidence.graph_origins.iter().any(|o| o.incomplete);
    let mut native_origins = Vec::new();
    let mut origin_indices = std::collections::BTreeMap::new();
    for (key, (_, retained)) in cache {
        if let Some(retained) = retained {
            origin_indices.insert(key, native_origins.len());
            native_origins.push(retained);
        }
    }
    let source_origin_indices = expanded
        .source_origin_keys
        .into_iter()
        .map(|key| key.and_then(|key| origin_indices.get(&key).copied()))
        .collect();
    Ok(PreparedContext {
        inputs: expanded.inputs,
        source_scopes: expanded.source_scopes,
        source_origin_indices,
        native_origins,
        evidence,
        effective_steps: expanded.request.steps,
        incomplete,
        origin_source_incomplete,
        findings,
    })
}
pub(crate) fn effective_steps(
    evidence: &PreparationEvidence,
    r: &Request,
) -> Result<Vec<physical::Step>, Error> {
    let expected = retained_steps(r)?;
    if expected != evidence.retained_steps {
        return Err(bad(
            "retained graph transfer prefix differs from native request",
        ));
    }
    let mut all = expected;
    all.extend(r.steps.iter().cloned());
    if all.len() > 64 {
        return Err(Error::TooLarge);
    }
    Ok(all)
}
fn retained_steps(r: &Request) -> Result<Vec<physical::Step>, Error> {
    let roots = r.recipe.materials.iter().filter_map(|m| {
        if let MaterialRecipe::SavedGraph { origin, .. } = m {
            Some(crate::material::capsule::NativeOriginRef::TopoGraph(origin))
        } else {
            None
        }
    });
    crate::material::capsule::native_origins_preflight(roots)?;
    fn collect(
        r: &Request,
        base: usize,
        step_base: usize,
        depth: usize,
    ) -> Result<(usize, Vec<physical::Step>), Error> {
        if depth > 32 {
            return Err(Error::TooLarge);
        }
        let mut n = 0usize;
        let mut out = vec![];
        for material in &r.recipe.materials {
            match material {
                MaterialRecipe::Prepared { strands, .. } => {
                    n = n.checked_add(strands.len()).ok_or(Error::TooLarge)?
                }
                MaterialRecipe::ExactPcr { .. } => n = n.checked_add(2).ok_or(Error::TooLarge)?,
                MaterialRecipe::SavedGraph { origin, .. } => {
                    let record = origin.validate_shallow()?;
                    let (count, steps) = collect(
                        &record.request.request,
                        base.checked_add(n).ok_or(Error::TooLarge)?,
                        step_base.checked_add(out.len()).ok_or(Error::TooLarge)?,
                        depth + 1,
                    )?;
                    n = n.checked_add(count).ok_or(Error::TooLarge)?;
                    out.extend(steps);
                }
            }
            if n > 64 || out.len() > 64 {
                return Err(Error::TooLarge);
            }
        }
        let map = (0..n)
            .map(|i| base.checked_add(i).ok_or(Error::TooLarge))
            .collect::<Result<Vec<_>, _>>()?;
        for step in &r.steps {
            out.push(remap_step(step, &map, step_base)?);
            if out.len() > 64 {
                return Err(Error::TooLarge);
            }
        }
        Ok((n, out))
    }
    let (_, mut full) = collect(r, 0, 0, 0)?;
    if r.steps.len() > full.len() {
        return Err(Error::TooLarge);
    }
    full.truncate(full.len() - r.steps.len());
    Ok(full)
}
fn check_origin_incompleteness(
    prepared: &PreparedContext,
    actual: &physical::Report,
) -> Result<(), Error> {
    if prepared.origin_source_incomplete.len() != prepared.evidence.graph_origins.len() {
        return Err(invalid("source completeness scope changed"));
    }
    for (origin, source_incomplete) in prepared
        .evidence
        .graph_origins
        .iter()
        .zip(&prepared.origin_source_incomplete)
    {
        let mut computed = *source_incomplete;
        for index in &origin.step_map {
            let step = actual
                .steps
                .get(*index)
                .filter(|_| *index < origin.after_steps)
                .ok_or_else(|| invalid("source transfer completeness map changed"))?;
            computed |= step.outcome() != physical::Outcome::Applied;
        }
        if computed != origin.incomplete {
            return Err(invalid(
                "native graph incomplete claim differs from fresh scoped source/transfer execution",
            ));
        }
    }
    Ok(())
}
pub(crate) fn execute(
    prepared: &PreparedContext,
    limits: &physical::Limits,
) -> Result<physical::Executed, Error> {
    let checkpoints = prepared
        .evidence
        .graph_origins
        .iter()
        .map(|o| physical::Checkpoint {
            after_steps: o.after_steps,
            strand_start: o.strand_start,
            strand_end: o.strand_end,
            activation_unmap: o.activation_unmap.iter().cloned().collect(),
            step_map: o.step_map.clone(),
        })
        .collect::<Vec<_>>();
    let fresh = physical::simulate_with_checkpoints(
        &prepared.evidence.substrates,
        &prepared.effective_steps,
        limits,
        &checkpoints,
    )
    .map_err(engine_error)?;
    check_origin_incompleteness(prepared, fresh.report())?;
    for (origin, product) in prepared
        .evidence
        .graph_origins
        .iter()
        .zip(fresh.checkpoints())
    {
        if hash(&(crate::cloning::topo::graph_history::ENGINE_VERSION, product))?
            != origin.selection.state_revision
        {
            return Err(invalid(
                "native graph prefix differs from retained exact state",
            ));
        }
    }
    Ok(fresh)
}
pub(crate) fn simulate_bound(
    inputs: &[Input],
    r: &Request,
    bindings: &[crate::cloning::topo::graph_history::SourceBinding],
) -> Result<TopoReport, Error> {
    let roles = bindings
        .iter()
        .filter(|b| b.origin.is_some())
        .map(|b| b.role.clone())
        .collect::<BTreeSet<_>>();
    source_preflight_inner(inputs, r, &roles)?;
    let prepared = prepare_bound(inputs, r, bindings, None, None)?;
    report_prepared(inputs, r, prepared)
}
fn report_prepared(
    inputs: &[Input],
    r: &Request,
    prepared: PreparedContext,
) -> Result<TopoReport, Error> {
    let physical = execute(&prepared, &r.engine_limits)?.into_report();
    let result = TopoReport {
        engine_version: ENGINE_VERSION.into(),
        batch_revision: batch_binding(inputs, r)?,
        sources: identities(&prepared.inputs),
        incomplete: prepared.incomplete
            || physical.unsupported_or_unknown
            || !physical.all_requested_steps_applied,
        findings: prepared.findings,
        preparation: prepared.evidence,
        physical,
    };
    bounded(&result, MAX_REPORT_BYTES)?;
    Ok(result)
}
pub fn simulate(inputs: &[Input], r: &Request) -> Result<TopoReport, Error> {
    source_preflight(inputs, r)?;
    simulate_preflighted(inputs, r)
}
fn simulate_preflighted(inputs: &[Input], r: &Request) -> Result<TopoReport, Error> {
    report_prepared(inputs, r, prepare(inputs, r)?)
}
pub fn product_revision(inputs: &[Input], r: &Request, s: &Selection) -> Result<String, Error> {
    let batch = batch_revision(inputs, r)?;
    if s.batch_revision.len() != 64 || s.batch_revision != batch || s.component >= 64 {
        return Err(invalid("selection does not bind current source/recipe"));
    }
    hash(&(ENGINE_VERSION, "selected-physical-component", &batch, s))
}
/// Fresh source reconstruction; imported Product/Report values are never authority.
pub fn product(inputs: &[Input], r: &Request, s: &Selection) -> Result<Product, Error> {
    product_revision(inputs, r, s)?;
    product_preflighted(inputs, r, s)
}
fn product_preflighted(inputs: &[Input], r: &Request, s: &Selection) -> Result<Product, Error> {
    product_native(inputs, r, s, None, None).map(|(product, _)| product)
}
fn product_native(
    inputs: &[Input],
    r: &Request,
    s: &Selection,
    capture: Option<crate::verification::ancestry::TraceBudget>,
    replay: Option<&mut crate::verification::replay_budget::ReplayBudget<'_>>,
) -> Result<
    (
        Product,
        Option<crate::verification::projection::GraphProjections>,
    ),
    Error,
> {
    let batch = batch_binding(inputs, r)?;
    if s.batch_revision != batch || s.component >= 64 {
        return Err(invalid("selection does not bind current source/recipe"));
    }
    let batch_revision = hash(&(ENGINE_VERSION, "selected-physical-component", &batch, s))?;
    let prepared = prepare_with_ancestry(inputs, r, capture.clone(), replay)?;
    let fresh = execute(&prepared, &r.engine_limits)?;
    let select = projection::Selection {
        source_batch_revision: fresh.report().batch_revision.clone(),
        component: s.component,
        reference_seed: s.reference_seed,
    };
    let projection = projection::project_fresh(
        &prepared.evidence.substrates,
        &prepared.effective_steps,
        &r.engine_limits,
        &select,
        &r.projection_options,
        fresh,
    )
    .map_err(projection_error)?;
    let reference = if let projection::Outcome::Proven(v) = &projection.outcome {
        let mut m = Molecule::dna(
            "TOPO paired-contour view",
            v.top.as_str(),
            match v.topology {
                projection::ContourTopology::Linear => Topology::Linear,
                projection::ContourTopology::Circular => Topology::Circular,
            },
        )?;
        m.metadata
            .insert("source:legacy:chemistry_unmapped".into(), "true".into());
        m.metadata
            .insert("topo:physical_graph_required".into(), "true".into());
        m.metadata
            .insert("topo:adapter_batch_revision".into(), batch.clone());
        m.metadata.insert(
            RECIPE_KEY.into(),
            serde_json::to_string(&(&batch, s, identities(&prepared.inputs)))
                .map_err(|_| Error::TooLarge)?,
        );
        Some(m)
    } else {
        None
    };
    let p = Product {
        engine_version: ENGINE_VERSION.into(),
        batch_revision,
        selection: s.clone(),
        sources: identities(&prepared.inputs),
        preparation: prepared.evidence,
        projection,
        annotations: None,
        reference,
    };
    bounded(&p, MAX_REPORT_BYTES)?;
    let (product, receipts) = annotate_fresh_product(&prepared.inputs, p, capture.is_some())?;
    Ok((
        product,
        graph_projections(
            prepared.inputs,
            prepared.source_scopes,
            prepared.source_origin_indices,
            prepared.native_origins,
            receipts,
        ),
    ))
}
/// Bind source occurrences recorded at expansion time. No later role/content
/// lookup reconstructs which retained graph supplied a source.
fn graph_projections(
    prepared: Vec<Input>,
    scopes: Vec<PreparedSourceScope>,
    origin_indices: Vec<Option<usize>>,
    retained: Vec<crate::material::run::RetainedReplay>,
    receipts: Option<Vec<crate::verification::projection::SourceProjection>>,
) -> Option<crate::verification::projection::GraphProjections> {
    use crate::verification::projection::{GraphProjectionSource, GraphProjections};
    let projections = receipts?;
    if prepared.len() != scopes.len()
        || prepared.len() != origin_indices.len()
        || crate::encoding::bounded_json::preflight(&prepared, 8 * 1024 * 1024).is_err()
    {
        return None;
    }
    let sources = prepared
        .into_iter()
        .zip(scopes)
        .zip(origin_indices)
        .map(|((input, scope), retained_index)| match scope {
            PreparedSourceScope::CurrentInput { index } => {
                GraphProjectionSource::CurrentInput { index }
            }
            PreparedSourceScope::RetainedGraph { role, scopes } => {
                GraphProjectionSource::RetainedGraph {
                    role,
                    scopes,
                    molecule: input.molecule,
                    retained_index,
                }
            }
        })
        .collect();
    Some(GraphProjections {
        sources,
        projections,
        retained,
    })
}
/// Annotation correspondence is attached to the same fresh canonical Product.
/// It never supplies native graph/terminal authority.
fn annotate_fresh_product(
    inputs: &[Input],
    mut product: Product,
    capture: bool,
) -> Result<
    (
        Product,
        Option<Vec<crate::verification::projection::SourceProjection>>,
    ),
    Error,
> {
    let mut receipts = None;
    if let (projection::Outcome::Proven(view), Some(reference)) =
        (&product.projection.outcome, &mut product.reference)
    {
        let (evidence, captured) = crate::cloning::topo::annotations::project_native(
            inputs,
            &product.preparation,
            view,
            reference,
            capture,
        )?;
        receipts = captured;
        reference.features = evidence.features.clone();
        reference.metadata.insert(
            "topo:annotation_engine".into(),
            ANNOTATION_ENGINE_VERSION.into(),
        );
        reference.metadata.insert(
            "topo:annotation_revision".into(),
            product.batch_revision.clone(),
        );
        reference.metadata.insert("topo:annotation_scope".into(),"Source-coordinate correspondence only; no inherited chemistry or biological validity claim.".into());
        reference.validate()?;
        product.annotations = Some(evidence);
    }
    bounded(&product, MAX_REPORT_BYTES)?;
    Ok((product, receipts))
}
pub fn target_revision(
    inputs: &[Input],
    r: &Request,
    s: &Selection,
    target: &Molecule,
    policy: &AcceptancePolicy,
) -> Result<String, Error> {
    target_preflight(target)?;
    hash(&(
        ENGINE_VERSION,
        "independent-target",
        product_revision(inputs, r, s)?,
        target.revision(),
        policy,
    ))
}
fn target_preflight(target: &Molecule) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[target])?;
    if target.sequence.len() > 100_000 || target.sequence.alphabet() != Alphabet::Dna {
        return Err(Error::TooLarge);
    }
    target.validate()
}
pub fn check_target(
    inputs: &[Input],
    r: &Request,
    s: &Selection,
    target: &Molecule,
    policy: &AcceptancePolicy,
) -> Result<TargetAssessment, Error> {
    source_preflight(inputs, r)?;
    target_preflight(target)?;
    let product = product(inputs, r, s)?;
    check_fresh_product(product, target, policy)
}
fn check_fresh_product(
    product: Product,
    target: &Molecule,
    policy: &AcceptancePolicy,
) -> Result<TargetAssessment, Error> {
    target_preflight(target)?;
    let target_revision = target.revision();
    let batch_revision = hash(&(
        ENGINE_VERSION,
        "independent-target",
        &product.batch_revision,
        &target_revision,
        policy,
    ))?;
    let outcome=match (&product.projection.outcome,&product.reference){(projection::Outcome::Proven(v),Some(m))=>{let mut reasons=vec![];if !policy.allow_nicks&&!v.nicks.is_empty(){reasons.push("selected contour retains nicks".into());}
if !policy.allow_unknown_terminal_chemistry&&!v.fully_known_terminal_chemistry{reasons.push("terminal chemistry is unknown".into());}
if !policy.allow_unapplied_steps&&!v.all_requested_steps_applied{reasons.push("requested step did not apply".into());}
if !policy.allow_other_components&&product.projection.physical.product.components.len()!=1{reasons.push("other physical components remain".into());}let dna=target.sequence.as_str();let same=|candidate:&str| candidate==dna||m.topology==Topology::Circular&&policy.allow_circular_rotation&&analysis::rotation_offset(candidate,dna).is_some();let sequence_matches=same(m.sequence.as_str())||policy.allow_reverse_complement&&same(m.sequence.reverse_complement()?.as_str());if m.topology!=target.topology||!sequence_matches{reasons.push("independent target DNA/topology differs under declared comparison policy".into());}
if reasons.is_empty(){TargetOutcome::AcceptedUnderDeclaredGraph}else{TargetOutcome::Rejected{reasons}}},_=>TargetOutcome::Withheld{reason:"No unique paired contour; full physical graph and withheld projection evidence retained. No flap/nick repair or target-derived path is attempted.".into()}};
    let result = TargetAssessment {
        engine_version: ENGINE_VERSION.into(),
        batch_revision,
        target_revision,
        policy: policy.clone(),
        product,
        outcome,
    };
    bounded(&result, MAX_REPORT_BYTES)?;
    Ok(result)
}
/// Heavy native audit: full generation/transfer/projection replay, not a cheap UI
/// validator. Serialized claims cannot provide physical input state.
pub fn replay_product(
    inputs: &[Input],
    r: &Request,
    s: &Selection,
    claim: &Product,
) -> Result<(), Error> {
    let fresh = product(inputs, r, s)?;
    // Native callers may construct arbitrarily nested feature trees; guard them
    // before bounded serialization, rather than relying on a JSON decoder limit.
    if let Some(m) = &claim.reference {
        collagene_model::bounds::model_preflight(&[m])?;
        if m.sequence.len() > 100_000 {
            return Err(Error::TooLarge);
        }
    }
    bounded(claim, MAX_REPORT_BYTES)?;
    if hash(&fresh)? != hash(claim)? {
        return Err(invalid("product differs from fresh reconstruction"));
    }
    Ok(())
}

/// Both declared vector/insert orientations, retaining the bounded engine's full
/// component graphs. This is not an exhaustive chemical reaction mixture.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct OrientationInspection {
    pub engine_version: String,
    pub batch_revision: String,
    pub sources: Vec<SourceIdentity>,
    pub preparation: PreparationEvidence,
    pub reference: physical::OrientationCandidates,
    /// Exact requests to use with fresh product materialization, in matching order.
    pub candidate_requests: Vec<Request>,
}
pub fn orientations_revision(
    inputs: &[Input],
    r: &Request,
    vector: [usize; 2],
    insert: [usize; 2],
    directional: bool,
) -> Result<String, Error> {
    source_preflight(inputs, r)?;
    if !r.steps.is_empty()
        || vector.iter().chain(&insert).any(|i| *i >= 64)
        || vector
            .iter()
            .chain(&insert)
            .copied()
            .collect::<BTreeSet<_>>()
            .len()
            != 4
    {
        return Err(bad(
            "orientation inspection needs four distinct bounded native strand handles and no new transfers",
        ));
    }
    hash(&(
        ENGINE_VERSION,
        "two-declared-orientations",
        identities(inputs),
        r,
        vector,
        insert,
        directional,
    ))
}
pub fn inspect_orientations(
    inputs: &[Input],
    r: &Request,
    vector: [usize; 2],
    insert: [usize; 2],
    directional: bool,
) -> Result<OrientationInspection, Error> {
    source_preflight(inputs, r)?;
    orientations_preflighted(inputs, r, vector, insert, directional, None)
}
fn orientations_preflighted(
    inputs: &[Input],
    r: &Request,
    vector: [usize; 2],
    insert: [usize; 2],
    directional: bool,
    replay: Option<&mut crate::verification::replay_budget::ReplayBudget<'_>>,
) -> Result<OrientationInspection, Error> {
    if !r.steps.is_empty() {
        return Err(bad(
            "orientation inspection requires an empty transfer list; orientation-specific requests are returned explicitly",
        ));
    }
    let prepared = prepare_with_ancestry(inputs, r, None, replay)?;
    if vector
        .iter()
        .chain(&insert)
        .any(|i| *i >= prepared.evidence.substrates.strands.len())
    {
        return Err(invalid("orientation source outside actual graph layout"));
    }
    let checkpoints = prepared
        .evidence
        .graph_origins
        .iter()
        .map(|o| physical::Checkpoint {
            after_steps: o.after_steps,
            strand_start: o.strand_start,
            strand_end: o.strand_end,
            activation_unmap: o.activation_unmap.iter().cloned().collect(),
            step_map: o.step_map.clone(),
        })
        .collect::<Vec<_>>();
    let reference = physical::orientations(
        &prepared.evidence.substrates,
        &prepared.evidence.retained_steps,
        vector,
        insert,
        directional,
        &r.engine_limits,
        &checkpoints,
    )
    .map_err(engine_error)?;
    let actual = reference
        .candidates
        .first()
        .ok_or_else(|| bad("native orientation prefix evidence missing"))?;
    check_origin_incompleteness(&prepared, actual)?;
    for (origin, p) in prepared
        .evidence
        .graph_origins
        .iter()
        .zip(&reference.checkpoints)
    {
        if hash(&(crate::cloning::topo::graph_history::ENGINE_VERSION, p))?
            != origin.selection.state_revision
        {
            return Err(bad(
                "orientation graph source differs from retained checkpoint",
            ));
        }
    }
    let sources = identities(&prepared.inputs);
    let preparation = prepared.evidence;
    let mut candidate_requests = Vec::new();
    for pending in &reference.pending_steps {
        let mut request = r.clone();
        request.steps = pending.clone();
        candidate_requests.push(request);
    }
    let batch_revision = hash(&(
        ENGINE_VERSION,
        "two-declared-orientations",
        identities(inputs),
        r,
        vector,
        insert,
        directional,
    ))?;
    let result = OrientationInspection {
        engine_version: ENGINE_VERSION.into(),
        batch_revision,
        sources,
        preparation,
        reference,
        candidate_requests,
    };
    bounded(&result, MAX_REPORT_BYTES)?;
    Ok(result)
}

// Only the native saved-run adapter calls this entry after obtaining real tracked
// state. All shared implementations remain private; public standalone APIs retain
// their original source guards.
pub(crate) fn inspect_reconstructed_task(
    inputs: &[Input],
    templates: &[(String, &crate::chemistry::tracked_duplex::TrackedDuplex)],
    task: &crate::cloning::topo::run_preparation::Task,
) -> Result<crate::cloning::topo::run_preparation::Reply, Error> {
    inspect_reconstructed_task_native(inputs, templates, task, None, None).map(|(reply, _)| reply)
}
pub(crate) fn inspect_reconstructed_task_with_projections(
    inputs: &[Input],
    templates: &[(String, &crate::chemistry::tracked_duplex::TrackedDuplex)],
    task: &crate::cloning::topo::run_preparation::Task,
    capture: Option<crate::verification::ancestry::TraceBudget>,
    replay: &mut crate::verification::replay_budget::ReplayBudget<'_>,
) -> Result<
    (
        crate::cloning::topo::run_preparation::Reply,
        Option<crate::verification::projection::GraphProjections>,
    ),
    Error,
> {
    inspect_reconstructed_task_native(inputs, templates, task, capture, Some(replay))
}
fn inspect_reconstructed_task_native(
    inputs: &[Input],
    templates: &[(String, &crate::chemistry::tracked_duplex::TrackedDuplex)],
    task: &crate::cloning::topo::run_preparation::Task,
    capture: Option<crate::verification::ancestry::TraceBudget>,
    mut replay: Option<&mut crate::verification::replay_budget::ReplayBudget<'_>>,
) -> Result<
    (
        crate::cloning::topo::run_preparation::Reply,
        Option<crate::verification::projection::GraphProjections>,
    ),
    Error,
> {
    use crate::cloning::topo::run_preparation::{Reply, Task};
    source_preflight_reconstructed(inputs, task.preparation(), templates)?;
    let mut receipts = None;
    let reply = match task {
        Task::Simulate { request } => {
            let prepared =
                prepare_with_ancestry(inputs, request, capture.clone(), replay.as_deref_mut())?;
            Reply::Simulated(report_prepared(inputs, request, prepared)?)
        }
        Task::Orientations {
            request,
            vector,
            insert,
            directional,
        } => Reply::Orientations(orientations_preflighted(
            inputs,
            request,
            *vector,
            *insert,
            *directional,
            replay.as_deref_mut(),
        )?),
        Task::Project { request, selection } => {
            let (product, captured) = product_native(
                inputs,
                request,
                selection,
                capture.clone(),
                replay.as_deref_mut(),
            )?;
            receipts = captured;
            Reply::Projected(product)
        }
        Task::Check {
            request,
            selection,
            target,
            policy,
        } => {
            let (fresh, captured) = product_native(
                inputs,
                request,
                selection,
                capture.clone(),
                replay.as_deref_mut(),
            )?;
            receipts = captured;
            let TargetAssessment {
                product,
                policy,
                target_revision,
                outcome,
                ..
            } = check_fresh_product(fresh, target, policy)?;
            Reply::Checked(TargetAssessment {
                engine_version: ENGINE_VERSION.into(),
                batch_revision: hash(&(
                    ENGINE_VERSION,
                    "independent-target",
                    &product.batch_revision,
                    &target_revision,
                    &policy,
                ))?,
                target_revision,
                policy,
                product,
                outcome,
            })
        }
        Task::ProcessActions { request } => {
            let report = if replay.is_some() {
                let (report, captured) =
                    crate::cloning::topo::actions::inspect_reconstructed_with_projections(
                        inputs,
                        request,
                        templates,
                        capture.clone(),
                        replay.as_deref_mut().expect("scoped native entry"),
                    )?;
                receipts = captured;
                report
            } else {
                crate::cloning::topo::actions::inspect_reconstructed(inputs, request, templates)?
            };
            Reply::ActionsProcessed(report)
        }
    };
    Ok((reply, receipts))
}
// Bounds/declaration check only. Never returns tracked material or bypasses native
// replay; caller must independently bind the supplied capsule and exact event.
pub(crate) fn declared_template_preflight(
    inputs: &[Input],
    r: &Request,
    role: &str,
) -> Result<(), Error> {
    text(role, 256)?;
    let mut roles = BTreeSet::new();
    roles.insert(role.to_string());
    source_preflight_inner(inputs, r, &roles)?;
    if inputs.iter().filter(|i| i.role == role).count() != 1 {
        return Err(invalid("saved template role must occur exactly once"));
    }
    let mut used = false;
    for material in &r.recipe.materials {
        match material {
            MaterialRecipe::SavedGraph { .. } => {}
            MaterialRecipe::Prepared { strands, .. }
                if strands
                    .iter()
                    .any(|s| s.segments.iter().any(|p| p.source_role == role)) =>
            {
                return Err(bad(
                    "saved template may only feed explicit new PCR; prepared vector ancestry is unsupported",
                ));
            }
            MaterialRecipe::ExactPcr { template_role, .. } if template_role == role => used = true,
            _ => (),
        }
    }
    if !used {
        return Err(invalid("saved template must feed explicit exact PCR"));
    }
    Ok(())
}

pub mod actions;
pub mod annotations;
pub mod graph_history;
pub mod preparation;
pub mod run_preparation;
pub mod worker;
