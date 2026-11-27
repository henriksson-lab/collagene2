//! Source-backed preparation declarations and fresh physical TOPO reconstruction.
//! No sequence-record metadata establishes activation, end chemistry or kit state.
use crate::{Alphabet, Error, Molecule, Sequence, Topology, analysis, chemistry, tailing};
pub use collagene_topo as physical;
use physical::projection;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const ENGINE_VERSION: &str = "topo-source-adapter-2";
pub const ANNOTATION_ENGINE_VERSION: &str = "topo-annotated-reference-1";
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
    SavedGraph {material_id:String,origin:crate::topo_graph_history::Capsule},
}}
dto! {pub struct GenerationRecipe {pub schema_version:u16,pub materials:Vec<MaterialRecipe>}}
dto! {pub struct Request {pub recipe:GenerationRecipe,pub transfers:Vec<physical::Transfer>,pub engine_limits:physical::Limits,pub projection_options:projection::Options}}
dto! {pub enum Origin {
    ReferenceCorrespondence(SourceInterval),
    /// Full freshly computed PCR strand before optional tailing. This identifies
    /// an ideal polymerized copy; it asserts neither inherited atoms nor methylation.
    ExactPcrCopy {material:usize,template_role:String,template_start:usize,template_length:usize,forward_tail_length:usize,reverse_tail_length:usize,reverse:bool},
    SynthesizedTail {material:usize,base:tailing::TailingBase},
}}
dto! {pub struct PreparedMap {pub strand:usize,pub strand_start:usize,pub length:usize,pub origin:Origin}}
dto! {pub struct ReferenceGap {pub source_role:String,pub start:usize,pub end:usize}}
dto! {pub struct PreparationEvidence {pub substrates:physical::Substrates,pub maps:Vec<PreparedMap>,pub not_represented:Vec<ReferenceGap>,pub assumptions:Vec<String>,pub retained_transfers:Vec<physical::Transfer>,pub layout:Vec<MaterialLayout>,pub graph_origins:Vec<GraphOriginEvidence>}}
dto! {pub struct TopoReport {pub engine_version:String,pub batch_revision:String,pub sources:Vec<SourceIdentity>,pub preparation:PreparationEvidence,pub physical:physical::Report,pub incomplete:bool,pub findings:Vec<String>}}
dto! {pub struct Selection {pub batch_revision:String,pub component:usize,pub reference_seed:physical::Call}}
dto! {pub struct Product {pub engine_version:String,pub batch_revision:String,pub selection:Selection,pub sources:Vec<SourceIdentity>,pub preparation:PreparationEvidence,pub projection:projection::ProjectionReport,pub annotations:Option<crate::topo_annotations::AnnotationEvidence>,pub reference:Option<Molecule>}}
dto! {#[derive(Default)] pub struct AcceptancePolicy {pub allow_nicks:bool,pub allow_unknown_terminal_chemistry:bool,pub allow_unapplied_transfers:bool,pub allow_other_components:bool,pub allow_circular_rotation:bool,pub allow_reverse_complement:bool}}
dto! {pub enum TargetOutcome {AcceptedUnderDeclaredGraph,Rejected{reasons:Vec<String>},Withheld{reason:String}}}
dto! {pub struct TargetAssessment {pub engine_version:String,pub batch_revision:String,pub target_revision:String,pub policy:AcceptancePolicy,pub product:Product,pub outcome:TargetOutcome}}
// Whole-graph evidence and a selected annotated contour answer distinct current
// questions. There is one Product document, with documentary annotations only.
dto! {pub struct MaterialLayout{pub owner:Option<String>,pub material:usize,pub source_strands:Vec<usize>,pub source_roles:Vec<(String,String)>}}
dto! {pub struct GraphOriginEvidence{pub selection:crate::topo_graph_history::Selection,pub strand_start:usize,pub strand_end:usize,pub after_transfers:usize,pub activation_unmap:Vec<(String,String)>,pub transfer_map:Vec<usize>,pub incomplete:bool}}
pub(crate) struct PreparedContext {
    pub inputs: Vec<Input>,
    pub evidence: PreparationEvidence,
    pub effective_transfers: Vec<physical::Transfer>,
    pub incomplete: bool,
    pub origin_source_incomplete: Vec<bool>,
    pub findings: Vec<String>,
}

fn bad(s: &str) -> Error {
    Error::Unsupported(format!("TOPO source adapter: {s}"))
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
fn engine_error(e: impl std::fmt::Display) -> Error {
    bad(&e.to_string())
}
fn source<'a>(inputs: &'a [Input], role: &str) -> Result<&'a Molecule, Error> {
    inputs
        .iter()
        .find(|x| x.role == role)
        .map(|x| &x.molecule)
        .ok_or_else(|| bad("unknown source role"))
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
        return Err(bad("invalid source interval"));
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
    let roots = r.recipe.materials.iter().filter_map(|m| {
        if let MaterialRecipe::SavedGraph { origin, .. } = m {
            Some(crate::run_material::NativeOriginRef::TopoGraph(origin))
        } else {
            None
        }
    });
    crate::run_material::native_origins_preflight(roots)?;
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
    if inputs.len() > 64
        || r.recipe.schema_version != 1
        || r.recipe.materials.is_empty()
        || r.recipe.materials.len() > 64
        || r.transfers.len() > 16
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
            return Err(bad("unique DNA source roles required"));
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
                    if !ids.insert(s.id.as_str())
                        || s.segments.is_empty()
                        || s.segments.len() > 1024
                    {
                        return Err(bad(
                            "unique strand ids and bounded nonempty segments required",
                        ));
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
                            return Err(bad(
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
                        return Err(bad("prepared pairing indices are material-local"));
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
            Some(crate::run_material::NativeOriginRef::TopoGraph(origin))
        } else {
            None
        }
    });
    crate::run_material::native_origins_preflight(roots)?;
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
    templates: &[(String, &crate::modified_duplex::TrackedDuplex)],
) -> Result<(), Error> {
    if inputs.len() > 32 || templates.len() > 32 || templates.len() > inputs.len() {
        return Err(Error::TooLarge);
    }
    let mut roles = BTreeSet::new();
    for (role, _) in templates {
        text(role, 256)?;
        if !roles.insert(role.clone()) || inputs.iter().filter(|i| i.role == *role).count() != 1 {
            return Err(bad(
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
            return Err(bad(
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
            return Err(bad(
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
                || key == crate::recombination::CHEMISTRY_UNMAPPED_KEY
        }) {
            return Err(bad(
                "reference-only, variant or unmapped reference ancestry cannot establish a physical PCR template",
            ));
        }
    }
    Ok(())
}
pub(crate) fn batch_revision_reconstructed(
    inputs: &[Input],
    r: &Request,
    templates: &[(String, &crate::modified_duplex::TrackedDuplex)],
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
    request: Request,
    bindings: Vec<crate::topo_graph_history::SourceBinding>,
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
    root_bindings: Option<&[crate::topo_graph_history::SourceBinding]>,
) -> Result<Expanded, Error> {
    let roots = r.recipe.materials.iter().filter_map(|m| {
        if let MaterialRecipe::SavedGraph { origin, .. } = m {
            Some(crate::run_material::NativeOriginRef::TopoGraph(origin))
        } else {
            None
        }
    });
    crate::run_material::native_origins_preflight(roots)?;
    let mut out = Expanded {
        inputs: vec![],
        request: Request {
            recipe: GenerationRecipe {
                schema_version: 1,
                materials: vec![],
            },
            transfers: vec![],
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
    bindings: &[crate::topo_graph_history::SourceBinding],
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
    bindings: Option<&[crate::topo_graph_history::SourceBinding]>,
    seen: &mut BTreeSet<String>,
    used: &mut BTreeSet<String>,
    counter: &mut usize,
    out: &mut Expanded,
) -> Result<Vec<usize>, Error> {
    if *counter > 128 {
        return Err(Error::TooLarge);
    }
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
        out.inputs.push(Input {
            role: role.clone(),
            molecule: input.molecule.clone(),
        });
        if let Some(binding) = bindings.and_then(|b| b.iter().find(|b| b.role == input.role)) {
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
                if record.engine_version != crate::topo_graph_history::ENGINE_VERSION
                    || record.scientific_engine != crate::verification::ENGINE_VERSION
                {
                    return Err(bad(
                        "historical graph model cannot provide current material",
                    ));
                }
                let transfers_before = out.request.transfers.len();
                let activation_namespace = origin.selection.record_hash.clone();
                let indices = expand_node(
                    &origin.source_snapshots,
                    &record.request.request,
                    Some(&activation_namespace),
                    Some(&record.request.bindings),
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
                let transfer_map =
                    (transfers_before..out.request.transfers.len()).collect::<Vec<_>>();
                out.origins.push(GraphOriginEvidence {
                    selection: origin.selection.clone(),
                    strand_start: start,
                    strand_end: end,
                    after_transfers: out.request.transfers.len(),
                    activation_unmap,
                    transfer_map,
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
    for t in &r.transfers {
        out.request.transfers.push(physical::Transfer {
            donor: *strand_map
                .get(t.donor)
                .ok_or_else(|| bad("transfer donor outside original graph namespace"))?,
            acceptor: *strand_map
                .get(t.acceptor)
                .ok_or_else(|| bad("transfer acceptor outside original graph namespace"))?,
            geometry: t.geometry,
        });
    }
    if out.inputs.len() > 64
        || out.request.recipe.materials.len() > 64
        || out.request.transfers.len() > 64
    {
        return Err(Error::TooLarge);
    }
    Ok(strand_map)
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
fn prepare_declared(inputs: &[Input], r: &Request) -> Result<PreparationEvidence, Error> {
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
                    f += analysis::find_matches(template, canonical.annealing.as_str(), false)?
                        .len();
                    b += analysis::find_matches(
                        template,
                        canonical.annealing.reverse_complement()?.as_str(),
                        false,
                    )?
                    .len();
                }
                let upper = f
                    .checked_mul(b)
                    .and_then(|x| x.checked_mul(template.sequence.len() + 4000))
                    .filter(|x| *x <= 1_000_000)
                    .ok_or(Error::TooLarge)?;
                let _ = upper;
                let possibilities = analysis::simulate_pcr_geometry(template, forward, reverse)?;
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
    physical::preflight(&substrates, &r.transfers, &r.engine_limits).map_err(engine_error)?;
    Ok(PreparationEvidence {
        substrates,
        maps,
        not_represented,
        assumptions,
        retained_transfers: vec![],
        layout: vec![],
        graph_origins: vec![],
    })
}
pub(crate) fn prepare(inputs: &[Input], r: &Request) -> Result<PreparedContext, Error> {
    prepare_bound(inputs, r, &[])
}
fn prepare_bound(
    inputs: &[Input],
    r: &Request,
    root_bindings: &[crate::topo_graph_history::SourceBinding],
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
            crate::run_material::MaterialOrigin::Run(c) => {
                crate::run_material::NativeOriginRef::Run(c)
            }
            crate::run_material::MaterialOrigin::Derived(c) => {
                crate::run_material::NativeOriginRef::Derived(c)
            }
        });
    crate::run_material::native_origins_preflight(roots)?;
    let mut cache = std::collections::BTreeMap::new();
    let mut incomplete = false;
    let mut conditional_roles = BTreeSet::new();
    let mut findings = vec![];
    for b in &expanded.bindings {
        if let Some(origin) = &b.origin {
            let input = source(&expanded.inputs, &b.role)?;
            let key = hash(&(&b.event, origin.selection()))?;
            if !cache.contains_key(&key) {
                cache.insert(key.clone(), origin.replay(input)?);
            }
            let fresh = &cache[&key];
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
            templates.push((b.role.clone(), cache[&key].tracked()));
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
    let mut evidence = prepare_declared(&expanded.inputs, &expanded.request)?;
    let prefix = expanded
        .request
        .transfers
        .len()
        .checked_sub(r.transfers.len())
        .ok_or(Error::TooLarge)?;
    evidence.retained_transfers = expanded.request.transfers[..prefix].to_vec();
    evidence.layout = expanded.layout;
    evidence.graph_origins = expanded.origins;
    let origin_source_incomplete=evidence.graph_origins.iter().map(|origin|evidence.maps.iter().any(|map|{
        origin.strand_start<=map.strand&&map.strand<origin.strand_end&&matches!(&map.origin,Origin::ExactPcrCopy{template_role,..}if conditional_roles.contains(template_role))
    })).collect();
    // These flags are compared with actual scoped execution before any result
    // is published; the retained record is never authority for completeness.
    incomplete |= evidence.graph_origins.iter().any(|o| o.incomplete);
    Ok(PreparedContext {
        inputs: expanded.inputs,
        evidence,
        effective_transfers: expanded.request.transfers,
        incomplete,
        origin_source_incomplete,
        findings,
    })
}
pub(crate) fn effective_transfers(
    evidence: &PreparationEvidence,
    r: &Request,
) -> Result<Vec<physical::Transfer>, Error> {
    let expected = retained_transfers(r)?;
    if expected != evidence.retained_transfers {
        return Err(bad(
            "retained graph transfer prefix differs from native request",
        ));
    }
    let mut all = expected;
    all.extend(r.transfers.iter().cloned());
    if all.len() > 64 {
        return Err(Error::TooLarge);
    }
    Ok(all)
}
fn retained_transfers(r: &Request) -> Result<Vec<physical::Transfer>, Error> {
    let roots = r.recipe.materials.iter().filter_map(|m| {
        if let MaterialRecipe::SavedGraph { origin, .. } = m {
            Some(crate::run_material::NativeOriginRef::TopoGraph(origin))
        } else {
            None
        }
    });
    crate::run_material::native_origins_preflight(roots)?;
    fn collect(
        r: &Request,
        base: usize,
        depth: usize,
    ) -> Result<(usize, Vec<physical::Transfer>), Error> {
        if depth > 32 {
            return Err(Error::TooLarge);
        }
        let mut n = 0usize;
        let mut out = vec![];
        for m in &r.recipe.materials {
            match m {
                MaterialRecipe::Prepared { strands, .. } => n += strands.len(),
                MaterialRecipe::ExactPcr { .. } => n += 2,
                MaterialRecipe::SavedGraph { origin, .. } => {
                    let record = origin.validate_shallow()?;
                    let (count, steps) = collect(&record.request.request, base + n, depth + 1)?;
                    n = n.checked_add(count).ok_or(Error::TooLarge)?;
                    out.extend(steps);
                }
            }
            if n > 64 || out.len() > 64 {
                return Err(Error::TooLarge);
            }
        }
        for t in &r.transfers {
            out.push(physical::Transfer {
                donor: base.checked_add(t.donor).ok_or(Error::TooLarge)?,
                acceptor: base.checked_add(t.acceptor).ok_or(Error::TooLarge)?,
                geometry: t.geometry,
            });
        }
        Ok((n, out))
    }
    let (_, mut full) = collect(r, 0, 0)?;
    if r.transfers.len() > full.len() {
        return Err(Error::TooLarge);
    }
    full.truncate(full.len() - r.transfers.len());
    Ok(full)
}
fn check_origin_incompleteness(
    prepared: &PreparedContext,
    actual: &physical::Report,
) -> Result<(), Error> {
    if prepared.origin_source_incomplete.len() != prepared.evidence.graph_origins.len() {
        return Err(bad("source completeness scope changed"));
    }
    for (origin, source_incomplete) in prepared
        .evidence
        .graph_origins
        .iter()
        .zip(&prepared.origin_source_incomplete)
    {
        let mut computed = *source_incomplete;
        for index in &origin.transfer_map {
            let step = actual
                .steps
                .get(*index)
                .filter(|_| *index < origin.after_transfers)
                .ok_or_else(|| bad("source transfer completeness map changed"))?;
            computed |= step.outcome != physical::Outcome::Applied;
        }
        if computed != origin.incomplete {
            return Err(bad(
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
            after_transfers: o.after_transfers,
            strand_start: o.strand_start,
            strand_end: o.strand_end,
            activation_unmap: o.activation_unmap.iter().cloned().collect(),
            transfer_map: o.transfer_map.clone(),
        })
        .collect::<Vec<_>>();
    let fresh = physical::simulate_with_checkpoints(
        &prepared.evidence.substrates,
        &prepared.effective_transfers,
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
        if hash(&(crate::topo_graph_history::ENGINE_VERSION, product))?
            != origin.selection.state_revision
        {
            return Err(bad("native graph prefix differs from retained exact state"));
        }
    }
    Ok(fresh)
}
pub(crate) fn simulate_bound(
    inputs: &[Input],
    r: &Request,
    bindings: &[crate::topo_graph_history::SourceBinding],
) -> Result<TopoReport, Error> {
    let roles = bindings
        .iter()
        .filter(|b| b.origin.is_some())
        .map(|b| b.role.clone())
        .collect::<BTreeSet<_>>();
    source_preflight_inner(inputs, r, &roles)?;
    let prepared = prepare_bound(inputs, r, bindings)?;
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
            || !physical.all_requested_transfers_applied,
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
        return Err(bad("selection does not bind current source/recipe"));
    }
    hash(&(ENGINE_VERSION, "selected-physical-component", &batch, s))
}
/// Fresh source reconstruction; imported Product/Report values are never authority.
pub fn product(inputs: &[Input], r: &Request, s: &Selection) -> Result<Product, Error> {
    product_revision(inputs, r, s)?;
    product_preflighted(inputs, r, s)
}
fn product_preflighted(inputs: &[Input], r: &Request, s: &Selection) -> Result<Product, Error> {
    let batch = batch_binding(inputs, r)?;
    if s.batch_revision != batch || s.component >= 64 {
        return Err(bad("selection does not bind current source/recipe"));
    }
    let batch_revision = hash(&(ENGINE_VERSION, "selected-physical-component", &batch, s))?;
    let prepared = prepare(inputs, r)?;
    let fresh = execute(&prepared, &r.engine_limits)?;
    let select = projection::Selection {
        source_batch_revision: fresh.report().batch_revision.clone(),
        component: s.component,
        reference_seed: s.reference_seed,
    };
    let projection = projection::project_fresh(
        &prepared.evidence.substrates,
        &prepared.effective_transfers,
        &r.engine_limits,
        &select,
        &r.projection_options,
        fresh,
    )
    .map_err(engine_error)?;
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
    annotate_fresh_product(&prepared.inputs, p)
}
/// Annotation correspondence is attached to the same fresh canonical Product.
/// It never supplies native graph/terminal authority.
fn annotate_fresh_product(inputs: &[Input], mut product: Product) -> Result<Product, Error> {
    if let (projection::Outcome::Proven(view), Some(reference)) =
        (&product.projection.outcome, &mut product.reference)
    {
        let evidence =
            crate::topo_annotations::project(inputs, &product.preparation, view, reference)?;
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
    Ok(product)
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
if !policy.allow_unapplied_transfers&&!v.all_requested_transfers_applied{reasons.push("requested transfer did not apply".into());}
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
        return Err(bad("product differs from fresh reconstruction"));
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
    if !r.transfers.is_empty()
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
    orientations_preflighted(inputs, r, vector, insert, directional)
}
fn orientations_preflighted(
    inputs: &[Input],
    r: &Request,
    vector: [usize; 2],
    insert: [usize; 2],
    directional: bool,
) -> Result<OrientationInspection, Error> {
    if !r.transfers.is_empty() {
        return Err(bad(
            "orientation inspection requires an empty transfer list; orientation-specific requests are returned explicitly",
        ));
    }
    let prepared = prepare(inputs, r)?;
    if vector
        .iter()
        .chain(&insert)
        .any(|i| *i >= prepared.evidence.substrates.strands.len())
    {
        return Err(bad("orientation source outside actual graph layout"));
    }
    let checkpoints = prepared
        .evidence
        .graph_origins
        .iter()
        .map(|o| physical::Checkpoint {
            after_transfers: o.after_transfers,
            strand_start: o.strand_start,
            strand_end: o.strand_end,
            activation_unmap: o.activation_unmap.iter().cloned().collect(),
            transfer_map: o.transfer_map.clone(),
        })
        .collect::<Vec<_>>();
    let reference = physical::orientations(
        &prepared.evidence.substrates,
        &prepared.evidence.retained_transfers,
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
        if hash(&(crate::topo_graph_history::ENGINE_VERSION, p))? != origin.selection.state_revision
        {
            return Err(bad(
                "orientation graph source differs from retained checkpoint",
            ));
        }
    }
    let sources = identities(&prepared.inputs);
    let preparation = prepared.evidence;
    let mut candidate_requests = Vec::new();
    for pending in &reference.pending_transfers {
        let mut request = r.clone();
        request.transfers = pending.clone();
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
    templates: &[(String, &crate::modified_duplex::TrackedDuplex)],
    task: &crate::topo_run_preparation::Task,
) -> Result<crate::topo_run_preparation::Reply, Error> {
    use crate::topo_run_preparation::{Reply, Task};
    source_preflight_reconstructed(inputs, task.preparation(), templates)?;
    Ok(match task {
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
        )?),
        Task::Project { request, selection } => {
            Reply::Projected(product_preflighted(inputs, request, selection)?)
        }
        Task::Check {
            request,
            selection,
            target,
            policy,
        } => {
            let TargetAssessment {
                product,
                policy,
                target_revision,
                outcome,
                ..
            } = check_fresh_product(
                product_preflighted(inputs, request, selection)?,
                target,
                policy,
            )?;
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
        Task::InspectActions { request } => Reply::ActionsInspected(
            crate::topo_actions::inspect_actions_reconstructed(inputs, request, templates)?,
        ),
        Task::ProcessActions { request } => Reply::ActionsProcessed(
            crate::topo_actions::inspect_reconstructed(inputs, request, templates)?,
        ),
    })
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
        return Err(bad("saved template role must occur exactly once"));
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
        return Err(bad("saved template must feed explicit exact PCR"));
    }
    Ok(())
}
