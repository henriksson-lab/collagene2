//! Expected junction sequence/frame inspection from fresh typed-plan reconstruction.
//! Coordinate correspondence and translation never establish fusion function,
//! experimental ligation, or reading-frame intent from feature labels.
use crate::{
    chemistry::{Duplex, DuplexSummary, End, Nick},
    cloning::assembly::recipe as assembly_recipe,
    planning::plan::{DesignPlan, Operation, OrientedSource, Source},
    sequence::translation::{self, InitiationPolicy, TranslationOptions, TranslationReport},
    verification::{self, evidence::VerificationReport},
};
use collagene_model::{Error, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const ENGINE_VERSION: &str = "assembly-junctions-1";
const MAX_JUNCTIONS: usize = 128;
const MAX_CDS: usize = 128;
const MAX_CDS_LENGTH: usize = 30_000;
const MAX_CDS_BASES: usize = 1_000_000;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)] pub struct Frame {pub reverse:bool,pub frame:u8}}
dto! {pub struct Options {pub flank_bases:usize,pub genetic_code:u8,pub frames:Vec<Frame>,pub max_junctions:usize,pub include_cds:bool}}
impl Default for Options {
    fn default() -> Self {
        Self {
            flank_bases: 30,
            genetic_code: 1,
            frames: [false, true]
                .into_iter()
                .flat_map(|reverse| (0..3).map(move |frame| Frame { reverse, frame }))
                .collect(),
            max_junctions: 128,
            include_cds: true,
        }
    }
}
impl Options {
    pub fn validate(&self) -> Result<(), Error> {
        if self.flank_bases == 0
            || self.flank_bases > 150
            || !matches!(self.genetic_code, 1 | 2 | 11)
            || self.frames.is_empty()
            || self.frames.len() > 6
            || self.max_junctions == 0
            || self.max_junctions > MAX_JUNCTIONS
        {
            return Err(bad(
                "invalid finite flank/frame/genetic-code/junction options",
            ));
        }
        let mut seen = BTreeSet::new();
        for f in &self.frames {
            if f.frame > 2 || !seen.insert((f.reverse, f.frame)) {
                return Err(bad(
                    "frames must be unique explicit strand/frame choices0..2",
                ));
            }
        }
        Ok(())
    }
}
dto! {#[derive(PartialEq,Eq)] pub struct Range {pub start:usize,pub end:usize}}
dto! {pub struct GoldenGatePart {pub path_index:usize,pub pool_fragment:usize,pub template_index:usize,pub template:Source,pub template_revision:String,pub template_length:usize,pub fragment_index:usize,pub top_start:usize,pub bottom_start:i64,pub reverse:bool,pub top_length:usize,pub bottom_length:usize,pub output_start:usize,pub output_end:usize,pub left_end:Option<End>,pub right_end:Option<End>}}
dto! {pub enum Method {Ligate,Overlap,ProfiledOverlap,GoldenGate{candidate_index:usize,path:Vec<GoldenGatePart>,enumeration_complete:bool,states_examined:usize,truncation_reasons:Vec<String>},HomologyAssembly,HomologyAnnealing}}
dto! {pub struct Window {pub parts:Vec<Range>,pub dna:String,pub boundary_offset:usize,pub left_bases:usize,pub right_bases:usize,pub clipped:bool}}
dto! {pub struct Codon {pub dna:String,pub amino_acid:char,/// Offset in the selected frame's oriented WINDOW, not the full product.
pub window_start:usize}}
dto! {pub struct FrameEvidence {pub spec:Frame,pub translation:TranslationReport,pub boundary_offset:usize,pub crossing_codon:Option<Codon>}}
dto! {pub struct CdsWindow {pub dna:String,pub protein:String,pub first_codon:usize,pub codon_count:usize,pub junction_coding_offset:usize,pub codon_start:u8,pub reverse:bool,pub genetic_code:u8,pub code_name:String,pub complete_under_model:bool}}
dto! {pub struct CdsEvidence {pub feature_index:usize,pub label:String,pub needs_review:bool,pub spans_boundary:bool,pub translation:Option<CdsWindow>,pub warning:Option<String>}}
dto! {pub struct JunctionPhysical {pub upstream_end:Option<End>,pub downstream_end:Option<End>,/// Full freshly reconstructed product nick state, not inferred bond-by-bond proof.
pub product_nicks:Vec<Nick>}}
dto! {pub struct Junction {pub index:usize,pub upstream:usize,pub downstream:usize,pub closing:bool,pub boundary:usize,pub homology:Vec<Range>,pub window:Window,pub frames:Vec<FrameEvidence>,pub cds:Vec<CdsEvidence>,pub cds_omitted:usize,pub physical:Option<JunctionPhysical>}}
dto! {pub struct Report {pub engine_version:String,pub context_revision:String,pub plan_revision:String,pub target_revision:String,pub options:Options,pub verification:VerificationReport,pub product_revision:Option<String>,pub product_length:Option<usize>,pub product_topology:Option<Topology>,pub method:Option<Method>,pub junctions:Vec<Junction>,pub complete:bool,pub omitted_junctions:usize,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("assembly junction inspection: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
/// Borrow originals before host clones. Every graph/source snapshot stays bounded
/// before canonical hashing; supported scientific engines retain their own limits.
pub fn source_preflight(plan: &DesignPlan, target: &Molecule) -> Result<(), Error> {
    if plan.steps.is_empty()
        || plan.steps.len() > 1000
        || plan.inputs.len() > 128
        || plan.name.len() > 1024
    {
        return Err(Error::TooLarge);
    }
    let mut models = plan.inputs.iter().map(|i| &i.molecule).collect::<Vec<_>>();
    models.push(target);
    collagene_model::bounds::model_preflight(&models)?;
    let mut bases = 0usize;
    for m in models {
        bases = bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= 8 << 20)
            .ok_or(Error::TooLarge)?;
        m.validate()?;
    }
    if plan
        .inputs
        .iter()
        .any(|i| i.id.is_empty() || i.id.len() > 128 || !hex(&i.revision))
        || !hex(&plan.target_revision)
    {
        return Err(bad("invalid original/target revision shape"));
    }
    verification::sources::validate_reference_steps(plan)?;
    assembly_recipe::hash(&(plan, target))?;
    Ok(())
}
pub fn inspection_revision(
    plan: &DesignPlan,
    target: &Molecule,
    options: &Options,
) -> Result<String, Error> {
    options.validate()?;
    source_preflight(plan, target)?;
    assembly_recipe::hash(&(ENGINE_VERSION, plan, target, options))
}
fn resolve<'a>(
    plan: &'a DesignPlan,
    products: &'a [Molecule],
    s: &Source,
) -> Result<&'a Molecule, Error> {
    match s {
        Source::Input(id) => plan
            .inputs
            .iter()
            .find(|i| &i.id == id)
            .map(|i| &i.molecule),
        Source::Step(i) => products.get(*i),
    }
    .ok_or_else(|| bad("prepared source is outside the actually reconstructed prefix"))
}
fn actual_physical(
    plan: &DesignPlan,
    products: &[Molecule],
    trace: &[Option<Duplex>],
    fragment: &OrientedSource,
) -> Result<DuplexSummary, Error> {
    match &fragment.source {
        Source::Step(i) => {
            let d = trace
                .get(*i)
                .and_then(Option::as_ref)
                .ok_or_else(|| bad("prepared source has no fresh physical trace"))?;
            if fragment.reverse {
                d.reverse()?.summary()
            } else {
                d.summary()
            }
        }
        Source::Input(_) => {
            let d = Duplex::from_molecule(resolve(plan, products, &fragment.source)?)?;
            if fragment.reverse {
                d.reverse()?.summary()
            } else {
                d.summary()
            }
        }
    }
}
struct Geometry {
    index: usize,
    upstream: usize,
    downstream: usize,
    boundary: usize,
    closing: bool,
    homology: Vec<Range>,
    physical: Option<JunctionPhysical>,
}
fn geometry(
    plan: &DesignPlan,
    verification: &VerificationReport,
    trace: &[Option<Duplex>],
    product: &Molecule,
    golden_gate: Option<&crate::cloning::golden_gate::GoldenGateReport>,
) -> Result<(Method, Vec<Geometry>), Error> {
    let n = product.sequence.len();
    let last = plan.steps.last().ok_or_else(|| bad("empty plan"))?;
    match last {
        Operation::GoldenGate {
            templates,
            options,
            candidate_index,
        }
        | Operation::SelectedGoldenGate {
            templates,
            options,
            candidate_index,
            ..
        } => {
            let pool = golden_gate
                .ok_or_else(|| bad("final Golden Gate has no fresh captured enumeration"))?;
            if pool.options != *options || pool.model != crate::cloning::golden_gate::MODEL {
                return Err(bad(
                    "fresh Golden Gate scope differs from the final operation",
                ));
            }
            let selected = pool
                .candidates
                .get(*candidate_index)
                .ok_or_else(|| bad("selected fresh Golden Gate path unavailable"))?;
            let final_state = trace
                .last()
                .and_then(Option::as_ref)
                .ok_or_else(|| bad("final Golden Gate has no actual physical state"))?;
            if selected.duplex != *final_state
                || selected.duplex.top != product.sequence
                || selected.duplex.topology != product.topology
                || selected.path.is_empty()
                || selected.path.len() > 32
                || final_state
                    .nicks
                    .len()
                    .checked_mul(selected.path.len())
                    .is_none_or(|n| n > 100_000)
            {
                return Err(bad(
                    "selected path, final duplex or physical output budget differs from actual product",
                ));
            }
            let mut path = Vec::with_capacity(selected.path.len());
            let mut offset = 0usize;
            for (index, choice) in selected.path.iter().enumerate() {
                let fragment = pool
                    .fragments
                    .get(choice.fragment)
                    .ok_or_else(|| bad("selected path fragment unavailable"))?;
                let source = templates
                    .get(fragment.source_index)
                    .ok_or_else(|| bad("selected path template unavailable"))?;
                let original = resolve(plan, &verification.products, source)?;
                if original.revision() != fragment.source_revision {
                    return Err(bad("actual digest template revision changed"));
                }
                let oriented = if choice.reverse {
                    fragment.duplex.reverse()?
                } else {
                    fragment.duplex.clone()
                };
                let end = offset
                    .checked_add(oriented.top.len())
                    .filter(|e| *e <= n)
                    .ok_or(Error::TooLarge)?;
                if product.sequence.as_str().get(offset..end) != Some(oriented.top.as_str()) {
                    return Err(bad(
                        "selected oriented path bases differ from actual final product",
                    ));
                }
                let (left, right) = oriented.ends()?;
                path.push(GoldenGatePart {
                    path_index: index,
                    pool_fragment: choice.fragment,
                    template_index: fragment.source_index,
                    template: source.clone(),
                    template_revision: fragment.source_revision.clone(),
                    template_length: original.sequence.len(),
                    fragment_index: fragment.fragment_index,
                    top_start: fragment.top_start,
                    bottom_start: fragment.bottom_start,
                    reverse: choice.reverse,
                    top_length: oriented.top.len(),
                    bottom_length: oriented.bottom.len(),
                    output_start: offset,
                    output_end: end,
                    left_end: Some(left),
                    right_end: Some(right),
                });
                offset = end;
            }
            if offset != n {
                return Err(bad("selected path does not cover the whole actual product"));
            }
            let mut out = Vec::new();
            for (i, part) in path.iter().enumerate() {
                if i + 1 < path.len() || options.circular {
                    let next = (i + 1) % path.len();
                    out.push(Geometry {
                        index: i,
                        upstream: i,
                        downstream: next,
                        boundary: if next == 0 { 0 } else { part.output_end },
                        closing: next == 0,
                        homology: vec![],
                        physical: Some(JunctionPhysical {
                            upstream_end: part.right_end.clone(),
                            downstream_end: path[next].left_end.clone(),
                            product_nicks: final_state.nicks.clone(),
                        }),
                    });
                }
            }
            Ok((
                Method::GoldenGate {
                    candidate_index: *candidate_index,
                    path,
                    enumeration_complete: pool.enumeration_complete,
                    states_examined: pool.states_examined,
                    truncation_reasons: pool.truncation_reasons.clone(),
                },
                out,
            ))
        }
        Operation::Ligate {
            fragments,
            circular,
            ..
        } => {
            if fragments.is_empty() || fragments.len() > 1000 {
                return Err(Error::TooLarge);
            }
            let mut source_nicks = 0usize;
            for fragment in fragments {
                if let Source::Step(i) = &fragment.source {
                    if let Some(Some(d)) = trace.get(*i) {
                        source_nicks = source_nicks
                            .checked_add(d.nicks.len())
                            .filter(|n| *n <= 100_000)
                            .ok_or(Error::TooLarge)?;
                    }
                }
            }
            let final_nicks = trace
                .get(plan.steps.len() - 1)
                .and_then(Option::as_ref)
                .map(|d| d.nicks.len())
                .unwrap_or(0);
            if final_nicks
                .checked_mul(fragments.len())
                .is_none_or(|n| n > 100_000)
            {
                return Err(bad(
                    "junction physical nick evidence exceeds the aggregate100000-copy output budget",
                ));
            }
            let summaries = fragments
                .iter()
                .map(|f| actual_physical(plan, &verification.products, trace, f))
                .collect::<Result<Vec<_>, _>>()?;
            let total = summaries.iter().try_fold(0usize, |sum, d| {
                sum.checked_add(d.top_length).ok_or(Error::TooLarge)
            })?;
            if total != n {
                return Err(bad(
                    "fresh oriented duplex lengths disagree with final product",
                ));
            }
            let final_state = trace
                .get(plan.steps.len() - 1)
                .and_then(Option::as_ref)
                .ok_or_else(|| bad("final ligation lacks fresh physical state"))?;
            let mut out = Vec::new();
            let mut end = 0usize;
            for (i, s) in summaries.iter().enumerate() {
                end = end.checked_add(s.top_length).ok_or(Error::TooLarge)?;
                if i + 1 < summaries.len() || *circular {
                    let downstream = (i + 1) % summaries.len();
                    out.push(Geometry {
                        index: i,
                        upstream: i,
                        downstream,
                        boundary: if downstream == 0 { 0 } else { end },
                        closing: downstream == 0,
                        homology: vec![],
                        physical: Some(JunctionPhysical {
                            upstream_end: s.ends.as_ref().map(|(_, r)| r.clone()),
                            downstream_end: summaries[downstream]
                                .ends
                                .as_ref()
                                .map(|(l, _)| l.clone()),
                            product_nicks: final_state.nicks.clone(),
                        }),
                    });
                }
            }
            Ok((Method::Ligate, out))
        }
        Operation::Overlap {
            fragments,
            overlaps,
            circular,
        } => overlap_geometry(
            plan,
            verification,
            product,
            &fragments
                .iter()
                .cloned()
                .map(|source| OrientedSource {
                    source,
                    reverse: false,
                })
                .collect::<Vec<_>>(),
            overlaps,
            *circular,
            Method::Overlap,
        ),
        Operation::HomologyAnnealing {
            fragments, policy, ..
        } => {
            let step = plan
                .steps
                .len()
                .checked_sub(1)
                .ok_or_else(|| bad("empty annealing plan"))?;
            let native = verification
                .physical_processing_reports
                .iter()
                .find_map(|entry| {
                    if entry.step_index == step {
                        match &entry.evidence {
                            verification::evidence::PhysicalProcessingEvidence::HomologyAnnealing(report) => {
                                Some(report)
                            }
                            _ => None,
                        }
                    } else {
                        None
                    }
                })
                .ok_or_else(|| bad("annealing physical strand evidence is unavailable"))?;
            native.validate_documentary()?;
            if assembly_recipe::hash(&native.options)? != assembly_recipe::hash(policy)?
                || native.reference != product.sequence
                || policy.topology != product.topology
            {
                return Err(bad("annealing junction reference/policy differs"));
            }
            overlap_geometry(
                plan,
                verification,
                product,
                fragments,
                &policy.overlaps,
                policy.topology == Topology::Circular,
                Method::HomologyAnnealing,
            )
        }
        Operation::HomologyAssembly {
            fragments, policy, ..
        } => {
            let step = plan
                .steps
                .len()
                .checked_sub(1)
                .ok_or_else(|| bad("empty homology plan"))?;
            let native = verification
                .physical_processing_reports
                .iter()
                .find_map(|entry| {
                    if entry.step_index != step {
                        return None;
                    }
                    match &entry.evidence {
                        verification::evidence::PhysicalProcessingEvidence::HomologyAssembly(
                            report,
                        ) => Some(report),
                        _ => None,
                    }
                })
                .ok_or_else(|| bad("native homology strand evidence is unavailable"))?;
            native.validate_documentary()?;
            if assembly_recipe::hash(&native.options)? != assembly_recipe::hash(policy)?
                || native.source_geometry_revisions.len() != fragments.len()
                || native.product.as_ref().is_none_or(|duplex| {
                    duplex.top != product.sequence || duplex.topology != product.topology
                })
            {
                return Err(bad(
                    "homology junction evidence does not match the actual final product/policy",
                ));
            }
            let physical = native
                .product
                .as_ref()
                .ok_or_else(|| bad("homology product withheld"))?;
            let mut output = Vec::with_capacity(native.junctions.len());
            for (index, junction) in native.junctions.iter().enumerate() {
                let start = usize::try_from(junction.homology_start)
                    .map_err(|_| bad("homology lies outside the displayed top reference"))?;
                let end = usize::try_from(junction.homology_end)
                    .map_err(|_| bad("homology lies outside the displayed top reference"))?;
                if start >= end
                    || end > product.sequence.len()
                    || junction.left_source >= fragments.len()
                    || junction.right_source >= fragments.len()
                {
                    return Err(bad(
                        "native homology coordinates exceed this displayed top reference",
                    ));
                }
                let closing = policy.topology == Topology::Circular && junction.right_source == 0;
                output.push(Geometry {
                    index,
                    upstream: junction.left_source,
                    downstream: junction.right_source,
                    closing,
                    boundary: if closing { 0 } else { end },
                    homology: vec![Range { start, end }],
                    physical: Some(JunctionPhysical {
                        upstream_end: None,
                        downstream_end: None,
                        product_nicks: physical.nicks.clone(),
                    }),
                });
            }
            Ok((Method::HomologyAssembly, output))
        }
        Operation::ProfiledOverlap {
            fragments,
            overlaps,
            circular,
            ..
        } => overlap_geometry(
            plan,
            verification,
            product,
            fragments,
            overlaps,
            *circular,
            Method::ProfiledOverlap,
        ),
        _ => Err(bad(
            "final operation is not supported by this junction adapter; fresh verification remains available",
        )),
    }
}
fn overlap_geometry(
    plan: &DesignPlan,
    verification: &VerificationReport,
    product: &Molecule,
    fragments: &[OrientedSource],
    overlaps: &[usize],
    circular: bool,
    method: Method,
) -> Result<(Method, Vec<Geometry>), Error> {
    if fragments.is_empty()
        || fragments.len() > 1000
        || overlaps.len() != fragments.len() - 1 + usize::from(circular)
    {
        return Err(bad("invalid final overlap operation shape"));
    }
    let lengths = fragments
        .iter()
        .map(|f| resolve(plan, &verification.products, &f.source).map(|m| m.sequence.len()))
        .collect::<Result<Vec<_>, _>>()?;
    let n = product.sequence.len();
    if overlaps.iter().any(|k| *k == 0 || *k > n) {
        return Err(bad(
            "shared homology exceeds one final-product turn; repeated or target-guided placement is withheld",
        ));
    }
    let mut expected = lengths
        .iter()
        .try_fold(0usize, |a, b| a.checked_add(*b).ok_or(Error::TooLarge))?;
    for k in overlaps {
        expected = expected
            .checked_sub(*k)
            .ok_or_else(|| bad("overlap extent underflow"))?;
    }
    if expected != n {
        return Err(bad(
            "prepared lengths/overlaps disagree with reconstructed product",
        ));
    }
    let mut end = lengths[0];
    let mut out = Vec::new();
    for (i, &k) in overlaps.iter().enumerate() {
        let closing = i + 1 == fragments.len();
        let (start, stop, boundary) = if closing {
            (0, k, 0)
        } else {
            let start = end
                .checked_sub(k)
                .ok_or_else(|| bad("homology position underflow"))?;
            if end > n {
                return Err(bad(
                    "homology interval cannot be placed without repeating product calls",
                ));
            }
            (start, end, end)
        };
        out.push(Geometry {
            index: i,
            upstream: i,
            downstream: if closing { 0 } else { i + 1 },
            closing,
            boundary,
            homology: vec![Range { start, end: stop }],
            physical: None,
        });
        if !closing {
            end = end.checked_add(lengths[i + 1] - k).ok_or(Error::TooLarge)?;
        }
    }
    Ok((method, out))
}
fn window(product: &Molecule, boundary: usize, flank: usize) -> Result<Window, Error> {
    let n = product.sequence.len();
    if n == 0 || boundary > n {
        return Err(bad("invalid product boundary"));
    }
    let (left, right, start) = if product.topology == Topology::Circular {
        if boundary >= n {
            return Err(bad("circular boundary must be normalized explicitly"));
        }
        let left = flank.min(n / 2);
        let right = flank.min(n - left);
        (left, right, (boundary + n - left) % n)
    } else {
        let left = flank.min(boundary);
        let right = flank.min(n - boundary);
        (left, right, boundary - left)
    };
    let count = left + right;
    let mut parts = vec![];
    let end = start + count;
    if end <= n {
        if count > 0 {
            parts.push(Range { start, end });
        }
    } else {
        parts.push(Range { start, end: n });
        parts.push(Range {
            start: 0,
            end: end - n,
        });
    }
    let mut dna = String::with_capacity(count);
    for r in &parts {
        dna.push_str(&product.sequence.as_str()[r.start..r.end]);
    }
    Ok(Window {
        parts,
        dna,
        boundary_offset: left,
        left_bases: left,
        right_bases: right,
        clipped: left < flank || right < flank,
    })
}
fn frames(w: &Window, options: &Options) -> Result<Vec<FrameEvidence>, Error> {
    let sequence = Sequence::dna(&w.dna)?;
    options
        .frames
        .iter()
        .map(|spec| {
            let opts = TranslationOptions {
                genetic_code: options.genetic_code,
                frame: spec.frame,
                reverse: spec.reverse,
                initiation: InitiationPolicy::Ordinary,
            };
            let translation = translation::translate(&sequence, &opts)?;
            let boundary = if spec.reverse {
                w.dna.len() - w.boundary_offset
            } else {
                w.boundary_offset
            };
            let oriented = if spec.reverse {
                sequence.reverse_complement()?
            } else {
                sequence.clone()
            };
            let frame = spec.frame as usize;
            let crossing_codon = if boundary > frame {
                let start = frame + (boundary - frame - 1) / 3 * 3;
                let end = start + 3;
                if start < boundary && end > boundary && end <= oriented.len() {
                    Some(Codon {
                        dna: oriented.as_str()[start..end].into(),
                        amino_acid: translation.protein.as_bytes()[(start - frame) / 3] as char,
                        window_start: start,
                    })
                } else {
                    None
                }
            } else {
                None
            };
            Ok(FrameEvidence {
                spec: spec.clone(),
                translation,
                boundary_offset: boundary,
                crossing_codon,
            })
        })
        .collect()
}
fn exact_cds_parts(location: &Location) -> Result<(Vec<(usize, usize)>, bool, usize), Error> {
    let mut root = location;
    let mut reverse = false;
    while let Location::Complement(inner) = root {
        reverse = !reverse;
        root = inner;
    }
    let mut stack = vec![root];
    let mut parts = vec![];
    let mut length = 0usize;
    while let Some(p) = stack.pop() {
        match p {
            Location::Range {
                start,
                end,
                before: false,
                after: false,
            } => {
                length = length
                    .checked_add(end - start)
                    .filter(|n| *n <= MAX_CDS_LENGTH)
                    .ok_or_else(|| bad("CDS exceeds the per-feature30000nt translation budget"))?;
                parts.push((*start, *end));
            }
            Location::Join(children) => {
                for p in children.iter().rev() {
                    stack.push(p)
                }
            }
            _ => {
                return Err(bad(
                    "CDS fuzzy/order/mixed-strand geometry cannot determine a local junction frame",
                ));
            }
        }
    }
    let mut sorted = parts.clone();
    sorted.sort_unstable();
    if sorted.windows(2).any(|w| w[0].1 > w[1].0) {
        return Err(bad(
            "CDS repeats reference calls; the junction has no unique coding offset",
        ));
    }
    if reverse {
        parts.reverse();
    }
    Ok((parts, reverse, length))
}
fn coding_offset(parts: &[(usize, usize)], reverse: bool, position: usize) -> Option<usize> {
    let mut offset = 0usize;
    for &(a, b) in parts {
        if a <= position && position < b {
            return Some(
                offset
                    + if reverse {
                        b - 1 - position
                    } else {
                        position - a
                    },
            );
        }
        offset += b - a;
    }
    None
}
fn local_cds(
    product: &Molecule,
    feature: &collagene_model::Feature,
    boundary: usize,
    revision: &str,
    options: &Options,
    used: &mut usize,
) -> Result<CdsWindow, Error> {
    let (parts, reverse, length) = exact_cds_parts(&feature.location)?;
    *used = used
        .checked_add(length)
        .filter(|n| *n <= MAX_CDS_BASES)
        .ok_or_else(|| bad("aggregate1Mnt CDS translation budget exhausted"))?;
    let n = product.sequence.len();
    if product.topology == Topology::Linear && (boundary == 0 || boundary == n) {
        return Err(bad(
            "linear terminal boundary has no two adjacent coding calls",
        ));
    }
    let left = (boundary + n - 1) % n;
    let right = boundary % n;
    let a = coding_offset(&parts, reverse, left)
        .ok_or_else(|| bad("CDS does not cover the upstream boundary call"))?;
    let b = coding_offset(&parts, reverse, right)
        .ok_or_else(|| bad("CDS does not cover the downstream boundary call"))?;
    let offset = if reverse {
        if a != b + 1 {
            return Err(bad(
                "product boundary is not adjacent in this CDS biological coordinate order",
            ));
        }
        a
    } else {
        if b != a + 1 {
            return Err(bad(
                "product boundary is not adjacent in this CDS biological coordinate order",
            ));
        }
        b
    };
    let full = translation::translate_cds_validated(product, feature, revision)?;
    let skip = (full.codon_start - 1) as usize;
    if offset < skip {
        return Err(bad(
            "junction falls within the skipped incomplete5prime CDS codon",
        ));
    }
    let coding = offset - skip;
    let center = coding / 3;
    let radius = (options.flank_bases / 3).max(1);
    let start = center
        .saturating_sub(radius)
        .min(full.translation.protein.len());
    let end = (center + radius + 1).min(full.translation.protein.len());
    let dna = feature.location.extract(&product.sequence)?;
    let a = skip + start * 3;
    let b = skip + end * 3;
    if b > dna.len() {
        return Err(bad("CDS codon window exceeds actual extracted calls"));
    }
    Ok(CdsWindow {
        dna: dna[a..b].into(),
        protein: full.translation.protein[start..end].into(),
        first_codon: start,
        codon_count: end - start,
        junction_coding_offset: coding,
        codon_start: full.codon_start,
        reverse,
        genetic_code: full.translation.options.genetic_code,
        code_name: full.translation.code_name,
        complete_under_model: full.complete_under_model,
    })
}
fn label(text: &str) -> String {
    if text.len() <= 512 {
        return text.into();
    }
    let mut end = 512;
    while !text.is_char_boundary(end) {
        end -= 1;
    }
    format!("{}…", &text[..end])
}
fn intersects(spans: &[(usize, usize, bool)], window: &Window) -> bool {
    spans
        .iter()
        .any(|&(a, b, _)| window.parts.iter().any(|r| a < r.end && r.start < b))
}
fn spans_boundary(
    spans: &[(usize, usize, bool)],
    boundary: usize,
    length: usize,
    topology: Topology,
) -> bool {
    if topology == Topology::Linear && (boundary == 0 || boundary == length) {
        return false;
    }
    let a = (boundary + length - 1) % length;
    let b = boundary % length;
    [a, b]
        .iter()
        .all(|p| spans.iter().any(|&(a, b, _)| a <= *p && *p < b))
}
pub fn inspect_plan(
    plan: &DesignPlan,
    target: &Molecule,
    options: &Options,
) -> Result<Report, Error> {
    let context_revision = inspection_revision(plan, target, options)?;
    let (verification, trace, golden_gate) =
        verification::engine::verify_with_golden_gate_trace(plan, target);
    let mut report=Report{engine_version:ENGINE_VERSION.into(),context_revision,plan_revision:assembly_recipe::hash(plan)?,target_revision:assembly_recipe::hash(target)?,options:options.clone(),verification,product_revision:None,product_length:None,product_topology:None,method:None,junctions:vec![],complete:false,omitted_junctions:0,warnings:vec!["Junction DNA is the expected reconstructed product, not measured material. Window-relative translations do not establish a CDS reading frame, fusion function or experimental success.".into(),"Frames0/1/2 skip calls at the beginning of each explicitly oriented WINDOW; they are not product-origin frames. CDS windows, when available, instead use explicit feature-directed coding geometry and the existing genetic-code engine.".into()]};
    if report.verification.products.len() != plan.steps.len() {
        report.warnings.push("Final operation did not reconstruct a product. Earlier preparation products are retained in verification evidence but are never presented as final junction DNA.".into());
        report.validate()?;
        return Ok(report);
    }
    let product = report
        .verification
        .products
        .last()
        .ok_or_else(|| bad("missing final reconstructed product"))?;
    collagene_model::bounds::model_preflight(&[product])?;
    let product_revision = assembly_recipe::hash(product)?;
    report.product_revision = Some(product_revision.clone());
    report.product_length = Some(product.sequence.len());
    report.product_topology = Some(product.topology);
    let (method, geometry) = match geometry(
        plan,
        &report.verification,
        &trace,
        product,
        golden_gate.as_ref(),
    ) {
        Ok(v) => v,
        Err(error) => {
            report.warnings.push(error.to_string());
            report.validate()?;
            return Ok(report);
        }
    };
    report.method = Some(method);
    report.complete = true;
    report.omitted_junctions = geometry.len().saturating_sub(options.max_junctions);
    if report.omitted_junctions > 0 {
        report.complete = false;
        report.warnings.push(format!(
            "{} junctions were withheld by the requested finite output limit.",
            report.omitted_junctions
        ));
    }
    let cds = if options.include_cds {
        product
            .features
            .iter()
            .enumerate()
            .filter(|(_, f)| f.kind.eq_ignore_ascii_case("CDS"))
            .map(|(i, f)| (i, f, f.location.spans()))
            .collect::<Vec<_>>()
    } else {
        vec![]
    };
    let unplaced = cds.iter().filter(|(_, _, spans)| spans.is_empty()).count();
    if unplaced > 0 {
        report.complete = false;
        report.warnings.push(format!("{unplaced} CDS annotations have no placeable source ranges; no junction relationship or coding frame is inferred."));
    }
    let mut cds_count = 0usize;
    let mut cds_bases = 0usize;
    for g in geometry.into_iter().take(options.max_junctions) {
        let w = window(product, g.boundary, options.flank_bases)?;
        let frame_evidence = frames(&w, options)?;
        let mut evidence = vec![];
        let mut omitted = 0usize;
        for (index, feature, spans) in &cds {
            if !intersects(spans, &w) {
                continue;
            }
            if cds_count >= MAX_CDS {
                omitted += 1;
                continue;
            }
            cds_count += 1;
            let spans_boundary =
                spans_boundary(spans, g.boundary, product.sequence.len(), product.topology);
            let (translation, warning) = if spans_boundary {
                match local_cds(
                    product,
                    feature,
                    g.boundary,
                    &product_revision,
                    options,
                    &mut cds_bases,
                ) {
                    Ok(t) => {
                        let warning = if !t.complete_under_model {
                            report.complete = false;
                            Some("Existing CDS engine reports partial/ambiguous/incomplete coding evidence; the local window does not establish fusion function.".into())
                        } else {
                            None
                        };
                        (Some(t), warning)
                    }
                    Err(error) => {
                        report.complete = false;
                        (None, Some(error.to_string()))
                    }
                }
            } else {
                (None,Some("CDS intersects the flank window but does not span both adjacent product calls at this boundary; no cross-junction coding frame is asserted.".into()))
            };
            let mut warning = warning;
            if feature.label().len() > 512 {
                warning = Some(format!(
                    "{} Feature label abbreviated; actual feature index identifies full native evidence.",
                    warning.unwrap_or_default()
                ));
            }
            evidence.push(CdsEvidence {
                feature_index: *index,
                label: label(feature.label()),
                needs_review: feature.needs_review,
                spans_boundary,
                translation,
                warning,
            });
        }
        if omitted > 0 {
            report.complete = false;
        }
        report.junctions.push(Junction {
            index: g.index,
            upstream: g.upstream,
            downstream: g.downstream,
            closing: g.closing,
            boundary: g.boundary,
            homology: g.homology,
            window: w,
            frames: frame_evidence,
            cds: evidence,
            cds_omitted: omitted,
            physical: g.physical,
        });
    }
    if report.junctions.iter().any(|j| j.window.clipped) {
        report.warnings.push("Some requested flanks exceed available unique product calls. Windows are clipped to one turn and never repeat short circular DNA to fill a display.".into());
    }
    if report.junctions.iter().any(|j| {
        j.homology
            .iter()
            .any(|r| r.end - r.start > j.window.dna.len())
    }) {
        report.warnings.push("Some complete homology intervals exceed the bounded display window. Exact product ranges are retained; the window is not represented as full homology sequence.".into());
    }
    report.validate()?;
    Ok(report)
}
impl Report {
    /// Readonly shape/binding validation; no scientific reconstruction.
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.context_revision)
            || !hex(&self.plan_revision)
            || !hex(&self.target_revision)
            || self.junctions.len() > self.options.max_junctions
            || self.verification.products.len() > 1000
            || self.verification.step_reports.len() > 1000
            || self.warnings.len() > 64
            || self.warnings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad("invalid bounded report shape"));
        }
        if self.verification.plan_revision != self.plan_revision
            || self.verification.target_revision != self.target_revision
        {
            return Err(bad("verification context binding changed"));
        }
        for m in &self.verification.products {
            collagene_model::bounds::model_preflight(&[m])?;
        }
        for report in &self.verification.step_reports {
            if report.step_index >= 1000 {
                return Err(Error::TooLarge);
            }
            if let Some(product) = &report.topo {
                crate::cloning::topo::worker::validate_projection(
                    product,
                    &product.batch_revision,
                )?;
            }
        }
        let presence = self.product_revision.is_some();
        if presence != self.product_length.is_some()
            || presence != self.product_topology.is_some()
            || self.product_revision.as_ref().is_some_and(|s| !hex(s))
            || !presence && (!self.junctions.is_empty() || self.method.is_some())
            || self.complete && (!presence || self.method.is_none() || self.omitted_junctions > 0)
        {
            return Err(bad("inconsistent product/method availability"));
        }
        if let Some(revision) = &self.product_revision {
            let product = self
                .verification
                .products
                .last()
                .ok_or_else(|| bad("product receipt has no actual reconstructed DNA"))?;
            if assembly_recipe::hash(product)? != *revision
                || Some(product.sequence.len()) != self.product_length
                || Some(product.topology) != self.product_topology
            {
                return Err(bad("final product receipt changed"));
            }
        }
        if let Some(Method::GoldenGate {
            candidate_index: _,
            path,
            enumeration_complete: _,
            states_examined,
            truncation_reasons,
        }) = &self.method
        {
            if path.is_empty()
                || path.len() > 32
                || *states_examined > 1_000_001
                || truncation_reasons.len() > 64
                || truncation_reasons.iter().any(|s| s.len() > 8192)
            {
                return Err(bad("invalid finite Golden Gate path scope"));
            }
            let product = self
                .verification
                .products
                .last()
                .ok_or_else(|| bad("missing actual Golden Gate product"))?;
            let mut end = 0usize;
            let mut used = BTreeSet::new();
            for (i, part) in path.iter().enumerate() {
                if part.path_index != i
                    || part.pool_fragment >= 128
                    || !used.insert(part.pool_fragment)
                    || part.template_index >= 32
                    || !hex(&part.template_revision)
                    || part.template_length == 0
                    || part.template_length > 8 << 20
                    || part.fragment_index >= 128
                    || part.top_start >= part.template_length
                    || part.bottom_start.unsigned_abs() > 32 << 20
                    || part.output_start != end
                    || part.output_end <= end
                    || part.output_end > product.sequence.len()
                    || part.top_length != part.output_end - part.output_start
                    || part.bottom_length > 16 << 20
                    || part.left_end.is_none()
                    || part.right_end.is_none()
                {
                    return Err(bad("Golden Gate path mapping changed"));
                }
                match &part.template {
                    Source::Input(id)
                        if !id.is_empty()
                            && id.len() <= 128
                            && !id.chars().any(char::is_control) => {}
                    Source::Step(step)
                        if *step < self.verification.products.len().saturating_sub(1) => {}
                    _ => {
                        return Err(bad(
                            "Golden Gate path refers beyond actual preparation prefix",
                        ));
                    }
                }
                if let Source::Step(step) = &part.template {
                    let source = &self.verification.products[*step];
                    if source.revision() != part.template_revision
                        || source.sequence.len() != part.template_length
                    {
                        return Err(bad("Golden Gate prepared template binding changed"));
                    }
                }
                end = part.output_end;
            }
            let total = path.len() - 1 + usize::from(product.topology == Topology::Circular);
            if end != product.sequence.len()
                || self.junctions.len() + self.omitted_junctions != total
            {
                return Err(bad("Golden Gate path/junction coverage changed"));
            }
            for (i, j) in self.junctions.iter().enumerate() {
                let next = (i + 1) % path.len();
                if j.upstream != i
                    || j.downstream != next
                    || j.closing != (next == 0)
                    || j.boundary != if next == 0 { 0 } else { path[i].output_end }
                    || !j.homology.is_empty()
                {
                    return Err(bad(
                        "Golden Gate junction differs from selected path positions",
                    ));
                }
                let physical = j
                    .physical
                    .as_ref()
                    .ok_or_else(|| bad("Golden Gate junction lacks actual ends"))?;
                if physical.upstream_end != path[i].right_end
                    || physical.downstream_end != path[next].left_end
                {
                    return Err(bad("Golden Gate junction ends differ from selected parts"));
                }
            }
        }
        let mut cds_count = 0usize;
        for (i, j) in self.junctions.iter().enumerate() {
            let n = self
                .product_length
                .ok_or_else(|| bad("junction has no product"))?;
            if j.index != i
                || j.upstream >= 1000
                || j.downstream >= 1000
                || j.boundary > n
                || self.product_topology == Some(Topology::Circular) && j.boundary >= n
                || j.homology.len() > 2
                || j.window.parts.len() > 2
                || j.window.dna.is_empty()
                || j.window.dna.len() > 2 * self.options.flank_bases
                || j.window.boundary_offset != j.window.left_bases
                || j.window.left_bases.checked_add(j.window.right_bases) != Some(j.window.dna.len())
                || j.window.left_bases > self.options.flank_bases
                || j.window.right_bases > self.options.flank_bases
                || j.frames.len() != self.options.frames.len()
                || j.cds.len() > MAX_CDS
            {
                return Err(bad("invalid junction/window shape"));
            }
            let product = self
                .verification
                .products
                .last()
                .ok_or_else(|| bad("missing window source"))?;
            let mut length = 0usize;
            let mut actual = String::new();
            for r in &j.window.parts {
                if r.start >= r.end || r.end > n {
                    return Err(bad("invalid actual window interval"));
                }
                length = length
                    .checked_add(r.end - r.start)
                    .filter(|v| *v <= 2 * self.options.flank_bases && *v <= n)
                    .ok_or(Error::TooLarge)?;
                actual.push_str(&product.sequence.as_str()[r.start..r.end]);
            }
            if length != j.window.dna.len() || actual != j.window.dna {
                return Err(bad("window DNA differs from actual final product"));
            }
            for r in &j.homology {
                if r.start >= r.end || r.end > n {
                    return Err(bad("homology exceeds one actual product turn"));
                }
            }
            for (frame, spec) in j.frames.iter().zip(&self.options.frames) {
                if frame.spec != *spec
                    || frame.translation.protein.len() > 100
                    || frame.translation.nucleotide_length != j.window.dna.len()
                    || frame.translation.options.reverse != spec.reverse
                    || frame.translation.options.frame != spec.frame
                    || frame.translation.options.genetic_code != self.options.genetic_code
                    || frame.boundary_offset
                        != if spec.reverse {
                            j.window.dna.len() - j.window.boundary_offset
                        } else {
                            j.window.boundary_offset
                        }
                {
                    return Err(bad("frame window binding changed"));
                }
                if frame.crossing_codon.as_ref().is_some_and(|c| {
                    c.dna.len() != 3
                        || c.window_start
                            .checked_add(3)
                            .is_none_or(|e| e > j.window.dna.len())
                        || !c.amino_acid.is_ascii()
                }) {
                    return Err(bad("invalid codon window"));
                }
            }
            cds_count += j.cds.len();
            if cds_count > MAX_CDS {
                return Err(Error::TooLarge);
            }
            for cds in &j.cds {
                if cds.label.len() > 515
                    || cds.warning.as_ref().is_some_and(|s| s.len() > 8192)
                    || cds.feature_index >= product.features.len()
                {
                    return Err(Error::TooLarge);
                }
                if let Some(t) = &cds.translation {
                    if t.dna.len() > 303
                        || t.protein.len() > 101
                        || Some(t.dna.len()) != t.codon_count.checked_mul(3)
                        || t.codon_count > 101
                        || t.first_codon > 10_000
                        || t.junction_coding_offset > MAX_CDS_LENGTH
                        || t.protein.len() != t.codon_count
                        || !matches!(t.codon_start, 1..=3)
                        || !matches!(t.genetic_code, 1 | 2 | 11)
                        || t.code_name.len() > 256
                        || !cds.spans_boundary
                    {
                        return Err(bad("invalid bounded CDS window"));
                    }
                }
            }
            if j.physical.as_ref().is_some_and(|p| {
                p.product_nicks.len() > 10_000
                    || p.upstream_end
                        .as_ref()
                        .is_some_and(|e| e.overhang.len() > 8 << 20)
                    || p.downstream_end
                        .as_ref()
                        .is_some_and(|e| e.overhang.len() > 8 << 20)
            }) {
                return Err(Error::TooLarge);
            }
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad(
                "inspection belongs to another plan/target/options context",
            ));
        }
        Ok(())
    }
}
