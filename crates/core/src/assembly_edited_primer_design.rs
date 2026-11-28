//! Whole overlap preparation proposals from original source contributions,
//! junction declarations and literal internal edits. No arbitrary target input.
use crate::{
    assembly_edit_intent::{self, Edit},
    assembly_junction_edit,
    assembly_recipe::{self, FragmentSpec, Preparation, RunPreparation},
    overlap_extension, overlap_extension_pcr,
    overlap_extension_primer_design::{self as extension, PairSearch},
    overlap_primer_design::{self, Fragment},
    overlap_profiles::ProfileOptions,
    primers::{self, DesignOptions, DesignReport},
    run_material::MaterialOrigin,
    verification::{self, Operation, PlanInput},
    workspace::Workspace,
};
use collagene_model::{Error, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ENGINE_VERSION: &str = "assembly-edited-primer-design-2";
const MAX_WORK: usize = 100_000_000;
const MAX_BYTES: usize = 32 << 20;
const MAX_ROWS: usize = 32;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request{pub intent:assembly_edit_intent::Request,pub overlap:usize,pub options:extension::Options,pub profile_options:ProfileOptions,pub finishing:assembly_recipe::OverlapFinishing}}
dto! {#[derive(PartialEq,Eq)] pub struct TailDeclaration{pub role:String,pub input_id:String,pub forward_tail:String,pub reverse_tail:String,pub outgoing_insert:String,pub homology:String,pub downstream_role:Option<String>}}
// Documentary source-view evidence for inspecting a proposal. This cannot be
// constructed into FreshMaterial or used as a physical template by the verifier.
dto! {pub struct SourceFrame{pub input_id:String,pub template_strand:crate::chemistry::Strand,pub state_revision:String,pub sequence:Sequence,pub topology:Topology}}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RowDesign {
    pub role: String,
    pub input_id: String,
    pub input_revision: String,
    pub fragment: Fragment,
    pub template_strand: crate::chemistry::Strand,
    pub edits: Vec<Edit>,
    pub designs: Vec<DesignReport>,
    pub search: Vec<PairSearch>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine_version: String,
    pub context_revision: String,
    pub source_revisions: Vec<String>,
    pub source_frames: Vec<SourceFrame>,
    pub request: Request,
    pub intent: assembly_edit_intent::Report,
    pub recipe: assembly_recipe::Recipe,
    pub tails: Vec<TailDeclaration>,
    pub rows: Vec<RowDesign>,
    pub preview: assembly_recipe::Preview,
    pub warnings: Vec<String>,
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("edited assembly primer planning: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn concrete(s: &str) -> bool {
    s.bytes().all(|b| b"ACGT".contains(&b))
}
struct Count(usize);
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| std::io::Error::other("bounded whole primer proposal"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn size(value: &impl Serialize) -> Result<(), Error> {
    serde_json::to_writer(Count(0), value).map_err(|_| Error::TooLarge)
}
fn shape(r: &Request) -> Result<(), Error> {
    let base = &r.intent.base.base;
    if base.fragments.is_empty()
        || base.fragments.len() > MAX_ROWS
        || base.name.len() > 1024
        || r.overlap == 0
        || r.overlap > 80
        || r.options.primer_options.max_pairs == 0
        || r.options.primer_options.max_pairs > 20
        || !r.options.outer_forward_tail.is_empty()
        || !r.options.outer_reverse_tail.is_empty()
        || !r.options.primer_options.forward_tail.is_empty()
        || !r.options.primer_options.reverse_tail.is_empty()
    {
        return Err(bad(
            "declare1..32 source roles, overlap1..80, bounded native pairs and empty global tails; derived assembly tails cannot silently replace declared extra DNA",
        ));
    }
    if base.fragments.iter().any(|f| {
        f.role.is_empty() || f.role.len() > 128 || f.input_id.is_empty() || f.input_id.len() > 128
    }) {
        return Err(bad(
            "bound all source and role identities before copying declarations",
        ));
    }
    r.profile_options.validate()?;
    r.finishing.validate_geometry(
        base.fragments.len(),
        &vec![r.overlap; base.fragments.len() - 1 + usize::from(base.circular)],
        base.circular,
        &r.profile_options,
    )?;
    overlap_primer_design::options(&overlap_primer_design::Request {
        name: base.name.clone(),
        fragments: base
            .fragments
            .iter()
            .map(frame_fragment)
            .collect::<Result<Vec<_>, _>>()?,
        overlap: r.overlap,
        circular: base.circular,
        primer_options: r.options.primer_options.clone(),
        profile_options: r.profile_options.clone(),
    })?;
    if r.options.homology_flank == 0
        || r.options.homology_flank > primers::MAX_OLIGO
        || r.options.max_oligo_length < r.options.primer_options.min_length
        || r.options.max_oligo_length > primers::MAX_OLIGO
        || r.options.min_paired_bases == 0
        || r.options.max_comparisons == 0
        || r.options.max_comparisons > overlap_extension::MAX_COMPARISONS
        || r.options.max_candidates == 0
        || r.options.max_candidates > overlap_extension::MAX_CANDIDATES
    {
        return Err(bad("invalid shared extension/total-oligo/search limits"));
    }
    crate::run_material::material_origins_preflight(
        base.material_bindings.iter().map(|binding| &binding.origin),
    )?;
    Ok(())
}
fn frame_fragment(contribution: &crate::assembly_intent::Contribution) -> Result<Fragment, Error> {
    match &contribution.selector {
        crate::assembly_intent::Selector::ReferenceRange { .. } => {
            contribution.reference_fragment()
        }
        crate::assembly_intent::Selector::MaterialStrandRange { start, end, .. } => Ok(Fragment {
            role: contribution.role.clone(),
            input_id: contribution.input_id.clone(),
            start: *start,
            end: *end,
            reverse: false,
        }),
    }
}
fn strand(contribution: &crate::assembly_intent::Contribution) -> crate::chemistry::Strand {
    match &contribution.selector {
        crate::assembly_intent::Selector::ReferenceRange { .. } => crate::chemistry::Strand::Top,
        crate::assembly_intent::Selector::MaterialStrandRange { strand, .. } => *strand,
    }
}
// These borrowed views are inspection/ranking coordinates, never native stock.
// Only pcr_source::Template::fresh can authorize physical saved-template work.
fn frame_views<'a>(
    inputs: &[(&str, &'a Molecule)],
    r: &Request,
    evidence: &'a [SourceFrame],
) -> Result<BTreeMap<String, crate::analysis::TemplateView<'a>>, Error> {
    if evidence.len() > MAX_ROWS
        || evidence
            .iter()
            .try_fold(0usize, |sum, frame| {
                sum.checked_add(frame.sequence.len())
                    .filter(|n| *n <= 16 << 20)
            })
            .is_none()
    {
        return Err(Error::TooLarge);
    }
    let mut used = BTreeSet::new();
    let mut result = BTreeMap::new();
    for contribution in &r.intent.base.base.fragments {
        let source = inputs
            .iter()
            .find(|(id, _)| *id == contribution.input_id)
            .ok_or_else(|| bad("missing original frame source"))?
            .1;
        let template_strand = strand(contribution);
        let binding = r
            .intent
            .base
            .base
            .material_bindings
            .iter()
            .find(|b| b.input_id == contribution.input_id);
        let (sequence, topology) = if binding.is_some() {
            let (index, frame) = evidence
                .iter()
                .enumerate()
                .find(|(_, frame)| {
                    frame.input_id == contribution.input_id
                        && frame.template_strand == template_strand
                })
                .ok_or_else(|| bad("native source frame evidence missing"))?;
            if !hex(&frame.state_revision)
                || frame.sequence.is_empty()
                || frame.sequence.alphabet() != collagene_model::Alphabet::Dna
                || frame.topology != source.topology
                || template_strand == crate::chemistry::Strand::Top
                    && frame.sequence != source.sequence
            {
                return Err(bad("invalid documentary source frame"));
            }
            if let crate::assembly_intent::Selector::MaterialStrandRange {
                expected_state_revision,
                ..
            } = &contribution.selector
            {
                if expected_state_revision != &frame.state_revision {
                    return Err(bad("source frame state differs from declared selector"));
                }
            }
            used.insert(index);
            (&frame.sequence, frame.topology)
        } else {
            if template_strand != crate::chemistry::Strand::Top {
                return Err(bad("bottom frame requires explicit native source"));
            }
            (&source.sequence, source.topology)
        };
        result.insert(
            contribution.role.clone(),
            crate::analysis::TemplateView {
                sequence,
                topology,
                template_strand,
            },
        );
    }
    if used.len() != evidence.len() {
        return Err(bad("duplicate or unused native frame evidence"));
    }
    Ok(result)
}
fn source_frames(
    resolved: &crate::assembly_intent::Resolved<'_>,
    r: &Request,
) -> Result<Vec<SourceFrame>, Error> {
    let mut result = vec![];
    let mut seen = BTreeSet::new();
    let mut bases = 0usize;
    for contribution in &r.intent.base.base.fragments {
        let Some(material) = resolved.material(&contribution.input_id) else {
            continue;
        };
        let template_strand = strand(contribution);
        if !seen.insert((
            contribution.input_id.clone(),
            template_strand == crate::chemistry::Strand::Bottom,
        )) {
            continue;
        }
        let duplex = &material.tracked().duplex;
        let sequence = match template_strand {
            crate::chemistry::Strand::Top => &duplex.top,
            crate::chemistry::Strand::Bottom => &duplex.bottom,
        };
        bases = bases
            .checked_add(sequence.len())
            .filter(|n| *n <= 16 << 20)
            .ok_or(Error::TooLarge)?;
        result.push(SourceFrame {
            input_id: contribution.input_id.clone(),
            template_strand,
            state_revision: crate::modified_duplex::summary(material.tracked())?.state_revision,
            sequence: sequence.clone(),
            topology: duplex.topology,
        });
    }
    Ok(result)
}
#[derive(Clone)]
struct RowGeometry {
    fragment: Fragment,
    edits: Vec<Edit>,
    edited_length: usize,
}
fn base_at(source: &crate::analysis::TemplateView<'_>, f: &Fragment, offset: usize) -> u8 {
    let local = if f.reverse {
        f.end - f.start - 1 - offset
    } else {
        offset
    };
    let raw = source.sequence.as_str().as_bytes()[(f.start + local) % source.sequence.len()];
    if f.reverse {
        collagene_model::complement(raw, collagene_model::Alphabet::Dna)
    } else {
        raw
    }
}
fn edited_prefix(
    source: &crate::analysis::TemplateView<'_>,
    row: &RowGeometry,
    count: usize,
) -> Result<String, Error> {
    let mut out = String::with_capacity(count);
    let mut cursor = 0usize;
    for e in &row.edits {
        for offset in cursor..e.start {
            if out.len() == count {
                return Ok(out);
            }
            out.push(base_at(source, &row.fragment, offset) as char);
        }
        for b in e.replacement.bytes() {
            if out.len() == count {
                return Ok(out);
            }
            out.push(b as char);
        }
        cursor = e.end;
    }
    for offset in cursor..row.fragment.end - row.fragment.start {
        if out.len() == count {
            return Ok(out);
        }
        out.push(base_at(source, &row.fragment, offset) as char);
    }
    if out.len() != count {
        return Err(bad(
            "downstream edited retained core cannot supply the requested complete homology prefix",
        ));
    }
    Ok(out)
}
fn geometry(
    sources: &BTreeMap<String, crate::analysis::TemplateView<'_>>,
    r: &Request,
) -> Result<(Vec<RowGeometry>, Vec<TailDeclaration>), Error> {
    shape(r)?;
    let base = &r.intent.base.base;
    let lengths = base
        .fragments
        .iter()
        .map(|f| {
            sources
                .get(f.role.as_str())
                .map(|m| m.sequence.len())
                .ok_or_else(|| bad("missing bound source"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let adjusted =
        assembly_junction_edit::intent_adjustment(base, &r.intent.base.junction_edits, &lengths)?;
    let mut rows = Vec::with_capacity(base.fragments.len());
    let mut inserts = vec![String::new(); base.fragments.len()];
    for j in &r.intent.base.junction_edits {
        let text = if j.inserted.bytes().all(|b| b.is_ascii_whitespace()) {
            String::new()
        } else {
            Sequence::dna(&j.inserted)?.as_str().to_owned()
        };
        if !concrete(&text) {
            return Err(bad(
                "ambiguous literal junction DNA cannot supply concrete ordered oligos",
            ));
        }
        if !base.circular && j.junction + 1 == base.fragments.len() && !text.is_empty() {
            return Err(bad(
                "a linear terminal insert requires an explicit extra-tail policy; the last core has no outgoing assembly tail",
            ));
        }
        inserts[j.junction] = text;
    }
    for (original, f) in base.fragments.iter().zip(adjusted) {
        let original = frame_fragment(original)?;
        let f = frame_fragment(&f)?;
        let source = &sources[f.role.as_str()];
        let turn = if source.topology == Topology::Circular && f.start < original.start {
            source.sequence.len()
        } else {
            0
        };
        let begin = if original.reverse {
            original
                .end
                .checked_sub(f.end.checked_add(turn).ok_or(Error::TooLarge)?)
        } else {
            f.start
                .checked_add(turn)
                .and_then(|n| n.checked_sub(original.start))
        }
        .ok_or_else(|| bad("retained contribution changed original coordinate direction"))?;
        let mut edits = Vec::new();
        for e in r.intent.edits.iter().filter(|e| e.role == f.role) {
            let mut local = e.clone();
            local.start = local
                .start
                .checked_sub(begin)
                .ok_or_else(|| bad("edit removed by incoming trim"))?;
            local.end = local
                .end
                .checked_sub(begin)
                .ok_or_else(|| bad("edit removed by incoming trim"))?;
            edits.push(local);
        }
        for e in r
            .intent
            .base
            .substitutions
            .iter()
            .filter(|e| e.role == f.role)
        {
            let start = e
                .offset
                .checked_sub(begin)
                .ok_or_else(|| bad("substitution removed by incoming trim"))?;
            edits.push(Edit {
                role: e.role.clone(),
                start,
                end: start + 1,
                expected: e.expected.clone(),
                replacement: e.replacement.clone(),
            });
        }
        edits.sort_by_key(|e| (e.start, e.end));
        let mut edited_length = f.end - f.start;
        for e in &edits {
            if e.end > f.end - f.start || !concrete(&e.expected) || !concrete(&e.replacement) {
                return Err(bad(
                    "internal edit exceeds the retained core or uses ambiguous oligo DNA",
                ));
            }
            edited_length = edited_length
                .checked_sub(e.end - e.start)
                .and_then(|n| n.checked_add(e.replacement.len()))
                .ok_or(Error::TooLarge)?;
        }
        if edited_length == 0 {
            return Err(bad(
                "a fully deleted source contribution requires an explicit role removal",
            ));
        }
        for offset in 0..f.end - f.start {
            if !b"ACGT".contains(&base_at(source, &f, offset)) {
                return Err(bad("selected PCR cores require concrete original calls"));
            }
        }
        rows.push(RowGeometry {
            fragment: f,
            edits,
            edited_length,
        });
    }
    let mut tails = Vec::with_capacity(rows.len());
    for (i, row) in rows.iter().enumerate() {
        let outgoing = i + 1 < rows.len() || base.circular;
        let (homology, downstream_role) = if outgoing {
            let next = &rows[(i + 1) % rows.len()];
            if next.edited_length <= r.overlap {
                return Err(bad(
                    "downstream edited core must exceed the shared homology",
                ));
            }
            (
                edited_prefix(&sources[next.fragment.role.as_str()], next, r.overlap)?,
                Some(next.fragment.role.clone()),
            )
        } else {
            (String::new(), None)
        };
        let added = format!("{}{}", inserts[i], homology);
        let tail = if added.is_empty() {
            String::new()
        } else {
            Sequence::dna(&added)?
                .reverse_complement()?
                .as_str()
                .to_owned()
        };
        if tail
            .len()
            .checked_add(r.options.primer_options.min_length)
            .is_none_or(|n| n > r.options.max_oligo_length)
        {
            return Err(bad(
                "complete junction insert/homology plus annealing domain exceeds total oligo length",
            ));
        }
        tails.push(TailDeclaration {
            role: row.fragment.role.clone(),
            input_id: row.fragment.input_id.clone(),
            forward_tail: if row.fragment.reverse {
                tail.clone()
            } else {
                String::new()
            },
            reverse_tail: if row.fragment.reverse {
                String::new()
            } else {
                tail
            },
            outgoing_insert: inserts[i].clone(),
            homology,
            downstream_role,
        });
    }
    Ok((rows, tails))
}
fn row_request(
    r: &Request,
    row: &RowGeometry,
    tail: &TailDeclaration,
    template_strand: crate::chemistry::Strand,
) -> extension::Request {
    let mut options = r.options.clone();
    options.outer_forward_tail = tail.forward_tail.clone();
    options.outer_reverse_tail = tail.reverse_tail.clone();
    extension::Request {
        template_strand,
        material_binding: None,
        fragment: row.fragment.clone(),
        edits: row.edits.clone(),
        options,
    }
}
fn exact_options(
    r: &Request,
    row: &RowGeometry,
    tail: &TailDeclaration,
) -> Result<DesignOptions, Error> {
    let mut o = r.options.primer_options.clone();
    let length = row.fragment.end - row.fragment.start;
    if length < o.min_length.checked_mul(2).ok_or(Error::TooLarge)? {
        return Err(bad(
            "ordinary PCR core cannot fit both nonoverlapping annealing domains",
        ));
    }
    o.target_start = row.fragment.start;
    o.target_end = row.fragment.end;
    o.forward_tail = tail.forward_tail.clone();
    o.reverse_tail = tail.reverse_tail.clone();
    o.max_length = o
        .max_length
        .min(length - o.min_length)
        .min(r.options.max_oligo_length - o.forward_tail.len().max(o.reverse_tail.len()));
    if o.max_length < o.min_length {
        return Err(bad(
            "actual assembly tail leaves no oligo length satisfying the native minimum",
        ));
    }
    Ok(o)
}
fn exact_work(
    source: crate::pcr_source::Template<'_>,
    o: &DesignOptions,
) -> Result<(usize, usize), Error> {
    let count = o.max_length - o.min_length + 1;
    let width = o.max_length + o.forward_tail.len().max(o.reverse_tail.len());
    let pair = count
        .checked_mul(count)
        .and_then(|n| n.checked_mul(width))
        .and_then(|n| n.checked_mul(width))
        .and_then(|n| n.checked_mul(2))
        .ok_or(Error::TooLarge)?;
    if pair > 20_000_000 {
        return Err(bad(
            "native exact-pair ranking exceeds its per-call bound; narrow lengths",
        ));
    }
    let metrics = count
        .checked_mul(width)
        .and_then(|n| n.checked_mul(width))
        .and_then(|n| n.checked_mul(12))
        .ok_or(Error::TooLarge)?;
    let scan = source
        .comparison_bases()?
        .checked_mul(o.max_length)
        .and_then(|n| n.checked_mul(8))
        .ok_or(Error::TooLarge)?;
    let material = (o.target_end - o.target_start)
        .checked_add(o.forward_tail.len() + o.reverse_tail.len())
        .and_then(|n| n.checked_mul(2))
        .ok_or(Error::TooLarge)?;
    Ok((
        pair.checked_add(metrics)
            .and_then(|n| n.checked_add(scan))
            .ok_or(Error::TooLarge)?,
        material,
    ))
}
fn reserved_scope(
    sources: &BTreeMap<String, crate::analysis::TemplateView<'_>>,
    r: &Request,
    rows: &[RowGeometry],
    tails: &[TailDeclaration],
) -> Result<(usize, usize), Error> {
    let mut work = 0usize;
    let mut material = 0usize;
    for (row, tail) in rows.iter().zip(tails) {
        let source = &sources[row.fragment.role.as_str()];
        let multiplicity = if row.edits.is_empty() { 1 } else { 3 };
        let scan = source
            .sequence
            .len()
            .checked_mul(r.options.primer_options.max_length)
            .and_then(|n| n.checked_mul(8 * multiplicity))
            .ok_or(Error::TooLarge)?;
        let copied = row
            .edited_length
            .checked_add(tail.forward_tail.len() + tail.reverse_tail.len())
            .and_then(|n| n.checked_mul(if row.edits.is_empty() { 2 } else { 6 }))
            .ok_or(Error::TooLarge)?;
        work = work
            .checked_add(scan)
            .and_then(|n| n.checked_add(copied))
            .filter(|n| *n <= MAX_WORK)
            .ok_or(Error::TooLarge)?;
        material = material.checked_add(copied).ok_or(Error::TooLarge)?;
    }
    // ProfiledOverlap scans every prepared core for homology/repeats. Reserve
    // the whole bounded pass once, not an independent100M allocation per row.
    let prepared = rows
        .iter()
        .zip(tails)
        .try_fold(0usize, |sum, (row, tail)| {
            sum.checked_add(row.edited_length)
                .and_then(|n| n.checked_add(tail.forward_tail.len() + tail.reverse_tail.len()))
                .ok_or(Error::TooLarge)
        })?;
    let junctions = rows.len() - 1 + usize::from(r.intent.base.base.circular);
    let profile = prepared
        .checked_mul(2)
        .and_then(|n| n.checked_mul(r.overlap))
        .and_then(|n| n.checked_add(r.overlap))
        .and_then(|n| n.checked_mul(junctions))
        .ok_or(Error::TooLarge)?
        .min(r.profile_options.max_comparisons);
    let finishing_work = match &r.finishing {
        assembly_recipe::OverlapFinishing::ReferenceOnly => 0,
        assembly_recipe::OverlapFinishing::HomologyAssembly { policy } => policy.max_comparisons,
        assembly_recipe::OverlapFinishing::HomologyAnnealing { policy } => policy.max_comparisons,
    };
    work = work
        .checked_add(finishing_work)
        .ok_or(Error::TooLarge)?
        .checked_add(profile.checked_mul(2).ok_or(Error::TooLarge)?)
        .filter(|n| *n <= MAX_WORK)
        .ok_or(Error::TooLarge)?;
    material = material
        .checked_mul(2)
        .filter(|n| *n <= 64 * 1024 * 1024)
        .ok_or(Error::TooLarge)?;
    Ok((work, material))
}
fn proposal_preflight(inputs: &[(&str, &Molecule)], r: &Request) -> Result<(), Error> {
    let views = frame_views(inputs, r, &[])?;
    let (rows, tails) = geometry(&views, r)?;
    let sources = inputs.iter().copied().collect::<BTreeMap<_, _>>();
    let (mut work, mut material) = reserved_scope(&views, r, &rows, &tails)?;
    for (row, tail) in rows.iter().zip(&tails) {
        let source =
            crate::pcr_source::Template::reference(sources[row.fragment.input_id.as_str()])?;
        let charge = if row.edits.is_empty() {
            exact_work(source, &exact_options(r, row, tail)?)?
        } else {
            extension::source_work_preflight(source, &row_request(r, row, tail, source.strand))?
        };
        work = work
            .checked_add(charge.0)
            .filter(|n| *n <= MAX_WORK)
            .ok_or(Error::TooLarge)?;
        material = material
            .checked_add(charge.1)
            .filter(|n| *n <= 64 * 1024 * 1024)
            .ok_or(Error::TooLarge)?;
    }
    size(&(inputs, r))
}
/// Attest saved sources documentarily before primer search or origin copying.
/// Native recipe verification, not this ownership check, supplies fresh chemistry.
pub fn source_preflight(
    inputs: &[(&str, &Molecule)],
    r: &Request,
    workspace: Option<&Workspace>,
) -> Result<(), Error> {
    shape(r)?;
    if inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    collagene_model::bounds::model_preflight(&inputs.iter().map(|(_, m)| *m).collect::<Vec<_>>())?;
    size(&(inputs, r))?;
    assembly_edit_intent::source_preflight(inputs, &r.intent, workspace)?;
    if r.intent.base.base.material_bindings.is_empty() {
        proposal_preflight(inputs, r)?;
    }
    Ok(())
}
pub fn preflight(
    inputs: &[PlanInput],
    r: &Request,
    workspace: Option<&Workspace>,
) -> Result<(), Error> {
    source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        r,
        workspace,
    )?;
    if inputs
        .iter()
        .any(|i| !hex(&i.revision) || i.revision != i.molecule.revision())
    {
        return Err(bad(
            "source input identity/revision differs from its exact original molecule",
        ));
    }
    Ok(())
}
fn context(inputs: &[PlanInput], r: &Request) -> Result<String, Error> {
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        verification::ENGINE_VERSION,
        inputs
            .iter()
            .map(|i| (&i.id, &i.revision))
            .collect::<Vec<_>>(),
        r,
    ))
}
pub fn design_revision(
    inputs: &[PlanInput],
    r: &Request,
    workspace: Option<&Workspace>,
) -> Result<String, Error> {
    preflight(inputs, r, workspace)?;
    context(inputs, r)
}
/// Derive literal tails without primer search or assembly. Actual saved strands
/// require fresh origin reconstruction; callers must run this native phase off the UI thread.
pub fn expected_tails(
    inputs: &[PlanInput],
    r: &Request,
    workspace: Option<&Workspace>,
) -> Result<Vec<TailDeclaration>, Error> {
    preflight(inputs, r, workspace)?;
    let borrowed = inputs
        .iter()
        .map(|i| (i.id.as_str(), &i.molecule))
        .collect::<Vec<_>>();
    if r.intent.base.base.material_bindings.is_empty() {
        return geometry(&frame_views(&borrowed, r, &[])?, r).map(|(_, tails)| tails);
    }
    let resolved = crate::assembly_intent::resolve(inputs, &r.intent.base.base, workspace)?;
    assembly_edit_intent::build_resolved(&resolved, &r.intent)?;
    let evidence = source_frames(&resolved, r)?;
    geometry(&frame_views(&borrowed, r, &evidence)?, r).map(|(_, tails)| tails)
}
pub fn tail_preflight(
    inputs: &[PlanInput],
    r: &Request,
    workspace: Option<&Workspace>,
) -> Result<Vec<TailDeclaration>, Error> {
    expected_tails(inputs, r, workspace)
}
fn ordinary(
    source: crate::pcr_source::Template<'_>,
    r: &Request,
    row: &RowGeometry,
    tail: &TailDeclaration,
    budget: &mut extension::Budget,
) -> Result<(Preparation, RowDesign), Error> {
    let options = exact_options(r, row, tail)?;
    let (work, material) = exact_work(source, &options)?;
    if !budget.reserve(work) || !budget.reserve_material(material) {
        return Err(Error::TooLarge);
    }
    let design = extension::propose(source, &options, "exact", &row.fragment.role)?;
    let (search, _) = extension::select(
        source,
        &design,
        "exact",
        options.target_start,
        options.target_end - options.target_start,
        budget,
    );
    let pair = design
        .pairs
        .get(search.selected_index)
        .ok_or_else(|| bad("native exact pair disappeared"))?;
    let preparation = Preparation::ExactPcr {
        forward: pair.forward.primer.clone(),
        reverse: pair.reverse.primer.clone(),
        template_strand: source.strand,
    };
    Ok((
        preparation,
        RowDesign {
            role: row.fragment.role.clone(),
            input_id: row.fragment.input_id.clone(),
            input_revision: source.reference.revision(),
            fragment: row.fragment.clone(),
            template_strand: source.strand,
            edits: vec![],
            designs: vec![design],
            search: vec![search],
        },
    ))
}
fn from_origin(preparation: Preparation, origin: &MaterialOrigin) -> Result<Preparation, Error> {
    let child = match preparation {
        Preparation::ExactPcr {
            forward,
            reverse,
            template_strand,
        } => RunPreparation::ExactPcr {
            forward,
            reverse,
            template_strand,
        },
        Preparation::OverlapExtensionPcr { request } => {
            RunPreparation::OverlapExtensionPcr { request }
        }
        _ => {
            return Err(bad(
                "whole primer proposal emitted an unsupported preparation method",
            ));
        }
    };
    Ok(match origin {
        MaterialOrigin::Run(origin) => Preparation::FromRun {
            origin: origin.clone(),
            preparation: child,
        },
        MaterialOrigin::Derived(origin) => Preparation::FromDerived {
            origin: origin.clone(),
            preparation: child,
        },
    })
}
pub fn design(
    inputs: &[PlanInput],
    r: &Request,
    workspace: Option<&Workspace>,
) -> Result<Report, Error> {
    let context_revision = design_revision(inputs, r, workspace)?;
    let borrowed = inputs
        .iter()
        .map(|i| (i.id.as_str(), &i.molecule))
        .collect::<Vec<_>>();
    let resolved = crate::assembly_intent::resolve(inputs, &r.intent.base.base, workspace)?;
    let intent = assembly_edit_intent::build_resolved(&resolved, &r.intent)?;
    let source_frames = source_frames(&resolved, r)?;
    let views = frame_views(&borrowed, r, &source_frames)?;
    let (geometry, tails) = geometry(&views, r)?;
    let sources = inputs
        .iter()
        .map(|i| (i.id.as_str(), &i.molecule))
        .collect::<BTreeMap<_, _>>();
    let (work, material) = reserved_scope(&views, r, &geometry, &tails)?;
    let mut budget = extension::Budget {
        used: work,
        limit: MAX_WORK,
        materialized: material,
    };
    // Intent is a separate source/literal construction before molecular proposals.
    let mut rows = Vec::with_capacity(geometry.len());
    let mut fragments = Vec::with_capacity(geometry.len());
    for (row, tail) in geometry.iter().zip(&tails) {
        let original = sources[row.fragment.input_id.as_str()];
        let contribution = r
            .intent
            .base
            .base
            .fragments
            .iter()
            .find(|c| c.role == row.fragment.role)
            .ok_or(Error::TooLarge)?;
        let source = if let Some(material) = resolved.material(&row.fragment.input_id) {
            crate::pcr_source::Template::fresh(original, material.tracked(), strand(contribution))?
        } else {
            crate::pcr_source::Template::reference(original)?
        };
        let (mut preparation, design) = if row.edits.is_empty() {
            ordinary(source, r, row, tail, &mut budget)?
        } else {
            let proposal = extension::compile_source_with_budget(
                source,
                &row_request(r, row, tail, source.strand),
                &mut budget,
            )?;
            (
                Preparation::OverlapExtensionPcr {
                    request: proposal.preparation,
                },
                RowDesign {
                    role: row.fragment.role.clone(),
                    input_id: row.fragment.input_id.clone(),
                    input_revision: proposal.source_revision,
                    fragment: row.fragment.clone(),
                    template_strand: source.strand,
                    edits: row.edits.clone(),
                    designs: proposal.designs,
                    search: proposal.search,
                },
            )
        };
        if let Some(origin) = r
            .intent
            .base
            .base
            .material_bindings
            .iter()
            .find(|binding| binding.input_id == row.fragment.input_id)
        {
            preparation = from_origin(preparation, &origin.origin)?;
        }
        fragments.push(FragmentSpec {
            generation: None,
            role: row.fragment.role.clone(),
            input_id: row.fragment.input_id.clone(),
            preparation,
            reverse: row.fragment.reverse,
        });
        rows.push(design);
    }
    let base = &r.intent.base.base;
    let recipe = assembly_recipe::Recipe {
        name: base.name.clone(),
        fragments,
        overlaps: vec![r.overlap; geometry.len() - 1 + usize::from(base.circular)],
        circular: base.circular,
        options: r.profile_options.clone(),
    };
    let preview = r.finishing.preview(inputs, &intent.target, &recipe)?;
    let report=Report{engine_version:ENGINE_VERSION.into(),context_revision,source_revisions:inputs.iter().map(|i|i.revision.clone()).collect(),source_frames,request:r.clone(),intent,recipe,tails,rows,preview,warnings:vec![
        "Every source contribution, junction trim/insert, SNP and internal edit is a literal declaration built independently of molecular reconstruction. Outgoing homology is taken from the next edited retained source core, never chosen from a successful reconstructed product or arbitrary target.".into(),
        "Ordinary rows use source-bound ranked exact PCR; edited rows use two PCR branches, physical reciprocal extension and outer PCR. One shared100M comparison-work/material ledger covers native primer/PCR/junction searches and the declared final whole-recipe molecular search reservation; no role receives an independent100M allowance. Literal-reference/annotation construction and closed-origin replay retain their separate native aggregate/copy/replay bounds;100M is not an absolute end-to-end runtime bound.".into(),
        "Per-row ranked calls and search notes are proposal evidence, not a separate molecular assessment. The one final preview replays declared exact saved origins and retains originating/terminal/PCR uncertainty. Component selection and reference agreement do not establish isolation or experimental success.".into(),
        format!("Shared charged comparison work{} of{}; materialized/reserved base units{} of64MiB. Native ranked lists and caps are retained per role; no global search optimality is asserted.",budget.used,MAX_WORK,budget.materialized),
    ]};
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
fn rc(text: &str) -> Result<String, Error> {
    if text.is_empty() {
        Ok(String::new())
    } else {
        Ok(Sequence::dna(text)?
            .reverse_complement()?
            .as_str()
            .to_owned())
    }
}
fn stored_segment(
    source: &crate::analysis::TemplateView<'_>,
    f: &Fragment,
    start: usize,
    end: usize,
) -> String {
    (start..end)
        .map(|j| source.sequence.as_str().as_bytes()[(f.start + j) % source.sequence.len()] as char)
        .collect()
}
fn edited_bridge(
    source: &crate::analysis::TemplateView<'_>,
    row: &RowGeometry,
) -> Result<(usize, usize, String), Error> {
    let start = row
        .edits
        .iter()
        .map(|e| e.start)
        .min()
        .ok_or_else(|| bad("missing internal bridge edit"))?;
    let end = row
        .edits
        .iter()
        .map(|e| e.end)
        .max()
        .ok_or(Error::TooLarge)?;
    let mut bridge = String::new();
    let mut cursor = start;
    for e in &row.edits {
        for j in cursor..e.start {
            bridge.push(base_at(source, &row.fragment, j) as char);
        }
        bridge.push_str(&e.replacement);
        cursor = e.end;
    }
    for j in cursor..end {
        bridge.push(base_at(source, &row.fragment, j) as char);
    }
    if row.fragment.reverse {
        Ok((
            row.fragment.end - row.fragment.start - end,
            row.fragment.end - row.fragment.start - start,
            rc(&bridge)?,
        ))
    } else {
        Ok((start, end, bridge))
    }
}
fn call_options(
    source: &crate::analysis::TemplateView<'_>,
    r: &Request,
    row: &RowGeometry,
    tail: &TailDeclaration,
) -> Result<Vec<DesignOptions>, Error> {
    if row.edits.is_empty() {
        return Ok(vec![exact_options(r, row, tail)?]);
    }
    let (s, e, insert) = edited_bridge(source, row)?;
    let len = row.fragment.end - row.fragment.start;
    let h = r.options.homology_flank;
    let inner = h.checked_add(insert.len()).ok_or(Error::TooLarge)?;
    if h >= s
        || h >= len - e
        || s < 2 * r.options.primer_options.min_length
        || len - e < 2 * r.options.primer_options.min_length
        || inner + r.options.primer_options.min_length > r.options.max_oligo_length
    {
        return Err(bad(
            "native unchanged flanks do not fit the explicit edited bridge",
        ));
    }
    let mut left = r.options.primer_options.clone();
    let mut right = left.clone();
    let mut outer = left.clone();
    left.target_start = row.fragment.start;
    left.target_end = row.fragment.start + s;
    right.target_start = (row.fragment.start + e) % source.sequence.len();
    right.target_end = right.target_start + len - e;
    outer.target_start = row.fragment.start;
    outer.target_end = row.fragment.end;
    left.reverse_tail = rc(&format!(
        "{insert}{}",
        stored_segment(source, &row.fragment, e, e + h)
    ))?;
    right.forward_tail = format!(
        "{}{insert}",
        stored_segment(source, &row.fragment, s - h, s)
    );
    outer.forward_tail = tail.forward_tail.clone();
    outer.reverse_tail = tail.reverse_tail.clone();
    left.max_length = left
        .max_length
        .min(s - left.min_length)
        .min(r.options.max_oligo_length - inner);
    right.max_length = right
        .max_length
        .min(len - e - right.min_length)
        .min(r.options.max_oligo_length - inner);
    outer.max_length = outer
        .max_length
        .min(s)
        .min(len - e)
        .min(r.options.max_oligo_length - tail.forward_tail.len().max(tail.reverse_tail.len()));
    Ok(vec![left, right, outer])
}
fn location_guard(
    location: &collagene_model::Location,
    n: usize,
    maximum: usize,
) -> Result<(), Error> {
    let mut pending = vec![(location, 0usize)];
    let mut count = 0usize;
    let mut bases = 0usize;
    while let Some((l, depth)) = pending.pop() {
        count += 1;
        if count > 8 || depth > 4 {
            return Err(Error::TooLarge);
        }
        match l {
            collagene_model::Location::Range {
                start,
                end,
                before: false,
                after: false,
            } if start < end && *end <= n => {
                bases = bases
                    .checked_add(end - start)
                    .filter(|n| *n <= maximum)
                    .ok_or(Error::TooLarge)?;
            }
            collagene_model::Location::Complement(inner) => pending.push((inner, depth + 1)),
            collagene_model::Location::Join(parts) if parts.len() <= 2 => {
                for p in parts {
                    pending.push((p, depth + 1));
                }
            }
            _ => return Err(bad("invalid bounded native primer location")),
        }
    }
    Ok(())
}
fn validate_design(
    source: &crate::analysis::TemplateView<'_>,
    source_revision: &str,
    d: &DesignReport,
    search: &PairSearch,
    o: &DesignOptions,
    call: &str,
    maximum: usize,
) -> Result<(), Error> {
    if d.source_revision != source_revision
        || d.template_strand != source.template_strand
        || d.options != *o
        || d.pairs.is_empty()
        || d.pairs.len() > o.max_pairs
        || d.evaluated_pairs > 4096
        || d.rejected_pairs > d.evaluated_pairs
        || d.warnings.len() > 32
        || d.warnings.iter().any(|s| s.len() > 8192)
        || search.call != call
        || search.candidates != d.pairs.len()
        || search.examined > search.candidates
        || search.rejected > search.examined
        || search.selected_index >= d.pairs.len()
        || search.notes.len() > 22
        || search.notes.iter().any(|s| s.len() > 4096)
        || search.unique_geometry
            && (search.examined != search.selected_index + 1
                || search.rejected + 1 != search.examined)
        || !search.unique_geometry
            && (search.selected_index != 0 || search.rejected != search.examined)
        || d.pairs.windows(2).any(|p| p[0].score > p[1].score)
    {
        return Err(bad(
            "native ranked call/search/source correspondence changed",
        ));
    }
    for pair in &d.pairs {
        if !pair.score.is_finite()
            || !pair.tm_difference_c.is_finite()
            || pair.forward.primer.annealing.len() + pair.reverse.primer.annealing.len()
                > o.target_end - o.target_start
            || pair.product_length
                != (o.target_end - o.target_start) + o.forward_tail.len() + o.reverse_tail.len()
        {
            return Err(bad("invalid native ranked pair geometry"));
        }
        for (candidate, reverse, tail) in [
            (&pair.forward, false, &o.forward_tail),
            (&pair.reverse, true, &o.reverse_tail),
        ] {
            let p = &candidate.primer;
            overlap_extension_pcr::primer(p)?;
            if p.annealing.len() < o.min_length
                || p.annealing.len() > o.max_length
                || p.tail != *tail
                || p.ordered_sequence().len() > maximum
                || candidate.analysis.name != p.name
                || candidate.analysis.annealing != p.annealing.as_str()
                || candidate.analysis.tail != p.tail
                || candidate.analysis.ordered_sequence != p.ordered_sequence()
            {
                return Err(bad(
                    "returned ordered oligo differs from its source-bound native call",
                ));
            }
            location_guard(
                &candidate.location,
                source.sequence.len(),
                p.annealing.len(),
            )?;
            let start = if reverse {
                o.target_end - p.annealing.len()
            } else {
                o.target_start
            };
            let actual = stored_segment(
                source,
                &Fragment {
                    role: String::new(),
                    input_id: String::new(),
                    start: 0,
                    end: source.sequence.len(),
                    reverse: false,
                },
                start,
                start + p.annealing.len(),
            );
            let actual = if reverse { rc(&actual)? } else { actual };
            if p.annealing.as_str() != actual
                || candidate.location.extract(&source.sequence)? != actual
            {
                return Err(bad(
                    "returned oligo no longer matches the declared original endpoint",
                ));
            }
        }
    }
    Ok(())
}
fn selected_preparation(
    source: &crate::analysis::TemplateView<'_>,
    r: &Request,
    row: &RowGeometry,
    d: &RowDesign,
) -> Result<Preparation, Error> {
    let selected = |i: usize| {
        d.designs
            .get(i)
            .and_then(|call| call.pairs.get(d.search.get(i)?.selected_index))
            .ok_or_else(|| bad("selected native pair is unavailable"))
    };
    if row.edits.is_empty() {
        let p = selected(0)?;
        return Ok(Preparation::ExactPcr {
            forward: p.forward.primer.clone(),
            reverse: p.reverse.primer.clone(),
            template_strand: d.template_strand,
        });
    }
    let (_, _, insert) = edited_bridge(source, row)?;
    let pair = |i| {
        let p = selected(i)?;
        Ok::<_, Error>(overlap_extension_pcr::PrimerPair {
            forward: p.forward.primer.clone(),
            reverse: p.reverse.primer.clone(),
            template_strand: if i == 2 {
                crate::chemistry::Strand::Top
            } else {
                d.template_strand
            },
        })
    };
    Ok(Preparation::OverlapExtensionPcr {
        request: overlap_extension_pcr::Request {
            left: pair(0)?,
            right: pair(1)?,
            outer: pair(2)?,
            extension: overlap_extension::Options {
                overlap: 2 * r.options.homology_flank + insert.len(),
                min_paired_bases: r.options.min_paired_bases,
                max_comparisons: r.options.max_comparisons,
                max_candidates: r.options.max_candidates,
            },
        },
    })
}
impl Report {
    /// Bounded documentary correspondence. No primer ranking, PCR, whole-intent
    /// rebuilding, closed-origin replay or molecular transformation runs here.
    pub fn validate(&self) -> Result<(), Error> {
        let n = self.request.intent.base.base.fragments.len();
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.context_revision)
            || n == 0
            || n > MAX_ROWS
            || self.rows.len() != n
            || self.tails.len() != n
            || self.recipe.fragments.len() != n
            || self.source_revisions.len() > MAX_ROWS
            || self.source_revisions.is_empty()
            || self.source_revisions.iter().any(|s| !hex(s))
            || self.warnings.len() > 32
            || self.warnings.iter().any(|s| s.len() > 8192)
            || self
                .rows
                .iter()
                .any(|row| row.designs.len() > 3 || row.search.len() > 3)
        {
            return Err(bad("invalid bounded whole-proposal report shape"));
        }
        shape(&self.request)?;
        size(self)?;
        let inputs = &self.preview.compiled.plan.inputs;
        if self.request.intent.base.base.material_bindings.is_empty() {
            proposal_preflight(
                &inputs
                    .iter()
                    .map(|i| (i.id.as_str(), &i.molecule))
                    .collect::<Vec<_>>(),
                &self.request,
            )?;
        }
        if inputs.iter().any(|i| i.revision != i.molecule.revision()) {
            return Err(bad("captured proposal source revision changed"));
        }
        if self.source_revisions
            != inputs
                .iter()
                .map(|i| i.revision.clone())
                .collect::<Vec<_>>()
            || context(inputs, &self.request)? != self.context_revision
        {
            return Err(bad("whole proposal source/request context changed"));
        }
        if assembly_recipe::hash(&self.intent.request)?
            != assembly_recipe::hash(&self.request.intent)?
        {
            return Err(bad("independent typed intent declaration changed"));
        }
        self.intent
            .validate_binding(&self.intent.context_revision)?;
        let borrowed = inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>();
        let views = frame_views(&borrowed, &self.request, &self.source_frames)?;
        let (geometry, tails) = geometry(&views, &self.request)?;
        if tails != self.tails {
            return Err(bad("literal source/edit/junction tail provenance changed"));
        }
        let sources = inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<BTreeMap<_, _>>();
        let base = &self.request.intent.base.base;
        if self.recipe.name != base.name
            || self.recipe.circular != base.circular
            || self.recipe.overlaps
                != vec![self.request.overlap; n - 1 + usize::from(base.circular)]
            || assembly_recipe::hash(&self.recipe.options)?
                != assembly_recipe::hash(&self.request.profile_options)?
        {
            return Err(bad(
                "whole recipe topology/junction/profile declaration changed",
            ));
        }
        for (((row, derived), tail), fragment) in self
            .rows
            .iter()
            .zip(&geometry)
            .zip(&tails)
            .zip(&self.recipe.fragments)
        {
            let original = sources[derived.fragment.input_id.as_str()];
            let source = &views[derived.fragment.role.as_str()];
            if row.role != derived.fragment.role
                || row.input_id != derived.fragment.input_id
                || row.input_revision != original.revision()
                || row.template_strand != source.template_strand
                || assembly_recipe::hash(&row.fragment)?
                    != assembly_recipe::hash(&derived.fragment)?
                || row.edits != derived.edits
                || fragment.role != row.role
                || fragment.input_id != row.input_id
                || fragment.reverse != row.fragment.reverse
            {
                return Err(bad(
                    "whole row original identity/adjusted interval/edit coordinates changed",
                ));
            }
            let expected = call_options(source, &self.request, derived, tail)?;
            if row.designs.len() != expected.len() || row.search.len() != expected.len() {
                return Err(bad("row native method evidence count changed"));
            }
            for (i, ((call, search), options)) in row
                .designs
                .iter()
                .zip(&row.search)
                .zip(&expected)
                .enumerate()
            {
                let label = if expected.len() == 1 {
                    "exact"
                } else {
                    ["left", "right", "outer"][i]
                };
                validate_design(
                    source,
                    &original.revision(),
                    call,
                    search,
                    options,
                    label,
                    self.request.options.max_oligo_length,
                )?;
            }
            let mut preparation = selected_preparation(source, &self.request, derived, row)?;
            if let Some(origin) = self
                .request
                .intent
                .base
                .base
                .material_bindings
                .iter()
                .find(|binding| binding.input_id == row.input_id)
            {
                preparation = from_origin(preparation, &origin.origin)?;
            }
            if assembly_recipe::hash(&preparation)? != assembly_recipe::hash(&fragment.preparation)?
            {
                return Err(bad(
                    "recipe differs from its exact selected oligos or closed-origin declaration",
                ));
            }
        }
        self.preview
            .validate_binding(&self.request.finishing.preview_revision(
                inputs,
                &self.intent.target,
                &self.recipe,
            )?)?;
        if self.preview.compiled.plan.target_revision != self.intent.target.revision() {
            return Err(bad("whole recipe replaced its independent target"));
        }
        let (mut expected_steps, expected_rows) = assembly_recipe::preparation_steps(
            inputs,
            &self.recipe.fragments,
            crate::assembly_recipe::CleavagePolicy::AssumedComplete,
        )?;
        let assembly_step = expected_steps.len();
        let fragments = expected_rows
            .iter()
            .map(|row| verification::OrientedSource {
                source: row.prepared_source.clone(),
                reverse: row.reverse,
            })
            .collect();
        expected_steps.push(match &self.request.finishing {
            assembly_recipe::OverlapFinishing::ReferenceOnly => Operation::ProfiledOverlap {
                fragments,
                overlaps: self.recipe.overlaps.clone(),
                circular: self.recipe.circular,
                options: self.recipe.options.clone(),
            },
            assembly_recipe::OverlapFinishing::HomologyAssembly { policy } => {
                Operation::HomologyAssembly {
                    fragments,
                    options: self.recipe.options.clone(),
                    policy: policy.clone(),
                }
            }
            assembly_recipe::OverlapFinishing::HomologyAnnealing { policy } => {
                Operation::HomologyAnnealing {
                    fragments,
                    options: self.recipe.options.clone(),
                    policy: policy.clone(),
                }
            }
        });
        if self.preview.compiled.assembly_step != assembly_step
            || self.preview.compiled.plan.name != self.recipe.name
            || assembly_recipe::hash(&expected_rows)?
                != assembly_recipe::hash(&self.preview.compiled.rows)?
            || assembly_recipe::hash(&expected_steps)?
                != assembly_recipe::hash(&self.preview.compiled.plan.steps)?
        {
            return Err(bad(
                "whole compiled graph differs from its declared original-source recipe",
            ));
        }
        Ok(())
    }
    /// Exact current workspace/source acceptance without repeating primer search
    /// or PCR/reconstructing material from the documentary preview.
    pub fn validate_sources(
        &self,
        inputs: &[PlanInput],
        expected: &str,
        workspace: Option<&Workspace>,
    ) -> Result<(), Error> {
        if design_revision(inputs, &self.request, workspace)? != expected {
            return Err(bad("current source/origin/options differ from proposal"));
        }
        self.validate_binding(expected)?;
        let retained = &self.preview.compiled.plan.inputs;
        if retained.len() != inputs.len()
            || retained
                .iter()
                .zip(inputs)
                .any(|(a, b)| a.id != b.id || a.revision != b.revision || a.molecule != b.molecule)
        {
            return Err(bad("proposal contains substituted source input"));
        }
        self.intent
            .validate_sources(inputs, &self.intent.context_revision, workspace)
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.context_revision != expected {
            return Err(bad("stale whole source/intent/options/origin context"));
        }
        self.validate()
    }
}
