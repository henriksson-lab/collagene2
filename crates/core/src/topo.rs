//! Source-backed preparation declarations and fresh physical TOPO reconstruction.
//! No sequence-record metadata establishes activation, end chemistry or kit state.
use crate::{Alphabet, Error, Molecule, Sequence, Topology, analysis, chemistry, tailing};
pub use collagene_topo as physical;
use physical::projection;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const ENGINE_VERSION: &str = "topo-source-adapter-1";
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
dto! {pub struct PreparationEvidence {pub substrates:physical::Substrates,pub maps:Vec<PreparedMap>,pub not_represented:Vec<ReferenceGap>,pub assumptions:Vec<String>}}
dto! {pub struct TopoReport {pub engine_version:String,pub batch_revision:String,pub sources:Vec<SourceIdentity>,pub preparation:PreparationEvidence,pub physical:physical::Report}}
dto! {pub struct Selection {pub batch_revision:String,pub component:usize,pub reference_seed:physical::Call}}
dto! {pub struct Product {pub engine_version:String,pub batch_revision:String,pub selection:Selection,pub sources:Vec<SourceIdentity>,pub preparation:PreparationEvidence,pub projection:projection::ProjectionReport,pub reference:Option<Molecule>}}
dto! {#[derive(Default)] pub struct AcceptancePolicy {pub allow_nicks:bool,pub allow_unknown_terminal_chemistry:bool,pub allow_unapplied_transfers:bool,pub allow_other_components:bool,pub allow_circular_rotation:bool,pub allow_reverse_complement:bool}}
dto! {pub enum TargetOutcome {AcceptedUnderDeclaredGraph,Rejected{reasons:Vec<String>},Withheld{reason:String}}}
dto! {pub struct TargetAssessment {pub engine_version:String,pub batch_revision:String,pub target_revision:String,pub policy:AcceptancePolicy,pub product:Product,pub outcome:TargetOutcome}}
// Additive envelopes: Product/TargetAssessment remain unchanged for old wire
// consumers. These fields describe correspondence, never chemical ancestry.
dto! {pub struct AnnotatedProduct {pub engine_version:String,pub batch_revision:String,pub physical:Product,pub annotations:Option<crate::topo_annotations::AnnotationEvidence>,pub reference:Option<Molecule>}}
dto! {pub struct AnnotatedAssessment {pub engine_version:String,pub batch_revision:String,pub target_revision:String,pub policy:AcceptancePolicy,pub product:AnnotatedProduct,pub outcome:TargetOutcome}}

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
struct AnnotationHash {
    hash: blake3::Hasher,
    bytes: usize,
}
impl Write for AnnotationHash {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|&n| n <= MAX_REPORT_BYTES)
            .ok_or_else(|| std::io::Error::other("bounded annotated TOPO binding"))?;
        self.hash.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn annotation_hash(value: &impl Serialize) -> Result<String, Error> {
    let mut writer = AnnotationHash {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(&mut writer, value).map_err(|_| Error::TooLarge)?;
    Ok(writer.hash.finalize().to_hex().to_string())
}
fn annotated_product_binding(physical_revision: &str) -> Result<String, Error> {
    annotation_hash(&(
        ANNOTATION_ENGINE_VERSION,
        "annotated-product",
        physical_revision,
    ))
}
fn annotated_target_binding(
    product_revision: &str,
    target_revision: &str,
    policy: &AcceptancePolicy,
) -> Result<String, Error> {
    annotation_hash(&(
        ANNOTATION_ENGINE_VERSION,
        "annotated-target",
        product_revision,
        target_revision,
        policy,
    ))
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
/// Borrowed aggregate guard before recursive validation, model hashes, sequence
/// construction or rich input cloning. Bounds apply to native callers as well.
pub fn source_preflight(inputs: &[Input], r: &Request) -> Result<(), Error> {
    source_preflight_inner(inputs, r, &BTreeSet::new())
}
fn source_preflight_inner(
    inputs: &[Input],
    r: &Request,
    reconstructed_templates: &BTreeSet<String>,
) -> Result<(), Error> {
    if inputs.is_empty()
        || inputs.len() > 32
        || r.recipe.schema_version != 1
        || r.recipe.materials.is_empty()
        || r.recipe.materials.len() > 16
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
pub(crate) fn prepare(inputs: &[Input], r: &Request) -> Result<PreparationEvidence, Error> {
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
    })
}
pub fn simulate(inputs: &[Input], r: &Request) -> Result<TopoReport, Error> {
    source_preflight(inputs, r)?;
    simulate_preflighted(inputs, r)
}
fn simulate_preflighted(inputs: &[Input], r: &Request) -> Result<TopoReport, Error> {
    let batch_revision = batch_binding(inputs, r)?;
    let preparation = prepare(inputs, r)?;
    let physical = physical::simulate(&preparation.substrates, &r.transfers, &r.engine_limits)
        .map_err(engine_error)?;
    let result = TopoReport {
        engine_version: ENGINE_VERSION.into(),
        batch_revision,
        sources: identities(inputs),
        preparation,
        physical,
    };
    bounded(&result, MAX_REPORT_BYTES)?;
    Ok(result)
}
pub fn product_revision(inputs: &[Input], r: &Request, s: &Selection) -> Result<String, Error> {
    let batch = batch_revision(inputs, r)?;
    if s.batch_revision.len() != 64 || s.batch_revision != batch || s.component >= 64 {
        return Err(bad("selection does not bind current source/recipe"));
    }
    hash(&(ENGINE_VERSION, "selected-physical-component", batch, s))
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
    let batch_revision = hash(&(ENGINE_VERSION, "selected-physical-component", batch, s))?;
    let fresh = simulate_preflighted(inputs, r)?;
    let select = projection::Selection {
        source_batch_revision: fresh.physical.batch_revision.clone(),
        component: s.component,
        reference_seed: s.reference_seed,
    };
    let projection = projection::project(
        &fresh.preparation.substrates,
        &r.transfers,
        &r.engine_limits,
        &select,
        &r.projection_options,
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
        m.metadata.insert(
            "topo:adapter_batch_revision".into(),
            fresh.batch_revision.clone(),
        );
        m.metadata.insert(
            RECIPE_KEY.into(),
            serde_json::to_string(&(r, s, &fresh.sources)).map_err(|_| Error::TooLarge)?,
        );
        Some(m)
    } else {
        None
    };
    let p = Product {
        engine_version: ENGINE_VERSION.into(),
        batch_revision,
        selection: s.clone(),
        sources: fresh.sources,
        preparation: fresh.preparation,
        projection,
        reference,
    };
    bounded(&p, MAX_REPORT_BYTES)?;
    Ok(p)
}
/// Borrowed source/request/selection binding only; no graph reconstruction or
/// annotation projection. Annotation changes already alter the physical source
/// revision, and this separate domain pins annotation-transform semantics.
pub fn annotated_product_revision(
    inputs: &[Input],
    request: &Request,
    selection: &Selection,
) -> Result<String, Error> {
    // Reuse the existing bounded streaming implementation of the legacy
    // physical product domain. It hashes borrowed inputs and dispatches no job.
    annotated_product_binding(&crate::topo_worker::projection_revision(
        inputs, request, selection,
    )?)
}
/// This private handoff accepts only a Product freshly built in this module.
/// External serialized Product values are never accepted as annotation authority.
fn annotate_fresh_product(inputs: &[Input], physical: Product) -> Result<AnnotatedProduct, Error> {
    let batch_revision = annotated_product_binding(&physical.batch_revision)?;
    let (annotations, reference) = match (&physical.projection.outcome, &physical.reference) {
        (projection::Outcome::Proven(view), Some(destination)) => {
            let evidence =
                crate::topo_annotations::project(inputs, &physical.preparation, view, destination)?;
            let mut reference = destination.clone();
            reference.features = evidence.features.clone();
            reference.metadata.insert(
                "topo:annotation_engine".into(),
                ANNOTATION_ENGINE_VERSION.into(),
            );
            reference
                .metadata
                .insert("topo:annotation_revision".into(), batch_revision.clone());
            reference.metadata.insert("topo:annotation_scope".into(),"Source-coordinate correspondence only; not inherited chemical ancestry or biological validity. Unmapped/review evidence remains in the full annotated product.".into());
            reference.validate()?;
            (Some(evidence), Some(reference))
        }
        (projection::Outcome::Withheld(_), None) => (None, None),
        _ => return Err(bad("fresh contour/reference presence disagrees")),
    };
    let result = AnnotatedProduct {
        engine_version: ANNOTATION_ENGINE_VERSION.into(),
        batch_revision,
        physical,
        annotations,
        reference,
    };
    bounded(&result, MAX_REPORT_BYTES)?;
    Ok(result)
}
/// Fresh physical reconstruction plus separately retained conservative annotation
/// correspondence. The original physical Product stays annotation-free unchanged.
pub fn annotated_product(
    inputs: &[Input],
    request: &Request,
    selection: &Selection,
) -> Result<AnnotatedProduct, Error> {
    let fresh = product(inputs, request, selection)?;
    annotate_fresh_product(inputs, fresh)
}
/// Independently bound target/policy request; this hashes borrowed data only.
pub fn annotated_target_revision(
    inputs: &[Input],
    request: &Request,
    selection: &Selection,
    target: &Molecule,
    policy: &AcceptancePolicy,
) -> Result<String, Error> {
    target_preflight(target)?;
    annotated_target_binding(
        &annotated_product_revision(inputs, request, selection)?,
        &annotation_hash(target)?,
        policy,
    )
}
/// Reconstruct/check the physical target exactly once. Annotation correspondence
/// enriches that same fresh graph, without changing its target-policy outcome.
pub fn check_target_annotated(
    inputs: &[Input],
    request: &Request,
    selection: &Selection,
    target: &Molecule,
    policy: &AcceptancePolicy,
) -> Result<AnnotatedAssessment, Error> {
    let TargetAssessment {
        product,
        policy,
        target_revision,
        outcome,
        ..
    } = check_target(inputs, request, selection, target, policy)?;
    let product = annotate_fresh_product(inputs, product)?;
    let batch_revision =
        annotated_target_binding(&product.batch_revision, &target_revision, &policy)?;
    let result = AnnotatedAssessment {
        engine_version: ANNOTATION_ENGINE_VERSION.into(),
        batch_revision,
        target_revision,
        policy,
        product,
        outcome,
    };
    bounded(&result, MAX_REPORT_BYTES)?;
    Ok(result)
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
    let preparation = prepare(inputs, r)?;
    if vector
        .iter()
        .chain(&insert)
        .any(|i| *i >= preparation.substrates.strands.len())
    {
        return Err(bad("orientation strand index outside prepared material"));
    }
    let reference = physical::orientations(
        &preparation.substrates,
        vector,
        insert,
        directional,
        &r.engine_limits,
    )
    .map_err(engine_error)?;
    let mut candidate_requests = Vec::new();
    for candidate in &reference.candidates {
        let mut request = r.clone();
        request.transfers = candidate.steps.iter().map(|s| s.transfer.clone()).collect();
        candidate_requests.push(request);
    }
    let sources = identities(inputs);
    let batch_revision = hash(&(
        ENGINE_VERSION,
        "two-declared-orientations",
        &sources,
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
        Task::Project { request, selection } => Reply::Projected(annotate_fresh_product(
            inputs,
            product_preflighted(inputs, request, selection)?,
        )?),
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
            let product = annotate_fresh_product(inputs, product)?;
            Reply::Checked(AnnotatedAssessment {
                engine_version: ANNOTATION_ENGINE_VERSION.into(),
                batch_revision: annotated_target_binding(
                    &product.batch_revision,
                    &target_revision,
                    &policy,
                )?,
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
