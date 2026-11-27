//! Independently declared reference contributions from explicit source domains.
//! Native source replay identifies actual input strands; no planned assembly or
//! desired target supplies their bases, chemistry or source authority.
use crate::{
    Alphabet, Error, Molecule, Sequence, Topology,
    annotations::{self, MappingSegment},
    assembly_recipe,
    chemistry::Strand,
    run_lineage::EventBinding,
    run_material::{self, MaterialOrigin},
    verification::PlanInput,
    workspace::Workspace,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE_VERSION: &str = "assembly-intent-2";
pub const MATERIAL_METADATA: &str = "assembly_material_sources:v1";
const MATERIAL_ENGINE: &str = "actual-strand-intent-1";
pub(crate) const MAX_BASES: usize = 8 << 20;
const MAX_ROWS: usize = 32;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)]pub enum Selector{
 ReferenceRange{start:usize,end:usize,reverse:bool},
 MaterialStrandRange{strand:Strand,start:usize,end:usize,expected_state_revision:String},
}}
dto! {#[derive(PartialEq,Eq)]pub struct Contribution{pub role:String,pub input_id:String,pub selector:Selector}}
dto! {pub struct MaterialBinding{pub input_id:String,pub event:EventBinding,pub origin:MaterialOrigin}}
dto! {pub struct Request{pub name:String,pub fragments:Vec<Contribution>,pub circular:bool,pub material_bindings:Vec<MaterialBinding>}}
dto! {#[derive(PartialEq,Eq)]pub struct StrandDomain{pub strand:Strand,pub start:i64,pub end:i64,pub length:usize,pub topology:Topology}}
dto! {#[derive(PartialEq,Eq)]pub enum ContributionPosition{Reference{position:usize,reverse:bool},Material{strand:Strand,offset:usize}}}
dto! {pub struct PlacedFragment{
 pub role:String,pub input_id:String,pub input_revision:String,pub selector:Selector,pub source_length:usize,
 pub source_event:Option<EventBinding>,pub source_state_revision:Option<String>,pub source_domain:Option<StrandDomain>,pub incomplete:bool,
 /// Full original oriented contribution, BEFORE trims, substitutions and indels.
 pub original_sequence:Sequence,pub retained_start:usize,pub retained_end:usize,
 /// Original source top→original contribution LOCAL coordinates, not edited target coordinates.
 pub top_mappings:Vec<MappingSegment>,pub output_start:usize,pub output_end:usize,
}}
dto! {pub struct Report{pub context_revision:String,pub request:Request,pub target:Molecule,pub fragments:Vec<PlacedFragment>,pub warnings:Vec<String>}}
dto! {pub struct MaterialDescriptor{pub input_id:String,pub event:EventBinding,pub selection:run_material::MaterialSelection,pub origin_revision:String,pub state_revision:String,pub top_length:usize,pub bottom_start:i64,pub bottom_length:usize,pub topology:Topology,pub incomplete:bool}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("assembly intent: {s}"))
}
pub(crate) fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn token(s: &str, n: usize) -> Result<(), Error> {
    if s.trim().is_empty() || s.len() > n || s.chars().any(char::is_control) {
        Err(bad("invalid bounded role/input/name"))
    } else {
        Ok(())
    }
}
impl Selector {
    pub fn length(&self) -> Result<usize, Error> {
        let (a, b) = match self {
            Self::ReferenceRange { start, end, .. }
            | Self::MaterialStrandRange { start, end, .. } => (*start, *end),
        };
        b.checked_sub(a)
            .filter(|n| *n > 0 && *n <= MAX_BASES)
            .ok_or_else(|| bad("nonempty bounded source interval required"))
    }
    pub fn reference_interval(&self) -> Result<(usize, usize, bool), Error> {
        self.length()?;
        match self {
            Self::ReferenceRange {
                start,
                end,
                reverse,
            } => Ok((*start, *end, *reverse)),
            Self::MaterialStrandRange { .. } => Err(bad(
                "actual material-strand domains are not original-reference PCR intervals; choose explicitly supported primers/preparation",
            )),
        }
    }
    pub fn is_reverse(&self) -> bool {
        match self {
            Self::ReferenceRange { reverse, .. } => *reverse,
            Self::MaterialStrandRange { strand, .. } => *strand == Strand::Bottom,
        }
    }
    /// Trim a selected window in contribution orientation; this does not prepare DNA.
    pub fn trim(&self, left: usize, right: usize, n: usize) -> Result<Self, Error> {
        let len = self.length()?;
        if n == 0 || left.checked_add(right).is_none_or(|k| k >= len) {
            return Err(bad(
                "trims must retain at least one original contribution base",
            ));
        }
        let mut result = self.clone();
        match &mut result {
            Self::ReferenceRange {
                start,
                end,
                reverse,
            } => {
                *start = start
                    .checked_add(if *reverse { right } else { left })
                    .ok_or(Error::TooLarge)?;
                *end -= if *reverse { left } else { right };
                if *start >= n {
                    *start -= n;
                    *end -= n;
                }
            }
            Self::MaterialStrandRange { start, end, .. } => {
                *start = start.checked_add(left).ok_or(Error::TooLarge)?;
                *end -= right;
                if *start >= n {
                    *start -= n;
                    *end -= n;
                }
            }
        }
        Ok(result)
    }
}
impl Contribution {
    /// A real domain conversion, never a material/reference compatibility adapter.
    pub fn reference_fragment(&self) -> Result<crate::overlap_primer_design::Fragment, Error> {
        token(&self.role, 128)?;
        token(&self.input_id, 128)?;
        let (start, end, reverse) = self.selector.reference_interval()?;
        Ok(crate::overlap_primer_design::Fragment {
            role: self.role.clone(),
            input_id: self.input_id.clone(),
            start,
            end,
            reverse,
        })
    }
}
impl PlacedFragment {
    pub fn source_position(&self, offset: usize) -> Result<ContributionPosition, Error> {
        if offset >= self.original_sequence.len() || self.source_length == 0 {
            return Err(bad("contribution offset outside original resolved strand"));
        }
        match &self.selector {
            Selector::ReferenceRange {
                start,
                end,
                reverse,
            } => Ok(ContributionPosition::Reference {
                position: if *reverse {
                    end - 1 - offset
                } else {
                    start.checked_add(offset).ok_or(Error::TooLarge)?
                } % self.source_length,
                reverse: *reverse,
            }),
            Selector::MaterialStrandRange { strand, start, .. } => {
                Ok(ContributionPosition::Material {
                    strand: *strand,
                    offset: start.checked_add(offset).ok_or(Error::TooLarge)? % self.source_length,
                })
            }
        }
    }
}
pub(crate) fn shape(r: &Request) -> Result<(), Error> {
    token(&r.name, 1024)?;
    if r.fragments.is_empty()
        || r.fragments.len() > MAX_ROWS
        || r.material_bindings.len() > MAX_ROWS
    {
        return Err(Error::TooLarge);
    }
    let mut roles = BTreeSet::new();
    let mut size = 0usize;
    for f in &r.fragments {
        token(&f.role, 128)?;
        token(&f.input_id, 128)?;
        if !roles.insert(&f.role) {
            return Err(bad("duplicate source role"));
        }
        size = size
            .checked_add(f.selector.length()?)
            .filter(|n| *n <= MAX_BASES)
            .ok_or(Error::TooLarge)?;
        if let Selector::MaterialStrandRange {
            expected_state_revision,
            ..
        } = &f.selector
        {
            if !hex(expected_state_revision) {
                return Err(bad("actual source state revision is required"));
            }
        }
    }
    let mut ids = BTreeSet::new();
    for b in &r.material_bindings {
        token(&b.input_id, 128)?;
        if !ids.insert(&b.input_id) || !r.fragments.iter().any(|f| f.input_id == b.input_id) {
            return Err(bad("duplicate or unused material source binding"));
        }
        crate::run_material::event_shape(&b.event)?;
        if b.origin.event() != &b.event {
            return Err(bad("material origin and exact source event differ"));
        }
    }
    for f in &r.fragments {
        if matches!(f.selector, Selector::MaterialStrandRange { .. })
            && !r.material_bindings.iter().any(|b| b.input_id == f.input_id)
        {
            return Err(bad("actual strand selector needs its exact native origin"));
        }
    }
    Ok(())
}
/// Borrowed context/ownership preflight; does not replay chemistry or copy source material.
pub fn preflight(
    inputs: &[(&str, &Molecule)],
    r: &Request,
    w: Option<&Workspace>,
) -> Result<(), Error> {
    shape(r)?;
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    crate::model::bounds::model_preflight(&inputs.iter().map(|(_, m)| *m).collect::<Vec<_>>())?;
    run_material::material_origins_preflight(r.material_bindings.iter().map(|b| &b.origin))?;
    assembly_recipe::hash(&(inputs, r))?;
    let mut seen = BTreeSet::new();
    let mut bases = 0usize;
    for (id, m) in inputs {
        token(id, 128)?;
        if !seen.insert(*id) || m.sequence.alphabet() != Alphabet::Dna {
            return Err(bad("distinct DNA input IDs required"));
        }
        m.validate()?;
        bases = bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= MAX_BASES)
            .ok_or(Error::TooLarge)?;
    }
    for f in &r.fragments {
        let m = inputs
            .iter()
            .find(|(id, _)| *id == f.input_id)
            .ok_or_else(|| bad("missing original input"))?
            .1;
        if let Selector::ReferenceRange { start, end, .. } = f.selector {
            range(start, end, m.sequence.len(), m.topology)?;
        }
    }
    if !r.material_bindings.is_empty() {
        let w =
            w.ok_or_else(|| bad("actual material strands require native workspace attestation"))?;
        run_material::validate_native_origins(
            w,
            r.material_bindings
                .iter()
                .map(|binding| match &binding.origin {
                    MaterialOrigin::Run(origin) => run_material::NativeOriginRef::Run(origin),
                    MaterialOrigin::Derived(origin) => {
                        run_material::NativeOriginRef::Derived(origin)
                    }
                }),
        )?;
        for b in &r.material_bindings {
            let m = inputs
                .iter()
                .find(|(id, _)| *id == b.input_id)
                .ok_or_else(|| bad("material source input missing"))?
                .1;
            if run_material::resolve(w, &b.event)? != m || m.revision() != b.event.revision {
                return Err(bad(
                    "native source event differs from the complete original input snapshot",
                ));
            }
        }
    }
    Ok(())
}
fn range(start: usize, end: usize, n: usize, t: Topology) -> Result<(), Error> {
    if start >= n
        || end.checked_sub(start).is_none_or(|k| k == 0 || k > n)
        || start.checked_add(n).is_none_or(|limit| end > limit)
        || t == Topology::Linear && end > n
    {
        return Err(bad(
            "interval exceeds the actual selected source strand; only circles permit one unrolled turn",
        ));
    }
    Ok(())
}
pub fn intent_revision(
    inputs: &[PlanInput],
    r: &Request,
    w: Option<&Workspace>,
) -> Result<String, Error> {
    preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        r,
        w,
    )?;
    for i in inputs {
        if !hex(&i.revision) || i.molecule.revision() != i.revision {
            return Err(bad("original source content changed"));
        }
    }
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        crate::verification::ENGINE_VERSION,
        inputs
            .iter()
            .map(|i| (&i.id, &i.revision))
            .collect::<Vec<_>>(),
        // Ownership was attested above. Exact events, closed source instructions
        // and the scientific model determine these source declarations; saving
        // an attachment or advancing an unrelated head cannot change target DNA
        // identity. Host request/adoption guards still bind the workspace scope.
        r,
    ))
}
fn slice(s: &Sequence, start: usize, end: usize) -> Result<Sequence, Error> {
    let n = s.len();
    let mut out = String::with_capacity(end - start);
    for i in start..end {
        out.push(s.as_str().as_bytes()[i % n] as char);
    }
    Sequence::dna(&out)
}
fn reference_maps(start: usize, end: usize, reverse: bool, n: usize) -> Vec<MappingSegment> {
    let mut result = vec![];
    let mut done = 0;
    let len = end - start;
    while done < len {
        let a = (start + done) % n;
        let k = (n - a).min(len - done);
        result.push(MappingSegment {
            source_start: a,
            source_end: a + k,
            destination_start: if reverse { len - done - k } else { done },
            reverse,
        });
        done += k;
    }
    result
}
fn bottom_maps(
    start: usize,
    end: usize,
    duplex: &crate::chemistry::Duplex,
) -> Result<Vec<MappingSegment>, Error> {
    bottom_domain_maps(
        start,
        end,
        duplex.top.len(),
        duplex.bottom_start,
        duplex.bottom.len(),
    )
}
fn bottom_domain_maps(
    start: usize,
    end: usize,
    t: usize,
    bottom_start: i64,
    b: usize,
) -> Result<Vec<MappingSegment>, Error> {
    let boundary = bottom_start
        .checked_add(i64::try_from(b).map_err(|_| Error::TooLarge)?)
        .ok_or(Error::TooLarge)?;
    let mut out = vec![];
    let mut done = 0;
    while done < end - start {
        let local = (start + done) % b;
        let k = (b - local).min(end - start - done);
        let hi = boundary - local as i64;
        let lo = hi - k as i64;
        let a = lo.max(0).min(t as i64);
        let z = hi.max(0).min(t as i64);
        if a < z {
            out.push(MappingSegment {
                source_start: a as usize,
                source_end: z as usize,
                destination_start: done + (hi - z) as usize,
                reverse: true,
            });
        }
        done += k;
    }
    Ok(out)
}
pub(crate) struct Built {
    pub target: Molecule,
    pub fragments: Vec<PlacedFragment>,
    pub warnings: Vec<String>,
}
pub(crate) struct Resolved<'a> {
    pub inputs: &'a [PlanInput],
    pub rows: Vec<PlacedFragment>,
    pub descriptors: Vec<MaterialDescriptor>,
    pub warnings: Vec<String>,
    pub context: String,
    materials: BTreeMap<String, run_material::FreshMaterial>,
    material_keys: BTreeMap<String, String>,
}
impl Resolved<'_> {
    /// Lend the exact freshly replayed source, shared by target construction and
    /// primer planning. No imported report or selected target supplies this state.
    pub(crate) fn material(&self, input_id: &str) -> Option<&run_material::FreshMaterial> {
        self.material_keys
            .get(input_id)
            .and_then(|key| self.materials.get(key))
    }
}
pub(crate) fn resolve<'a>(
    inputs: &'a [PlanInput],
    r: &Request,
    w: Option<&Workspace>,
) -> Result<Resolved<'a>, Error> {
    let context = intent_revision(inputs, r, w)?;
    let mut cache = BTreeMap::new();
    let mut material_keys = BTreeMap::new();
    let mut material = 0usize;
    let mut descriptors = vec![];
    let mut warnings=vec!["Independently declared reference only. Native source replay identifies actual input strand calls; no planned digest/PCR/assembly product or desired target selects these bases. This is not a PCR/digest/assembly prediction, physical material or experimental success.".into()];
    for b in &r.material_bindings {
        let key = assembly_recipe::hash(&(&b.event, b.origin.selection()))?;
        if !cache.contains_key(&key) {
            let source = &inputs
                .iter()
                .find(|i| i.id == b.input_id)
                .ok_or_else(|| bad("missing material source"))?
                .molecule;
            let fresh = b.origin.replay(source)?;
            material = material
                .checked_add(fresh.tracked().duplex.top.len())
                .and_then(|n| n.checked_add(fresh.tracked().duplex.bottom.len()))
                .filter(|n| *n <= 16 * 1024 * 1024)
                .ok_or(Error::TooLarge)?;
            cache.insert(key.clone(), fresh);
        }
        let fresh = &cache[&key];
        material_keys.insert(b.input_id.clone(), key);
        let d = &fresh.tracked().duplex;
        let state = crate::modified_duplex::summary(fresh.tracked())?.state_revision;
        descriptors.push(MaterialDescriptor {
            input_id: b.input_id.clone(),
            event: b.event.clone(),
            selection: b.origin.selection(),
            origin_revision: assembly_recipe::hash(&b.origin)?,
            state_revision: state,
            top_length: d.top.len(),
            bottom_start: d.bottom_start,
            bottom_length: d.bottom.len(),
            topology: d.topology,
            incomplete: fresh.incomplete(),
        });
        if fresh.incomplete() {
            warnings.push(format!("Input {} retains incomplete native-source reconstruction; source selection is not experimental isolation.",b.input_id));
        }
        for finding in fresh.findings().iter().take(2) {
            warnings.push(format!(
                "Input {} origin: {}",
                b.input_id,
                finding.chars().take(2048).collect::<String>()
            ));
        }
    }
    let mut rows = vec![];
    for f in &r.fragments {
        let i = inputs
            .iter()
            .find(|i| i.id == f.input_id)
            .ok_or_else(|| bad("original input absent"))?;
        let (sequence, length, event, state, domain, incomplete, maps) = match &f.selector {
            Selector::ReferenceRange {
                start,
                end,
                reverse,
            } => {
                let s = slice(&i.molecule.sequence, *start, *end)?;
                let s = if *reverse { s.reverse_complement()? } else { s };
                (
                    s,
                    i.molecule.sequence.len(),
                    None,
                    None,
                    None,
                    false,
                    reference_maps(*start, *end, *reverse, i.molecule.sequence.len()),
                )
            }
            Selector::MaterialStrandRange {
                strand,
                start,
                end,
                expected_state_revision,
            } => {
                let b = r
                    .material_bindings
                    .iter()
                    .find(|b| b.input_id == f.input_id)
                    .ok_or_else(|| bad("missing native binding"))?;
                let key = assembly_recipe::hash(&(&b.event, b.origin.selection()))?;
                let fresh = &cache[&key];
                let d = &fresh.tracked().duplex;
                let descriptor = descriptors
                    .iter()
                    .find(|d| d.input_id == f.input_id)
                    .ok_or(Error::TooLarge)?;
                if &descriptor.state_revision != expected_state_revision {
                    return Err(bad("actual native strand/modification state changed"));
                }
                let s = if *strand == Strand::Top {
                    &d.top
                } else {
                    &d.bottom
                };
                range(*start, *end, s.len(), d.topology)?;
                let offset = if *strand == Strand::Top {
                    0
                } else {
                    d.bottom_start
                };
                let maps = if *strand == Strand::Top {
                    reference_maps(*start, *end, false, d.top.len())
                } else {
                    bottom_maps(*start, *end, d)?
                };
                (
                    slice(s, *start, *end)?,
                    s.len(),
                    Some(b.event.clone()),
                    Some(descriptor.state_revision.clone()),
                    Some(StrandDomain {
                        strand: *strand,
                        start: offset,
                        end: offset.checked_add(s.len() as i64).ok_or(Error::TooLarge)?,
                        length: s.len(),
                        topology: d.topology,
                    }),
                    fresh.incomplete(),
                    maps,
                )
            }
        };
        let len = sequence.len();
        rows.push(PlacedFragment {
            role: f.role.clone(),
            input_id: f.input_id.clone(),
            input_revision: i.revision.clone(),
            selector: f.selector.clone(),
            source_length: length,
            source_event: event,
            source_state_revision: state,
            source_domain: domain,
            incomplete,
            original_sequence: sequence,
            retained_start: 0,
            retained_end: len,
            top_mappings: maps,
            output_start: 0,
            output_end: len,
        });
    }
    Ok(Resolved {
        inputs,
        rows,
        descriptors,
        warnings,
        context,
        materials: cache,
        material_keys,
    })
}
/// Clip original-local mappings to retained original offsets and place in target.
pub(crate) fn clipped_maps(row: &PlacedFragment) -> Result<Vec<MappingSegment>, Error> {
    let mut maps = vec![];
    for m in &row.top_mappings {
        let len = m.source_end - m.source_start;
        let lo = m.destination_start.max(row.retained_start);
        let hi = (m.destination_start + len).min(row.retained_end);
        if lo >= hi {
            continue;
        }
        let a = lo - m.destination_start;
        let z = hi - m.destination_start;
        maps.push(MappingSegment {
            source_start: if m.reverse {
                m.source_end - z
            } else {
                m.source_start + a
            },
            source_end: if m.reverse {
                m.source_end - a
            } else {
                m.source_start + z
            },
            destination_start: row.output_start + lo - row.retained_start,
            reverse: m.reverse,
        });
    }
    Ok(maps)
}
#[derive(Serialize)]
struct PlacementDescriptor<'a> {
    role: &'a str,
    input_id: &'a str,
    input_revision: &'a str,
    selector: &'a Selector,
    source_length: usize,
    source_event: &'a Option<EventBinding>,
    source_state_revision: &'a Option<String>,
    source_domain: &'a Option<StrandDomain>,
    incomplete: bool,
    original_sequence_revision: String,
    retained_start: usize,
    retained_end: usize,
    top_mappings: &'a [MappingSegment],
    output_start: usize,
    output_end: usize,
}
pub(crate) fn placements_json(rows: &[PlacedFragment]) -> Result<String, Error> {
    let v = rows
        .iter()
        .map(|r| PlacementDescriptor {
            role: &r.role,
            input_id: &r.input_id,
            input_revision: &r.input_revision,
            selector: &r.selector,
            source_length: r.source_length,
            source_event: &r.source_event,
            source_state_revision: &r.source_state_revision,
            source_domain: &r.source_domain,
            incomplete: r.incomplete,
            original_sequence_revision: blake3::hash(r.original_sequence.as_str().as_bytes())
                .to_hex()
                .to_string(),
            retained_start: r.retained_start,
            retained_end: r.retained_end,
            top_mappings: &r.top_mappings,
            output_start: r.output_start,
            output_end: r.output_end,
        })
        .collect::<Vec<_>>();
    Ok(serde_json::to_string(&v).map_err(|_| Error::TooLarge)?)
}
pub(crate) fn request_json(r: &Request) -> Result<String, Error> {
    #[derive(Serialize)]
    struct Declaration<'a> {
        name: &'a str,
        fragments: &'a [Contribution],
        circular: bool,
        material_sources: Vec<(
            &'a str,
            &'a EventBinding,
            run_material::MaterialSelection,
            String,
        )>,
    }
    serde_json::to_string(&Declaration {
        name: &r.name,
        fragments: &r.fragments,
        circular: r.circular,
        material_sources: r
            .material_bindings
            .iter()
            .map(|b| {
                Ok((
                    b.input_id.as_str(),
                    &b.event,
                    b.origin.selection(),
                    assembly_recipe::hash(&b.origin)?,
                ))
            })
            .collect::<Result<Vec<_>, Error>>()?,
    })
    .map_err(|_| Error::TooLarge)
}
pub(crate) fn build_layout(
    resolved: &Resolved<'_>,
    r: &Request,
    layout: &[(usize, usize)],
    context: &str,
) -> Result<(Molecule, Vec<PlacedFragment>, Vec<String>), Error> {
    if layout.len() != resolved.rows.len() {
        return Err(bad("layout/source count differs"));
    }
    let mut dna = String::new();
    let mut rows = resolved.rows.clone();
    for (row, (a, b)) in rows.iter_mut().zip(layout) {
        if a >= b || *b > row.original_sequence.len() {
            return Err(bad("retained layout outside original contribution"));
        }
        row.retained_start = *a;
        row.retained_end = *b;
        row.output_start = dna.len();
        dna.push_str(&row.original_sequence.as_str()[*a..*b]);
        row.output_end = dna.len();
    }
    let mut target = Molecule::dna(
        &r.name,
        &dna,
        if r.circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )?;
    let mut warnings = resolved.warnings.clone();
    for row in &rows {
        let input = resolved
            .inputs
            .iter()
            .find(|i| i.id == row.input_id)
            .ok_or(Error::TooLarge)?;
        if input.molecule.features.len() > 4096usize.saturating_sub(target.features.len()) {
            warnings.push(format!("Role {}: annotations withheld by the 4096-feature budget; original source annotations remain unchanged.", row.role));
            continue;
        }
        let mapped = match annotations::project_features_with_topology(
            &input.molecule,
            &clipped_maps(row)?,
            target.sequence.len(),
            target.topology,
        ) {
            Ok(mapped) => mapped,
            Err(error) => {
                warnings.push(format!("Role {}: annotation projection withheld ({error}); original source evidence remains unchanged.", row.role));
                continue;
            }
        };
        let mut uncertain = false;
        for m in &row.top_mappings {
            let count = m.source_end - m.source_start;
            let original = &input.molecule.sequence.as_str()[m.source_start..m.source_end];
            let actual =
                &row.original_sequence.as_str()[m.destination_start..m.destination_start + count];
            let expected = if m.reverse {
                Sequence::dna(original)?
                    .reverse_complement()?
                    .as_str()
                    .to_owned()
            } else {
                original.to_owned()
            };
            uncertain |= actual != expected;
        }
        let mut features = mapped.features;
        for feature in &mut features {
            if uncertain {
                feature.needs_review = true;
                feature
                    .qualifiers
                    .retain(|(k, _)| k != "translation" && k != "codon_start");
                feature.qualifiers.push(("note".into(),Some("Actual strand calls differ from displayed-top ambiguity; review source sequence/function assertions.".into())));
            }
        }
        let mut budget = crate::model::bounds::ModelBudget::default();
        if budget
            .molecule(&target)
            .and_then(|()| budget.features(&features))
            .is_err()
        {
            warnings.push(format!("Role {}: projected annotations exceed the aggregate native-feature/qualifier budget; original source evidence remains unchanged.", row.role));
            continue;
        }
        target.features.extend(features);
        if !mapped.unmapped_features.is_empty() || mapped.dropped_features > 0 {
            warnings.push(format!("Role {}: {} annotations lack a complete source-top mapping; bottom-only calls receive no invented annotation.",row.role,mapped.unmapped_features.len()+mapped.dropped_features));
        }
        for warning in mapped.warnings.iter().take(8) {
            warnings.push(format!(
                "Role {}: {}",
                row.role,
                warning.chars().take(1024).collect::<String>()
            ));
        }
    }

    target.metadata.insert(
        "duplex_model".into(),
        "Reference-only assembly design intent; physical material state unavailable".into(),
    );
    target
        .metadata
        .insert("assembly:intent:reference_only".into(), "true".into());
    target
        .metadata
        .insert("assembly_intent:engine".into(), ENGINE_VERSION.into());
    target
        .metadata
        .insert("assembly_intent:context_revision".into(), context.into());
    target
        .metadata
        .insert("assembly_intent:request".into(), request_json(r)?);
    target.metadata.insert(
        "assembly_intent:contributions".into(),
        placements_json(&rows)?,
    );
    target.metadata.insert(
        "assembly_intent:sequence_blake3".into(),
        blake3::hash(target.sequence.as_str().as_bytes())
            .to_hex()
            .to_string(),
    );
    if !resolved.descriptors.is_empty() {
        let manifest = serde_json::to_string(&resolved.descriptors).map_err(|_| Error::TooLarge)?;
        target.metadata.insert(
            "assembly_material_sources:engine".into(),
            MATERIAL_ENGINE.into(),
        );
        target.metadata.insert(
            "assembly_material_sources:revision".into(),
            blake3::hash(manifest.as_bytes()).to_hex().to_string(),
        );
        target.metadata.insert(MATERIAL_METADATA.into(), manifest);
    }
    if dna.bytes().any(|b| !b"ACGT".contains(&b)) {
        warnings.push("IUPAC ambiguity is retained as requested reference notation, without establishing concrete primer calls or unique material.".into());
    }
    crate::model::bounds::model_preflight(&[&target])?;
    Ok((target, rows, warnings))
}
pub fn build(inputs: &[PlanInput], r: &Request, w: Option<&Workspace>) -> Result<Report, Error> {
    let resolved = resolve(inputs, r, w)?;
    let layout = resolved
        .rows
        .iter()
        .map(|r| (0, r.original_sequence.len()))
        .collect::<Vec<_>>();
    let (target, fragments, warnings) = build_layout(&resolved, r, &layout, &resolved.context)?;
    let report = Report {
        context_revision: resolved.context,
        request: r.clone(),
        target,
        fragments,
        warnings,
    };
    report.validate()?;
    Ok(report)
}
/// One current documentary descriptor parser. Hash/shape validation proves
/// neither native parent ownership nor material state. Sharing MUST attest each
/// exact event in the target's real native history ancestry and the selected
/// native Run/Derived owner before following any reference.
pub fn material_source_descriptors(m: &Molecule) -> Result<Vec<MaterialDescriptor>, Error> {
    let present = m
        .metadata
        .keys()
        .any(|k| k.starts_with("assembly_material_sources:"));
    if !present {
        return Ok(vec![]);
    }
    let bytes = m
        .metadata
        .get(MATERIAL_METADATA)
        .ok_or_else(|| bad("partial material source descriptor"))?;
    if bytes.len() > 256 * 1024
        || m.metadata
            .get("assembly_material_sources:engine")
            .map(String::as_str)
            != Some(MATERIAL_ENGINE)
        || m.metadata.get("assembly_material_sources:revision")
            != Some(&blake3::hash(bytes.as_bytes()).to_hex().to_string())
    {
        return Err(bad("material source descriptor changed"));
    }
    let records: Vec<MaterialDescriptor> =
        serde_json::from_str(bytes).map_err(|_| bad("invalid material source descriptor"))?;
    if records.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    for r in &records {
        token(&r.input_id, 128)?;
        crate::run_material::event_shape(&r.event)?;
        if r.selection.event() != &r.event
            || !hex(&r.origin_revision)
            || !hex(&r.state_revision)
            || r.top_length == 0
            || r.bottom_length == 0
            || r.top_length > crate::model::MAX_SEQUENCE_LEN
            || r.bottom_length > crate::model::MAX_SEQUENCE_LEN
        {
            return Err(bad("invalid material source descriptor"));
        }
        r.selection.validate()?;
    }
    Ok(records)
}
pub(crate) fn validate_material_descriptors(
    m: &Molecule,
    request: &Request,
    rows: &[PlacedFragment],
) -> Result<(), Error> {
    let descriptors = material_source_descriptors(m)?;
    if descriptors.len() != request.material_bindings.len() {
        return Err(bad("material descriptor count differs from exact bindings"));
    }
    for (b, d) in request.material_bindings.iter().zip(&descriptors) {
        if d.input_id != b.input_id
            || d.event != b.event
            || d.selection != b.origin.selection()
            || d.origin_revision != assembly_recipe::hash(&b.origin)?
        {
            return Err(bad(
                "material descriptor differs from exact native source binding",
            ));
        }
        for row in rows.iter().filter(|r| {
            r.input_id == b.input_id && matches!(r.selector, Selector::MaterialStrandRange { .. })
        }) {
            let domain = row
                .source_domain
                .as_ref()
                .ok_or_else(|| bad("missing material strand domain"))?;
            let (start, length) = if domain.strand == Strand::Top {
                (0, d.top_length)
            } else {
                (d.bottom_start, d.bottom_length)
            };
            let (strand, a, b) = match row.selector {
                Selector::MaterialStrandRange {
                    strand, start, end, ..
                } => (strand, start, end),
                _ => unreachable!(),
            };
            range(a, b, length, d.topology)?;
            let maps = if strand == Strand::Top {
                reference_maps(a, b, false, d.top_length)
            } else {
                bottom_domain_maps(a, b, d.top_length, d.bottom_start, d.bottom_length)?
            };
            if row.top_mappings != maps {
                return Err(bad("material source-top correspondence changed"));
            }
            if row.source_state_revision.as_ref() != Some(&d.state_revision)
                || row.incomplete != d.incomplete
                || domain.start != start
                || domain.length != length
                || domain.topology != d.topology
            {
                return Err(bad(
                    "material descriptor differs from resolved source domain",
                ));
            }
        }
    }
    Ok(())
}
pub(crate) fn validate_rows(
    rows: &[PlacedFragment],
    request: &Request,
    target: &Molecule,
) -> Result<(), Error> {
    if rows.len() != request.fragments.len() || rows.len() > MAX_ROWS {
        return Err(bad("placement count changed"));
    }
    let mut offset = 0;
    for (row, f) in rows.iter().zip(&request.fragments) {
        if row.role != f.role
            || row.input_id != f.input_id
            || row.selector != f.selector
            || !hex(&row.input_revision)
            || row.original_sequence.alphabet() != Alphabet::Dna
            || row.original_sequence.len() != f.selector.length()?
            || row.source_length == 0
            || row.retained_start >= row.retained_end
            || row.retained_end > row.original_sequence.len()
            || row.output_start != offset
            || row.output_end <= offset
            || row.output_end > target.sequence.len()
        {
            return Err(bad("original contribution/retained placement changed"));
        }
        if let Selector::MaterialStrandRange {
            strand,
            expected_state_revision,
            ..
        } = &row.selector
        {
            let d = row
                .source_domain
                .as_ref()
                .ok_or_else(|| bad("missing actual strand domain"))?;
            if d.strand != *strand
                || d.length != row.source_length
                || d.end.checked_sub(d.start) != Some(d.length as i64)
                || row.source_state_revision.as_ref() != Some(expected_state_revision)
                || row.source_event.as_ref()
                    != request
                        .material_bindings
                        .iter()
                        .find(|b| b.input_id == row.input_id)
                        .map(|b| &b.event)
            {
                return Err(bad("actual source state/domain/event changed"));
            }
        } else if row.source_domain.is_some()
            || row.source_event.is_some()
            || row.source_state_revision.is_some()
            || row.incomplete
        {
            return Err(bad("reference range cannot claim native material evidence"));
        }
        if let Selector::ReferenceRange {
            start,
            end,
            reverse,
        } = row.selector
        {
            if row.top_mappings != reference_maps(start, end, reverse, row.source_length) {
                return Err(bad("reference source-top correspondence changed"));
            }
        }
        if row.top_mappings.len() > 4
            || row.top_mappings.iter().any(|m| {
                m.source_start >= m.source_end
                    || m.destination_start
                        .checked_add(m.source_end - m.source_start)
                        .is_none_or(|n| n > row.original_sequence.len())
            })
        {
            return Err(bad("invalid original source-top correspondence"));
        }
        offset = row.output_end;
    }
    if offset != target.sequence.len() {
        return Err(bad("placements do not cover intended DNA"));
    }
    Ok(())
}
impl Report {
    pub fn validate(&self) -> Result<(), Error> {
        shape(&self.request)?;
        crate::model::bounds::model_preflight(&[&self.target])?;
        self.target.validate()?;
        validate_rows(&self.fragments, &self.request, &self.target)?;
        if !hex(&self.context_revision)
            || self.target.name != self.request.name
            || self.target.topology
                != if self.request.circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                }
            || self.target.sequence.len() > MAX_BASES
            || self.warnings.len() > 400
            || self.warnings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad("invalid intended report"));
        }
        for row in &self.fragments {
            if row.retained_start != 0
                || row.retained_end != row.original_sequence.len()
                || self
                    .target
                    .sequence
                    .as_str()
                    .get(row.output_start..row.output_end)
                    != Some(row.original_sequence.as_str())
            {
                return Err(bad(
                    "plain target differs from original resolved contribution",
                ));
            }
        }
        if self
            .target
            .metadata
            .get("assembly_intent:engine")
            .map(String::as_str)
            != Some(ENGINE_VERSION)
            || self.target.metadata.get("assembly_intent:context_revision")
                != Some(&self.context_revision)
            || self.target.metadata.get("assembly_intent:request")
                != Some(&request_json(&self.request)?)
            || self.target.metadata.get("assembly_intent:contributions")
                != Some(&placements_json(&self.fragments)?)
            || self.target.metadata.get("assembly_intent:sequence_blake3")
                != Some(
                    &blake3::hash(self.target.sequence.as_str().as_bytes())
                        .to_hex()
                        .to_string(),
                )
            || self.target.metadata.get("duplex_model").map(String::as_str)
                != Some(
                    "Reference-only assembly design intent; physical material state unavailable",
                )
        {
            return Err(bad("intended reference receipt changed"));
        }
        validate_material_descriptors(&self.target, &self.request, &self.fragments)?;
        report_bound(self)?;
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad("source/selector context changed"));
        }
        Ok(())
    }
}

pub(crate) fn report_bound(value: &impl Serialize) -> Result<(), Error> {
    struct Count(usize);
    impl std::io::Write for Count {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.0 = self
                .0
                .checked_add(b.len())
                .filter(|n| *n <= 32 * 1024 * 1024 - 65536)
                .ok_or_else(|| {
                    std::io::Error::other("intent report exceeds reserved worker bound")
                })?;
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(Count(0), value).map_err(|_| Error::TooLarge)
}

pub(crate) fn reference_sequence(m: &Molecule, selector: &Selector) -> Result<Sequence, Error> {
    let (start, end, reverse) = selector.reference_interval()?;
    range(start, end, m.sequence.len(), m.topology)?;
    let s = slice(&m.sequence, start, end)?;
    if reverse {
        s.reverse_complement()
    } else {
        Ok(s)
    }
}
pub(crate) fn validate_reference_rows(
    inputs: &[PlanInput],
    rows: &[PlacedFragment],
) -> Result<(), Error> {
    for row in rows {
        let input = inputs
            .iter()
            .find(|i| i.id == row.input_id)
            .ok_or_else(|| bad("source input absent"))?;
        if input.revision != row.input_revision {
            return Err(bad("original source revision changed"));
        }
        if matches!(row.selector, Selector::ReferenceRange { .. })
            && (row.source_length != input.molecule.sequence.len()
                || reference_sequence(&input.molecule, &row.selector)? != row.original_sequence)
        {
            return Err(bad("original reference contribution changed"));
        }
    }
    Ok(())
}
dto! {
/// Documentary source shape from an actual native inspection. This DTO is not
/// a duplex, has no nucleotide calls, and cannot restore physical material.
#[derive(PartialEq,Eq)]pub struct MaterialDomain{pub top_length:usize,pub bottom_start:i64,pub bottom_length:usize,pub topology:Topology,pub state_revision:String}}
dto! {#[derive(PartialEq,Eq)]pub struct DigestGeometry{pub top_start:usize,pub top_length:usize,pub bottom_start:i64,pub bottom_length:usize}}
pub fn digest_material_selector(
    source: &MaterialDomain,
    fragment: &DigestGeometry,
    reverse: bool,
) -> Result<Selector, Error> {
    let max = crate::model::MAX_SEQUENCE_LEN;
    if source.top_length == 0
        || source.bottom_length == 0
        || source.top_length > max
        || source.bottom_length > max
        || source.bottom_start.unsigned_abs() > 2 * max as u64
        || !hex(&source.state_revision)
    {
        return Err(bad("invalid documentary source strand domain"));
    }
    if source.topology == Topology::Circular
        && (source.bottom_start != 0 || source.bottom_length != source.top_length)
    {
        return Err(bad(
            "circular source must have actual aligned equal-length strands",
        ));
    }
    let (strand, start, end) = if !reverse {
        let start = if source.topology == Topology::Circular {
            fragment.top_start % source.top_length
        } else {
            fragment.top_start
        };
        let end = start
            .checked_add(fragment.top_length)
            .ok_or(Error::TooLarge)?;
        range(start, end, source.top_length, source.topology)?;
        (Strand::Top, start, end)
    } else {
        if fragment.bottom_length == 0 || fragment.bottom_length > source.bottom_length {
            return Err(bad(
                "selected bottom fragment exceeds actual source bottom strand",
            ));
        }
        let end = fragment
            .bottom_start
            .checked_add(i64::try_from(fragment.bottom_length).map_err(|_| Error::TooLarge)?)
            .ok_or(Error::TooLarge)?;
        let source_end = source
            .bottom_start
            .checked_add(source.bottom_length as i64)
            .ok_or(Error::TooLarge)?;
        let start = if source.topology == Topology::Circular {
            usize::try_from((source_end - end).rem_euclid(source.bottom_length as i64))
                .map_err(|_| Error::TooLarge)?
        } else {
            if fragment.bottom_start < source.bottom_start || end > source_end {
                return Err(bad(
                    "fragment bottom interval lies outside the actual source strand",
                ));
            }
            usize::try_from(source_end - end).map_err(|_| Error::TooLarge)?
        };
        let end = start
            .checked_add(fragment.bottom_length)
            .ok_or(Error::TooLarge)?;
        range(start, end, source.bottom_length, source.topology)?;
        (Strand::Bottom, start, end)
    };
    Ok(Selector::MaterialStrandRange {
        strand,
        start,
        end,
        expected_state_revision: source.state_revision.clone(),
    })
}
