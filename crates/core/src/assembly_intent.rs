//! Explicit intended reference construction, independent of all preparation and
//! assembly simulators. A reference contribution is not physical DNA ancestry.
use crate::{
    Alphabet, Error, Location, Molecule, Topology,
    annotations::{self, MappingSegment},
    assembly_recipe,
    overlap_primer_design::Fragment,
    verification::PlanInput,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ENGINE_VERSION: &str = "assembly-intent-1";
const MAX_ROWS: usize = 32;
const MAX_BASES: usize = 8 << 20;
const MAX_FEATURES: usize = 4096;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request {pub name:String,pub fragments:Vec<Fragment>,pub circular:bool}}
dto! {pub struct PlacedFragment {pub role:String,pub input_id:String,pub input_revision:String,pub start:usize,pub end:usize,pub reverse:bool,pub output_start:usize,pub output_end:usize}}
dto! {pub struct Report {pub context_revision:String,pub request:Request,pub target:Molecule,pub fragments:Vec<PlacedFragment>,pub warnings:Vec<String>}}
struct AnnotationBudget(usize);
impl Write for AnnotationBudget {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= 30 * 1024 * 1024)
            .ok_or_else(|| std::io::Error::other("intended annotation output budget"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn abbreviated(text: &str) -> String {
    if text.len() <= 4096 {
        return text.into();
    }
    let mut end = 4000;
    while !text.is_char_boundary(end) {
        end -= 1;
    }
    format!(
        "{}… (abbreviated; full original input evidence retained)",
        &text[..end]
    )
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("assembly intent: {s}"))
}
fn token(s: &str, cap: usize) -> Result<(), Error> {
    if s.trim().is_empty() || s.len() > cap || s.chars().any(char::is_control) {
        Err(bad("invalid explicit bounded name/source/role"))
    } else {
        Ok(())
    }
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn request_shape(request: &Request) -> Result<(), Error> {
    token(&request.name, 1024)?;
    if request.fragments.is_empty() || request.fragments.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    let mut roles = BTreeSet::new();
    let mut output = 0usize;
    for f in &request.fragments {
        token(&f.role, 128)?;
        token(&f.input_id, 128)?;
        if !roles.insert(f.role.as_str()) {
            return Err(bad(
                "distinct stable roles are required even when a source is reused",
            ));
        }
        let length = f
            .end
            .checked_sub(f.start)
            .filter(|n| *n > 0)
            .ok_or_else(|| bad("interval must be nonempty and end must exceed start"))?;
        output = output
            .checked_add(length)
            .filter(|n| *n <= MAX_BASES)
            .ok_or(Error::TooLarge)?;
    }
    Ok(())
}
/// Borrowed original-document guards. Physical markers may restrict preparation,
/// but do not prevent reading a valid reference for explicit design intent.
pub fn preflight(inputs: &[(&str, &Molecule)], request: &Request) -> Result<(), Error> {
    request_shape(request)?;
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    let models = inputs.iter().map(|(_, m)| *m).collect::<Vec<_>>();
    crate::model::bounds::model_preflight(&models)?;
    let mut lookup = BTreeMap::new();
    let mut bases = 0usize;
    for &(id, m) in inputs {
        token(id, 128)?;
        if lookup.insert(id, m).is_some() {
            return Err(bad("duplicate original input ID"));
        }
        if m.sequence.alphabet() != Alphabet::Dna {
            return Err(bad(
                "intended cloning reference requires DNA; RNA/protein conversion is not inferred",
            ));
        }
        m.validate()?;
        bases = bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= MAX_BASES)
            .ok_or(Error::TooLarge)?;
    }
    for f in &request.fragments {
        let m = lookup
            .get(f.input_id.as_str())
            .copied()
            .ok_or_else(|| bad("requested contribution has no bound original input"))?;
        let n = m.sequence.len();
        let length = f.end - f.start;
        if f.start >= n
            || length > n
            || f.end > f.start.checked_add(n).ok_or(Error::TooLarge)?
            || m.topology == Topology::Linear && f.end > n
        {
            return Err(bad(
                "invalid interval: one unrolled circular turn at most, with linear ends inside the actual source",
            ));
        }
    }
    assembly_recipe::hash(&(inputs, request))?;
    Ok(())
}
pub fn intent_revision(inputs: &[PlanInput], request: &Request) -> Result<String, Error> {
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        request,
    )?;
    for i in inputs {
        if !hex(&i.revision) || assembly_recipe::hash(&i.molecule)? != i.revision {
            return Err(bad("original source revision changed"));
        }
    }
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        inputs
            .iter()
            .map(|i| (&i.id, &i.revision))
            .collect::<Vec<_>>(),
        request,
    ))
}
fn parts(f: &Fragment, n: usize) -> Vec<(usize, usize)> {
    let first_end = f.end.min(n);
    let mut parts = vec![(f.start, first_end)];
    if f.end > n {
        parts.push((0, f.end - n));
    }
    parts
}
fn location(f: &Fragment, n: usize) -> Location {
    let p = parts(f, n)
        .into_iter()
        .map(|(a, b)| Location::range(a, b))
        .collect::<Vec<_>>();
    let location = if p.len() == 1 {
        p.into_iter().next().unwrap()
    } else {
        Location::Join(p)
    };
    if f.reverse {
        Location::Complement(Box::new(location))
    } else {
        location
    }
}
fn maps(f: &Fragment, n: usize, offset: usize) -> Vec<MappingSegment> {
    let parts = parts(f, n);
    let length = f.end - f.start;
    let mut consumed = 0usize;
    parts
        .into_iter()
        .map(|(a, b)| {
            let count = b - a;
            let destination_start = offset
                + if f.reverse {
                    length - consumed - count
                } else {
                    consumed
                };
            consumed += count;
            MappingSegment {
                source_start: a,
                source_end: b,
                destination_start,
                reverse: f.reverse,
            }
        })
        .collect()
}
pub fn build(inputs: &[PlanInput], request: &Request) -> Result<Report, Error> {
    let context_revision = intent_revision(inputs, request)?;
    let lookup = inputs
        .iter()
        .map(|i| (i.id.as_str(), i))
        .collect::<BTreeMap<_, _>>();
    let length = request.fragments.iter().map(|f| f.end - f.start).sum();
    let mut dna = String::with_capacity(length);
    let mut fragments = Vec::with_capacity(request.fragments.len());
    let mut warnings=vec!["Explicit intended reference only: concatenated selected source regions, order and orientation. This is not a PCR/digest/assembly prediction, physical material, source-atom ancestry or experimental result.".into(),"This reference is frozen separately before primer/preparation/assembly planning. Original source documents and event bindings remain separate; changes to source revisions or selected contributions invalidate the intent context.".into()];
    for f in &request.fragments {
        let input = lookup
            .get(f.input_id.as_str())
            .copied()
            .ok_or_else(|| bad("missing original input"))?;
        let output_start = dna.len();
        dna.push_str(
            &location(f, input.molecule.sequence.len()).extract(&input.molecule.sequence)?,
        );
        let output_end = dna.len();
        fragments.push(PlacedFragment {
            role: f.role.clone(),
            input_id: f.input_id.clone(),
            input_revision: input.revision.clone(),
            start: f.start,
            end: f.end,
            reverse: f.reverse,
            output_start,
            output_end,
        });
    }
    let mut target = Molecule::dna(
        &request.name,
        &dna,
        if request.circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
    )?;
    target
        .metadata
        .insert("assembly_intent:engine".into(), ENGINE_VERSION.into());
    target.metadata.insert(
        "assembly_intent:context_revision".into(),
        context_revision.clone(),
    );
    target
        .metadata
        .insert("assembly:intent:reference_only".into(), "true".into());
    target.metadata.insert(
        "duplex_model".into(),
        "Reference-only assembly design intent; physical material state unavailable".into(),
    );
    target.metadata.insert(
        "assembly_intent:contributions".into(),
        serde_json::to_string(&fragments).map_err(|e| bad(&e.to_string()))?,
    );
    target.metadata.insert("assembly_intent:scope".into(),"Operator-declared intended reference from source-coordinate contributions, independent of all preparation/assembly simulation. No materialization or experimental claim.".into());
    target.metadata.insert(
        "assembly_intent:sequence_blake3".into(),
        blake3::hash(target.sequence.as_str().as_bytes())
            .to_hex()
            .to_string(),
    );
    for (f, placement) in request.fragments.iter().zip(&fragments) {
        let source = &lookup
            .get(f.input_id.as_str())
            .ok_or_else(|| bad("missing annotation original"))?
            .molecule;
        if source.features.len() > MAX_FEATURES.saturating_sub(target.features.len()) {
            warnings.push(format!("Role {}: annotations withheld by the4096-feature intended-reference budget. All {} original features remain in original input {} at revision {}; no annotation coordinates were guessed.",f.role,source.features.len(),f.input_id,placement.input_revision));
            continue;
        }
        match annotations::project_features_with_topology(source,&maps(f,source.sequence.len(),placement.output_start),length,target.topology){Ok(projected)=>{let reviewed=projected.features.iter().filter(|f|f.needs_review).count();let mut budget=crate::model::bounds::ModelBudget::default();if budget.molecule(&target).and_then(|()|budget.features(&projected.features)).is_err() || serde_json::to_writer(&mut AnnotationBudget(0),&(&target,&projected.features)).is_err(){warnings.push(format!("Role {}: projected annotations exceed the aggregate native-feature/qualifier budget. Original input {} at revision {} retains all source evidence; intended DNA is unchanged.",f.role,f.input_id,placement.input_revision));continue;}target.features.extend(projected.features);if reviewed>0||!projected.unmapped_features.is_empty()||projected.dropped_features>0{warnings.push(format!("Role {}: {} projected features need review; {} original annotations are unprojectable and {} lie outside the chosen contribution. Original features/raw evidence remain in input {} at revision {}; correspondence does not establish inherited function or chemistry.",f.role,reviewed,projected.unmapped_features.len(),projected.dropped_features,f.input_id,placement.input_revision));}for warning in projected.warnings.iter().take(8){warnings.push(format!("Role {}: {}",f.role,warning.chars().take(1024).collect::<String>()));}if projected.warnings.len()>8{warnings.push(format!("Role {}: additional annotation findings omitted; review flags and original input evidence remain available.",f.role));}},Err(error)=>warnings.push(format!("Role {}: annotation projection withheld ({error}); original annotations remain in input {} at revision {}, with no fictitious current coordinates.",f.role,f.input_id,placement.input_revision))}
    }
    if dna.bytes().any(|b| !b"ACGT".contains(&b)) {
        warnings.push("IUPAC ambiguity is retained as requested reference notation. This does not provide concrete primer/preparation calls or establish a unique sequence.".into());
    }
    target.validate()?;
    crate::model::bounds::model_preflight(&[&target])?;
    let warnings = warnings.into_iter().map(|s| abbreviated(&s)).collect();
    let report = Report {
        context_revision,
        request: request.clone(),
        target,
        fragments,
        warnings,
    };
    report.validate()?;
    Ok(report)
}
impl Report {
    /// Readonly shape/context checks. A decoded intended reference is not
    /// assembly verification or authority to replace any source document.
    pub fn validate(&self) -> Result<(), Error> {
        request_shape(&self.request)?;
        crate::model::bounds::model_preflight(&[&self.target])?;
        self.target.validate()?;
        if !hex(&self.context_revision)
            || self.fragments.len() != self.request.fragments.len()
            || self.target.sequence.alphabet() != Alphabet::Dna
            || self.target.name != self.request.name
            || self.target.topology
                != if self.request.circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                }
            || self.target.sequence.len() > MAX_BASES
            || self.warnings.len() > MAX_ROWS * 10 + 8
            || self.warnings.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("invalid intended-reference report shape"));
        }
        let mut offset = 0usize;
        for (f, p) in self.request.fragments.iter().zip(&self.fragments) {
            if f.role != p.role
                || f.input_id != p.input_id
                || f.start != p.start
                || f.end != p.end
                || f.reverse != p.reverse
                || !hex(&p.input_revision)
                || p.output_start != offset
                || p.output_end.checked_sub(p.output_start) != f.end.checked_sub(f.start)
            {
                return Err(bad("source contribution mapping changed"));
            }
            offset = p.output_end;
        }
        if offset != self.target.sequence.len() {
            return Err(bad(
                "contribution boundaries do not cover the actual intended reference",
            ));
        }
        if self.target.metadata.get("assembly_intent:context_revision")
            != Some(&self.context_revision)
            || self
                .target
                .metadata
                .get("assembly_intent:engine")
                .map(String::as_str)
                != Some(ENGINE_VERSION)
            || self
                .target
                .metadata
                .get("assembly:intent:reference_only")
                .map(String::as_str)
                != Some("true")
            || self.target.metadata.get("duplex_model").map(String::as_str)
                != Some(
                    "Reference-only assembly design intent; physical material state unavailable",
                )
        {
            return Err(bad("intended-reference scope/context guard changed"));
        }
        if self.target.metadata.get("assembly_intent:contributions")
            != Some(&serde_json::to_string(&self.fragments).map_err(|e| bad(&e.to_string()))?)
        {
            return Err(bad("retained source contribution evidence changed"));
        }
        if self
            .target
            .metadata
            .get("assembly_intent:sequence_blake3")
            .map(String::as_str)
            != Some(
                blake3::hash(self.target.sequence.as_str().as_bytes())
                    .to_hex()
                    .as_str(),
            )
        {
            return Err(bad(
                "intended reference calls disagree with their readonly sequence receipt",
            ));
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad(
                "intent belongs to changed sources or a different contribution request",
            ));
        }
        Ok(())
    }
}
