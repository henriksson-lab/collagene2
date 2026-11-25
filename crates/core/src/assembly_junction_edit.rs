//! Coherent source-bound junction edit intent and explicit PCR oligo proposals.
//! Intended DNA is declared before, and independently of, molecular reconstruction.
use crate::{
    Error, Feature, Location, Molecule, Sequence, Topology,
    assembly_intent::{self, PlacedFragment},
    assembly_recipe::{self, FragmentSpec, Preparation},
    overlap_primer_design::{self, Fragment},
    overlap_profiles::ProfileOptions,
    primers::{self, DesignOptions, DesignReport},
    verification::PlanInput,
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const ENGINE_VERSION: &str = "assembly-junction-edit-1";
const MAX_ROWS: usize = 32;
const MAX_BASES: usize = 8 << 20;
const REFERENCE_MODEL: &str =
    "Reference-only junction edit design intent; physical material state unavailable";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)] pub struct Edit {pub junction:usize,pub remove_left:usize,pub remove_right:usize,pub inserted:String}}
dto! {pub struct Intent {pub context_revision:String,pub base_request:assembly_intent::Request,pub edits:Vec<Edit>,pub adjusted_request:assembly_intent::Request,pub target:Molecule,pub fragments:Vec<PlacedFragment>,pub warnings:Vec<String>}}
dto! {pub struct Request {pub base:assembly_intent::Request,pub edits:Vec<Edit>,pub overlap:usize,pub primer_options:DesignOptions,pub profile_options:ProfileOptions}}
dto! {pub struct Report {pub context_revision:String,pub request:Request,pub intent:Intent,pub recipe:assembly_recipe::Recipe,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))] pub designs:Vec<DesignReport>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("assembly junction edit: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn inserts(base: &assembly_intent::Request, edits: &[Edit]) -> Result<Vec<String>, Error> {
    let n = base.fragments.len();
    if n == 0
        || n > MAX_ROWS
        || edits.len() > MAX_ROWS
        || base.name.trim().is_empty()
        || base.name.len() > 1024
        || base.name.chars().any(char::is_control)
    {
        return Err(bad(
            "requires1..32 source rows, bounded name and at most32 unique edits",
        ));
    }
    let count = n - 1 + usize::from(base.circular);
    let mut seen = BTreeSet::new();
    let mut total = 0usize;
    let mut result = vec![String::new(); n];
    for e in edits {
        if e.junction >= count
            || !seen.insert(e.junction)
            || e.inserted.len() > 200
            || !e.inserted.is_ascii()
        {
            return Err(bad(
                "invalid or repeated junction, or inserted text exceeds bounded DNA profile",
            ));
        }
        let text = if e.inserted.bytes().all(|b| b.is_ascii_whitespace()) {
            String::new()
        } else {
            Sequence::dna(&e.inserted)?.as_str().to_owned()
        };
        if text.len() > 100 {
            return Err(bad(
                "each explicit junction insert is limited to100 DNA bases",
            ));
        }
        total = total
            .checked_add(text.len())
            .filter(|n| *n <= 4096)
            .ok_or(Error::TooLarge)?;
        result[e.junction] = text;
    }
    Ok(result)
}
fn adjustment(
    base: &assembly_intent::Request,
    edits: &[Edit],
    lengths: &[usize],
) -> Result<assembly_intent::Request, Error> {
    inserts(base, edits)?;
    if lengths.len() != base.fragments.len() {
        return Err(bad("source lengths disagree with rows"));
    }
    let mut left = vec![0usize; lengths.len()];
    let mut right = left.clone();
    for e in edits {
        right[e.junction] = e.remove_left;
        left[(e.junction + 1) % lengths.len()] = e.remove_right;
    }
    let mut seen = BTreeSet::new();
    for (f, n) in base.fragments.iter().zip(lengths) {
        if f.role.is_empty()
            || f.role.len() > 128
            || f.input_id.is_empty()
            || f.input_id.len() > 128
            || !seen.insert(&f.role)
            || f.start >= *n
            || f.end.checked_sub(f.start).is_none_or(|v| v == 0 || v > *n)
            || f.start.checked_add(*n).is_none_or(|end| f.end > end)
        {
            return Err(bad(
                "invalid bounded source row identity or original interval",
            ));
        }
    }
    let mut adjusted = base.clone();
    for (i, f) in adjusted.fragments.iter_mut().enumerate() {
        let length = f
            .end
            .checked_sub(f.start)
            .filter(|v| *v > 0 && *v <= lengths[i])
            .ok_or_else(|| bad("invalid original source interval"))?;
        if left[i].checked_add(right[i]).is_none_or(|n| n >= length) {
            return Err(bad(
                "combined incoming/outgoing trims must retain at least one original core base",
            ));
        }
        if f.reverse {
            f.start = f.start.checked_add(right[i]).ok_or(Error::TooLarge)?;
            f.end -= left[i];
        } else {
            f.start = f.start.checked_add(left[i]).ok_or(Error::TooLarge)?;
            f.end -= right[i];
        }
        if lengths[i] == 0 {
            return Err(bad("empty source"));
        }
        if f.start >= lengths[i] {
            f.start -= lengths[i];
            f.end -= lengths[i];
        }
    }
    Ok(adjusted)
}
struct Count(usize);
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= 32 << 20)
            .ok_or_else(|| std::io::Error::other("annotation byte cap"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn lengths(
    inputs: &[(&str, &Molecule)],
    base: &assembly_intent::Request,
) -> Result<Vec<usize>, Error> {
    base.fragments
        .iter()
        .map(|f| {
            inputs
                .iter()
                .find(|(id, _)| *id == f.input_id)
                .map(|(_, m)| m.sequence.len())
                .ok_or_else(|| bad("source role lacks its explicitly bound input"))
        })
        .collect()
}
pub fn intent_preflight(
    inputs: &[(&str, &Molecule)],
    base: &assembly_intent::Request,
    edits: &[Edit],
) -> Result<(), Error> {
    let insert = inserts(base, edits)?;
    assembly_intent::preflight(inputs, base)?;
    let adjusted = adjustment(base, edits, &lengths(inputs, base)?)?;
    assembly_intent::preflight(inputs, &adjusted)?;
    let target = adjusted
        .fragments
        .iter()
        .try_fold(0usize, |sum, f| {
            sum.checked_add(f.end - f.start).ok_or(Error::TooLarge)
        })?
        .checked_add(insert.iter().map(String::len).sum::<usize>())
        .ok_or(Error::TooLarge)?;
    let original = inputs.iter().try_fold(0usize, |sum, (_, m)| {
        sum.checked_add(m.sequence.len()).ok_or(Error::TooLarge)
    })?;
    if original.checked_add(target).is_none_or(|n| n > MAX_BASES) {
        return Err(Error::TooLarge);
    }
    let mut annotations = Count(0);
    for (_, m) in inputs {
        serde_json::to_writer(&mut annotations, &m.features).map_err(|_| Error::TooLarge)?;
    }
    if annotations
        .0
        .checked_add(target)
        .and_then(|n| n.checked_mul(1 + insert.iter().filter(|s| !s.is_empty()).count()))
        .is_none_or(|n| n > 128 << 20)
    {
        return Err(Error::TooLarge);
    }
    assembly_recipe::hash(&(inputs, base, edits))?;
    Ok(())
}
fn pins(inputs: &[PlanInput]) -> Result<Vec<(&str, &str)>, Error> {
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    for i in inputs {
        if !hex(&i.revision) || i.molecule.revision() != i.revision {
            return Err(bad("original input revision changed"));
        }
    }
    Ok(inputs
        .iter()
        .map(|i| (i.id.as_str(), i.revision.as_str()))
        .collect())
}
pub fn intent_revision(
    inputs: &[PlanInput],
    base: &assembly_intent::Request,
    edits: &[Edit],
) -> Result<String, Error> {
    intent_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        base,
        edits,
    )?;
    assembly_recipe::hash(&(ENGINE_VERSION, "intent", pins(inputs)?, base, edits))
}
fn json(v: &impl Serialize) -> Result<String, Error> {
    assembly_recipe::hash(v)?;
    serde_json::to_string(v).map_err(|e| bad(&e.to_string()))
}
pub fn build_intent(
    inputs: &[PlanInput],
    base: &assembly_intent::Request,
    edits: &[Edit],
) -> Result<Intent, Error> {
    let context_revision = intent_revision(inputs, base, edits)?;
    let lens = lengths(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        base,
    )?;
    let adjusted_request = adjustment(base, edits, &lens)?;
    let insert = inserts(base, edits)?;
    let initial = assembly_intent::build(inputs, &adjusted_request)?;
    let mut target = initial.target;
    let mut warnings = initial.warnings;
    for (i, text) in insert.iter().enumerate().rev() {
        if !text.is_empty() {
            let boundary = initial.fragments[i].output_end;
            target = target.replace(boundary, boundary, text)?;
        }
    }
    let mut fragments = initial.fragments;
    let mut offset = 0usize;
    for (i, f) in fragments.iter_mut().enumerate() {
        f.output_start = offset;
        let core = f.end - f.start;
        let synthetic_start = offset + core;
        offset = synthetic_start + insert[i].len();
        f.output_end = offset;
        if !insert[i].is_empty() {
            target.features.push(Feature{kind:"misc_feature".into(),location:Location::range(synthetic_start,offset),qualifiers:vec![("label".into(),Some(format!("Declared junction {} insert",i+1))),("note".into(),Some("Operator-declared synthetic reference bases, not inherited source coordinates or physical synthesis evidence.".into()))],needs_review:true});
        }
    }
    target
        .metadata
        .retain(|k, _| !k.starts_with("assembly_intent:"));
    for (key, value) in [
        ("duplex_model", REFERENCE_MODEL.into()),
        ("assembly_junction_edit:reference_only", "true".into()),
        ("assembly_junction_edit:engine", ENGINE_VERSION.into()),
        (
            "assembly_junction_edit:context_revision",
            context_revision.clone(),
        ),
        ("assembly_junction_edit:base_request", json(base)?),
        ("assembly_junction_edit:edits", json(&edits)?),
        (
            "assembly_junction_edit:adjusted_request",
            json(&adjusted_request)?,
        ),
        ("assembly_junction_edit:contributions", json(&fragments)?),
        ("assembly_junction_edit:source_lengths", json(&lens)?),
        ("assembly_junction_edit:inputs", json(&pins(inputs)?)?),
        (
            "assembly_junction_edit:sequence_blake3",
            blake3::hash(target.sequence.as_str().as_bytes())
                .to_hex()
                .to_string(),
        ),
    ] {
        target.metadata.insert(key.into(), value);
    }
    warnings.push("Source trims and literal inserted DNA define this intended reference independently of PCR or assembly simulation. Synthetic suffixes belong to the upstream row; original sources remain unchanged.".into());
    warnings.push("A changed intended reference invalidates prior primers, product checks and order evidence. Reconstruct the new explicit preparation recipe before treating it as computationally consistent.".into());
    let result = Intent {
        context_revision: context_revision.clone(),
        base_request: base.clone(),
        edits: edits.to_vec(),
        adjusted_request,
        target,
        fragments,
        warnings,
    };
    result.validate_binding(&context_revision)?;
    Ok(result)
}
fn primer_request(
    request: &Request,
    adjusted: &assembly_intent::Request,
) -> overlap_primer_design::Request {
    overlap_primer_design::Request {
        name: adjusted.name.clone(),
        fragments: adjusted.fragments.clone(),
        overlap: request.overlap,
        circular: adjusted.circular,
        primer_options: request.primer_options.clone(),
        profile_options: request.profile_options.clone(),
    }
}
pub fn preflight(inputs: &[(&str, &Molecule)], r: &Request) -> Result<(), Error> {
    intent_preflight(inputs, &r.base, &r.edits)?;
    if !r.primer_options.forward_tail.is_empty() || !r.primer_options.reverse_tail.is_empty() {
        return Err(bad(
            "existing primer-option tails must be empty; user bases are not silently discarded",
        ));
    }
    let adjusted = adjustment(&r.base, &r.edits, &lengths(inputs, &r.base)?)?;
    let insert = inserts(&r.base, &r.edits)?;
    overlap_primer_design::options(&primer_request(r, &adjusted))?;
    let count = r.primer_options.max_length - r.primer_options.min_length + 1;
    let a = r.primer_options.max_length;
    let mut work = 0usize;
    let mut scans = 0usize;
    for (i, f) in adjusted.fragments.iter().enumerate() {
        let source = inputs
            .iter()
            .find(|(id, _)| *id == f.input_id)
            .map(|(_, m)| *m)
            .ok_or_else(|| bad("missing original"))?;
        overlap_primer_design::pcr_source(source)?;
        let core = f.end - f.start;
        if core <= r.overlap || core < r.primer_options.min_length * 2 {
            return Err(bad(
                "adjusted core must exceed homology and retain both nonoverlapping annealing domains",
            ));
        }
        for j in 0..core {
            if !b"ACGT".contains(
                &source.sequence.as_str().as_bytes()[(f.start + j) % source.sequence.len()],
            ) {
                return Err(bad(
                    "ambiguous core cannot provide concrete oligo proposals",
                ));
            }
        }
        if !insert[i].bytes().all(|b| b"ACGT".contains(&b)) {
            return Err(bad(
                "ambiguous literal insert is retained as intent, but cannot provide concrete oligo proposals",
            ));
        }
        let outgoing = i + 1 < adjusted.fragments.len() || adjusted.circular;
        let b = a + if outgoing {
            insert[i].len() + r.overlap
        } else {
            0
        };
        if b > primers::MAX_OLIGO {
            return Err(bad(
                "complete designed oligo including insert and homology must fit the200-base profile",
            ));
        }
        let pair = count
            .checked_mul(count)
            .and_then(|n| n.checked_mul(a))
            .and_then(|n| n.checked_mul(b))
            .and_then(|n| n.checked_mul(2))
            .ok_or(Error::TooLarge)?;
        if pair > 20_000_000 {
            return Err(bad(
                "narrow primer lengths: per-row pair work exceeds its native limit",
            ));
        }
        let pair_charge = count
            .checked_mul(count)
            .and_then(|n| n.checked_mul(a + b))
            .and_then(|n| n.checked_mul(a.max(b)))
            .ok_or(Error::TooLarge)?;
        let metrics = count
            .checked_mul(a * a + b * b)
            .and_then(|n| n.checked_mul(8))
            .ok_or(Error::TooLarge)?;
        work = work
            .checked_add(pair_charge)
            .and_then(|n| n.checked_add(metrics))
            .filter(|n| *n <= 100_000_000)
            .ok_or(Error::TooLarge)?;
        scans = scans
            .checked_add(
                source
                    .sequence
                    .len()
                    .checked_mul(a)
                    .and_then(|n| n.checked_mul(8))
                    .ok_or(Error::TooLarge)?,
            )
            .filter(|n| *n <= 100_000_000)
            .ok_or(Error::TooLarge)?;
    }
    assembly_recipe::hash(&(inputs, r))?;
    Ok(())
}
pub fn design_revision(inputs: &[PlanInput], r: &Request) -> Result<String, Error> {
    preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        r,
    )?;
    assembly_recipe::hash(&(ENGINE_VERSION, "design", pins(inputs)?, r))
}
fn prefix(m: &Molecule, f: &Fragment, h: usize) -> Result<String, Error> {
    let n = m.sequence.len();
    let start = if f.reverse { f.end - h } else { f.start };
    let mut bases = String::with_capacity(h);
    for j in 0..h {
        bases.push(m.sequence.as_str().as_bytes()[(start + j) % n] as char);
    }
    if f.reverse {
        Ok(Sequence::dna(&bases)?.reverse_complement()?.as_str().into())
    } else {
        Ok(bases)
    }
}
pub fn design(inputs: &[PlanInput], r: &Request) -> Result<Report, Error> {
    let context_revision = design_revision(inputs, r)?;
    let intent = build_intent(inputs, &r.base, &r.edits)?;
    let insert = inserts(&r.base, &r.edits)?;
    let mut designs = Vec::new();
    let mut fragments = Vec::new();
    for (i, f) in intent.adjusted_request.fragments.iter().enumerate() {
        let m = &inputs
            .iter()
            .find(|s| s.id == f.input_id)
            .ok_or_else(|| bad("missing original"))?
            .molecule;
        let mut options = r.primer_options.clone();
        options.target_start = f.start;
        options.target_end = f.end;
        options.max_pairs = 1;
        if i + 1 < intent.fragments.len() || r.base.circular {
            let next = &intent.adjusted_request.fragments[(i + 1) % intent.fragments.len()];
            let next_source = &inputs
                .iter()
                .find(|s| s.id == next.input_id)
                .ok_or_else(|| bad("missing downstream original"))?
                .molecule;
            let added = format!("{}{}", insert[i], prefix(next_source, next, r.overlap)?);
            let tail = Sequence::dna(&added)?
                .reverse_complement()?
                .as_str()
                .to_owned();
            if f.reverse {
                options.forward_tail = tail;
            } else {
                options.reverse_tail = tail;
            }
        }
        let mut d = primers::design_primers(m, &options)?;
        let pair=d.pairs.first_mut().ok_or_else(||bad("no anchored pair passes the declared filters with the complete inserted/homology tails"))?;
        pair.forward.primer.name = format!("{}-forward", f.role);
        pair.reverse.primer.name = format!("{}-reverse", f.role);
        pair.forward.analysis.name = pair.forward.primer.name.clone();
        pair.reverse.analysis.name = pair.reverse.primer.name.clone();
        fragments.push(FragmentSpec {
            role: f.role.clone(),
            input_id: f.input_id.clone(),
            preparation: Preparation::ExactPcr {
                forward: pair.forward.primer.clone(),
                reverse: pair.reverse.primer.clone(),
            },
            reverse: f.reverse,
        });
        designs.push(d);
    }
    let recipe = assembly_recipe::Recipe {
        name: r.base.name.clone(),
        fragments,
        overlaps: vec![r.overlap; intent.fragments.len() - 1 + usize::from(r.base.circular)],
        circular: r.base.circular,
        options: r.profile_options.clone(),
    };
    assembly_recipe::preflight(inputs, &intent.target, &recipe)?;
    let report=Report{context_revision:context_revision.clone(),request:r.clone(),intent,recipe,designs,warnings:vec!["One anchored in-house ranked pair per row, evaluated with the full actual insert/homology tails. This is a proposal, not unique amplification, mixture specificity, laboratory success or a verified construct.".into(),"Freshly reconstruct the returned recipe against this independently declared intent target. Manual or saved oligos are never overwritten without explicit application.".into()]};
    report.validate_binding(&context_revision)?;
    Ok(report)
}
fn warning_guard(w: &[String]) -> Result<(), Error> {
    if w.len() > 1024 || w.iter().any(|s| s.len() > 8192) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
impl Intent {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !hex(expected)
            || self.context_revision != expected
            || self.fragments.len() != self.base_request.fragments.len()
        {
            return Err(bad("intent context or row count changed"));
        }
        let insert = inserts(&self.base_request, &self.edits)?;
        crate::model::bounds::model_preflight(&[&self.target])?;
        self.target.validate()?;
        warning_guard(&self.warnings)?;
        let raw = self
            .target
            .metadata
            .get("assembly_junction_edit:source_lengths")
            .filter(|s| s.len() < 2048)
            .ok_or_else(|| bad("missing source lengths"))?;
        let lengths: Vec<usize> =
            serde_json::from_str(raw).map_err(|_| bad("invalid source lengths"))?;
        if lengths.iter().any(|n| *n == 0 || *n > MAX_BASES) {
            return Err(bad("invalid source bounds"));
        }
        let adjusted = adjustment(&self.base_request, &self.edits, &lengths)?;
        if assembly_recipe::hash(&adjusted)? != assembly_recipe::hash(&self.adjusted_request)? {
            return Err(bad("adjusted source coordinates changed"));
        }
        let mut offset = 0usize;
        for (i, (f, p)) in self
            .adjusted_request
            .fragments
            .iter()
            .zip(&self.fragments)
            .enumerate()
        {
            if p.role != f.role
                || p.input_id != f.input_id
                || p.start != f.start
                || p.end != f.end
                || p.reverse != f.reverse
                || !hex(&p.input_revision)
                || p.output_start != offset
                || p.output_end.checked_sub(offset) != Some(f.end - f.start + insert[i].len())
            {
                return Err(bad("actual source/synthetic contribution mapping changed"));
            }
            offset = p.output_end;
        }
        if offset != self.target.sequence.len()
            || self.target.name != self.base_request.name
            || self.target.topology
                != if self.base_request.circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                }
        {
            return Err(bad("intended reference dimensions changed"));
        }
        let raw = self
            .target
            .metadata
            .get("assembly_junction_edit:inputs")
            .filter(|s| s.len() < 16_384)
            .ok_or_else(|| bad("missing source revision pins"))?;
        let pins: Vec<(String, String)> =
            serde_json::from_str(raw).map_err(|_| bad("invalid source pins"))?;
        if pins.is_empty()
            || pins.len() > MAX_ROWS
            || pins
                .iter()
                .any(|(id, rev)| id.is_empty() || id.len() > 128 || !hex(rev))
            || pins.iter().map(|(id, _)| id).collect::<BTreeSet<_>>().len() != pins.len()
            || self.fragments.iter().any(|f| {
                !pins
                    .iter()
                    .any(|(id, rev)| id == &f.input_id && rev == &f.input_revision)
            })
            || assembly_recipe::hash(&(
                ENGINE_VERSION,
                "intent",
                &pins,
                &self.base_request,
                &self.edits,
            ))? != self.context_revision
        {
            return Err(bad("intent source/request revision pins changed"));
        }
        for (key, value) in [
            ("duplex_model", REFERENCE_MODEL.into()),
            ("assembly_junction_edit:reference_only", "true".into()),
            ("assembly_junction_edit:engine", ENGINE_VERSION.into()),
            (
                "assembly_junction_edit:context_revision",
                self.context_revision.clone(),
            ),
            (
                "assembly_junction_edit:base_request",
                json(&self.base_request)?,
            ),
            ("assembly_junction_edit:edits", json(&self.edits)?),
            (
                "assembly_junction_edit:adjusted_request",
                json(&self.adjusted_request)?,
            ),
            (
                "assembly_junction_edit:contributions",
                json(&self.fragments)?,
            ),
            (
                "assembly_junction_edit:sequence_blake3",
                blake3::hash(self.target.sequence.as_str().as_bytes())
                    .to_hex()
                    .to_string(),
            ),
        ] {
            if self.target.metadata.get(key) != Some(&value) {
                return Err(bad("intended target evidence changed"));
            }
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
fn location_guard(l: &Location) -> Result<(), Error> {
    let mut stack = vec![(l, 0usize)];
    let mut count = 0usize;
    while let Some((l, depth)) = stack.pop() {
        count += 1;
        if depth > 32 || count > 1024 {
            return Err(Error::TooLarge);
        }
        match l {
            Location::Complement(inner) => stack.push((inner, depth + 1)),
            Location::Join(parts) | Location::Order(parts) => {
                if parts.len() > 1024 {
                    return Err(Error::TooLarge);
                }
                for p in parts {
                    stack.push((p, depth + 1));
                }
            }
            _ => {}
        }
    }
    Ok(())
}
impl Report {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !hex(expected)
            || self.context_revision != expected
            || self.designs.len() != self.intent.fragments.len()
            || self.recipe.fragments.len() != self.designs.len()
        {
            return Err(bad("proposal context or row count changed"));
        }
        self.intent
            .validate_binding(&self.intent.context_revision)?;
        let raw = self
            .intent
            .target
            .metadata
            .get("assembly_junction_edit:inputs")
            .filter(|s| s.len() < 16_384)
            .ok_or_else(|| bad("missing design source pins"))?;
        let pins: Vec<(String, String)> =
            serde_json::from_str(raw).map_err(|_| bad("invalid design source pins"))?;
        if assembly_recipe::hash(&(ENGINE_VERSION, "design", &pins, &self.request))?
            != self.context_revision
        {
            return Err(bad("design request/source context changed"));
        }

        if assembly_recipe::hash(&self.request.base)?
            != assembly_recipe::hash(&self.intent.base_request)?
            || self.request.edits != self.intent.edits
            || self.recipe.name != self.request.base.name
            || self.recipe.circular != self.request.base.circular
            || self.recipe.overlaps
                != vec![
                    self.request.overlap;
                    self.designs.len() - 1 + usize::from(self.request.base.circular)
                ]
            || self.recipe.options != self.request.profile_options
        {
            return Err(bad("proposal intent/recipe binding changed"));
        }
        if !self.request.primer_options.forward_tail.is_empty()
            || !self.request.primer_options.reverse_tail.is_empty()
        {
            return Err(bad("unexpected original option tails"));
        }
        overlap_primer_design::options(&primer_request(
            &self.request,
            &self.intent.adjusted_request,
        ))?;
        let insert = inserts(&self.request.base, &self.request.edits)?;
        for (i, ((f, row), d)) in self
            .intent
            .adjusted_request
            .fragments
            .iter()
            .zip(&self.recipe.fragments)
            .zip(&self.designs)
            .enumerate()
        {
            if f.role != row.role
                || f.input_id != row.input_id
                || f.reverse != row.reverse
                || d.source_revision != self.intent.fragments[i].input_revision
                || d.pairs.len() != 1
                || d.options.max_pairs != 1
                || d.options.target_start != f.start
                || d.options.target_end != f.end
            {
                return Err(bad("ranked pair source/interval changed"));
            }
            let Preparation::ExactPcr { forward, reverse } = &row.preparation else {
                return Err(bad("junction proposals require exact original-source PCR"));
            };
            let pair = &d.pairs[0];
            let mut options = self.request.primer_options.clone();
            options.target_start = f.start;
            options.target_end = f.end;
            options.max_pairs = 1;
            options.forward_tail = forward.tail.clone();
            options.reverse_tail = reverse.tail.clone();
            if options != d.options
                || assembly_recipe::hash(forward)? != assembly_recipe::hash(&pair.forward.primer)?
                || assembly_recipe::hash(reverse)? != assembly_recipe::hash(&pair.reverse.primer)?
            {
                return Err(bad(
                    "recipe oligos disagree with the actually ranked pair/options",
                ));
            }
            let added = if i + 1 < self.designs.len() || self.request.base.circular {
                let next = &self.intent.fragments[(i + 1) % self.designs.len()];
                let end = next
                    .output_start
                    .checked_add(self.request.overlap)
                    .filter(|e| *e <= next.output_end - insert[(i + 1) % insert.len()].len())
                    .ok_or_else(|| bad("homology exceeds downstream original core"))?;
                format!(
                    "{}{}",
                    insert[i],
                    &self.intent.target.sequence.as_str()[next.output_start..end]
                )
            } else {
                String::new()
            };
            let expected_tail = if added.is_empty() {
                String::new()
            } else {
                Sequence::dna(&added)?.reverse_complement()?.as_str().into()
            };
            let (tail, other) = if f.reverse {
                (&forward.tail, &reverse.tail)
            } else {
                (&reverse.tail, &forward.tail)
            };
            if tail != &expected_tail
                || !other.is_empty()
                || pair.product_length != f.end - f.start + forward.tail.len() + reverse.tail.len()
            {
                return Err(bad("ranked complete-tail geometry changed"));
            }
            for (c, p) in [(&pair.forward, forward), (&pair.reverse, reverse)] {
                location_guard(&c.location)?;
                if p.name.len() > 256
                    || p.annealing.len() + p.tail.len() > primers::MAX_OLIGO
                    || c.analysis.name != p.name
                    || c.analysis.annealing != p.annealing.as_str()
                    || c.analysis.tail != p.tail
                    || c.analysis.ordered_sequence != p.ordered_sequence()
                {
                    return Err(bad(
                        "actual primer metrics do not describe the returned complete oligo",
                    ));
                }
            }
            warning_guard(&d.warnings)?;
        }
        warning_guard(&self.warnings)?;
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
