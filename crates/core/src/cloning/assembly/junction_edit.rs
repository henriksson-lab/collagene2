//! Coherent source-bound junction edit intent and explicit PCR oligo proposals.
//! Intended DNA is declared before, and independently of, molecular reconstruction.
use crate::{
    cloning::{
        assembly::{
            intent::{self as assembly_intent, PlacedFragment},
            recipe::{self as assembly_recipe, FragmentSpec, Preparation},
        },
        overlap::{
            primer_design::{self as overlap_primer_design, Fragment},
            profiles::ProfileOptions,
        },
        pcr::primers::{self, DesignOptions, DesignReport},
    },
    planning::plan::PlanInput,
};
use collagene_model::{Error, Feature, Location, Molecule, Sequence};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const ENGINE_VERSION: &str = "assembly-junction-edit-4";
const MAX_ROWS: usize = 32;
const REFERENCE_MODEL: &str =
    "Reference-only junction edit design intent; physical material state unavailable";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(PartialEq,Eq)] pub struct Edit {pub junction:usize,pub remove_left:usize,pub remove_right:usize,pub inserted:String}}
dto! {pub struct Intent {pub context_revision:String,pub source_context_revision:String,pub base_request:assembly_intent::Request,pub edits:Vec<Edit>,pub adjusted_fragments:Vec<assembly_intent::Contribution>,pub target:Molecule,pub fragments:Vec<PlacedFragment>,pub warnings:Vec<String>}}
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
pub(crate) fn inserts(
    base: &assembly_intent::Request,
    edits: &[Edit],
) -> Result<Vec<String>, Error> {
    insertions(base, edits, false)
}
/// Reference intent can declare the final outgoing suffix of linear DNA. The
/// primer/assembly design paths retain their narrower internal-junction scope.
pub(crate) fn intent_inserts(
    base: &assembly_intent::Request,
    edits: &[Edit],
) -> Result<Vec<String>, Error> {
    insertions(base, edits, true)
}
fn insertions(
    base: &assembly_intent::Request,
    edits: &[Edit],
    allow_linear_terminal: bool,
) -> Result<Vec<String>, Error> {
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
    let count = n - 1 + usize::from(base.circular || allow_linear_terminal);
    let mut seen = BTreeSet::new();
    let mut total = 0usize;
    let mut result = vec![String::new(); n];
    for e in edits {
        if e.junction >= count
            || (!base.circular && e.junction == n - 1 && e.remove_right != 0)
            || !seen.insert(e.junction)
            || e.inserted.len() > 400
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
        if text.len() > 200 {
            return Err(bad(
                "each explicit junction insert is limited to200 DNA bases",
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
pub(crate) fn adjustment(
    base: &assembly_intent::Request,
    edits: &[Edit],
    lengths: &[usize],
) -> Result<Vec<assembly_intent::Contribution>, Error> {
    adjusted_sources(base, edits, lengths, false)
}
pub(crate) fn intent_adjustment(
    base: &assembly_intent::Request,
    edits: &[Edit],
    lengths: &[usize],
) -> Result<Vec<assembly_intent::Contribution>, Error> {
    adjusted_sources(base, edits, lengths, true)
}
fn trims(
    base: &assembly_intent::Request,
    edits: &[Edit],
    terminal: bool,
) -> Result<Vec<(usize, usize)>, Error> {
    insertions(base, edits, terminal)?;
    let mut a = vec![(0, 0); base.fragments.len()];
    for e in edits {
        a[e.junction].1 = e.remove_left;
        if e.junction + 1 < a.len() {
            a[e.junction + 1].0 = e.remove_right;
        } else if base.circular {
            a[0].0 = e.remove_right;
        }
    }
    Ok(a)
}
fn adjusted_sources(
    base: &assembly_intent::Request,
    edits: &[Edit],
    lengths: &[usize],
    terminal: bool,
) -> Result<Vec<assembly_intent::Contribution>, Error> {
    let trim = trims(base, edits, terminal)?;
    if lengths.len() != base.fragments.len() {
        return Err(bad("source domain count differs"));
    }
    base.fragments
        .iter()
        .zip(lengths)
        .zip(trim)
        .map(|((f, n), (a, b))| {
            let selector = f.selector.trim(a, b, *n)?;
            Ok(assembly_intent::Contribution {
                role: f.role.clone(),
                input_id: f.input_id.clone(),
                selector,
            })
        })
        .collect()
}
fn lengths(
    inputs: &[(&str, &Molecule)],
    base: &assembly_intent::Request,
) -> Result<Vec<usize>, Error> {
    base.fragments
        .iter()
        .map(|f| {
            f.reference_fragment()?;
            inputs
                .iter()
                .find(|(id, _)| *id == f.input_id)
                .map(|(_, m)| m.sequence.len())
                .ok_or_else(|| bad("missing original reference"))
        })
        .collect()
}
fn pins(inputs: &[PlanInput]) -> Result<Vec<(String, String)>, Error> {
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    Ok(inputs
        .iter()
        .map(|i| (i.id.clone(), i.revision.clone()))
        .collect())
}
fn json<T: Serialize>(value: &T) -> Result<String, Error> {
    serde_json::to_string(value).map_err(|_| Error::TooLarge)
}
fn context(source: &str, base: &assembly_intent::Request, edits: &[Edit]) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, source, base, edits))
}
pub fn intent_preflight(
    inputs: &[(&str, &Molecule)],
    base: &assembly_intent::Request,
    edits: &[Edit],
    workspace: Option<&crate::documents::workspace::Workspace>,
) -> Result<(), Error> {
    assembly_intent::preflight(inputs, base, workspace)?;
    let trim = trims(base, edits, true)?;
    for (f, (a, b)) in base.fragments.iter().zip(trim) {
        if a.checked_add(b)
            .is_none_or(|k| k >= f.selector.length().unwrap_or(0))
        {
            return Err(bad(
                "junction trims remove the entire original contribution",
            ));
        }
    }
    Ok(())
}
pub fn intent_revision(
    inputs: &[PlanInput],
    base: &assembly_intent::Request,
    edits: &[Edit],
    workspace: Option<&crate::documents::workspace::Workspace>,
) -> Result<String, Error> {
    intent_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        base,
        edits,
        workspace,
    )?;
    context(
        &assembly_intent::intent_revision(inputs, base, workspace)?,
        base,
        edits,
    )
}
pub(crate) fn retained_layout(
    base: &assembly_intent::Request,
    edits: &[Edit],
) -> Result<Vec<(usize, usize)>, Error> {
    let trim = trims(base, edits, true)?;
    base.fragments
        .iter()
        .zip(trim)
        .map(|(f, (a, b))| {
            let len = f.selector.length()?;
            if a.checked_add(b).is_none_or(|n| n >= len) {
                return Err(bad("trims remove contribution"));
            }
            Ok((a, len - b))
        })
        .collect()
}
pub(crate) fn build_resolved(
    resolved: &assembly_intent::Resolved<'_>,
    base: &assembly_intent::Request,
    edits: &[Edit],
) -> Result<assembly_intent::Built, Error> {
    let context_revision = context(&resolved.context, base, edits)?;
    let layout = retained_layout(base, edits)?;
    let (target, mut fragments, mut warnings) =
        assembly_intent::build_layout(resolved, base, &layout, &resolved.context)?;
    let mut target = target;
    let insert = intent_inserts(base, edits)?;
    for (i, text) in insert.iter().enumerate().rev() {
        if !text.is_empty() {
            target = target.replace(fragments[i].output_end, fragments[i].output_end, text)?;
        }
    }
    let mut offset = 0;
    for (i, f) in fragments.iter_mut().enumerate() {
        let core = f.retained_end - f.retained_start;
        f.output_start = offset;
        let synthetic = offset + core;
        offset = synthetic + insert[i].len();
        f.output_end = offset;
        if !insert[i].is_empty() {
            target.features.push(Feature{kind:"misc_feature".into(),location:Location::range(synthetic,offset),qualifiers:vec![("label".into(),Some(format!("Declared junction {} insert",i+1))),("note".into(),Some("Operator-declared reference bases, not inherited material or synthesis proof.".into()))],needs_review:true});
        }
    }
    let adjusted_fragments = intent_adjustment(
        base,
        edits,
        &resolved
            .rows
            .iter()
            .map(|r| r.source_length)
            .collect::<Vec<_>>(),
    )?;
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
        (
            "assembly_junction_edit:source_context_revision",
            resolved.context.clone(),
        ),
        (
            "assembly_junction_edit:base_request",
            assembly_intent::request_json(base)?,
        ),
        ("assembly_junction_edit:edits", json(&edits)?),
        (
            "assembly_junction_edit:adjusted_fragments",
            json(&adjusted_fragments)?,
        ),
        (
            "assembly_junction_edit:contributions",
            assembly_intent::placements_json(&fragments)?,
        ),
        (
            "assembly_junction_edit:inputs",
            json(&pins(resolved.inputs)?)?,
        ),
        (
            "assembly_junction_edit:sequence_blake3",
            blake3::hash(target.sequence.as_str().as_bytes())
                .to_hex()
                .to_string(),
        ),
    ] {
        target.metadata.insert(key.into(), value);
    }
    warnings.push("Junction edits affect only the independently declared reference. Original selected strand calls, event/state pins and original-local mappings remain unchanged provenance.".into());
    Ok(assembly_intent::Built {
        target,
        fragments,
        warnings,
    })
}
pub fn build_intent(
    inputs: &[PlanInput],
    base: &assembly_intent::Request,
    edits: &[Edit],
    workspace: Option<&crate::documents::workspace::Workspace>,
) -> Result<Intent, Error> {
    intent_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        base,
        edits,
        workspace,
    )?;
    let resolved = assembly_intent::resolve(inputs, base, workspace)?;
    report_resolved(&resolved, base, edits)
}
pub(crate) fn report_resolved(
    resolved: &assembly_intent::Resolved<'_>,
    base: &assembly_intent::Request,
    edits: &[Edit],
) -> Result<Intent, Error> {
    let built = build_resolved(resolved, base, edits)?;
    let report = Intent {
        context_revision: context(&resolved.context, base, edits)?,
        source_context_revision: resolved.context.clone(),
        base_request: base.clone(),
        edits: edits.to_vec(),
        adjusted_fragments: intent_adjustment(
            base,
            edits,
            &resolved
                .rows
                .iter()
                .map(|r| r.source_length)
                .collect::<Vec<_>>(),
        )?,
        target: built.target,
        fragments: built.fragments,
        warnings: built.warnings,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
fn primer_request(request: &Request, adjusted: &[Fragment]) -> overlap_primer_design::Request {
    overlap_primer_design::Request {
        name: request.base.name.clone(),
        fragments: adjusted.to_vec(),
        overlap: request.overlap,
        circular: request.base.circular,
        primer_options: request.primer_options.clone(),
        profile_options: request.profile_options.clone(),
    }
}
pub fn preflight(inputs: &[(&str, &Molecule)], r: &Request) -> Result<(), Error> {
    // An explicit terminal intention is not a modeled primer tail/homology
    // design. Reject it before any sequence-dependent preparation work.
    inserts(&r.base, &r.edits)?;
    intent_preflight(inputs, &r.base, &r.edits, None)?;
    if !r.primer_options.forward_tail.is_empty() || !r.primer_options.reverse_tail.is_empty() {
        return Err(bad(
            "existing primer-option tails must be empty; user bases are not silently discarded",
        ));
    }
    let adjusted = adjustment(&r.base, &r.edits, &lengths(inputs, &r.base)?)?
        .iter()
        .map(|f| f.reference_fragment())
        .collect::<Result<Vec<_>, _>>()?;
    let insert = inserts(&r.base, &r.edits)?;
    overlap_primer_design::options(&primer_request(r, &adjusted))?;
    let count = r.primer_options.max_length - r.primer_options.min_length + 1;
    let a = r.primer_options.max_length;
    let mut work = 0usize;
    let mut scans = 0usize;
    for (i, f) in adjusted.iter().enumerate() {
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
        let outgoing = i + 1 < adjusted.len() || r.base.circular;
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
    crate::encoding::bounded_json::hash(&(inputs, r))?;
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
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, "design", pins(inputs)?, r))
}
pub(crate) fn prefix(m: &Molecule, f: &Fragment, h: usize) -> Result<String, Error> {
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
    let intent = build_intent(inputs, &r.base, &r.edits, None)?;
    let insert = inserts(&r.base, &r.edits)?;
    let mut designs = Vec::new();
    let mut fragments = Vec::new();
    let adjusted = intent
        .adjusted_fragments
        .iter()
        .map(|f| f.reference_fragment())
        .collect::<Result<Vec<_>, _>>()?;
    for (i, f) in adjusted.iter().enumerate() {
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
            let next = &adjusted[(i + 1) % intent.fragments.len()];
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
            generation: None,
            role: f.role.clone(),
            input_id: f.input_id.clone(),
            preparation: Preparation::ExactPcr {
                forward: pair.forward.primer.clone(),
                reverse: pair.reverse.primer.clone(),
                template_strand: crate::chemistry::Strand::Top,
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
        assembly_intent::shape(&self.base_request)?;
        if self.context_revision != expected
            || !hex(&self.source_context_revision)
            || context(
                &self.source_context_revision,
                &self.base_request,
                &self.edits,
            )? != expected
        {
            return Err(bad("junction source/context changed"));
        }
        collagene_model::bounds::model_preflight(&[&self.target])?;
        self.target.validate()?;
        warning_guard(&self.warnings)?;
        assembly_intent::validate_rows(&self.fragments, &self.base_request, &self.target)?;
        let layout = retained_layout(&self.base_request, &self.edits)?;
        let insert = intent_inserts(&self.base_request, &self.edits)?;
        let adjusted = intent_adjustment(
            &self.base_request,
            &self.edits,
            &self
                .fragments
                .iter()
                .map(|f| f.source_length)
                .collect::<Vec<_>>(),
        )?;
        if crate::encoding::bounded_json::hash(&adjusted)?
            != crate::encoding::bounded_json::hash(&self.adjusted_fragments)?
        {
            return Err(bad("adjusted fragment selectors changed"));
        }
        for (i, (f, (a, b))) in self.fragments.iter().zip(layout).enumerate() {
            if f.retained_start != a
                || f.retained_end != b
                || self
                    .target
                    .sequence
                    .as_str()
                    .get(f.output_start..f.output_end)
                    != Some(
                        format!("{}{}", &f.original_sequence.as_str()[a..b], insert[i]).as_str(),
                    )
            {
                return Err(bad("literal original-strand/junction target changed"));
            }
        }
        for (key, value) in [
            ("duplex_model", REFERENCE_MODEL.into()),
            ("assembly_junction_edit:engine", ENGINE_VERSION.into()),
            (
                "assembly_junction_edit:context_revision",
                self.context_revision.clone(),
            ),
            (
                "assembly_junction_edit:source_context_revision",
                self.source_context_revision.clone(),
            ),
            (
                "assembly_junction_edit:base_request",
                assembly_intent::request_json(&self.base_request)?,
            ),
            ("assembly_junction_edit:edits", json(&self.edits)?),
            (
                "assembly_junction_edit:adjusted_fragments",
                json(&self.adjusted_fragments)?,
            ),
            (
                "assembly_junction_edit:contributions",
                assembly_intent::placements_json(&self.fragments)?,
            ),
            (
                "assembly_junction_edit:sequence_blake3",
                blake3::hash(self.target.sequence.as_str().as_bytes())
                    .to_hex()
                    .to_string(),
            ),
        ] {
            if self.target.metadata.get(key) != Some(&value) {
                return Err(bad("junction intended-reference evidence changed"));
            }
        }
        assembly_intent::validate_material_descriptors(
            &self.target,
            &self.base_request,
            &self.fragments,
        )?;
        assembly_intent::report_bound(self)
    }
}
pub(crate) fn location_guard(l: &Location) -> Result<(), Error> {
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
        if crate::encoding::bounded_json::hash(&(ENGINE_VERSION, "design", &pins, &self.request))?
            != self.context_revision
        {
            return Err(bad("design request/source context changed"));
        }

        if crate::encoding::bounded_json::hash(&self.request.base)?
            != crate::encoding::bounded_json::hash(&self.intent.base_request)?
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
        let adjusted = self
            .intent
            .adjusted_fragments
            .iter()
            .map(|f| f.reference_fragment())
            .collect::<Result<Vec<_>, _>>()?;
        overlap_primer_design::options(&primer_request(&self.request, &adjusted))?;
        let insert = inserts(&self.request.base, &self.request.edits)?;
        for (i, ((f, row), d)) in adjusted
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
            let Preparation::ExactPcr {
                forward,
                reverse,
                template_strand: crate::chemistry::Strand::Top,
                ..
            } = &row.preparation
            else {
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
                || crate::encoding::bounded_json::hash(forward)?
                    != crate::encoding::bounded_json::hash(&pair.forward.primer)?
                || crate::encoding::bounded_json::hash(reverse)?
                    != crate::encoding::bounded_json::hash(&pair.reverse.primer)?
            {
                return Err(bad(
                    "recipe oligos disagree with the actually ranked pair/options",
                ));
            }
            let added = if i + 1 < self.designs.len() || self.request.base.circular {
                let next = &self.intent.fragments[(i + 1) % self.designs.len()];
                let end = next
                    .retained_start
                    .checked_add(self.request.overlap)
                    .filter(|e| *e <= next.retained_end)
                    .ok_or_else(|| bad("homology exceeds downstream original core"))?;
                format!(
                    "{}{}",
                    insert[i],
                    &next.original_sequence.as_str()[next.retained_start..end]
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
        crate::encoding::bounded_json::hash(self)?;
        Ok(())
    }
}
