//! Selective proposals preserve every unselected preparation literally. Intended
//! DNA stays source/edit-declared, independent of the later complete recipe check.
use crate::{
    Error, Molecule, Sequence, assembly_intent,
    assembly_junction_edit::{self as junction, Edit, Intent},
    assembly_recipe::{self, FragmentSpec, Preparation},
    overlap_primer_design,
    overlap_profiles::ProfileOptions,
    primers::{self, DesignOptions, DesignReport},
    verification::PlanInput,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const ENGINE_VERSION: &str = "selective-junction-edit-1";
const MAX_ROWS: usize = 32;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub enum Row {Retain(FragmentSpec),DesignPcr{role:String,input_id:String,reverse:bool}}}
dto! {pub struct Request {pub base:assembly_intent::Request,pub edits:Vec<Edit>,pub rows:Vec<Row>,pub overlaps:Vec<usize>,pub primer_options:DesignOptions,pub profile_options:ProfileOptions}}
dto! {pub struct Report {pub context_revision:String,pub request:Request,pub intent:Intent,pub recipe:assembly_recipe::Recipe,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))] pub designs:Vec<Option<DesignReport>>,pub retained_rows:Vec<usize>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("selective junction edit: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn raw_preparation(f: &FragmentSpec) -> Result<(), Error> {
    if f.role.len() > 128 || f.input_id.len() > 128 {
        return Err(Error::TooLarge);
    }
    let primer = |p: &crate::analysis::Primer| -> Result<(), Error> {
        if p.name.len() > 256 || p.annealing.len() > 1000 || p.tail.len() > 1000 {
            return Err(Error::TooLarge);
        }
        Ok(())
    };
    match &f.preparation {
        Preparation::ExactPcr { forward, reverse } => {
            primer(forward)?;
            primer(reverse)?;
        }
        Preparation::PcrCatalogDigest {
            forward,
            reverse,
            selection,
            ..
        } => {
            primer(forward)?;
            primer(reverse)?;
            selection.validate()?;
        }
        Preparation::CatalogDigest { selection, .. } => {
            selection.validate()?;
        }
        Preparation::Digest { enzymes, .. } => {
            if enzymes.len() > 32 || enzymes.iter().any(|s| s.len() > 128) {
                return Err(Error::TooLarge);
            }
        }
        _ => {}
    }
    Ok(())
}
fn shape(r: &Request) -> Result<(), Error> {
    let n = r.base.fragments.len();
    if n == 0
        || n > MAX_ROWS
        || r.rows.len() != n
        || r.overlaps.len() != n - 1 + usize::from(r.base.circular)
        || r.overlaps.iter().any(|h| *h == 0 || *h > 10_000)
    {
        return Err(bad(
            "rows/order and explicit positive junction homologies must match the base intent",
        ));
    }
    r.profile_options.validate()?;
    let o = &r.primer_options;
    if [
        o.target_tm_c,
        o.max_tm_difference_c,
        o.min_gc_fraction,
        o.max_gc_fraction,
        o.conditions.sodium_molar,
        o.conditions.total_strand_molar,
    ]
    .iter()
    .any(|v| !v.is_finite())
    {
        return Err(bad(
            "nonfinite primer configuration cannot form a canonical request",
        ));
    }

    let total = r
        .overlaps
        .iter()
        .try_fold(0usize, |s, h| s.checked_add(*h).ok_or(Error::TooLarge))?;
    if total > 1_000_000 || total > r.profile_options.max_comparisons {
        return Err(Error::TooLarge);
    }
    for (row, f) in r.rows.iter().zip(&r.base.fragments) {
        let (role, input, reverse) = match row {
            Row::Retain(spec) => {
                raw_preparation(spec)?;
                (&spec.role, &spec.input_id, spec.reverse)
            }
            Row::DesignPcr {
                role,
                input_id,
                reverse,
            } => (role, input_id, *reverse),
        };
        if role != &f.role || input != &f.input_id || reverse != f.reverse {
            return Err(bad(
                "row order, role, original input or orientation differs from base intent",
            ));
        }
    }
    if r.primer_options.forward_tail.len() > 1000 || r.primer_options.reverse_tail.len() > 1000 {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn lens(
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
                .ok_or_else(|| bad("unknown original source"))
        })
        .collect()
}
fn options(r: &Request, adjusted: &assembly_intent::Request) -> Result<(), Error> {
    if r.rows.iter().any(|r| matches!(r, Row::DesignPcr { .. })) {
        if !r.primer_options.forward_tail.is_empty() || !r.primer_options.reverse_tail.is_empty() {
            return Err(bad(
                "selected PCR designs require empty option-tail slots; declared user bases cannot be silently replaced",
            ));
        }
        overlap_primer_design::options(&overlap_primer_design::Request {
            name: adjusted.name.clone(),
            fragments: adjusted.fragments.clone(),
            overlap: 1,
            circular: adjusted.circular,
            primer_options: r.primer_options.clone(),
            profile_options: r.profile_options.clone(),
        })?;
    }
    Ok(())
}
/// Guard originals and retained fields before copying; only selected rows require
/// a usable PCR source, concrete core/domain geometry or primer design work.
pub fn preflight(inputs: &[(&str, &Molecule)], r: &Request) -> Result<(), Error> {
    shape(r)?;
    junction::intent_preflight(inputs, &r.base, &r.edits)?;
    let adjusted = junction::adjustment(&r.base, &r.edits, &lens(inputs, &r.base)?)?;
    let insert = junction::inserts(&r.base, &r.edits)?;
    options(r, &adjusted)?;
    let retained = r
        .rows
        .iter()
        .filter_map(|row| {
            if let Row::Retain(f) = row {
                Some(f.clone())
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    if !retained.is_empty() {
        assembly_recipe::preparation_rows_preflight(inputs, &r.base.name, &retained)?;
    }
    let mut work = 0usize;
    let mut scans = 0usize;
    for (i, (row, f)) in r.rows.iter().zip(&adjusted.fragments).enumerate() {
        if !matches!(row, Row::DesignPcr { .. }) {
            continue;
        }
        let source = inputs
            .iter()
            .find(|(id, _)| *id == f.input_id)
            .map(|(_, m)| *m)
            .ok_or_else(|| bad("missing selected PCR original"))?;
        overlap_primer_design::pcr_source(source)?;
        let core = f.end - f.start;
        if core < r.primer_options.min_length * 2 {
            return Err(bad(
                "selected core must fit both nonoverlapping anchored annealing domains",
            ));
        }
        for j in 0..core {
            if !b"ACGT".contains(
                &source.sequence.as_str().as_bytes()[(f.start + j) % source.sequence.len()],
            ) {
                return Err(bad("selected PCR core contains ambiguous calls"));
            }
        }
        let tail = if i < r.overlaps.len() {
            let next = &adjusted.fragments[(i + 1) % adjusted.fragments.len()];
            if r.overlaps[i] > next.end - next.start {
                return Err(bad(
                    "selected overlap exceeds the next adjusted original core",
                ));
            }
            if r.primer_options
                .max_length
                .checked_add(insert[i].len())
                .and_then(|n| n.checked_add(r.overlaps[i]))
                .is_none_or(|n| n > primers::MAX_OLIGO)
            {
                return Err(bad("complete selected oligo exceeds200 bases"));
            }
            let next_source = inputs
                .iter()
                .find(|(id, _)| *id == next.input_id)
                .map(|(_, m)| *m)
                .ok_or_else(|| bad("missing downstream original"))?;
            if !junction::prefix(next_source, next, r.overlaps[i])?
                .bytes()
                .chain(insert[i].bytes())
                .all(|b| b"ACGT".contains(&b))
            {
                return Err(bad(
                    "selected insert and homology must provide concrete DNA",
                ));
            }
            insert[i].len() + r.overlaps[i]
        } else {
            0
        };
        let a = r.primer_options.max_length;
        let b = a.checked_add(tail).ok_or(Error::TooLarge)?;
        if b > primers::MAX_OLIGO {
            return Err(bad(
                "complete selected oligo including insert and homology exceeds200 bases",
            ));
        }
        let count = r.primer_options.max_length - r.primer_options.min_length + 1;
        let pairs = count.checked_mul(count).ok_or(Error::TooLarge)?;
        if pairs
            .checked_mul(a)
            .and_then(|n| n.checked_mul(b))
            .and_then(|n| n.checked_mul(2))
            .is_none_or(|n| n > 20_000_000)
        {
            return Err(bad(
                "selected pair workload exceeds its native profile; narrow lengths",
            ));
        }
        let pairwork = pairs
            .checked_mul(a + b)
            .and_then(|n| n.checked_mul(a.max(b)))
            .ok_or(Error::TooLarge)?;
        let metrics = count
            .checked_mul(a * a + b * b)
            .and_then(|n| n.checked_mul(8))
            .ok_or(Error::TooLarge)?;
        work = work
            .checked_add(pairwork)
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
fn pins(inputs: &[PlanInput]) -> Result<Vec<(&str, &str)>, Error> {
    if inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    for i in inputs {
        if !hex(&i.revision) || i.revision != i.molecule.revision() {
            return Err(bad("original source revision changed"));
        }
    }
    Ok(inputs
        .iter()
        .map(|i| (i.id.as_str(), i.revision.as_str()))
        .collect())
}
pub fn design_revision(inputs: &[PlanInput], r: &Request) -> Result<String, Error> {
    if inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        r,
    )?;
    assembly_recipe::hash(&(ENGINE_VERSION, pins(inputs)?, r))
}
pub fn design(inputs: &[PlanInput], r: &Request) -> Result<Report, Error> {
    let context_revision = design_revision(inputs, r)?;
    let intent = junction::build_intent(inputs, &r.base, &r.edits)?;
    let inserts = junction::inserts(&r.base, &r.edits)?;
    let mut fragments = Vec::with_capacity(r.rows.len());
    let mut designs = Vec::with_capacity(r.rows.len());
    let mut retained_rows = Vec::new();
    for (i, (row, f)) in r
        .rows
        .iter()
        .zip(&intent.adjusted_request.fragments)
        .enumerate()
    {
        match row {
            Row::Retain(spec) => {
                fragments.push(spec.clone());
                designs.push(None);
                retained_rows.push(i);
            }
            Row::DesignPcr { .. } => {
                let source = &inputs
                    .iter()
                    .find(|s| s.id == f.input_id)
                    .ok_or_else(|| bad("missing selected source"))?
                    .molecule;
                let mut options = r.primer_options.clone();
                options.target_start = f.start;
                options.target_end = f.end;
                options.max_pairs = 1;
                if i < r.overlaps.len() {
                    let next = &intent.adjusted_request.fragments[(i + 1) % intent.fragments.len()];
                    let next_source = &inputs
                        .iter()
                        .find(|s| s.id == next.input_id)
                        .ok_or_else(|| bad("missing downstream source"))?
                        .molecule;
                    let added = format!(
                        "{}{}",
                        inserts[i],
                        junction::prefix(next_source, next, r.overlaps[i])?
                    );
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
                let mut report = primers::design_primers(source, &options)?;
                let pair=report.pairs.first_mut().ok_or_else(||bad("no selected anchored pair passes the full-tail filters; retained rows were not changed"))?;
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
                designs.push(Some(report));
            }
        }
    }
    let recipe = assembly_recipe::Recipe {
        name: r.base.name.clone(),
        fragments,
        overlaps: r.overlaps.clone(),
        circular: r.base.circular,
        options: r.profile_options.clone(),
    };
    assembly_recipe::source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        &intent.target,
        &recipe,
    )?;
    let result=Report{context_revision:context_revision.clone(),request:r.clone(),intent,recipe,designs,retained_rows,warnings:vec!["Only explicitly selected PCR rows were regenerated. Every retained preparation, oligo name/tail, enzyme pin and repair choice is unchanged.".into(),"Retained preparations can be incompatible with the newly declared intent or overlaps. Freshly verify the complete recipe against this independent intent; a selective proposal establishes no positive preparation, assembly, mixture-specificity or experimental verdict.".into()]};
    result.validate_binding(&context_revision)?;
    Ok(result)
}
fn warnings(values: &[String]) -> Result<(), Error> {
    if values.len() > 1024 || values.iter().any(|s| s.len() > 8192) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
impl Report {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        shape(&self.request)?;
        if !hex(expected)
            || self.context_revision != expected
            || self.designs.len() != self.request.rows.len()
            || self.recipe.fragments.len() != self.request.rows.len()
        {
            return Err(bad("selective report context or row shape changed"));
        }
        self.intent
            .validate_binding(&self.intent.context_revision)?;
        if assembly_recipe::hash(&self.request.base)?
            != assembly_recipe::hash(&self.intent.base_request)?
            || self.request.edits != self.intent.edits
            || self.recipe.name != self.request.base.name
            || self.recipe.circular != self.request.base.circular
            || self.recipe.overlaps != self.request.overlaps
            || self.recipe.options != self.request.profile_options
        {
            return Err(bad("independent intent or recipe options changed"));
        }
        let raw = self
            .intent
            .target
            .metadata
            .get("assembly_junction_edit:inputs")
            .filter(|s| s.len() < 16_384)
            .ok_or_else(|| bad("missing original source pins"))?;
        let pins: Vec<(String, String)> =
            serde_json::from_str(raw).map_err(|_| bad("invalid source pins"))?;
        if assembly_recipe::hash(&(ENGINE_VERSION, &pins, &self.request))? != self.context_revision
        {
            return Err(bad("request/source revision binding changed"));
        }
        options(&self.request, &self.intent.adjusted_request)?;
        let inserts = junction::inserts(&self.request.base, &self.request.edits)?;
        let mut retained = Vec::new();
        for (i, ((request, row), design)) in self
            .request
            .rows
            .iter()
            .zip(&self.recipe.fragments)
            .zip(&self.designs)
            .enumerate()
        {
            match request {
                Row::Retain(original) => {
                    if design.is_some()
                        || assembly_recipe::hash(original)? != assembly_recipe::hash(row)?
                    {
                        return Err(bad("retained preparation changed"));
                    }
                    retained.push(i);
                }
                Row::DesignPcr { .. } => {
                    let f = &self.intent.adjusted_request.fragments[i];
                    let d = design
                        .as_ref()
                        .ok_or_else(|| bad("selected design evidence missing"))?;
                    if d.pairs.len() != 1
                        || d.source_revision != self.intent.fragments[i].input_revision
                        || d.options.max_pairs != 1
                        || d.options.target_start != f.start
                        || d.options.target_end != f.end
                        || row.role != f.role
                        || row.input_id != f.input_id
                        || row.reverse != f.reverse
                    {
                        return Err(bad("selected design no longer binds its source/interval"));
                    }
                    let Preparation::ExactPcr { forward, reverse } = &row.preparation else {
                        return Err(bad("selected row must remain original-source exact PCR"));
                    };
                    let pair = &d.pairs[0];
                    let mut options = self.request.primer_options.clone();
                    options.target_start = f.start;
                    options.target_end = f.end;
                    options.max_pairs = 1;
                    options.forward_tail = forward.tail.clone();
                    options.reverse_tail = reverse.tail.clone();
                    if options != d.options
                        || assembly_recipe::hash(forward)?
                            != assembly_recipe::hash(&pair.forward.primer)?
                        || assembly_recipe::hash(reverse)?
                            != assembly_recipe::hash(&pair.reverse.primer)?
                    {
                        return Err(bad(
                            "selected recipe primers differ from actual ranked pair",
                        ));
                    }
                    if i < self.request.overlaps.len()
                        && self
                            .request
                            .primer_options
                            .max_length
                            .checked_add(inserts[i].len())
                            .and_then(|n| n.checked_add(self.request.overlaps[i]))
                            .is_none_or(|n| n > primers::MAX_OLIGO)
                    {
                        return Err(bad("selected tail exceeds its bounded ordered profile"));
                    }
                    let added = if i < self.request.overlaps.len() {
                        let ni = (i + 1) % self.intent.fragments.len();
                        let next = &self.intent.fragments[ni];
                        let end = next
                            .output_start
                            .checked_add(self.request.overlaps[i])
                            .filter(|e| *e <= next.output_end - inserts[ni].len())
                            .ok_or_else(|| bad("overlap exceeds downstream original core"))?;
                        format!(
                            "{}{}",
                            inserts[i],
                            &self.intent.target.sequence.as_str()[next.output_start..end]
                        )
                    } else {
                        String::new()
                    };
                    let tail = if added.is_empty() {
                        String::new()
                    } else {
                        Sequence::dna(&added)?
                            .reverse_complement()?
                            .as_str()
                            .to_owned()
                    };
                    let (actual, empty) = if f.reverse {
                        (&forward.tail, &reverse.tail)
                    } else {
                        (&reverse.tail, &forward.tail)
                    };
                    if actual != &tail
                        || !empty.is_empty()
                        || pair.product_length
                            != f.end - f.start + forward.tail.len() + reverse.tail.len()
                    {
                        return Err(bad("selected actual tail geometry changed"));
                    }
                    for (candidate, p) in [(&pair.forward, forward), (&pair.reverse, reverse)] {
                        junction::location_guard(&candidate.location)?;
                        if p.name.len() > 256
                            || p.annealing.len() + p.tail.len() > 200
                            || p.annealing.alphabet() != crate::Alphabet::Dna
                            || !p
                                .annealing
                                .as_str()
                                .bytes()
                                .chain(p.tail.bytes())
                                .all(|b| b"ACGT".contains(&b))
                            || candidate.analysis.name != p.name
                            || candidate.analysis.annealing != p.annealing.as_str()
                            || candidate.analysis.tail != p.tail
                            || candidate.analysis.ordered_sequence != p.ordered_sequence()
                        {
                            return Err(bad(
                                "selected metrics no longer describe the full returned oligo",
                            ));
                        }
                    }
                    warnings(&d.warnings)?;
                }
            }
        }
        if self.retained_rows != retained
            || retained.iter().collect::<BTreeSet<_>>().len() != retained.len()
        {
            return Err(bad("retained row index receipt changed"));
        }
        warnings(&self.warnings)?;
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
