//! Source-bound signed fusion-placement proposals. Search preserves a frozen
//! independent reference; only explicit materialization invokes PCR/digestion.
use crate::{
    Error, Molecule, Sequence,
    assembly_junction_edit::{self as junction, Intent},
    assembly_recipe::{self, Preparation},
    catalog_digest,
    chemistry::{End, Polarity, Terminus},
    golden_gate_junction_edit::{self as gate, Row},
    golden_gate_primer_design as primers, ligation_recipe,
    verification::PlanInput,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const ENGINE_VERSION: &str = "golden-gate-fusion-design-1";
pub const MAX_ROWS: usize = 32;
pub const MAX_SHIFT: i32 = 128;
pub const MAX_STATES: usize = 100_000;
pub const MAX_CANDIDATES: usize = 256;
pub const MAX_SEARCH_WORK: usize = 20_000_000;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Window {pub junction:usize,pub min_shift:i32,pub max_shift:i32}}
dto! {pub struct Constraints {pub reject_palindromes:bool,pub distinct_reverse_complements:bool,pub minimum_hamming:usize}}
impl Default for Constraints {
    fn default() -> Self {
        Self {
            reject_palindromes: true,
            distinct_reverse_complements: true,
            minimum_hamming: 1,
        }
    }
}
dto! {pub struct Limits {pub max_states:usize,pub max_candidates:usize}}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_states: 10_000,
            max_candidates: 64,
        }
    }
}
dto! {pub struct Request {pub base:gate::Request,pub windows:Vec<Window>,pub constraints:Constraints,pub limits:Limits}}
dto! {pub struct Candidate {pub id:String,pub shifts:Vec<i32>,pub fusions:Vec<String>}}
dto! {pub struct Selection {pub context_revision:String,pub candidate_id:String}}
dto! {pub struct Report {
    pub context_revision:String,pub request:Request,pub original_intent:Intent,
    pub terminal_fusions:Vec<String>,pub candidates:Vec<Candidate>,pub states_examined:usize,
    pub search_complete:bool,pub limits:Vec<String>,pub warnings:Vec<String>,
}}
// Intended joining pairs and literal, unpaired physical competitors. Linear
// terminal keys in Report are assembly-top keys; the right physical 5′ end
// therefore has the reverse-complement orientation.
dto! {#[derive(PartialEq,Eq)] pub struct EndScope {pub intended_junctions:Vec<String>,pub physical_competitors:Vec<String>}}
pub(crate) fn candidate_end_scope(
    report: &Report,
    candidate: &Candidate,
) -> Result<EndScope, Error> {
    Ok(EndScope {
        intended_junctions: candidate.fusions.clone(),
        physical_competitors: if report.request.base.base.circular {
            Vec::new()
        } else {
            if report.terminal_fusions.len() != 2 {
                return Err(bad("linear terminal end count changed"));
            }
            vec![
                report.terminal_fusions[0].clone(),
                rc(&report.terminal_fusions[1])?,
            ]
        },
    })
}
/// Read-only scope extraction, never a fresh scientific receipt or a default
/// candidate selection. Imported reports must still be recomputed for authority.
pub fn end_scope(report: &Report, candidate_id: &str) -> Result<EndScope, Error> {
    report.validate_binding(&report.context_revision)?;
    let candidate = report
        .candidates
        .iter()
        .find(|c| c.id == candidate_id)
        .ok_or_else(|| bad("explicit candidate is not retained"))?;
    candidate_end_scope(report, candidate)
}
dto! {pub struct Placement {pub role:String,pub input_id:String,pub start:usize,pub end:usize,pub reverse:bool,pub prefix:String,pub suffix:String,pub oriented_top:String}}
dto! {pub struct Materialized {
    pub context_revision:String,pub parent_context:String,pub request:Request,pub selection:Selection,
    pub candidate:Candidate,pub original_intent:Intent,pub terminal_fusions:Vec<String>,
    pub recipe:ligation_recipe::Recipe,pub designs:Vec<Option<primers::Row>>,pub placements:Vec<Placement>,pub warnings:Vec<String>,
}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("Golden Gate fusion placement: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn dna(s: &str) -> bool {
    s.bytes().all(|b| b"ACGT".contains(&b))
}
fn rc(s: &str) -> Result<String, Error> {
    if s.is_empty() {
        Ok(String::new())
    } else {
        Ok(Sequence::dna(s)?.reverse_complement()?.as_str().into())
    }
}
fn shape(r: &Request) -> Result<usize, Error> {
    let k = gate::dimensions(&r.base)?;
    let n = r.base.rows.len();
    let count = n - 1 + usize::from(r.base.base.circular);
    if n > MAX_ROWS
        || r.windows.len() > count
        || r.limits.max_states == 0
        || r.limits.max_states > MAX_STATES
        || r.limits.max_candidates == 0
        || r.limits.max_candidates > MAX_CANDIDATES
        || r.constraints.minimum_hamming > k
    {
        return Err(bad("invalid bounded search policy"));
    }
    let mut seen = BTreeSet::new();
    for w in &r.windows {
        if w.junction >= count
            || !seen.insert(w.junction)
            || w.min_shift > w.max_shift
            || w.min_shift < -MAX_SHIFT
            || w.max_shift > MAX_SHIFT
        {
            return Err(bad("invalid or duplicate signed placement window"));
        }
    }
    assembly_recipe::hash(r)?;
    Ok(k)
}
/// Borrowed bounds before input copies, hashing, source scans or native preparation.
pub fn preflight(inputs: &[(&str, &Molecule)], r: &Request) -> Result<(), Error> {
    shape(r)?;
    // Fusion preparation does not model independently declared outer termini.
    junction::inserts(&r.base.base, &r.base.edits)?;
    junction::intent_preflight(inputs, &r.base.base, &r.base.edits)?;
    let catalog = catalog_digest::builtin_catalog()?;
    let mut retained = Vec::new();
    let mut work = 0usize;
    for row in &r.base.rows {
        match row {
            Row::Retain(row) => {
                gate::retained_pin(row, &catalog, r.base.enzyme)?;
                retained.push(row.source.clone());
            }
            Row::DesignPcr {
                input_id,
                forward_length,
                reverse_length,
                ..
            } => {
                let m = inputs
                    .iter()
                    .find(|(id, _)| *id == input_id)
                    .map(|(_, m)| *m)
                    .ok_or_else(|| bad("missing explicitly bound original"))?;
                primers::physical_source(m)?;
                if !dna(m.sequence.as_str()) {
                    return Err(bad("selected PCR source requires concrete original DNA"));
                }
                work = work
                    .checked_add(
                        m.sequence
                            .len()
                            .checked_add(800)
                            .and_then(|n| n.checked_mul((forward_length + reverse_length + 32) * 8))
                            .ok_or(Error::TooLarge)?,
                    )
                    .filter(|n| *n <= 100_000_000)
                    .ok_or(Error::TooLarge)?;
            }
        }
    }
    if !retained.is_empty() {
        assembly_recipe::preparation_rows_preflight(inputs, &r.base.base.name, &retained)?;
    }
    assembly_recipe::hash(&(inputs, r))?;
    Ok(())
}
fn pins(inputs: &[PlanInput]) -> Result<Vec<(&str, &str)>, Error> {
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    for i in inputs {
        if !hex(&i.revision) || i.revision != i.molecule.revision() {
            return Err(bad("original input revision changed"));
        }
    }
    Ok(inputs
        .iter()
        .map(|i| (i.id.as_str(), i.revision.as_str()))
        .collect())
}
pub fn inspect_revision(inputs: &[PlanInput], r: &Request) -> Result<String, Error> {
    preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        r,
    )?;
    assembly_recipe::hash(&(ENGINE_VERSION, pins(inputs)?, r))
}
fn from_intent(intent: &Intent, r: &Request) -> Result<String, Error> {
    intent.validate_binding(&intent.context_revision)?;
    if assembly_recipe::hash(&intent.base_request)? != assembly_recipe::hash(&r.base.base)?
        || intent.edits != r.base.edits
    {
        return Err(bad("original independent intent changed"));
    }
    let raw = intent
        .target
        .metadata
        .get("assembly_junction_edit:inputs")
        .filter(|s| s.len() <= 16_384)
        .ok_or_else(|| bad("missing original revision pins"))?;
    let pins: Vec<(String, String)> =
        serde_json::from_str(raw).map_err(|_| bad("invalid source pins"))?;
    if pins.is_empty()
        || pins.len() > MAX_ROWS
        || pins
            .iter()
            .any(|(id, rev)| id.is_empty() || id.len() > 128 || !hex(rev))
    {
        return Err(bad("invalid bounded source pins"));
    }
    assembly_recipe::hash(&(ENGINE_VERSION, pins, r))
}
#[derive(Clone)]
struct Geometry {
    start: i64,
    end: i64,
    core_start: i64,
    core_end: i64,
}
fn geometries(
    intent: &Intent,
    r: &Request,
    shifts: &[i32],
    k: usize,
) -> Result<Option<Vec<Geometry>>, Error> {
    let n = r.base.rows.len();
    let count = n - 1 + usize::from(r.base.base.circular);
    if shifts.len() != count {
        return Err(bad("wrong ordered junction count"));
    }
    for (i, shift) in shifts.iter().enumerate() {
        let bounds = r
            .windows
            .iter()
            .find(|w| w.junction == i)
            .map(|w| (w.min_shift, w.max_shift))
            .unwrap_or((0, 0));
        if *shift < bounds.0 || *shift > bounds.1 {
            return Err(bad("candidate exceeds declared placement windows"));
        }
    }
    let mut result = Vec::with_capacity(n);
    for (i, (f, row)) in intent
        .adjusted_request
        .fragments
        .iter()
        .zip(&r.base.rows)
        .enumerate()
    {
        let old = &intent.fragments[i];
        let s = old.output_start as i64;
        let e = old.output_end as i64;
        let core_end = s + (f.end - f.start) as i64;
        let incoming = if i > 0 {
            shifts[i - 1]
        } else if r.base.base.circular {
            shifts[n - 1]
        } else {
            0
        };
        let outgoing = if i + 1 < n || r.base.base.circular {
            shifts[i]
        } else {
            0
        };
        let start = s + i64::from(incoming);
        let end = e + i64::from(outgoing);
        let cs = start.max(s);
        let ce = end.min(core_end);
        if end <= start || ce <= cs || cs - start > 200 || end - ce > 200 {
            return Ok(None);
        }
        if !r.base.base.circular && (start < 0 || end > intent.target.sequence.len() as i64) {
            return Ok(None);
        }
        match row {
            Row::Retain(_) => {
                if start != s || end != e {
                    return Ok(None);
                }
            }
            Row::DesignPcr {
                forward_length,
                reverse_length,
                ..
            } => {
                let core = (ce - cs) as usize;
                let annealed = core.checked_sub(if f.reverse { k } else { 0 });
                if core <= k || annealed.is_none_or(|a| a < forward_length + reverse_length) {
                    return Ok(None);
                }
                let tail = r.base.clamp.len()
                    + r.base.enzyme.profile().recognition.len()
                    + r.base.spacer.len();
                let prefix = (cs - start) as usize;
                let suffix = (end - ce) as usize;
                let (ft, rt) = if f.reverse {
                    (k + suffix, k + prefix)
                } else {
                    (prefix, k + suffix)
                };
                if forward_length + tail + ft > 200 || reverse_length + tail + rt > 200 {
                    return Ok(None);
                }
            }
        }
        result.push(Geometry {
            start,
            end,
            core_start: cs,
            core_end: ce,
        });
    }
    Ok(Some(result))
}
fn slice(intent: &Intent, start: i64, end: i64) -> Result<String, Error> {
    let bytes = intent.target.sequence.as_str().as_bytes();
    let len = end
        .checked_sub(start)
        .filter(|n| *n >= 0 && *n <= bytes.len() as i64)
        .ok_or_else(|| bad("placement exceeds one intended turn"))? as usize;
    if !intent.base_request.circular && (start < 0 || end > bytes.len() as i64) {
        return Err(bad("linear placement crosses a free end"));
    }
    if bytes.is_empty() {
        return Err(bad("empty independent target"));
    }
    let start = start.rem_euclid(bytes.len() as i64) as usize;
    Ok((0..len)
        .map(|i| bytes[(start + i) % bytes.len()] as char)
        .collect())
}
fn fusions(intent: &Intent, r: &Request, g: &[Geometry], k: usize) -> Result<Vec<String>, Error> {
    (0..r.base.rows.len() - 1 + usize::from(r.base.base.circular))
        .map(|i| {
            let next = &g[(i + 1) % g.len()];
            slice(intent, next.start, next.start + k as i64)
        })
        .collect()
}
fn fixed(
    inputs: &[PlanInput],
    r: &Request,
    intent: &Intent,
    k: usize,
) -> Result<(Vec<Option<(End, End)>>, Vec<String>), Error> {
    let rows = r
        .base
        .rows
        .iter()
        .enumerate()
        .filter_map(|(i, row)| {
            if let Row::Retain(row) = row {
                Some((i, row.clone()))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    let mut ends = vec![None; r.base.rows.len()];
    if !rows.is_empty() {
        let recipe = ligation_recipe::Recipe {
            name: r.base.base.name.clone(),
            fragments: rows.iter().map(|(_, r)| r.clone()).collect(),
            circular: false,
            allow_nicked: r.base.allow_nicked,
        };
        // Only unchanged manual preparations run here, once outside the search loop.
        // The final-only target assessment of this partial recipe is not consumed.
        let preview = ligation_recipe::preview(inputs, &intent.target, &recipe)?;
        for ((original, _), physical) in rows.iter().zip(&preview.physical) {
            let pair = physical
                .as_ref()
                .and_then(|p| p.ends.clone())
                .ok_or_else(|| bad("retained preparation has no established fixed ends"))?;
            if [&pair.0, &pair.1].iter().any(|e| {
                e.polarity != Polarity::FivePrime || e.overhang.len() != k || !dna(&e.overhang)
            }) {
                return Err(bad(
                    "retained preparation is outside the chosen cohesive-end geometry; it is not silently repaired",
                ));
            }
            ends[*original] = Some(pair);
        }
        if rows.iter().any(|(i, _)| ends[*i].is_none()) {
            return Err(bad(
                "retained preparation failed before fixed ends were established",
            ));
        }
    }
    let terminal = if r.base.base.circular {
        vec![]
    } else {
        let left = if let Some((left, _)) = &ends[0] {
            left.overhang.clone()
        } else {
            slice(intent, 0, k as i64)?
        };
        let last = ends.len() - 1;
        let right = if let Some((_, right)) = &ends[last] {
            rc(&right.overhang)?
        } else {
            r.base
                .terminal_right_overhang
                .clone()
                .ok_or_else(|| bad("linear PCR row needs its explicit terminal fusion"))?
        };
        vec![left, right]
    };
    Ok((ends, terminal))
}
fn compatible(keys: &[String], terminals: &[String], policy: &Constraints) -> Result<bool, Error> {
    let all = keys.iter().chain(terminals).collect::<Vec<_>>();
    for (i, key) in all.iter().enumerate() {
        let reverse = rc(key)?;
        if policy.reject_palindromes && key.as_str() == reverse {
            return Ok(false);
        }
        for other in all.iter().take(i) {
            let other_reverse = rc(other)?;
            if policy.distinct_reverse_complements
                && (key.as_str() == other.as_str() || key.as_str() == other_reverse)
            {
                return Ok(false);
            }
            let distance = key
                .bytes()
                .zip(other.bytes())
                .filter(|(a, b)| a != b)
                .count()
                .min(
                    key.bytes()
                        .zip(other_reverse.bytes())
                        .filter(|(a, b)| a != b)
                        .count(),
                );
            if distance < policy.minimum_hamming {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
fn identity(
    context: &str,
    shifts: &[i32],
    fusions: &[String],
    terminal: &[String],
) -> Result<String, Error> {
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        "candidate",
        context,
        shifts,
        fusions,
        terminal,
    ))
}
fn advance(values: &mut [usize], domains: &[Vec<i32>]) -> bool {
    for i in (0..values.len()).rev() {
        if values[i] + 1 < domains[i].len() {
            values[i] += 1;
            return true;
        }
        values[i] = 0;
    }
    false
}
pub fn inspect(inputs: &[PlanInput], r: &Request) -> Result<Report, Error> {
    let context_revision = inspect_revision(inputs, r)?;
    let k = shape(r)?;
    let intent = junction::build_intent(inputs, &r.base.base, &r.base.edits)?;
    if !dna(intent.target.sequence.as_str()) {
        return Err(bad(
            "fusion proposals require concrete independent intended DNA",
        ));
    }
    let (fixed, terminal) = fixed(inputs, r, &intent, k)?;
    let count = r.base.rows.len() - 1 + usize::from(r.base.base.circular);
    let mut domains = vec![vec![0]; count];
    for w in &r.windows {
        domains[w.junction] = (w.min_shift..=w.max_shift).collect();
    }
    for (i, domain) in domains.iter_mut().enumerate() {
        if matches!(&r.base.rows[i], Row::Retain(_))
            || matches!(&r.base.rows[(i + 1) % r.base.rows.len()], Row::Retain(_))
        {
            domain.retain(|shift| *shift == 0);
        }
        domain.sort_by_key(|shift| (shift.unsigned_abs(), *shift < 0));
    }
    let impossible = domains.iter().any(Vec::is_empty);
    let mut indices = vec![0; count];
    let mut candidates = Vec::new();
    let mut states = 0;
    let mut work = 0usize;
    let mut complete = impossible;
    let mut limits = Vec::new();
    while !impossible {
        let shifts = domains
            .iter()
            .zip(&indices)
            .map(|(domain, i)| domain[*i])
            .collect::<Vec<_>>();
        let keys = count + terminal.len();
        let charge = keys
            .checked_mul(keys + 1)
            .and_then(|n| n.checked_mul(k * 4))
            .and_then(|n| n.checked_add(r.base.rows.len() * 32))
            .ok_or(Error::TooLarge)?;
        if states >= r.limits.max_states
            || work.checked_add(charge).is_none_or(|n| n > MAX_SEARCH_WORK)
        {
            limits.push(
                "Declared assignment/comparison work limit reached; unexamined placements remain."
                    .into(),
            );
            break;
        }
        states += 1;
        work += charge;
        if let Some(geometry) = geometries(&intent, r, &shifts, k)? {
            let keys = fusions(&intent, r, &geometry, k)?;
            let mut manual_ok = true;
            for (i, ends) in fixed.iter().enumerate() {
                if let Some((left, right)) = ends {
                    let incoming = if i > 0 {
                        Some(&keys[i - 1])
                    } else if r.base.base.circular {
                        keys.last()
                    } else {
                        None
                    };
                    let outgoing = keys.get(i);
                    if incoming.is_some_and(|key| key != &left.overhang)
                        || outgoing
                            .is_some_and(|key| rc(key).ok().as_ref() != Some(&right.overhang))
                    {
                        manual_ok = false;
                        break;
                    }
                }
            }
            if manual_ok && compatible(&keys, &terminal, &r.constraints)? {
                candidates.push(Candidate {
                    id: identity(&context_revision, &shifts, &keys, &terminal)?,
                    shifts: shifts.clone(),
                    fusions: keys,
                });
            }
        }
        let more = advance(&mut indices, &domains);
        if !more {
            complete = true;
            break;
        }
        if candidates.len() >= r.limits.max_candidates {
            limits.push(
                "Candidate retention limit reached; further placements were not searched.".into(),
            );
            break;
        }
    }
    let report=Report {context_revision:context_revision.clone(),request:r.clone(),original_intent:intent,terminal_fusions:terminal,candidates,states_examined:states,search_complete:complete,limits,
        warnings:vec!["Signed placements preserve independently declared intended DNA. Copied flanks are explicit primer-introduced bases, not inherited source-coordinate annotations.".into(),"Palindrome, reverse-complement and Hamming rules are declared mathematical constraints, not measured ligation fidelity, reaction yield or exhaustive mixture uniqueness. Linear free cohesive ends are included as possible competitors.".into(),"Search does not establish selected-row PCR uniqueness or cleavage. Materialize one explicit candidate, then freshly verify its complete physical recipe against the unchanged independent target.".into(),"Unchanged manual rows safely fix both adjacent placement domains to zero. Remaining domains visit zero/nearest shifts first; an incompatible declared window can therefore have no assignments without exhaustive molecular claims.".into()]};
    report.validate_binding(&context_revision)?;
    Ok(report)
}
pub fn materialize_revision(
    inputs: &[PlanInput],
    r: &Request,
    s: &Selection,
) -> Result<String, Error> {
    let parent = inspect_revision(inputs, r)?;
    if s.context_revision != parent || !hex(&s.candidate_id) {
        return Err(bad("selection belongs to another source/request context"));
    }
    assembly_recipe::hash(&(ENGINE_VERSION, "materialize", parent, s))
}
fn placements(
    intent: &Intent,
    r: &Request,
    shifts: &[i32],
    k: usize,
) -> Result<Vec<Placement>, Error> {
    let geometry = geometries(intent, r, shifts, k)?
        .ok_or_else(|| bad("placement no longer preserves its original source cores"))?;
    intent
        .adjusted_request
        .fragments
        .iter()
        .zip(geometry)
        .enumerate()
        .map(|(i, (original, g))| {
            let old = &intent.fragments[i];
            let left = (g.core_start - old.output_start as i64) as usize;
            let right = (old.output_start + (original.end - original.start)) as i64 - g.core_end;
            let mut start = original.start
                + if original.reverse {
                    right as usize
                } else {
                    left
                };
            let mut end = original.end
                - if original.reverse {
                    left
                } else {
                    right as usize
                };
            // Existing independent intent stores original lengths in bounded metadata.
            let raw = intent
                .target
                .metadata
                .get("assembly_junction_edit:source_lengths")
                .ok_or_else(|| bad("missing original source lengths"))?;
            let lengths: Vec<usize> =
                serde_json::from_str(raw).map_err(|_| bad("invalid source lengths"))?;
            let n = *lengths.get(i).ok_or_else(|| bad("missing source length"))?;
            if start >= n {
                start -= n;
                end -= n;
            }
            Ok(Placement {
                role: original.role.clone(),
                input_id: original.input_id.clone(),
                start,
                end,
                reverse: original.reverse,
                prefix: slice(intent, g.start, g.core_start)?,
                suffix: slice(intent, g.core_end, g.end)?,
                oriented_top: slice(intent, g.start, g.end)?,
            })
        })
        .collect()
}
pub fn materialize(
    inputs: &[PlanInput],
    r: &Request,
    s: &Selection,
) -> Result<Materialized, Error> {
    let context_revision = materialize_revision(inputs, r, s)?;
    let report = inspect(inputs, r)?;
    let candidate = report
        .candidates
        .iter()
        .find(|c| c.id == s.candidate_id)
        .cloned()
        .ok_or_else(|| bad("explicit candidate was not retained by this exact bounded search"))?;
    let k = shape(r)?;
    let placed = placements(&report.original_intent, r, &candidate.shifts, k)?;
    let catalog = catalog_digest::builtin_catalog()?;
    let selection = catalog_digest::select(
        &catalog,
        &[r.base.enzyme.profile().name.to_ascii_lowercase()],
    )?;
    let enzymes = selection.resolved_enzymes(&catalog)?;
    let specs = placed
        .iter()
        .enumerate()
        .filter_map(|(i, p)| {
            if let Row::DesignPcr {
                forward_length,
                reverse_length,
                ..
            } = &r.base.rows[i]
            {
                Some(primers::Fragment {
                    role: p.role.clone(),
                    input_id: p.input_id.clone(),
                    start: p.start,
                    end: p.end,
                    reverse: p.reverse,
                    forward_length: *forward_length,
                    reverse_length: *reverse_length,
                    right_overhang: candidate
                        .fusions
                        .get(i)
                        .cloned()
                        .unwrap_or_else(|| report.terminal_fusions[1].clone()),
                })
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    let template = primers::Request {
        name: r.base.base.name.clone(),
        fragments: specs,
        enzyme: r.base.enzyme,
        clamp: r.base.clamp.clone(),
        spacer: r.base.spacer.clone(),
    };
    let mut designs = Vec::new();
    let mut fragments = Vec::new();
    for (i, (row, p)) in r.base.rows.iter().zip(&placed).enumerate() {
        match row {
            Row::Retain(row) => {
                fragments.push(row.clone());
                designs.push(None);
            }
            Row::DesignPcr { .. } => {
                let original = inputs
                    .iter()
                    .find(|v| v.id == p.input_id)
                    .ok_or_else(|| bad("selected source unavailable"))?;
                let spec = template
                    .fragments
                    .iter()
                    .find(|f| f.role == p.role)
                    .ok_or_else(|| bad("selected role unavailable"))?;
                let (fragment, evidence) = primers::design_row_with_flanks(
                    original, &template, spec, &p.prefix, &p.suffix, &selection, &enzymes,
                )?;
                if evidence.oriented_top != p.oriented_top {
                    return Err(bad(
                        "fresh physical preparation differs from the frozen independent placement",
                    ));
                }
                fragments.push(ligation_recipe::Fragment {
                    source: fragment,
                    treatment: None,
                });
                designs.push(Some(evidence));
            }
        }
        let _ = i;
    }
    let recipe = ligation_recipe::Recipe {
        name: r.base.base.name.clone(),
        fragments,
        circular: r.base.base.circular,
        allow_nicked: r.base.allow_nicked,
    };
    ligation_recipe::preflight(inputs, &report.original_intent.target, &recipe)?;
    let result = Materialized {
        context_revision: context_revision.clone(),
        parent_context: report.context_revision,
        request: r.clone(),
        selection: s.clone(),
        candidate,
        original_intent: report.original_intent,
        terminal_fusions: report.terminal_fusions,
        recipe,
        designs,
        placements: placed,
        warnings: report.warnings,
    };
    result.validate_binding(&context_revision)?;
    Ok(result)
}
impl Report {
    /// Structural binding only. Imported proposals never authorize native chemistry.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let k = shape(&self.request)?;
        if !hex(expected)
            || self.context_revision != expected
            || from_intent(&self.original_intent, &self.request)? != expected
            || self.candidates.len() > self.request.limits.max_candidates
            || self.states_examined > self.request.limits.max_states
            || self.states_examined == 0 && !self.candidates.is_empty()
            || self.search_complete != self.limits.is_empty()
            || self.terminal_fusions.len()
                != if self.request.base.base.circular {
                    0
                } else {
                    2
                }
            || self
                .terminal_fusions
                .iter()
                .any(|f| f.len() != k || !dna(f))
        {
            return Err(bad("search report binding/counts changed"));
        }
        let mut ids = BTreeSet::new();
        for c in &self.candidates {
            let g = geometries(&self.original_intent, &self.request, &c.shifts, k)?
                .ok_or_else(|| bad("reported placement is invalid"))?;
            if c.fusions != fusions(&self.original_intent, &self.request, &g, k)?
                || c.id != identity(expected, &c.shifts, &c.fusions, &self.terminal_fusions)?
                || !ids.insert(&c.id)
                || !compatible(
                    &c.fusions,
                    &self.terminal_fusions,
                    &self.request.constraints,
                )?
            {
                return Err(bad("candidate identity/fusion/constraint changed"));
            }
        }
        texts(&self.limits, 8)?;
        texts(&self.warnings, 16)?;
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
fn texts(values: &[String], max: usize) -> Result<(), Error> {
    if values.len() > max || values.iter().any(|s| s.len() > 4096) {
        Err(Error::TooLarge)
    } else {
        Ok(())
    }
}
impl Materialized {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let k = shape(&self.request)?;
        let parent = from_intent(&self.original_intent, &self.request)?;
        if !hex(expected)
            || self.context_revision != expected
            || self.parent_context != parent
            || self.selection.context_revision != parent
            || self.selection.candidate_id != self.candidate.id
            || self.candidate.id
                != identity(
                    &parent,
                    &self.candidate.shifts,
                    &self.candidate.fusions,
                    &self.terminal_fusions,
                )?
            || expected
                != assembly_recipe::hash(&(
                    ENGINE_VERSION,
                    "materialize",
                    &parent,
                    &self.selection,
                ))?
            || self.recipe.name != self.request.base.base.name
            || self.recipe.circular != self.request.base.base.circular
            || self.recipe.allow_nicked != self.request.base.allow_nicked
            || self.designs.len() != self.request.base.rows.len()
            || self.placements.len() != self.designs.len()
            || self.recipe.fragments.len() != self.designs.len()
            || self.terminal_fusions.len()
                != if self.request.base.base.circular {
                    0
                } else {
                    2
                }
        {
            return Err(bad("materialized binding/recipe changed"));
        }
        let expected_placements = placements(
            &self.original_intent,
            &self.request,
            &self.candidate.shifts,
            k,
        )?;
        if assembly_recipe::hash(&expected_placements)? != assembly_recipe::hash(&self.placements)?
        {
            return Err(bad(
                "materialized original-source intervals or copied flanks changed",
            ));
        }
        let g = geometries(
            &self.original_intent,
            &self.request,
            &self.candidate.shifts,
            k,
        )?
        .ok_or_else(|| bad("invalid placed rows"))?;
        if self.candidate.fusions != fusions(&self.original_intent, &self.request, &g, k)?
            || !compatible(
                &self.candidate.fusions,
                &self.terminal_fusions,
                &self.request.constraints,
            )?
        {
            return Err(bad("materialized fusion set changed"));
        }
        let catalog = catalog_digest::builtin_catalog()?;
        let pin = catalog_digest::select(
            &catalog,
            &[self.request.base.enzyme.profile().name.to_ascii_lowercase()],
        )?;
        for (i, (((original, row), design), placement)) in self
            .request
            .base
            .rows
            .iter()
            .zip(&self.recipe.fragments)
            .zip(&self.designs)
            .zip(&self.placements)
            .enumerate()
        {
            match original {
                Row::Retain(original) => {
                    if design.is_some()
                        || assembly_recipe::hash(original)? != assembly_recipe::hash(row)?
                    {
                        return Err(bad("manual preparation changed"));
                    }
                }
                Row::DesignPcr {
                    forward_length,
                    reverse_length,
                    ..
                } => {
                    let d = design
                        .as_ref()
                        .ok_or_else(|| bad("missing selected PCR evidence"))?;
                    let Preparation::PcrCatalogDigest {
                        forward,
                        reverse,
                        selection,
                        fragment_index,
                        repair,
                    } = &row.source.preparation
                    else {
                        return Err(bad("selected preparation is not actual pinned PCR/digest"));
                    };
                    let fusion = self
                        .candidate
                        .fusions
                        .get(i)
                        .or_else(|| self.terminal_fusions.get(1))
                        .ok_or_else(|| bad("missing right fusion"))?;
                    if d.role != placement.role
                        || d.input_id != placement.input_id
                        || row.source.role != placement.role
                        || row.source.input_id != placement.input_id
                        || row.source.reverse != placement.reverse
                        || row.treatment.is_some()
                        || repair.is_some()
                        || d.fragment_index != *fragment_index
                        || d.forward.annealing.len() != *forward_length
                        || d.reverse.annealing.len() != *reverse_length
                        || d.oriented_top != placement.oriented_top
                        || d.forward.tail != d.forward_tail.sequence()
                        || d.reverse.tail != d.reverse_tail.sequence()
                        || assembly_recipe::hash(forward)? != assembly_recipe::hash(&d.forward)?
                        || assembly_recipe::hash(reverse)? != assembly_recipe::hash(&d.reverse)?
                        || assembly_recipe::hash(selection)? != assembly_recipe::hash(&pin)?
                        || assembly_recipe::hash(&d.selection)? != assembly_recipe::hash(&pin)?
                        || d.left.polarity != Polarity::FivePrime
                        || d.right.polarity != Polarity::FivePrime
                        || d.left.chemistry != Terminus::CLEAVED
                        || d.right.chemistry != Terminus::CLEAVED
                        || d.left.overhang != placement.oriented_top[..k]
                        || d.right.overhang != rc(fusion)?
                    {
                        return Err(bad("actual PCR/digest row changed"));
                    }
                    if forward.ordered_sequence().len() > 200
                        || reverse.ordered_sequence().len() > 200
                    {
                        return Err(Error::TooLarge);
                    }
                }
            }
        }
        // An independent letter check of the declared partition, including circle
        // origin movement. This never derives intended DNA from physical products.
        let actual = self
            .placements
            .iter()
            .map(|p| p.oriented_top.as_str())
            .collect::<String>();
        let first = g[0].start;
        let expected_letters = slice(
            &self.original_intent,
            first,
            first + self.original_intent.target.sequence.len() as i64,
        )?;
        if actual != expected_letters {
            return Err(bad("placement changed independent intended DNA"));
        }
        texts(&self.warnings, 16)?;
        assembly_recipe::hash(self)?;
        Ok(())
    }
}
