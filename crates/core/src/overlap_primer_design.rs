//! Deterministic anchored PCR primer selection with explicit overlap tails.
//! Original source intervals determine all sequence; the independent target is
//! bound for later verification and never selects domains, tails or endpoints.
use crate::{
    Alphabet, Error, Location, Molecule, Sequence, Topology,
    assembly_recipe::{self, FragmentSpec, Preparation, Recipe},
    chemistry,
    overlap_profiles::ProfileOptions,
    primers::{self, DesignOptions, DesignReport},
    verification::PlanInput,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE_VERSION: &str = "overlap-primer-design-1";
const MAX_ROWS: usize = 32;
const MAX_BASES: usize = 8 << 20;
const MAX_WORK: usize = 100_000_000;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Fragment {pub role:String,pub input_id:String,pub start:usize,pub end:usize,pub reverse:bool}}
dto! {pub struct Request {pub name:String,pub fragments:Vec<Fragment>,pub overlap:usize,pub circular:bool,pub primer_options:DesignOptions,pub profile_options:ProfileOptions}}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub context_revision: String,
    pub source_revisions: Vec<String>,
    pub request: Request,
    pub recipe: Recipe,
    pub designs: Vec<DesignReport>,
    pub warnings: Vec<String>,
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("overlap primer planning: {s}"))
}
fn token(s: &str, cap: usize) -> Result<(), Error> {
    if s.trim().is_empty() || s.len() > cap || s.chars().any(char::is_control) {
        Err(bad("invalid bounded source/role/name"))
    } else {
        Ok(())
    }
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
pub(crate) fn options(request: &Request) -> Result<(), Error> {
    if request.fragments.is_empty()
        || request.fragments.len() > MAX_ROWS
        || request.overlap == 0
        || request.overlap > 80
    {
        return Err(bad("requires1..32 rows and explicit overlap1..80"));
    }
    token(&request.name, 1024)?;
    request.profile_options.validate()?;
    let junctions = request.fragments.len() - 1 + usize::from(request.circular);
    if request
        .overlap
        .checked_mul(junctions)
        .is_none_or(|n| n > request.profile_options.max_comparisons)
    {
        return Err(bad(
            "explicit overlaps exceed the later profile comparison budget",
        ));
    }

    let o = &request.primer_options;
    if !o.conditions.sodium_molar.is_finite()
        || !(0.01..=1.0).contains(&o.conditions.sodium_molar)
        || !o.conditions.total_strand_molar.is_finite()
        || !(1e-12..=1e-2).contains(&o.conditions.total_strand_molar)
    {
        return Err(bad(
            "existing Tm profile requires sodium0.01..1M and totalstrands1pM..10mM",
        ));
    }
    if !o.forward_tail.is_empty() || !o.reverse_tail.is_empty() {
        return Err(bad(
            "preexisting primer-option tails must be empty; no user sequence is silently replaced",
        ));
    }
    // These are the existing engine's public filter bounds. Target endpoints and
    // output count are deliberately replaced by each explicit row and one pair.
    if o.min_length < 8
        || o.max_length > 80
        || o.min_length > o.max_length
        || o.max_length - o.min_length > 63
        || !o.target_tm_c.is_finite()
        || !(-50.0..=100.0).contains(&o.target_tm_c)
        || !o.max_tm_difference_c.is_finite()
        || !(0.0..=100.0).contains(&o.max_tm_difference_c)
        || !o.min_gc_fraction.is_finite()
        || !o.max_gc_fraction.is_finite()
        || !(0.0..=1.0).contains(&o.min_gc_fraction)
        || !(o.min_gc_fraction..=1.0).contains(&o.max_gc_fraction)
        || o.max_homopolymer == 0
        || o.max_homopolymer > 200
        || o.max_three_prime_complementarity > 200
    {
        return Err(bad(
            "invalid primer length/Tm/GC/complementarity filter bounds",
        ));
    }
    Ok(())
}
pub(crate) fn pcr_source(m: &Molecule) -> Result<(), Error> {
    let raw_marker = m.original_genbank.as_ref().is_some_and(|raw| {
        raw.lines()
            .filter_map(|line| {
                line.strip_prefix("COMMENT     ")
                    .or_else(|| line.strip_prefix("CC   "))
            })
            .any(|value| {
                matches!(
                    value.trim_end(),
                    "Collagene physical product v1: unverified sequence-only evidence"
                        | "Collagene reference-only v1: unverified sequence-only evidence"
                )
            })
    });
    if raw_marker
        || chemistry::has_unmapped_physical_evidence(m)
        || m.metadata.keys().any(|k| {
            k == crate::recombination::REFERENCE_ONLY_KEY
                || k == crate::recombination::RECIPE_KEY
                || k == "source:recombination:reference_only"
                || k == "source:recombination:unverified_recipe"
                || k == "source:legacy:chemistry_unmapped"
                || k == "source:snapgene:chemistry_unmapped"
                || k == "tracked_modifications"
                || k == "source:declared_termini"
                || k.starts_with("topo:")
                || k.starts_with("variant:")
                || k.starts_with("restriction_inspection:")
                || k.starts_with("source:restriction_inspection:")
                || k.starts_with("restriction_catalog:")
                || k.starts_with("source:restriction_catalog:")
        })
    {
        return Err(bad(
            "physical/reference-only/variant precursor evidence cannot supply a new PCR substrate; reconstruct the original physical graph",
        ));
    }
    Ok(())
}
/// Borrowed source/interval/filter guards before any original model copies or
/// primer selection. Concrete cores do not establish unique binding elsewhere.
pub fn preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    request: &Request,
) -> Result<(), Error> {
    options(request)?;
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    let mut models = inputs.iter().map(|(_, m)| *m).collect::<Vec<_>>();
    models.push(target);
    crate::model::bounds::model_preflight(&models)?;
    if target.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("independent target must be DNA"));
    }
    target.validate()?;
    let mut lookup = BTreeMap::new();
    let mut bases = target.sequence.len();
    for &(id, m) in inputs {
        token(id, 128)?;
        if lookup.insert(id, m).is_some() {
            return Err(bad("duplicate original input ID"));
        }
        if m.sequence.alphabet() != Alphabet::Dna {
            return Err(bad("PCR sources must be DNA"));
        }
        m.validate()?;
        bases = bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= MAX_BASES)
            .ok_or(Error::TooLarge)?;
    }
    let mut roles = BTreeSet::new();
    let mut work = 0usize;
    let mut scans = 0usize;
    let count = request.primer_options.max_length - request.primer_options.min_length + 1;
    for (index, f) in request.fragments.iter().enumerate() {
        token(&f.role, 128)?;
        token(&f.input_id, 128)?;
        if !roles.insert(f.role.as_str()) {
            return Err(bad(
                "distinct recipe roles are required, even when an original is reused",
            ));
        }
        let m = lookup
            .get(f.input_id.as_str())
            .copied()
            .ok_or_else(|| bad("requested row has no explicitly bound original"))?;
        pcr_source(m)?;
        let n = m.sequence.len();
        let length = f
            .end
            .checked_sub(f.start)
            .filter(|length| *length > 0 && *length <= n)
            .ok_or_else(|| {
                bad("requested interval must contain at most one original-source turn")
            })?;
        if f.start >= n
            || f.end > f.start.checked_add(n).ok_or(Error::TooLarge)?
            || m.topology == Topology::Linear && f.end > n
        {
            return Err(bad(
                "invalid source interval; only circular sources may cross the origin",
            ));
        }
        if length <= request.overlap
            || length
                < request
                    .primer_options
                    .min_length
                    .checked_mul(2)
                    .ok_or(Error::TooLarge)?
        {
            return Err(bad(
                "core must exceed the shared homology and fit both nonoverlapping anchored annealing domains",
            ));
        }
        for j in 0..length {
            if !b"ACGT".contains(&m.sequence.as_str().as_bytes()[(f.start + j) % n]) {
                return Err(bad(
                    "requested core contains ambiguous calls; concrete primer/tail sequence is withheld",
                ));
            }
        }
        let outgoing = index + 1 < request.fragments.len() || request.circular;
        let a = request.primer_options.max_length;
        let b = a + if outgoing { request.overlap } else { 0 };
        let pair = count
            .checked_mul(count)
            .and_then(|n| n.checked_mul(a))
            .and_then(|n| n.checked_mul(b))
            .and_then(|n| n.checked_mul(2))
            .ok_or(Error::TooLarge)?;
        if pair > 20_000_000 {
            return Err(bad(
                "requested length interval/tails exceed the existing per-row primer-pair work bound; narrow lengths",
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
            .filter(|n| *n <= MAX_WORK)
            .ok_or_else(|| {
                bad("aggregate primer structure/pair planning work exceeds100million units")
            })?;
        scans = scans
            .checked_add(
                n.checked_mul(a)
                    .and_then(|n| n.checked_mul(8))
                    .ok_or(Error::TooLarge)?,
            )
            .filter(|n| *n <= MAX_WORK)
            .ok_or_else(|| {
                bad("aggregate later exact-PCR scan reservation exceeds100million comparisons")
            })?;
    }
    assembly_recipe::hash(&(inputs, target, request))?;
    Ok(())
}
pub fn design_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    request: &Request,
) -> Result<String, Error> {
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        target,
        request,
    )?;
    if inputs.iter().any(|i| !hex(&i.revision)) {
        return Err(bad("explicit canonical source revisions required"));
    }
    let target_revision = assembly_recipe::bound_revisions(inputs, target)?;
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        inputs
            .iter()
            .map(|i| (&i.id, &i.revision))
            .collect::<Vec<_>>(),
        target_revision,
        request,
    ))
}
fn interval(start: usize, length: usize, n: usize) -> Location {
    let start = start % n;
    let first = length.min(n - start);
    if first == length {
        Location::range(start, start + length)
    } else {
        Location::Join(vec![
            Location::range(start, n),
            Location::range(0, length - first),
        ])
    }
}
fn prefix(m: &Molecule, f: &Fragment, h: usize) -> Result<String, Error> {
    let l = if f.reverse {
        Location::Complement(Box::new(interval(f.end - h, h, m.sequence.len())))
    } else {
        interval(f.start, h, m.sequence.len())
    };
    l.extract(&m.sequence)
}
pub fn design(inputs: &[PlanInput], target: &Molecule, request: &Request) -> Result<Report, Error> {
    let context_revision = design_revision(inputs, target, request)?;
    let originals = inputs
        .iter()
        .map(|i| (i.id.as_str(), &i.molecule))
        .collect::<BTreeMap<_, _>>();
    let mut designs = Vec::with_capacity(request.fragments.len());
    let mut fragments = Vec::with_capacity(request.fragments.len());
    for (index, f) in request.fragments.iter().enumerate() {
        let m = originals
            .get(f.input_id.as_str())
            .copied()
            .ok_or_else(|| bad("missing original"))?;
        let mut o = request.primer_options.clone();
        o.target_start = f.start;
        o.target_end = f.end;
        o.max_pairs = 1;
        if index + 1 < request.fragments.len() || request.circular {
            let next = &request.fragments[(index + 1) % request.fragments.len()];
            let next_source = originals
                .get(next.input_id.as_str())
                .copied()
                .ok_or_else(|| bad("missing downstream original"))?;
            let tail = Sequence::dna(&prefix(next_source, next, request.overlap)?)?
                .reverse_complement()?
                .as_str()
                .to_owned();
            if f.reverse {
                o.forward_tail = tail;
            } else {
                o.reverse_tail = tail;
            }
        }
        let mut report = primers::design_primers(m, &o)?;
        let pair=report.pairs.first_mut().ok_or_else(||bad("no anchored primer pair passes the explicit filters; adjust criteria or source interval deliberately"))?;
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
        designs.push(report);
    }
    let recipe = Recipe {
        name: request.name.clone(),
        fragments,
        overlaps: vec![
            request.overlap;
            request.fragments.len() - 1 + usize::from(request.circular)
        ],
        circular: request.circular,
        options: request.profile_options.clone(),
    };
    assembly_recipe::preflight(inputs, target, &recipe)?;
    let result=Report{context_revision,source_revisions:inputs.iter().map(|i|i.revision.clone()).collect(),request:request.clone(),recipe,designs,warnings:vec!["Deterministic anchored endpoint/heuristic primer selection only. All cores, overlap tails and orientations come from explicitly selected original-source intervals; the independent target never chooses sequence.".into(),"Each upstream oriented core receives the next oriented core's prefix at its right end. The last linear core has no synthetic tail; circular closure extends toward the first prefix and retains the requested first origin.".into(),"One highest-ranked pair per row is returned. Ranking and exact interval geometry do not establish unique PCR, mixture specificity, reagent compatibility or experimental success. Freshly verify the generated recipe against the independently supplied target, including every alternative PCR geometry.".into(),"Request primer-option target endpoints, max_pairs and empty tail slots are replaced explicitly by row endpoints, max_pairs1 and calculated overlaps; all other thermodynamic/filter settings are preserved.".into()]};
    result.validate()?;
    Ok(result)
}
fn location_guard(location: &Location) -> Result<(), Error> {
    let mut stack = vec![(location, 0usize)];
    let mut nodes = 0usize;
    while let Some((l, depth)) = stack.pop() {
        nodes += 1;
        if nodes > 16384 || depth > 32 {
            return Err(Error::TooLarge);
        }
        match l {
            Location::Range {
                start,
                end,
                before: false,
                after: false,
            } if start < end && *end <= MAX_BASES => {}
            Location::Complement(inner) => stack.push((inner, depth + 1)),
            Location::Join(parts) if parts.len() <= 16384 => {
                for p in parts {
                    stack.push((p, depth + 1))
                }
            }
            _ => return Err(bad("unsupported returned primer coordinate geometry")),
        }
    }
    Ok(())
}
impl Report {
    /// Readonly returned-recipe consistency, never design/PCR replay or authority.
    pub fn validate(&self) -> Result<(), Error> {
        options(&self.request)?;
        let n = self.request.fragments.len();
        if !hex(&self.context_revision)
            || self.source_revisions.is_empty()
            || self.source_revisions.len() > MAX_ROWS
            || self.source_revisions.iter().any(|s| !hex(s))
            || self.designs.len() != n
            || self.recipe.fragments.len() != n
            || self.recipe.name != self.request.name
            || self.recipe.circular != self.request.circular
            || self.recipe.overlaps
                != vec![self.request.overlap; n - 1 + usize::from(self.request.circular)]
            || assembly_recipe::hash(&self.recipe.options)?
                != assembly_recipe::hash(&self.request.profile_options)?
            || self.warnings.len() > 32
            || self.warnings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad("invalid report/request/recipe shape"));
        }
        let mut roles = BTreeSet::new();
        for (index, ((f, row), d)) in self
            .request
            .fragments
            .iter()
            .zip(&self.recipe.fragments)
            .zip(&self.designs)
            .enumerate()
        {
            token(&f.role, 128)?;
            token(&f.input_id, 128)?;
            if !roles.insert(f.role.as_str())
                || f.end
                    .checked_sub(f.start)
                    .is_none_or(|n| n <= self.request.overlap || n > MAX_BASES)
                || f.start >= MAX_BASES
                || f.end > 2 * MAX_BASES
            {
                return Err(bad("invalid returned source interval/role"));
            }
            if f.role != row.role
                || f.input_id != row.input_id
                || f.reverse != row.reverse
                || d.pairs.len() != 1
                || d.options.max_pairs != 1
                || d.options.target_start != f.start
                || d.options.target_end != f.end
                || !self.source_revisions.contains(&d.source_revision)
                || d.warnings.len() > 32
                || d.warnings.iter().any(|s| s.len() > 8192)
            {
                return Err(bad("returned original/interval/pair mapping changed"));
            }
            let Preparation::ExactPcr {
                forward,
                reverse,
                template_strand: crate::chemistry::Strand::Top,
                ..
            } = &row.preparation
            else {
                return Err(bad(
                    "overlap planning must return actual exact-PCR preparations",
                ));
            };
            let pair = &d.pairs[0];
            for p in [forward, reverse] {
                if p.name.len() > 256
                    || p.annealing.len() + p.tail.len() > 200
                    || p.annealing.alphabet() != Alphabet::Dna
                    || !p.ordered_sequence().bytes().all(|b| b"ACGT".contains(&b))
                {
                    return Err(bad("invalid concrete returned primer"));
                }
            }
            if assembly_recipe::hash(forward)? != assembly_recipe::hash(&pair.forward.primer)?
                || assembly_recipe::hash(reverse)? != assembly_recipe::hash(&pair.reverse.primer)?
                || forward.name != format!("{}-forward", f.role)
                || reverse.name != format!("{}-reverse", f.role)
                || pair.forward.analysis.name != forward.name
                || pair.reverse.analysis.name != reverse.name
                || d.options.forward_tail != forward.tail
                || d.options.reverse_tail != reverse.tail
            {
                return Err(bad("returned recipe oligos disagree with ranked pair"));
            }
            for candidate in [&pair.forward, &pair.reverse] {
                if candidate.analysis.annealing != candidate.primer.annealing.as_str()
                    || candidate.analysis.tail != candidate.primer.tail
                    || candidate.analysis.ordered_sequence != candidate.primer.ordered_sequence()
                {
                    return Err(bad(
                        "returned primer analysis fields disagree with the actual oligo",
                    ));
                }
            }
            let core_length = f.end - f.start;
            if pair.forward.primer.annealing.len() + pair.reverse.primer.annealing.len()
                > core_length
                || core_length
                    .checked_add(forward.tail.len())
                    .and_then(|n| n.checked_add(reverse.tail.len()))
                    != Some(pair.product_length)
            {
                return Err(bad(
                    "returned anchored pair length differs from requested core/tails",
                ));
            }
            let outgoing = index + 1 < n || self.request.circular;
            let (tail, empty) = if f.reverse {
                (&forward.tail, &reverse.tail)
            } else {
                (&reverse.tail, &forward.tail)
            };
            if !empty.is_empty() || tail.len() != if outgoing { self.request.overlap } else { 0 } {
                return Err(bad("returned synthetic-tail side/length changed"));
            }
            let mut expected = self.request.primer_options.clone();
            expected.target_start = f.start;
            expected.target_end = f.end;
            expected.max_pairs = 1;
            expected.forward_tail = forward.tail.clone();
            expected.reverse_tail = reverse.tail.clone();
            if expected != d.options {
                return Err(bad("returned thermodynamic/filter settings changed"));
            }
            location_guard(&pair.forward.location)?;
            location_guard(&pair.reverse.location)?;
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad(
                "design belongs to another source/target/request context",
            ));
        }
        Ok(())
    }
}
