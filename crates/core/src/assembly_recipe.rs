//! Compile explicit preparation rows into existing typed operations. Fresh
//! reconstruction remains the verifier's job; this module never slices a source
//! document and calls that physical preparation, or searches toward a target.
use crate::{
    Alphabet, Error, Molecule, Topology,
    analysis::Primer,
    overlap_profiles::{self, ProfileOptions, ProfileReport},
    verification::{
        self, DesignPlan, Operation, OrientedSource, PlanInput, Source, VerificationReport,
    },
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ENGINE_VERSION: &str = "overlap-recipe-1";
pub const MAX_ROWS: usize = 128;
const MAX_BYTES: usize = 32 << 20;
const MAX_INPUT_BASES: usize = 8 << 20;
const MAX_ROLE_BASES: usize = 16 << 20;
const MAX_PREPARATION_WORK: usize = 100_000_000;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub enum Preparation {WholeLinear,Digest{enzymes:Vec<String>,fragment_index:usize,repair:Option<crate::chemistry::EndRepairMethod>},ExactPcr{forward:Primer,reverse:Primer}}}
dto! {pub struct FragmentSpec {pub role:String,pub input_id:String,pub preparation:Preparation,pub reverse:bool}}
dto! {pub struct Recipe {pub name:String,pub fragments:Vec<FragmentSpec>,pub overlaps:Vec<usize>,pub circular:bool,pub options:ProfileOptions}}
dto! {pub struct PreparedRow {pub row_index:usize,pub role:String,pub input_id:String,pub input_revision:String,pub preparation_step:Option<usize>,pub repair_step:Option<usize>,pub prepared_source:Source,pub reverse:bool}}
dto! {pub struct Compiled {pub plan:DesignPlan,pub rows:Vec<PreparedRow>,pub assembly_step:usize}}
dto! {pub struct Preview {pub context_revision:String,pub compiled:Compiled,pub verification:VerificationReport,
/// Fresh prepared documents in recipe order, before optional orientation.
/// On a failed preparation only the available leading row prefix is returned.
pub prepared:Vec<Molecule>,pub profile:Option<ProfileReport>,pub preparation_complete:bool}}
fn bad(s: &str) -> Error {
    Error::Format(format!("overlap recipe: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn token(s: &str, cap: usize) -> Result<(), Error> {
    if s.is_empty() || s.len() > cap || s.chars().any(char::is_control) {
        Err(bad("invalid bounded source/role/name"))
    } else {
        Ok(())
    }
}
struct Hash {
    hash: blake3::Hasher,
    bytes: usize,
}
impl Write for Hash {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(b.len())
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| std::io::Error::other("bounded assembly recipe"))?;
        self.hash.update(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn hash(v: &impl Serialize) -> Result<String, Error> {
    let mut w = Hash {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(&mut w, v).map_err(|_| Error::TooLarge)?;
    Ok(w.hash.finalize().to_hex().to_string())
}
fn primer(p: &Primer) -> Result<(), Error> {
    token(&p.name, 256)?;
    if p.annealing.alphabet() != Alphabet::Dna
        || p.annealing.is_empty()
        || p.annealing.len() > 1000
        || p.tail.len() > 1000
        || !p
            .annealing
            .as_str()
            .bytes()
            .chain(p.tail.bytes())
            .all(|b| b"ACGT".contains(&b))
    {
        return Err(bad(
            "exact-PCR primers require concrete DNA, annealing1..1000 and tail0..1000 bases",
        ));
    }
    Ok(())
}
/// Borrowed guards for hosts before cloning source documents or hashing revisions.
/// Input IDs identify originals; distinct row roles may reuse the same original.
/// Document/event identity is captured by the host and is never inferred here.
pub fn source_preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<(), Error> {
    if inputs.is_empty()
        || inputs.len() > MAX_ROWS
        || recipe.fragments.is_empty()
        || recipe.fragments.len() > MAX_ROWS
    {
        return Err(Error::TooLarge);
    }
    token(&recipe.name, 1024)?;
    recipe.options.validate()?;
    if recipe.overlaps.len() != recipe.fragments.len() - 1 + usize::from(recipe.circular)
        || recipe.overlaps.iter().any(|n| *n == 0 || *n > 10_000)
    {
        return Err(bad(
            "one explicit positive overlap1..10000 is required per junction",
        ));
    }
    let total_homology = recipe
        .overlaps
        .iter()
        .try_fold(0usize, |a, b| a.checked_add(*b).ok_or(Error::TooLarge))?;
    if total_homology > 1_000_000 || total_homology > recipe.options.max_comparisons {
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
            return Err(bad("source preparation requires DNA"));
        }
        m.validate()?;
        bases = bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= MAX_INPUT_BASES)
            .ok_or(Error::TooLarge)?;
    }
    let mut roles = BTreeSet::new();
    let mut role_bases = 0usize;
    let mut work = 0usize;
    let catalog = crate::analysis::starter_enzymes();
    for row in &recipe.fragments {
        token(&row.role, 128)?;
        token(&row.input_id, 128)?;
        if !roles.insert(row.role.as_str()) {
            return Err(bad(
                "duplicate row role ID; reused originals require distinct role identities",
            ));
        }
        let m = lookup
            .get(row.input_id.as_str())
            .copied()
            .ok_or_else(|| bad("row refers to an unknown original input"))?;
        role_bases = role_bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= MAX_ROLE_BASES)
            .ok_or(Error::TooLarge)?;
        let units = match &row.preparation {
            Preparation::WholeLinear => {
                if m.topology != Topology::Linear {
                    return Err(bad(
                        "WholeLinear cannot materialize a circular source; specify actual PCR or digest preparation",
                    ));
                }
                0
            }
            Preparation::ExactPcr { forward, reverse } => {
                primer(forward)?;
                primer(reverse)?;
                m.sequence
                    .len()
                    .checked_mul(forward.annealing.len() + reverse.annealing.len())
                    .and_then(|n| n.checked_mul(4))
                    .ok_or(Error::TooLarge)?
            }
            Preparation::Digest {
                enzymes,
                fragment_index,
                ..
            } => {
                if enzymes.is_empty() || enzymes.len() > 32 || *fragment_index > 100_000 {
                    return Err(Error::TooLarge);
                }
                let mut names = BTreeSet::new();
                let mut width = 0usize;
                for name in enzymes {
                    token(name, 128)?;
                    if !names.insert(name) {
                        return Err(bad("duplicate digest enzyme"));
                    }
                    let enzyme=catalog.iter().find(|e|e.name.eq_ignore_ascii_case(name)).ok_or_else(||bad("digest enzyme is not supported by the existing typed Digest operation"))?;
                    width = width
                        .checked_add(enzyme.recognition.len() * 2)
                        .ok_or(Error::TooLarge)?;
                }
                m.sequence.len().checked_mul(width).ok_or(Error::TooLarge)?
            }
        };
        work = work
            .checked_add(units)
            .filter(|n| *n <= MAX_PREPARATION_WORK)
            .ok_or(Error::TooLarge)?;
    }
    hash(recipe)?;
    Ok(())
}
pub fn preflight(inputs: &[PlanInput], target: &Molecule, recipe: &Recipe) -> Result<(), Error> {
    source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        target,
        recipe,
    )?;
    if inputs.iter().any(|i| !hex(&i.revision)) {
        return Err(bad("source revisions must be explicit canonical hashes"));
    }
    Ok(())
}
fn bound_revisions(inputs: &[PlanInput], target: &Molecule) -> Result<String, Error> {
    let target_revision = hash(target)?;
    for input in inputs {
        let actual = hash(&input.molecule)?;
        if actual != input.revision {
            return Err(bad("original input revision changed"));
        }
    }
    Ok(target_revision)
}
pub fn preview_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<String, Error> {
    preflight(inputs, target, recipe)?;
    bound_revisions(inputs, target)?;
    hash(&(ENGINE_VERSION, "preview", inputs, target, recipe))
}
/// Pure compilation: a structurally valid plan is not a verified assembly.
pub fn compile(
    inputs: &[PlanInput],
    target: &Molecule,
    recipe: &Recipe,
) -> Result<Compiled, Error> {
    preflight(inputs, target, recipe)?;
    let target_revision = bound_revisions(inputs, target)?;
    let mut steps = Vec::new();
    let mut rows = Vec::with_capacity(recipe.fragments.len());
    let mut fragments = Vec::with_capacity(recipe.fragments.len());
    for (index, row) in recipe.fragments.iter().enumerate() {
        let original = inputs
            .iter()
            .find(|i| i.id == row.input_id)
            .ok_or_else(|| bad("missing original"))?;
        let input = Source::Input(row.input_id.clone());
        let operation = match &row.preparation {
            Preparation::WholeLinear => None,
            Preparation::Digest {
                enzymes,
                fragment_index,
                ..
            } => Some(Operation::Digest {
                template: input.clone(),
                enzymes: enzymes.clone(),
                fragment_index: *fragment_index,
            }),
            Preparation::ExactPcr { forward, reverse } => Some(Operation::Pcr {
                template: input.clone(),
                forward: forward.clone(),
                reverse: reverse.clone(),
            }),
        };
        let preparation_step = operation.as_ref().map(|_| steps.len());
        let mut prepared_source = if let Some(operation) = operation {
            let step = steps.len();
            steps.push(operation);
            Source::Step(step)
        } else {
            input
        };
        let repair_step = if let Preparation::Digest {
            repair: Some(method),
            ..
        } = &row.preparation
        {
            let index = steps.len();
            steps.push(Operation::RepairEnds {
                template: prepared_source,
                method: *method,
            });
            prepared_source = Source::Step(index);
            Some(index)
        } else {
            None
        };
        rows.push(PreparedRow {
            row_index: index,
            role: row.role.clone(),
            input_id: row.input_id.clone(),
            input_revision: original.revision.clone(),
            preparation_step,
            repair_step,
            prepared_source: prepared_source.clone(),
            reverse: row.reverse,
        });
        fragments.push(OrientedSource {
            source: prepared_source,
            reverse: row.reverse,
        });
    }
    let assembly_step = steps.len();
    steps.push(Operation::ProfiledOverlap {
        fragments,
        overlaps: recipe.overlaps.clone(),
        circular: recipe.circular,
        options: recipe.options.clone(),
    });
    let compiled = Compiled {
        plan: DesignPlan {
            name: recipe.name.clone(),
            target_revision,
            inputs: inputs.to_vec(),
            steps,
        },
        rows,
        assembly_step,
    };
    hash(&compiled)?;
    Ok(compiled)
}
pub fn preview(inputs: &[PlanInput], target: &Molecule, recipe: &Recipe) -> Result<Preview, Error> {
    let context_revision = preview_revision(inputs, target, recipe)?;
    let compiled = compile(inputs, target, recipe)?;
    // Exactly the existing verifier prepares PCR/digest products and projects
    // annotations. Neither the target nor a preview chooses a fragment or geometry.
    let mut verification = verification::verify(&compiled.plan, target);
    let mut available = Vec::with_capacity(compiled.rows.len());
    for row in &compiled.rows {
        let m = match &row.prepared_source {
            Source::Input(id) => compiled
                .plan
                .inputs
                .iter()
                .find(|i| &i.id == id)
                .map(|i| &i.molecule),
            Source::Step(i) => verification.products.get(*i),
        };
        if let Some(m) = m {
            available.push(m)
        } else {
            break;
        }
    }
    let preparation_complete = available.len() == compiled.rows.len();
    crate::model::bounds::model_preflight(&available)?;
    let mut bases = 0usize;
    for m in &available {
        bases = bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= MAX_ROLE_BASES)
            .ok_or(Error::TooLarge)?;
    }
    hash(&available)?;
    let prepared = available.into_iter().cloned().collect::<Vec<_>>();
    // An advisory junction report remains useful after a target mismatch. Physical
    // validity/blocked/staggered substrates remain the verifier's separate result.
    let profile = if preparation_complete {
        match overlap_profiles::check_oriented_fragments(
            &prepared,
            &compiled.rows.iter().map(|r| r.reverse).collect::<Vec<_>>(),
            &recipe.overlaps,
            recipe.circular,
            &recipe.options,
        ) {
            Ok(report) => Some(report),
            Err(error) => {
                verification.findings.push(format!(
                    "Recipe advisory junction profile withheld: {error}"
                ));
                None
            }
        }
    } else {
        None
    };
    let result = Preview {
        context_revision,
        compiled,
        verification,
        prepared,
        profile,
        preparation_complete,
    };
    result.validate()?;
    Ok(result)
}
impl Preview {
    /// Readonly shape/binding checks only. No reconstruction or authority restoration.
    pub fn validate(&self) -> Result<(), Error> {
        if !hex(&self.context_revision)
            || self.compiled.rows.is_empty()
            || self.compiled.rows.len() > MAX_ROWS
            || self.prepared.len() > self.compiled.rows.len()
            || self.compiled.plan.inputs.len() > MAX_ROWS
            || self.verification.products.len() > 2 * MAX_ROWS + 1
            || self.verification.step_reports.len() > 2 * MAX_ROWS + 1
            || self.verification.annotation_reports.len() > 2 * MAX_ROWS + 1
            || self.compiled.plan.steps.len() > 2 * MAX_ROWS + 1
            || self.compiled.assembly_step.checked_add(1) != Some(self.compiled.plan.steps.len())
        {
            return Err(bad("invalid preview shape"));
        }
        if self.preparation_complete != (self.prepared.len() == self.compiled.rows.len())
            || self.profile.is_some() && !self.preparation_complete
            || self.verification.target_revision != self.compiled.plan.target_revision
        {
            return Err(bad("preview plan/evidence binding mismatch"));
        }
        let models = self
            .compiled
            .plan
            .inputs
            .iter()
            .map(|i| &i.molecule)
            .chain(self.prepared.iter())
            .chain(self.verification.products.iter())
            .collect::<Vec<_>>();
        for model in models {
            crate::model::bounds::model_preflight(&[model])?;
        }
        for report in &self.verification.annotation_reports {
            if report.step_index >= self.compiled.plan.steps.len()
                || report.evidence.sources.len() > MAX_ROWS
            {
                return Err(Error::TooLarge);
            }
            let mut budget = crate::model::bounds::ModelBudget::default();
            budget.features(&report.evidence.features)?;
            for role in &report.evidence.sources {
                if role.unmapped.len() > 4096 {
                    return Err(Error::TooLarge);
                }
                for unmapped in &role.unmapped {
                    budget.features(std::slice::from_ref(&unmapped.feature))?;
                }
            }
        }
        if self.verification.plan_revision != hash(&self.compiled.plan)? {
            return Err(bad("preview plan revision changed"));
        }
        for (index, row) in self.compiled.rows.iter().enumerate() {
            token(&row.role, 128)?;
            token(&row.input_id, 128)?;
            if row.row_index != index || !hex(&row.input_revision) {
                return Err(bad("invalid preview row identity"));
            }
            let input = self
                .compiled
                .plan
                .inputs
                .iter()
                .find(|i| i.id == row.input_id)
                .ok_or_else(|| bad("preview row lost its source"))?;
            if input.revision != row.input_revision {
                return Err(bad("preview row revision mismatch"));
            }
            match (
                &row.prepared_source,
                row.repair_step.or(row.preparation_step),
            ) {
                (Source::Input(id), None) if id == &row.input_id => {}
                (Source::Step(i), Some(j)) if *i == j && j < self.compiled.assembly_step => {}
                _ => return Err(bad("invalid prepared-source mapping")),
            }
        }
        let Some(Operation::ProfiledOverlap {
            fragments,
            overlaps,
            circular,
            options,
        }) = self.compiled.plan.steps.last()
        else {
            return Err(bad("preview has no final profiled overlap operation"));
        };
        if fragments.len() != self.compiled.rows.len()
            || self.compiled.plan.steps[..self.compiled.assembly_step]
                .iter()
                .any(|step| {
                    !matches!(
                        step,
                        Operation::Pcr { .. }
                            | Operation::Digest { .. }
                            | Operation::RepairEnds { .. }
                    )
                })
        {
            return Err(bad("preview operation shape changed"));
        }
        for (fragment, row) in fragments.iter().zip(&self.compiled.rows) {
            let same = match (&fragment.source, &row.prepared_source) {
                (Source::Input(a), Source::Input(b)) => a == b,
                (Source::Step(a), Source::Step(b)) => a == b,
                _ => false,
            };
            if !same || fragment.reverse != row.reverse {
                return Err(bad("final overlap source mapping changed"));
            }
        }
        if let Some(profile) = &self.profile {
            profile.validate()?;
            if &profile.overlaps != overlaps
                || profile.circular != *circular
                || hash(&profile.options)? != hash(options)?
            {
                return Err(bad("profile recipe options changed"));
            }
            if profile.source_revisions
                != self
                    .prepared
                    .iter()
                    .map(hash)
                    .collect::<Result<Vec<_>, _>>()?
                || profile.orientations
                    != self
                        .compiled
                        .rows
                        .iter()
                        .map(|r| r.reverse)
                        .collect::<Vec<_>>()
            {
                return Err(bad("profile row bindings changed"));
            }
        }
        hash(self)?;
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad(
                "preview belongs to a different source/target/recipe context",
            ));
        }
        Ok(())
    }
}
