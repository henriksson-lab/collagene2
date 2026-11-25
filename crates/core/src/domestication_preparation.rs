//! Original-template physical preparation for explicitly checked domestication.
//! A reference-only recoded target never becomes an invented PCR template.
use crate::{
    Error, Molecule, Topology, assembly_recipe, domestication,
    mutagenesis::{self, InversePcrOptions, InversePcrPreparation},
    verification::{self, Operation, PlanInput, Source},
};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io};

pub const ENGINE_VERSION: &str = "domestication-preparation-1";
pub const INPUT_ID: &str = "template";
const MAX_SOURCE_BASES: usize = 8 << 20;
const MAX_WIRE_BYTES: usize = 32 << 20;
macro_rules! dto {($($tokens:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($tokens)*};}
dto! {#[derive(Copy)] pub enum Route {
    InversePcr {edit_start:usize,edit_end:usize,forward_length:usize,reverse_length:usize,reverse_tail_bases:usize},
}}
dto! {pub struct Request {pub design:domestication::Request,pub selection:domestication::Selection,pub route:Route}}
dto! {pub struct Prepared {
    pub context_revision:String,pub request:Request,pub source_revision:String,pub target:Molecule,
    pub check:domestication::CheckReport,pub preparation:InversePcrPreparation,
    pub verification:verification::VerificationReport,pub warnings:Vec<String>,
}}
fn bad(message: &str) -> Error {
    Error::Unsupported(format!("domestication preparation: {message}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}
struct Bytes(usize);
impl io::Write for Bytes {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(bytes.len())
            .filter(|size| *size <= MAX_WIRE_BYTES)
            .ok_or_else(|| io::Error::other("domestication preparation wire budget"))?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn dimensions(
    input: &PlanInput,
    request: &Request,
) -> Result<(usize, usize, usize, usize, usize), Error> {
    let Route::InversePcr {
        edit_start,
        edit_end,
        forward_length,
        reverse_length,
        reverse_tail_bases,
    } = request.route;
    let n = input.molecule.sequence.len();
    if n == 0 || n > MAX_SOURCE_BASES || input.molecule.topology != Topology::Circular {
        return Err(bad(
            "this explicit physical route requires a bounded circular DNA original; other routes remain unsupported",
        ));
    }
    let span = edit_end
        .checked_sub(edit_start)
        .filter(|span| *span > 0 && *span < n)
        .ok_or_else(|| bad("replace a nonempty interval shorter than one original turn"))?;
    if edit_start >= n
        || edit_start.checked_add(n).is_none_or(|end| edit_end > end)
        || !(8..=80).contains(&forward_length)
        || !(8..=80).contains(&reverse_length)
        || forward_length + reverse_length > n - span
    {
        return Err(bad(
            "choose original-source coordinates and nonoverlapping8–80 base annealing domains",
        ));
    }
    if reverse_tail_bases > span
        || reverse_tail_bases > 200 - reverse_length
        || span - reverse_tail_bases > 200 - forward_length
    {
        return Err(bad(
            "every substituted base must fit the explicitly split actual primer tails, with each complete ordered primer at most200 bases; distant edits require another physical route or reviewed rounds",
        ));
    }
    if request.selection.patches.is_empty() {
        return Err(bad(
            "no substitutions were selected for a physical preparation",
        ));
    }
    let mut seen = BTreeSet::new();
    for patch in &request.selection.patches {
        if patch.position >= n {
            return Err(bad(
                "a selected substitution lies outside the original source",
            ));
        }
        let offset = if patch.position >= edit_start {
            patch.position - edit_start
        } else {
            patch.position + n - edit_start
        };
        if patch.position >= n || offset >= span || !seen.insert(patch.position) {
            return Err(bad(
                "each selected substitution must occur exactly once inside the explicit replacement interval",
            ));
        }
        if patch.expected.len() != 1
            || patch.replacement.len() != 1
            || !patch.expected.bytes().all(|base| b"ACGT".contains(&base))
            || !patch
                .replacement
                .bytes()
                .all(|base| b"ACGT".contains(&base))
            || patch.expected == patch.replacement
            || input.molecule.sequence.as_str().as_bytes()[patch.position]
                != patch.expected.as_bytes()[0]
        {
            return Err(bad(
                "selected substitutions must name the actual concrete original base and a distinct concrete replacement",
            ));
        }
    }
    Ok((
        edit_start,
        edit_end,
        forward_length,
        reverse_length,
        reverse_tail_bases,
    ))
}
/// Cheap bounded identity/geometry guards; no genetic checking or PCR is replayed.
pub fn preflight(input: &PlanInput, request: &Request) -> Result<(), Error> {
    if input.id != INPUT_ID
        || input.revision != request.design.source_revision
        || !digest(&input.revision)
    {
        return Err(bad(
            "the explicit template input and design must bind the same original source revision",
        ));
    }
    crate::model::bounds::model_preflight(&[&input.molecule])?;
    if input.molecule.sequence.len() > MAX_SOURCE_BASES {
        return Err(Error::TooLarge);
    }
    input.molecule.validate()?;
    if input.molecule.revision() != input.revision {
        return Err(bad("original template content changed"));
    }
    // Native domestication preflight owns catalog, feature, patch and search bounds.
    domestication::check_revision(&input.molecule, &request.design, &request.selection)?;
    dimensions(input, request)?;
    serde_json::to_writer(&mut Bytes(0), &(input, request)).map_err(|_| Error::TooLarge)?;
    Ok(())
}
pub fn revision(input: &PlanInput, request: &Request) -> Result<String, Error> {
    preflight(input, request)?;
    assembly_recipe::hash(&(ENGINE_VERSION, input, request))
}
fn options(input: &PlanInput, request: &Request) -> Result<InversePcrOptions, Error> {
    let (edit_start, edit_end, forward_length, reverse_length, reverse_tail_bases) =
        dimensions(input, request)?;
    let n = input.molecule.sequence.len();
    let original = input.molecule.sequence.as_str().as_bytes();
    let mut inserted = (edit_start..edit_end)
        .map(|position| original[position % n])
        .collect::<Vec<_>>();
    for patch in &request.selection.patches {
        let offset = if patch.position >= edit_start {
            patch.position - edit_start
        } else {
            patch.position + n - edit_start
        };
        inserted[offset] = patch.replacement.as_bytes()[0];
    }
    Ok(InversePcrOptions {
        edit_start,
        edit_end,
        inserted: String::from_utf8(inserted).map_err(|_| bad("non-ASCII original DNA"))?,
        forward_length,
        reverse_length,
        reverse_tail_bases,
    })
}
fn primer_region(start: usize, count: usize, n: usize, reverse: bool) -> crate::Location {
    let first = count.min(n - start);
    let mut parts = vec![crate::Location::range(start, start + first)];
    if first < count {
        parts.push(crate::Location::range(0, count - first));
    }
    let region = if parts.len() == 1 {
        parts.remove(0)
    } else {
        crate::Location::Join(parts)
    };
    if reverse {
        crate::Location::Complement(Box::new(region))
    } else {
        region
    }
}
/// Recheck literal substitutions, freeze the independent target, then compile and
/// verify actual original-template PCR → phosphorylation → sealed circular ligation.
pub fn materialize(input: &PlanInput, request: &Request) -> Result<Prepared, Error> {
    let context_revision = revision(input, request)?;
    let checked = domestication::check(&input.molecule, &request.design, &request.selection)?;
    let target = checked.target().clone(); // Frozen before any molecular reconstruction.
    let check = checked.report().clone();
    let route = options(input, request)?;
    let preparation = mutagenesis::compile_inverse_pcr(&input.molecule, &target, &route)?;
    if preparation.plan.inputs.len() != 1
        || preparation.plan.inputs[0].id != INPUT_ID
        || preparation.plan.inputs[0].revision != input.revision
    {
        return Err(bad(
            "the physical compiler did not retain the explicit original template",
        ));
    }
    let verification = verification::verify(&preparation.plan, &target);
    let prepared=Prepared{context_revision:context_revision.clone(),request:request.clone(),source_revision:input.revision.clone(),target,check,preparation,verification,
        warnings:vec![
            "Every selected substitution is supplied by actual primer tails replacing one explicit original-template interval. Unchanged intervening bases are also resupplied by those tails; no recoded reference is used as material.".into(),
            "The plan models exact PCR, terminal phosphorylation and sealed circular ligation. It does not model template clearance, enzyme kinetics, transformation, cellular repair or experimental success.".into(),
            "This single-interval route is bounded by complete200-base oligos and nonoverlapping annealing domains. Unsupported scattered edits remain reference-only until another explicit route or reviewed multi-round plan is supplied.".into(),
        ]};
    serde_json::to_writer(&mut Bytes(0), &prepared).map_err(|_| Error::TooLarge)?;
    prepared.validate_binding(&context_revision)?;
    Ok(prepared)
}
impl Prepared {
    /// Standalone wire/binding guard, not a replay or a successful science badge.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.context_revision != expected
            || !digest(expected)
            || !digest(&self.source_revision)
            || self.preparation.plan.inputs.len() != 1
        {
            return Err(bad("prepared source/context is malformed"));
        }
        let input = &self.preparation.plan.inputs[0];
        if revision(input, &self.request)? != expected
            || input.revision != self.source_revision
            || self.preparation.source_revision != self.source_revision
            || self.target.revision() != self.preparation.target_revision
            || self.preparation.plan.target_revision != self.preparation.target_revision
        {
            return Err(bad(
                "prepared source, request or independent target binding changed",
            ));
        }
        self.target.validate()?;
        if self.target.topology != Topology::Circular
            || self.target.sequence.alphabet() != crate::Alphabet::Dna
            || self.target.sequence.len() != input.molecule.sequence.len()
        {
            return Err(bad(
                "the independent target must retain the original substitution-only reference geometry",
            ));
        }
        let mut patches = self.request.selection.patches.iter().peekable();
        for (position, (original, target)) in input
            .molecule
            .sequence
            .as_str()
            .bytes()
            .zip(self.target.sequence.as_str().bytes())
            .enumerate()
        {
            let expected = if patches
                .peek()
                .is_some_and(|patch| patch.position == position)
            {
                patches
                    .next()
                    .expect("peeked literal substitution")
                    .replacement
                    .as_bytes()[0]
            } else {
                original
            };
            if target != expected {
                return Err(bad(
                    "the independent target contains letters outside the exact original-source substitutions",
                ));
            }
        }
        let original_check = domestication::check_revision(
            &input.molecule,
            &self.request.design,
            &self.request.selection,
        )?;
        self.check.validate_binding(&original_check)?;
        if self.check.source_revision != self.source_revision
            || self.check.target_revision != self.target.revision()
            || self.check.selection.candidate_id != self.request.selection.candidate_id
        {
            return Err(bad(
                "checked substitutions belong to another original or target",
            ));
        }
        let desired = options(input, &self.request)?;
        let actual = &self.preparation.options;
        if actual.edit_start != desired.edit_start
            || actual.edit_end != desired.edit_end
            || actual.inserted != desired.inserted
            || actual.forward_length != desired.forward_length
            || actual.reverse_length != desired.reverse_length
            || actual.reverse_tail_bases != desired.reverse_tail_bases
        {
            return Err(bad(
                "physical preparation route no longer encodes the exact substitutions",
            ));
        }
        let n = input.molecule.sequence.len();
        let forward_location =
            primer_region(desired.edit_end % n, desired.forward_length, n, false);
        let reverse_location = primer_region(
            (desired.edit_start + n - desired.reverse_length) % n,
            desired.reverse_length,
            n,
            true,
        );
        let forward_domain = forward_location.extract(&input.molecule.sequence)?;
        let reverse_domain = reverse_location.extract(&input.molecule.sequence)?;
        let reverse_tail = if desired.reverse_tail_bases == 0 {
            String::new()
        } else {
            crate::Sequence::dna(&desired.inserted[..desired.reverse_tail_bases])?
                .reverse_complement()?
                .as_str()
                .to_string()
        };
        if self.preparation.forward_location != forward_location
            || self.preparation.reverse_location != reverse_location
            || self.preparation.retained_template_length
                != n - (desired.edit_end - desired.edit_start)
            || self.preparation.forward.name != "Mutagenesis forward"
            || self.preparation.reverse.name != "Mutagenesis reverse"
            || self.preparation.forward.annealing.as_str() != forward_domain
            || self.preparation.reverse.annealing.as_str() != reverse_domain
            || self.preparation.forward.tail != desired.inserted[desired.reverse_tail_bases..]
            || self.preparation.reverse.tail != reverse_tail
            || self.preparation.plan.name
                != format!("Inverse-PCR replacement in {}", input.molecule.name)
        {
            return Err(bad(
                "actual primer domains, locations, tails or retained-template geometry differ from the explicit original-source route",
            ));
        }
        let steps = &self.preparation.plan.steps;
        if steps.len() != 3 {
            return Err(bad(
                "inverse-PCR preparation must retain all three actual operations",
            ));
        }
        if !matches!(&steps[0],Operation::Pcr{template:Source::Input(id),forward,reverse}
            if id==INPUT_ID
                && forward.name==self.preparation.forward.name
                && forward.annealing==self.preparation.forward.annealing
                && forward.tail==self.preparation.forward.tail
                && reverse.name==self.preparation.reverse.name
                && reverse.annealing==self.preparation.reverse.annealing
                && reverse.tail==self.preparation.reverse.tail)
            || !matches!(
                &steps[1],
                Operation::TreatEnds {
                    template: Source::Step(0),
                    treatment: crate::chemistry::EndTreatment::Phosphorylate
                }
            )
            || !matches!(&steps[2],Operation::Ligate{fragments,circular:true,allow_nicked:false}
                if fragments.len()==1&&matches!(fragments[0].source,Source::Step(1))&&!fragments[0].reverse)
        {
            return Err(bad("inverse-PCR physical operation chain changed"));
        }
        for primer in [&self.preparation.forward, &self.preparation.reverse] {
            if primer.ordered_sequence().len() > 200 {
                return Err(bad("complete physical oligo exceeds200 bases"));
            }
        }
        if self.verification.plan_revision != verification::plan_revision(&self.preparation.plan)
            || self.verification.target_revision != self.target.revision()
            || self.verification.products.len() > 3
        {
            return Err(bad(
                "verification describes a different plan or independent target",
            ));
        }
        if self.warnings.len() > 64 || self.warnings.iter().any(|warning| warning.len() > 65536) {
            return Err(Error::TooLarge);
        }
        serde_json::to_writer(&mut Bytes(0), self).map_err(|_| Error::TooLarge)?;
        Ok(())
    }
}
