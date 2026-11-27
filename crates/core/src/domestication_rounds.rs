//! Sequential declarations from original-coordinate patches. Intended intermediate
//! frames supply primer literals only; physical PCR always consumes the preceding
//! actual plan step. No intermediate reference is installed as a template stock.
use crate::{
    Error, Sequence, Topology,
    chemistry::EndTreatment,
    domestication::Patch,
    domestication_preparation::{InverseRound, RoundEvidence},
    mutagenesis::{self, InversePcrOptions},
    verification::{Operation, OrientedSource, PlanInput, Source},
};
const MAX_ROUNDS: usize = 32;
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("sequential domestication: {s}"))
}
fn contains(round: &InverseRound, position: usize, n: usize) -> bool {
    let position = if position < round.edit_start {
        position + n
    } else {
        position
    };
    position >= round.edit_start && position < round.edit_end
}
pub(crate) fn preflight(
    input: &PlanInput,
    patches: &[Patch],
    rounds: &[InverseRound],
) -> Result<(), Error> {
    let n = input.molecule.sequence.len();
    if input.molecule.topology != Topology::Circular
        || n == 0
        || !(2..=MAX_ROUNDS).contains(&rounds.len())
        || patches.is_empty()
    {
        return Err(bad(
            "declare2..32 rounds over one circular original and actual selected substitutions",
        ));
    }
    let mut work = 0usize;
    let mut material = 0usize;
    for round in rounds {
        let span = round
            .edit_end
            .checked_sub(round.edit_start)
            .filter(|s| *s > 0 && *s < n)
            .ok_or_else(|| bad("each round replaces a nonempty proper original interval"))?;
        if round.edit_start >= n
            || round.edit_end > round.edit_start + n
            || !(8..=80).contains(&round.forward_length)
            || !(8..=80).contains(&round.reverse_length)
            || round.forward_length + round.reverse_length > n - span
            || round.reverse_tail_bases > span
            || span - round.reverse_tail_bases + round.forward_length > 200
            || round.reverse_tail_bases + round.reverse_length > 200
        {
            return Err(bad(
                "original-coordinate round or complete oligo lengths exceed supported bounds",
            ));
        }
        if !patches.iter().any(|p| contains(round, p.position, n)) {
            return Err(bad(
                "every round must carry at least one selected substitution",
            ));
        }
        work = work
            .checked_add(
                n.checked_mul(8 * (round.forward_length + round.reverse_length))
                    .ok_or(Error::TooLarge)?,
            )
            .filter(|w| *w <= 100_000_000)
            .ok_or(Error::TooLarge)?;
        material = material
            .checked_add(n.checked_mul(8).ok_or(Error::TooLarge)?)
            .filter(|m| *m <= 64 * 1024 * 1024)
            .ok_or(Error::TooLarge)?;
    }
    let mut seen = std::collections::BTreeSet::new();
    for patch in patches {
        if patch.position >= n
            || patch.expected.len() != 1
            || patch.replacement.len() != 1
            || !patch.expected.bytes().all(|b| b"ACGT".contains(&b))
            || !patch.replacement.bytes().all(|b| b"ACGT".contains(&b))
            || patch.expected == patch.replacement
            || input.molecule.sequence.as_str().as_bytes()[patch.position]
                != patch.expected.as_bytes()[0]
            || !seen.insert(patch.position)
            || rounds
                .iter()
                .filter(|r| contains(r, patch.position, n))
                .count()
                != 1
        {
            return Err(bad(
                "assign each distinct exact original substitution to exactly one round",
            ));
        }
    }
    Ok(())
}
pub(crate) fn compile(
    input: &PlanInput,
    patches: &[Patch],
    rounds: &[InverseRound],
    material: Option<&crate::run_material::MaterialOrigin>,
) -> Result<(Vec<Operation>, Vec<RoundEvidence>), Error> {
    preflight(input, patches, rounds)?;
    let n = input.molecule.sequence.len();
    let mut intended = input.molecule.sequence.as_str().as_bytes().to_vec();
    let mut origin = 0usize;
    let mut steps = Vec::with_capacity(3 * rounds.len() + usize::from(material.is_some()));
    let mut template =
        crate::domestication_preparation::append_origin(&input.id, material, &mut steps);
    let mut evidence = Vec::with_capacity(rounds.len());
    for (index, round) in rounds.iter().enumerate() {
        // This sequence is explicitly an intended coordinate frame for oligo
        // declarations. It is never a Molecule or a physical plan input.
        let text = (0..n)
            .map(|i| intended[(i + origin) % n] as char)
            .collect::<String>();
        let frame = Sequence::dna(&text)?;
        let edit_start = (round.edit_start + n - origin) % n;
        let span = round.edit_end - round.edit_start;
        let mut inserted = (0..span)
            .map(|i| intended[(round.edit_start + i) % n])
            .collect::<Vec<_>>();
        for patch in patches.iter().filter(|p| contains(round, p.position, n)) {
            let offset = (patch.position + n - round.edit_start) % n;
            inserted[offset] = patch.replacement.as_bytes()[0];
        }
        let options = InversePcrOptions {
            edit_start,
            edit_end: edit_start + span,
            inserted: String::from_utf8(inserted)
                .map_err(|_| bad("non-ASCII literal substitution"))?,
            forward_length: round.forward_length,
            reverse_length: round.reverse_length,
            reverse_tail_bases: round.reverse_tail_bases,
        };
        let parts = mutagenesis::inverse_parts(&frame, &options)?;
        let mut forward = parts.forward;
        let mut reverse = parts.reverse;
        forward.name = format!("Domestication round{} forward", index + 1);
        reverse.name = format!("Domestication round{} reverse", index + 1);
        let pcr = steps.len();
        steps.push(Operation::Pcr {
            template: template.clone(),
            forward: forward.clone(),
            reverse: reverse.clone(),
            template_strand: crate::chemistry::Strand::Top,
        });
        steps.push(Operation::TreatEnds {
            template: Source::Step(pcr),
            treatment: EndTreatment::Phosphorylate,
        });
        steps.push(Operation::Ligate {
            fragments: vec![OrientedSource {
                source: Source::Step(pcr + 1),
                reverse: false,
            }],
            circular: true,
            allow_nicked: false,
            options: Default::default(),
        });
        evidence.push(RoundEvidence {
            round: round.clone(),
            template: template.clone(),
            source_origin: origin,
            options,
            forward,
            reverse,
            forward_location: parts.forward_location,
            reverse_location: parts.reverse_location,
            retained_template_length: parts.retained,
        });
        template = Source::Step(pcr + 2);
        for patch in patches.iter().filter(|p| contains(round, p.position, n)) {
            intended[patch.position] = patch.replacement.as_bytes()[0];
        }
        // PCR starts with the forward5′ tail (insertion suffix), then retained
        // template and reverse-tail complement. Keep that actual product origin.
        origin = (round.edit_start + round.reverse_tail_bases) % n;
    }
    Ok((steps, evidence))
}
