//! Narrow host transition for applying whole-assembly primer suggestions.
//! This guard does not turn a suggestion into verified physical material.
use crate::assembly::draft::model::Draft;
use crate::cloning::preparation::cloning_methods_model::Mode;
use crate::cloning::preparation::cloning_methods_model::OverlapInput;
use crate::cloning::preparation::cloning_methods_model::OverlapPreparation;

fn exact(row: &OverlapInput) -> Result<(), String> {
    for (name, domain, tail) in [
        (&row.forward_name, &row.forward, &row.forward_tail),
        (&row.reverse_name, &row.reverse_primer, &row.reverse_tail),
    ] {
        if name.trim().is_empty() || name.len() > 256 || name.chars().any(char::is_control) {
            return Err("Each proposed PCR primer requires its own bounded name.".into());
        }
        let p = collagene_core::sequence::analysis::Primer::new(name.as_str(), domain, tail)
            .map_err(|e| e.to_string())?;
        if p.annealing.is_empty()
            || p.annealing.len() + p.tail.len() > collagene_core::cloning::pcr::primers::MAX_OLIGO
            || !p
                .annealing
                .as_str()
                .bytes()
                .chain(p.tail.bytes())
                .all(|b| b"ACGT".contains(&b))
        {
            return Err("Proposed oligos require nonempty concrete annealing DNA within the native length profile.".into());
        }
    }
    Ok(())
}

pub(crate) fn validate(before: &Draft, after: &Draft, current: &Draft) -> Result<usize, String> {
    before.validate()?;
    current.validate()?;
    after.validate()?;
    if current != before {
        return Err(
            "The assembly draft changed. Propose primers from its current declaration again."
                .into(),
        );
    }
    let n = before.slots.len();
    if before.methods.mode != Mode::Overlap || !(2..=12).contains(&n) {
        return Err("Whole-assembly primer proposals require an overlap construct with two to twelve source contributions.".into());
    }
    let roles = before.methods.role_ids(n)?;
    let overlaps = after
        .methods
        .overlaps
        .split(',')
        .map(|v| {
            v.trim()
                .parse::<usize>()
                .map_err(|_| "Invalid proposed homology length.".to_owned())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let count = n - 1 + usize::from(before.circular);
    if overlaps.len() != count
        || overlaps
            .iter()
            .any(|h| !(1..=80).contains(h) || *h != overlaps[0])
    {
        return Err(
            "The proposed assembly requires one explicit uniform homology length per junction."
                .into(),
        );
    }
    let mut expected = before.clone();
    expected.methods.overlaps = after.methods.overlaps.clone();
    for (i, role) in roles.iter().enumerate() {
        let original = &before.methods.rows[i];
        if !matches!(
            original.preparation,
            OverlapPreparation::WholeLinear
                | OverlapPreparation::ExactPcr
                | OverlapPreparation::MismatchPcr
                | OverlapPreparation::OverlapExtensionPcr
        ) || before.slots[i].anneal_bottom.is_some()
        {
            return Err("This proposal cannot replace an active digest, paired-oligo or declared-duplex route. Keep that preparation or explicitly choose a PCR route first.".into());
        }
        let edited = before.edits.iter().any(|e| e.role == role.role)
            || before.substitutions.iter().any(|e| e.role == role.role);
        let row = &after.methods.rows[i];
        let dest = &mut expected.methods.rows[i];
        dest.reverse = row.reverse;
        if edited {
            if row.preparation != OverlapPreparation::OverlapExtensionPcr {
                return Err("Internally edited contributions require their explicit overlap-extension PCR preparation.".into());
            }
            row.overlap_extension_pcr
                .as_ref()
                .ok_or("Proposed six-oligo preparation is missing.")?
                .request()?;
            dest.preparation = OverlapPreparation::OverlapExtensionPcr;
            dest.overlap_extension_pcr = row.overlap_extension_pcr.clone();
        } else {
            if row.preparation != OverlapPreparation::ExactPcr {
                return Err(
                    "An unchanged contribution requires the proposed exact-PCR pair.".into(),
                );
            }
            exact(row)?;
            dest.preparation = OverlapPreparation::ExactPcr;
            dest.template_strand = row.template_strand;
            dest.forward_name = row.forward_name.clone();
            dest.forward = row.forward.clone();
            dest.forward_tail = row.forward_tail.clone();
            dest.reverse_name = row.reverse_name.clone();
            dest.reverse_primer = row.reverse_primer.clone();
            dest.reverse_tail = row.reverse_tail.clone();
        }
    }
    expected.sync_rows();
    if expected.alternatives.is_some() {
        expected = crate::assembly::draft::assembly_alternatives::save_active(&expected)?;
    }
    if expected != *after {
        return Err("Primer suggestions may change only active PCR oligos and the requested overlap lengths. Sources, intent, saved origins and other settings must remain unchanged.".into());
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cloning::preparation::cloning_methods_model::OverlapExtensionPcrInput;

    fn transition() -> (Draft, Draft) {
        let mut before = Draft::default();
        before.methods.mode = Mode::Overlap;
        before.circular = false;
        let roles = before.methods.role_ids(2).unwrap();
        before
            .edits
            .push(collagene_core::cloning::assembly::edit_intent::Edit {
                role: roles[1].role.clone(),
                start: 16,
                end: 18,
                expected: "GG".into(),
                replacement: "TT".into(),
            });
        before.methods.source_roles = Some(roles);
        before.methods.rows[1].preparation = OverlapPreparation::WholeLinear;
        before.methods.rows[1].forward_tail = "inactive legacy tail".into();
        before.sync_rows();
        let mut after = before.clone();
        after.methods.overlaps = "12".into();
        let first = &mut after.methods.rows[0];
        first.preparation = OverlapPreparation::ExactPcr;
        first.forward_name = "Vector F".into();
        first.forward = "AGTCCGAT".into();
        first.reverse_name = "Vector R".into();
        first.reverse_primer = "TACGTAGC".into();
        first.reverse_tail = "ACGTACGTACGT".into();
        let mut raw = OverlapExtensionPcrInput::default();
        raw.left.forward.annealing = "AGTCCGAT".into();
        raw.left.reverse.annealing = "TACGTAGC".into();
        raw.right.forward.annealing = "TTCAGACC".into();
        raw.right.reverse.annealing = "TGACGATC".into();
        raw.outer.forward.annealing = "AGTCCGAT".into();
        raw.outer.reverse.annealing = "TGACGATC".into();
        after.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcr;
        after.methods.rows[1].overlap_extension_pcr = Some(raw);
        after.sync_rows();
        (before, after)
    }

    #[test]
    fn mixed_pair_and_six_oligos_change_only_active_fields_and_uniform_homology() {
        let (before, after) = transition();
        assert_eq!(validate(&before, &after, &before).unwrap(), 2);
        assert_eq!(after.methods.rows[1].forward_tail, "inactive legacy tail");
        assert!(after.slots == before.slots);
        assert_eq!(after.edits, before.edits);
        let mut without_edit_route = after.clone();
        without_edit_route.methods.rows[1].preparation = OverlapPreparation::ExactPcr;
        assert!(validate(&before, &without_edit_route, &before).is_err());
    }

    #[test]
    fn stale_drafts_changed_targets_and_inactive_fields_cannot_apply() {
        let (before, after) = transition();
        let mut current = before.clone();
        current.active_slot = 1;
        assert!(validate(&before, &after, &current).is_err());
        for field in 0..5 {
            let mut changed = after.clone();
            match field {
                0 => changed.name = "different intended target".into(),
                1 => changed.edits[0].replacement = "AA".into(),
                2 => changed.methods.rows[1].forward_tail.clear(),
                3 => changed.methods.comparisons = "10".into(),
                _ => changed.methods.overlaps = "12,12".into(),
            }
            assert!(validate(&before, &changed, &before).is_err());
        }
    }

    #[test]
    fn active_digest_or_declared_material_is_not_silently_replaced() {
        let (mut before, after) = transition();
        for prep in [
            OverlapPreparation::Digest,
            OverlapPreparation::CatalogDigest,
            OverlapPreparation::DeclaredDuplex,
        ] {
            before.methods.rows[0].preparation = prep;
            assert!(validate(&before, &after, &before).is_err());
        }
    }
}
