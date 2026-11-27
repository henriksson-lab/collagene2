//! Guard the narrow draft transition which installs suggested preparation oligos.
//! No source, edit intent, historical origin or target is replaced here.
use crate::assembly_draft_model::Draft;
#[cfg(test)]
use crate::cloning_methods_model::OverlapPreparation;

pub(crate) fn validate(before: &Draft, after: &Draft, current: &Draft) -> Result<usize, String> {
    before.validate()?;
    current.validate()?;
    after.validate()?;
    if current != before {
        return Err("The draft changed while primer proposals were being reviewed. Propose from the current draft again.".into());
    }
    let index = before.active_slot;
    let next = after
        .methods
        .rows
        .get(index)
        .ok_or("The proposed fragment is unavailable.")?;
    let original = before
        .methods
        .rows
        .get(index)
        .ok_or("The original fragment is unavailable.")?;
    let preparation =
        crate::assembly_extension_primer_model::proposed_preparation(original.preparation);
    if next.preparation != preparation {
        return Err(
            "The proposal must install the explicit overlap-extension PCR preparation.".into(),
        );
    }
    next.overlap_extension_pcr
        .as_ref()
        .ok_or("The six proposed oligos are missing.")?
        .request()?;
    // Captured source molecules and origin capsules are shared immutable values;
    // the bounded copy preserves all original event identities and declarations.
    let mut expected = before.clone();
    let row = expected
        .methods
        .rows
        .get_mut(index)
        .ok_or("The original fragment is unavailable.")?;
    row.preparation = preparation;
    row.overlap_extension_pcr = next.overlap_extension_pcr.clone();
    expected.sync_rows();
    if expected != *after {
        return Err("A primer proposal may change only this fragment's preparation oligos and their draft schema. Original sources, edits and other settings must remain unchanged.".into());
    }
    Ok(index)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cloning_methods_model::OverlapExtensionPcrInput;

    fn transition() -> (Draft, Draft) {
        let mut before = Draft::default();
        before.active_slot = 1;
        before.methods.rows[1].forward_tail = "ACGT".into();
        let mut raw = OverlapExtensionPcrInput::default();
        raw.left.forward.annealing = "AGTCCGAT".into();
        raw.left.reverse.annealing = "TACGTAGC".into();
        raw.left.reverse.tail = "GGTCTGAAAA".into();
        raw.right.forward.annealing = "TTCAGACC".into();
        raw.right.forward.tail = "GCTACGTATT".into();
        raw.right.reverse.annealing = "TGACGATC".into();
        raw.outer.forward.annealing = "AGTCCGAT".into();
        raw.outer.reverse.annealing = "TGACGATC".into();
        raw.overlap = "18".into();
        raw.min_paired_bases = "8".into();
        let mut after = before.clone();
        after.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcr;
        after.methods.rows[1].overlap_extension_pcr = Some(raw);
        after.sync_rows();
        (before, after)
    }

    #[test]
    fn suggested_oligos_preserve_independent_draft_and_inactive_settings() {
        let (before, after) = transition();
        assert_eq!(validate(&before, &after, &before).unwrap(), 1);
        assert_eq!(after.version, crate::assembly_draft_model::CURRENT_VERSION);
        assert_eq!(after.methods.rows[1].forward_tail, "ACGT");
        assert!(after.slots == before.slots);
        assert!(after.methods.rows[0] == before.methods.rows[0]);
    }

    #[test]
    fn upstream_suggestions_preserve_pinned_downstream_digest_and_reject_dropped_or_changed_selection()
     {
        for current in [
            OverlapPreparation::CatalogDigest,
            OverlapPreparation::PcrCatalogDigest,
            OverlapPreparation::MismatchPcrCatalogDigest,
            OverlapPreparation::OverlapExtensionPcrCatalogDigest,
        ] {
            let (mut before, mut after) = transition();
            let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
            let id = catalog
                .edition()
                .enzymes
                .iter()
                .find(|e| e.name == "EcoRI")
                .unwrap()
                .id
                .clone();
            let pin = collagene_core::catalog_digest::select(&catalog, &[id]).unwrap();
            before.methods.rows[1].preparation = current;
            before.methods.rows[1].catalog_selection = Some(pin.clone());
            before.methods.rows[1].fragment = "3".into();
            before.sync_rows();
            after.methods.rows[1].preparation =
                OverlapPreparation::OverlapExtensionPcrCatalogDigest;
            after.methods.rows[1].catalog_selection = Some(pin);
            after.methods.rows[1].fragment = "3".into();
            after.sync_rows();
            assert_eq!(validate(&before, &after, &before).unwrap(), 1);
            assert_eq!(after.version, crate::assembly_draft_model::CURRENT_VERSION);
            let mut dropped = after.clone();
            dropped.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcr;
            assert!(validate(&before, &dropped, &before).is_err());
            let mut changed = after.clone();
            changed.methods.rows[1].fragment = "2".into();
            assert!(validate(&before, &changed, &before).is_err());
            let mut changed = after.clone();
            changed.methods.rows[1].catalog_selection = None;
            assert!(validate(&before, &changed, &before).is_err());
        }
    }

    #[test]
    fn stale_or_expanded_transitions_cannot_replace_intent_other_rows_or_settings() {
        let (before, after) = transition();
        let mut current = before.clone();
        current.active_slot = 0;
        assert!(validate(&before, &after, &current).is_err());
        let mut changed = after.clone();
        changed.name = "Invented replacement target".into();
        assert!(validate(&before, &changed, &before).is_err());
        let mut changed = after.clone();
        changed.methods.rows[0].forward_tail = "TTTT".into();
        assert!(validate(&before, &changed, &before).is_err());
        let mut changed = after.clone();
        changed.methods.rows[1].forward_tail.clear();
        assert!(validate(&before, &changed, &before).is_err());
        let mut changed = after.clone();
        changed.methods.rows[1]
            .overlap_extension_pcr
            .as_mut()
            .unwrap()
            .outer
            .forward
            .annealing
            .clear();
        assert!(validate(&before, &changed, &before).is_err());
    }
}
