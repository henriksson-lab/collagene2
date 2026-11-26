//! Shared restrictions for reference-only exchange profiles which cannot carry
//! native physical evidence. This does not restrict reference alphabets/topology.
use crate::{Alphabet, Error, Molecule, chemistry};
pub(crate) fn reference_export_guard(m: &Molecule) -> Result<(), Error> {
    let producer = m.original_genbank.as_ref().is_some_and(|raw| {
        raw.lines()
            .filter_map(|s| {
                s.strip_prefix("COMMENT     ")
                    .or_else(|| s.strip_prefix("CC   "))
            })
            .any(|s| {
                matches!(
                    s.trim_end(),
                    "Collagene physical product v1: unverified sequence-only evidence"
                        | "Collagene reference-only v1: unverified sequence-only evidence"
                )
            })
    });
    let marked = m.metadata.keys().any(|key| {
        matches!(
            key.as_str(),
            "tracked_modifications"
                | "source:declared_termini"
                | "source:legacy:chemistry_unmapped"
                | "source:snapgene:chemistry_unmapped"
                | "topo:physical_graph_required"
                | "source:recombination:reference_only"
                | "source:recombination:unverified_recipe"
                | "source:restriction_inspection:recipe"
        ) || key == crate::recombination::REFERENCE_ONLY_KEY
            || key == crate::recombination::RECIPE_KEY
            || key == crate::restriction_inspection::RECIPE_KEY
            || key.starts_with("topo:")
            || key.starts_with("variant:")
            || key.starts_with("source:restriction_inspection:")
            || key.starts_with("restriction_catalog:")
            || key.starts_with("source:restriction_catalog:")
    });
    let protein_state = m.sequence.alphabet() == Alphabet::Protein
        && m.metadata.keys().any(|k| {
            k.eq_ignore_ascii_case("protein:modifications")
                || k.eq_ignore_ascii_case("protein:modified")
                || k.eq_ignore_ascii_case("modifications")
        });
    if producer || marked || protein_state || chemistry::has_unmapped_physical_evidence(m) {
        return Err(Error::Unsupported("Physical, damaged, single/mixed-strand, modified or reference-only evidence cannot be preserved by this reference export profile. Retain a native archive; accepting omissions cannot erase this restriction.".into()));
    }
    Ok(())
}
