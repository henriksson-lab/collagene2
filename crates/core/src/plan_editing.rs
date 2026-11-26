//! Atomic edits of a typed plan's dependency graph. Original step identities are
//! remapped explicitly; molecule content, targets and science are never evaluated.
use crate::{
    Error,
    verification::{Operation, Source, operation_sources, operation_sources_mut},
};
use std::collections::BTreeSet;
pub const MAX_INPUTS: usize = 1000;
pub const MAX_STEPS: usize = 1000;
pub const MAX_SOURCES_PER_STEP: usize = 1000;
fn bad(s: &str) -> Error {
    Error::Format(format!("Plan editing: {s}"))
}
fn identifiers(ids: &[String]) -> Result<BTreeSet<&str>, Error> {
    if ids.len() > MAX_INPUTS {
        return Err(Error::TooLarge);
    }
    let mut found = BTreeSet::new();
    for id in ids {
        if id.trim().is_empty() || id.len() > 128 || !found.insert(id.as_str()) {
            return Err(bad("input IDs must be nonempty, bounded and unique"));
        }
    }
    Ok(found)
}
// Bound the variable-size reference vectors before the shared visitor allocates
// its borrowed references. Other operation settings remain editable declarations.
fn reference_count(op: &Operation) -> usize {
    match op {
        Operation::Topo { sources, .. }
        | Operation::CheckedTopo { sources, .. }
        | Operation::ProcessedTopo { sources, .. }
        | Operation::TopoActions { sources, .. }
        | Operation::MultiSiteRecombineReference { sources, .. } => sources.len(),
        Operation::Overlap { fragments, .. }
        | Operation::GoldenGate {
            templates: fragments,
            ..
        }
        | Operation::SelectedGoldenGate {
            templates: fragments,
            ..
        } => fragments.len(),
        Operation::Ligate { fragments, .. }
        | Operation::LigateWithOutcomes { fragments, .. }
        | Operation::ProfiledOverlap { fragments, .. }
        | Operation::HomologyAssembly { fragments, .. }
        | Operation::HomologyAnnealing { fragments, .. } => fragments.len(),
        Operation::RecombineReference { .. }
        | Operation::AnnealOligos { .. }
        | Operation::OverlapExtension { .. } => 2,
        Operation::ResectThreePrime { .. }
        | Operation::ReconstructRunProduct { .. }
        | Operation::ReconstructDigestFragment { .. }
        | Operation::ReconstructMaterial { .. }
        | Operation::Pcr { .. }
        | Operation::PcrWithMismatches { .. }
        | Operation::Digest { .. }
        | Operation::CatalogDigest { .. }
        | Operation::CheckedDigest { .. }
        | Operation::DeclareSourceModifications { .. }
        | Operation::DeclareSourceTermini { .. }
        | Operation::TreatEnds { .. }
        | Operation::RepairEnds { .. }
        | Operation::Tail3Prime { .. } => 1,
    }
}
fn references<'a>(op: &'a Operation) -> Result<Vec<&'a Source>, Error> {
    if reference_count(op) > MAX_SOURCES_PER_STEP {
        return Err(Error::TooLarge);
    }
    Ok(operation_sources(op))
}
fn validate_operation(
    ids: &BTreeSet<&str>,
    operation: &Operation,
    index: usize,
) -> Result<(), Error> {
    for source in references(operation)? {
        match source {
            Source::Input(id) => {
                if id.trim().is_empty() || id.len() > 128 || !ids.contains(id.as_str()) {
                    return Err(bad("step references an unavailable original input ID"));
                }
            }
            Source::Step(prior) => {
                if *prior >= index {
                    return Err(bad(
                        "step dependencies must refer strictly to earlier steps",
                    ));
                }
            }
        }
    }
    Ok(())
}
/// Graph validation only. Empty editable plans are valid here; the scientific
/// verifier separately requires a complete, meaningful cloning plan.
pub fn validate_steps(input_ids: &[String], steps: &[Operation]) -> Result<(), Error> {
    if steps.len() > MAX_STEPS {
        return Err(Error::TooLarge);
    }
    let ids = identifiers(input_ids)?;
    for (index, operation) in steps.iter().enumerate() {
        validate_operation(&ids, operation, index)?;
    }
    Ok(())
}
/// `to` is the moved step's final zero-based index. Dependencies retain their
/// original producing-step identity; moving a consumer before its producer fails.
pub fn moved(
    input_ids: &[String],
    steps: &[Operation],
    from: usize,
    to: usize,
) -> Result<Vec<Operation>, Error> {
    validate_steps(input_ids, steps)?;
    if from >= steps.len() || to >= steps.len() {
        return Err(bad("move index is outside the step list"));
    }
    let mut order = (0..steps.len()).collect::<Vec<_>>();
    let old = order.remove(from);
    order.insert(to, old);
    let mut mapping = vec![0; steps.len()];
    for (new, old) in order.iter().copied().enumerate() {
        mapping[old] = new;
    }
    // Refuse before cloning any operation when the proposed order breaks a DAG edge.
    for (new, old) in order.iter().copied().enumerate() {
        for source in references(&steps[old])? {
            if let Source::Step(prior) = source {
                if mapping[*prior] >= new {
                    return Err(bad("proposed order places a step before its dependency"));
                }
            }
        }
    }
    let mut result = Vec::with_capacity(steps.len());
    for old in order {
        let mut operation = steps[old].clone();
        for source in operation_sources_mut(&mut operation) {
            if let Source::Step(prior) = source {
                *prior = mapping[*prior];
            }
        }
        result.push(operation);
    }
    Ok(result)
}
/// Refuses removal of a referenced producer. Remaining dependency identities are
/// retained with shifted indexes; dependent subgraphs are never silently deleted.
pub fn removed(
    input_ids: &[String],
    steps: &[Operation],
    index: usize,
) -> Result<Vec<Operation>, Error> {
    validate_steps(input_ids, steps)?;
    if index >= steps.len() {
        return Err(bad("remove index is outside the step list"));
    }
    for operation in steps {
        if references(operation)?
            .into_iter()
            .any(|s| matches!(s,Source::Step(i) if *i==index))
        {
            return Err(bad("step is still referenced by a later operation"));
        }
    }
    let mut result = Vec::with_capacity(steps.len() - 1);
    for (old, operation) in steps.iter().enumerate() {
        if old == index {
            continue;
        }
        let mut operation = operation.clone();
        for source in operation_sources_mut(&mut operation) {
            if let Source::Step(prior) = source {
                if *prior > index {
                    *prior -= 1;
                }
            }
        }
        result.push(operation);
    }
    Ok(result)
}
/// Existing steps use their old indexes. A new operation's Step references also
/// use PRE-insertion indexes and must already precede the insertion position.
pub fn inserted(
    input_ids: &[String],
    steps: &[Operation],
    index: usize,
    operation: Operation,
) -> Result<Vec<Operation>, Error> {
    validate_steps(input_ids, steps)?;
    if index > steps.len() || steps.len() >= MAX_STEPS {
        return Err(bad("insert position or step capacity is invalid"));
    }
    validate_operation(&identifiers(input_ids)?, &operation, index)?;
    let mut result = Vec::with_capacity(steps.len() + 1);
    let mut new = Some(operation);
    for position in 0..=steps.len() {
        if position == index {
            result.push(new.take().expect("one validated insertion position"));
        }
        if let Some(old) = steps.get(position) {
            let mut old = old.clone();
            for source in operation_sources_mut(&mut old) {
                if let Source::Step(prior) = source {
                    if *prior >= index {
                        *prior += 1;
                    }
                }
            }
            result.push(old);
        }
    }
    Ok(result)
}
/// Replacement retains this step's identity for later consumers. It does not
/// claim that changed preparation settings still produce a suitable molecule.
pub fn replaced(
    input_ids: &[String],
    steps: &[Operation],
    index: usize,
    operation: Operation,
) -> Result<Vec<Operation>, Error> {
    if steps.len() > MAX_STEPS || index >= steps.len() {
        return Err(bad("replace index or step capacity is invalid"));
    }
    let ids = identifiers(input_ids)?;
    // Validate the result before cloning. The replaced declaration can repair
    // its own stale/invalid dependency, while every other node stays intact.
    for (i, existing) in steps.iter().enumerate() {
        validate_operation(&ids, if i == index { &operation } else { existing }, i)?;
    }
    let mut result = steps.to_vec();
    result[index] = operation;
    Ok(result)
}
