//! Source/event-bound adapter for native overlap-extension oligo proposals.
//! Applying a proposal changes one preparation declaration, never design intent
//! or source documents. Ordinary recipe verification remains authoritative.
use crate::assembly::draft::model::Draft;
use crate::assembly::review::assembly_review::IntendedReference;
use crate::cloning::preparation::cloning_methods_model::ExtensionPrimerInput;
use crate::cloning::preparation::cloning_methods_model::ExtensionPrimerPairInput;
use crate::cloning::preparation::cloning_methods_model::OverlapExtensionPcrInput;
use crate::cloning::preparation::cloning_methods_model::OverlapPreparation;
use collagene_core::{
    cloning::{
        assembly::edit_intent as assembly_edit_intent, overlap::extension::primers as native,
    },
    documents::workspace::{RevisionEvent, Workspace},
    material::lineage::EventBinding,
    sequence::analysis::Primer,
    verification::evidence::VerificationReport,
};
use collagene_model::{Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::{io::Write, sync::Arc};
const ENGINE_VERSION: &str = "assembly-extension-primer-proposal-2";
const MAX_REPORT: usize = 32 << 20;
const MAX_CONTEXT: usize = 128 << 20;
struct Digest {
    hash: blake3::Hasher,
    bytes: usize,
    max: usize,
}
impl Write for Digest {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(b.len())
            .filter(|n| *n <= self.max)
            .ok_or_else(|| std::io::Error::other("overlap-extension proposal bound"))?;
        self.hash.update(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn receipt(value: &impl Serialize, max: usize) -> Result<String, String> {
    let mut writer = Digest {
        hash: blake3::Hasher::new(),
        bytes: 0,
        max,
    };
    serde_json::to_writer(&mut writer, value).map_err(|e| e.to_string())?;
    Ok(writer.hash.finalize().to_hex().to_string())
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Capture {
    pub context_revision: String,
    pub selected_index: usize,
    pub source: EventBinding,
    pub role: String,
    pub input_id: String,
    pub intent_context_revision: String,
    pub target_revision: String,
    pub draft_revision: String,
    pub options_revision: String,
    pub design_revision: String,
}
#[derive(Clone)]
pub struct Proposal {
    pub capture: Capture,
    pub preparation: collagene_core::cloning::overlap::extension::pcr::Request,
    pub verification: VerificationReport,
    pub warnings: Vec<String>,
    worker: Arc<native::Report>,
    worker_receipt: String,
    public_receipt: String,
}
/// Keep a declared downstream digestion when replacing the upstream PCR route.
/// Catalog pins and fragment numbers remain untouched and require fresh inspection.
pub(crate) fn proposed_preparation(current: OverlapPreparation) -> OverlapPreparation {
    match current {
        OverlapPreparation::CatalogDigest
        | OverlapPreparation::PcrCatalogDigest
        | OverlapPreparation::MismatchPcrCatalogDigest
        | OverlapPreparation::OverlapExtensionPcrCatalogDigest => {
            OverlapPreparation::OverlapExtensionPcrCatalogDigest
        }
        _ => OverlapPreparation::OverlapExtensionPcr,
    }
}

/// Retain the active preparation's outer-oligo sequences using the same case
/// and whitespace normalization as ordinary Primer::new. A blank OE editor has
/// no explicit tails and does not require its inner oligos to be completed.
pub fn current_outer_tails(draft: &Draft, index: usize) -> Result<(String, String), String> {
    let row = draft
        .methods
        .rows
        .get(index)
        .ok_or("Select a source contribution.")?;
    let (forward, reverse) = if matches!(
        row.preparation,
        OverlapPreparation::OverlapExtensionPcr
            | OverlapPreparation::OverlapExtensionPcrCatalogDigest
    ) {
        let Some(raw) = row.overlap_extension_pcr.as_ref() else {
            return Ok((String::new(), String::new()));
        };
        raw.validate_raw()?;
        (
            raw.outer.forward.tail.as_str(),
            raw.outer.reverse.tail.as_str(),
        )
    } else {
        (row.forward_tail.as_str(), row.reverse_tail.as_str())
    };
    let normalize = |tail: &str| -> Result<String, String> {
        if tail.len() > 1000 {
            return Err("Existing PCR tails exceed the supported bound.".into());
        }
        let normalized = Primer::new("Retained outer tail", "A", tail)
            .map_err(|e| e.to_string())?
            .tail;
        if !normalized.bytes().all(|base| b"ACGT".contains(&base)) {
            return Err("Outer-PCR tails require concrete DNA; complete or clear the partial tail explicitly.".into());
        }
        Ok(normalized)
    };
    Ok((normalize(forward)?, normalize(reverse)?))
}
/// Attest exact retained events using borrowed history/model records. Equal DNA
/// in another document, a different event, or producer metadata is insufficient.
pub fn attest(draft: &Draft, index: usize, frame: &Workspace) -> Result<EventBinding, String> {
    draft.validate()?;
    if frame.documents.len() > 1000
        || frame.history.len() > collagene_core::documents::workspace::MAX_HISTORY_EVENTS
    {
        return Err("Workspace source history exceeds its supported bounds.".into());
    }
    for slot in &draft.slots {
        let primary = slot
            .captured
            .as_ref()
            .ok_or("Pin all original source contributions first.")?;
        for captured in std::iter::once(primary).chain(slot.anneal_bottom.iter()) {
            let pin = &captured.source.binding;
            if frame
                .documents
                .iter()
                .filter(|d| d.id == pin.document_id)
                .count()
                != 1
            {
                return Err("The exact source document is unavailable or duplicated.".into());
            }
            let mut events = frame.history.iter().filter(|e| e.id == pin.event_id);
            let event = events
                .next()
                .ok_or("The exact historical source event is unavailable.")?;
            if events.next().is_some()
                || event.document_id != pin.document_id
                || event.revision != pin.revision
                || event.parents.len() > 128
                || event.operation.len() > 4096
                || event.parents.iter().any(|id| id.len() > 80)
            {
                return Err("The historical source identity changed.".into());
            }
            let canonical = RevisionEvent::new(
                event.document_id.clone(),
                event.revision.clone(),
                event.parents.clone(),
                event.operation.clone(),
            )
            .map_err(|e| e.to_string())?;
            let molecule = frame
                .revisions
                .get(&pin.revision)
                .ok_or("The original source revision is unavailable.")?;
            if canonical.id != event.id
                || molecule != &captured.molecule
                || molecule.revision() != pin.revision
            {
                return Err("The retained source event differs from its pinned molecule.".into());
            }
        }
    }
    draft
        .slots
        .get(index)
        .and_then(|s| s.captured.as_ref())
        .map(|c| c.source.binding.clone())
        .ok_or_else(|| "Select a pinned original source contribution.".into())
}
fn independent_intent(
    draft: &Draft,
    intent: &IntendedReference,
    inputs: &[collagene_core::planning::plan::PlanInput],
    base: collagene_core::cloning::assembly::intent::Request,
    workspace: Option<&Workspace>,
) -> Result<String, String> {
    let declared = draft.edit_request(base);
    match intent {
        IntendedReference::Edits(report) => {
            if declared.edits.is_empty() {
                return Err("Declare an internal contribution edit first.".into());
            }
            let expected = assembly_edit_intent::intent_revision(inputs, &declared, workspace)
                .map_err(|e| e.to_string())?;
            report
                .validate_sources(inputs, &expected, workspace)
                .map_err(|e| e.to_string())?;
            if receipt(&report.request, MAX_CONTEXT)? != receipt(&declared, MAX_CONTEXT)? {
                return Err("The frozen edit declaration changed.".into());
            }
            Ok(expected)
        }
        IntendedReference::Substitutions(report) => {
            if !declared.edits.is_empty() || declared.base.substitutions.is_empty() {
                return Err("The frozen substitution declaration differs from this draft.".into());
            }
            let expected = collagene_core::cloning::assembly::substitution_intent::intent_revision(
                inputs,
                &declared.base,
                workspace,
            )
            .map_err(|e| e.to_string())?;
            report
                .validate_sources(inputs, &expected, workspace)
                .map_err(|e| e.to_string())?;
            if receipt(&report.request, MAX_CONTEXT)? != receipt(&declared.base, MAX_CONTEXT)? {
                return Err("The frozen substitution declaration changed.".into());
            }
            Ok(expected)
        }
        _ => Err(
            "Declare internal edits or substitutions before proposing overlap-extension primers."
                .into(),
        ),
    }
}
pub fn request(
    draft: &Draft,
    intent: &IntendedReference,
    index: usize,
    current_source: &EventBinding,
    options: &native::Options,
    workspace: Option<&Workspace>,
) -> Result<(Capture, Molecule, native::Request), String> {
    draft.validate()?;
    if draft.alternatives.is_some() {
        return Err("Choose an ordinary construct before preparing an edited contribution.".into());
    }
    if draft.methods.mode
        == crate::cloning::preparation::cloning_methods_model::Mode::MultiSiteReference
    {
        return Err("Multi-parent reference recombination uses whole original parents; it does not consume row PCR preparations.".into());
    }
    if index != draft.active_slot {
        return Err("The selected source contribution changed.".into());
    }
    let slot = draft
        .slots
        .get(index)
        .ok_or("Select a source contribution.")?;
    let row = draft
        .methods
        .rows
        .get(index)
        .ok_or("The selected preparation is unavailable.")?;
    if row.preparation == OverlapPreparation::Digest {
        return Err("Pin the digest enzymes using the catalog before proposing an upstream PCR preparation; manual starter-enzyme digestion cannot be silently discarded.".into());
    }
    if row.preparation == OverlapPreparation::AnnealedOligos || slot.anneal_bottom.is_some() {
        return Err(
            "A paired-oligo contribution requires its existing two-source preparation.".into(),
        );
    }
    let captured = slot
        .captured
        .as_ref()
        .ok_or("Pin the original source contribution first.")?;
    if captured.source.binding != *current_source
        || captured.molecule.revision() != current_source.revision
    {
        return Err("The current source event differs from the pinned contribution.".into());
    }
    if !options.primer_options.forward_tail.is_empty()
        || !options.primer_options.reverse_tail.is_empty()
        || options.outer_forward_tail.len() > collagene_core::cloning::pcr::primers::MAX_OLIGO
        || options.outer_reverse_tail.len() > collagene_core::cloning::pcr::primers::MAX_OLIGO
    {
        return Err("Use bounded explicit outer tails; native primer-call tails are derived from the declared edit.".into());
    }
    let tails = current_outer_tails(draft, index)?;
    if options.outer_forward_tail != tails.0 || options.outer_reverse_tail != tails.1 {
        return Err("Existing outer-PCR tails changed; review their explicit settings before proposing primers.".into());
    }
    let (inputs, base, _) = draft.inputs_request(workspace)?;
    let original = base
        .fragments
        .get(index)
        .ok_or("The exact source role is unavailable.")?
        .clone();
    let material_binding = base
        .material_bindings
        .iter()
        .find(|b| b.input_id == original.input_id)
        .map(|b| native::BoundOrigin {
            event: b.event.clone(),
            origin: b.origin.clone(),
        });
    let intent_context_revision = independent_intent(draft, intent, &inputs, base, workspace)?;
    // PlacedFragment retains the actual source core after explicit junction
    // trims, before internal edits. Its source coordinates may be normalized
    // across the origin; subtracting its start from the untrimmed start is wrong.
    let retained = intent
        .fragments()
        .get(index)
        .ok_or("The retained intended source core is unavailable.")?;
    if retained.role != original.role
        || retained.input_id != original.input_id
        || retained.input_revision != current_source.revision
        || retained.selector != original.selector
    {
        return Err("The intended retained core differs from its original source role.".into());
    }
    let right_trim = retained
        .original_sequence
        .len()
        .checked_sub(retained.retained_end)
        .ok_or("The retained source window exceeds its original contribution.")?;
    let adjusted = retained
        .selector
        .trim(retained.retained_start, right_trim, retained.source_length)
        .map_err(|e| e.to_string())?;
    let (start, end, reverse, template_strand) = match adjusted {
        collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
            start,
            end,
            reverse,
        } => (start, end, reverse, collagene_core::chemistry::Strand::Top),
        collagene_core::cloning::assembly::intent::Selector::MaterialStrandRange {
            start,
            end,
            strand,
            expected_state_revision,
        } => {
            if retained.source_state_revision.as_deref() != Some(expected_state_revision.as_str())
                || retained.source_event.as_ref() != Some(current_source)
            {
                return Err(
                    "The retained actual strand differs from the exact source state/event.".into(),
                );
            }
            (start, end, false, strand)
        }
    };
    let fragment = collagene_core::cloning::overlap::primer_design::Fragment {
        role: retained.role.clone(),
        input_id: retained.input_id.clone(),
        start,
        end,
        reverse,
    };
    let incoming = index
        .checked_sub(1)
        .or_else(|| draft.circular.then(|| draft.slots.len() - 1));
    let trim = incoming
        .and_then(|junction| {
            draft
                .junction_edits
                .iter()
                .find(|edit| edit.junction == junction)
        })
        .map_or(0, |edit| edit.remove_right);
    let mut edits = draft
        .edits
        .iter()
        .filter(|e| e.role == fragment.role)
        .cloned()
        .collect::<Vec<_>>();
    for e in draft
        .substitutions
        .iter()
        .filter(|e| e.role == fragment.role)
    {
        edits.push(assembly_edit_intent::Edit {
            role: e.role.clone(),
            start: e.offset,
            end: e
                .offset
                .checked_add(1)
                .ok_or("Substitution coordinates overflow.")?,
            expected: e.expected.clone(),
            replacement: e.replacement.clone(),
        });
    }
    // Local edit declarations use the original contribution's construct
    // orientation. Only the declared incoming trim shifts that coordinate zero.
    for edit in &mut edits {
        edit.start = edit
            .start
            .checked_sub(trim)
            .ok_or("An internal edit lies outside the retained incoming source core.")?;
        edit.end = edit
            .end
            .checked_sub(trim)
            .ok_or("An internal edit lies outside the retained incoming source core.")?;
    }
    edits.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));
    if edits.is_empty() {
        return Err("The selected contribution has no declared internal edit.".into());
    }
    let native = native::Request {
        template_strand,
        material_binding,
        fragment,
        edits,
        options: options.clone(),
    };
    native::preflight(&captured.molecule, &native, workspace).map_err(|e| e.to_string())?;
    let mut capture = Capture {
        context_revision: String::new(),
        selected_index: index,
        source: current_source.clone(),
        role: native.fragment.role.clone(),
        input_id: native.fragment.input_id.clone(),
        intent_context_revision,
        target_revision: intent.target().revision(),
        draft_revision: receipt(draft, MAX_CONTEXT)?,
        options_revision: receipt(options, MAX_REPORT)?,
        design_revision: native::design_revision(&captured.molecule, &native, workspace)
            .map_err(|e| e.to_string())?,
    };
    capture.context_revision = receipt(&(ENGINE_VERSION, &capture), MAX_REPORT)?;
    Ok((capture, captured.molecule.clone(), native))
}
fn check_expected(
    intent: &IntendedReference,
    index: usize,
    request: &native::Request,
    expected: &Molecule,
) -> Result<(), String> {
    let retained = intent
        .fragments()
        .get(index)
        .ok_or("Retained contribution unavailable")?;
    let source = retained
        .original_sequence
        .as_str()
        .get(retained.retained_start..retained.retained_end)
        .ok_or("Retained contribution bounds changed")?;
    let mut text = String::new();
    let mut cursor = 0;
    for edit in &request.edits {
        if edit.start < cursor
            || edit.end < edit.start
            || edit.end > source.len()
            || &source[edit.start..edit.end] != edit.expected
        {
            return Err("Literal edit differs from retained source calls".into());
        }
        text.push_str(&source[cursor..edit.start]);
        text.push_str(&edit.replacement);
        cursor = edit.end;
    }
    text.push_str(&source[cursor..]);
    let sequence = Sequence::dna(&text).map_err(|e| e.to_string())?;
    let physical = if request.fragment.reverse {
        sequence.reverse_complement().map_err(|e| e.to_string())?
    } else {
        sequence
    };
    let reverse_tail = if request.options.outer_reverse_tail.is_empty() {
        String::new()
    } else {
        Sequence::dna(&request.options.outer_reverse_tail)
            .and_then(|s| s.reverse_complement())
            .map_err(|e| e.to_string())?
            .as_str()
            .to_owned()
    };
    let dna = format!(
        "{}{}{}",
        request.options.outer_forward_tail,
        physical.as_str(),
        reverse_tail.as_str()
    );
    if expected.topology != Topology::Linear || expected.sequence.as_str() != dna {
        return Err("Proposed preparation differs from independently edited retained contribution and explicit tails".into());
    }
    Ok(())
}
pub fn from_design(
    draft: &Draft,
    intent: &IntendedReference,
    index: usize,
    current_source: &EventBinding,
    options: &native::Options,
    report: native::Report,
    workspace: Option<&Workspace>,
) -> Result<Proposal, String> {
    let worker_receipt = receipt(&report, MAX_REPORT)?;
    let (capture, source, native) =
        request(draft, intent, index, current_source, options, workspace)?;
    report
        .validate_sources(&source, &capture.design_revision, workspace)
        .map_err(|e| e.to_string())?;
    if report.engine_version != native::ENGINE_VERSION
        || report.source_revision != source.revision()
        || receipt(&report.request, MAX_REPORT)? != receipt(&native, MAX_REPORT)?
    {
        return Err("The native primer result belongs to another source or edit request.".into());
    }
    check_expected(intent, index, &native, &report.expected)?;
    let checked = crate::assembly::review::assembly_product_model::Snapshot::new(
        report.plan.clone(),
        report.expected.clone(),
        report.verification.clone(),
    );
    if let Some(error) = checked.prefix_error {
        return Err(format!("The diagnostic worker receipt is unbound: {error}"));
    }
    let preparation = report.preparation.clone();
    let verification = report.verification.clone();
    let warnings = report.warnings.clone();
    let public_receipt = receipt(
        &(&capture, &preparation, &verification, &warnings),
        MAX_CONTEXT,
    )?;
    Ok(Proposal {
        capture,
        preparation,
        verification,
        warnings,
        worker: Arc::new(report),
        worker_receipt,
        public_receipt,
    })
}
fn raw_primer(primer: &Primer) -> ExtensionPrimerInput {
    ExtensionPrimerInput {
        name: primer.name.clone(),
        annealing: primer.annealing.as_str().into(),
        tail: primer.tail.clone(),
    }
}
fn raw_pair(
    pair: &collagene_core::cloning::overlap::extension::pcr::PrimerPair,
) -> ExtensionPrimerPairInput {
    ExtensionPrimerPairInput {
        template_strand: pair.template_strand,
        forward: raw_primer(&pair.forward),
        reverse: raw_primer(&pair.reverse),
    }
}
impl Proposal {
    pub fn validate_binding(
        &self,
        draft: &Draft,
        intent: &IntendedReference,
        current_source: &EventBinding,
        options: &native::Options,
        workspace: Option<&Workspace>,
    ) -> Result<(), String> {
        if receipt(self.worker.as_ref(), MAX_REPORT)? != self.worker_receipt
            || receipt(
                &(
                    &self.capture,
                    &self.preparation,
                    &self.verification,
                    &self.warnings,
                ),
                MAX_CONTEXT,
            )? != self.public_receipt
        {
            return Err("The retained primer proposal or native worker receipt changed.".into());
        }
        let (capture, _, native) = request(
            draft,
            intent,
            self.capture.selected_index,
            current_source,
            options,
            workspace,
        )?;
        if capture != self.capture
            || receipt(&self.worker.request, MAX_REPORT)? != receipt(&native, MAX_REPORT)?
        {
            return Err(
                "The primer proposal is stale for the current draft, intent, source or settings."
                    .into(),
            );
        }
        self.worker
            .validate_binding(&capture.design_revision)
            .map_err(|e| e.to_string())
    }
    pub fn apply(
        &self,
        draft: &Draft,
        intent: &IntendedReference,
        current_source: &EventBinding,
        options: &native::Options,
        workspace: Option<&Workspace>,
    ) -> Result<Draft, String> {
        self.validate_binding(draft, intent, current_source, options, workspace)?;
        let mut next = draft.clone();
        let row = next
            .methods
            .rows
            .get_mut(self.capture.selected_index)
            .ok_or("The selected source contribution is unavailable.")?;
        let o = &self.preparation.extension;
        let raw = OverlapExtensionPcrInput {
            left: raw_pair(&self.preparation.left),
            right: raw_pair(&self.preparation.right),
            outer: raw_pair(&self.preparation.outer),
            overlap: o.overlap.to_string(),
            min_paired_bases: o.min_paired_bases.to_string(),
            max_comparisons: o.max_comparisons.to_string(),
            max_candidates: o.max_candidates.to_string(),
        };
        raw.validate_raw()?;
        row.preparation = proposed_preparation(row.preparation);
        row.overlap_extension_pcr = Some(raw);
        next.sync_rows();
        next.validate()?;
        Ok(next)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::draft::model::CapturedSource;
    use crate::assembly::edit::assembly_edit_model;
    use crate::cloning::preparation::cloning_methods_model::Mode;
    use crate::cloning::preparation::cloning_methods_model::Treatment;
    use crate::sources::identity::SourceInfo;
    use collagene_core::{
        chemistry::EndRepairMethod,
        material::lineage::InputEventBinding,
        planning::plan::{DesignPlan, Operation, PlanInput, Source},
        verification::evidence::Status,
    };
    use collagene_model::Alphabet;
    const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
    fn event(frame: &Workspace, index: usize) -> EventBinding {
        let doc = &frame.documents[index];
        EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        }
    }
    fn captured(molecule: Molecule, binding: EventBinding) -> Arc<CapturedSource> {
        Arc::new(CapturedSource {
            source: SourceInfo {
                binding,
                name: molecule.name.clone(),
                length: molecule.sequence.len(),
                alphabet: Alphabet::Dna,
                topology: molecule.topology,
            },
            molecule,
        })
    }
    fn intent(draft: &Draft, workspace: Option<&Workspace>) -> IntendedReference {
        let (inputs, base, _) = draft.inputs_request(workspace).unwrap();
        IntendedReference::Edits(
            assembly_edit_intent::build(&inputs, &draft.edit_request(base), workspace).unwrap(),
        )
    }
    fn options(draft: &Draft) -> native::Options {
        let mut options = native::Options::default();
        options.primer_options.min_length = 8;
        options.primer_options.max_length = 8;
        options.primer_options.max_pairs = 1;
        options.primer_options.max_tm_difference_c = 100.;
        options.primer_options.min_gc_fraction = 0.;
        options.primer_options.max_gc_fraction = 1.;
        options.primer_options.max_homopolymer = 200;
        options.primer_options.max_three_prime_complementarity = 200;
        options.homology_flank = 4;
        options.min_paired_bases = 8;
        let tails = current_outer_tails(draft, draft.active_slot).unwrap();
        options.outer_forward_tail = tails.0;
        options.outer_reverse_tail = tails.1;
        options
    }
    fn fixture(reverse_wrapped: bool) -> (Draft, IntendedReference, Workspace, native::Options) {
        let vector =
            Molecule::dna("Other original source", "CCGGAATTGGCC", Topology::Linear).unwrap();
        let source = Molecule::dna(
            "Original edited source",
            ORIGINAL,
            if reverse_wrapped {
                Topology::Circular
            } else {
                Topology::Linear
            },
        )
        .unwrap();
        let frame = Workspace::from_documents(
            "Exact source events",
            vec![vector.clone(), source.clone()],
            0,
        )
        .unwrap();
        let mut draft = Draft::default();
        draft.circular = false;
        draft.active_slot = 1;
        draft.methods.mode = Mode::Overlap;
        draft.methods.rows[0].preparation = OverlapPreparation::WholeLinear;
        draft.slots[0].captured = Some(captured(vector, event(&frame, 0)));
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                end: value,
                ..
            } = &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 12;
        }
        draft.slots[1].captured = Some(captured(source, event(&frame, 1)));
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                start: value,
                ..
            } = &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = if reverse_wrapped { 32 } else { 0 };
        }
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                end: value,
                ..
            } = &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = if reverse_wrapped { 70 } else { 40 };
        }
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                reverse: value,
                ..
            } = &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = reverse_wrapped;
        }
        draft.methods.rows[1].preparation = OverlapPreparation::ExactPcr;
        draft.methods.rows[1].repair = Some(EndRepairMethod::Blunt);
        draft.methods.rows[1].treatment = Treatment::Phosphorylate;
        draft.methods.rows[1].forward_name = "Previous forward".into();
        draft.methods.rows[1].forward = "AAAA".into();
        draft.methods.rows[1].reverse_name = "Previous reverse".into();
        draft.methods.rows[1].reverse_primer = "TTTT".into();
        draft.methods.rows[1].fragment = "7".into();
        draft.sync_rows();
        let mut draft = assembly_edit_model::declare(&draft, 1, "17", "19", "T", None).unwrap();
        let role = draft.methods.role_ids(2).unwrap()[1].role.clone();
        let expected = assembly_edit_model::source_text(&draft, 1, 19, 20, None).unwrap();
        draft.substitutions.push(
            collagene_core::cloning::assembly::substitution_intent::Substitution {
                role,
                offset: 19,
                replacement: if expected == "A" { "C" } else { "A" }.into(),
                expected,
            },
        );
        draft.sync_rows();
        draft.validate().unwrap();
        let intent = intent(&draft, None);
        let options = options(&draft);
        (draft, intent, frame, options)
    }
    fn proposal(
        draft: &Draft,
        intent: &IntendedReference,
        frame: &Workspace,
        options: &native::Options,
    ) -> Proposal {
        let source = attest(draft, draft.active_slot, frame).unwrap();
        let (_, molecule, request) = request(
            draft,
            intent,
            draft.active_slot,
            &source,
            options,
            Some(frame),
        )
        .unwrap();
        let report = native::design(&molecule, &request, Some(frame)).unwrap();
        from_design(
            draft,
            intent,
            draft.active_slot,
            &source,
            options,
            report,
            Some(frame),
        )
        .unwrap()
    }
    #[test]
    fn manual_digest_cannot_be_silently_removed_by_a_primer_proposal() {
        let (mut draft, _, frame, _) = fixture(false);
        draft.methods.rows[1].preparation = OverlapPreparation::Digest;
        draft.methods.rows[1].enzymes = "EcoRI".into();
        draft.sync_rows();
        let intent = intent(&draft, None);
        let options = options(&draft);
        let before = receipt(&draft, MAX_CONTEXT).unwrap();
        let message = request(&draft, &intent, 1, &event(&frame, 1), &options, None)
            .err()
            .unwrap();
        assert!(message.contains("catalog"), "{message}");
        assert_eq!(receipt(&draft, MAX_CONTEXT).unwrap(), before);
    }

    #[test]
    fn proposed_internal_edits_keep_downstream_catalog_request_and_independent_target() {
        let (mut draft, _, frame, _) = fixture(false);
        let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
        let id = catalog
            .edition()
            .enzymes
            .iter()
            .find(|e| e.name == "EcoRI")
            .unwrap()
            .id
            .clone();
        let pin = collagene_core::cloning::restriction::catalog::select(&catalog, &[id]).unwrap();
        draft.methods.rows[1].preparation = OverlapPreparation::PcrCatalogDigest;
        draft.methods.rows[1].catalog_selection = Some(pin.clone());
        draft.methods.rows[1].forward_tail = "CCGAATTC".into();
        draft.methods.rows[1].reverse_tail = "GGGAATTC".into();
        draft.sync_rows();
        let intent = intent(&draft, None);
        let target = intent.target().revision();
        let options = options(&draft);
        let p = proposal(&draft, &intent, &frame, &options);
        let after = p
            .apply(&draft, &intent, &p.capture.source, &options, None)
            .unwrap();
        assert_eq!(
            after.version,
            crate::assembly::draft::model::CURRENT_VERSION
        );
        assert!(
            after.methods.rows[1].preparation
                == OverlapPreparation::OverlapExtensionPcrCatalogDigest
        );
        assert_eq!(after.methods.rows[1].catalog_selection.as_ref(), Some(&pin));
        assert_eq!(after.methods.rows[1].fragment, "7");
        assert_eq!(after.methods.rows[1].repair, draft.methods.rows[1].repair);
        assert_eq!(intent.target().revision(), target);
        crate::assembly::primer_design::assembly_extension_primer_commit::validate(
            &draft, &after, &draft,
        )
        .unwrap();
        // Subsequent proposal requests retain the six-oligo editor's explicit tails,
        // even when unrelated legacy single-pair fields contain other text.
        let mut changed = after.clone();
        changed.methods.rows[1].forward_tail = "TTTT".into();
        assert_eq!(
            current_outer_tails(&changed, 1).unwrap(),
            ("CCGAATTC".into(), "GGGAATTC".into())
        );
    }

    #[test]
    fn mixed_internal_edits_propose_six_oligos_and_change_only_selected_preparation() {
        let (draft, intent, frame, options) = fixture(false);
        assert_eq!(
            intent.target().sequence.as_str(),
            "CCGGAATTGGCCAGTCCGATGCTACGTATAATTCAGACCTGTAGATCGTCA"
        );
        let before = receipt(&draft, MAX_CONTEXT).unwrap();
        let target = intent.target().revision();
        let p = proposal(&draft, &intent, &frame, &options);
        assert_eq!(p.worker.request.edits.len(), 2);
        assert_eq!(
            p.worker.expected.sequence.as_str(),
            "AGTCCGATGCTACGTATAATTCAGACCTGTAGATCGTCA"
        );
        let after = p
            .apply(&draft, &intent, &p.capture.source, &options, None)
            .unwrap();
        let mut expected = draft.clone();
        expected.version = crate::assembly::draft::model::CURRENT_VERSION;
        expected.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcr;
        expected.methods.rows[1].overlap_extension_pcr =
            after.methods.rows[1].overlap_extension_pcr.clone();
        assert!(after == expected);
        assert_eq!(receipt(&draft, MAX_CONTEXT).unwrap(), before);
        assert_eq!(intent.target().revision(), target);
        assert_eq!(
            after.slots[1]
                .captured
                .as_ref()
                .unwrap()
                .molecule
                .sequence
                .as_str(),
            ORIGINAL
        );
        let raw = after.methods.rows[1]
            .overlap_extension_pcr
            .as_ref()
            .unwrap();
        assert_eq!(
            receipt(&raw.request().unwrap(), MAX_REPORT).unwrap(),
            receipt(&p.preparation, MAX_REPORT).unwrap()
        );
        assert_eq!(raw.left.forward.name, p.preparation.left.forward.name);
        assert_eq!(
            raw.right.reverse.annealing,
            p.preparation.right.reverse.annealing.as_str()
        );
        assert_eq!(raw.outer.forward.tail, p.preparation.outer.forward.tail);
        assert_eq!(raw.outer.reverse.tail, p.preparation.outer.reverse.tail);
    }
    #[test]
    fn partial_reversed_wrapped_source_and_literal_outer_tails_stay_bound() {
        let (mut draft, _, frame, _) = fixture(true);
        draft.methods.rows[1].forward_tail = "AAG".into();
        draft.methods.rows[1].reverse_tail = "CCT".into();
        draft.sync_rows();
        let intent = intent(&draft, None);
        let options = options(&draft);
        let p = proposal(&draft, &intent, &frame, &options);
        assert_eq!(p.worker.request.fragment.start, 32);
        assert_eq!(p.worker.request.fragment.end, 70);
        assert!(p.worker.request.fragment.reverse);
        assert!(p.worker.expected.sequence.as_str().starts_with("AAG"));
        assert!(p.worker.expected.sequence.as_str().ends_with("AGG"));
        let after = p
            .apply(&draft, &intent, &p.capture.source, &options, None)
            .unwrap();
        assert!(after.slots == draft.slots);
        assert!(after.slots[1].selector.reference_interval().unwrap().2);
        assert_eq!(
            after.methods.rows[1]
                .overlap_extension_pcr
                .as_ref()
                .unwrap()
                .outer
                .forward
                .tail,
            "AAG"
        );
        assert_eq!(
            after.methods.rows[1]
                .overlap_extension_pcr
                .as_ref()
                .unwrap()
                .outer
                .reverse
                .tail,
            "CCT"
        );
        let mut wrong = options.clone();
        wrong.outer_forward_tail.clear();
        assert!(request(&draft, &intent, 1, &p.capture.source, &wrong, None).is_err());
    }
    #[test]
    fn incoming_junction_trim_prepares_only_the_retained_edited_core() {
        let (mut draft, _, _, _) = fixture(false);
        let vector = draft.slots[0].captured.as_ref().unwrap().molecule.clone();
        let original = Molecule::dna(
            "Prefix-bearing original",
            &format!("GTT{ORIGINAL}"),
            Topology::Linear,
        )
        .unwrap();
        let frame = Workspace::from_documents(
            "Trimmed original events",
            vec![vector.clone(), original.clone()],
            0,
        )
        .unwrap();
        draft.slots[0].captured = Some(captured(vector, event(&frame, 0)));
        draft.slots[1].captured = Some(captured(original.clone(), event(&frame, 1)));
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                start: value,
                ..
            } = &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 0;
        }
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                end: value,
                ..
            } = &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 43;
        }
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                reverse: value,
                ..
            } = &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = false;
        }
        draft.edits.clear();
        draft.substitutions.clear();
        draft.junction_edits = vec![collagene_core::cloning::assembly::junction_edit::Edit {
            junction: 0,
            remove_left: 0,
            remove_right: 2,
            inserted: String::new(),
        }];
        draft.sync_rows();
        let draft = assembly_edit_model::declare(&draft, 1, "20", "24", "TT", None).unwrap();
        let intended = intent(&draft, None);
        let options = options(&draft);
        let target_revision = intended.target().revision();
        assert_eq!(
            intended.target().sequence.as_str(),
            "CCGGAATTGGCCTAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA"
        );
        let p = proposal(&draft, &intended, &frame, &options);
        assert_eq!(p.worker.request.fragment.start, 2);
        assert_eq!(p.worker.request.fragment.end, 43);
        assert_eq!(p.worker.request.edits[0].start, 17);
        assert_eq!(p.worker.request.edits[0].end, 21);
        assert_eq!(p.worker.request.edits[0].expected, "GGAC");
        assert_eq!(
            p.worker.expected.sequence.as_str(),
            "TAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA"
        );
        assert_eq!(p.worker.expected.sequence.len(), 39);
        let after = p
            .apply(&draft, &intended, &p.capture.source, &options, None)
            .unwrap();
        assert!(after.slots == draft.slots);
        assert_eq!(after.slots[1].captured.as_ref().unwrap().molecule, original);
        assert!(after.edits == draft.edits && after.junction_edits == draft.junction_edits);
        assert_eq!(intended.target().revision(), target_revision);
        assert_eq!(
            after.methods.rows[1]
                .overlap_extension_pcr
                .as_ref()
                .unwrap()
                .request()
                .unwrap()
                .outer
                .forward
                .ordered_sequence(),
            p.preparation.outer.forward.ordered_sequence()
        );
    }
    #[test]
    fn incoming_trim_uses_construct_coordinates_with_reverse_and_origin_normalization() {
        for (
            start,
            end,
            reverse,
            trim,
            edit_start,
            edit_end,
            snp_offset,
            retained_start,
            retained_end,
        ) in [
            (32, 70, true, 2, 18, 20, 21, 32, 68),
            (38, 78, false, 3, 19, 21, 22, 1, 38),
        ] {
            let (mut draft, _, _, _) = fixture(false);
            let vector = draft.slots[0].captured.as_ref().unwrap().molecule.clone();
            let source = Molecule::dna("Wrapped original", ORIGINAL, Topology::Circular).unwrap();
            let frame = Workspace::from_documents(
                "Wrapped trim events",
                vec![vector.clone(), source.clone()],
                0,
            )
            .unwrap();
            draft.slots[0].captured = Some(captured(vector, event(&frame, 0)));
            draft.slots[1].captured = Some(captured(source, event(&frame, 1)));
            {
                let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    start: value,
                    ..
                } = &mut draft.slots[1].selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = start;
            }
            {
                let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    end: value,
                    ..
                } = &mut draft.slots[1].selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = end;
            }
            {
                let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    reverse: value,
                    ..
                } = &mut draft.slots[1].selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = reverse;
            }
            draft.edits.clear();
            draft.substitutions.clear();
            draft.junction_edits = vec![collagene_core::cloning::assembly::junction_edit::Edit {
                junction: 0,
                remove_left: 0,
                remove_right: trim,
                inserted: String::new(),
            }];
            draft.sync_rows();
            let mut draft = assembly_edit_model::declare(
                &draft,
                1,
                &(edit_start + 1).to_string(),
                &(edit_end + 1).to_string(),
                "T",
                None,
            )
            .unwrap();
            let expected =
                assembly_edit_model::source_text(&draft, 1, snp_offset, snp_offset + 1, None)
                    .unwrap();
            draft.substitutions.push(
                collagene_core::cloning::assembly::substitution_intent::Substitution {
                    role: draft.methods.role_ids(2).unwrap()[1].role.clone(),
                    offset: snp_offset,
                    replacement: if expected == "A" { "C" } else { "A" }.into(),
                    expected,
                },
            );
            draft.sync_rows();
            let intended = intent(&draft, None);
            let options = options(&draft);
            let source = attest(&draft, 1, &frame).unwrap();
            let (_, _, native) = request(&draft, &intended, 1, &source, &options, None).unwrap();
            assert_eq!(
                (
                    native.fragment.start,
                    native.fragment.end,
                    native.fragment.reverse
                ),
                (retained_start, retained_end, reverse)
            );
            assert_eq!((native.edits[0].start, native.edits[0].end), (16, 18));
            assert_eq!((native.edits[1].start, native.edits[1].end), (19, 20));
        }
    }
    #[test]
    fn blank_oe_editor_and_unfinished_inner_oligos_allow_native_proposals_with_retained_tails() {
        let (mut draft, _, frame, _) = fixture(false);
        draft.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcr;
        draft.methods.rows[1].overlap_extension_pcr = None;
        draft.methods.rows[1].forward_tail = "Inactive old PCR tail".into();
        draft.sync_rows();
        assert_eq!(
            current_outer_tails(&draft, 1).unwrap(),
            (String::new(), String::new())
        );
        let intended = intent(&draft, None);
        let empty_options = options(&draft);
        let p = proposal(&draft, &intended, &frame, &empty_options);
        let after = p
            .apply(&draft, &intended, &p.capture.source, &empty_options, None)
            .unwrap();
        assert!(after.methods.rows[1].overlap_extension_pcr.is_some());
        assert!(
            after.methods.rows[1]
                .overlap_extension_pcr
                .as_ref()
                .unwrap()
                .outer
                .forward
                .tail
                .is_empty()
        );
        assert_eq!(after.methods.rows[1].forward_tail, "Inactive old PCR tail");

        let mut raw = OverlapExtensionPcrInput::default();
        raw.left.forward.annealing = "unfinished domain".into();
        raw.right.reverse.tail = "unfinished inner tail".into();
        raw.outer.forward.annealing.clear();
        raw.outer.forward.tail = " a a g ".into();
        raw.outer.reverse.tail = "c c t\n".into();
        draft.methods.rows[1].overlap_extension_pcr = Some(raw);
        draft.sync_rows();
        assert_eq!(
            current_outer_tails(&draft, 1).unwrap(),
            ("AAG".into(), "CCT".into())
        );
        let intended = intent(&draft, None);
        let retained_options = options(&draft);
        let p = proposal(&draft, &intended, &frame, &retained_options);
        let after = p
            .apply(
                &draft,
                &intended,
                &p.capture.source,
                &retained_options,
                None,
            )
            .unwrap();
        assert_eq!(
            after.methods.rows[1]
                .overlap_extension_pcr
                .as_ref()
                .unwrap()
                .outer
                .forward
                .tail,
            "AAG"
        );
        assert_eq!(
            after.methods.rows[1]
                .overlap_extension_pcr
                .as_ref()
                .unwrap()
                .outer
                .reverse
                .tail,
            "CCT"
        );
        draft.methods.rows[1]
            .overlap_extension_pcr
            .as_mut()
            .unwrap()
            .outer
            .forward
            .tail = "ACN".into();
        assert!(current_outer_tails(&draft, 1).is_err());
    }
    #[test]
    fn stale_source_intent_options_and_mutated_public_proposal_cannot_apply() {
        let (draft, intended, frame, options) = fixture(false);
        let p = proposal(&draft, &intended, &frame, &options);
        let mut multisite = draft.clone();
        multisite.methods.mode = Mode::MultiSiteReference;
        assert!(request(&multisite, &intended, 1, &p.capture.source, &options, None).is_err());
        assert!(
            p.apply(&multisite, &intended, &p.capture.source, &options, None)
                .is_err()
        );
        let mut wrong = p.capture.source.clone();
        wrong.document_id.push_str("other");
        assert!(p.apply(&draft, &intended, &wrong, &options, None).is_err());
        let mut changed = draft.clone();
        changed.methods.rows[0].forward_name = "Unrelated settings changed".into();
        assert!(
            p.apply(&changed, &intended, &p.capture.source, &options, None)
                .is_err()
        );
        let mut changed = draft.clone();
        changed.active_slot = 0;
        assert!(
            p.apply(&changed, &intended, &p.capture.source, &options, None)
                .is_err()
        );
        let mut changed = options.clone();
        changed.homology_flank += 1;
        assert!(
            p.apply(&draft, &intended, &p.capture.source, &changed, None)
                .is_err()
        );
        let mut changed = intended.clone();
        let IntendedReference::Edits(report) = &mut changed else {
            unreachable!()
        };
        report.target.name.push_str(" altered");
        assert!(
            p.apply(&draft, &changed, &p.capture.source, &options, None)
                .is_err()
        );
        let mut changed = p.clone();
        changed.preparation.left.forward.tail.push('A');
        assert!(
            changed
                .apply(&draft, &intended, &p.capture.source, &options, None)
                .is_err()
        );
        let mut changed = p.clone();
        changed.capture.role.push_str("other");
        assert!(
            changed
                .apply(&draft, &intended, &p.capture.source, &options, None)
                .is_err()
        );
        let mut bad_frame = frame.clone();
        let id = draft.slots[0]
            .captured
            .as_ref()
            .unwrap()
            .source
            .binding
            .event_id
            .clone();
        bad_frame
            .history
            .iter_mut()
            .find(|e| e.id == id)
            .unwrap()
            .operation
            .push_str(" forged");
        assert!(attest(&draft, 1, &bad_frame).is_err());
        let mut report = p.worker.as_ref().clone();
        report.expected.name.push_str(" forged");
        assert!(
            from_design(
                &draft,
                &intended,
                1,
                &p.capture.source,
                &options,
                report,
                None
            )
            .is_err()
        );
    }
    #[test]
    fn source_metadata_is_diagnostic_and_saved_origin_survives_for_fresh_recipe_replay() {
        let (mut draft, _, _, _) = fixture(false);
        let other = draft.slots[0].captured.as_ref().unwrap().molecule.clone();
        let original = draft.slots[1].captured.as_ref().unwrap().molecule.clone();
        let target = Molecule::dna(
            "Separate original-stock PCR requirement",
            ORIGINAL,
            Topology::Linear,
        )
        .unwrap();
        let frame = Workspace::from_documents(
            "Saved stock origin",
            vec![other, original.clone(), target.clone()],
            0,
        )
        .unwrap();
        let plan = DesignPlan {
            name: "Original stock amplification".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "stock".into(),
                revision: original.revision(),
                molecule: original,
            }],
            steps: vec![Operation::Pcr {
                template_strand: collagene_core::chemistry::Strand::Top,
                template: Source::Input("stock".into()),
                forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            }],
        };
        let commit = collagene_core::material::outputs::finalize(
            &frame,
            &plan,
            event(&frame, 2),
            vec![InputEventBinding {
                input_id: "stock".into(),
                event: event(&frame, 1),
            }],
            &collagene_core::material::outputs::Options {
                name: "Saved original stock".into(),
                save_primers: false,
                show_prepared: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let selection = collagene_core::material::run::Selection {
            run_id: commit.receipt.run.run_id.clone(),
            step_index: 0,
            event: commit.receipt.construct.clone(),
        };
        let origin =
            collagene_core::material::capsule::capsule(&commit.workspace, &selection).unwrap();
        draft.slots[0].captured = Some(captured(
            commit.workspace.revisions[&frame.documents[0].head].clone(),
            event(&frame, 0),
        ));
        draft.slots[1].captured = Some(captured(
            commit.workspace.revisions[&selection.event.revision].clone(),
            selection.event.clone(),
        ));
        draft.methods.rows[1].run_origin = Some(
            crate::history::run_origin_model::PinnedOrigin::from_material(
                collagene_core::material::run::MaterialOrigin::Run(origin.clone()),
            )
            .unwrap(),
        );
        draft.slots[1].material_origin = draft.methods.rows[1].run_origin.clone();
        draft.sync_rows();
        let intended = intent(&draft, Some(&commit.workspace));
        let options = options(&draft);
        let p = proposal(&draft, &intended, &commit.workspace, &options);
        assert_eq!(
            p.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            p.verification.findings
        );
        let after = p
            .apply(
                &draft,
                &intended,
                &p.capture.source,
                &options,
                Some(&commit.workspace),
            )
            .unwrap();
        assert_eq!(
            receipt(
                after.methods.rows[1].run_origin.as_ref().unwrap(),
                MAX_CONTEXT
            )
            .unwrap(),
            receipt(
                draft.methods.rows[1].run_origin.as_ref().unwrap(),
                MAX_CONTEXT
            )
            .unwrap()
        );
        assert_eq!(
            after.slots[1].captured.as_ref().unwrap().source.binding,
            selection.event
        );
        let (inputs, _, _) = after.inputs_request(Some(&commit.workspace)).unwrap();
        let (_, fragments) = crate::cloning::preparation::cloning_preparations::prepared_rows(
            &after.methods,
            inputs.len(),
        )
        .unwrap();
        let compiled = collagene_core::cloning::assembly::recipe::compile_preparations(
            &inputs,
            &p.worker.expected,
            "Actual saved material preparation",
            &fragments,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap();
        let physical = DesignPlan {
            name: "Fresh selected physical preparation".into(),
            target_revision: p.worker.expected.revision(),
            inputs,
            steps: compiled.steps,
        };
        let checked = collagene_core::verification::engine::verify(&physical, &p.worker.expected);
        assert_eq!(
            checked.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            checked.findings
        );
    }
    #[test]
    fn applying_native_oligos_preserves_bottom_frame_and_outer_top_frame() {
        let bottom = collagene_core::cloning::overlap::extension::pcr::PrimerPair {
            template_strand: collagene_core::chemistry::Strand::Bottom,
            forward: Primer::new("Bottom extension endpoint", "GGGGTGAC", "").unwrap(),
            reverse: Primer::new("Bottom opposite endpoint", "AGTCCGAT", "").unwrap(),
        };
        let raw = raw_pair(&bottom);
        assert_eq!(
            raw.template_strand,
            collagene_core::chemistry::Strand::Bottom
        );
        assert_eq!(raw.forward.annealing, "GGGGTGAC");
        let outer = collagene_core::cloning::overlap::extension::pcr::PrimerPair {
            template_strand: collagene_core::chemistry::Strand::Top,
            ..bottom
        };
        assert_eq!(
            raw_pair(&outer).template_strand,
            collagene_core::chemistry::Strand::Top
        );
    }
}
