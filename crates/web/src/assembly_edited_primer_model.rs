//! Atomic declaration adapter for all original-source overlap-assembly primers.
//! Independent typed intent remains unchanged; native reports are diagnostic
//! until the ordinary recipe is freshly reconstructed against that intent.
use crate::{
    assembly_draft_model::{self, Draft},
    assembly_review::IntendedReference,
    cloning_methods_model::{
        ExtensionPrimerInput, ExtensionPrimerPairInput, Mode, OverlapExtensionPcrInput,
        OverlapPreparation,
    },
};
#[cfg(test)]
use collagene_core::Molecule;
use collagene_core::{
    analysis::Primer,
    assembly_edit_intent, assembly_edited_primer_design as native, overlap_extension_primer_design,
    run_lineage::EventBinding,
    verification::{PlanInput, VerificationReport},
    workspace::Workspace,
};
use serde::{Deserialize, Serialize};
use std::{io::Write, sync::Arc};
const ENGINE_VERSION: &str = "assembly-edited-primer-proposal-1";
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
            .ok_or_else(|| std::io::Error::other("bounded whole-assembly primer receipt"))?;
        self.hash.update(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn receipt(v: &impl Serialize, max: usize) -> Result<String, String> {
    let mut w = Digest {
        hash: blake3::Hasher::new(),
        bytes: 0,
        max,
    };
    serde_json::to_writer(&mut w, v).map_err(|e| e.to_string())?;
    Ok(w.hash.finalize().to_hex().to_string())
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Options {
    pub overlap: usize,
    pub preparation: overlap_extension_primer_design::Options,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            overlap: 20,
            preparation: overlap_extension_primer_design::Options::default(),
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Capture {
    pub context_revision: String,
    pub sources: Vec<EventBinding>,
    pub intent_kind: String,
    pub intent_context_revision: String,
    pub target_revision: String,
    pub draft_revision: String,
    pub options_revision: String,
    pub design_revision: String,
}
#[derive(Clone)]
pub struct Proposal {
    pub capture: Capture,
    pub recipe: collagene_core::assembly_recipe::Recipe,
    pub verification: VerificationReport,
    pub warnings: Vec<String>,
    worker: Arc<native::Report>,
    worker_receipt: String,
    public_receipt: String,
}
/// Contribution navigation is view state; all declarations remain receipt-bound.
pub(crate) fn context_draft(mut draft: Draft) -> Draft {
    draft.active_slot = 0;
    draft
}
pub fn eligible(draft: &Draft) -> bool {
    draft.methods.mode == Mode::Overlap
        && (2..=12).contains(&draft.slots.len())
        && draft
            .slots
            .iter()
            .all(|slot| slot.captured.is_some() && slot.anneal_bottom.is_none())
        && draft.methods.rows.iter().all(|row| {
            matches!(
                row.preparation,
                OverlapPreparation::WholeLinear
                    | OverlapPreparation::ExactPcr
                    | OverlapPreparation::MismatchPcr
                    | OverlapPreparation::OverlapExtensionPcr
            )
        })
}
pub fn attest(draft: &Draft, frame: &Workspace) -> Result<Vec<EventBinding>, String> {
    crate::assembly_extension_primer_model::attest(draft, draft.active_slot, frame)?;
    draft
        .slots
        .iter()
        .map(|slot| {
            slot.captured
                .as_ref()
                .map(|c| c.source.binding.clone())
                .ok_or_else(|| "Pin every original source contribution first.".into())
        })
        .collect()
}
fn independent(
    draft: &Draft,
    intent: &IntendedReference,
    inputs: &[PlanInput],
    base: collagene_core::assembly_intent::Request,
    workspace: Option<&Workspace>,
) -> Result<(String, String, assembly_edit_intent::Request), String> {
    let declared = assembly_draft_model::edit_request(draft, base);
    let (kind, context) = match intent {
        IntendedReference::Regions(report) => {
            if !draft.edits.is_empty()
                || !draft.substitutions.is_empty()
                || !draft.junction_edits.is_empty()
            {
                return Err(
                    "The frozen plain-region intent differs from this edited draft.".into(),
                );
            }
            let context = collagene_core::assembly_intent::intent_revision(
                inputs,
                &declared.base.base,
                workspace,
            )
            .map_err(|e| e.to_string())?;
            report
                .validate_binding(&context)
                .map_err(|e| e.to_string())?;
            if receipt(&report.request, MAX_CONTEXT)? != receipt(&declared.base.base, MAX_CONTEXT)?
            {
                return Err("The frozen original source declaration changed".into());
            }
            ("regions", context)
        }
        IntendedReference::Junctions(report) => {
            if !draft.edits.is_empty() || !draft.substitutions.is_empty() {
                return Err(
                    "The frozen junction intent does not include these internal edits.".into(),
                );
            }
            let context = collagene_core::assembly_junction_edit::intent_revision(
                inputs,
                &declared.base.base,
                &declared.base.junction_edits,
                workspace,
            )
            .map_err(|e| e.to_string())?;
            report
                .validate_binding(&context)
                .map_err(|e| e.to_string())?;
            if receipt(&report.base_request, MAX_CONTEXT)?
                != receipt(&declared.base.base, MAX_CONTEXT)?
                || receipt(&report.edits, MAX_CONTEXT)?
                    != receipt(&declared.base.junction_edits, MAX_CONTEXT)?
            {
                return Err("The frozen source or junction declaration changed".into());
            }
            ("junctions", context)
        }
        IntendedReference::Substitutions(report) => {
            if !draft.edits.is_empty() {
                return Err(
                    "The frozen substitution intent does not include these replacement edits."
                        .into(),
                );
            }
            let context = collagene_core::assembly_substitution_intent::intent_revision(
                inputs,
                &declared.base,
                workspace,
            )
            .map_err(|e| e.to_string())?;
            report
                .validate_sources(inputs, &context, workspace)
                .map_err(|e| e.to_string())?;
            if receipt(&report.request, MAX_CONTEXT)? != receipt(&declared.base, MAX_CONTEXT)? {
                return Err("The frozen source or substitution declaration changed".into());
            }
            ("substitutions", context)
        }
        IntendedReference::Edits(report) => {
            let context = assembly_edit_intent::intent_revision(inputs, &declared, workspace)
                .map_err(|e| e.to_string())?;
            report
                .validate_sources(inputs, &context, workspace)
                .map_err(|e| e.to_string())?;
            if receipt(&report.request, MAX_CONTEXT)? != receipt(&declared, MAX_CONTEXT)? {
                return Err("The original internal-edit declaration changed.".into());
            }
            ("edits", context)
        }
    };
    Ok((kind.into(), context, declared))
}
fn normalized_tail(text: &str) -> Result<String, String> {
    if text.len() > 1000 {
        return Err("An existing active PCR tail exceeds its field bound.".into());
    }
    let tail = Primer::new("Existing active outer tail", "A", text)
        .map_err(|e| e.to_string())?
        .tail;
    if !tail.bytes().all(|b| b"ACGT".contains(&b)) {
        return Err(
            "Complete or clear ambiguous active outer tails before proposing assembly primers."
                .into(),
        );
    }
    Ok(tail)
}
fn active_tails(
    row: &crate::cloning_methods_model::OverlapInput,
) -> Result<(String, String), String> {
    match row.preparation{
        OverlapPreparation::WholeLinear=>Ok((String::new(),String::new())),
        OverlapPreparation::ExactPcr|OverlapPreparation::MismatchPcr=>Ok((normalized_tail(&row.forward_tail)?,normalized_tail(&row.reverse_tail)?)),
        OverlapPreparation::OverlapExtensionPcr=>match &row.overlap_extension_pcr{
            None=>Ok((String::new(),String::new())),
            Some(raw)=>{raw.validate_raw()?;Ok((normalized_tail(&raw.outer.forward.tail)?,normalized_tail(&raw.outer.reverse.tail)?))},
        },
        _=>Err("This active preparation has a digest, source declaration or paired-oligo route. Keep its manual workflow, or explicitly choose a plain PCR preparation before proposing all overlap primers.".into()),
    }
}
fn tail_guard(draft: &Draft, tails: &[native::TailDeclaration]) -> Result<(), String> {
    let roles = draft.methods.role_ids(draft.slots.len())?;
    if tails.len() != roles.len() {
        return Err("The native tail declaration omitted a source contribution.".into());
    }
    for (index, ((role, row), proposed)) in
        roles.iter().zip(&draft.methods.rows).zip(tails).enumerate()
    {
        if role.role != proposed.role || role.input_id != proposed.input_id {
            return Err("The native tail order differs from the pinned source roles.".into());
        }
        let (forward, reverse) = active_tails(row)?;
        for (name, existing, derived) in [
            ("forward", forward, proposed.forward_tail.as_str()),
            ("reverse", reverse, proposed.reverse_tail.as_str()),
        ] {
            if !existing.is_empty() && existing != derived {
                return Err(format!(
                    "Fragment {} has an explicit {name} outer tail that differs from the assembly-derived tail. Keep its manual preparation or review and clear that tail before proposing all primers.",
                    index + 1
                ));
            }
        }
    }
    Ok(())
}
pub fn request(
    draft: &Draft,
    intent: &IntendedReference,
    current_sources: &[EventBinding],
    options: &Options,
    workspace: Option<&Workspace>,
) -> Result<(Capture, Vec<PlanInput>, native::Request), String> {
    draft.validate()?;
    if !eligible(draft) {
        return Err("Whole-assembly proposals require 2–12 overlap contributions with plain PCR/whole-source preparations. Active digest, declared-chemistry and paired-oligo routes must be kept or changed explicitly.".into());
    }
    if current_sources.len() != draft.slots.len() {
        return Err("The exact original-source event roster changed.".into());
    }
    for (slot, event) in draft.slots.iter().zip(current_sources) {
        let source = slot.captured.as_ref().unwrap();
        if source.source.binding != *event || source.molecule.revision() != event.revision {
            return Err("An original source event differs from its pinned contribution.".into());
        }
    }
    if !options.preparation.outer_forward_tail.is_empty()
        || !options.preparation.outer_reverse_tail.is_empty()
        || !options.preparation.primer_options.forward_tail.is_empty()
        || !options.preparation.primer_options.reverse_tail.is_empty()
    {
        return Err("Whole-assembly tails are derived from declared junction inserts and edited source cores; global tail overrides must be empty.".into());
    }
    crate::assembly_alternatives::require_active_choice_match(draft)?;
    if let Some(workspace) = workspace {
        attest(draft, workspace)?;
    }
    let (inputs, base, _) = assembly_draft_model::inputs_request(draft, workspace)?;
    let (intent_kind, intent_context_revision, declared) =
        independent(draft, intent, &inputs, base, workspace)?;
    let native = native::Request {
        intent: declared,
        overlap: options.overlap,
        options: options.preparation.clone(),
        profile_options: draft.methods.profile_options()?,
    };
    native::preflight(&inputs, &native, workspace).map_err(|e| e.to_string())?;
    if native.intent.base.base.material_bindings.is_empty() {
        tail_guard(
            draft,
            &native::expected_tails(&inputs, &native, workspace).map_err(|e| e.to_string())?,
        )?;
    }
    let mut capture = Capture {
        context_revision: String::new(),
        sources: current_sources.to_vec(),
        intent_kind,
        intent_context_revision,
        target_revision: intent.target().revision(),
        draft_revision: receipt(&context_draft(draft.clone()), MAX_CONTEXT)?,
        options_revision: receipt(options, MAX_REPORT)?,
        design_revision: native::design_revision(&inputs, &native, workspace)
            .map_err(|e| e.to_string())?,
    };
    capture.context_revision = receipt(&(ENGINE_VERSION, &capture), MAX_CONTEXT)?;
    Ok((capture, inputs, native))
}
pub fn from_design(
    draft: &Draft,
    intent: &IntendedReference,
    current_sources: &[EventBinding],
    options: &Options,
    report: native::Report,
    workspace: Option<&Workspace>,
) -> Result<Proposal, String> {
    let worker_receipt = receipt(&report, MAX_REPORT)?;
    let (capture, inputs, request) = request(draft, intent, current_sources, options, workspace)?;
    report
        .validate_binding(&capture.design_revision)
        .map_err(|e| e.to_string())?;
    if report.engine_version != native::ENGINE_VERSION
        || report.source_revisions
            != inputs
                .iter()
                .map(|input| input.revision.clone())
                .collect::<Vec<_>>()
        || receipt(&report.request, MAX_CONTEXT)? != receipt(&request, MAX_CONTEXT)?
    {
        return Err(
            "The whole-assembly worker result differs from its original source/request pins."
                .into(),
        );
    }
    // The native API normalizes all kinds to the general edit vocabulary. Keep
    // that receipt intact and validate it separately; metadata from the original
    // SNP/plain/junction target is never retagged as a native verification receipt.
    let canonical_context =
        assembly_edit_intent::intent_revision(&inputs, &request.intent, workspace)
            .map_err(|e| e.to_string())?;
    report
        .intent
        .validate_sources(&inputs, &canonical_context, workspace)
        .map_err(|e| e.to_string())?;
    if report.intent.target.sequence != intent.target().sequence
        || report.intent.target.topology != intent.target().topology
    {
        return Err("The proposed recipe differs from the frozen independent intended DNA.".into());
    }
    report
        .validate_sources(&inputs, &capture.design_revision, workspace)
        .map_err(|e| e.to_string())?;
    tail_guard(draft, &report.tails)?;
    let checked = crate::assembly_product_model::Snapshot::new(
        report.preview.compiled.plan.clone(),
        report.intent.target.clone(),
        report.preview.verification.clone(),
    );
    if let Some(error) = checked.prefix_error {
        return Err(format!(
            "The whole-assembly diagnostic receipt is unbound: {error}"
        ));
    }
    let recipe = report.recipe.clone();
    let verification = report.preview.verification.clone();
    let warnings = report.warnings.clone();
    let public_receipt = receipt(&(&capture, &recipe, &verification, &warnings), MAX_CONTEXT)?;
    Ok(Proposal {
        capture,
        recipe,
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
fn raw_pair(pair: &collagene_core::overlap_extension_pcr::PrimerPair) -> ExtensionPrimerPairInput {
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
        current_sources: &[EventBinding],
        options: &Options,
        workspace: Option<&Workspace>,
    ) -> Result<(), String> {
        if receipt(self.worker.as_ref(), MAX_REPORT)? != self.worker_receipt
            || receipt(
                &(
                    &self.capture,
                    &self.recipe,
                    &self.verification,
                    &self.warnings,
                ),
                MAX_CONTEXT,
            )? != self.public_receipt
        {
            return Err("The retained whole-assembly primer proposal changed.".into());
        }
        let (capture, _, request) = request(draft, intent, current_sources, options, workspace)?;
        if capture != self.capture
            || receipt(&self.worker.request, MAX_CONTEXT)? != receipt(&request, MAX_CONTEXT)?
        {
            return Err("The whole-assembly proposal is stale for the current draft, intent, sources or settings.".into());
        }
        self.worker
            .validate_binding(&capture.design_revision)
            .map_err(|e| e.to_string())
    }
    pub fn apply(
        &self,
        draft: &Draft,
        intent: &IntendedReference,
        current_sources: &[EventBinding],
        options: &Options,
        workspace: Option<&Workspace>,
    ) -> Result<Draft, String> {
        use collagene_core::assembly_recipe::{Preparation, RunPreparation};
        self.validate_binding(draft, intent, current_sources, options, workspace)?;
        if self.verification.status == collagene_core::verification::Status::Invalid {
            return Err("The proposed full recipe is inconsistent with its declared reference. Review its findings before changing preparations.".into());
        }
        if self.recipe.fragments.len() != draft.methods.rows.len() {
            return Err("The proposed recipe source count changed.".into());
        }
        let roles = draft.methods.role_ids(draft.slots.len())?;
        let mut next = draft.clone();
        for (index, fragment) in self.recipe.fragments.iter().enumerate() {
            if fragment.role != roles[index].role || fragment.input_id != roles[index].input_id {
                return Err(
                    "The proposed row differs from its original source role or orientation.".into(),
                );
            }
            let prep = match &fragment.preparation {
                Preparation::FromRun {
                    origin,
                    preparation,
                } => {
                    let saved = draft.methods.rows[index]
                        .run_origin
                        .as_ref()
                        .ok_or("An unrequested saved origin appeared in the proposed recipe.")?;
                    if receipt(
                        &collagene_core::run_material::MaterialOrigin::Run(origin.clone()),
                        MAX_CONTEXT,
                    )? != receipt(saved.as_ref(), MAX_CONTEXT)?
                    {
                        return Err("The proposed recipe changed a saved material origin.".into());
                    }
                    match preparation {
                        RunPreparation::ExactPcr {
                            template_strand,
                            forward,
                            reverse,
                        } => Preparation::ExactPcr {
                            template_strand: *template_strand,
                            forward: forward.clone(),
                            reverse: reverse.clone(),
                        },
                        RunPreparation::OverlapExtensionPcr { request } => {
                            Preparation::OverlapExtensionPcr {
                                request: request.clone(),
                            }
                        }
                        _ => {
                            return Err(
                                "The proposed preparation is outside the supported PCR vocabulary."
                                    .into(),
                            );
                        }
                    }
                }
                Preparation::FromDerived {
                    origin,
                    preparation,
                } => {
                    let saved = draft.methods.rows[index]
                        .run_origin
                        .as_ref()
                        .ok_or("An unrequested saved origin appeared in the proposed recipe.")?;
                    if receipt(
                        &collagene_core::run_material::MaterialOrigin::Derived(origin.clone()),
                        MAX_CONTEXT,
                    )? != receipt(saved.as_ref(), MAX_CONTEXT)?
                    {
                        return Err("The proposed recipe changed a saved material origin.".into());
                    }
                    match preparation {
                        RunPreparation::ExactPcr {
                            template_strand,
                            forward,
                            reverse,
                        } => Preparation::ExactPcr {
                            template_strand: *template_strand,
                            forward: forward.clone(),
                            reverse: reverse.clone(),
                        },
                        RunPreparation::OverlapExtensionPcr { request } => {
                            Preparation::OverlapExtensionPcr {
                                request: request.clone(),
                            }
                        }
                        _ => {
                            return Err(
                                "The proposed preparation is outside the supported PCR vocabulary."
                                    .into(),
                            );
                        }
                    }
                }
                other => {
                    if draft.methods.rows[index].run_origin.is_some() {
                        return Err(
                            "The native whole-assembly recipe dropped a saved material origin."
                                .into(),
                        );
                    }
                    other.clone()
                }
            };
            let edited = draft.edits.iter().any(|e| e.role == fragment.role)
                || draft.substitutions.iter().any(|e| e.role == fragment.role);
            let row = &mut next.methods.rows[index];
            row.reverse = fragment.reverse;
            match prep {
                Preparation::ExactPcr {
                    template_strand,
                    forward,
                    reverse,
                } if !edited => {
                    row.template_strand = template_strand;
                    row.preparation = OverlapPreparation::ExactPcr;
                    row.forward_name = forward.name;
                    row.forward = forward.annealing.as_str().into();
                    row.forward_tail = forward.tail;
                    row.reverse_name = reverse.name;
                    row.reverse_primer = reverse.annealing.as_str().into();
                    row.reverse_tail = reverse.tail;
                }
                Preparation::OverlapExtensionPcr { request } if edited => {
                    let o = &request.extension;
                    let raw = OverlapExtensionPcrInput {
                        left: raw_pair(&request.left),
                        right: raw_pair(&request.right),
                        outer: raw_pair(&request.outer),
                        overlap: o.overlap.to_string(),
                        min_paired_bases: o.min_paired_bases.to_string(),
                        max_comparisons: o.max_comparisons.to_string(),
                        max_candidates: o.max_candidates.to_string(),
                    };
                    raw.request()?;
                    row.preparation = OverlapPreparation::OverlapExtensionPcr;
                    row.overlap_extension_pcr = Some(raw);
                }
                _ => return Err(
                    "Edited roles need six-oligo OE preparation and unedited roles need exact PCR."
                        .into(),
                ),
            }
        }
        next.methods.overlaps = std::iter::repeat_n(
            options.overlap.to_string(),
            draft.slots.len() - 1 + usize::from(draft.circular),
        )
        .collect::<Vec<_>>()
        .join(",");
        next.sync_rows();
        next.validate()?;
        if next.alternatives.is_some() {
            next = crate::assembly_alternatives::save_active(&next)?;
        }
        Ok(next)
    }
    pub fn rows(&self) -> &[native::RowDesign] {
        &self.worker.rows
    }
    pub fn source_frame(&self, input_id: &str) -> Option<&native::SourceFrame> {
        self.worker
            .source_frames
            .iter()
            .find(|frame| frame.input_id == input_id)
    }
    pub fn tails(&self) -> &[native::TailDeclaration] {
        &self.worker.tails
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        assembly_draft_model::CapturedSource, assembly_edit_model, source_identity::SourceInfo,
    };
    use collagene_core::{Alphabet, Topology, chemistry::EndRepairMethod, verification::Status};
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
    fn intent(draft: &Draft, workspace: &Workspace) -> IntendedReference {
        let (inputs, base, _) =
            assembly_draft_model::inputs_request(draft, Some(workspace)).unwrap();
        IntendedReference::Edits(
            assembly_edit_intent::build(
                &inputs,
                &assembly_draft_model::edit_request(draft, base),
                Some(workspace),
            )
            .unwrap(),
        )
    }
    fn fixture() -> (Draft, Workspace, Options) {
        let vector = Molecule::dna("Vector", "CCGGAATTGGCCATACCGTA", Topology::Linear).unwrap();
        let insert = Molecule::dna("Original edited stock", ORIGINAL, Topology::Linear).unwrap();
        let frame = Workspace::from_documents(
            "Original exact events",
            vec![vector.clone(), insert.clone()],
            0,
        )
        .unwrap();
        let mut draft = Draft::default();
        draft.circular = false;
        draft.active_slot = 1;
        draft.methods.mode = Mode::Overlap;
        draft.slots[0].captured = Some(captured(vector, event(&frame, 0)));
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 20;
        }
        draft.slots[1].captured = Some(captured(insert, event(&frame, 1)));
        {
            let collagene_core::assembly_intent::Selector::ReferenceRange { end: value, .. } =
                &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 40;
        }
        draft.methods.rows[0].preparation = OverlapPreparation::WholeLinear;
        draft.methods.rows[0].forward_tail = "inactive non-DNA".into();
        draft.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcr;
        draft.methods.rows[1].overlap_extension_pcr = None;
        draft.methods.rows[1].repair = Some(EndRepairMethod::Blunt);
        draft.methods.rows[1].treatment = crate::cloning_methods_model::Treatment::Phosphorylate;
        draft.methods.rows[1].fragment = "7".into();
        draft.sync_rows();
        let draft = assembly_edit_model::declare(&draft, 1, "17", "21", "TT", None).unwrap();
        let mut options = Options::default();
        let p = &mut options.preparation.primer_options;
        p.min_length = 8;
        p.max_length = 8;
        p.max_pairs = 1;
        p.max_tm_difference_c = 100.;
        p.min_gc_fraction = 0.;
        p.max_gc_fraction = 1.;
        p.max_homopolymer = 200;
        p.max_three_prime_complementarity = 200;
        options.preparation.homology_flank = 4;
        options.preparation.min_paired_bases = 8;
        (draft, frame, options)
    }
    fn propose(
        draft: &Draft,
        intended: &IntendedReference,
        frame: &Workspace,
        options: &Options,
    ) -> Proposal {
        let events = attest(draft, frame).unwrap();
        let (_, inputs, native) = request(draft, intended, &events, options, Some(frame)).unwrap();
        let report = native::design(&inputs, &native, Some(frame)).unwrap();
        from_design(draft, intended, &events, options, report, Some(frame)).unwrap()
    }
    #[test]
    fn explicit_library_construct_proposal_changes_only_its_overrides_and_survives_native_archive()
    {
        let (mut original, mut frame, options) = fixture();
        original.edits.clear();
        original.substitutions.clear();
        let first = crate::assembly_alternatives::enable(&original, "edit-primers").unwrap();
        let first_id = first.alternatives.as_ref().unwrap().active_row.clone();
        let choices = first.alternatives.as_ref().unwrap().rows[0].choices.clone();
        let mut second = crate::assembly_alternatives::add_row(&first, choices).unwrap();
        second = assembly_edit_model::declare(&second, 1, "17", "21", "TT", None).unwrap();
        let first_saved = second.alternatives.as_ref().unwrap().rows[0].clone();
        let frozen = intent(&second, &frame);
        let target = frozen.target().revision();
        let proposal = propose(&second, &frozen, &frame, &options);
        assert_eq!(
            proposal.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            proposal.verification.findings
        );
        let applied = proposal
            .apply(
                &second,
                &frozen,
                &proposal.capture.sources,
                &options,
                Some(&frame),
            )
            .unwrap();
        assert!(applied.alternatives.as_ref().unwrap().rows[0] == first_saved);
        assert_eq!(frozen.target().revision(), target);
        assert!(matches!(
            applied.methods.rows[1].preparation,
            OverlapPreparation::OverlapExtensionPcr
        ));
        crate::assembly_edited_primer_commit::validate(&second, &applied, &second).unwrap();
        let first_again = crate::assembly_alternatives::open_row(&applied, &first_id).unwrap();
        assert!(
            proposal
                .apply(
                    &first_again,
                    &frozen,
                    &proposal.capture.sources,
                    &options,
                    Some(&frame)
                )
                .is_err()
        );
        let bytes = assembly_draft_model::encode(&applied).unwrap();
        let hash = blake3::hash(&bytes).to_hex().to_string();
        frame.attachments.insert(
            hash.clone(),
            collagene_core::workspace::WorkspaceAttachment {
                kind: assembly_draft_model::ATTACHMENT_KIND.into(),
                name: "Row declarations".into(),
                media_type: "application/json".into(),
                bytes,
            },
        );
        let reopened = collagene_core::workspace::load_workspace(
            &collagene_core::workspace::save_workspace(&frame).unwrap(),
        )
        .unwrap();
        let draft = assembly_draft_model::decode(&reopened.attachments[&hash].bytes).unwrap();
        assert!(draft == applied);
        assert!(draft.alternatives.as_ref().unwrap().rows[0] == first_saved);
        let mut changed = second.clone();
        changed.edits[0].replacement = "AA".into();
        assert!(
            proposal
                .apply(
                    &changed,
                    &frozen,
                    &proposal.capture.sources,
                    &options,
                    Some(&frame)
                )
                .is_err()
        );
    }
    #[test]
    fn whole_proposal_uses_edited_downstream_homology_and_applies_only_active_fields() {
        let (draft, frame, options) = fixture();
        let intended = intent(&draft, &frame);
        let original_target = intended.target().revision();
        let p = propose(&draft, &intended, &frame, &options);
        assert_eq!(
            p.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            p.verification.findings
        );
        assert_eq!(p.tails()[0].homology, "AGTCCGATGCTACGTATTTT");
        assert_ne!(p.tails()[0].homology, &ORIGINAL[..20]);
        assert_eq!(p.rows()[0].designs.len(), 1);
        assert_eq!(p.rows()[1].designs.len(), 3);
        let after = p
            .apply(
                &draft,
                &intended,
                &p.capture.sources,
                &options,
                Some(&frame),
            )
            .unwrap();
        assert_eq!(
            crate::assembly_edited_primer_commit::validate(&draft, &after, &draft).unwrap(),
            2
        );
        assert!(after.methods.rows[0].preparation == OverlapPreparation::ExactPcr);
        assert!(after.methods.rows[1].preparation == OverlapPreparation::OverlapExtensionPcr);
        let raw = after.methods.rows[1]
            .overlap_extension_pcr
            .as_ref()
            .unwrap();
        let typed = raw.request().unwrap();
        assert!(!typed.left.forward.annealing.is_empty());
        assert!(!typed.left.reverse.annealing.is_empty());
        assert!(!typed.right.forward.annealing.is_empty());
        assert!(!typed.right.reverse.annealing.is_empty());
        assert!(!typed.outer.forward.annealing.is_empty());
        assert!(!typed.outer.reverse.annealing.is_empty());
        assert_eq!(after.methods.overlaps, "20");
        assert_eq!(after.methods.rows[1].repair, draft.methods.rows[1].repair);
        assert!(after.methods.rows[1].treatment == draft.methods.rows[1].treatment);
        assert_eq!(after.methods.rows[1].fragment, "7");
        assert!(after.slots == draft.slots);
        assert_eq!(intended.target().revision(), original_target);
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
    }
    #[test]
    fn contribution_navigation_preserves_whole_proposal_but_scientific_changes_stale_it() {
        let (mut draft, frame, options) = fixture();
        draft.active_slot = 0;
        let intended = intent(&draft, &frame);
        let p = propose(&draft, &intended, &frame, &options);
        let mut navigated = draft.clone();
        navigated.active_slot = 1;
        p.validate_binding(
            &navigated,
            &intended,
            &p.capture.sources,
            &options,
            Some(&frame),
        )
        .unwrap();
        let after = p
            .apply(
                &navigated,
                &intended,
                &p.capture.sources,
                &options,
                Some(&frame),
            )
            .unwrap();
        assert_eq!(after.active_slot, 1);
        crate::assembly_edited_primer_commit::validate(&navigated, &after, &navigated).unwrap();
        navigated.methods.rows[0]
            .forward_name
            .push_str(" changed declaration");
        assert!(
            p.apply(
                &navigated,
                &intended,
                &p.capture.sources,
                &options,
                Some(&frame)
            )
            .is_err()
        );
    }
    #[test]
    fn guards_preserve_explicit_tails_source_events_and_sealed_public_receipts() {
        let (draft, frame, options) = fixture();
        let intended = intent(&draft, &frame);
        let p = propose(&draft, &intended, &frame, &options);
        for prep in [
            OverlapPreparation::Digest,
            OverlapPreparation::CatalogDigest,
            OverlapPreparation::PcrCatalogDigest,
            OverlapPreparation::OverlapExtensionPcrCatalogDigest,
            OverlapPreparation::DeclaredDuplex,
            OverlapPreparation::AnnealedOligos,
        ] {
            let mut changed = draft.clone();
            changed.methods.rows[0].preparation = prep;
            changed.sync_rows();
            assert!(
                request(
                    &changed,
                    &intended,
                    &p.capture.sources,
                    &options,
                    Some(&frame)
                )
                .is_err()
            );
        }
        let mut changed = draft.clone();
        changed.methods.rows[0].forward_name.push_str("stale");
        assert!(
            p.apply(
                &changed,
                &intended,
                &p.capture.sources,
                &options,
                Some(&frame)
            )
            .is_err()
        );
        let mut events = p.capture.sources.clone();
        events[0].event_id.push('x');
        assert!(
            p.apply(&draft, &intended, &events, &options, Some(&frame))
                .is_err()
        );
        let mut other_options = options.clone();
        other_options.overlap -= 1;
        assert!(
            p.apply(
                &draft,
                &intended,
                &p.capture.sources,
                &other_options,
                Some(&frame)
            )
            .is_err()
        );
        let mut changed = p.clone();
        changed.warnings.push("tampered".into());
        assert!(
            changed
                .apply(
                    &draft,
                    &intended,
                    &p.capture.sources,
                    &options,
                    Some(&frame)
                )
                .is_err()
        );
        let mut changed = p.clone();
        changed.recipe.fragments[0].reverse = true;
        assert!(
            changed
                .apply(
                    &draft,
                    &intended,
                    &p.capture.sources,
                    &options,
                    Some(&frame)
                )
                .is_err()
        );
        let mut bad_frame = frame.clone();
        bad_frame.history[0].operation.push_str("forged");
        assert!(attest(&draft, &bad_frame).is_err());
        let mut explicit = draft.clone();
        explicit.methods.rows[0].preparation = OverlapPreparation::ExactPcr;
        explicit.methods.rows[0].forward_tail = "A".into();
        explicit.sync_rows();
        assert!(
            request(
                &explicit,
                &intent(&explicit, &frame),
                &p.capture.sources,
                &options,
                Some(&frame)
            )
            .is_err()
        );
        explicit.methods.rows[0].forward_tail.clear();
        explicit.methods.rows[0].reverse_tail = p.tails()[0]
            .reverse_tail
            .to_lowercase()
            .chars()
            .map(|c| format!("{c} "))
            .collect();
        assert!(
            request(
                &explicit,
                &intent(&explicit, &frame),
                &p.capture.sources,
                &options,
                Some(&frame)
            )
            .is_ok()
        );
    }
    #[test]
    fn substitution_reference_keeps_original_metadata_receipt_while_canonical_check_is_separate() {
        let (mut draft, frame, options) = fixture();
        draft.edits.clear();
        draft
            .substitutions
            .push(collagene_core::assembly_substitution_intent::Substitution {
                role: draft.methods.role_ids(2).unwrap()[1].role.clone(),
                offset: 16,
                expected: "G".into(),
                replacement: "T".into(),
            });
        draft.sync_rows();
        let (inputs, base, _) = assembly_draft_model::inputs_request(&draft, None).unwrap();
        let declared = assembly_draft_model::edit_request(&draft, base);
        let intended = IntendedReference::Substitutions(
            collagene_core::assembly_substitution_intent::build(&inputs, &declared.base, None)
                .unwrap(),
        );
        let frozen = intended.target().revision();
        let p = propose(&draft, &intended, &frame, &options);
        assert_eq!(p.capture.intent_kind, "substitutions");
        assert_eq!(p.capture.target_revision, frozen);
        assert_ne!(p.worker.intent.target.revision(), frozen);
        let after = p
            .apply(
                &draft,
                &intended,
                &p.capture.sources,
                &options,
                Some(&frame),
            )
            .unwrap();
        assert_eq!(after.substitutions, draft.substitutions);
        assert_eq!(intended.target().revision(), frozen);
        crate::assembly_edited_primer_commit::validate(&draft, &after, &draft).unwrap();
        let mut forged = p.worker.as_ref().clone();
        forged.intent.target.name.push_str("forged");
        assert!(
            from_design(
                &draft,
                &intended,
                &p.capture.sources,
                &options,
                forged,
                Some(&frame)
            )
            .is_err()
        );
    }
    #[test]
    fn saved_origin_is_kept_and_replayed_in_complete_assembly_check() {
        use collagene_core::{
            run_lineage::InputEventBinding,
            verification::{DesignPlan, Operation, Source},
        };
        let (mut draft, _, options) = fixture();
        let vector = draft.slots[0].captured.as_ref().unwrap().molecule.clone();
        let stock = draft.slots[1].captured.as_ref().unwrap().molecule.clone();
        let target = Molecule::dna(
            "Independent stock amplification requirement",
            ORIGINAL,
            Topology::Linear,
        )
        .unwrap();
        let frame = Workspace::from_documents(
            "Saved stock",
            vec![vector, stock.clone(), target.clone()],
            0,
        )
        .unwrap();
        let plan = DesignPlan {
            name: "Actual stock PCR".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "stock".into(),
                revision: stock.revision(),
                molecule: stock,
            }],
            steps: vec![Operation::Pcr {
                template_strand: collagene_core::chemistry::Strand::Top,
                template: Source::Input("stock".into()),
                forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            }],
        };
        let saved = collagene_core::cloning_outputs::finalize(
            &frame,
            &plan,
            event(&frame, 2),
            vec![InputEventBinding {
                input_id: "stock".into(),
                event: event(&frame, 1),
            }],
            &collagene_core::cloning_outputs::Options {
                name: "Saved original stock".into(),
                save_primers: false,
                show_prepared: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let selection = collagene_core::run_material::Selection {
            run_id: saved.receipt.run.run_id.clone(),
            step_index: 0,
            event: saved.receipt.construct.clone(),
        };
        let origin = collagene_core::run_material::capsule(&saved.workspace, &selection).unwrap();
        draft.slots[0].captured = Some(captured(
            saved.workspace.revisions[&frame.documents[0].head].clone(),
            event(&frame, 0),
        ));
        draft.slots[1].captured = Some(captured(
            saved.workspace.revisions[&selection.event.revision].clone(),
            selection.event.clone(),
        ));
        draft.methods.rows[1].run_origin = Some(
            crate::run_origin_model::PinnedOrigin::from_material(
                collagene_core::run_material::MaterialOrigin::Run(origin),
            )
            .unwrap(),
        );
        draft.slots[1].material_origin = draft.methods.rows[1].run_origin.clone();
        draft.sync_rows();
        let intended = intent(&draft, &saved.workspace);
        let p = propose(&draft, &intended, &saved.workspace, &options);
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
                &p.capture.sources,
                &options,
                Some(&saved.workspace),
            )
            .unwrap();
        assert_eq!(
            after.methods.rows[1].run_origin,
            draft.methods.rows[1].run_origin
        );
        assert!(matches!(
            p.recipe.fragments[1].preparation,
            collagene_core::assembly_recipe::Preparation::FromRun { .. }
        ));
        crate::assembly_edited_primer_commit::validate(&draft, &after, &draft).unwrap();
    }
}
