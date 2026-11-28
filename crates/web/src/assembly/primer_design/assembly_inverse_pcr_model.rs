//! Proposal adapter for one independently declared edit in one whole circular
//! source. No target is inferred from primers or changed to fit a PCR product.
//!
//! This pure adapter compares the caller's current source-event identity with the
//! pinned draft. The host must still attest that exact event in its workspace and
//! use ordinary atomic assembly staging before installing the returned candidate.
use crate::assembly::draft::model::Draft;
use collagene_core::{
    cloning::{
        assembly::{
            edit_intent as assembly_edit_intent,
            substitution_intent as assembly_substitution_intent,
        },
        pcr::{
            mutagenesis::{self, InversePcrDesign, InversePcrOptions, InversePcrPreparation},
            mutation_intent,
        },
    },
    documents::{transfer::Candidate, workspace::Workspace},
    material::lineage::EventBinding,
    planning::plan::{Operation, PlanInput, Source},
    verification::evidence::VerificationReport,
};
use collagene_model::{Molecule, Topology};
use serde::{Deserialize, Serialize};

const ENGINE_VERSION: &str = "assembly-inverse-pcr-proposal-1";

mod sealed {
    pub trait Sealed {}
    impl Sealed for collagene_core::cloning::assembly::edit_intent::Report {}
    impl Sealed for collagene_core::cloning::assembly::substitution_intent::Report {}
    impl Sealed for crate::assembly::review::assembly_review::IntendedReference {}
}

/// The two independently source-bound intent kinds supported by this bridge.
/// A frozen target is never converted between report kinds or rewritten to match
/// prepared DNA. This sealed interface also lets the dock retain its typed enum.
pub trait IntendedEdit: sealed::Sealed {
    fn intended_target(&self) -> &Molecule;
    fn validate_original(
        &self,
        inputs: &[PlanInput],
        base: &assembly_substitution_intent::Request,
        edits: &[assembly_edit_intent::Edit],
        workspace: Option<&Workspace>,
    ) -> Result<String, String>;
}
impl IntendedEdit for assembly_edit_intent::Report {
    fn intended_target(&self) -> &Molecule {
        &self.target
    }
    fn validate_original(
        &self,
        inputs: &[PlanInput],
        base: &assembly_substitution_intent::Request,
        edits: &[assembly_edit_intent::Edit],
        workspace: Option<&Workspace>,
    ) -> Result<String, String> {
        if edits.len() != 1 || !base.substitutions.is_empty() {
            return Err(
                "The general-edit report does not describe the declared single edit.".into(),
            );
        }
        let request = assembly_edit_intent::Request {
            base: base.clone(),
            edits: edits.to_vec(),
        };
        let expected = assembly_edit_intent::intent_revision(inputs, &request, workspace)
            .map_err(|e| e.to_string())?;
        self.validate_sources(inputs, &expected, workspace)
            .map_err(|e| e.to_string())?;
        Ok(expected)
    }
}
impl IntendedEdit for assembly_substitution_intent::Report {
    fn intended_target(&self) -> &Molecule {
        &self.target
    }
    fn validate_original(
        &self,
        inputs: &[PlanInput],
        base: &assembly_substitution_intent::Request,
        edits: &[assembly_edit_intent::Edit],
        workspace: Option<&Workspace>,
    ) -> Result<String, String> {
        if !edits.is_empty() || base.substitutions.len() != 1 {
            return Err(
                "The substitution report does not describe the declared single-base edit.".into(),
            );
        }
        let expected = assembly_substitution_intent::intent_revision(inputs, base, workspace)
            .map_err(|e| e.to_string())?;
        self.validate_binding(&expected)
            .map_err(|e| e.to_string())?;
        let fresh = assembly_substitution_intent::build(inputs, base, workspace)
            .map_err(|e| e.to_string())?;
        if receipt(self)? != receipt(&fresh)? {
            return Err(
                "The substitution intent differs from its exact original source declaration."
                    .into(),
            );
        }
        Ok(expected)
    }
}
impl IntendedEdit for crate::assembly::review::assembly_review::IntendedReference {
    fn intended_target(&self) -> &Molecule {
        self.target()
    }
    fn validate_original(
        &self,
        inputs: &[PlanInput],
        base: &assembly_substitution_intent::Request,
        edits: &[assembly_edit_intent::Edit],
        workspace: Option<&Workspace>,
    ) -> Result<String, String> {
        match self {
            Self::Edits(intent) => intent.validate_original(inputs, base, edits, workspace),
            Self::Substitutions(intent) => intent.validate_original(inputs, base, edits, workspace),
            _ => Err("Declare one internal edit before preparing it by inverse PCR.".into()),
        }
    }
}
fn literal_edit(draft: &Draft) -> Result<assembly_edit_intent::Edit, String> {
    match (draft.edits.as_slice(),draft.substitutions.as_slice()) {
        ([edit],[])=>Ok(edit.clone()),
        ([],[edit])=>Ok(assembly_edit_intent::Edit {
            role:edit.role.clone(),start:edit.offset,
            end:edit.offset.checked_add(1).ok_or("Substitution coordinates overflow.")?,
            expected:edit.expected.clone(),replacement:edit.replacement.clone(),
        }),
        _=>Err("Inverse-PCR proposals require exactly one internal edit or one point substitution; mixed or multiple declarations need a staged preparation.".into()),
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Options {
    pub forward_length: usize,
    pub reverse_length: usize,
    pub reverse_tail_bases: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            forward_length: 20,
            reverse_length: 20,
            reverse_tail_bases: 0,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Capture {
    pub context_revision: String,
    pub source: EventBinding,
    pub input_id: String,
    pub origin: Option<crate::history::run_origin_model::PinnedOrigin>,
    pub intent_context_revision: String,
    pub target_revision: String,
    pub options: Options,
}
#[derive(Clone, Debug, Serialize)]
pub struct Proposal {
    pub capture: Capture,
    pub preparation: InversePcrPreparation,
    pub candidate: Candidate,
    /// Original worker-plan diagnostic only. Its fixed "template" identity and
    /// product tags remain unchanged; this is not a receipt for the rebound
    /// candidate. Opening/finalizing that candidate requires fresh verification.
    pub verification: VerificationReport,
    worker_design: InversePcrDesign,
    worker_receipt: String,
}

fn context(
    event: &EventBinding,
    input_id: &str,
    intent: &str,
    target: &str,
    options: &Options,
    origin: Option<&crate::history::run_origin_model::PinnedOrigin>,
) -> Result<String, String> {
    receipt(&(
        ENGINE_VERSION,
        event,
        input_id,
        intent,
        target,
        options,
        origin,
    ))
}

/// Rebuild independent intent from the original draft before preparing any
/// oligos. Accepting a current event pin does not forge historical workspace proof.
pub fn capture<I: IntendedEdit + ?Sized>(
    draft: &Draft,
    intent: &I,
    current_source: &EventBinding,
    frozen_target: &Molecule,
    options: &Options,
    workspace: Option<&Workspace>,
) -> Result<Capture, String> {
    draft.validate()?;
    collagene_model::bounds::model_preflight(&[frozen_target]).map_err(|e| e.to_string())?;
    frozen_target.validate().map_err(|e| e.to_string())?;
    if frozen_target.sequence.len() > 8 << 20 {
        return Err("The frozen target exceeds the ordinary assembly size bound.".into());
    }
    if draft.alternatives.is_some()
        || draft.slots.len() != 1
        || !draft.circular
        || !draft.junction_edits.is_empty()
    {
        return Err("Inverse-PCR proposals require one whole circular source and one independent internal edit or point substitution, without junction edits.".into());
    }
    let edit = literal_edit(draft)?;
    let slot = &draft.slots[0];
    let source = slot
        .captured
        .as_ref()
        .ok_or("Pin the original circular source first.")?;
    if source.source.binding != *current_source
        || source.source.binding.revision != source.molecule.revision()
    {
        return Err(
            "The current original-source event differs from the pinned contribution.".into(),
        );
    }
    if source.molecule.topology != Topology::Circular
        || slot.selector
            != (collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                start: 0,
                end: source.molecule.sequence.len(),
                reverse: false,
            })
        || slot.digest_geometry.is_some()
        || slot.anneal_bottom.is_some()
    {
        return Err("Inverse-PCR proposals currently require the complete forward original circular reference, without alternative source geometry. Saved native origins must match the exact source event.".into());
    }
    let origin = slot.material_origin.clone();
    if let Some(pin) = &origin {
        pin.validate_documentary()?;
        if pin.event() != current_source {
            return Err("Saved material belongs to another source event.".into());
        }
    }
    if !(8..=80).contains(&options.forward_length)
        || !(8..=80).contains(&options.reverse_length)
        || options.reverse_tail_bases > edit.replacement.len()
    {
        return Err("Choose 8–80 base primer domains and a valid insertion-tail split.".into());
    }
    let roles = draft.methods.role_ids(1)?;
    let role = &roles[0];
    let (inputs, base, _) = draft.inputs_request(workspace)?;
    let base = assembly_substitution_intent::Request {
        base,
        junction_edits: draft.junction_edits.clone(),
        substitutions: draft.substitutions.clone(),
    };
    // Preserve and rebuild the original report kind, including its full target
    // metadata. Normalizing the literal below is only for native PCR options.
    let expected = intent.validate_original(&inputs, &base, &draft.edits, workspace)?;
    if frozen_target.revision() != intent.intended_target().revision() {
        return Err("The frozen intended target changed; review the edit again.".into());
    }
    let independently_edited = mutation_intent::build(
        &source.molecule,
        &mutation_intent::Request {
            name: draft.name.clone(),
            edit_start: edit.start,
            edit_end: edit.end,
            inserted: edit.replacement.clone(),
        },
    )
    .map_err(|e| e.to_string())?;
    if independently_edited.target.sequence != frozen_target.sequence
        || independently_edited.target.topology != frozen_target.topology
    {
        return Err("The declared mutation does not describe the frozen dock target.".into());
    }
    let target_revision = frozen_target.revision();
    Ok(Capture {
        context_revision: context(
            current_source,
            &role.input_id,
            &expected,
            &target_revision,
            options,
            origin.as_ref(),
        )?,
        source: current_source.clone(),
        input_id: role.input_id.clone(),
        origin,
        intent_context_revision: expected,
        target_revision,
        options: options.clone(),
    })
}

/// Capture a bounded request for the existing inverse-PCR worker job. The
/// source-event guard is checked independently from the eventual response.
pub fn request<I: IntendedEdit + ?Sized>(
    draft: &Draft,
    intent: &I,
    current_source: &EventBinding,
    frozen_target: &Molecule,
    options: &Options,
    workspace: Option<&Workspace>,
) -> Result<(Capture, Molecule, InversePcrOptions), String> {
    let captured = capture(
        draft,
        intent,
        current_source,
        frozen_target,
        options,
        workspace,
    )?;
    let edit = literal_edit(draft)?;
    let native = InversePcrOptions {
        edit_start: edit.start,
        edit_end: edit.end,
        inserted: edit.replacement.clone(),
        forward_length: options.forward_length,
        reverse_length: options.reverse_length,
        reverse_tail_bases: options.reverse_tail_bases,
    };
    Ok((
        captured,
        draft.slots[0].captured.as_ref().unwrap().molecule.clone(),
        native,
    ))
}

struct ReceiptWriter {
    hash: blake3::Hasher,
    bytes: usize,
    max: usize,
}
impl std::io::Write for ReceiptWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|n| *n <= self.max)
            .ok_or_else(|| std::io::Error::other("inverse-PCR proposal receipt size"))?;
        self.hash.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn receipt_with_limit(value: &impl Serialize, max: usize) -> Result<String, String> {
    let mut writer = ReceiptWriter {
        hash: blake3::Hasher::new(),
        bytes: 0,
        max,
    };
    serde_json::to_writer(&mut writer, value).map_err(|e| e.to_string())?;
    Ok(writer.hash.finalize().to_hex().to_string())
}
fn receipt(value: &impl Serialize) -> Result<String, String> {
    receipt_with_limit(value, 128 << 20)
}
fn preparation_from_design(report: &InversePcrDesign) -> InversePcrPreparation {
    InversePcrPreparation {
        context_revision: report.context_revision.clone(),
        source_revision: report.source_revision.clone(),
        target_revision: report.target_revision.clone(),
        options: report.options.clone(),
        forward: report.forward.clone(),
        reverse: report.reverse.clone(),
        forward_location: report.forward_location.clone(),
        reverse_location: report.reverse_location.clone(),
        retained_template_length: report.retained_template_length,
        plan: report.plan.clone(),
        warnings: report.warnings.clone(),
    }
}

/// Validate an existing worker result without rerunning molecular reconstruction
/// on the UI thread. Native deterministic compilation checks the exact source,
/// target, options, ordered oligos and steps; the shared product-snapshot adapter
/// checks the original worker report's outer pins, engine and product tags.
///
/// Renaming an input changes its plan revision. Consequently the original worker
/// verification stays diagnostic evidence for the original "template" plan. It
/// is never rewritten or installed as evidence for the stable-role candidate.
pub fn from_design<I: IntendedEdit + ?Sized>(
    draft: &Draft,
    intent: &I,
    current_source: &EventBinding,
    frozen_target: &Molecule,
    options: &Options,
    report: InversePcrDesign,
    workspace: Option<&Workspace>,
) -> Result<Proposal, String> {
    // Bound rich response fields before cloning or checking their correspondence.
    let worker_receipt = receipt_with_limit(&report, 32 << 20)?;
    let (captured, source, native) = request(
        draft,
        intent,
        current_source,
        frozen_target,
        options,
        workspace,
    )?;
    let mut preparation = mutagenesis::compile_inverse_pcr(
        &source,
        frozen_target,
        captured.origin.as_ref().map(|p| p.as_ref()),
        &native,
    )
    .map_err(|e| e.to_string())?;
    if receipt(&preparation)? != receipt(&preparation_from_design(&report))? {
        return Err("The worker inverse-PCR plan, oligos or options differ from the pinned source declaration.".into());
    }
    let checked = crate::assembly::review::assembly_product_model::Snapshot::new(
        report.plan.clone(),
        frozen_target.clone(),
        report.verification.clone(),
    );
    if let Some(error) = checked.prefix_error {
        return Err(format!(
            "The original worker reconstruction receipt is unbound: {error}"
        ));
    }
    // Rebind only the input identifier; do not relabel worker receipt/product pins.
    if preparation.plan.inputs.len() != 1 || preparation.plan.inputs[0].id != "template" {
        return Err("The native inverse-PCR input shape changed; inspect the adapter.".into());
    }
    match preparation.plan.steps.first_mut() {
        Some(Operation::Pcr{template:Source::Input(id),..}) if captured.origin.is_none()&&id=="template" => *id=captured.input_id.clone(),
        Some(Operation::ReconstructRunProduct{template:Source::Input(id),..}|Operation::ReconstructMaterial{template:Source::Input(id),..}) if captured.origin.is_some()&&id=="template"=>*id=captured.input_id.clone(),
        _=>return Err("The native inverse-PCR preparation no longer starts with the expected original-source PCR.".into()),
    }
    preparation.plan.inputs[0].id = captured.input_id.clone();
    let candidate = Candidate::new(preparation.plan.clone());
    candidate.validate().map_err(|e| e.to_string())?;
    Ok(Proposal {
        capture: captured,
        preparation,
        candidate,
        verification: report.verification.clone(),
        worker_design: report,
        worker_receipt,
    })
}

impl Proposal {
    /// Recheck the capture, exact native compilation and retained worker receipt.
    /// This rejects edits to public proposal fields without rerunning physics.
    /// Host staging and fresh verification still authorize the rebound candidate.
    pub fn validate_binding<I: IntendedEdit + ?Sized>(
        &self,
        draft: &Draft,
        intent: &I,
        current_source: &EventBinding,
        frozen_target: &Molecule,
        options: &Options,
        workspace: Option<&Workspace>,
    ) -> Result<(), String> {
        if receipt_with_limit(&self.worker_design, 32 << 20)? != self.worker_receipt {
            return Err("The retained worker receipt changed.".into());
        }
        let fresh = from_design(
            draft,
            intent,
            current_source,
            frozen_target,
            options,
            self.worker_design.clone(),
            workspace,
        )?;
        if receipt(self)? != receipt(&fresh)? {
            return Err(
                "The inverse-PCR proposal differs from its source-bound worker receipt.".into(),
            );
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::draft::model::CapturedSource;
    use crate::assembly::edit::assembly_edit_model;
    use crate::cloning::preparation::cloning_methods_model::Mode;
    use crate::sources::identity::SourceInfo;
    use collagene_core::{
        documents::workspace::Workspace,
        sequence::analysis,
        verification::{self, evidence::Status},
    };
    use collagene_model::Alphabet;
    use std::sync::Arc;
    /// Native fixture/reference path: the expensive molecular reconstruction runs
    /// once. The web UI uses request → existing worker → from_design instead.
    fn propose<I: IntendedEdit + ?Sized>(
        draft: &Draft,
        intent: &I,
        current_source: &EventBinding,
        frozen_target: &Molecule,
        options: &Options,
        workspace: Option<&Workspace>,
    ) -> Result<Proposal, String> {
        let (captured, source, native) = request(
            draft,
            intent,
            current_source,
            frozen_target,
            options,
            workspace,
        )?;
        let report = mutagenesis::design_inverse_pcr(
            &source,
            frozen_target,
            captured.origin.as_ref().map(|p| p.as_ref()),
            &native,
        )
        .map_err(|e| e.to_string())?;
        from_design(
            draft,
            intent,
            current_source,
            frozen_target,
            options,
            report,
            workspace,
        )
    }
    fn fixture(
        start: &str,
        end: &str,
        inserted: &str,
    ) -> (Draft, assembly_edit_intent::Report, EventBinding, Options) {
        let source = Molecule::dna(
            "Original parent",
            "AGTCCGATGCTACGTAGGACTTCAGACCTGTA",
            Topology::Circular,
        )
        .unwrap();
        let frame = Workspace::from_documents("Original events", vec![source.clone()], 0).unwrap();
        let doc = &frame.documents[0];
        let event = EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        };
        let capture = Arc::new(CapturedSource {
            source: SourceInfo {
                binding: event.clone(),
                name: source.name.clone(),
                length: source.sequence.len(),
                alphabet: Alphabet::Dna,
                topology: Topology::Circular,
            },
            molecule: source.clone(),
        });
        let mut draft = Draft::default();
        draft.slots.truncate(1);
        draft.methods.rows.truncate(1);
        draft.methods.mode = Mode::Overlap;
        draft.slots[0].captured = Some(capture);
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                start: value,
                ..
            } = &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 0;
        }
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                end: value,
                ..
            } = &mut draft.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = source.sequence.len();
        }
        draft.sync_rows();
        let draft = assembly_edit_model::declare(&draft, 0, start, end, inserted, None).unwrap();
        let role = &draft.methods.role_ids(1).unwrap()[0];
        let inputs = vec![PlanInput {
            id: role.input_id.clone(),
            revision: source.revision(),
            molecule: source,
        }];
        let intent = assembly_edit_intent::build(
            &inputs,
            &assembly_edit_intent::Request {
                base: collagene_core::cloning::assembly::substitution_intent::Request {
                    base: collagene_core::cloning::assembly::intent::Request {
                        material_bindings: Vec::new(),
                        name: draft.name.clone(),
                        circular: true,
                        fragments: vec![collagene_core::cloning::assembly::intent::Contribution {
                            role: role.role.clone(),
                            input_id: role.input_id.clone(),
                            selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                                start: 0,
                                end: 32,
                                reverse: false,
                            },
                        }],
                    },
                    junction_edits: vec![],
                    substitutions: vec![],
                },
                edits: draft.edits.clone(),
            },
            None,
        )
        .unwrap();
        (
            draft,
            intent,
            event,
            Options {
                forward_length: 8,
                reverse_length: 8,
                reverse_tail_bases: 0,
            },
        )
    }
    #[test]
    fn qualifier_aware_source_edit_can_prepare_inverse_pcr_without_reusing_stale_annotation_events()
    {
        use collagene_model::{Feature, Location};
        let (mut draft, old_intent, _, options) = fixture("4", "6", "GGA");
        let mut capture = draft.slots[0].captured.as_ref().unwrap().as_ref().clone();
        capture.molecule.features.push(Feature {
            kind: "tRNA".into(),
            location: Location::range(12, 15),
            qualifiers: vec![(
                "anticodon".into(),
                Some("(pos:13..15,aa:Arg,seq:acg)".into()),
            )],
            needs_review: false,
        });
        let original = capture.molecule.clone();
        let frame = Workspace::from_documents("Source annotation event", vec![original.clone()], 0)
            .unwrap();
        let document = &frame.documents[0];
        let event = EventBinding {
            document_id: document.id.clone(),
            event_id: document.head_event.clone().unwrap(),
            revision: document.head.clone(),
        };
        capture.source.binding = event.clone();
        draft.slots[0].captured = Some(Arc::new(capture));
        assert!(
            request(
                &draft,
                &old_intent,
                &event,
                &old_intent.target,
                &options,
                None
            )
            .is_err()
        );
        let (inputs, base, _) = draft.inputs_request(None).unwrap();
        let intent = assembly_edit_intent::build(&inputs, &draft.edit_request(base), None).unwrap();
        assert_eq!(
            intent.target.sequence.as_str(),
            "AGTGGAGATGCTACGTAGGACTTCAGACCTGTA"
        );
        assert_eq!(intent.target.features[0].location, Location::range(13, 16));
        assert_eq!(
            intent.target.features[0].qualifiers,
            vec![(
                "anticodon".into(),
                Some("(pos:14..16,aa:Arg,seq:acg)".into())
            )]
        );
        let proposal = propose(&draft, &intent, &event, &intent.target, &options, None).unwrap();
        assert_eq!(
            proposal.verification.status,
            Status::ConsistentUnderAssumptions
        );
        assert_eq!(proposal.candidate.plan.inputs[0].molecule, original);
        assert_eq!(draft.slots[0].captured.as_ref().unwrap().molecule, original);
    }
    #[test]
    fn original_source_inverse_pcr_reconstructs_independent_insertion_deletion_and_replacement() {
        for (start, end, inserted) in [("11", "11", "AACG"), ("10", "14", ""), ("4", "6", "GGA")] {
            let (draft, intent, event, mut options) = fixture(start, end, inserted);
            for split in 0..=inserted.len() {
                options.reverse_tail_bases = split;
                let proposal =
                    propose(&draft, &intent, &event, &intent.target, &options, None).unwrap();
                assert_eq!(
                    proposal.verification.status,
                    Status::ConsistentUnderAssumptions,
                    "{:?}",
                    proposal.verification.findings
                );
                assert_eq!(
                    proposal.candidate.plan.inputs[0].molecule,
                    draft.slots[0].captured.as_ref().unwrap().molecule
                );
                assert_eq!(
                    proposal.candidate.plan.inputs[0].id,
                    proposal.capture.input_id
                );
                assert_eq!(
                    proposal.candidate.plan.target_revision,
                    intent.target.revision()
                );
                assert_eq!(
                    proposal.verification.target_revision,
                    intent.target.revision()
                );
                // The diagnostic keeps the worker's original plan identity. It
                // cannot become a receipt for the stable-role candidate by rename.
                assert_eq!(
                    proposal.verification.plan_revision,
                    collagene_core::planning::plan::DesignPlan::revision(
                        &proposal.worker_design.plan
                    )
                );
                assert_ne!(
                    proposal.verification.plan_revision,
                    proposal.candidate.plan_revision
                );
                let rebound_check =
                    verification::engine::verify(&proposal.candidate.plan, &intent.target);
                assert_eq!(rebound_check.status, Status::ConsistentUnderAssumptions);
                assert_eq!(
                    rebound_check.plan_revision,
                    proposal.candidate.plan_revision
                );
                assert_eq!(proposal.verification.products.len(), 3);
                let actual = proposal.verification.products.last().unwrap();
                assert_eq!(actual.topology, Topology::Circular);
                assert!(
                    analysis::rotation_offset(
                        actual.sequence.as_str(),
                        intent.target.sequence.as_str()
                    )
                    .is_some()
                );
                proposal
                    .validate_binding(&draft, &intent, &event, &intent.target, &options, None)
                    .unwrap();
            }
        }
    }
    #[test]
    fn source_event_wrong_target_and_forged_receipt_are_rejected() {
        let (draft, intent, event, options) = fixture("11", "11", "AACG");
        let mut stale = event.clone();
        stale.event_id = "0".repeat(64);
        assert!(propose(&draft, &intent, &stale, &intent.target, &options, None).is_err());
        let mut wrong = intent.target.clone();
        wrong.sequence = collagene_model::Sequence::dna(&"A".repeat(wrong.sequence.len())).unwrap();
        assert!(propose(&draft, &intent, &event, &wrong, &options, None).is_err());
        wrong = intent.target.clone();
        wrong.metadata.insert("changed".into(), "yes".into());
        assert!(propose(&draft, &intent, &event, &wrong, &options, None).is_err());
        let mut forged = intent.clone();
        forged.applied[0].output_start += 1;
        assert!(propose(&draft, &forged, &event, &forged.target, &options, None).is_err());
        let mut proposal =
            propose(&draft, &intent, &event, &intent.target, &options, None).unwrap();
        proposal.verification.findings.push("forged receipt".into());
        assert!(
            proposal
                .validate_binding(&draft, &intent, &event, &intent.target, &options, None)
                .is_err()
        );
    }
    #[test]
    fn worker_response_preserves_original_receipts_and_rejects_wrong_native_fields_or_products() {
        let (draft, intent, event, options) = fixture("11", "11", "AACG");
        let (captured, source, native) =
            request(&draft, &intent, &event, &intent.target, &options, None).unwrap();
        assert_eq!(captured.target_revision, intent.target.revision());
        let report = mutagenesis::design_inverse_pcr(
            &source,
            &intent.target,
            draft.methods.rows[0]
                .run_origin
                .as_ref()
                .map(|p| p.as_ref()),
            &native,
        )
        .unwrap();
        let original = serde_json::to_value(&report.verification).unwrap();
        let proposal = from_design(
            &draft,
            &intent,
            &event,
            &intent.target,
            &options,
            report.clone(),
            None,
        )
        .unwrap();
        assert_eq!(
            serde_json::to_value(&proposal.verification).unwrap(),
            original
        );
        let mut wrong = report.clone();
        wrong.forward =
            collagene_core::sequence::analysis::Primer::new("Forged", "AAAAAAAA", "").unwrap();
        assert!(
            from_design(
                &draft,
                &intent,
                &event,
                &intent.target,
                &options,
                wrong,
                None
            )
            .is_err()
        );
        let mut wrong = report.clone();
        wrong.options.reverse_tail_bases = 1;
        assert!(
            from_design(
                &draft,
                &intent,
                &event,
                &intent.target,
                &options,
                wrong,
                None
            )
            .is_err()
        );
        let mut wrong = report.clone();
        wrong.verification.target_revision = "0".repeat(64);
        assert!(
            from_design(
                &draft,
                &intent,
                &event,
                &intent.target,
                &options,
                wrong,
                None
            )
            .is_err()
        );
        let mut wrong = report.clone();
        wrong.verification.products[0]
            .metadata
            .insert("plan_revision".into(), "0".repeat(64));
        assert!(
            from_design(
                &draft,
                &intent,
                &event,
                &intent.target,
                &options,
                wrong,
                None
            )
            .is_err()
        );
        let mut wrong = report;
        wrong.verification.engine_version = "other-engine".into();
        assert!(
            from_design(
                &draft,
                &intent,
                &event,
                &intent.target,
                &options,
                wrong,
                None
            )
            .is_err()
        );
        let mut changed = proposal;
        changed.candidate.plan.steps.clear();
        assert!(
            changed
                .validate_binding(&draft, &intent, &event, &intent.target, &options, None)
                .is_err()
        );
    }
    #[test]
    fn unsupported_declarations_never_become_guessed_preparations() {
        let (draft, intent, event, options) = fixture("11", "11", "AACG");
        let mut partial = draft.clone();
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                start: value,
                ..
            } = &mut partial.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 1;
        }
        partial.sync_rows();
        assert!(propose(&partial, &intent, &event, &intent.target, &options, None).is_err());
        let mut reverse = draft.clone();
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                reverse: value,
                ..
            } = &mut reverse.slots[0].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = true;
        }
        reverse.sync_rows();
        assert!(propose(&reverse, &intent, &event, &intent.target, &options, None).is_err());
        let mut multiple = draft.clone();
        multiple.edits.push(assembly_edit_intent::Edit {
            role: multiple.edits[0].role.clone(),
            start: 20,
            end: 21,
            expected: "T".into(),
            replacement: "C".into(),
        });
        multiple.sync_rows();
        assert!(propose(&multiple, &intent, &event, &intent.target, &options, None).is_err());
        let mut mixture = draft.clone();
        mixture.substitutions.push(
            collagene_core::cloning::assembly::substitution_intent::Substitution {
                role: mixture.edits[0].role.clone(),
                offset: 20,
                expected: "T".into(),
                replacement: "C".into(),
            },
        );
        mixture.sync_rows();
        assert!(propose(&mixture, &intent, &event, &intent.target, &options, None).is_err());
        let mut multifragment = draft.clone();
        multifragment.slots.push(multifragment.slots[0].clone());
        multifragment.methods.source_roles = None;
        multifragment.sync_rows();
        assert!(
            propose(
                &multifragment,
                &intent,
                &event,
                &intent.target,
                &options,
                None
            )
            .is_err()
        );
    }
    #[test]
    fn existing_binary_worker_uses_native_context_and_returns_bound_preparation() {
        use collagene_core::protocol::{
            self,
            request::Request,
            schema::{AnalysisJob, AnalysisResult},
        };
        let (draft, intent, event, options) = fixture("11", "11", "AACG");
        let (capture, source, native) =
            request(&draft, &intent, &event, &intent.target, &options, None).unwrap();
        let original = source.clone();
        let request = Request {
            id: 721,
            revision: mutagenesis::inverse_pcr_revision(
                &source,
                &intent.target,
                capture.origin.as_ref().map(|p| p.as_ref()),
                &native,
            )
            .unwrap(),
            job: AnalysisJob::InversePcr {
                molecule: source,
                origin: capture.origin.as_ref().map(|p| p.as_ref().clone()),
                target: intent.target.clone(),
                options: native,
            },
        };
        let raw = protocol::execution::process_request(
            &protocol::request::Request::encode(&request).unwrap(),
        )
        .unwrap();
        let response = protocol::response::Response::decode(&raw).unwrap();
        assert_eq!(response.id, 721);
        assert_eq!(response.revision, request.revision);
        let AnalysisResult::InversePcr(report) = response.result else {
            panic!("inverse-PCR reply");
        };
        let proposal = from_design(
            &draft,
            &intent,
            &event,
            &intent.target,
            &options,
            report,
            None,
        )
        .unwrap();
        let fresh = verification::engine::verify(&proposal.candidate.plan, &intent.target);
        assert_eq!(
            fresh.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            fresh.findings
        );
        assert_eq!(proposal.candidate.plan.inputs[0].molecule, original);
        assert_eq!(proposal.capture.target_revision, intent.target.revision());
        let mut stale = request;
        stale.revision = intent.target.revision();
        assert_ne!(stale.revision, event.revision);
        assert!(
            protocol::execution::process_request(
                &protocol::request::Request::encode(&stale).unwrap()
            )
            .is_err()
        );
    }
    #[test]
    fn point_substitution_keeps_its_original_intent_through_binary_worker_and_rebound_review() {
        use collagene_core::protocol::{
            self,
            request::Request,
            schema::{AnalysisJob, AnalysisResult},
        };
        let (mut draft, _, event, mut options) = fixture("11", "11", "AACG");
        draft.edits.clear();
        draft.sync_rows();
        let draft =
            crate::assembly::edit::assembly_substitution_model::declare(&draft, 0, "11", "A", None)
                .unwrap();
        let source = draft.slots[0].captured.as_ref().unwrap().molecule.clone();
        let role = draft.methods.role_ids(1).unwrap().remove(0);
        let input = PlanInput {
            id: role.input_id.clone(),
            revision: source.revision(),
            molecule: source.clone(),
        };
        let request = assembly_substitution_intent::Request {
            base: collagene_core::cloning::assembly::intent::Request {
                material_bindings: Vec::new(),
                name: draft.name.clone(),
                circular: true,
                fragments: vec![collagene_core::cloning::assembly::intent::Contribution {
                    role: role.role.clone(),
                    input_id: role.input_id.clone(),
                    selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                        start: 0,
                        end: source.sequence.len(),
                        reverse: false,
                    },
                }],
            },
            junction_edits: vec![],
            substitutions: draft.substitutions.clone(),
        };
        let intent =
            assembly_substitution_intent::build(std::slice::from_ref(&input), &request, None)
                .unwrap();
        let target_revision = intent.target.revision();
        assert_eq!(
            intent.target.sequence.as_str(),
            "AGTCCGATGCAACGTAGGACTTCAGACCTGTA"
        );
        let typed = crate::assembly::review::assembly_review::IntendedReference::Substitutions(
            intent.clone(),
        );
        for split in 0..=1 {
            options.reverse_tail_bases = split;
            let (capture, original, native) =
                super::request(&draft, &typed, &event, &intent.target, &options, None).unwrap();
            assert_eq!(
                (native.edit_start, native.edit_end, native.inserted.as_str()),
                (10, 11, "A")
            );
            assert_eq!(capture.intent_context_revision, intent.context_revision);
            let binary = Request {
                id: 723 + split as u32,
                revision: mutagenesis::inverse_pcr_revision(
                    &original,
                    &intent.target,
                    capture.origin.as_ref().map(|p| p.as_ref()),
                    &native,
                )
                .unwrap(),
                job: AnalysisJob::InversePcr {
                    molecule: original,
                    origin: capture.origin.as_ref().map(|p| p.as_ref().clone()),
                    target: intent.target.clone(),
                    options: native,
                },
            };
            let response = protocol::response::Response::decode(
                &protocol::execution::process_request(
                    &protocol::request::Request::encode(&binary).unwrap(),
                )
                .unwrap(),
            )
            .unwrap();
            assert_eq!(response.id, binary.id);
            assert_eq!(response.revision, binary.revision);
            let AnalysisResult::InversePcr(report) = response.result else {
                panic!("inverse-PCR SNP reply")
            };
            let original_worker_plan = report.verification.plan_revision.clone();
            let proposal = from_design(
                &draft,
                &typed,
                &event,
                &intent.target,
                &options,
                report,
                None,
            )
            .unwrap();
            assert_eq!(proposal.verification.plan_revision, original_worker_plan);
            assert_ne!(
                proposal.verification.plan_revision,
                proposal.candidate.plan_revision
            );
            assert_eq!(proposal.candidate.plan.inputs[0].molecule, source);
            assert_eq!(proposal.candidate.plan.inputs[0].id, role.input_id);
            assert_eq!(proposal.candidate.plan.target_revision, target_revision);
            proposal
                .validate_binding(&draft, &typed, &event, &intent.target, &options, None)
                .unwrap();
            // The bridge has not converted the substitution target into another
            // design-intent schema or changed its metadata to fit a PCR result.
            assert_eq!(typed.target().revision(), target_revision);
            assert_eq!(
                super::capture(&draft, &intent, &event, &intent.target, &options, None)
                    .unwrap()
                    .intent_context_revision,
                intent.context_revision
            );
            let fresh = verification::engine::verify(&proposal.candidate.plan, &intent.target);
            assert_eq!(
                fresh.status,
                Status::ConsistentUnderAssumptions,
                "{:?}",
                fresh.findings
            );
            let actual = fresh.products.last().unwrap();
            assert!(
                analysis::rotation_offset(
                    actual.sequence.as_str(),
                    intent.target.sequence.as_str()
                )
                .is_some()
            );
        }
        let mut mixed = draft.clone();
        mixed.edits.push(assembly_edit_intent::Edit {
            role: role.role.clone(),
            start: 20,
            end: 21,
            expected: "T".into(),
            replacement: "C".into(),
        });
        mixed.sync_rows();
        assert!(super::request(&mixed, &typed, &event, &intent.target, &options, None).is_err());
        let multiple =
            crate::assembly::edit::assembly_substitution_model::declare(&draft, 0, "21", "C", None)
                .unwrap();
        assert!(super::request(&multiple, &typed, &event, &intent.target, &options, None).is_err());
        let mut forged = intent.clone();
        forged.applied[0].output_position += 1;
        assert!(super::request(&draft, &forged, &event, &intent.target, &options, None).is_err());
        let mut stale = event;
        stale.event_id = "0".repeat(64);
        assert!(super::request(&draft, &typed, &stale, &intent.target, &options, None).is_err());
    }
    #[test]
    fn saved_run_and_derived_circle_proposals_retain_exact_origin_after_restart() {
        use collagene_core::{
            chemistry::Terminus,
            material::{
                lineage::{self as run_lineage, InputEventBinding},
                run::{self as run_material, MaterialOrigin},
            },
            planning::plan::OrientedSource,
        };
        let (base, _, _, options) = fixture("11", "11", "AACG");
        let circle = base.slots[0].captured.as_ref().unwrap().molecule.clone();
        let mut linear = circle.clone();
        linear.topology = Topology::Linear;
        let frame =
            Workspace::from_documents("Native circles", vec![linear.clone(), circle.clone()], 1)
                .unwrap();
        let pin = |i: usize| {
            let d = &frame.documents[i];
            EventBinding {
                document_id: d.id.clone(),
                event_id: d.head_event.clone().unwrap(),
                revision: d.head.clone(),
            }
        };
        let plan = collagene_core::planning::plan::DesignPlan {
            name: "Explicitly close circle".into(),
            target_revision: circle.revision(),
            inputs: vec![PlanInput {
                id: "stock".into(),
                revision: linear.revision(),
                molecule: linear,
            }],
            steps: vec![
                Operation::DeclareSourceTermini {
                    template: Source::Input("stock".into()),
                    left: Terminus::CLEAVED,
                    right: Terminus::CLEAVED,
                },
                Operation::Ligate {
                    fragments: vec![OrientedSource {
                        source: Source::Step(0),
                        reverse: false,
                    }],
                    circular: true,
                    allow_nicked: false,
                    options: Default::default(),
                },
            ],
        };
        let kept = run_lineage::record_run(
            &frame,
            &plan,
            pin(1),
            vec![InputEventBinding {
                input_id: "stock".into(),
                event: pin(0),
            }],
            false,
        )
        .unwrap();
        let row = kept.receipt.products.last().unwrap();
        let event = EventBinding {
            document_id: row.document_id.clone(),
            event_id: row.event_id.clone(),
            revision: row.revision.clone(),
        };
        let run = collagene_core::material::capsule::capsule(
            &kept.workspace,
            &run_material::Selection {
                run_id: kept.receipt.run_id,
                step_index: 1,
                event: event.clone(),
            },
        )
        .unwrap();
        let origin = MaterialOrigin::Run(run);
        let mut saved = kept.workspace;
        saved.active_document = event.document_id.clone();
        let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
        let selection =
            collagene_core::cloning::restriction::catalog::select(&catalog, &["ecori".into()])
                .unwrap();
        let derived = collagene_core::material::derivation::install(
            &saved,
            &collagene_core::material::derivation::Request {
                sources: vec![collagene_core::material::derivation::BoundSource {
                    role: collagene_core::material::derivation::SourceRole::Template,
                    event: event.clone(),
                    origin: collagene_core::material::derivation::SourceOrigin::Saved(
                        origin.clone(),
                    ),
                }],
                name: "Derived uncut circle".into(),
                save_primers: false,
                derivation: collagene_core::material::derivation::Derivation::CatalogDigest {
                    conditions: None,
                    selection,
                    options: Default::default(),
                    fragment_index: 0,
                },
            },
        )
        .unwrap();
        let derived_frame =
            collagene_core::documents::workspace::load_workspace(&derived.archive).unwrap();
        let derived_event = derived.receipt.product.clone();
        let derived_origin = MaterialOrigin::Derived(
            collagene_core::material::derivation::capsule(&derived_frame, &derived_event).unwrap(),
        );
        for (workspace, event, origin) in [
            (saved, event, origin),
            (derived_frame, derived_event, derived_origin),
        ] {
            let reopened = collagene_core::documents::workspace::load_workspace(
                &collagene_core::documents::workspace::save_workspace(&workspace).unwrap(),
            )
            .unwrap();
            origin.validate_for_workspace(&reopened, &event).unwrap();
            let before = serde_json::to_vec(&reopened).unwrap();
            let molecule = reopened.revisions[&event.revision].clone();
            let mut draft = base.clone();
            draft.slots[0].captured = Some(Arc::new(CapturedSource {
                source: SourceInfo {
                    binding: event.clone(),
                    name: molecule.name.clone(),
                    length: molecule.sequence.len(),
                    alphabet: Alphabet::Dna,
                    topology: molecule.topology,
                },
                molecule,
            }));
            draft.methods.rows[0].run_origin = Some(
                crate::history::run_origin_model::PinnedOrigin::from_material(origin.clone())
                    .unwrap(),
            );
            draft.slots[0].material_origin = draft.methods.rows[0].run_origin.clone();
            draft.sync_rows();
            let (inputs, base_request, _) = draft.inputs_request(Some(&reopened)).unwrap();
            let intent = assembly_edit_intent::build(
                &inputs,
                &draft.edit_request(base_request),
                Some(&reopened),
            )
            .unwrap();
            assert!(propose(&draft, &intent, &event, &intent.target, &options, None).is_err());
            let proposal = propose(
                &draft,
                &intent,
                &event,
                &intent.target,
                &options,
                Some(&reopened),
            )
            .unwrap();
            assert!(proposal.capture.origin.is_some());
            assert_eq!(proposal.candidate.plan.steps.len(), 4);
            assert!(
                matches!(&proposal.candidate.plan.steps[0],Operation::ReconstructRunProduct{template:Source::Input(id),..}|Operation::ReconstructMaterial{template:Source::Input(id),..} if id==&proposal.capture.input_id)
            );
            assert!(matches!(
                proposal.candidate.plan.steps[1],
                Operation::Pcr {
                    template: Source::Step(0),
                    ..
                }
            ));
            assert_eq!(
                verification::engine::verify(&proposal.candidate.plan, &intent.target).status,
                Status::ConsistentUnderAssumptions
            );
            let mut alias = event.clone();
            alias.event_id = "equal-DNA-wrong-event".into();
            assert!(
                propose(
                    &draft,
                    &intent,
                    &alias,
                    &intent.target,
                    &options,
                    Some(&reopened)
                )
                .is_err()
            );
            let mut changed = draft.clone();
            changed.methods.rows[0].run_origin = None;
            changed.slots[0].material_origin = changed.methods.rows[0].run_origin.clone();
            assert!(
                proposal
                    .validate_binding(
                        &changed,
                        &intent,
                        &event,
                        &intent.target,
                        &options,
                        Some(&reopened)
                    )
                    .is_err()
            );
            let mut wrong = intent.target.clone();
            wrong.sequence = collagene_model::Sequence::dna("TTTTCCCCAAAAGGGG").unwrap();
            assert!(propose(&draft, &intent, &event, &wrong, &options, Some(&reopened)).is_err());
            assert_eq!(serde_json::to_vec(&reopened).unwrap(), before);
        }
    }
}
